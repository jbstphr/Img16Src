// I16_Ut5c.cpp : quad-laps' treatment (CPU-only)
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut5.h"
#include "BS_Ut7.h" // ad-hoc service CPU-transforms
#include "BS_Ut5c.h"

//// miscelanneous commons ////  
static int ixQuadLapsEnabled = 0; int BS5c_QuadLapsEnabled() { return ixQuadLapsEnabled; } // -- >0:is enabled

static struct {
  int rLapOrgThr;  // applies at beginning; if |LapOrg[i]|>ThrLapOrg then 
                   //  "constraint LapOrg[i] by ThrLapOrg", "set all QuaLaps#e[i] to 0"
  int rCorrMetMn,  // applies at finalizing;
      rCorrMetGa;  //  "constraint |CorrMet| within |Mn|..|Mn+Ga|"
} ixQL = {
  300,            // LapOrgThr
   99,            // CorrMetMn|minimum
  222             // CorrMetGa|gamut
};

void BS5c_QuadLapsLoad( GIniFile *AF ) { // -- load estimation parameters
  ixQuadLapsEnabled = atoi( AF->SGetStr( "nEnabled", "0" ) );
  if (!ixQuadLapsEnabled) return;
  const char *cS; char qS[256];
  cS = AF->DGetStr( "LapOrg", NULL ); // treatment of LapOrg at beginning
  if (!cS || !sscanf( cS, "%d", &ixQL.rLapOrgThr )) {
    sprintf( qS, "%d;Thr", ixQL.rLapOrgThr ); AF->SetStr( "LapOrg", qS ); };
  cS = AF->DGetStr( "CorrMet", NULL ); // treatment of CorrMet when finalizing
  if (!cS || !sscanf( cS, "%d %d", &ixQL.rCorrMetMn, &ixQL.rCorrMetGa )) {
    sprintf( qS, "%d %d;Mn Ga", ixQL.rCorrMetMn, ixQL.rCorrMetGa ); AF->SetStr( "CorrMet", qS ); };
}

tBS_Val **BS5c_QuadLapsAlloc( int AW, int AH, tBS_Val **ABufs, // -- pointer to allocated buffers
  int ARO, int ARQ ) { // implementation-specific params
  BS_DatasAlloc( cBS5c_QuadLaps, AW*AH, ABufs ); return ABufs; }

#if 1 /// production !!!!
// input & output:
//   inp Lap:LapOrg, Laps:8 quad-laps (assume anisotropically blurred)
//   out Lap:~"mx correlation" between LapOrg & Laps#e
// elementary process-flow:
//   step 1) Encode inp LapOrg and Laps#e into few bits so that the integrals
// (along intervals of const sign) would be ~small (i.e. within +|-"Length*10").
//   step 2) Couple-pack each encoded Laps#e together with encoded LapOrg (this
// is easy with few bits).
//   step 3) For each couple-pack, integrate relevant metrics along the
// intervals of const sign.
//   step 4) For each pixel, choose max of 8 metrics. This stage outputs the max
// among aliased "correlation metrics" between LapOrg and Laps#e.

struct tQuadWr {
  char     rDbgPfx[256], *rDbgSfx; // debug-dump's strings
  int      rW, rH; // frame's dimension
  tBS_Val *rLap; // as the input it could be LapMid|LapAda|cLapMM; it is resulting laplacian-like mark
  tBS_Val *rLaps[cBS5c_QuadLaps]; // input quad-laps (mostly unchanged)
  tBS_Val *rEncs[9]; // encoded frame-buffers (memory gap for skew)

  // -- init & release
  void DbgSet( const char *ADbgPfx ) { // debug-trace
    if (ADbgPfx) rDbgSfx = rDbgPfx+sprintf( rDbgPfx, "%s_", ADbgPfx ); }
  void Init( int AW, int AH, tBS_Val **ALaps, bool ACUDA, const tBS_Val *ALapOrg, const tBS_Val* ) { // fills "self"
    rW = AW; rH = AH; const int cN = AW*AH; rLap = new tBS_Val[cN];
    if (ACUDA) BC_Move( -cN*cBS_Sz, rLap, (tBS_Val*)ALapOrg );
          else BS_Cpy( cN, ALapOrg, rLap );
    memcpy( rLaps, ALaps, cBS5c_QuadLaps*sizeof(void*) );
    BS_DatasAlloc( mcASz(rEncs), cN, rEncs ); }
  void Release() { // releases "self" (except the results and org inps)
    BS_DatasRelease( mcASz(rEncs), rEncs );
    BS_DatasRelease( cBS5c_QuadLaps, rLaps ); } // caller allocated these buffers by BS5c_QuadLapsAlloc()

  // -- encode inp quad-laps to be short & simply countable;
  // -- then pack pixels of Lap & Laps#i into single BS_Val.
  static void itppr_LapEnc( tBS_SpawnThreadRec *AP ) { // neg_pixel~>"-1", pos_pixel~>"+1"
    tQuadWr *cThis = (tQuadWr*)AP->rP[0]; const int cN = cThis->rW*cThis->rH;
    const int cEncodeShr = 6, cEncodeMx = 3; for (int qe=AP->rFst;qe<AP->rLst;qe++) {
      tBS_Val *qLap = cThis->rEncs[qe]; for (int qi=0;qi<cN;qi++) { int qV = qLap[qi]; if (!qV) continue;
      int qA = abs(qV)>>cEncodeShr; if (qA>cEncodeMx) qA = cEncodeMx; qLap[qi] = qV<0?-qA-1:qA+1; };
    }; }
  static void itppr_LapPck( tBS_SpawnThreadRec *AP ) { // magnitude & sign by "Laps#e", low bits keep "Lap"
    tQuadWr *cThis = (tQuadWr*)AP->rP[0]; const int cN = cThis->rW*cThis->rH;
    const tBS_Val *cLap = cThis->rEncs[8]; for (int qe=AP->rFst;qe<AP->rLst;qe++) {
      tBS_Val *qLap = cThis->rEncs[qe]; for (int qi=0;qi<cN;qi++) {
      int qV = qLap[qi]; if (qV) qLap[qi] = ((qV<<8)&~255)|(cLap[qi]&255);
    }; }; }
  void iEncAndPckLaps() { // encode and pack
    BS_SpawnThreads( 0, 0, 9, itppr_LapEnc, this ); // encode inp quad-laps and LapOrg
    BS_SpawnThreads( 0, 0, 8, itppr_LapPck, this ); } // pack encoded quad-laps with LapOrg

  // -- count correlations between Lap & Laps#e along relevant directions;
  // -- then estimate correlation metric for each pixel (~as simple maximum).
  static int icpl_CorrelCount( int AL, tBS_Val *ADat ) { // count correlation-metric along a line
    int qL = 0, qZq = 0, qZo = 0; ADat[AL] = 0; // simplify loop's condition
    for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; 
      if (qV>0) { if (qL<0) goto lab1; qL++; qZq += abs(qV>>8); qZo += (char)qV; continue; };
      if (qV<0) { if (qL>0) goto lab1; qL--; qZq += abs(qV>>8); qZo += (char)qV; continue; };
      if (!qL) continue; lab1:
      qL = abs( qL ); int qZoa = abs(qZo);
      int //qEmpZqMn = qZq>>2, qEmpZqMx = qZq+(qZq>>1), // empirical constraints
          //qCorr = (qZoa<qEmpZqMn)||(qZoa>qEmpZqMx)?0:(qZoa+qZq)*mcRng(qZo,-1,1);
          qCorr = (qZoa+qZq)*mcRng(qZo,-1,1);
      BS_VFil( qL, ADat+qi-qL, qCorr );
      qL = mcRng( qV, -1, 1 ); qZq = abs(qV>>8); qZo = (char)qV;
    }; return AL; } // for

  void iCountCorrel( tBS_Val **AEnc8, tBS_Val *ARes, tBS_Val **ATLimbo8 ) { // estimate correlation-metric of Ecns#e
    const int cN = rW*rH; for (int qi=0;qi<8;qi++) BS_Cpy( cN, AEnc8[qi], ATLimbo8[qi] );
    BS_CPLA( icpl_CorrelCount, 8, cBS_AllDirsT, rW, rH, ATLimbo8 ); // count correlation-metric along a line
    //strcpy( rDbgSfx, "bqPacks" ); BS_DatasWrt( rDbgPfx, cN, qLimbo, 0, 8 );
    BS_CPB( mccpb_MxA<8>, 8|cCPBNDD, cN, ATLimbo8, ARes ); } // mx of correclations

  // -- process
  void QuadLapsProcess() {
    const int cW = rW, cH = rH, cN = cW*cH;
    { // copy the inputs to encode-buffers
      for (int qe=0;qe<8;qe++) BS_Cpy( cN, rLaps[qe], rEncs[qe] );
      tBS_Val *qLapOrg = rLap, *qLapOrgEnc = rEncs[8]; // set Encs#0..7 to zero when 
      const int cLapOrgThr = ixQL.rLapOrgThr;          // |Encs#8]|>ThrLapOrg; constraint Encs#8; 
      for (int qi=0;qi<cN;qi++) { int qV = qLapOrg[qi];
        if (abs(qV)<cLapOrgThr) qLapOrgEnc[qi] = qV; else {
          qLapOrgEnc[qi] = qLapOrg[qi] = qV<0?-cLapOrgThr:cLapOrgThr;
          for (int qe=0;qe<8;qe++) rLaps[qe][qi] = 0;
        };
      };
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "aInps" ); BS_DatasWrt( rDbgPfx, cN, rEncs, 0, 9 ); };
    };
    { // encode and pack
      iEncAndPckLaps(); 
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "bPcks" ); BS_DatasWrt( rDbgPfx, cN, rEncs, 0, 9 ); }; };
    // estimate correlation-metric
    tBS_Val    *qMetric = rLaps[7], // output estimated plain metric to ...
               *qMetSkw = rEncs[8]; // output estimated skewed metric to ...
    tBS_Val    *qTLimbo[8]; BS_DatasAlloc( 8, cN, qTLimbo ); // allocate temp_limbo
    { // estimate metric of plain Ecns
      iCountCorrel( rEncs, qMetric, qTLimbo ); // count the correlations
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "caLimbo" ); BS_DatasWrt( rDbgPfx, cN, qTLimbo ); };
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "cbMetric.raw" ); BS_Write( rDbgPfx, cN, qMetric ); }; };
    { // estimate metric of skewed Ecns
      #if 1
        const int cShearM = 7, cShears[][2] = { {-cShearM*2,-cShearM*2}, {-cShearM,-cShearM} };
      #else
        const int cShearM = 5, cShears[][2] = { {-cShearM*3,-cShearM*3}, {-cShearM*2,-cShearM*2}, {-cShearM,-cShearM} };
      #endif
      const int cWskw = BS_Align( cW+BS7_SkewShear(cShearM,cH)*2+1 ), // skew buffer's gap by X
                cHskw = BS_Align( cH+BS7_SkewShear(cShearM,cW)*2+1 ); // skew buffer's gap by Y
      tBS_Val  *qEnks[8]; BS_DatasAlloc( 8, cWskw*cHskw, qEnks ); // allocate temp_encs
      for (int qskw=0;qskw<mcASz(cShears);qskw++) { // estimate metric of skewed Ecns
        const int cShearH = cShears[qskw][0], cShearV = cShears[qskw][1];
        for (int qe=0;qe<8;qe++) { BS_Cpy( cN, rEncs[qe], qEnks[qe] );
          BS7_SkewH( cShearH, false, qEnks[qe], cW, cH ); BS7_SkewV( cShearV, false, qEnks[qe], cW, cH ); }; // skew
        if (0 && rDbgSfx) { sprintf( rDbgSfx, "d%daPcks", qskw ); BS_DatasWrt( rDbgPfx, cN, qEnks ); };
        iCountCorrel( qEnks, qMetSkw, qTLimbo ); // count the correlations
        if (0 && rDbgSfx) { sprintf( rDbgSfx, "d%dbLimb", qskw ); BS_DatasWrt( rDbgPfx, cN, qTLimbo ); }
        BS7_SkewV( -cShearV, false, qMetSkw, cW, cH ); BS7_SkewH( -cShearH, false, qMetSkw, cW, cH ); // deskew
        if (1 && rDbgSfx) { sprintf( rDbgSfx, "d%dcMetr.raw", qskw ); BS_Write( rDbgPfx, cN, qMetSkw ); };
        for (int qi=0;qi<cN;qi++) { // accumulate the max metric
          int qV = qMetSkw[qi]; if (abs(qMetric[qi])<abs(qV)) qMetric[qi] = qV; }; };
      BS_DatasRelease( 8, qEnks ); }; // release temp_encs
    BS_DatasRelease( 8, qTLimbo ); // release temp_limbo
    { // constraint combined metric
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "eaMetCns.raw" ); BS_Write( rDbgPfx, cN, qMetric ); };
      const int cCorrMetMn = ixQL.rCorrMetMn, cCorrMetGa = ixQL.rCorrMetGa;
      for (int qi=0;qi<cN;qi++) { int qV = qMetric[qi], qA = abs(qV)-cCorrMetMn;
        if (qA<=0) { qMetric[qi] = 0; continue; }; if (qA>cCorrMetGa) qA = cCorrMetGa;
        qMetric[qi] = qV<0?-qA:qA; };
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "ebMetCns.raw" ); BS_Write( rDbgPfx, cN, qMetric ); }; };
    { // result as "Met0+Met1*LapOrg"
      tBS_Val  *qMet0 = qMetric, *qMet1 = qMetSkw; BS_Cpy( cN, qMet0, qMet1 ); BS_Abs( cN, qMet1 );
      const int cRad0 =  1; BS2_MovAve2Dp( cRad0%100, qMet0, rLaps, cRad0/100 );
      const int cRad1 =102; BS2_MovAve2Dp( cRad1%100, qMet1, rLaps, cRad1/100 );
      if (1 && rDbgSfx) {
        strcpy( rDbgSfx, "fInp.raw" ); BS_Write( rDbgPfx, cN, rLap );
        strcpy( rDbgSfx, "fMet0.raw" ); BS_Write( rDbgPfx, cN, qMet0 );
        strcpy( rDbgSfx, "fMet1.raw" ); BS_Write( rDbgPfx, cN, qMet1 ); }
      for (int qi=0;qi<cN;qi++) rLap[qi] = qMet0[qi]*2+(rLap[qi]*qMetSkw[qi]>>5);
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "fOut.raw" ); BS_Write( rDbgPfx, cN, rLap ); } };

      #if 0  // !!!! test blur
      mcShiftRgt( cN, rLap, 5 ); mcShiftLft( cN, rLap, 4 ); 
      BS_Rng( cN, rLap, -99, 99 );
      BS2_MovAve2Dp( 3, rLap, rEncs, 0 );
      #endif
  }; // end of QuadLapsProcess
}; // end of tQuadWr

#else //// test ???

struct tQuadWr { // class for  test 
  char     rDbgPfx[256], *rDbgSfx; // debug-dump's strings
  int      rW, rH; // frame's dimension
  tBS_Val *rLap, // as the input it could be LapMid|LapAda|cLapMM; it is resulting laplacian-like mark
          *rMsk;
  tBS_Val *rLaps[cBS5c_QuadLaps]; // input quad-laps (mostly unchanged)
  // -- init & release
  void DbgSet( const char *ADbgPfx ) { // debug-trace
    if (ADbgPfx) rDbgSfx = rDbgPfx+sprintf( rDbgPfx, "%s_", ADbgPfx ); }
  void Init( int AW, int AH, tBS_Val **ALaps, bool ACUDA, const tBS_Val *ALapOrg, const tBS_Val *AMsk ) { // fills "self"
    rW = AW; rH = AH; const int cN = AW*AH;
    rLap = new tBS_Val[cN]; rMsk = new tBS_Val[cN];
    if (ACUDA) {
      BC_Move( -cN*cBS_Sz, rLap, (tBS_Val*)ALapOrg ); BC_Move( -cN*cBS_Sz, rMsk, (tBS_Val*)AMsk );
    } else {
      BS_Cpy( cN, ALapOrg, rLap ); BS_Cpy( cN, AMsk, rMsk );
    }
    memcpy( rLaps, ALaps, cBS5c_QuadLaps*sizeof(void*) );
  }
  void Release() { // releases "self" (except the results and org inps)
    delete rMsk;
    BS_DatasRelease( cBS5c_QuadLaps, rLaps ); } // caller allocated these buffers by BS5c_QuadLapsAlloc()

  void QuadLapsProcess() {
    const int cW = rW, cH = rH, cN = cW*cH;
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "aInps" ); BS_DatasWrt( rDbgPfx, cN, rLaps, 0, 8 ); };
    BS_CPB( mccpb_MxA<8>, 8|cCPBNDD, cN, rLaps, rLap ); // mx of correclations
    BS_Abs( cN, rLap );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "bMax.raw" ); BS_Write( rDbgPfx, cN, rLap ); };
    BS_Mul( cN, rLap, rMsk, -9 ); 
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "cMul.raw" ); BS_Write( rDbgPfx, cN, rMsk ); };
  }

}; // end of tQuadWr test

#endif

static tQuadWr *jxQuWrSpawn1; static void jQuWrSpawn1( int ATag ) { // could be inline or spawned (via the only exemplar of tAiWr)
  if (ATag!=1) return; jxQuWrSpawn1->QuadLapsProcess(); }

void BS5c_QuadLaps( int AW, int AH, tBS_Val **ALaps, bool ACUDA, // forks quad-laps thread
  const char *ADbgPfx, const tBS_Val *ALapOrg, const tBS_Val *AMsk ) { 
  static tQuadWr qQuWr; memset( &qQuWr, 0, sizeof(tQuadWr) ); // this static locates in DSeg
  qQuWr.DbgSet( ADbgPfx ); 
  qQuWr.Init( AW, AH, ALaps, ACUDA, ALapOrg, AMsk );
  jxQuWrSpawn1 = &qQuWr; BS_Spawn1( jQuWrSpawn1, 2 ); } //  inline execution "jQuWrSpawn1( 1 )"

tBS_Val *BS5c_QuadLapsWait() { // wait for forked quad-laps thread
  BS_Spawn1Wait( 2 ); jxQuWrSpawn1->Release(); return jxQuWrSpawn1->rLap; }
