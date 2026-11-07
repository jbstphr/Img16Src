// in might be includes into "stdafx.h", into "kernel.cu"
#ifndef I16_INCV_H_INCLUDED
#define I16_INCV_H_INCLUDED

  #define mcASz( AnArray ) ((int)(sizeof(AnArray)/sizeof(AnArray[0])))
  #define mcSwap( tN, A, B ) { tN T=(tN)(A); A=(tN)(B); B=T; }
  #define mcRng( V, Min, Max ) (((V)<(Min))?(Min):(((V)>(Max))?(Max):(V)))
  #define mcRg( V ) mcRng( V, -0x7FFF, 0x7FFF )
  #define mcIpl(X,IX,Y0,Y1) ((Y0)+((Y1)-(Y0))*(X)/(IX)) // linear interpolation
  #define mcVi( V ) (int)(size_t)(V) // pointer_container to int
  #define mcIv( I ) (void*)(size_t)(I) // int to pointer_container
  #define strend( S ) (S+strlen(S)) // EoS
  #define cTPS "gc12" // project's message_prefix
  typedef short tBS_Val; const int cBS_Sz = sizeof(tBS_Val); // working signed pixels
  typedef unsigned short tBS_Card; // rather use (instead of int) for cardinal's dbg output
  #ifdef _WIN32
    typedef __int64 tLL;
    typedef unsigned __int64 tLLU;
  #else
    typedef long long tLL;
    typedef unsigned long long tLLU;
    #define __min(a,b) (a<b?a:b)
    #define __max(a,b) (a<b?b:a)
    #define __fastcall
    #define __cdecl
  #endif

  const float cPi = 3.14159265359f,
              cPitch45 = 1.4142135624f, // pixel's distance for 45 deg
              cPitch22 = 1.0823922024f; // pixel's distance for 22.5 deg 
  
  // == the Commons == dispatched in VAPrcFlu.cpp
  void WriteLn( const char *AMsg ); // WriteLn() to ~a console
  void ErrMsg( const char *AMsg ); // WriteLn( "ERR " Msg )
  extern int xDbgDmp; // extern char *xDbgPfx, *xDbgSfx; // dbg_output_label is 'strcat(DbgSfx,...)+printf(DbgPfx)'
  extern int xN, xModeId; // the parameters of I16_Init, I16_Invoke
  extern int xPars[]; // defined by _ParSet
  const int  cParCUDA       = 0, // 0:"CPU", 1:"CUDA"
             cParIniSetDflt = 1; // 1:IniFile.SGetStr overwrites certains section<s>
  
  #ifdef __CUDACC__ // HostAndDevice for NVCC
    #define HSTnDEV __host__ __device__
  #else // a module's intrinsic for cl
    #define HSTnDEV static
  #endif

  // 1D-HermiteLike-spline
  const int cHerLikFo = 0xFFFF, // fill by nearest known origin (or explicit short value within this mask)
            cHerLikFs = 0x10000; // sharp spline
  inline HSTnDEV bool mcHerLikBuf4( int AFill, int AS, int AL, const tBS_Val *ASrc, int Ao, int *ABuf4 ) { // -- fill the buffer
    // Fill contains cHerLikF...; S is stride; L is length
    if (Ao<1) { ABuf4[0] = Ao<0 ? ((AFill&cHerLikFo)==cHerLikFo? ASrc[0] : (short)AFill)
                                : (ASrc[0]+ASrc[AS])>>1; return false; };
    if (Ao>=AL-2) { int qo = (AL-1)*AS;
      ABuf4[0] = Ao>=AL-1 ? ((AFill&cHerLikFo)==cHerLikFo ? ASrc[qo] : (short)AFill)
                          : (ASrc[qo]+ASrc[qo-AS])>>1; return false; };
    Ao = (Ao-1)*AS; for (int qi=0;qi<4;qi++) { ABuf4[qi] = ASrc[Ao]; Ao += AS; }; return true; }
  // 1D-HermiteLike-splines; usually tA=int; 2nd stage of bicubic_spline uses tA=float
  template <class tA> inline HSTnDEV float mcHerLikB3( float AT, const tA *A4 ) { // -- wiki's bicubic spline; T is offset against A[1]
    tA m0 = A4[1]-A4[0], m1 = A4[3]-A4[2]; return (((2*A4[1]+m0-2*A4[2]+m1)*AT-3*A4[1]+3*A4[2]-2*m0-m1)*AT+m0)*AT+A4[1]; }
  template <class tA> inline HSTnDEV float mcHerLikPB( float AT, const tA *A4 ) { // -- from http://astronomy.swin.edu.au/~pbourke/
    float qR = (float)(6*A4[1]+(((-1*A4[0]+3*A4[1]-3*A4[2]+A4[3])*AT+3*A4[0]-6*A4[1]+3*A4[2])*AT-2*A4[0]-3*A4[1]+6*A4[2]-A4[3])*AT); return qR/6; }
  template <class tA> inline HSTnDEV float mcHerLikStdSpline( int AFill, float AT, const tA *A4 ) { // -- select HermiteLike-spline
    return AFill&cHerLikFs? mcHerLikB3( AT, A4 ) : mcHerLikPB( AT, A4 ); } // 'bicubic' is sharp; 'Bourke' is blur
    // as debug return (A4[1]+A4[2])/2
#endif // I16_INCV_H_INCLUDED

 