// maintenance-like services
// by S.B.
#include "stdafx.h"
#include "BS_GIni.h"

char xGIni_DLLPath[192], *xGIni_DLLName;

#ifdef _WIN32
  void Gini_StartUps( HMODULE AModule ) { // -- CWD etc at startup
    GetModuleFileName( AModule, xGIni_DLLPath, mcASz(xGIni_DLLPath) );
    xGIni_DLLName = strrchr( xGIni_DLLPath, '\\' )+1; }
#else // ifLinux
  #include <dlfcn.h>
  void Gini_StartUps() { // -- CWD etc at startup
    Dl_info qInfo; // pass the address of the current function or any local symbol
    if (!dladdr( (void*)Gini_StartUps, &qInfo )) {
      ErrMsg( "Failed to find module info" ); return; };
    strncpy( xGIni_DLLPath, qInfo.dli_fname, 191 );
    xGIni_DLLName = strrchr( xGIni_DLLPath, '/' )+1; }
#endif

#ifdef GINIFILE_WINDOWS
  GIniFile::GIniFile() { rFName[0] = rSName[0] = 0; } GIniFile::~GIniFile() {}
  void GIniFile::SetFName( const char *AFileName ) { strcpy( rFName, AFileName ); }
  void GIniFile::SetSName( const char *ASectName ) { if (ASectName) strcpy( rSName, ASectName ); }
  int GIniFile::GetStr( const char *AKeyName, char *AResult ) {
    int qOk = GetPrivateProfileStringA( rSName, AKeyName, "", AResult, 256, rFName );
    if (!qOk) return -1;
    if (AResult[0]=='\\' && AResult[1]<=' ') { // it allows to specify empty lines using back slash
      strcpy( AResult, &(AResult[1]) ); qOk--; }; return qOk; }

  void GIniFile::SetStr( const char *AKeyName, const char *AStr ) {
    if (!AStr || !AStr[0]) AStr = "\\";
    WritePrivateProfileStringA( rSName, AKeyName, AStr, rFName ); }
#else

  GIniFile::GIniFile() { rFName[0] = rSName[0] = 0; iSectionReSet(); }

  GIniFile::~GIniFile() { iSectionWrite(); }

  int GIniFile::GetStr( const char *AKeyName, char *AResult ) {
    for (int qi=0;qi<rN;qi++) if (!strcasecmp( AKeyName, rK[qi] )) {
      strcpy( AResult, rV[qi] ); return (int)strlen(AResult); };
    return -1; }

  void GIniFile::SetStr( const char *AKeyName, const char *AStr ) {
    rChg = true;
    for (int qi=0;qi<rN;qi++) if (!strcasecmp( AKeyName, rK[qi] )) {
      strcpy( rV[qi], AStr ); return; }; if (rN>=cN) return;
    strcpy( rK[rN], AKeyName ); strcpy( rV[rN], AStr ); rN++; }

  void GIniFile::SetFName( const char *AFileName ) { 
    iSectionWrite(); strcpy( rFName, AFileName ); iSectionReSet(); }

  void GIniFile::SetSName( const char *ASectName ) { iSectionWrite();
    if (!ASectName) return; strcpy( rSName, ASectName ); iSectionReSet(); iSectionRead();  }

  void GIniFile::iSectionReSet() {
    rN = 0; rChg = false; for (int qi=0;qi<cN;qi++) rV[qi][0] = rK[qi][0] = 0; }

  static void istruprANSI( char *AStr ) { char qCh = *AStr;
    while (qCh) { if (qCh>='a' && qCh<='z') *AStr = qCh-('a'-'A'); qCh = *++AStr; }; }

  void GIniFile::iSectionRead() {
    char      qS[256], qSName[52]; 
    FILE     *qF = fopen( rFName, "rt" ); if (!qF) return;
    bool      qOn = !rSName[0]; if (!qOn) {
      sprintf( qSName, "[%s]", rSName ); istruprANSI( qSName ); };
    while (fgets( qS, 252, qF )) {
      int     qL = (int)strlen( qS ); if (qL<3) continue;
      while (qS[qL-1]<=' ') if (--qL<3) break; if (qL<3) continue; qS[qL] = 0; // trim ending spaces
      char   *qStr = qS; while (*qStr<=' ') qStr++; qL -= (int)(qStr-qS); // trim trailing spaces
      if (qL<3) continue; 
      if (qOn) { if (qS[0]=='[') break; // from while
        char *qVal = strchr( qStr, '=' ); if (!qVal || (int)(qVal-qStr)<1) continue;
        char *qKey = qVal-1; while (*qKey<=' ') qKey--; qKey[1] = 0; qKey = qS;
        strcpy( rK[rN], qKey ); // then will compare by 'stricmp'
        qVal++; while (*qVal && (*qVal<=' ')) qVal++;
        if ((qVal[0]=='\\' && (!qVal[1]||qVal[1]==';'))
           || qVal[0]==';') qVal[0] = 0; // single_backslash:=empty_string
        strcpy( rV[rN++], qVal ); // from for i
        if (rN>=cN) break; // from while
      } else if (qStr[0]=='[') { istruprANSI( qStr ); // uppercase for compare
        if (!strcmp( qSName, qStr )) qOn = true;
      };
    }; // while
    fclose( qF ); }

  void GIniFile::iSectionWrite() {
    if (!rChg) return; rChg = false;
    char   qSName[52]; sprintf( qSName, "[%s]", rSName ); istruprANSI( qSName );
    const int cBuf = 1<<16; char *qBuf = new char[cBuf], *qBuU = new char[cBuf];
    int    qSz = 0;
    FILE  *qF = fopen( rFName, "rb" ); if (qF) {  
      qSz = (int)fread( qBuf, 1, cBuf, qF ); fclose( qF ); };
    char *qSeB = qBuU, *qSeE = qBuU; if (qSz>0) {
      memcpy( qBuU, qBuf, qSz ); qBuU[qSz] = 0; istruprANSI( qBuU );
      qSeB = strstr( qBuU, qSName ); if (qSeB) {
        qSeE = strstr( qSeB+1, "\n[" ); if (qSeE) qSeE++; }; // && qSeE[-1]=='\r') qSeE--; };
      if (!qSeE) { qSeE = qBuU+qSz; // keep EOF/EOT if already happened
        if (qSeE[-1]==4 || qSeE[-1]==26) qSeE--; };
      if (!qSeB || qSeB>qSeE) qSeB = qSeE; }; // if Sz>0
    char  *qCh = qBuU; // if (qSeB>qBuU) qCh += sprintf( qCh, "\n" );
    qCh += sprintf( qCh, "%s\n", qSName );
    for (int qi=0;qi<rN;qi++) qCh += sprintf( qCh, "%s=%s\n", rK[qi], rV[qi] );
    if (qSeE<qBuU+qSz-2) qCh += sprintf( qCh, "\n" ); // voluntarily add empty line
    qF = fopen( rFName, "wb" ); if (qF) {
      fwrite( qBuf, 1, (int)(qSeB-qBuU), qF );
      fwrite( qBuU, 1, (int)(qCh-qBuU), qF );
      int qPoE = (int)(qSeE-qBuU); fwrite( qBuf+qPoE, 1, qSz-qPoE, qF );
      fclose( qF ); }; // if F<>nil
    delete qBuf; delete qBuU; }
#endif // not GINIFILE_WINDOWS

const char *GIniFile::DGetStr( const char *AKName, const char *ADeflt ) {
  return GetStr(AKName,rResult)>=0?rResult:ADeflt; } // ret ini-string
const char *GIniFile::SGetStr( const char *AKName, const char *ADeflt ) {
if (GetStr(AKName,rResult)>=0) return rResult; SetStr( AKName, ADeflt ); return ADeflt; } // ret/set ini-string

#define mcscanf(S,F,R) sscanf( S, F F F F F F F F F F F F F F F F, R+0, R+1, R+2, R+3, R+4, R+5, R+6, R+7, R+8, R+9, R+10, R+11, R+12, R+13, R+14, R+15 )
int Gini_sscanf( const char *AStr, int *ARes ) { memset( ARes, 0, sizeof(ARes[0])*16 ); return AStr?mcscanf( AStr, "%d", ARes ):-1; }
int Gini_sscanf( const char *AStr, float *ARes ) { memset( ARes, 0, sizeof(ARes[0])*16 ); return AStr?mcscanf( AStr, "%f", ARes ):-1; }

//// Text output on 16-bit frames via binary bitmap ////
void GIni_Lab8x8( int AW, int AH, tBS_Val *ACanvas, const char *AStr, int AV ) {
  // github.com/dhepper/font8x8/blob/master/font8x8_basic.h#L23
  typedef unsigned char tChar8x8[8]; const tChar8x8 font8x8_basic0[] = {
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // 20 space
    {0x18,0x3C,0x3C,0x18,0x18,0x00,0x18,0x00}, // 21 !
    {0x36,0x36,0x00,0x00,0x00,0x00,0x00,0x00}, // 22 "
    {0x36,0x36,0x7F,0x36,0x7F,0x36,0x36,0x00}, // 23 #
    {0x0C,0x3E,0x03,0x1E,0x30,0x1F,0x0C,0x00}, // 24 $
    {0x00,0x63,0x33,0x18,0x0C,0x66,0x63,0x00}, // 25 %
    {0x1C,0x36,0x1C,0x6E,0x3B,0x33,0x6E,0x00}, // 26 &
    {0x06,0x06,0x03,0x00,0x00,0x00,0x00,0x00}, // 27 '
    {0x18,0x0C,0x06,0x06,0x06,0x0C,0x18,0x00}, // 28 (
    {0x06,0x0C,0x18,0x18,0x18,0x0C,0x06,0x00}, // 29 )
    {0x00,0x66,0x3C,0xFF,0x3C,0x66,0x00,0x00}, // 2A *
    {0x00,0x0C,0x0C,0x3F,0x0C,0x0C,0x00,0x00}, // 2B +
    {0x00,0x00,0x00,0x00,0x00,0x0C,0x0C,0x06}, // 2C ,
    {0x00,0x00,0x00,0x3F,0x00,0x00,0x00,0x00}, // 2D -
    {0x00,0x00,0x00,0x00,0x00,0x0C,0x0C,0x00}, // 2E .
    {0x60,0x30,0x18,0x0C,0x06,0x03,0x01,0x00}, // 2F /
    {0x3E,0x63,0x73,0x7B,0x6F,0x67,0x3E,0x00}, // 30 0
    {0x0C,0x0E,0x0C,0x0C,0x0C,0x0C,0x3F,0x00}, // 31 1
    {0x1E,0x33,0x30,0x1C,0x06,0x33,0x3F,0x00}, // 32 2
    {0x1E,0x33,0x30,0x1C,0x30,0x33,0x1E,0x00}, // 33 3
    {0x38,0x3C,0x36,0x33,0x7F,0x30,0x78,0x00}, // 34 4
    {0x3F,0x03,0x1F,0x30,0x30,0x33,0x1E,0x00}, // 35 5
    {0x1C,0x06,0x03,0x1F,0x33,0x33,0x1E,0x00}, // 36 6
    {0x3F,0x33,0x30,0x18,0x0C,0x0C,0x0C,0x00}, // 37 7
    {0x1E,0x33,0x33,0x1E,0x33,0x33,0x1E,0x00}, // 38 8
    {0x1E,0x33,0x33,0x3E,0x30,0x18,0x0E,0x00}, // 39 9
    {0x00,0x0C,0x0C,0x00,0x00,0x0C,0x0C,0x00}, // 3A :
    {0x00,0x0C,0x0C,0x00,0x00,0x0C,0x0C,0x06}, // 3B ;
    {0x18,0x0C,0x06,0x03,0x06,0x0C,0x18,0x00}, // 3C <
    {0x00,0x00,0x3F,0x00,0x00,0x3F,0x00,0x00}, // 3D =
    {0x06,0x0C,0x18,0x30,0x18,0x0C,0x06,0x00}, // 3E >
    {0x1E,0x33,0x30,0x18,0x0C,0x00,0x0C,0x00}, // 3F ?
    {0x3E,0x63,0x7B,0x7B,0x7B,0x03,0x1E,0x00}, // 40 @
    {0x0C,0x1E,0x33,0x33,0x3F,0x33,0x33,0x00}, // 41 A
    {0x3F,0x66,0x66,0x3E,0x66,0x66,0x3F,0x00}, // 42 B
    {0x3C,0x66,0x03,0x03,0x03,0x66,0x3C,0x00}, // 43 C
    {0x1F,0x36,0x66,0x66,0x66,0x36,0x1F,0x00}, // 44 D
    {0x7F,0x46,0x16,0x1E,0x16,0x46,0x7F,0x00}, // 45 E
    {0x7F,0x46,0x16,0x1E,0x16,0x06,0x0F,0x00}, // 46 F
    {0x3C,0x66,0x03,0x03,0x73,0x66,0x7C,0x00}, // 47 G
    {0x33,0x33,0x33,0x3F,0x33,0x33,0x33,0x00}, // 48 H
    {0x1E,0x0C,0x0C,0x0C,0x0C,0x0C,0x1E,0x00}, // 49 I
    {0x78,0x30,0x30,0x30,0x33,0x33,0x1E,0x00}, // 4A J
    {0x67,0x66,0x36,0x1E,0x36,0x66,0x67,0x00}, // 4B K
    {0x0F,0x06,0x06,0x06,0x46,0x66,0x7F,0x00}, // 4C L
    {0x63,0x77,0x7F,0x7F,0x6B,0x63,0x63,0x00}, // 4D M
    {0x63,0x67,0x6F,0x7B,0x73,0x63,0x63,0x00}, // 4E N
    {0x1C,0x36,0x63,0x63,0x63,0x36,0x1C,0x00}, // 4F O
    {0x3F,0x66,0x66,0x3E,0x06,0x06,0x0F,0x00}, // 50 P
    {0x1E,0x33,0x33,0x33,0x3B,0x1E,0x38,0x00}, // 51 Q
    {0x3F,0x66,0x66,0x3E,0x36,0x66,0x67,0x00}, // 52 R
    {0x1E,0x33,0x07,0x0E,0x38,0x33,0x1E,0x00}, // 53 S
    {0x3F,0x2D,0x0C,0x0C,0x0C,0x0C,0x1E,0x00}, // 54 T
    {0x33,0x33,0x33,0x33,0x33,0x33,0x3F,0x00}, // 55 U
    {0x33,0x33,0x33,0x33,0x33,0x1E,0x0C,0x00}, // 56 V
    {0x63,0x63,0x63,0x6B,0x7F,0x77,0x63,0x00}, // 57 W
    {0x63,0x63,0x36,0x1C,0x1C,0x36,0x63,0x00}, // 58 X
    {0x33,0x33,0x33,0x1E,0x0C,0x0C,0x1E,0x00}, // 59 Y
    {0x7F,0x63,0x31,0x18,0x4C,0x66,0x7F,0x00}, // 5A Z
    {0x1E,0x06,0x06,0x06,0x06,0x06,0x1E,0x00}, // 5B [
    {0x03,0x06,0x0C,0x18,0x30,0x60,0x40,0x00}, // 5C backslash
    {0x1E,0x18,0x18,0x18,0x18,0x18,0x1E,0x00}, // 5D ]
    {0x08,0x1C,0x36,0x63,0x00,0x00,0x00,0x00}, // 5E ^
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xFF}, // 5F _
    {0x0C,0x0C,0x18,0x00,0x00,0x00,0x00,0x00}, // 60 `
  };
  const tChar8x8 font8x8_basic1[] = {
    {0x38,0x0C,0x0C,0x07,0x0C,0x0C,0x38,0x00}, // 7B {
    {0x18,0x18,0x18,0x00,0x18,0x18,0x18,0x00}, // 7C |
    {0x07,0x0C,0x0C,0x38,0x0C,0x0C,0x07,0x00}, // 7D }
    {0x6E,0x3B,0x00,0x00,0x00,0x00,0x00,0x00}, // 7E ~
  };
  const int cX0 = -6; int qX = cX0, qY = 2; do {
    unsigned char qCh = *AStr++; if (!qCh) return; 
    if (qCh=='\n') { qY += 10; if (qY>=AH-8) return; qX = cX0; continue; };
    qX += 8; if (qX>=AW-8 || qCh<' ' || qCh>'~') continue; const tChar8x8 *cCh;
    if (qCh>='{') cCh = font8x8_basic1+qCh-'{'; else {
      if (qCh>='a' && qCh<='z') qCh -= 'a'-'A'; cCh = font8x8_basic0+qCh-' '; };
    tBS_Val *qCanvas = ACanvas+qY*AW+qX;
    for (int qi=0;qi<8;qi++) { unsigned char qc = (*cCh)[qi];
      for (int qj=0;qj<8;qj++) if (qc & (1<<qj)) qCanvas[qj]= AV;
      qCanvas += AW; };
  } while (1); }
