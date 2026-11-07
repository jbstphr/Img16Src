// I16Prc.cpp : Defines the entry point for the console application.
// Refer ReadMe.txt, VaPrcFlu.rtf.
#include "stdafx.h"
#include "../I16_Cmn/I16_Intf.h"

char *iFNam( char *AFPath ) { // file's name at the end of a path
  char *qF0 = strrchr( AFPath, '\\' ), *qF1 = strrchr( AFPath, '/' );
  if (qF0<qF1) qF0 = qF1;return qF0?qF0+1:AFPath; }

#ifdef _WIN32 // Windows
  #include <IO.h>
  struct dirent { char d_name[256]; }; // to use PathMatchSpecA
  struct DIR { intptr_t rFFH; struct _finddata_t rFindData; struct dirent rRes; };
  DIR *opendir( const char *AName ) {
    DIR *qDir = new DIR; memset( qDir, 0, sizeof(DIR) );
    char qS[256]; sprintf( qS, "%s\\*.*", AName );
    qDir->rFFH = _findfirst( qS, &qDir->rFindData );
    if (qDir->rFFH==-1) { delete qDir; return NULL; };
    strcpy( qDir->rRes.d_name, qDir->rFindData.name );
    qDir->rFindData.name[0] = 0; // indicate first entry
    return qDir; }
  void closedir( DIR *ADir ) {
    if (ADir) { _findclose( ADir->rFFH ); delete ADir; }; }
  struct dirent *readdir( DIR *ADir ) {
    if (!ADir->rFindData.name[0]) { ADir->rFindData.name[0] = 1; return &ADir->rRes; };
    if (_findnext( ADir->rFFH, &ADir->rFindData )) return NULL;
    strcpy( ADir->rRes.d_name, ADir->rFindData.name ); return &ADir->rRes; }
  //#define strncasecmp strnicmp // compare two strings of limited length ignoring case
  //#define strcasecmp stricmp   // compare two strings ignoring case
  void iGetModuleFileName( char *APath ) { APath[0] = 0; GetModuleFileNameA( NULL, APath, 192 ); }
  char *iGetDLLFileNameExtPtr( char *AExePath ) {
    char *qdlleptr = strrchr( AExePath, '.' );  if (!qdlleptr) qdlleptr = AExePath+strlen( AExePath );
    strcpy( qdlleptr, ".dll" ); return qdlleptr; }
  int iGetTickCount() { return (int)GetTickCount(); }
  #define cSlash "\\"
#else // Linux
  #include <dirent.h>
  void iGetModuleFileName( char *APath ) { APath[readlink("/proc/self/exe",APath,192)] = 0; }
  char *iGetDLLFileNameExtPtr( char *AExePath ) {
    char *qdlleptr = strrchr( AExePath, '.' ); if (!qdlleptr) qdlleptr = AExePath+strlen( AExePath );
    strcpy( qdlleptr, ".so" ); return qdlleptr; }
  int iGetTickCount() { struct timespec qNow; clock_gettime( CLOCK_BOOTTIME, &qNow );
    return (int)(qNow.tv_sec*1000L+qNow.tv_nsec/1000000L); }
  #define cSlash "/"
#endif

void iWriteBMP( const char *AFName, int AW, int AH, const tI16_Out *ASrc, int AZer, int AGam ) {
  typedef unsigned int  tES_DWORD; typedef unsigned int tES_LONG;
  typedef unsigned short tES_WORD; typedef unsigned char tES_BYTE;
  #pragma pack(2)
  typedef struct { tES_DWORD biSize; tES_LONG biWidth; tES_LONG biHeight; tES_WORD biPlanes;
    tES_WORD biBitCount; tES_DWORD biCompression; tES_DWORD biSizeImage; tES_LONG biXPelsPerMeter;
    tES_LONG biYPelsPerMeter; tES_DWORD biClrUsed; tES_DWORD biClrImportant; } tES_BitMapInfoHeader;
  typedef struct { tES_BYTE rgbBlue; tES_BYTE rgbGreen; tES_BYTE rgbRed; tES_BYTE rgbReserved; } tES_RgbQuad;
  typedef struct { tES_BitMapInfoHeader rBMI; tES_RgbQuad rPal[256]; } tES_BMP;
  typedef struct { tES_WORD bfType; tES_DWORD bfSize; tES_WORD bfReserved1; tES_WORD bfReserved2; tES_DWORD bfOffBits; } tES_BMPFileHdr;
  tES_BMP   qBMP;  memset( &qBMP, 0, sizeof(tES_BMP) );
  int       qStride = AW; if (qStride&3) qStride = (qStride+4)&~3;
  qBMP.rBMI.biSize = sizeof(tES_BitMapInfoHeader); qBMP.rBMI.biWidth = AW; qBMP.rBMI.biHeight = AH; 
  qBMP.rBMI.biPlanes = 1; qBMP.rBMI.biBitCount = 8; qBMP.rBMI.biSizeImage = qStride*AH;
  qBMP.rBMI.biClrUsed = qBMP.rBMI.biClrImportant = 256;
  qBMP.rBMI.biXPelsPerMeter = qBMP.rBMI.biYPelsPerMeter = (long)(1000/0.20f);
  for (int qi=0;qi<256;qi++) qBMP.rPal[qi].rgbRed = qBMP.rPal[qi].rgbGreen = qBMP.rPal[qi].rgbBlue = (tES_BYTE)qi;
  tES_BYTE *qByteBuf = new tES_BYTE[qBMP.rBMI.biSizeImage], *qDst = qByteBuf+qBMP.rBMI.biSizeImage-qStride; // bitmap's-data
  for (int qy=0;qy<AH;qy++) {
    for (int qx=0;qx<AW;qx++) { int qV = (ASrc[qx]-AZer)*256/AGam;
      if (qV<=0) qDst[qx] = 0; else if (qV>=255) qDst[qx] = 255; else qDst[qx] = (tES_BYTE)qV;
    }; // for x
    ASrc += AW; qDst -= qStride;
  }; // for y
  // write file
  char      qS[256]; strcpy( qS, AFName ); strcpy( strchr(qS,'.'), ".bmp"  );
  tES_BMPFileHdr qFHdr; memset( &qFHdr, 0, sizeof(tES_BMPFileHdr) );
  qFHdr.bfType = 0x4D42; qFHdr.bfOffBits = sizeof(tES_BMPFileHdr)+sizeof(tES_BMP);
  qFHdr.bfSize = qFHdr.bfOffBits+qBMP.rBMI.biSizeImage;
  FILE     *qF = fopen( qS, "wb" ); if (qF) {
    fwrite( &qFHdr, 1, sizeof(tES_BMPFileHdr), qF );
    fwrite( &qBMP,1, sizeof(tES_BMP), qF );
    fwrite( qByteBuf, 1, qBMP.rBMI.biSizeImage, qF ); fclose( qF ); };
  delete [] qByteBuf; }

void iDisplayVPFRes( int ARes, const char *AMsg ) { // print out result-code
  // "application-meaning" of specific output of VAPrcFlu.dll
  //   <0: an error,
  //   =0: VaCal didn't change input data,
  //   >0: adhoc flags (just to trace),
  //   >1 shl 18: out dimension <changed against the inp>.
  if (ARes<0) printf( "ERR " ); printf( "%s", AMsg ); printf( ":=" );
  if (ARes<=0) printf( "%d", ARes );
  else if (ARes>1<<19) printf( "%dx%d", ARes&0x7FFF, ARes>>16 );
  else printf( "%X", ARes ); printf( "\n" ); fflush( stdout ); }

int xW, xH; // frames' dimension from the command-line
unsigned iInterpretVPFRes( int ANonFatMn, int ARes ) { // retcode by VaCal.Init,Invoke 
  // "behavioral-meaning" of a result-code for this files' iterator
  //   0: err fatal,
  //   1: err non-fatal (ignore this specific output only),
  //   2: the caller to proceed as is,
  //   'W or (H shl 16)' otherwise.
  if (ARes<ANonFatMn) return 0; if (ARes<0) return 1; if (!ARes) return 2;
  if (ARes>1<<19) return ARes; return xW|(xH<<16); }

static bool iis_name_to_process( const char *AFName ) {
  const char *cFExt = strrchr( AFName+1, '.' ); if (!cFExt) return false;
  const int c_poss[2] = {3,5}; for (int qi=0;qi<2;qi++) { int c_p = c_poss[qi];
    if (cFExt-AFName>c_p && cFExt[-c_p]=='_') { bool qSuit = false;
    do { c_p--; char qc = cFExt[-c_p];
      if ((qc<'0'||qc>'9') && (qc<'A'||qc>'F')) { qSuit = true; break; };
    } while (c_p>1);
    if (!qSuit) return false; }; }; return true; } // for i; if FExt[-p]='_'

int main( int argc, char* argv[] ) {
  static const int
    cFl_SavBMP = 1, // write *.BMP
    cFl_OuNRAW = 2, // don't write file *.RAW
    cFl_OuNadj = 4, // don't write Adj in RAW's name
    cFl_Debugs = 15<<4;
  char qS[256]; // ad-hoc string-buffer
  #ifdef _WIN32  // for debug command-line when called by ~hand-made CreateProcess(...)
    char *qChDbg = qS+sprintf( qS, "I16.exe" ); for (int qa=0;qa<__min(argc,7);qa++) {
      const char *qarg = argv[qa]; int qL = (int)strlen( qarg ); if (qL>22) qarg += qL-15;
      qChDbg += sprintf( qChDbg, " %s", qarg ); }; OutputDebugString( qS );
  #endif
  if (argc<6) { iGetModuleFileName( qS ); printf(
    "Usage\t:\\>%s WDT HGT WildCard.raw ModeId Adj <Flg>\n", iFNam(qS) ); return -99; }; // if argc<
  xW = atoi( argv[1] ); xH = atoi( argv[2] ); // dimension
  // This console-file-iterator invokes DbgDmp(<>""), which (in some VAPrcFlu.dll) changes
  // current directory. Hence (for relative FPath) the wildcard-buffer here saves initial
  // CWD (where the Names[] are enlisted).
  const char *cFPath = argv[3]; char qFPath[192]; // wildcard-buffer
  if (cFPath[0]=='\\' || cFPath[0]=='/' || cFPath[1]==':') strcpy( qFPath, cFPath ); else
    sprintf( qFPath+strlen(getcwd(qFPath,128)), cSlash"%s", cFPath );
  const int
    cModeId = atoi(argv[4]), // refer about 'ModeId' in .\I16_Test\ReadMe.txt
    cAdj = strtol(argv[5],NULL,16), // refer about 'Adjustments' in .\I16_Test\ReadMe.txt
    cFlg = argc>6?strtol( argv[6], NULL, 16 ) : 0, // ad-hoc flags for local use
    cNames = 512, cName = 32;
  char *qFNam = iFNam( qFPath ); if (qFNam>qFPath) qFNam[-1] = 0; // assign file-name's pointer 
  char  qNames[cNames][cName]; int qCnt = 0; // <sorted> files to process
  DIR  *qDP = opendir( qFNam>qFPath?qFPath:"./" ); if (qDP) { // POSIX file-iterator
    struct dirent *qdirp = readdir( qDP ); while (qdirp) {
      if (qdirp->d_name[0]=='.') goto labContinue;
      #ifdef _WIN32
        if (!PathMatchSpecA( qdirp->d_name, qFNam )) goto labContinue;
      #else
        if (fnmatch( qFNam, qdirp->d_name, 0 )) goto labContinue;
      #endif  
      { int qLen = (int)strlen( qdirp->d_name ); if (qLen>=cName) goto labContinue;
        if (!iis_name_to_process(qdirp->d_name)) goto labContinue;
        { int qIns = 0; while (qIns<qCnt) { if (strcmp(qNames[qIns],qdirp->d_name)>0) break; qIns++; };
          memmove( qNames+qIns+1, qNames+qIns, sizeof(qNames[0])*(qCnt-qIns) ); // ascending sort
          memcpy( qNames[qIns], qdirp->d_name, (qLen+1)*sizeof(char) ); }; // of Ins:=
      if (++qCnt>=cNames) break; }; // of Len:=
      labContinue: qdirp = readdir( qDP ); // findnext
    }; closedir( qDP ); // while
  }; // if DP
  if (!qCnt) { printf( "No files by wildcard %s\n", cFPath ); return -97; };
  tI16_Intf qIntf; { // I16_MapDLL calls loadlibrary in ANSI-coding
    // try *.dll co-existing with *.exe, then try upper directory
    char *qdllpath = qS; iGetModuleFileName( qdllpath );
    char *qdlleptr = iGetDLLFileNameExtPtr( qdllpath );
    //char *qdlleptr = strrchr( qdllpath, '.' );
    //if (!qdlleptr) qdlleptr = qdllpath+strlen( qdllpath ); strcpy( qdlleptr, ".dll" );
    struct stat qStt; if (stat(qdllpath,&qStt)==-1) { char *qdllnam = NULL;
      while (--qdlleptr>qdllpath) if (*qdlleptr=='\\' || *qdlleptr=='/') {
        if (!qdllnam) { qdllnam = qdlleptr; } else { strcpy( qdlleptr, qdllnam ); break; }; };
    };
    if (!I16_MapDLL( qdllpath, &qIntf )) return -96; // ^^ map the DLL into process' memory
  };
  // process
  char     *qParStr0 = qIntf.rParStr( 0 ); // ^^ aux_str0 := files_dir without trailing slash
  if (qParStr0) strcpy( qParStr0, qFPath );
  int       qDbgDmp = ((cFlg&cFl_Debugs)>>4)-1; if (qDbgDmp>8) qDbgDmp = 8;
  qIntf.rDbgDmp( "" ); // always allow debug-trace in I16_Init
  int       qRes = qIntf.rInit( xW, xH, cModeId ); // ^^ init DLL's resources
  sprintf( qS, "I16_Init(%dx%d,%d)", xW, xH, cModeId ); iDisplayVPFRes( qRes, qS );
  unsigned  qRet = iInterpretVPFRes( 0, qRes ); if (qRet<1<<19) return -95;
  // allocate fluoro-process buffers
  unsigned  qNinp = xW*xH; // input buffer's size
  tI16_Inp *qInp = new tI16_Inp[xW*xH]; memset( qInp, 0, sizeof(tI16_Inp)*xW*xH );
  qRet = (qRet&0x7FFF)*(qRet>>16); // allocate not less than input buffer
  tI16_Out *qOut = new tI16_Out[__max(qNinp,qRet)]; memset( qOut, 0, sizeof(tI16_Out)*qRet );
  int       qTicks = 0, qC = 0, qAppRetCode = 0;
  if (qDbgDmp<0) qIntf.rDbgDmp( NULL ); else
  if (qDbgDmp>0) qIntf.rDbgDmpReset();
  if (qFNam>qFPath) qFNam[-1] = '/';
  // invocation's loop
  for (int qi=0;qi<qCnt;qi++) {
    if (qi<qDbgDmp) { // ^^ enable debug-dump
      char qDbgPfx[40]; sprintf( qDbgPfx, "%c%02d_%04X", 'a'+qi, cModeId, cAdj );
      qIntf.rDbgDmp( qDbgPfx );
    } else if (qi==qDbgDmp && qi>0) 
      qIntf.rDbgDmp( "" );
    strcpy( qFNam, qNames[qi] ); // read input file
    FILE    *qF = fopen( qFPath, "rb" );
    if (fread( qInp, sizeof(tI16_Inp), qNinp, qF )!=qNinp)
      printf( "WRN\tfread(%s)!=%d", qFNam, qNinp );
    fclose( qF );
    if (qParStr0) strcpy( qParStr0, qFNam ); // ^^ aux_str0 := file_name
    int     qTick0 = iGetTickCount();
    qRes = qIntf.rInvoke( xW, xH, qInp, xW, qOut, cAdj ); // ^^ invoke process
    qRet = iInterpretVPFRes( -63, qRes ); // _Invoke:=-63..-1 denotes "err non-fatal"
    if (qRet<2 || qi<10 || !(qi%10)) {
      sprintf( qS, "\t%s(%X)", qFNam, cAdj ); iDisplayVPFRes( qRes, qS ); };
    switch (qRet) {
      case 0: qAppRetCode = -89; break; // err fatal
      case 1: continue; }; // err non-fatal
    qTicks += abs(iGetTickCount()-qTick0); qC++; // encount timing
    if (qRet>2)  { // _Invoke returned an image
      int   qWo = qRet&0x7FFF, qHo = qRet>>16;
      if (!(cFlg&cFl_OuNRAW)) { // output *_??_????.raw
        char *qOutSfx = strrchr( qFNam, '.' ); qOutSfx += sprintf( qOutSfx, "_%02d", cModeId );
        if (!(cFlg&cFl_OuNadj)) qOutSfx += sprintf( qOutSfx,"_%04X", cAdj ); qOutSfx += sprintf( qOutSfx, ".raw" );
        qF = fopen( qFPath, "wb" );
        if (qF) { fwrite( qOut, sizeof(tI16_Out), qWo*qHo, qF ); fclose( qF ); }
           else printf( "ERR\tcreate(%s):=%d\n", qFPath, errno );
      }; // output *_??_????.raw
      // Adj=0 indicates "free-style summary"; a VPF.dll returns 0 if it has nothing to do
      if (cFlg&cFl_SavBMP) iWriteBMP( qFPath, qWo, qHo, qOut, 0, 1<<14 ); // output *.bmp
    }; // _Invoke returned an image
  }; // for i
  if (qDbgDmp<0) qIntf.rDbgDmp( "" ); // always allow debug-trace at finalization
  if (qC && !qAppRetCode) { // finalization 
    printf( "%d I16_Invoke(%04dx%04d,%X) took %dms,%.1ffps\n", // print the timing
       qC, xW, xH, cAdj, qTicks, qTicks>0?1000.0*qC/qTicks:0.0f ); fflush( stdout );
    if (qParStr0) *qParStr0 = 0; // ^^ call Invoke(Src=nil) for "single summary"
    qRes = qIntf.rInvoke( xW, xH, NULL, xW, qOut, 0 );
    iDisplayVPFRes( qRes, "I16_Sum" );
    qRet = iInterpretVPFRes( -63, qRes ); // _Invoke:=-63..-1 denotes "err non-fatal"
    if (!qRet) qAppRetCode = -88; else
    if (qRet>2) { // _Invoke returned an image
      int   qWo = qRet&0x7FFF, qHo = qRet>>16;
      sprintf( qFNam, "Z%04dx%04d.raw", qWo, qHo );
      FILE *qF = fopen( qFPath, "wb" );
      if (qF) { fwrite( qOut, sizeof(tI16_Out), qWo*qHo, qF ); fclose( qF ); }
         else printf( "ERR create(%s):=%d\n", qFPath, errno );
    }; // _Invoke returned an image
  }; // finalization
  delete []qInp; delete []qOut; fflush( stdout );// deallocate fluoro-process buffers
  qIntf.rRelease(); I16_ReleaseDLL( &qIntf ); // ^^ release resources and memory
  return qAppRetCode; }