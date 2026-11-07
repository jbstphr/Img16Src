#ifndef BS_GINIFILE_H_INCLUDED
#define BS_GINIFILE_H_INCLUDED

extern char xGIni_DLLPath[], *xGIni_DLLName;

#ifdef _WIN32
  void Gini_StartUps( HMODULE AModule ); // CWD etc at startup
#else // ifLinux
  void Gini_StartUps(); // CWD etc at startup
#endif

//// A wrapper around INI-files' (or registry's) API ////
// "Key=":denotes absent ini-string, "Key=\":explicit backslash denotes empty ini-string.
// #define GINIFILE_WINDOWS // use Windows API
class GIniFile { public:
  char rFName[192], rSName[48]; GIniFile(); ~GIniFile();
  void SetFName( const char *AFileName ); // sets name of ini file
  void SetSName( const char *ASectName ); // sets name of section
  int  GetStr( const char *AKeyName, char *AResult ); // returns ini-value and its length
  void SetStr( const char *AKeyName, const char *AStr ); // writes ini-value
  #ifndef GINIFILE_WINDOWS // hand-made INI-file
    protected: // internal procedures
    static const int cN = 48; // max #entries in a section
    int   rN; bool rChg; // #entries in cur section; fact of SetIniStr
    char  rV[cN][256], // values in uppercase
          rK[cN][48]; // keys in uppercase
    void iSectionReSet(); 
    void iSectionRead();
    void iSectionWrite();
  #endif
  char rResult[256];
  public:
  const char *DGetStr( const char *AKName, const char *ADeflt ); // ret existing ini-string or default one
  const char *SGetStr( const char *AKName, const char *ADeflt ); // ret ini-string (set if doesn't exist)
};

int Gini_sscanf( const char *AStr, int *ARes ); // ret #ints in the string
int Gini_sscanf( const char *AStr, float *ARes ); // ret #floats in the string

//// Text output on 16-bit frames via binary bitmap ////
void GIni_Lab8x8( int AW, int AH, tBS_Val *ACanvas, const char *AStr, int AV );

#endif // #ifndef BS_GINIFILE_H_INCLUDED
