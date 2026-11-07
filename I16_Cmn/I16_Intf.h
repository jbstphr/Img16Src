/** *************************************************************** **
 **                        I16 Process library                      **
 **  edition : May-2026                                 by S.B.     **
 ** *************************************************************** **/
#ifndef I16_INTF_H_INCLUDED
#define I16_INTF_H_INCLUDED
  typedef short tI16_Inp; // most typical inp pixels 
  typedef short tI16_Out; // most typical out pixels
  #ifdef INTERNAL_I16_LIBRARY // procedures' declarations ('as is' in the DLL)
    #ifdef _WIN32
      #define I16_DECLSPEC __declspec(dllexport) 
    #else
      #define I16_DECLSPEC
    #endif
  #else // for use in applications
    // for dynamic link include "VAPrcFlu.cpp" into project
    typedef struct {
      #ifdef _WIN32
        HINSTANCE rHDLL; 
      #else
        void *rHDLL; 
      #endif
      int(*rInit)(int,int,int); void(*rRelease)();
      void(*rDbgDmpReset)(); void(*rDbgDmp)(const char*); 
      int(*rInvoke)(int,int,const tI16_Inp*,int,tI16_Out*,int);
      void(*rParSet)(int,int); int(*rParGet)(int);
      char*(*rParStr)(int); } tI16_Intf;
    bool I16_MapDLL( const char *ADllPath, tI16_Intf *AIntf );
    void I16_ReleaseDLL( tI16_Intf *AIntf );
    #ifdef _WIN32
      #define I16_DECLSPEC __declspec(dllimport)
      // #pragma comment( lib, "VAPrcFlu.lib" )
    #else
      #define I16_DECLSPEC
    #endif
  #endif

  // static link
  // refer paragraphs about ~"presentation enhacmenent API" in a project's documentation
  I16_DECLSPEC int I16_Init( int AW, int AH, int AModeId ); // inits and allocates resources
    // <0: an error (e.g. wrong dim of existing CalData)
    // =0: further _Invoke will have nothing to do (e.g. empty CalDir)
    //     the caller would better just call _Release now
    // >0: Ok (e.g. #of occupied frames)
    // >1 shl 18: maximum dimension "(Hmx shl 16) or Wmx"
  I16_DECLSPEC void I16_Release(); // releases resources
  I16_DECLSPEC void I16_DbgDmpReset(); // resets debug-dump
  I16_DECLSPEC void I16_DbgDmp( const char *ASubDir ); // nil:minimalistic trace, "":verbose trace, "xyz":dump images
  I16_DECLSPEC int I16_Invoke( int AW, int AH, const tI16_Inp *ASrc, int AWSrc, tI16_Out *ADst, int AAdj ); // invokes process
    // <0: an error (bit "1 shl 31")
    // =0: the caller has nothing to do with this Dst
    // >0: processed (~"1=Ok", "flags", etc )
    // >=1 shl 18: changed dimension "(Hout shl 16) or Wout"
    // Src=nil,WSrc=0: invoke for single summary
  I16_DECLSPEC void I16_ParSet( int AParIdx, int AParVal ); // auxilliary int parameters' assign
  I16_DECLSPEC int I16_ParGet( int AParIdx ); // auxilliary parameters' read
  I16_DECLSPEC char *I16_ParStr( int AParIdx ); // 1D-array for auxiliary ASCIIz-pars "char[256]"
    // This API might be absent in a specfic library. When it is exported then:
    // string #0 for Init: specifies ~CalDir or preset's ini-profile,
    //         -"- Invoke: optionally specifies ~current filename (as meta-information).
#endif // ifdef I16_INTF_H_INCLUDED
