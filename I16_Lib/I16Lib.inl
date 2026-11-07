// -- Include into I16Lib.cpp.
// -- Frequently changing declarations for jUnsharp..., I16_Init.

//// depend on tUnSharp in BS_Ut5.h ////
static void jUnSharpLoad( GIniFile *AF, tUnSharp *AUS, tUnSharpAdjs *AUSA,
  const tAdjMnMx *AAdjMnMx ) { int qI[16]; char qS[256]; // load 'tUnSharp', init 'tUnSharpAdjs'
  if (Gini_sscanf( AF->DGetStr("nEnabled",NULL), qI )>0) { AUS->rEnabled = qI[0]>0; } else {
    sprintf( qS, "%d", (int)AUS->rEnabled ); AF->SetStr( "nEnabled", qS ); };
  if (!AUS->rEnabled) return;
  // fields of tUnSharp
  memset( AUSA, 0, sizeof(tUnSharpAdjs) ); // mainly USA.N := 0
  jjUSkey( AF, 1, AUS->rDeR, "nDeG1", "Thr" );
  jjUSkey( AF, 4, AUS->fDeG, "fDeG4", "Gamma Kgr KnP KSnP" );
  jjUSkey( AF, 1, AUS->rDeR, "nDeR1", "RadR" );
  jjUSkey( AF, 3, AUS->fDeR, "fDeR3", "KR KU Kmf" );
  jjUSkey( AF, 4, AUS->fDeD, "fDeD4", "A+a* A+l* L+a* L+l*" );
  jjUSkey( AF, 4, AUS->rDeN, "nDeN4", "RadN LgN RadN2 ZN2" );
  jjUSkey( AF, 2, AUS->fDeN, "fDeN2", "KN KN2" );
  jjUSkey( AF, 3, AUS->rDeS, "nDeS3", "RadSnP r0SnP r1SnP" );
  jjUSkey( AF, 3, AUS->rMgn, "nMgn3", "Msk Bev Cor" );
  jjUSkey( AF, 6, AUS->rA,   "nA6",   "RadHF RadMM TGMM ZMM Lim LimP" );
  jjUSkey( AF, 5, AUS->fA,   "fA5",   "KHF KMM KAO KCPU DnSpDi" );
  jjUSkey( AF, 3, AUS->fB,   "fB3",   "GammaSpa CLA Att" );
  AUS->fB[0] = iiGammaByInp( AUS->fB[0] ); AUS->fB[1] = iiGammaByInp( AUS->fB[1] );
  jjUSkey( AF, 3, AUS->rC,   "nC3",   "RadBas RadBig RadMM" );
  jjUSkey( AF, 3, AUS->rC,   "fC3",   "RefBig DelBig KMM" );
  jjUSkey( AF, 2, AUS->rD,   "nD2",   "Blk Blr" );
  jjUSkey( AF, 3, AUS->fD,   "fD3",   "BodyOne ExtMn ExtMx" );
  jjUSkey( AF, 4, AUS->rE,   "nE4",   "Blk Blr Lim Zero" );
  jjUSkey( AF, 4, AUS->fF,   "fF4",   "Att EmLthn EmLedg MxLbias" );
  const int ci1 = mcASz(AUS->fFK); for (int qi1=0;qi1<ci1;qi1++) { sprintf( qS, "fFK%c3", qi1+'a' );
    const char *cShint[ci1] = {"GSh +|*| D Dcla Lap","GDi","Edg"}; jjUSkey( AF, 3, AUS->fFK[qi1], qS, cShint[qi1] ); };
  const int ci2 = mcASz(AUS->rFC); for (int qi2=0;qi2<ci2;qi2++) { sprintf( qS, "nFC%c5", qi2+'a' ); 
    const char *cShint[ci2] = {"Thn","Thk","GDi","Edg"}; jjUSkey( AF, 5, AUS->rFC[qi2], qS, cShint[qi2] ); };
  jjUSkey( AF, 2, AUS->rH,   "rH2",   "Lim LimP" );
  jjUSkey( AF, 4, AUS->rJ,   "nJ4",   "RadMM RadMMO RadENm Lim" );
  jjUSkey( AF, 2, AUS->fJ,   "fJ2",   "KMM KENm" );
  jjUSkey( AF, 3, AUS->rL,   "nL3",   "RadHess RadMM Lim" );
  jjUSkey( AF, 3, AUS->fL,   "fL3",   "KEdg KRdg KMM" );
  jjUSkey( AF, 3, AUS->fN,   "fN3",   "Thk GDi Edg" );
  jjUSkey( AF, 4, AUS->fSA,  "fSA4",  "Body +* GSpa GCLA HiF LapMid" );
  jjUSkey( AF, 3, AUS->fSB,  "fSB3",  "Lap +* CLAInc LapAda LapAdaEff" ); const int cPrecs_fU[] = {2,4,4};
  jjUSkey( AF, 3, AUS->fU,   "fU3",   "BlOuDG Air0 Air1", cPrecs_fU );
  //  float qn[3] = {0.4f,0.75f,0.7f}, // Tag=1
  //        qx[3] = {2.6f,1.22f,1.6f}; // Tag=15
  float qn[3]; memcpy( qn, AAdjMnMx->rn, sizeof(qn) ); // jjUSkey( AF, 3, qn, "fBCSn", "B* C* S*" );
  float qx[3]; memcpy( qx, AAdjMnMx->rx, sizeof(qx) ); // jjUSkey( AF, 3, qx, "fBCSx", "B* C* S*" ); 
  // BSCindex #0:"Brightness", #1:"Contrast", #2:"Sharpness"
  jjadj( AUSA, AUS, 0, true, qx[0], qn[0], &AUS->fB[0] );  // Gamma0
  jjadj( AUSA, AUS, 0, true, qx[0], qn[0], &AUS->fB[1] );  // Gamma1
  jjadj( AUSA, AUS, 1, true, qn[1], qx[1], &AUS->fDeR[2] ); // des_Kmf
  jjadj( AUSA, AUS, 1, true, qn[1], qx[1], &AUS->fA[0] );   // a_KHF
  jjadj( AUSA, AUS, 1, false,qx[1], qn[1], &AUS->rD[1] );   // d_Blr
  jjadj( AUSA, AUS, 1, true, qn[1], qx[1], &AUS->fSB[0] );  // +Zb0*CLAInc
  jjadj( AUSA, AUS, 2, true, qn[2], qx[2], &AUS->fSA[2] ); // +Za2*HiFrem
  jjadj( AUSA, AUS, 2, true, qn[2], qx[2], &AUS->fSA[3] ); // +Za3*LapMid
  jjadj( AUSA, AUS, 2, true, qn[2], qx[2], &AUS->fSB[1] ); // +Zb1*LapAda
  jjadj( AUSA, AUS, 2, true, qn[2], qx[2], &AUS->fSB[2] ); // +Zb2*LapAdaEff
}

static const tUnSharp cUS01 = { true, // pano x2717x1372/111.raw ~#1
  {28}, // des_Thr
  {0, 0, 1.5f, 0.6f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {8}, // des_RadR
  {0.7f, 1.0f, 2.2f}, // des_KR des_KU des_Kmf
  {0, 0.1f, 3.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {9, 108, 7, 22}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.1f, 0.6f}, // des_KN des_KN2
  {104, -333, 111}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, 0, 0, 48, -1}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {3.0f, 5.0f, 0, 1.0f, 0.0f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {2.0f, 3.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 3, 6}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {5, 22}, // d_Blk d_Blr
  {0.0f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {6, 3, 99, 0}, // e_Blk e_Blr e_Lim e_Zero
  {2.0f, 0.33f, 0.33f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 1.0f, 0.0f},  // f_K#a ~GSh
   { 0.0f, 1.0f, 1.0f},  //    #b ~GDi <opt>
   { 0.0f,-0.5f, 1.0f}}, //    #c ~Thn <opt>
  {{   0,2222,  0, 99,  1},  // f_C#a mskThn
   {-555,1111,123,  0,  3},  //    #b mskThk
   {   0, 1777, 99, 55,103},  //   #c mskGDi
   {  99, 987,  0,  0,  1}}, //    #d mskEdg
  {111, 333}, // h_Lim h_LimP
  {8, 0, 0, 0}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.0f, 0.0f}, // j_KMM j_ENm
  {2, 5, 333}, // l_RadHess l_RadMM l_Lim
  {2.8f, 0.03f, 1.0f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.25f, 0.25f, 12.0f, 2.2f}, // s_Za
  {0.66f, 0.66f, 4.0f}, // s_Zb
  {0.2f, 3E-3f, 1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static tUnSharp cUS03 = { true, // cephalo-SAA x2289x2592/01.raw ~#3
  {2}, // des_Thr
  {0, 0, 0, 0.8f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {5}, // des_RadR
  {0.6f, -1.0f, 0.5f}, // des_KR des_KU des_Kmf
  {-0.7f, 0.0f, 1.0f, -0.8f}, // A+a* A+l* L+a* L+l*
  {0, 0, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {0.0f, 0.0f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {3, 9, -1, 0, 2222, 0}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.0f, 0.5f, 0.0f, 1.0f, 0.0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {7.0f, 5.0f, 8.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 5, 9}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {6, 22}, // d_Blk d_Blr
  {0.25f, 0.66f, 0.77f}, // d_BodyOne d_ExtMn d_ExtMx
  {8, 1, 156, 111}, // e_Blk e_Blr e_Lim e_Zero
  {0.3f, 0.04f, 0.015f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 0.0f, 0.0f},  // f_K#a ~GSh
   { 0.0f, 0.0f, 0.0f},  //    #b ~GDi <opt>
   { 0.0f, 0.0f, 0.0f}}, //    #c ~Thn <opt>
  {{  33, 333,  0,  0,  2},  // f_C#a mskThn
   {   0, 555,  0, 99,104},  //    #b mskThk
   {   0, 444,  0,  0,102},  //    #c mskGDi
   {   0,  55,  0,  0,  3}}, //    #d mskEdg
  {-33, 66}, // h_Lim h_LimP
  {8, 5, 0, 33}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.5f, 0.0f}, // j_KMM j_ENm
  {2, 4, 66}, // l_RadHess l_RadMM l_Lim
  {1.0f, 0.5f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {1.0f, 0.0f, 33.0f, -1.0f}, // s_Za
  {1.1f, 1.0f, 17.0f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS05 = { true, // low Poisson x1610x2220/e21.raw
  {12}, // des_Thr
  {1.5f, 0, 2.0f, 0.48f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {7}, // des_RadR
  {0.5f, -1, 3.3f}, // des_KR des_KU des_Kmf
  {0.33f, 0.0f, 2.0f, -1.1f}, // A+a* A+l* L+a* L+l*
  {10, 7, 7, 99}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.2f, 1.6f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {-25, 0, 2}, // mgn_Msk mgn_Bev mgn_Cor "144 for 2.0 810x1100"
  {1, 3, 8, 22, 66, 33}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.0f, 0.0f, 33.3f, 1.25f, 0.55f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {1.0f, 2.0f, 0.5f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 3, 4}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.3f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {6, 33}, // d_Blk d_Blr
  {0.0f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {7, 7, 128, 0}, // e_Blk e_Blr e_Lim e_Zero
  {3.0f, 0.07f, 0.07f, 1.01f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 0.0f, 1.0f, 0.0f},  // f_K#a ~GSh
   { 1.0f, 2.0f, 1.0f},  //    #b ~GDi <opt>
   { 1.0f, 1.0f,-1.5f}}, //    #c ~Thn <opt>
  {{-111,1111,  0,123,  2},  // f_C#a mskThn
   {   0, 888, 99,  0,104},  //    #b mskThk
   {-555,3333,222,  0,102},  //    #c mskGDi
   {   0, 400, 55,  0,  1}}, //    #d mskEdg
  {333, 0}, // h_Lim h_LimP
  {8, 8, 1, 0}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.7f, 0.25f}, // j_KMM j_ENm
  {2, 4, 99}, // l_RadHess l_RadMM l_Lim
  {0.5f, 0.15f, 0.7f}, // l_KEdg l_KRdg l_KMM
  {0.15f, 0.15f, 0.15f}, // n_Z
  {0.5f, 0.5f, 30.0f, 1.0f}, // s_Za
  {0.44f, 0.3f, 2.6f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS06 = { true, // molar XSen x1610x2220
  {11}, // des_Thr
  {1.5f, 0, -0.7f, 0.2f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {11}, // des_RadR
  {0.27f, 0.05f, 1.0f}, // des_KR des_KU des_Kmf
  {-1.0f, 0.0f, 2.0f, -1.2f}, // A+a* A+l* L+a* L+l*
  {9, 107, 10, 44}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.5f, 1.0f}, // des_KN des_KN2
  {-206, 1444, 99}, // des_RadSnP des_r0SnP des_r1SnP
  {-25, 0, 2}, // mgn_Msk mgn_Bev mgn_Cor
  {2, 5, 4, 111, 33, 0}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {0.7f, 3.3f, 111.1f, 1.4f, 0.55f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi=0.75f
  {1.0f, 1.8f, 0.7f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 7, 0}, // c_RadBas c_RadBig c_RadMM
  {88.0f, 0.1f, 0.0f}, // c_RefBig c_DelBig c_KMM
  {7, 15}, // d_Blk d_Blr
  {0.0f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {7, 7, 128, 0}, // e_Blk e_Blr e_Lim e_Zero
  {3.0f, -0.05f, -0.02f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1, 1,-0.5f}, // f_K#a ~GSh
   {-1, 1, 1},    //    #b ~GDi <opt>
   { 1, 1, 0}},   //    #c ~Thn <opt>
  {{-111,6789,  0,123,  3},  // f_C#a mskThn
   {   0,1234, 99,  0,104},  //    #b mskThk
   {-555,3333,222,  0,103},  //    #c mskGDi
   {  33, 456, 55,  0,  2}}, //    #d mskEdg
  {333, 0}, // h_Lim h_LimP
  {8, 2, 1, 44}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.2f, 0.2f}, // j_KMM j_ENm
  {2, 9, 0}, // l_RadHess l_RadMM l_Lim
  {0.1f, 0.03f, 2.0f}, // l_KEdg l_KRdg l_KMM
  {0.1f, 0.4f, 0.0f}, // n_Z
  {0.3f, 0.7f, 17.0, 0.6f}, // s_Za
  {0.27f, 0.5f, 11.0f}, // s_Zb
  {0.1f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS10 = { true, // lims
  {22}, // des_Thr
  {0, 0, 2.0f, 0.28f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {7}, // des_RadR
  {0.33f, -1, 3.0f}, // des_KR des_KU des_Kmf
  {0.0f, 0.0f, 0.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {6, 7, 4, 11}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.1f, 1.1f}, // des_KN des_KN2
  {-104, -123, 123}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 8, 3, 1, 17, 34}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {4.0f, 2.0f, 0.0f, 1.14f, 0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {2.0f, 4.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 3, 5}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {6, 22}, // d_Blk d_Blr
  {0.15f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {0, 0, 0, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.0f, 0.05f, 0.50f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 0.0f, 1.0f, 1.0f},  // f_K#a ~GSh
   {-6.0f, 1.0f, 0.0f},  //    #b ~GDi <opt>
   { 1.0f, 0.0f, 0.0f}}, //    #c ~Thn <opt>
  {{  789, 2345,  0,222,  1},  // f_C#a mskThn
   { 3333,  999, 55, 55,103},  //    #b mskThk
   {-3333,    0,  0,  0,102},  //    #c mskGDi
   {    0, 2678,  0,  0,  1}}, //    #d mskEdg
  {22,77}, // h_Lim h_LimP
  {14, 0, 0, 99}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.0f, 0.0f}, // j_KMM j_ENm
  {2, 3, 77}, // l_RadHess l_RadMM l_Lim
  {5.0f, 0.0f, 0.6f}, // l_KEdg l_KRdg l_KMM
  {0.0f, 0.0f, 0.0f}, // n_Z
  {0.3f, 0.4f, 28.0f, 0.7f}, // s_Za
  {0.7f, 0.0f, 9.9f}, // s_Zb
  {0.05f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS20 = { true, // torso lat (between skul lat & torso fr)
  {9}, // des_Thr
  {0, 0, 1.0f, 1.0f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {7}, // des_RadR
  {0.3f, -1.0f, 1.0f}, // des_KR des_KU des_Kmf
  {-0.5f, 0, 0.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {0, 0, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {0.0f, 0.0f}, // des_KN des_KN2
  {-103, -123, 123}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 8, 3, 0, 22, 44}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.0f, 2.0f, 0.0f, 1.4f, 0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {5.0f, 7.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {4, 6, 9}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {5, 33}, // d_Blk d_Blr; 5020:20, 5010:30
  {0.15f, 0.9f, 0.1f}, // d_BodyOne d_ExtMn'~>' d_ExtMx'~<'
  {7, 2, 888, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.3f, 0.01f, 0.1f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 1.0f,-1.0f},  // f_K#a ~GSh
   { 0.0f, 0.0f, 0.0f},  //    #b ~GDi <opt>
   { 0.0f, 0.0f, 0.0f}}, //    #c ~Thn <opt>
  {{   0,2345,  0,222,103},  // f_C#a mskThn
   { 999,   0,  0,  0,105},  //    #b mskThk
   { 111, 777,111,  0,104},  //    #c mskGDi
   { 111,1789,  0,  0,  2}}, //    #d mskEdg
  {-123, 99}, // h_Lim h_LimP
  {11, 0, 0, 44}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.0f, 0.1f}, // j_KMM j_ENm
  {2, 3, 33}, // l_RadHess lRadMM l_Lim
  {0.4f, 0.03f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.0f, 0.6f, 33.0f, 0.5f}, // s_Za
  {1.2f, 0.2f, 0.8f}, // s_Zb
  {0.05f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS22 = { true, // torso fr
  {9}, // des_Thr
  {0, 0, 1.0f, 0}, // des_Gamma des_Kgr des_KnP des_KSnP
  {5}, // des_RadR
  {0.7f, -1.0f, 2.0f}, // des_KR des_KU des_Kmf
  {0.3f, 0.03f, 3.0f, -0.85f}, // A+a* A+l* L+a* L+l*
  {6, 107, 7, 6}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.2f, 1.2f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 4, -1, 0, 22, 44}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {2.5f, 0.5f, 0, 1.0f, 0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {2.9f, 4.4f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {4, 6, 8}, // c_RadBas c_RadBig c_RadMM
  {31.0f, 0.3f, 0.5f}, // c_RefBig c_DelBig c_KMM; 5020:31, 5010:...
  {6, 20}, // d_Blk d_Blr; 5020:20, 5010:30
  {0.33f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {7, 2, 777, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.3f, 0.01f, 0.2f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 1.0f, 1.0f},  // f_K#a ~GSh
   { 1.0f, 1.0f, 0.0f},  //    #b ~GDi <opt>
   {-1.0f, 2.0f, 0.5f}}, //    #c ~Thn <opt>
  {{   0, 777,111,  0,  1},  // f_C#a mskThn
   {   0, 999,  0,  0,103},  //    #b mskThk
   {   0, 777,111,  0,102},  //    #c mskGDi
   {   0, 333,  0,  0,  2}}, //    #d mskEdg
  {99, 66}, // h_Lim h_LimP
  {7, 4, 0, 77}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.8f, 0.0f}, // j_KMM j_ENm
  {3, 4, 333}, // l_RadHess l_RadMM l_Lim
  {0.3f, 0.2f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.3f, 0.3f, 33.0f, 0.5f}, // s_Za
  {0.9f, 0.3f, 4.6f}, // s_Zb
  {0.05f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static tUnSharp cUS30 = { true, // skull lat (~= torso lat)
  {9}, // des_Thr
  {0, 0, 2.0f, 0}, // des_Gamma des_Kgr des_KnP0 des_KSnP
  {5}, // des_RadR
  {1.0f, 0.05f, 1.5f}, // des_KR des_KU des_Kmf
  {1.0f, 0.1f, 3.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {6, 107, 9, 9}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.4f, 1.2f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, -1, 0, 33, 55}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {3.0f, 0.5f, 0, 1.0f, 0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {3.33f, 4.44f, 4.44f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 3, 2}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.3f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {6, 11}, // d_Blk d_Blr; smaller Blr ~="20" or "-1" here
  {0.2f, 0.5f, 0.77f}, // d_BodyOne d_ExtMn d_ExtMx
  {6, 3, 256, 0}, // e_Blk e_Blr e_Lim e_Zero
  {1.0f, 0.05f, 0.3f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 0.5f, 1.0f,-0.3f},  // f_K#a ~GSh
   { 0.0f, 0.0f, 0.0f},  //    #b ~GDi <opt>
   { 0.0f, 0.0f, 0.0f}}, //    #c ~Thn <opt>
  {{   0,2222,  0,256,  2},  // f_C#a mskThn
   { 456,2345,  0,  0,  3},  //    #b mskThk
   {   0,2222,111, 99,102},  //    #c mskGDi
   { 222,1234,  0,  0,  2}}, //    #d mskEdg
  {123, 0}, // h_Lim h_LimP
  {9, 5, 1, 0}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.8f, 0.1f}, // j_KMM j_ENm
  {2, 5, 333}, // l_RadHess l_RadMM l_Lim
  {1.0f, 0.1f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.7f, 0.0f, 33.0f, 1.0f}, // s_Za
  {0.9f, 0.3f, 4.0f}, // s_Zb
  {0.05f,3E-3f,1E-6f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static tUnSharp cUS32 = { true, // skull fr
  {7}, // des_Thr
  {0, 0, 3.0f, 0}, // des_Gamma des_KnP des_KSnP
  {6}, // des_RadR
  {0.7f, 1.0f, 2.0f}, // des_KR des_KU des_Kmf
  {0.0f, 0.1f, 0.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {6, 107, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.1f, 1.1f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 4, -1, 0, 22, 33}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.5f, 0.5f, 0, 1.0f, 0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {7.0f, 4.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 4, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.44f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {6, 22}, // d_Blk d_Blr; smaller Blr ~="20" or "-1" here 
  {0.2f, 0.5f, 0.55f}, // d_BodyOne d_ExtMn d_ExtMx
  {7, 2, 888, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.5f, 0.01f, 0.1f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 0.5f, 1.0f,-0.5f},  // f_K#a ~GSh
   { 0.0f, 0.0f, 0.0f},  //    #b ~GDi <opt>
   { 0.0f, 0.0f, 0.0f}}, //    #c ~Thn <opt>
  {{   0,1234,222, 55,  1},  // f_C#a mskThn
   {   0,3456,  0,  0,102},  //    #b mskThk
   {   0,2345,256,  0,103},  //    #c mskGDi
   {   0, 999,  0,  0,  3}}, //    #d mskEdg
  {77, 88}, // h_Lim h_LimP
  {7, 4, 2, 88}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {2.0f, 0.05f}, // j_KMM j_ENm
  {2, 3, 888}, // l_RadHess l_RadMM l_Lim
  {0.5f, 0.05f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.8f, 0.0f, 33.0f, 4.0f}, // s_Za
  {1.23f, 0.4f, 7.4f}, // s_Zb
  {0.05f,3E-3f,1E-6f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS40 = { true, // chest,lung
  {27}, // des_Thr
  {5.0f, 0, 0, 0}, // des_Gamma des_Kgr des_KnP des_KSnP
  {6}, // des_RadR
  {1.8f, -1.0f, 1.0f}, // des_KR des_KU des_Kmf
  {-1.0f, 0, 0.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {0, 0, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {0.0f, 0.0f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {2, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, -1, 0, 99, 0}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {0.0f, 0.5f, 0, 1.0f, 0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {1.2f, 1.0f, 0.4f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 4, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.3f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {5, 33}, // d_Blk d_Blr
  {0.15f, 0.22f, 0.22f}, // d_BodyOne d_ExtMn d_ExtMx
  {0, 0, 0, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.0f, 0.0f, 0.0f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 0.6f, 0.3f},  // f_K#a ~GSh
   { 1.0f,-1.0f, 0.0f},  //    #b ~GDi <opt>
   { 0.0f, 1.0f,-1.0f}}, //    #c ~Thn <opt>
  {{  333,-567,  0,111,  1}, // f_C#a mskThn
   {  234,1234,  0, 55,102}, //    #b mskThk
   {-1234, 444,111,  0,103}, //    #c mskGDi
   {  999,   0,  0,  0,  1}},//    #d mskEdg
  {88, 0}, // h_Lim h_LimP
  {5, 0, 0, 0}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.0f, 0.0f}, // j_KMM j_ENm
  {1, 3, 55}, // l_RadHess l_RadMM l_Lim
  {0.05f, 0.01f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.7f, 0.3f, 8.0f, 0.5f}, // s_Za
  {1.0f, 0.2f, 1.5f}, // s_Zb
  {0.05f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS70 = { true, // pipe x1472x2172ndt/2050.raw ~#70
  {111}, // des_Thr
  {0.0f, 0, 0.66f, 0}, // des_Gamma des_Kgr des_KnP 
  {8}, // des_RadR
  {0.5f, 0.05f, 1.2f}, // des_KR des_KU des_Kmf
  {-0.5f, 0.0f, 2.0f, -1.2f}, // A+a* A+l* L+a* L+l*
  {11, 7, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.0f, 1.0f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 7, -1, 0, 33, 111}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {2.0f, 0.5f, 17.0f, -0.8f, 0.0f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {2.2f, 1.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 3, 6}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.4f}, // c_RefBig c_DelBig c_KMM
  {5, 44}, // d_Blk d_Blr
  {0.25f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {7, 7, 128, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.5f, 0.0f, 0.0f, 1.15f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1, 1,-2.5f},  // f_K#a ~GSh
   { 0, 0, 0},     //    #b ~GDi <opt>
   { 0, 0, 0}},    //    #c ~Thn <opt>
  {{   0, 777,222,  0,  2},  // f_C#a mskThn
   { 333, 777,  0,123,102},  //    #b mskThk
   { 444, 999, 55,  0,103},  //    #c mskGDi
   {   0, 333,  0,  0,  2}}, //    #d mskEdg
  {1,111}, // h_Lim h_LimP
  {11, 4, 0, 0}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {3.0f, 0.0f}, // j_KMM j_ENm
  {0, 5, 0}, // l_RadHess l_RadMM l_Lim
  {0.0f, 0.05f, 1.0f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.0f, 0.6f, 11.1f, 0.5f}, // s_Za
  {0.8f, 0.0f, 5.0f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS71 = { true, // PCB x2892x1920ndt/b03.raw 71 0888
  {-16}, // des_Thr
  {0, 0.7f, 0.7f, 0}, // des_Gamma des_Kgr des_KnP des_KSnP
  {4}, // des_RadR
  {0.3f, 0.1f, 1.0f}, // des_KR des_KU des_Kmf
  {-0.6f, 0.1f, 0.8f, -0.8f}, // A+a* A+l* L+a* L+l*
  {5, 106, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.2f, 1.0f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 5, -1, 15, 0, 0}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.0f, 0.5f, 0.0f, 1.25f, 0.0f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {1.0f, 1.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 1, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 1.0f}, // c_RefBig c_DelBig c_KMM
  {6, 33}, // d_Blk d_Blr
  {0.0f, 0.72f, 0.42f}, // d_BodyOne d_ExtMn d_ExtMx
  {6, 6, 111, 0}, // e_Blk e_Blr e_Lim e_Zero 
  {3.0f, -0.01f, -0.01f, 1.2f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 2.0f, 0.0f, 1.0f},  // f_K#a ~GSh
   { 1.0f,-0.5f, 0.0f},  //    #b ~GDi <opt>
   { 1.0f, 0.0f, 1.0f}}, //    #c ~Thn <opt>
  {{ 555,1234,  0,  0,102},  // f_C#a mskThn
   { 777,333,   0,  0,103},  //    #b mskThk
   { 444,1234, 99,  0,  3},  //    #c mskGDi
   { 345, 777,  0,  0,  2}}, //    #d mskEdg
  {111, 0}, // h_Lim h_LimP
  {6, 6, 0, 33}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.0f, 0.0f}, // j_KMM j_ENm
  {2, 5, 44}, // l_RadHess l_RadMM l_Lim
  {0.1f, 0.005f, 0.5f}, // l_KEdg l_KRdg l_KMM
  {0, 0.5f, 0.1f}, // n_Z
  {0.0f, 0.6f, 12.0f, 0.0f}, // s_Za
  {0.1f, 0.5f, 8.0f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

#if 0
static const tUnSharp cUS90 = { true, // MRI x0256/5o.raw ~#90 home adj
  {-333}, // des_Thr
  {0, 0.2f, 0.2f, 0.15f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {3}, // des_RadR
  {0.7f, 3.0f, 0.25f}, // des_KR des_KU des_Kmf
  {-0.5f, 0.0f, 2.0f, -0.8f}, // A+a* A+l* L+a* L+l*
  {5, 108, 4, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.2f, 1.2f}, // des_KN des_KN2
  {-103, -2222, 2222}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, 999, 234, 888, 345}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.0f, 10.0f, 0.0f, 1.1f, 0.0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {0.6f, 1.5f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 2, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {5, 55}, // d_Blk d_Blr
  {1.0f, 0.2f, 0.2f}, // d_BodyOne d_ExtMn d_ExtMx
  {0, 0, 0, 0}, // e_Blk e_Blr e_Lim e_Zero 
  {0.0f, 0.12f, -0.4f, 1.15f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 1.0f,-1.0f},  // f_K#a ~Msk
   {-1.0f,-1.0f, 2.2f},  //    #b ~GDi
   { 0.0f, 0.0f, 0.0f}}, //    #c ~Thn
  {{  222, 3000,  0,  0, 0},  // f_C#a mskThn
   {  222, 3000,  0,  0, 3},  //    #b mskThk
   {  222, 2000,  0,  0, 2},  //    #c mskGDi
   { -333,  999,  0,  0, 2}}, //    #d mskEdg
  {2345, 0}, // h_Lim h_LimP
  {3, 3, 1, 345}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {2.0f, 0.1f}, // j_KMM j_ENm
  {1, 3, 999}, // l_RadHess l_RadMM l_Lim
  {0.5f, 0.77f, 1.0f}, // l_KEdg l_KRdg l_KMM
  {-0.1f, 0.95f, -0.1f}, // n_Z
  {0.33f, 0.22f, 7.0f, 0.2f}, // s_Za
  {0.2f, 0.2f, 7.0f}, // s_Zb
  {0.05f,3E-4f,1E-4f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};
#else
static const tUnSharp cUS90 = { true, // MRI x0256/5o.raw ~#90 home adj
  {-333}, // des_Thr
  {0, 0.2f, 0.2f, 0.15f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {3}, // des_RadR
  {0.7f, 3.0f, 0.25f}, // des_KR des_KU des_Kmf
  {-0.5f, 0.0f, 2.0f, -0.8f}, // A+a* A+l* L+a* L+l*
  {5, 108, 4, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.2f, 1.2f}, // des_KN des_KN2
  {-103, -2222, 2222}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, 999, 234, 888, 345}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.0f, 10.0f, 0.0f, 1.1f, 0.0}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {0.6f, 1.5f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 2, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {5, 55}, // d_Blk d_Blr
  {1.0f, 0.2f, 0.2f}, // d_BodyOne d_ExtMn d_ExtMx
  {0, 0, 0, 0}, // e_Blk e_Blr e_Lim e_Zero 
  {0.0f, 0.12f, -0.4f, 1.15f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 1.0f,-1.0f},  // f_K#a ~Msk
   {-1.0f,-1.0f, 2.2f},  //    #b ~GDi
   { 0.0f, 0.0f, 0.0f}}, //    #c ~Thn
  {{  222, 3000,  0,  0, 0},  // f_C#a mskThn
   {  222, 3000,  0,  0, 3},  //    #b mskThk
   {  222, 2000,  0,  0, 2},  //    #c mskGDi
   { -333,  999,  0,  0, 2}}, //    #d mskEdg
  {2345, 0}, // h_Lim h_LimP
  {3, 3, 1, 345}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {2.0f, 0.1f}, // j_KMM j_ENm
  {1, 3, 999}, // l_RadHess l_RadMM l_Lim
  {0.5f, 0.77f, 1.0f}, // l_KEdg l_KRdg l_KMM
  {-0.1f, 0.95f, -0.1f}, // n_Z
  {0.33f, 0.22f, 7.0f, 0.2f}, // s_Za
  {0.2f, 0.2f, 7.0f}, // s_Zb
  {0.05f,3E-4f,1E-4f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};
#endif

static const tUnSharp cUS100 = { true, // x1024ifr/0000.raw ~#100
  {22}, // des_Thr
  {0.83f, 0, 1.4f, 0}, // des_Gamma des_Kgr des_KnP des_KSnP
  {3}, // des_RadR
  {1.0f, 0.05f, 1.0f}, // des_KR des_KU des_Kmf
  {0, 0, 0, 0}, // A+a* A+l* L+a* L+l*
  {4, 9, 12, 44}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.0f, 1.0f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, -1, 0, 555, 0}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {2.0f, 0.5f, 0, 1.0f, 0.0f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {0.6f, 0.8f, 0.2f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {2, 3, 4}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 1.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {6, 6}, // d_Blk d_Blr
  {0.0f, 0.2f, 0.2f}, // d_BodyOne d_ExtMn d_ExtMx
  {5, 2, 88, 0}, // e_Blk e_Blr e_Lim e_Zero
  {5.0f, 0.0f, 0.0f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 2.0f, 1.0f},  // f_K#a ~GSh
   {-1.0f, 3.0f,-1.0f},  //    #b ~GDi <opt>
   {-1.0f, 2.0f,-1.0f}}, //    #c ~Thn <opt>
  {{-111,3456,  0,123,  2},  // f_C#a mskThn
   {   0, 1234, 0,  0,102},  //    #b mskThk
   {1234,2345,111,  0,102},  //    #c mskGDi
   {   0, 400, 55,  0,  1}}, //    #d mskEdg
  {123, 0}, // h_Lim h_LimP
  {5, 4, 0, 0}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.6f, 0.0f}, // j_KMM j_ENm
  {1, 0, 123}, // l_RadHess l_RadMM l_Lim
  {0.5f, 0.2f, 0.0f}, // l_KEdg l_KRdg l_KMM
  {0.25f, 0.25f, 0.25f}, // n_Z
  {0.7f, 0.3f, 3.0f, 0.33f}, // s_Za
  {0.05f, 0.0f, 3.00f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS102 = { true, // x0324x256ifr/01210.raw ~#102
  {199}, // des_Thr
  {0, 0, 1.2f, 0}, // des_Gamma des_Kgr des_KnP des_KSnP
  {2}, // des_RadR
  {1.0f, 0.05f, 1.5f}, // des_KR des_KU des_Kmf
  {-0.5f, 0.0f, 1.0f, -1.0f}, // A+a* A+l* L+a* L+l*
  {3, 106, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.2f, 0.0f}, // des_KN des_KN2
  {0, 0, 0}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 2, -1, 0, 300, 500}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {1.2f, 0.5f, 0, -0.9f, 0.0f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {1.5f, 1.0f, 0.5f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 2, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 0.5f}, // c_RefBig c_DelBig c_KMM
  {5, 15}, // d_Blk d_Blr
  {0.0f, 0.2f, 0.2f}, // d_BodyOne d_ExtMn d_ExtMx
  {0, 0, 0, 0}, // e_Blk e_Blr e_Lim e_Zero
  {5.0f, 0.0f, 0.0f, 0.0f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1.0f, 2.0f, 1.0f},  // f_K#a ~GSh
   {-1.0f, 3.0f,-1.0f},  //    #b ~GDi <opt>
   {-1.0f, 2.0f,-1.0f}}, //    #c ~Thn <opt>
  {{-111,3456,  0,123,  2},  // f_C#a mskThn
   {   0,1234,  0,  0,  3},  //    #b mskThk
   {   0,2345,111,  0,  4},  //    #c mskGDi
   {   0, 400, 55,  0,  1}}, //    #d mskEdg
  {999, 0}, // h_Lim h_LimP
  {5, 3, 0, 99}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {1.0f, 0.0f}, // j_KMM j_ENm
  {1, 2, 234}, // l_RadHess l_RadMM l_Lim
  {2.2f, 0.0f, 0.5f}, // l_KEdg l_KRdg l_KMM`
  {0.25f, 0.25f, 0.25f}, // n_Z
  {0.5f, 0.5f, 3.5f, 0.1f}, // s_Za
  {0.33f, 0.1f, 2.0f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};

static const tUnSharp cUS200 = { true, // x1280x0960bar/001.raw ~#200
  {-267}, // des_Thr
  {0.0f, 3.0f, -0.75, 0.3f}, // des_Gamma des_Kgr des_KnP des_KSnP
  {5}, // des_RadR
  {0.6f, 1.5f, 2.0f}, // des_KR des_KU des_Kmf
  {-1.0f, 0.0f, 3.0f, 0.0f}, // A+a* A+l* L+a* L+l*
  {5, 105, 0, 0}, // des_RadN des_LgN des_RadN2 des_ZN2
  {1.27f, 0.0f}, // des_KN des_KN2
  {103, 111, -222}, // des_RadSnP des_r0SnP des_r1SnP
  {0, 0, 0}, // mgn_Msk mgn_Bev mgn_Cor
  {1, 3, 0, 111, 300, 500}, // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  {2.0f, 0.5f, 2.0f, 1.15f, 0.0f}, // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  {2.0f, 1.0f, 1.0f}, // b_GammaSpa b_GammaCLA b_GammaAtt
  {1, 2, 3}, // c_RadBas c_RadBig c_RadMM
  {0.0f, 0.0f, 1.0f}, // c_RefBig c_DelBig c_KMM
  {5, 22}, // d_Blk d_Blr
  {0.3f, 0.2f, 0.2f}, // d_BodyOne d_ExtMn d_ExtMx
  {0, 0, 0, 0}, // e_Blk e_Blr e_Lim e_Zero
  {0.0f, 0.0f, 0.2f, 1.1f}, // f_Att f_EmLthn f_EmLedg f_MxLbias
  {{ 1, 1,-0.2f},   // f_K#a ~GSh
   { 1,-1, 0.3f},   //    #b ~GDi <opt>
   { 0, 0, 0} },    //    #c ~Thn <opt>
  {{111 , 999, 0, 0, 1},  // f_C#a mskThn
   { 0,  1234, 0, 0, 2},  //    #b mskThk
   { 0,  2345, 0, 0, 1},  //    #c mskGDi
   { 0,  1234, 0, 0, 0}}, //    #d mskEdg
  {999, 0}, // h_Lim h_LimP
  {4, 2, 0, 99}, // j_RadMM j_RadMMO j_RadENm j_Lim
  {2.0f, 0.0f}, // j_KMM j_ENm
  {0, 4, 44}, // l_RadHess l_RadMM l_Lim
  {0.0f, 0.00f, 0.77f}, // l_KEdg l_KRdg l_KMM
  {0.25f, 0.25f, 0.25f}, // n_Z
  {0.5f, 0.0f, 6.6f, -0.5f}, // s_Za
  {0.15f, 0.1f, 2.2f}, // s_Zb
  {0.2f,3E-3f,1E-5f}, // u_BlOuDG u_Air0Mn u_Air1Mn
};
