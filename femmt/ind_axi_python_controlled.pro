// ----------------------
// Files and Directories
// ----------------------
Include "Parameter.pro";
Include "BH.pro";
ExtGmsh = ".pos";
DirRes = "res/";


// ----------------------
// All Variables - remove or create in python
// ----------------------
Nb_max_iter = 100;
relaxation_factor = 1.;
stop_criterion = 1e-8;
SymFactor = 1. ; //half inductor with axisymmetry
Flag_Circuit = Flag_ImposedVoltage;


NbrTurns = 1; // does not make sense for this way of homogenisation
//Len = 1000*(2*Pi*(Xw1+Xw2)/2)*NbrTurns ; // average length
AreaCond = Pi*Rc^2;
//Dex = 2.2*Rc;
//Dey = Dex; // only valid for squared packing


//AreaCell = 1; //Dex*Dey;  // what is the meaning of this factor?
//Rdc = Len/SigmaCu/AreaCond; // fill factor is missing


// ----------------------
// Physical numbers
// ----------------------
OUTBND  = 1111;
AIR    = 1000;
IRON   = 2000;
iCOND = 4000;
istrandedCOND = 6000; // sTill


// ----------------------
// Groups
// ----------------------
Group{
  // ----------------------
  // Physical Domains
  // ----------------------
  Air  = Region[{AIR}];
  Iron = Region[{IRON}];
  DomainCC = Region[{Air, Iron}] ; // Non Conducting Domain
  OuterBoundary = Region[{OUTBND}]; // including symmetry


  // Current Conducting Domains
  Winding =  Region[{}] ;
  StrandedWinding =  Region[{}] ; // sTill


  nbturns = NbrCond/SymFactor; // number of turns
  For iF In {1:nbturns}
      Turn~{iF} = Region[{(iCOND+iF-1)}] ;
      Winding  += Region[{(iCOND+iF-1)}] ;
  EndFor

  nbturns = NbrstrandedCond/SymFactor; // number of turns // sTill
  For isF In {1:nbturns} // sTill
      TurnStrand~{isF} = Region[{(istrandedCOND+isF-1)}] ; // sTill
      StrandedWinding  += Region[{(istrandedCOND+isF-1)}] ; // sTill
  EndFor // sTill

  DomainC = Region[{Winding}] ;
  DomainS = Region[{StrandedWinding}] ;
  DomainCC += Region[{DomainS}] ; // sTill

  If(Flag_NL)
    Domain_Lin      = Region[{Air, Winding, StrandedWinding}]; // sTill
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Iron}];
  Else
    Domain_Lin      = Region[{Air, Winding, StrandedWinding, Iron}]; // sTill
    Domain_Lin_NoJs = Region[{Air, Iron}];
    Domain_NonLin   = Region[{}];
  EndIf

  Domain = Region[{DomainC, DomainCC}] ;

  DomainDummy = Region[ 12345 ] ; // Dummy region number for postpro with functions


  // ----------------------
  // Circuit Domains
  // ----------------------
  Input = # 12345 ;
  iZH    = 10000 ;
  iLH    = 20000 ; // Zskin in homog. coil
  iLHp   = 30000 ;


  Zh~{1} = Region[{(iZH+1)}];
  Lh~{1} = Region[{(iLH+1)}];


  Resistance_Cir  = Region[{}];
  If(Flag_HomogenisedModel==1)
    Resistance_Cir += Region[{Zh~{1}}];
  EndIf

  Inductance_Cir  = Region[{}];
  Capacitance1_Cir = Region[ {} ] ;
  Capacitance2_Cir = Region[ {} ] ;
  Capacitance_Cir = Region[ {Capacitance1_Cir, Capacitance2_Cir} ] ;

  SourceV_Cir = Region[ {Input} ] ;
  SourceI_Cir = Region[ {} ] ;

  DomainZ_Cir = Region[ {Resistance_Cir, Inductance_Cir, Capacitance_Cir} ] ;

  DomainSource_Cir = Region[ {SourceV_Cir, SourceI_Cir} ] ;
  DomainZt_Cir = Region[ {DomainZ_Cir, DomainSource_Cir} ] ;
}


// ----------------------
// Excitation
// ----------------------
Function {

  CoefGeo = 2*Pi*SymFactor ; // axisymmetry + symmetry factor

  sigma[#{Winding, StrandedWinding}] = SigmaCu ;  // sTill
  sigma[#{Air, Iron}] = 0.;
  rho[] = 1/sigma[];

  nu[#{Air}] = nu0;

  If (!Flag_NL)
    nu[#{Iron}]   = nu0/mur;
  Else
    //nu[ #{Iron} ] = nu_3kW[$1] ;
    //h[ #{Iron} ]  = h_3kW[$1];
    //dhdb_NL[ #{Iron} ]= dhdb_3kW_NL[$1] ;
    //dhdb[ #{Iron} ]   = dhdb_3kW[$1] ;


    //nu[ #{Iron} ] = nu_N95[$1] ;
    nu[ #{Iron} ] = nu~{Core_Material}[$1] ;
    h[ #{Iron} ]  = h~{Core_Material}[$1];
    dhdb_NL[ #{Iron} ]= dhdb_95_NL[$1] ;
    dhdb[ #{Iron} ]   = dhdb~{Core_Material}[$1] ;

  EndIf


  FSinusoidal[] = Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Fct_Src[] = FSinusoidal[];



  // new files, finer mesh, max X=8
  la = Fill ;
  file_ZSkinRe  = Sprintf("pre/coeff/pI_RS_la%.2g_%.2glayer.dat", la, NbrLayers);
  file_ZSkinIm  = Sprintf("pre/coeff/qI_RS_la%.2g_%.2glayer.dat", la, NbrLayers);
  file_NuProxRe = Sprintf("pre/coeff/qB_RS_la%.2g_%.2glayer.dat", la, NbrLayers);
  file_NuProxIm = Sprintf("pre/coeff/pB_RS_la%.2g_%.2glayer.dat", la, NbrLayers);

  // new files, finer mesh, max X=8
  // la = 0.65; // fill factor
  // file_ZSkinRe  = Sprintf("coeff/pI_RS_la%.2g.dat", la);
  // file_ZSkinIm  = Sprintf("coeff/qI_RS_la%.2g.dat", la);
  // file_NuProxRe = Sprintf("coeff/qB_RS_la%.2g.dat", la);
  // file_NuProxIm = Sprintf("coeff/pB_RS_la%.2g.dat", la);

  // Homogenization coefficients: round conductor & square packing
  // Frequency domain
  skin_rhor_list() = ListFromFile[ file_ZSkinRe ];
  skin_rhoi_list() = ListFromFile[ file_ZSkinIm ];
  prox_nur_list()  = ListFromFile[ file_NuProxRe ];
  prox_nui_list()  = ListFromFile[ file_NuProxIm ];

  skin_rhor[] = InterpolationLinear[$1]{ skin_rhor_list() };
  skin_rhoi[] = InterpolationLinear[$1]{ skin_rhoi_list() };

  prox_nur[]  = InterpolationLinear[$1]{ prox_nur_list() } ;
  prox_nui[]  = InterpolationLinear[$1]{ prox_nui_list() } ;


  nu[Winding] = nu0;

  // Formula from Paper:
  nu[StrandedWinding] = nu0*Complex[prox_nur[Rr], prox_nui[Rr]*Fill*Rr^2/2];


  // used if Flag_HomogenisedModel==1
  // Skin effect => complex impedance
  // Formula from Paper:
  //Zskin[] = 1/SymFactor*Complex[ skin_rhor[Rr]*Rdc, 2*Pi*Freq*skin_rhoi[Rr]*mu0*Len/(8*Pi*Fill)] ;

  // Auxiliary functions for post-processing
  nuOm[#{Air}] = -nu[]*Complex[0.,1.];
  nuOm[#{Iron}] = -nu[$1]*Complex[0.,1.];
  nuOm[#{Winding}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ];
  nuOm[#{StrandedWinding}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill

  kkk[] =  SymFactor * skin_rhor[Rr] / Fill / SigmaCu ;




  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  Ns[] = 1; // NbrCond/SymFactor ;
  // Sc[] =  SurfaceArea[]{iCOND};   // sTill

  //If (Flag_HomogenisedModel==1)
    // Accounting for eddy currents (homogenization)
    //Resistance[Zh~{1}] = Zskin[] ;
  //EndIf

  // List of nodes related to circuit
  N1() = {1:nbturns};   // Node 1 for each turn
  N2() = {2:nbturns+1}; // Node 2 for each turn


  Sc[] =  SurfaceArea[]{istrandedCOND};   // sTill


}


// ----------------------
// definition of constraints
// ----------------------
Constraint {

  { Name MVP_2D ;
    Case {
      { Region OuterBoundary ; Type Assign ;  Value 0. ; }
    }
  }

  // Massive/stranded conductor constraints
  { Name Current_2D ;
    Case {
      If(Flag_Circuit==0)
        { Region Winding ; Value Val_EE; TimeFunction Fct_Src[] ; }
      EndIf

      If(Flag_Circuit==0 && Flag_HomogenisedModel)
        { Region StrandedWinding ; Value Val_EE; TimeFunction Fct_Src[] ; }
      EndIf

    }
  }

  { Name Voltage_2D ;
    Case{
    }
  }

  { Name Voltage_Cir ;
    Case {
      If(Flag_Circuit && Flag_ImposedVoltage)
        { Region Input ; Value Val_EE; TimeFunction Fct_Src[] ; }
      EndIf
    }
  }
  { Name Current_Cir ;
    Case {
      If(Flag_Circuit && !Flag_ImposedVoltage)
        { Region Input ; Value -Val_EE; TimeFunction Fct_Src[] ; }
      EndIf
    }
  }

  { Name ElectricalCircuit ; Type Network ;
    // Common to fine and homogenised models
    Case Circuit1 {
      If(Flag_HomogenisedModel==0)
        { Region Input ;  Branch {N1(0), N2(nbturns-1)} ; }
      EndIf

      If(Flag_HomogenisedModel==1)
          { Region Input  ; Branch {777, N2(nbturns-1)} ; }
          { Region Zh~{1} ; Branch {777, N1(0)}; } // Complex impedance: Zskin
      EndIf

      For k In {0:nbturns-1} // list indexes start at 0
        { Region Turn~{k+1} ; Branch {N1(k), N2(k)} ; }
      EndFor
    }
  }
}


// ----------------------
// Include solver file(s)
// ----------------------
Include "solver.pro"


// ----------------------
// call of the chosen problem formulation
// ----------------------
Resolution {

  { Name Analysis ;
    System {
      If(Flag_HomogenisedModel==0)
        { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue ; Frequency Freq ; }
      EndIf

      If(Flag_HomogenisedModel==1)
        { Name A ; NameOfFormulation MagDyn_a_Homog ; Type ComplexValue ; Frequency Freq ; }
      EndIf
    }

    Operation {
      CreateDir[DirRes];

      If(!Flag_NL)
          Generate[A] ; Solve[A] ;
          Else
          IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor]{
              GenerateJac[A] ; SolveJac[A] ;
          }
      EndIf
      SaveSolution[A] ;

      If(Flag_HomogenisedModel==0)
          PostOperation[Map_local] ;
          PostOperation[Get_global] ;
      Else
          PostOperation[Map_local_homog] ;
          PostOperation[Get_global_homog] ;
      EndIf
    }// Operation
  }
}// Reslution


// ----------------------
// calculation of varoius post-processing quantities with the help of the mag. vec.pot.
// ----------------------
PostProcessing {

  { Name MagDyn_a ; NameOfFormulation MagDyn_a ; NameOfSystem A;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }

      { Name ur ; Value { Term { [ {ur} ] ; In Domain  ; Jacobian Vol  ;} } }

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name Magb ; Value { Term { [ Norm[ {d a} ] ]; In Domain ; Jacobian Vol ; } } }
      { Name h ; Value { Term { [ nu[{d a}]*{d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name j ; Value { Term { [ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ; In DomainC ; Jacobian Vol ; } } }
      { Name jz ; Value {
          Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
          //Term { [ CompZ[ {ir}/AreaCond ] ] ; In DomainS ; Jacobian Vol ; }
        } }

      { Name J_rms ; Value {
          Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
          //Term { [ CompZ[ {ir}/AreaCond ] ] ; In DomainS ; Jacobian Vol ; }
      } }

      //{ Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell ] ;
      //      In Domain ; Jacobian Vol  ; Integration II ; } } }


      { Name j2F ; Value { Integral { // Joule losses
            [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
            In DomainC ; Jacobian Vol ; Integration II ; } } }

      { Name b2F ; Value { Integral { // Magnetic Energy
            [ CoefGeo*nu[{d a}]*SquNorm[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }

      { Name SoF ; Value {
          Integral {
            [ CoefGeo * Complex[ sigma[] * SquNorm[(Dt[{a}]+{ur}/CoefGeo)], nu[{d a}]*SquNorm[{d a}] ] ] ; // added {d a} in nu[...]
            In Domain ; Jacobian Vol ; Integration II ; }
        } }//Complex power

      { Name U ; Value {
          Term { [ {U} ]   ; In DomainC ; }
          Term { [ {Us} ]   ; In DomainS ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } }
      { Name I ; Value {
          Term { [ {I} ]   ; In DomainC ; }
          Term { [ {Is} ]   ; In DomainS ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        } }

      { Name MagEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }


      //{ Name Flux ; Value {
      //    Integral { [ Lz*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
      //      In Inds  ; Jacobian Vol ; Integration II ; } } }

      //{ Name Inductance_from_Flux ; Value { Term { Type Global; [ $Flux * 1/Val_EE ] ;
      //  In DomainDummy ; } } }

      { Name Inductance_from_MagEnergy ;
        Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE*Val_EE) ] ; In Domain ; } } }
     }
  }


  { Name MagDyn_a_Homog ; NameOfFormulation MagDyn_a_Homog ;
    PostQuantity {
      { Name a   ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az  ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name Magb ; Value { Term { [ Norm[ {d a} ] ]; In Domain ; Jacobian Vol ; } } }

      { Name j  ; Value { Term { [ -1/AreaCell*{ir} ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz ; Value { Term { [ CompZ[ -1/AreaCell*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }


      { Name j2H ; Value { Integral { // Joule losses
            [ CoefGeo*( Re[{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell*{ir}]) ] ;
            In DomainS ; Jacobian Vol ; Integration II ; } } }

      { Name nuOm ; Value { Term { [ nuOm[] ] ; In DomainS ; Jacobian Vol ; } } }

      //{ Name j2H ; Value { Integral { // Joule losses
      //      [ CoefGeo*(Norm[Re[{d a}*Conj[nuOm[]*{d a}]]] + kkk[]*SquNorm[{ir}]) ] ;
      //      In DomainS ; Jacobian Vol ; Integration II ; } } }

      //{ Name j2F ; Value { Integral { // Joule losses
      //      [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
      //      In DomainC ; Jacobian Vol ; Integration II ; } } }

      { Name j2Hprox ; Value { Integral { // Joule losses
            [ CoefGeo*Re[{d a}*Conj[nuOm[]*{d a}]] ] ;
            In DomainS ; Jacobian Vol ; Integration II ; } } }

      { Name j2Hskin ; Value { Integral { // Joule losses
            [ CoefGeo*kkk[]*SquNorm[-1/AreaCell*{ir}] ] ;
            In DomainS ; Jacobian Vol ; Integration II ; } } }

      { Name SoH ; Value { Integral { // Complex power = Active power +j * Reactive power => S = P+j*Q
            [ CoefGeo * ({d a}*Conj[nuOm[{d a}]*{d a}] + kkk[]*SquNorm[-1/AreaCell*{ir}]) ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } } //Complex power

      { Name MagEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }

      { Name Inductance_from_MagEnergy ;
        Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE*Val_EE) ] ; In Domain ; } } }


      { Name U ; Value {
          Term { [ {Us} ]  ; In DomainS ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } }
      { Name I ; Value {
          Term { [ {Is} ]  ; In DomainS ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        } }

      { Name conjI ; Value {
          Term { [ Conj[{Is}] ]  ; In DomainS ; }
          Term { [ Conj[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
      { Name squnormI ; Value {
          Term { [ SquNorm[{Is}] ]  ; In DomainS ; }
          Term { [ SquNorm[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
      { Name reI ; Value {
          Term { [ Re[{Is}] ]  ; In DomainS ; }
          Term { [ Re[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
      { Name imI ; Value {
          Term { [ Im[{Is}] ]  ; In DomainS ; }
          Term { [ Im[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
    }
  }
  // sTill
}


// ----------------------
// output of post-processing quantities
// printing in gmsh with a .pos file or saving in a .dat file
// ----------------------

// -- No Homogenisation --
PostOperation Map_local UsingPost MagDyn_a {
  Print[ raz,OnElementsOf Domain,  File StrCat[DirRes, "a", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ur,OnElementsOf Domain,  File StrCat[DirRes, "ur", ExtGmsh], LastTimeStepOnly ] ;

  //Print[ b,  OnElementsOf Domain,  File StrCat[DirRes, "b", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ h,  OnElementsOf Domain,  File StrCat[DirRes, "h", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ Magb,  OnElementsOf Domain,  File StrCat[DirRes, "Magb", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ MagEnergy,  OnElementsOf Domain,  File StrCat[DirRes, "MagEnergy", ExtGmsh],  LastTimeStepOnly ] ;

  Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "jz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ J_rms, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "J_rms", ExtGmsh], LastTimeStepOnly ] ;
  Print[ j2F, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "j2F", ExtGmsh], LastTimeStepOnly ] ;

  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File >Sprintf("res/I_f%g.dat", Freq), LastTimeStepOnly];
    Print[ U, OnRegion Winding, Format TimeTable, File >Sprintf("res/U_f%g.dat", Freq), LastTimeStepOnly];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("res/I_f%g.dat", Freq), LastTimeStepOnly];
    Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("res/U_f%g.dat", Freq), LastTimeStepOnly];
  EndIf

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[DirRes,"ME.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];
  //Print[ MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/MagEnergy.dat") ];
  Print[ Inductance_from_MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/Inductance.dat") ];
}

PostOperation Get_global UsingPost MagDyn_a {
  Print[ j2F[ Winding ], OnGlobal, Format TimeTable, File > Sprintf("res/j2F_iron.dat")] ;// Joule losses
  Print[ SoF[ Domain ], OnGlobal, Format TimeTable,  File > Sprintf("res/SF_iron.dat")] ; // Complex power
}


// -- Homogenisation --
  PostOperation Map_local_homog UsingPost MagDyn_a_Homog {
  // Magnetic Vector Potential and Magnetic Fields
  Print[ raz, OnElementsOf Domain, File StrCat[DirRes,"aH",ExtGmsh] ] ;
  Print[ b,   OnElementsOf Domain, File StrCat[DirRes,"bH",ExtGmsh] ] ;
  Print[ Magb,   OnElementsOf Domain, File StrCat[DirRes,"MagbH",ExtGmsh] ] ;

  // Current Density
  Print[ j,   OnElementsOf DomainS, File StrCat[DirRes,"j",ExtGmsh] ] ;
  Print[ jz,   OnElementsOf DomainS, File StrCat[DirRes,"jz",ExtGmsh] ] ;

  // Losses
  Print[ j2H,   OnElementsOf DomainS, File StrCat[DirRes,"jH",ExtGmsh] ] ;
  Print[ j2Hprox,   OnElementsOf DomainS, File StrCat[DirRes,"jHprox",ExtGmsh] ] ;
  Print[ j2Hskin,   OnElementsOf DomainS, File StrCat[DirRes,"jHskin",ExtGmsh] ] ;
  Print[ nuOm,   OnElementsOf DomainS, File StrCat[DirRes,"nuOm",ExtGmsh] ] ;

  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  }

  PostOperation Get_global_homog UsingPost MagDyn_a_Homog {
  // Complex power: S = P + i*Q, P=active power, Q=reactive power

  Print[ SoH[Domain], OnGlobal, Format TimeTable, File Sprintf("res/SH_f%g.dat", Freq) ] ;
  // Print[ j2H[StrandedWinding], OnGlobal, Format Table, File Sprintf("res/j2H_f%g.dat", Freq) ] ;
  Print[ j2H[StrandedWinding], OnGlobal, Format TimeTable, File > StrCat["res/j2H.dat"] ] ;
  //Print[ Inductance_from_MagEnergy, OnGlobal, Format Table, File Sprintf("res/Inductance.dat") ];
  //Print[ Inductance_from_MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/Inductance.dat") ];

  If(!Flag_Circuit)
    Print[ U, OnRegion StrandedWinding, Format Table, File Sprintf("res/U_f%g.dat", Freq) ] ;
    Print[ I, OnRegion StrandedWinding, Format Table, File Sprintf("res/I_f%g.dat", Freq) ] ;
  Else
    Print[ U, OnRegion Input, Format Table, File Sprintf("res/U_f%g.dat", Freq) ];
    Print[ I, OnRegion Input, Format Table, File Sprintf("res/I_f%g.dat", Freq) ];
  EndIf
  }


// ---- Not Used ----
PostOperation Get_allTS UsingPost MagDyn_a {
  // Print[ Inductance_from_MagEnergy, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
  //  File StrCat[DirRes,"Inductance"];

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File Sprintf("res/I_f%g.dat", Freq) ];
    Print[ U, OnRegion Winding, Format TimeTable, File Sprintf("res/U_f%g.dat", Freq) ];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File Sprintf("res/I_f%g.dat", Freq),
      SendToServer StrCat[po,"I [A]"]{0}, Color "Pink"];
    Print[ U, OnRegion Input, Format TimeTable, File Sprintf("res/U_f%g.dat", Freq),
      SendToServer StrCat[po,"V [V]"]{0}, Color "LightGreen"];
    Print[ I, OnRegion Turn~{nbturns} , Format TimeTable, File Sprintf("res/Iw_f%g.dat", Freq) ];
  EndIf

  Print[ j2F[Winding], OnGlobal, Format TimeTable, File Sprintf("res/jl_f%g.dat", Freq) ] ; // Joule losses
}