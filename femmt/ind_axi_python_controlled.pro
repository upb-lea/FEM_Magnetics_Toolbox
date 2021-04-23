// ----------------------
// Included Files
// ----------------------
Include "Parameter.pro";
Include "BH.pro";
// Include "BH_ind.pro";


// ----------------------
// Directories
// ----------------------
ExtGmsh = ".pos";
DirRes = "res/";
po = "{Output/";


// ----------------------
// All Variables - remove or create in python
// ----------------------
// new files, finer mesh, max X=8
// file_ZSkinRe  = Sprintf("coeff/pI_RS_la%.2g.dat", la);
// file_ZSkinIm  = Sprintf("coeff/qI_RS_la%.2g.dat", la);
// file_NuProxRe = Sprintf("coeff/qB_RS_la%.2g.dat", la);
// file_NuProxIm = Sprintf("coeff/pB_RS_la%.2g.dat", la);
Zh~{1} = 100;
Flag_FD = 1;
Flag_IronCore = 1;
Nb_max_iter = 100;
iter_max = Nb_max_iter;
relaxation_factor = 1.;
stop_criterion = 1e-8;
SymFactor = 1. ; //half inductor with axisymmetry
Rc = Sqrt[1/Pi]*1e-3; // needs review Till: Conductor radius
NbrTurns = 1;
Len = (2*Pi*(Xw1+Xw2)/2)*NbrTurns ;
AreaCond = Pi*Rc^2;
Dex = 2.2*Rc;
Dey = Dex;
//AreaCell = Dex*Dey;
AreaCell = 0.008*0.028;
Rdc = Len/SigmaCu/AreaCond;
Fill = 0.9;
la = Fill; // fill factor


// ----------------------
// Declaration of further quantities from Parameter.pro file
// ----------------------
Flag_imposedFreq = !Flag_imposedRr;
Omega  = 2*Pi*Freq;
Period = 1./Freq;
Flag_Circuit = Flag_ImposedVoltage;


// ----------------------
// Physical numbers
// ----------------------
OUTBND  = 1111;
AIR    = 1000;
IRON   = 2000;
iCOND = 4000;


// ----------------------
// Groups
// ----------------------
Group{
  // ----------------------
  // Physical Domains
  // ----------------------
  Air  = Region[{AIR}];

  If(Flag_IronCore)
    Iron = Region[{IRON}];
  Else
    Iron = Region[{}];
    Air  += Region[{IRON}];
  EndIf

  OuterBoundary = Region[{OUTBND}]; // including symmetry
  Winding =  Region[{}] ;
  DomainCC = Region[{Air, Iron}] ;
  nbturns = NbrCond/SymFactor; // number of turns
  For iF In {1:nbturns}
      Turn~{iF} = Region[{(iCOND+iF-1)}] ;
      Winding  += Region[{(iCOND+iF-1)}] ;
  EndFor

  DomainC = Region[{Winding}] ;
  DomainS = Region[{}] ;
  //?  DomainCC += Region[{DomainS}] ;

  If(Flag_NL)
    Domain_Lin      = Region[{Air, Winding}];
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Iron}];
  Else
    Domain_Lin      = Region[{Air, Winding, Iron}];
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

  Resistance_Cir  = Region[{}];
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

  sigma[#{Winding}] = SigmaCu ;
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

//Sprintf("res/U_f%g.dat", Core_Material)

    //nu[ #{Iron} ] = nu_N95[$1] ;
    nu[ #{Iron} ] = nu~{Core_Material}[$1] ;
    h[ #{Iron} ]  = h_95[$1];
    dhdb_NL[ #{Iron} ]= dhdb_95_NL[$1] ;
    dhdb[ #{Iron} ]   = dhdb_95[$1] ;

  EndIf


  FSinusoidal[] = Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Fct_Src[] = FSinusoidal[];


  nu[Winding] = nu0;


  //If(Flag_FD) // used if Flag_HomogenisedModel==1
    // Skin effect => complex impedance
  //  Zskin[] = 1/SymFactor*Complex[ skin_rhor[Rr]*Rdc, 2*Pi*Freq*skin_rhoi[Rr]*mu0*Len/(8*Pi*Fill)] ;
  //EndIf
  // Auxiliary functions for post-processing
  nuOm[#{Air}] = -nu[]*Complex[0.,1.];
  nuOm[#{Iron}] = -nu[$1]*Complex[0.,1.];
  nuOm[#{Winding}] = Complex[ Omega * Im[nu[]], -Re[nu[]] ];

  //kkk[] =  SymFactor * skin_rhor[Rr] / Fill /SigmaCu ;




  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  Ns[] = 1; // NbrCond/SymFactor ;
  Sc[] =  SurfaceArea[]{iCOND};

  // List of nodes related to circuit
  N1() = {1:nbturns};   // Node 1 for each turn
  N2() = {2:nbturns+1}; // Node 2 for each turn
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
      { Region Input ;  Branch {N1(0), N2(nbturns-1)} ; }

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

    If(Flag_FD) // Frequency domain
        { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue ; Frequency Freq ; }
    EndIf
    }


  //  Operation {
  //    CreateDir[DirRes];
  //    InitSolution[A]; SaveSolution[A];
  //
  //    If(Flag_FD) // Frequency domain
  //      Generate[A] ; Solve[A] ; SaveSolution[A] ;
  //      PostOperation [Map_local];
  //      PostOperation [Get_global];
  //    EndIf
  //  }


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

      PostOperation[Map_local] ;
      PostOperation[Get_global] ;
    }
  }
}


// ----------------------
// calculation of varoius post-processing quantities with the help of the mag. vec.pot.
// ----------------------
PostProcessing {

  { Name MagDyn_a ; NameOfFormulation MagDyn_a ; NameOfSystem A;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }

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

      { Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell ] ;
            In Domain ; Jacobian Vol  ; Integration II ; } } }


      { Name j2F ; Value { Integral { // Joule losses
            [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
            In DomainC ; Jacobian Vol ; Integration II ; } } }

      { Name b2F ; Value { Integral { // Magnetic Energy
            [ CoefGeo*nu[{d a}]*SquNorm[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }

      { Name SoF ; Value {
          Integral {
            [ CoefGeo * Complex[ sigma[] * SquNorm[(Dt[{a}]+{ur}/CoefGeo)], nu[]*SquNorm[{d a}] ] ] ;
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
}


// ----------------------
// output of post-processing quantities
// printing in gmsh with a .pos file or saving in a .dat file
// ----------------------
PostOperation Map_local UsingPost MagDyn_a {
  Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "jz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ J_rms, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "J_rms", ExtGmsh], LastTimeStepOnly ] ;
  Print[ j2F, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "j2F", ExtGmsh], LastTimeStepOnly ] ;
  Print[ b,  OnElementsOf Domain,  File StrCat[DirRes, "b", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ h,  OnElementsOf Domain,  File StrCat[DirRes, "h", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ Magb,  OnElementsOf Domain,  File StrCat[DirRes, "Magb", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ MagEnergy,  OnElementsOf Domain,  File StrCat[DirRes, "MagEnergy", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ raz,OnElementsOf Domain,  File StrCat[DirRes, "a", ExtGmsh], LastTimeStepOnly ] ;

  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File >Sprintf("res/I_f%g.dat", Freq),
      SendToServer StrCat[po,"I [A]"], Color "Pink", LastTimeStepOnly];
    Print[ U, OnRegion Winding, Format TimeTable, File >Sprintf("res/U_f%g.dat", Freq),
      SendToServer StrCat[po,"V [V]"], Color "LightGreen", LastTimeStepOnly];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("res/I_f%g.dat", Freq),
      SendToServer StrCat[po,"I [A]"], Color "Pink", LastTimeStepOnly];
    Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("res/U_f%g.dat", Freq),
      SendToServer StrCat[po,"V [V]"], Color "LightGreen", LastTimeStepOnly];
  EndIf

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[DirRes,"ME.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];


  //Print[ MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/MagEnergy.dat") ];
  Print[ Inductance_from_MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/Inductance.dat") ];

}

PostOperation Get_global UsingPost MagDyn_a {
  Print[ j2F[ Winding ], OnGlobal, Format TimeTable, File > Sprintf("res/j2F_iron%g.dat", Flag_IronCore)] ;// Joule losses
  Print[ SoF[ Domain ], OnGlobal, Format TimeTable,  File > Sprintf("res/SF_iron%g.dat", Flag_IronCore)] ; // Complex power
}



// ---- Not Used ----
PostOperation Get_allTS UsingPost MagDyn_a {
  //Print[ Inductance_from_MagEnergy, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
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
