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


// ----------------------
// Physical numbers
// ----------------------
OUTBND  = 1111;
AIR    = 1000;
IRON   = 2000;
iCOND1 = 4000;
istrandedCOND1 = 6000; // sTill
If(Flag_Transformer)// xfmr
  iCOND2 = 5000;
  istrandedCOND2 = 7000; // sTill
EndIf

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
  Winding1 =  Region[{}] ;
  Winding2 =  Region[{}] ;

  StrandedWinding1 =  Region[{}] ; // sTill
  StrandedWinding2 =  Region[{}] ; // sTill


  // Inductor + Transformer
  nbturns1 = NbrCond1/SymFactor; // Primary
  For iF In {1:nbturns1}
      Turn1~{iF} = Region[{(iCOND1+iF-1)}] ;
      Winding1  += Region[{(iCOND1+iF-1)}] ;
  EndFor
  For isF In {1:nbturns1}
      TurnStrand1~{isF} = Region[{(istrandedCOND1+isF-1)}] ;
      StrandedWinding1  += Region[{(istrandedCOND1+isF-1)}] ;
  EndFor


  // Transformer
  If(Flag_Transformer) // xfmr
    nbturns2 = NbrCond2/SymFactor; // Secondary
    For iF In {1:nbturns2}
        Turn2~{iF} = Region[{(iCOND2+iF-1)}] ;
        Winding2  += Region[{(iCOND2+iF-1)}] ;
    EndFor
    For isF In {1:nbturns2}
        TurnStrand2~{isF} = Region[{(istrandedCOND2+isF-1)}] ;
        StrandedWinding2  += Region[{(istrandedCOND2+isF-1)}] ;
    EndFor
  EndIf


  DomainC = Region[{Winding1, Winding2}] ;
  DomainS = Region[{StrandedWinding1, StrandedWinding2}] ;
  DomainCC += Region[{DomainS}] ; // sTill

  If(Flag_NL)
    Domain_Lin      = Region[{Air, Winding1, Winding2, StrandedWinding1, StrandedWinding2}]; // sTill
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Iron}];
  Else
    Domain_Lin      = Region[{Air, Iron, Winding1, Winding2, StrandedWinding1, StrandedWinding2}]; // sTill
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
  /*
  If(Flag_HomogenisedModel==1)
    Resistance_Cir += Region[{Zh~{1}}];
  EndIf
  */

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


  AreaCell[#{StrandedWinding1}] = AreaCell1;
  If(Flag_Transformer)
    AreaCell[#{StrandedWinding2}] = AreaCell2;
  EndIf
  AreaCell[#{Air, Iron, Winding1, Winding2}] = 1.;


  CoefGeo = 2*Pi*SymFactor ; // axisymmetry + symmetry factor

  sigma[#{Winding1, Winding2}] = SigmaCu ;
  sigma[#{Air, Iron}] = 0.;
  //rho[] = 1/sigma[];

  nu[#{Air}] = nu0;
  nu[#{Winding1, Winding2}] = nu0;

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

  // Excitation Current
  FSinusoidal1[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_1}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Fct_Src1[] = FSinusoidal1[];

  If(Flag_Transformer)
    FSinusoidal2[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_2}; //Complex_MH[1, 0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
    Fct_Src2[] = FSinusoidal2[];
  EndIf

  // Auxiliary functions for post-processing
  nuOm[#{Air}] = -nu[]*Complex[0.,1.];
  nuOm[#{Iron}] = -nu[$1]*Complex[0.,1.];
  nuOm[#{Winding1, Winding2}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ];


  // Skin Coefficient - will be multiplied with current "ir", which is zero except in windings
  kkk[#{Iron, Air, Winding1, Winding2}] =  1 ; // choose arbitrary value

  If(Flag_HomogenisedModel1)
    // Homogenization coefficients
    // Primary (Inductor + Transformer)
    file_ZSkinRe_1  = Sprintf("pre/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    file_ZSkinIm_1  = Sprintf("pre/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    file_NuProxRe_1= Sprintf("pre/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    file_NuProxIm_1 = Sprintf("pre/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    skin_rhor_list_1() = ListFromFile[ file_ZSkinRe_1 ];
    skin_rhoi_list_1() = ListFromFile[ file_ZSkinIm_1 ];
    prox_nur_list_1()  = ListFromFile[ file_NuProxRe_1 ];
    prox_nui_list_1()  = ListFromFile[ file_NuProxIm_1 ];
    skin_rhor_1[] = InterpolationLinear[$1]{ skin_rhor_list_1() };
    skin_rhoi_1[] = InterpolationLinear[$1]{ skin_rhoi_list_1() };
    prox_nur_1[]  = InterpolationLinear[$1]{ prox_nur_list_1() } ;
    prox_nui_1[]  = InterpolationLinear[$1]{ prox_nui_list_1() } ;
    nu[#{StrandedWinding1}] = nu0*Complex[prox_nur_1[Rr1], prox_nui_1[Rr1]*Fill1*Rr1^2/2];
    nuOm[#{StrandedWinding1}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill
    kkk[#{StrandedWinding1}] =  SymFactor * skin_rhor_1[Rr1] / SigmaCu / Fill1 ;
  EndIf

  If(Flag_Transformer)
    If(Flag_HomogenisedModel2)
      // Secondary
      file_ZSkinRe_2  = Sprintf("pre/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_ZSkinIm_2  = Sprintf("pre/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_NuProxRe_2= Sprintf("pre/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_NuProxIm_2 = Sprintf("pre/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      skin_rhor_list_2() = ListFromFile[ file_ZSkinRe_2 ];
      skin_rhoi_list_2() = ListFromFile[ file_ZSkinIm_2 ];
      prox_nur_list_2()  = ListFromFile[ file_NuProxRe_2 ];
      prox_nui_list_2()  = ListFromFile[ file_NuProxIm_2 ];
      skin_rhor_2[] = InterpolationLinear[$1]{ skin_rhor_list_2() };
      skin_rhoi_2[] = InterpolationLinear[$1]{ skin_rhoi_list_2() };
      prox_nur_2[]  = InterpolationLinear[$1]{ prox_nur_list_2() } ;
      prox_nui_2[]  = InterpolationLinear[$1]{ prox_nui_list_2() } ;
      // Formula from Paper:
      nu[#{StrandedWinding2}] = nu0*Complex[prox_nur_2[Rr2], prox_nui_2[Rr2]*Fill2*Rr2^2/2];
      nuOm[#{StrandedWinding2}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill
      kkk[#{StrandedWinding2}] =  SymFactor * skin_rhor_2[Rr2] / SigmaCu / Fill2 ;
    EndIf
  EndIf

  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  // List of nodes related to circuit
  // Inductor
  // Primary
  N1_1() = {1:nbturns1};   // Node 1 for each turn
  N1_2() = {2:nbturns1+1}; // Node 2 for each turn


  // Transformer
  If(Flag_Transformer)
    // Secondary
    N2_1() = {1:nbturns2};   // Node 1 for each turn
    N2_2() = {2:nbturns2+1}; // Node 2 for each turn
  EndIf

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
      // Inductor + Transformer
      If(Val_EE_1!=0)
          If(Flag_Circuit==0 && Flag_HomogenisedModel1==0)
            { Region Winding1 ; Value Val_EE_1; TimeFunction Fct_Src1[] ; }
          EndIf
          If(Flag_Circuit==0 && Flag_HomogenisedModel1==1)
            { Region StrandedWinding1 ; Value Val_EE_1; TimeFunction Fct_Src1[] ; }
          EndIf
      EndIf
      // Transformer
      If(Flag_Transformer)
        If(Val_EE_2!=0)
            If(Flag_Circuit==0 && Flag_HomogenisedModel2==0)
              { Region Winding2 ; Value Val_EE_2; TimeFunction Fct_Src2[] ; }
            EndIf
            If(Flag_Circuit==0 && Flag_HomogenisedModel2==1)
              { Region StrandedWinding2 ; Value Val_EE_2; TimeFunction Fct_Src2[] ; }
            EndIf
        EndIf
      EndIf
    }
  }


  // ---- Only For Voltage Excitation ----
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

  /*
  { Name ElectricalCircuit ; Type Network ;
    // Common to fine and homogenised models
    Case Circuit1 {
      If(Flag_HomogenisedModel1==0)
        { Region Input ;  Branch {N1_1(0), N1_2(nbturns1-1)} ; }
      EndIf

      If(Flag_HomogenisedModel1==1)
          { Region Input  ; Branch {777, N2(nbturns-1)} ; }
          { Region Zh~{1} ; Branch {777, N1(0)}; } // Complex impedance: Zskin
      EndIf

      For k In {0:nbturns1-1} // list indexes start at 0
        { Region Turn1~{k+1} ; Branch {N1_1(k), N1_2(k)} ; }
      EndFor
    }
  }
  // ---- Only For Voltage Excitation ----
  */
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
      { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue ; Frequency Freq ; }
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


      PostOperation[Map_local] ;
      PostOperation[Get_global] ;

    }// Operation
  }
}// Resolution


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

      //{ Name normal ; Value { Term { [ Normal[] ] ; In Domain ; Jacobian Vol ; } } }
      //{ Name tangent ; Value { Term { [ Tangent[] ] ; In Domain ; Jacobian Vol ; } } }

      /*
      //{ Name tmp ; Value {
      //    Term { [  nu[{d a}] * {d a}  ]; In Domain ; Jacobian Vol ; }
      //} }


      // Error Estimations for adaptive mesh refinement
      { Name integrated_error ; Value {
          Integral { [ SquNorm[ -1 * { d tmp } - sigma[]*(Dt[{a}]+{ur}/CoefGeo)] +
                       SquNorm[ -1 * { d tmp } - sigma[]*(Dt[{a}]+{ur}/CoefGeo)] +
                       SquNorm[ Cross [ Normal[], nu[{d a}] * {d a} ] ]
                     ];
          In Domain ; Jacobian Vol ; Integration II ;
      } } }

      { Name errora ; Value {
          Term { [   Sqrt [ SquNorm[  ElementVol[] * ( Cross [ Normal[], nu[{d a}] * {d a} ] ) ]
                     ] ];
          In Domain ; Jacobian Vol ;
      } } }

      { Name error ; Value {
          Term { [   Sqrt [ SquNorm[ ElementVol[] * ( -1 * { d tmp } - sigma[]*(Dt[{a}]+{ur}/CoefGeo) ) ] +
                            //SquNorm[ ElementVol[] * ( sigma[]*( Dt[{a}] + {ur} ) * Normal[] ) ] +
                            SquNorm[  ElementVol[] * ( Cross [ Normal[], nu[{d a}] * {d a} ] ) ]
                     ] ];
          In Domain ; Jacobian Vol ;
      } } }

      { Name errorb ; Value {
          Term { [   Sqrt [ SquNorm[ ( -1 * { d tmp } - sigma[]*(Dt[{a}]+{ur}/CoefGeo) ) ] +
                            //SquNorm[ ( sigma[]*( Dt[{a}] + {ur} ) * Normal[] ) ] +
                            SquNorm[ ( Cross [ Normal[], nu[{d a}] * {d a} ] ) ]
                     ] ];
          In Domain ; Jacobian Vol ; Integration II ;
      } } }

      { Name error1 ; Value {
          Term { [ Norm[ -1 * { d tmp } - sigma[]*(Dt[{a}]+{ur}/CoefGeo)] ];
          In Domain ; Jacobian Vol ; Integration II ;
      } } }

      //{ Name error2 ; Value {
      //    Term { [ Norm[ {Div jc} ] ];
      //    In DomainC ; Jacobian Vol ; Integration II ;
      //} } }

      { Name jc ; Value { Term { [ {jc} ] ; In Domain ; Jacobian Vol  ;} } }

      { Name error4 ; Value {
          Term { [ ElementVol[] * Norm[ Cross [ Normal[], nu[{d a}] * {d a} ] ] ];
          In Domain ; Jacobian Vol ; Integration II ;
      } } }

      { Name error5 ; Value {
          Term { [ Norm[ sigma[]*( Dt[{a}] + {ur} ) * Normal[] ] ];
          In Domain ; Jacobian Vol ; Integration II ;
      } } }
      */

      { Name j ; Value {
        Term { [ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ; In DomainC ; Jacobian Vol ; }
        Term { [ -1/AreaCell[]*{ir} ] ; In DomainS ; Jacobian Vol ; }
      } }
      { Name jz ; Value {
          Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
          Term { [ CompZ[ -1/AreaCell[]*{ir} ] ] ; In DomainS ; Jacobian Vol ; }
      } }

      { Name e ; Value {
          Term { [ -1*(Dt[{a}]+{ur}/CoefGeo) ] ; In Domain ; Jacobian Vol ; }
      } }

      { Name MagEz ; Value {
          Term { [ Norm [ -1*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In Domain ; Jacobian Vol ; }
      } }


      { Name J_rms ; Value {
          Term { [ Norm[ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ] ; In DomainC ; Jacobian Vol ; }
          Term { [ Norm[ CompZ[ {ir}/AreaCell[] ] ] ] ; In DomainS ; Jacobian Vol ; }
          Term { [ Norm[ 0 ] ] ; In Domain_Lin_NoJs ; Jacobian Vol ; } //non_lin case is missing
      } }

      //{ Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell[] ] ;
      //      In Domain ; Jacobian Vol  ; Integration II ; } } }


      If(Freq==0.0)
          { Name j2F ; Value { Integral { // Joule losses
             [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2F ; Value { Integral { // Joule losses
             [ 0.5*CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      EndIf


      // together
      If(Freq==0.0)
           { Name j2H ; Value { Integral { // Joule losses
             [ CoefGeo*( Re[{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
             In DomainS ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2H ; Value { Integral { // Joule losses
             [ 0.5*CoefGeo*( Re[{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 added by Till
             In DomainS ; Jacobian Vol ; Integration II ; } } }
      EndIf


      { Name j2Hprox ; Value { Integral { // Joule losses
            [ 0.5*CoefGeo*Re[{d a}*Conj[nuOm[]*{d a}]] ] ;// 0.5 added by Till
            In DomainS ; Jacobian Vol ; Integration II ; } } }

      { Name j2Hskin ; Value { Integral { // Joule losses
            [ 0.5*CoefGeo*kkk[]*SquNorm[-1/AreaCell[]*{ir}] ] ;// 0.5 added by Till
            In DomainS ; Jacobian Vol ; Integration II ; } } }


      { Name b2F ; Value { Integral { // Magnetic Energy
            [ CoefGeo*nu[{d a}]*SquNorm[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }


      { Name SoF ; Value {
          Integral {
            [ CoefGeo * Complex[ sigma[] * SquNorm[(Dt[{a}]+{ur}/CoefGeo)], nu[{d a}]*SquNorm[{d a}] ] ] ; // added {d a} in nu[...]
            In DomainC ; Jacobian Vol ; Integration II ; }
        } }//Complex power
        // together: Domain -> DomainC



      { Name SoH ; Value { Integral { // Complex power = Active power +j * Reactive power => S = P+j*Q
            [ CoefGeo * ({d a}*Conj[nuOm[{d a}]*{d a}] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
            In DomainC ; Jacobian Vol ; Integration II ; } } } //Complex power
            // xfmr changed Domain to DomainC
            // to prevent from div by zero error in "Air" and "Core" domains


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

      { Name ElectEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }


      { Name Flux ; Value {
          Integral { [ NbrCond1/AreaCell1 * Norm[{a}] ] ;
            In Domain  ; Jacobian Vol ; Integration II ; } } }

//      { Name Flux ; Value {
//          Integral { [ SymmetryFactor*Lz*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
//            In Inds  ; Jacobian Vol ; Integration I1 ; } } }

/*
      { Name Inductance_from_Flux ;
        Value { Term { Type Global; [ $Flux / Val_EE_1 ] ; In Domain ; } } } // xfmr

//      { Name Inductance_from_Flux ; Value {
//          Integral { [ NbrCond1/AreaCell1 * CompZ[{a}] / Val_EE_1 ] ;
//            In Winding1  ; Jacobian Vol ; Integration II ; } } }
*/

      { Name Inductance_from_MagEnergy ;
        Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE_1*Val_EE_1) ] ; In Domain ; } } } // xfmr

    If(Val_EE_1!=0)
         { Name L_11 ;
            Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE_1*Val_EE_1) ] ; In Domain ; } } } // xfmr
    EndIf

     If(Flag_Transformer)// xfmr

     //    { Name L_12 ;
     //       Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE_1*Val_EE_2) ] ; In Domain ; } } } // xfmr

     //    { Name L_21 ;
     //      Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE_2*Val_EE_1) ] ; In Domain ; } } } // xfmr

        If(Val_EE_2!=0)
            { Name L_22 ;
               Value { Term { Type Global; [ 2 * $MagEnergy * 1/(Val_EE_2*Val_EE_2) ] ; In Domain ; } } } // xfmr
        EndIf
     EndIf

    }
  }
}


// ----------------------
// output of post-processing quantities
// printing in gmsh with a .pos file or saving in a .dat file
// ----------------------

// -- No Homogenisation --
PostOperation Map_local UsingPost MagDyn_a {
  //Print[ raz,OnElementsOf Domain,  File StrCat[DirRes, "a", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ur,OnElementsOf Domain,  File StrCat[DirRes, "ur", ExtGmsh], LastTimeStepOnly ] ;

  //Print[ e,  OnElementsOf Domain,  File StrCat[DirRes, "e", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ MagEz,  OnElementsOf Domain,  File StrCat[DirRes, "MagEz", ExtGmsh],  LastTimeStepOnly ] ;

  //Print[ b,  OnElementsOf Domain,  File StrCat[DirRes, "b", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ h,  OnElementsOf Domain,  File StrCat[DirRes, "h", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ Magb,  OnElementsOf Domain,  File StrCat[DirRes, "Magb", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ Flux,  OnElementsOf Domain,  File StrCat[DirRes, "Flux_field", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ MagEnergy,  OnElementsOf Domain,  File StrCat[DirRes, "MagEnergy", ExtGmsh],  LastTimeStepOnly ] ;
  /*
  Print[ error, OnElementsOf Region[{Domain}], File StrCat[DirRes, "error", ExtGmsh], LastTimeStepOnly ] ;
  Print[ error1, OnElementsOf Region[{Domain}], File StrCat[DirRes, "error1", ExtGmsh], LastTimeStepOnly ] ;
  Print[ jc, OnElementsOf Region[{Domain}], File StrCat[DirRes, "jc", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ error2, OnElementsOf Region[{Domain}], File StrCat[DirRes, "error2", ExtGmsh], LastTimeStepOnly ] ;
  Print[ error4, OnElementsOf Region[{Domain}], File StrCat[DirRes, "error4", ExtGmsh], LastTimeStepOnly ] ;
  Print[ error5, OnElementsOf Region[{Domain}], File StrCat[DirRes, "error5", ExtGmsh], LastTimeStepOnly ] ;
  Print[ normal, OnElementsOf Region[{Domain}], File StrCat[DirRes, "normal", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ tangent, OnElementsOf Region[{Domain}], File StrCat[DirRes, "tangent", ExtGmsh], LastTimeStepOnly ] ;
  */
  //Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "jz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ j, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "j", ExtGmsh], LastTimeStepOnly ] ;
  Print[ J_rms, OnElementsOf Region[{Domain}], File StrCat[DirRes, "J_rms", ExtGmsh], LastTimeStepOnly ] ;
  Print[ j2F, OnElementsOf Region[{DomainC}], File StrCat[DirRes, "j2F", ExtGmsh], LastTimeStepOnly ] ;

  // together Losses
  Print[ j2H,   OnElementsOf DomainS, File StrCat[DirRes,"jH",ExtGmsh] ] ;
  //try Print[ j2Hprox,   OnElementsOf DomainS, File StrCat[DirRes,"jHprox",ExtGmsh] ] ;
  //try Print[ j2Hskin,   OnElementsOf DomainS, File StrCat[DirRes,"jHskin",ExtGmsh] ] ;
  //Print[ nuOm,   OnElementsOf DomainS, File StrCat[DirRes,"nuOm",ExtGmsh] ] ;


  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding1, Format TimeTable, File >Sprintf("res/I_f%g.dat", Freq), LastTimeStepOnly];
    Print[ U, OnRegion Winding1, Format TimeTable, File >Sprintf("res/U_f%g.dat", Freq), LastTimeStepOnly];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("res/I_f%g.dat", Freq), LastTimeStepOnly];
    Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("res/U_f%g.dat", Freq), LastTimeStepOnly];
  EndIf

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[DirRes,"ME.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];
  //Print[ MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/MagEnergy.dat") ];
  //Print[ Inductance_from_MagEnergy, OnRegion Domain, Format Table, File Sprintf("res/Inductance.dat") ];

  If(Val_EE_1!=0)
    //Print[ L_11, OnRegion Domain, Format Table, File Sprintf("res/L_11.dat") ];
    Print[ L_11, OnRegion Domain, Format TimeTable, File > Sprintf("res/L_11.dat")] ;
  EndIf
  If(Flag_Transformer)// xfmr
      //Print[ L_12, OnRegion Domain, Format Table, File Sprintf("res/L_12.dat") ];
      //Print[ L_21, OnRegion Domain, Format Table, File Sprintf("res/L_21.dat") ];
      If(Val_EE_2!=0)
        //Print[ L_22, OnRegion Domain, Format Table, File Sprintf("res/L_22.dat")] ;
        Print[ L_22, OnRegion Domain, Format TimeTable, File > Sprintf("res/L_22.dat")] ;

      EndIf
  EndIf
  //Print[ Inductance_from_Flux, OnRegion Domain, Format Table, File Sprintf("res/Inductance_from_Flux.dat") ];
  //Print[ Flux[Winding1], OnGlobal, Format Table, File Sprintf("res/Flux.dat") ];
}

PostOperation Get_global UsingPost MagDyn_a {
  Print[ j2F[ DomainC ], OnGlobal, Format TimeTable, File > Sprintf("res/j2F_iron.dat")] ;// Joule losses
  //try Print[ SoF[ DomainC ], OnGlobal, Format TimeTable,  File > Sprintf("res/SF_iron.dat")] ; // Complex power
  Print[ j2H[ DomainS ], OnGlobal, Format TimeTable, File > StrCat["res/j2H.dat"] ] ;
  //Print[ SoH[ DomainS ], OnGlobal, Format TimeTable, File Sprintf("res/SH_f%g.dat", Freq) ] ;

}

