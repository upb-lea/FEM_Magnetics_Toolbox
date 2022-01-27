// ----------------------
// Files and Directories
// ----------------------
Include "../Parameter.pro";
Include "BH.pro";
Include "mu_imag.pro";
ExtGmsh = ".pos";



// ----------------------
// All Variables - remove or create in python
// ----------------------
Nb_max_iter             = 20;
relaxation_factor       = 1.;
stop_criterion          = 1e-8;
Flag_Circuit            = Flag_ImposedVoltage;
// ----------------------
// half inductor with axisymmetry
// 1 means full zylinder
SymFactor               = 1. ;
CoefGeo                 = 2*Pi*SymFactor ; // axisymmetry +/* symmetry factor
e_0                     = 8.8541878128e-12;


Flag_Conducting_Core    = 1;
sigma_core              = e_r_imag * 2*Pi*Freq * e_0;

// ----------------------
// Physical numbers
// ----------------------
OUTBND              = 1111;
AIR                 = 1000;
AIR_EXT             = 1001;
IRON                = 2000;
iCOND1              = 4000;
istrandedCOND1      = 6000;
If(Flag_Transformer)
  iCOND2            = 5000;
  istrandedCOND2    = 7000;
EndIf



// ----------------------
// Groups
// ----------------------
Group{
  // ----------------------
  // Physical Domains
  // ----------------------
  Air  = Region[{AIR, AIR_EXT}];
  Iron = Region[{IRON}];

  // Non Conducting Domain:
  DomainCC = Region[{Air}];
  If(!Flag_Conducting_Core)
    DomainCC += Region[{Iron}];
  EndIf

  // Boundary Conditions
  // including symmetry
  OuterBoundary = Region[{OUTBND}];

  // Current Conducting Domains
  Winding1 =  Region[{}] ;
  Winding2 =  Region[{}] ;
  StrandedWinding1 =  Region[{}] ;
  StrandedWinding2 =  Region[{}] ;


  // Primary (Inductor + Transformer)
  nbturns1 = NbrCond1/SymFactor;
  For iF In {1:nbturns1}
      Turn1~{iF} = Region[{(iCOND1+iF-1)}] ;
      Winding1  += Region[{(iCOND1+iF-1)}] ;
  EndFor
  For isF In {1:nbturns1}
      TurnStrand1~{isF} = Region[{(istrandedCOND1+isF-1)}] ;
      StrandedWinding1  += Region[{(istrandedCOND1+isF-1)}] ;
  EndFor


  // Secondary (Transformer)
  If(Flag_Transformer)
    nbturns2 = NbrCond2/SymFactor;
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
  If(Flag_Conducting_Core)
    DomainC += Region[{Iron}] ;
  EndIf
  DomainS = Region[{StrandedWinding1, StrandedWinding2}] ;
  DomainCC += Region[{DomainS}] ;

  If(Flag_NL)
    Domain_Lin      = Region[{Air, Winding1, Winding2, StrandedWinding1, StrandedWinding2}];
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Iron}];
  Else
    Domain_Lin      = Region[{Air, Iron, Winding1, Winding2, StrandedWinding1, StrandedWinding2}];
    Domain_Lin_NoJs = Region[{Air, Iron}];
    Domain_NonLin   = Region[{}];
  EndIf

  Domain = Region[{DomainC, DomainCC}] ;


  DomainCond1 = Region[{Winding1, StrandedWinding1}];
  DomainCond2 = Region[{Winding2, StrandedWinding2}];

  // Dummy region number for postpro with functions
  DomainDummy = Region[ 12345 ] ;


  // ----------------------
  // Circuit Domains
  // ----------------------
  Input  = # 12345 ;
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

  // Strand sizes
  AreaCell[#{StrandedWinding1}] = AreaCell1;
  If(Flag_Transformer)
    AreaCell[#{StrandedWinding2}] = AreaCell2;
  EndIf
  // in non-stranded domains, def. AreaCell to 1 (neutral element of mutiplication)
  AreaCell[#{Air, Iron, Winding1, Winding2}] = 1.;


  // Material Properties

  // sigma: conductivity (= imaginary part of complex permitivity)
  //rho[] = 1/sigma[];
  sigma[#{Winding1, Winding2}] = SigmaCu ;
  If(Flag_Conducting_Core)
    sigma[#{Iron}] = sigma_core;
  EndIf
  If(!Flag_Conducting_Core)
    sigma[#{Air, Iron}] = 0.;
  EndIf

  // nu: reluctivity
  // nu = 1/mu
  nu[#{Air}] = nu0;
  nu[#{Winding1, Winding2}] = nu0;

  // Hysteresis Loss
  // Imaginary Part Of Permeability
  // Liste von Lukas hinterlegen
  mu_imag[ #{Iron} ] = mu0 * f_N95_mu_imag[$1, $2];
  //mu_imag[#{Iron}] = mu0 * mu_r_imag;

  If(!Flag_NL)
    nu[#{Iron}]   = nu0/mur;
  Else
    //nu[ #{Iron} ] = nu_3kW[$1] ;
    //h[ #{Iron} ]  = h_3kW[$1];
    //dhdb_NL[ #{Iron} ]= dhdb_3kW_NL[$1] ;
    //dhdb[ #{Iron} ]   = dhdb_3kW[$1] ;


    //nu[ #{Iron} ] = nu_N95[$1] ;
    nu[ #{Iron} ] = nu~{Core_Material}[$1] ;
    h[ #{Iron} ]  = h~{Core_Material}[$1];
    //dhdb_NL[ #{Iron} ]= dhdb_95_NL[$1] ;
    dhdb_NL[ #{Iron} ]= dhdb_95100_NL[$1] ;
    dhdb[ #{Iron} ]   = dhdb~{Core_Material}[$1] ;
  EndIf

  // Excitation Current
  FSinusoidal1[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_1}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Fct_Src1[] = FSinusoidal1[];
  Sign1 = (Phase_1==Pi) ? -1 : 1;

  If(Flag_Transformer)
    FSinusoidal2[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_2}; //Complex_MH[1, 0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
    Fct_Src2[] = FSinusoidal2[];
    Sign2 = (Phase_2==Pi) ? -1 : 1;
  EndIf

  // Auxiliary functions for post-processing
  nuOm[#{Air}] = -nu[]*Complex[0.,1.];
  nuOm[#{Iron}] = -nu[$1]*Complex[0.,1.];
  //nuOm[#{Winding1, Winding2}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ];


  // Skin Coefficient - will be multiplied with current "ir", which is zero except in windings
  kkk[#{Iron, Air, Winding1, Winding2}] =  1 ; // choose arbitrary value

  If(Flag_HomogenisedModel1)
    // Homogenization coefficients
    // Primary (Inductor + Transformer)
    file_ZSkinRe_1  = Sprintf("Strands_Coefficents/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    file_ZSkinIm_1  = Sprintf("Strands_Coefficents/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    file_NuProxRe_1= Sprintf("Strands_Coefficents/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
    file_NuProxIm_1 = Sprintf("Strands_Coefficents/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill1, NbrLayers1);
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
      file_ZSkinRe_2  = Sprintf("Strands_Coefficents/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_ZSkinIm_2  = Sprintf("Strands_Coefficents/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_NuProxRe_2= Sprintf("Strands_Coefficents/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_NuProxIm_2 = Sprintf("Strands_Coefficents/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
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
  // Imprinted Currents
  { Name Current_2D ;
    Case {
      // Inductor + Transformer
      //If(Val_EE_1!=0)
      If(1)
          If(Flag_Circuit==0 && Flag_HomogenisedModel1==0)
            { Region Winding1 ; Value Val_EE_1; TimeFunction Fct_Src1[] ; }
          EndIf
          If(Flag_Circuit==0 && Flag_HomogenisedModel1==1)
            { Region StrandedWinding1 ; Value Val_EE_1; TimeFunction Fct_Src1[] ; }
          EndIf
      EndIf
      // Transformer
      If(Flag_Transformer)
        //If(Val_EE_2!=0)
        If(1)
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
      If(Flag_Conducting_Core)
        { Region Iron ; Value 0; }
      EndIf
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
      CreateDir[DirResFields];
      CreateDir[DirResVals];
      CreateDir[DirResValsPrimary];
      CreateDir[DirResValsSecondary];
      CreateDir[DirResCirc];

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

      // Potentials
      { Name a ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name ur ; Value { Term { [ {ur} ] ; In Domain  ; Jacobian Vol  ;} } }

      // Electrical Field
      { Name e ; Value { Term { [ -1*(Dt[{a}]+{ur}/CoefGeo) ] ; In Domain ; Jacobian Vol ; } } }
      { Name MagEz ; Value { Term { [ Norm [ -1*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In Domain ; Jacobian Vol ; } } }

      // Magnetic Field
      { Name h ; Value { Term { [ nu[{d a}]*{d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name Magh ; Value { Term { [ Norm[ nu[{d a}]*{d a} ] ] ; In Domain ; Jacobian Vol ; } } }

      // Magnetic Flux Density
      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name b_pol ; Value { Term { [ Norm [ Re [ Cart2Pol[ {d a} ] ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name im_b_pol ; Value { Term { [ Norm [ Im [ Cart2Pol[ {d a} ] ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name Magb ; Value { Term { [ Norm[ {d a} ] ]; In Domain ; Jacobian Vol ; } } }

      // Core Loss
      // for piecewise linear currents
      // iGSE Integral explicitely solved
      // needs the result of Magb at peak current to evaluate the peak flux density
      // (Norm[{d a}]*2) is delta_B


      If(Flag_Generalized_Steinmetz_loss)
        { Name piGSE ; Value { Integral { [ Freq * ki * (Norm[{d a}]*2)^(beta-alpha) * (
                                        ((Norm[{d a}]*2 / t_rise )^alpha) * t_rise +
                                        ((Norm[{d a}]*2 / t_fall )^alpha) * t_fall
                                        // 10 abschnitte reinbauen
                                        // python überprüfung + vorfaktoren zu NULL
                                   ) ] ; In Iron ; Jacobian Vol ; Integration II ;} } }
      EndIf

      If(Flag_Steinmetz_loss)
        { Name pSE ; Value { Integral { [ CoefGeo * ki * Freq^alpha * (Norm[{d a}])^beta
                                     ] ; In Iron ; Jacobian Vol ; Integration II ;} } }

        { Name pSE_density ; Value { Integral { [ CoefGeo* ki * Freq^alpha * (Norm[{d a}])^beta
                                     ] ; In Iron ; Jacobian Vol ; Integration II ;} } }
      EndIf

      // Permeability
      { Name mur ; Value { Term { [ 1 / Norm[ nu[{d a}] / mu0 ] ] ; In Domain ; Jacobian Vol ; } } }

      // Magnetic Flux Density
      { Name p_hyst ; Value { Integral {
        [ 0.5 * CoefGeo * 2*Pi*Freq * mu_imag[{d a}, Freq] * SquNorm[nu[{d a}]*{d a}] ] ;
        In Iron ; Jacobian Vol ; Integration II ;} } }
      { Name p_hyst_density ; Value { Integral {
        [ 0.5 * CoefGeo/ElementVol[] * 2*Pi*Freq * mu_imag[{d a}, Freq] * SquNorm[nu[{d a}]*{d a}] ] ;
        In Iron ; Jacobian Vol ; Integration II ;} } }

      // Energy
      { Name MagEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }
      { Name ElectEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }
      //{ Name b2F ; Value { Integral { // Magnetic Energy
      //      [ CoefGeo*nu[{d a}]*SquNorm[{d a}] ] ;
      //     In Domain ; Jacobian Vol ; Integration II ; } } }
      //{ Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell[] ] ;
      //      In Domain ; Jacobian Vol  ; Integration II ; } } }

      // Current Density
      { Name j ; Value {
            Term { [ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ; In DomainC ; Jacobian Vol ; }
            Term { [ -1/AreaCell[]*{ir} ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz ; Value {
            Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
            Term { [ CompZ[ -1/AreaCell[]*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }
      { Name J_rms ; Value {
            Term { [ Norm[ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ] ; In DomainC ; Jacobian Vol ; }
            Term { [ Norm[ CompZ[ {ir}/AreaCell[] ] ] ] ; In DomainS ; Jacobian Vol ; }
            Term { [ Norm[ 0 ] ] ; In Domain_Lin_NoJs ; Jacobian Vol ; } } } //non_lin case is missing

      // Ohmic Loss
      // DomainC (Solid Conductors)
      If(Freq==0.0)
          { Name j2F ; Value { Integral {
             [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)]] ;
             In DomainC ; Jacobian Vol ; Integration II ; } } }
          { Name j2F_density ; Value { Integral {
             [ CoefGeo/ElementVol[]*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }

      Else
           { Name j2F ; Value { Integral {
             [ 0.5*CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }

           { Name j2F_density ; Value { Integral {
             [ 0.5*CoefGeo/ElementVol[]*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      EndIf

      // DomainS (Stranded Conductors)
      If(Freq==0.0)
           { Name j2H ; Value { Integral {
             [ CoefGeo*( Re[{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
             In DomainS ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2H ; Value { Integral {
             [ 0.5*CoefGeo*( Re[{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 added
             In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2H_density ; Value { Integral {
             [ 0.5*CoefGeo/ElementVol[]*( Re[{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 added
             In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2Hprox ; Value { Integral {
            [ 0.5*CoefGeo*Re[{d a}*Conj[nuOm[]*{d a}]] ] ;// 0.5 added by Till
            In DomainS ; Jacobian Vol ; Integration II ; } } }
           { Name j2Hskin ; Value { Integral {
            [ 0.5*CoefGeo*kkk[]*SquNorm[-1/AreaCell[]*{ir}] ] ;// 0.5 added by Till
            In DomainS ; Jacobian Vol ; Integration II ; } } }
      EndIf

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

      // Voltage (Voltage_i = dFlux_Linkage_i / dt)
      { Name Voltage_1 ; Value {
        Integral { [ CoefGeo / AreaCell1 * CompZ[Dt[{a}]] ]; In DomainCond1; Jacobian Vol; Integration II; } } }
      If(Flag_Transformer)
        { Name Voltage_2 ; Value {
          Integral { [ CoefGeo / AreaCell2 * CompZ[Dt[{a}]] ]; In DomainCond2; Jacobian Vol; Integration II; } } }
      EndIf

      // Flux (Linkage)
      { Name Flux_Linkage_1 ; Value {
        Integral { [ CoefGeo / AreaCell1 * CompZ[{a}] ]; In DomainCond1; Jacobian Vol; Integration II; } } }
      If(Flag_Transformer)
        { Name Flux_Linkage_2 ; Value {
          Integral { [ CoefGeo / AreaCell2 * CompZ[{a}] ]; In DomainCond2; Jacobian Vol; Integration II; } } }
      EndIf

      // (Self) Inductances
      If(Val_EE_1!=0)
        { Name L_11 ; Value { Integral {
          [ Sign1 * CoefGeo / AreaCell1 * CompZ[{a}] / Val_EE_1 ]; In DomainCond1; Jacobian Vol; Integration II; } } }
        { Name L_11_from_MagEnergy ; Value { Integral {
          [ 2 * CoefGeo*nu[{d a}]*({d a}*{d a})/2 / (Val_EE_1*Val_EE_1) ]; In Domain; Jacobian Vol; Integration II; } } }
      EndIf

      If(Flag_Transformer)
        If(Val_EE_2!=0)
          { Name L_22 ; Value { Integral {
            [ Sign2 * CoefGeo / AreaCell2 * CompZ[{a}] / Val_EE_2 ]; In DomainCond2; Jacobian Vol; Integration II; } } }
          { Name L_22_from_MagEnergy ; Value { Integral {
            [ 2 * CoefGeo*nu[{d a}]*({d a}*{d a})/2 / (Val_EE_2*Val_EE_2) ]; In Domain; Jacobian Vol; Integration II; } } }
        EndIf
     EndIf


      // Circuit Quantities
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
    }
  }
}


// ----------------------
// output of post-processing quantities
// printing fields in gmsh with a .pos file or saving results in a .dat file
// ----------------------


// === Field Quantities ===

PostOperation Map_local UsingPost MagDyn_a {

  // Potentials
  //Print[ raz,OnElementsOf Domain,  File StrCat[DirResFields, "raz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ a,OnElementsOf Domain,  File StrCat[DirResFields, "a", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ur,OnElementsOf Domain,  File StrCat[DirResFields, "ur", ExtGmsh], LastTimeStepOnly ] ;

  // Electrical Field
  //Print[ e,  OnElementsOf Domain,  File StrCat[DirResFields, "e", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ MagEz,  OnElementsOf Domain,  File StrCat[DirResFields, "MagEz", ExtGmsh],  LastTimeStepOnly ] ;

  // Magnetic Field
  //Print[ h,  OnElementsOf Domain,  File StrCat[DirResFields, "h", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ Magh,  OnElementsOf Domain,  File StrCat[DirResFields, "Magh", ExtGmsh],  LastTimeStepOnly ] ;

  // Core Loss Density
  If(Flag_Generalized_Steinmetz_loss)
    Print[ piGSE,  OnElementsOf Domain,  File StrCat[DirResFields, "piGSE", ExtGmsh],  LastTimeStepOnly ] ;
  EndIf

  If(Flag_Steinmetz_loss)
    Print[ pSE,  OnElementsOf Domain,  File StrCat[DirResFields, "pSE", ExtGmsh],  LastTimeStepOnly ] ;
    Print[ pSE_density,  OnElementsOf Domain,  File StrCat[DirResFields, "pSE_density", ExtGmsh],  LastTimeStepOnly ] ;
  EndIf

  //Print[ mur,  OnElementsOf Domain,  File StrCat[DirResFields, "mur", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ p_hyst,  OnElementsOf Domain,  File StrCat[DirResFields, "p_hyst", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ p_hyst_density,  OnElementsOf Domain,  File StrCat[DirResFields, "p_hyst_density", ExtGmsh],  LastTimeStepOnly ] ;

  // Magnetic Flux (Density)
  Print[ b,  OnElementsOf Domain,  File StrCat[DirResFields, "b", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ b_pol,  OnElementsOf Domain,  File StrCat[DirResFields, "b_pol", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ im_b_pol,  OnElementsOf Domain,  File StrCat[DirResFields, "im_b_pol", ExtGmsh],  LastTimeStepOnly ] ;
  If(Flag_show_standard_fields)
    Print[ Magb,  OnElementsOf Domain,  File StrCat[DirResFields, "Magb", ExtGmsh],  LastTimeStepOnly] ;
  EndIf
  //  , StoreInVariable $Magb maybe use this for Core Loss

  // Energy
  //Print[ MagEnergy,  OnElementsOf Domain,  File StrCat[DirResFields, "MagEnergy", ExtGmsh],  LastTimeStepOnly ] ;

  // Current Density
  //Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirResFields, "jz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ j, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirResFields, "j", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ J_rms, OnElementsOf Region[{Domain}], File StrCat[DirResFields, "J_rms", ExtGmsh], LastTimeStepOnly ] ;

  // Ohmic Loss
  If(Flag_show_standard_fields)
    Print[ j2F, OnElementsOf Region[{DomainC}], File StrCat[DirResFields, "j2F", ExtGmsh], LastTimeStepOnly ] ;
  EndIf
  //Print[ j2F_density, OnElementsOf Region[{DomainC}], File StrCat[DirResFields, "j2F_density", ExtGmsh], LastTimeStepOnly ] ;
  If(Flag_show_standard_fields)
    Print[ j2H,   OnElementsOf DomainS, File StrCat[DirResFields,"jH",ExtGmsh] ] ;
  EndIf
  //Print[ j2H_density,   OnElementsOf DomainS, File StrCat[DirResFields,"jH_density",ExtGmsh] ] ;

  //Print[ j2Hprox,   OnElementsOf DomainS, File StrCat[DirResFields,"jHprox",ExtGmsh] ] ;
  //Print[ j2Hskin,   OnElementsOf DomainS, File StrCat[DirResFields,"jHskin",ExtGmsh] ] ;


  // Settings
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "results/option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

}


// === Global/Integrated Quantities ===

PostOperation Get_global UsingPost MagDyn_a {

  // Losses

  // Core
  Print[ j2F[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"CoreEddyCurrentLosses.dat"]] ;

  // Windings Total
  //Print[ SoF[ DomainC ], OnGlobal, Format TimeTable,  File > Sprintf("results/SF_iron.dat")] ; // Complex power
  Print[ j2H[ DomainS ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"LossesStrandedWindings.dat"] ] ;
  Print[ j2F[ Winding1 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2F_1.dat"]] ;
  Print[ j2F[ Winding2 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2F_2.dat"]] ;

  Print[ j2H[ StrandedWinding1 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_1.dat"] ] ;
  //Print[ j2H[ StrandedWinding1 ], OnGlobal, Format Table];
  Print[ j2H[ StrandedWinding2 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_2.dat"] ] ;
  //Print[ j2H[ StrandedWinding2 ], OnGlobal, Format Table];

  // Single Turns



  If(Flag_HomogenisedModel1) // Differentiate fine und hom
    For isF In {1:nbturns1}
      Print[ j2H[ TurnStrand1~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsPrimary,"Losses_turn_%g.dat"], isF] ] ;
    EndFor
  Else
    For isF In {1:nbturns1}
      Print[ j2F[ Turn1~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsPrimary,"Losses_turn_%g.dat"], isF] ] ;
    EndFor
  EndIf

  If(Flag_Transformer)
    If(Flag_HomogenisedModel2) // Differentiate fine und hom
      For isF In {1:nbturns2}
        Print[ j2H[ TurnStrand2~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsSecondary,"Losses_turn_%g.dat"], isF] ] ;
      EndFor
    Else
      For isF In {1:nbturns2}
        Print[ j2F[ Turn2~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsSecondary,"Losses_turn_%g.dat"], isF] ] ;
      EndFor
    EndIf
  EndIf



  //Print[ j2Hskin[StrandedWinding1],   OnGlobal , Format Table];
  //Print[ j2Hprox[StrandedWinding1],   OnGlobal , Format Table];
  //Print[ j2Hskin[StrandedWinding2],   OnGlobal , Format Table];
  //Print[ j2Hprox[StrandedWinding2],   OnGlobal , Format Table];

  //Print[ SoH[ DomainS ], OnGlobal, Format TimeTable, File Sprintf("results/SH_f%g.dat", Freq) ] ;

  If(Flag_Generalized_Steinmetz_loss)
    Print[ piGSE[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"piGSE.dat"]] ;// Core losses
    Print[ piGSE[ Iron ], OnGlobal, Format Table];
  EndIf

  If(Flag_Steinmetz_loss)
    Print[ pSE[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"pSE.dat"]] ;// Core losses
    Print[ pSE[ Iron ], OnGlobal, Format Table];
  EndIf

  Print[ p_hyst[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"p_hyst.dat"]] ;// Core losses



  // Stored Energy
  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
   File > StrCat[DirResVals,"ME.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];
  //Print[ MagEnergy, OnRegion Domain, Format Table, File Sprintf("results/MagEnergy.dat") ];

  // Flux (Linkage)
  Print[ Flux_Linkage_1[DomainCond1], OnGlobal, Format Table, File > StrCat[DirResVals,"Flux_Linkage_1.dat"]];
  Print[ Flux_Linkage_1[DomainCond1], OnGlobal, Format Table];
  If(Flag_Transformer)
    Print[ Flux_Linkage_2[DomainCond2], OnGlobal, Format Table, File > StrCat[DirResVals,"Flux_Linkage_2.dat"]];
    Print[ Flux_Linkage_2[DomainCond2], OnGlobal, Format Table];
  EndIf

  // Inductances
  If(Val_EE_1!=0)
    Print[ L_11[DomainCond1], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_11.dat"]] ;
    Print[ L_11_from_MagEnergy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_11_from_MagEnergy.dat"]] ;
  EndIf
  If(Flag_Transformer)
      If(Val_EE_2!=0)
        Print[ L_22[DomainCond2], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_22.dat"]] ;
        Print[ L_22_from_MagEnergy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_22_from_MagEnergy.dat"]] ;
      EndIf
  EndIf

  // Voltage
  Print[ Voltage_1[DomainCond1], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Voltage_1.dat"]];
  If(Flag_Transformer)
    Print[ Voltage_2[DomainCond2], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Voltage_2.dat"]];
  EndIf

  // Circuit Quantities
  If(!Flag_Circuit)
    Print[ I, OnRegion Winding1, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_f%g.dat"], Freq] , LastTimeStepOnly];
    Print[ U, OnRegion Winding1, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_f%g.dat"], Freq] , LastTimeStepOnly];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("results/I_f%g.dat", Freq), LastTimeStepOnly];
    Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("results/U_f%g.dat", Freq), LastTimeStepOnly];
  EndIf

}

