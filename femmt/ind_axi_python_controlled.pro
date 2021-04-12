// Time  domain + fine model:
// gmsh  ind_axi.geo -setnumber Flag_HomogenisedModel 0 -2 -o fine.msh -v 3
// getdp ind_axi -msh fine.msh -setnumber Flag_HomogenisedModel 0 -sn Flag_FD 0 -sn Flag_SrcType 2 -sn NbT 1 -sol Analysis

// Resistance (Ohms): p 'SF_s0.dat' u 1:2 w l lw 2, 'SH.dat' u 1:2 w p ps 2 pt 6 lw 2
// Inductance (mH):   p 'SF_s0.dat' u 1:($3*1e3) w l lw 2, 'SH.dat' u 1:($3*1e3) w p ps 2 pt 6 lw 2

// ====================================================
// Added constants by Till
Fill = 0.9;
mfem = "{1Analysis parameters/";
col1 = "AliceBlue";
col2 = "Blue";
col3 = "Ivory";

// coordinates of rectangular winding window
Xw1 = 11e-3 ;
Xw2 = Xw1 + 8.e-3 ;
Yw1 = -25.5e-3/2 ;
Yw2 =  25.5e-3/2 ;

SigmaCu = 6e7;
mu0 = 4.e-7 * Pi ;
nu0 = 1./mu0;
Rc = Sqrt[1/Pi]*1e-3; // needs review Till: Conductor radius
NbrCond = 6;
NbrTurns = 1;
Len = (2*Pi*(Xw1+Xw2)/2)*NbrTurns ;
AreaCond = Pi*Rc^2;
Dex = 2.2*Rc;
Dey = Dex;
//AreaCell = Dex*Dey;
AreaCell = 0.008*0.028;
Rdc = Len/SigmaCu/AreaCond;
Flag_HalfModel = 0;
Flag_HomogenisedModel = 0;
Zh~{1} = 100;


//----------------------------------
// Physical numbers
//----------------------------------
OUTBND  = 1111;

AIR    = 1000;
AIRGAP = 1100;
IRON   = 2000;
//Till INSULATION = 3000;

iCOND = 4000;

ALLCOND = 5000;
// ====================================================

Include "ind_axi_dat.pro"
Include "BH.pro"

la = Fill; // fill factor

// new files, finer mesh, max X=8
file_ZSkinRe  = Sprintf("coeff/pI_RS_la%.2g.dat", la);
file_ZSkinIm  = Sprintf("coeff/qI_RS_la%.2g.dat", la);
file_NuProxRe = Sprintf("coeff/qB_RS_la%.2g.dat", la);
file_NuProxIm = Sprintf("coeff/pB_RS_la%.2g.dat", la);

// Till removed time domain

DirRes = "res/";
po = "{Output/";

DefineConstant[
  visu = {0, Choices{0, 1}, AutoCheck 0,
    Name StrCat[mfem,"Visu/Real-time visualization"], Highlight "LightPink"}

  Flag_FD = { 1, Choices{0,1}, Name StrCat[mfem,"2Frequency domain analysis?"], Highlight "AliceBlue" }

// Till removed time domain

  Flag_IronCore    = {1, Choices{0,1}, Name StrCat[mfem, "3Core/Iron core?"], Highlight Str[col1]}

  Flag_NL  = {0, Choices{0,1}, Name StrCat[mfem, "3Core/Nonlinear bh-curve?"],
              Highlight Str[col1], Visible Flag_IronCore, ReadOnly Flag_FD}

  Nb_max_iter = {100, Name StrCat[mfem, "3Core/Nonlinear solver/Max. num. iterations"],
    Visible Flag_NL, Highlight Str[col3]}
  iter_max = Nb_max_iter
  relaxation_factor = {1., Name StrCat[mfem, "3Core/Nonlinear solver/Relaxation factor"],
    Visible Flag_NL, Highlight Str[col3]}
  stop_criterion = {1e-8, Name StrCat[mfem, "3Core/Nonlinear solver/Stopping criterion"],
    Visible Flag_NL, Highlight Str[col3]} //1e-4

  Flag_ImposedVoltage = {0, Choices{0,1}, Name StrCat[mfem,"Source/001Imposed Voltage?"]}
  Flag_Circuit = {Flag_ImposedVoltage, Choices{0,1}, Name StrCat[mfem,"Source/002Use circuit"],
    ReadOnly (Flag_ImposedVoltage==1)}

  Flag_imposedRr   = {0, Choices{0,1},
    Name StrCat[mfem, "Source/1Imposed reduced frequency X"], Highlight Str[col1]}
  Flag_imposedFreq = !Flag_imposedRr

// Till  removed time domain

  ExtGmsh = ".pos"
  ExtGnuplot = ".dat"
];

If(Flag_imposedRr) // Reduced frequency
  DefineConstant[
    Rr = {4, Min 0.1, Max 5, Step 0.1, Name StrCat[mfem,"Source/1Reduced frequency"], ReadOnly 0, Highlight "Ivory"}
    delta = {Rc/Rr, Name StrCat[mfem,"Source/2Skin depth [m]"], ReadOnly 1, Highlight "LightGrey"}
    Freq  = {1/(delta*delta*mu0*SigmaCu*Pi), Name StrCat[mfem,"Source/3Frequency [Hz]"], ReadOnly 1, Highlight "LightGrey"}
  ];
Else
  DefineConstant[
    Freq = { 200000, Min 0.1, Max 500e3, Name StrCat[mfem,"Source/3Frequency [Hz]"], ReadOnly 0, Highlight "Ivory"}
    delta = {1/Sqrt[mu0*SigmaCu*Freq*Pi], Name StrCat[mfem,"Source/2Skin depth [m]"], ReadOnly 1, Highlight "LightGrey"}
    Rr = {Rc/delta, Name StrCat[mfem,"Source/1Reduced frequency"], ReadOnly 1, Highlight "LightGrey"}
  ];
EndIf


DefineConstant[
  Omega  = 2*Pi*Freq,
  Period = 1./Freq,

  Vdc = 50 // amplitude of PWM voltage
  Val_EE = { Vdc, Name StrCat[mfem,"Source/003source amplitude"],  Highlight "Ivory", Visible 1}

  //------------------------------------------------------

// Till removed time domain
];


Group{
  Air  = Region[{AIR, AIRGAP}];
  Insulation = Region[{INSULATION}];

  If(Flag_IronCore)
    Iron = Region[{IRON}];
  Else
    Iron = Region[{}];
    Air  += Region[{IRON}];
  EndIf

  OuterBoundary = Region[{OUTBND}]; // including symmetry

  Winding =  Region[{}] ;

  DomainCC = Region[{Air, Insulation, Iron}] ;

  SymFactor = Flag_HalfModel ? 2.:1. ; //half inductor with axisymmetry

  nbturns = (Flag_HomogenisedModel==0) ? NbrCond/SymFactor : 1  ; // number of turns

  If (Flag_HomogenisedModel==0) // Fine case
    For iF In {1:nbturns}
      Turn~{iF} = Region[{(iCOND+iF-1)}] ;
      Winding  += Region[{(iCOND+iF-1)}] ;
    EndFor

    DomainC = Region[{Winding}] ;
    DomainS = Region[{}] ;
  EndIf

  If (Flag_HomogenisedModel==1) //Homogenised case
    Turn~{1} = Region[{iCOND}] ;
    Winding = Region[{iCOND}];
    DomainC = Region[{}] ;
    DomainS = Region[{Winding}] ;
  EndIf

  //DomainCC += Region[{DomainS}] ;

  If(Flag_NL)
    Domain_Lin      = Region[{Air, Insulation, Winding}];
    Domain_Lin_NoJs = Region[{Air, Insulation}];
    Domain_NonLin   = Region[{Iron}];
  Else
    Domain_Lin      = Region[{Air, Insulation, Winding, Iron}];
    Domain_Lin_NoJs = Region[{Air, Insulation, Iron}];
    Domain_NonLin   = Region[{}];
  EndIf

  Domain = Region[{DomainC, DomainCC}] ;

  //--------------------------------------------------------
  //--------------------------------------------------------

  // Groups related to source circuit
  Input = # 12345 ;

   // Groups related to the circuit
  Input = # 12345 ;
  iZH    = 10000 ;
  iLH    = 20000 ; // Zskin in homog. coil
  iLHp   = 30000 ;

// Till removed time domain

  Resistance_Cir  = Region[{}];
  If(Flag_FD && Flag_HomogenisedModel==1) // Frequency domain
    Resistance_Cir += Region[{Zh~{1}}];
  EndIf

// Till removed time domain

  Inductance_Cir  = Region[{}];

// Till removed time domain

  Capacitance1_Cir = Region[ {} ] ;
  Capacitance2_Cir = Region[ {} ] ;
  Capacitance_Cir = Region[ {Capacitance1_Cir, Capacitance2_Cir} ] ;

  SourceV_Cir = Region[ {Input} ] ;
  SourceI_Cir = Region[ {} ] ;

  DomainZ_Cir = Region[ {Resistance_Cir, Inductance_Cir, Capacitance_Cir} ] ;

  DomainSource_Cir = Region[ {SourceV_Cir, SourceI_Cir} ] ;
  DomainZt_Cir = Region[ {DomainZ_Cir, DomainSource_Cir} ] ;

}


Function {

  CoefGeo = 2*Pi*SymFactor ; // axisymmetry + symmetry factor

  sigma[#{Winding}] = SigmaCu ;
  sigma[#{Air, Insulation, Iron}] = 0.;
  rho[] = 1/sigma[];

  nu[#{Air, Insulation}] = nu0;

  If (!Flag_NL)
    nu[#{Iron}]   = nu0/1000;
  Else
    nu[ #{Iron} ] = nu_3kW[$1] ;
    h[ #{Iron} ]  = h_3kW[$1];
    dhdb_NL[ #{Iron} ]= dhdb_3kW_NL[$1] ;
    dhdb[ #{Iron} ]   = dhdb_3kW[$1] ;
  EndIf


  //==================================================================
  If(Flag_FD) // Frequency domain
    FSinusoidal[] = Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Else // Time domain
    FSinusoidal[] = Complex_MH[0,-1]{Freq} ; //Sin F_Sin_wt_p[]{2*Pi*Freq, 0};
  EndIf
  //------------------------------------------------------

  Fct_Src[] = FSinusoidal[];

  //==================================================================

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


  If(Flag_HomogenisedModel==0)
    nu[Winding] = nu0 ;
  Else
    //Proximity effect
    nu[Winding] = nu0*Complex[prox_nur[Rr], prox_nui[Rr]*Fill*Rr^2/2];
  EndIf
  If(Flag_FD) // used if Flag_HomogenisedModel==1
    // Skin effect => complex impedance
    Zskin[] = 1/SymFactor*Complex[ skin_rhor[Rr]*Rdc, 2*Pi*Freq*skin_rhoi[Rr]*mu0*Len/(8*Pi*Fill)] ;
  EndIf
  //==================================================================
  // Auxiliary functions for post-processing
  nuOm[#{Air, Insulation}] = -nu[]*Complex[0.,1.];
  nuOm[#{Iron}] = -nu[$1]*Complex[0.,1.];
  nuOm[#{Winding}] = Complex[ Omega * Im[nu[]], -Re[nu[]] ];

  kkk[] =  SymFactor * skin_rhor[Rr] / Fill /SigmaCu ;
  //==================================================================

// Till removed time domain

  //==================================================================

// Till removed time domain

  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  Ns[] = (Flag_HomogenisedModel==0) ? 1 : NbrCond/SymFactor ;
  Sc[] =  SurfaceArea[]{iCOND};

  If (Flag_HomogenisedModel==1)
    // Accounting for eddy currents (homogenization)
    If(Flag_FD)
      Resistance[Zh~{1}] = Zskin[] ;
    EndIf

// Till removed time domain
  EndIf


  // List of nodes related to circuit
  N1() = {1:nbturns};   // Node 1 for each turn
  N2() = {2:nbturns+1}; // Node 2 for each turn

}


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
      If(Flag_HomogenisedModel==0)
        { Region Input ;  Branch {N1(0), N2(nbturns-1)} ; }
      EndIf
      If(Flag_HomogenisedModel==1 && Flag_FD)
          { Region Input  ; Branch {777, N2(nbturns-1)} ; }
          { Region Zh~{1} ; Branch {777, N1(0)}; } // Complex impedance: Zskin
      EndIf

      // Till removed time domain

      For k In {0:nbturns-1} // list indexes start at 0
        { Region Turn~{k+1} ; Branch {N1(k), N2(k)} ; }
      EndFor
    }
  }

}

//-----------------------------------------------------------------------------

Jacobian {
  { Name Vol ; Case { { Region All ; Jacobian VolAxiSqu ; } } }
  { Name Sur ; Case { { Region All ; Jacobian SurAxi ; } } }
}

Integration {
  { Name II ; Case {
      { Type Gauss ; Case {
          { GeoElement Triangle ;    NumberOfPoints 4 ; }
          { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
        } }
    } }
}

//-----------------------------------------------------------------------------

FunctionSpace {

  { Name Hcurl_a_2D ; Type Form1P ; // Split for TD + homog: subspace isolating unknowns
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ All, Not Winding ] ; }
      { Name seh ; NameOfCoef aeh ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ Winding ] ; }
    }
    SubSpace {
      { Name aH ; NameOfBasisFunction {seh}; } // Subspace only used in TD homog
    }

    Constraint {
      { NameOfCoef ae  ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
     }
  }

  { Name Hregion_u_2D ; Type Form1P ; // Gradient of Electric scalar potential (2D)
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC ; Entity DomainC ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef I ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS ; Entity DomainS ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


  // For circuit equations
  { Name Hregion_Z ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainZt_Cir ; Entity DomainZt_Cir ; }
    }
    GlobalQuantity {
      { Name Iz ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Uz ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }

// Till removed time domain

}


Formulation {
  { Name MagDyn_a ; Type FemEquation ;
    Quantity {
       { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

       { Name ur ; Type Local  ; NameOfSpace Hregion_u_2D  ; }
       { Name I  ; Type Global ; NameOfSpace Hregion_u_2D[I] ; }
       { Name U  ; Type Global ; NameOfSpace Hregion_u_2D[U] ; }

       { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
       { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
       { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

       { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
       { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Domain_Lin ; Jacobian Vol ; Integration II ; }
      If(Flag_NL)
        Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
        Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{ur}/CoefGeo , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{ur}/CoefGeo , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{I}, {U} ] ; In DomainC ; }

      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ CoefGeo*Ns[]/Sc[] * Dof{a}, {ir} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }

      Galerkin { [ Ns[]/Sc[] / sigma[] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
        In DomainS ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{Us}/CoefGeo , {Is} ] ; In DomainS ; }

      If(Flag_Circuit)
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }

        GlobalTerm { [ Dof{Uz}                      , {Iz} ] ; In Inductance_Cir ; }
        GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ] ; In Inductance_Cir ; }

        GlobalTerm { [ Dof{Iz}        , {Iz} ] ; In Capacitance1_Cir ; }
        GlobalTerm { NeverDt[ Dof{Iz} , {Iz} ] ; In Capacitance2_Cir ; }
        GlobalTerm { DtDof [ Capacitance[] * Dof{Uz} , {Iz} ] ; In Capacitance_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {I};  Loop {U};  Equation {I};  In DomainC ; }
          { Node {Is}; Loop {Us}; Equation {Us}; In DomainS ; }
          { Node {Iz}; Loop {Uz}; Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }

  { Name MagDyn_a_Homog ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

      { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
      { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {

      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Domain_Lin ; Jacobian Vol ; Integration II ; }

      If(Flag_NL)
        Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
        Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      Galerkin { [ -1/AreaCell * Dof{ir}, {a} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ 1/AreaCell * Dof{a}, {ir} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{Us}/CoefGeo, {Is} ]     ; In DomainS ; }

      // Circuit equations
      If(Flag_Circuit)
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }

        GlobalTerm { [ Dof{Uz}                      , {Iz} ] ; In Inductance_Cir ; }
        GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ] ; In Inductance_Cir ; }

        GlobalTerm { [ Dof{Iz}        , {Iz} ] ; In Capacitance1_Cir ; }
        GlobalTerm { NeverDt[ Dof{Iz} , {Iz} ] ; In Capacitance2_Cir ; }
        GlobalTerm { DtDof [ Capacitance[] * Dof{Uz} , {Iz} ] ; In Capacitance_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {Is};  Loop {Us};  Equation {Us}; In DomainS ; }
          { Node {Iz};  Loop {Uz};  Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }

}


Resolution {

  { Name Analysis ;
    System {
      If(Flag_HomogenisedModel==0)
        If(Flag_FD) // Frequency domain
          { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue ; Frequency Freq ; }
// Till removed time domain
        EndIf
      EndIf
      If(Flag_HomogenisedModel==1)
        If(Flag_FD) // Frequency domain
          { Name A ; NameOfFormulation MagDyn_a_Homog ; Type ComplexValue ; Frequency Freq ; }
// Till removed time domain
        EndIf
      EndIf
    }
    Operation {
      CreateDir[DirRes];
      InitSolution[A]; SaveSolution[A];

      If(Flag_FD) // Frequency domain
        Generate[A] ; Solve[A] ; SaveSolution[A] ;
        If(Flag_HomogenisedModel==0)
          PostOperation [Map_local];
          PostOperation [Get_global];
        Else
          PostOperation [Map_local_Homog];
          PostOperation [Get_global_Homog];
        EndIf

// Till removed time domain
      EndIf
    }
  }


}// Resolution

PostProcessing {

  { Name MagDyn_a ; NameOfFormulation MagDyn_a ; NameOfSystem A;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name h ; Value { Term { [ nu[{d a}]*{d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name j ; Value { Term { [ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ; In DomainC ; Jacobian Vol ; } } }
      { Name jz ; Value {
          Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
          Term { [ CompZ[ {ir}/AreaCond ] ] ; In DomainS ; Jacobian Vol ; }
        } }

      { Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell ] ;
            In Domain ; Jacobian Vol  ; Integration II ; } } }


      { Name j2F ; Value { Integral { // Joue losses
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
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2 ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }
     }
  }


  { Name MagDyn_a_Homog ; NameOfFormulation MagDyn_a_Homog ;
    PostQuantity {
      { Name a   ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az  ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }

      { Name j  ; Value { Term { [ -1/AreaCell*{ir} ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz ; Value { Term { [ CompZ[ -1/AreaCell*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz_homo ; Value { Term { [  AreaCell*{ir}  ] ; In DomainS ; Jacobian Vol ; } } }


      { Name j2H ; Value { Integral { // Joule losses
            [ CoefGeo*(Re[{d a}*Conj[nuOm[]*{d a}]]+kkk[]*SquNorm[-1/AreaCell*{ir}]) ] ;
            In DomainS ; Jacobian Vol ; Integration II ; } } }

      { Name SoH ; Value { Integral { // Complex power = Active power +j * Reactive power => S = P+j*Q
            [ CoefGeo * ({d a}*Conj[nuOm[{d a}]*{d a}] + kkk[]*SquNorm[-1/AreaCell*{ir}]) ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } } //Complex power

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

 // Till removed time domain

}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

PostOperation Map_local UsingPost MagDyn_a {
  Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "jz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ J_rms, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "J_rms", ExtGmsh], LastTimeStepOnly ] ;
  Print[ j2F, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "j2F", ExtGmsh], LastTimeStepOnly ] ;
  Print[ b,  OnElementsOf Domain,  File StrCat[DirRes, "b", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ MagEnergy,  OnElementsOf Domain,  File StrCat[DirRes, "MagEnergy", ExtGmsh],  LastTimeStepOnly ] ;
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
}

PostOperation Get_global UsingPost MagDyn_a {
  Print[ j2F[ Winding ], OnGlobal, Format TimeTable, File > Sprintf("res/j2F_iron%g.dat", Flag_IronCore)] ;// Joule losses
  Print[ SoF[ Domain ], OnGlobal, Format TimeTable,  File > Sprintf("res/SF_iron%g.dat", Flag_IronCore)] ; // Complex power
}


PostOperation Get_allTS UsingPost MagDyn_a {
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



//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

PostOperation Map_local_Homog UsingPost MagDyn_a_Homog {
  //Print[ jz,   OnElementsOf DomainS, File StrCat[DirRes,"jH",ExtGmsh] ] ;
  Print[ j2H,   OnElementsOf DomainS, File StrCat[DirRes,"jH",ExtGmsh] ] ;
  Print[ b,   OnElementsOf Domain, File StrCat[DirRes,"bH",ExtGmsh] ] ;
  Print[ raz, OnElementsOf Domain, File StrCat[DirRes,"aH",ExtGmsh] ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
}

PostOperation Get_global_Homog UsingPost MagDyn_a_Homog {
  // Complex power: S = P + i*Q, P=active power, Q=reactive power

  Print[ SoH[Domain], OnGlobal, Format TimeTable,
    File Sprintf("res/SH_f%g.dat", Freq) ] ;
  Print[ j2H[Winding], OnGlobal, Format Table,
    File Sprintf("res/j2H_f%g.dat", Freq) ] ;

  If(!Flag_Circuit)
    Print[ U, OnRegion Winding, Format Table, File Sprintf("res/U_f%g.dat", Freq) ] ;
    Print[ I, OnRegion Winding, Format Table, File Sprintf("res/I_f%g.dat", Freq) ] ;
  Else
    Print[ U, OnRegion Input, Format Table, File Sprintf("res/U_f%g.dat", Freq) ];
    Print[ I, OnRegion Input, Format Table, File Sprintf("res/I_f%g.dat", Freq) ];
  EndIf

}

// Till removed time domain



// This is only for Onelab
DefineConstant[
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 1, Closed 1}
  C_ = {"-solve -v 4 -v2 -bin", Name "GetDP/9ComputeCommand", Visible 1}
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 1}
];
