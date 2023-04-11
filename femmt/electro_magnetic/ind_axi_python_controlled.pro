// ----------------------
// Files and Directories
// ----------------------
Include "Parameter.pro";
Include "postquantities.pro";
Include "core_materials_temp.pro";
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
// 1 means full cylinder
SymFactor               = 1. ;
CoefGeo                 = 2*Pi*SymFactor ; // axisymmetry +/* symmetry factor */
n_windings = Number_of_Windings;  //added by Othman

// ----------------------
// Physical numbers
// ----------------------
OUTBND              = 111111;
AIR                 = 110000;
AIR_EXT             = 110001;
IRON                = 120000;
iCOND~{1}           = 130000;
istrandedCOND~{1}   = 140000;

If(Flag_Transformer)
   For n In {2:n_windings}
       iCOND~{n} = iCOND~{1} + 1000*(n-1);
       istrandedCOND~{n} = istrandedCOND~{1} + 1000*(n-1);
   EndFor
EndIf


/*
If(Flag_Transformer)
  iCOND2            = 131000;
  istrandedCOND2    = 141000;
EndIf
If(Flag_Three_Transformer)
  iCOND3            = 132000;
  istrandedCOND3    = 142000;
EndIf
*/




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
  /*
  Winding1 =  Region[{}] ;
  Winding2 =  Region[{}] ;
  Winding3 =  Region[{}] ;
  StrandedWinding1 =  Region[{}] ;
  StrandedWinding2 =  Region[{}] ;
  StrandedWinding3 =  Region[{}] ;

  Winding~{1} =  Region[{}] ;
  /*
  For n In {1:n_windings} // loop over each winding //added by Othman
      Winding~{n} = Region[{}]; // create a region for the winding
      StrandedWinding~{n} = Region[{}]; // create a region for the stranded winding
  EndFor
  */

  If(Flag_Transformer)
    For n In {2:n_windings} // loop over each winding //added by Othman
        Winding~{n} = Region[{}]; // create a region for the winding
        StrandedWinding~{n} = Region[{}]; // create a region for the stranded winding
    EndFor
  EndIf


  /*
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
  */
  /*
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
  */
  /*
  // Tertiary(Transformer)
  If(Flag_Three_Transformer)
    nbturns3 = NbrCond3/SymFactor;
    For iF In {1:nbturns3}
        Turn3~{iF} = Region[{(iCOND3+iF-1)}] ;
        Winding3  += Region[{(iCOND3+iF-1)}] ;
    EndFor
    For isF In {1:nbturns3}
        TurnStrand3~{isF} = Region[{(istrandedCOND3+isF-1)}] ;
        StrandedWinding3  += Region[{(istrandedCOND3+isF-1)}] ;
    EndFor
  EndIf
  */
  /*
  For n In {1:n_windings} //added by Othman
      nbturns~{n} = NbrCond~{n}/SymFactor;
      For iF In {1:nbturns~{n}}
          Turn~{n}~{iF} = Region[{(iCOND~{n}+iF-1)}] ;
          Winding~{n}  += Region[{(iCOND~{n}+iF-1)}] ;
      EndFor
      For isF In {1:nbturns~{n}}
          TurnStrand~{n}~{isF} = Region[{(istrandedCOND~{n}+isF-1)}] ;
          StrandedWinding~{n}~{isF} += Region[{(istrandedCOND~{n}+isF-1)}] ;
      EndFor
  EndFor
  */

  //Inductor
  nbturns~{1} = NbrCond~{1}/SymFactor;
  For iF In {1:nbturns~{1}}
      Turn~{1}~{iF} = Region[{(iCOND~{1}+iF-1)}] ;
      Winding~{1}  += Region[{(iCOND~{1}+iF-1)}] ;
  EndFor
  For isF In {1:nbturns~{1}}
      TurnStrand~{1}~{isF} = Region[{(istrandedCOND~{1}+isF-1)}] ;
      StrandedWinding~{1}  += Region[{(istrandedCOND~{1}+isF-1)}] ;
  EndFor
  //Transformer
  If (Flag_Transformer)
    For n In {2:n_windings} //added by Othman
        nbturns~{n} = NbrCond~{n}/SymFactor;
        For iF In {1:nbturns~{n}}
            Turn~{n}~{iF} = Region[{(iCOND~{n}+iF-1)}] ;
            Winding~{n}  += Region[{(iCOND~{n}+iF-1)}] ;
        EndFor
        For isF In {1:nbturns~{n}}
            TurnStrand~{n}~{isF} = Region[{(istrandedCOND~{n}+isF-1)}] ;
            StrandedWinding~{n}~{isF} += Region[{(istrandedCOND~{n}+isF-1)}] ;
        EndFor
    EndFor
  EndIf

  /*
  For n In {1:n_windings}  //added by Othman
      DomainC = Region[{Winding~{n}}] ;
  EndFor
  */
  DomainC  = Region[{Winding~{1}}];
  If(Flag_Transformer)
     For n In {2:n_windings}  //added by Othman
        DomainC += Region[{Winding~{n}}] ;
     EndFor
  EndIf




  //DomainC           = Region[{Winding1, Winding2, Winding3}] ;
  If(Flag_Conducting_Core)
    DomainC         += Region[{Iron}] ;
  EndIf
  //DomainS           = Region[{StrandedWinding1, StrandedWinding2, StrandedWinding3}] ;
  /*
  For n In {1:n_windings} //added by Othman
      DomainS = Region[{StrandedWinding~{n}}] ;
  EndFor
  */


  DomainS = Region[{StrandedWinding~{1}}];
  If(Flag_Transformer)
    For n In {2:n_windings} //added by Othman
        DomainS += Region[{StrandedWinding~{n}}] ;
    EndFor
  EndIf

  DomainCC          += Region[{DomainS}] ;
  /*
  If(Flag_NL)
    For n In {1:n_windings}
        Domain_Lin = Region[{Air, Winding~{n}, StrandedWinding~{n}}];
    EndFor
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Iron}];
  Else
    For n In {1:n_windings}
        Domain_Lin = Region[{Air, Iron,Winding~{n}, StrandedWinding~{n}}];
    EndFor
    Domain_Lin_NoJs = Region[{Air, Iron}];
    Domain_NonLin   = Region[{}];
  EndIf
  */



  If(Flag_NL)
    Domain_Lin      = Region[{Air, Winding~{1}, StrandedWinding~{1}}];
    If (Flag_Transformer)
        For n In {2:n_windings}
            Domain_Lin += Region[{Winding~{n}, StrandedWinding~{n}}];
        EndFor
    EndIf
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Iron}];
  Else
    Domain_Lin      = Region[{Air, Iron, Winding~{1}, StrandedWinding~{1}}];
    If (Flag_Transformer)
        For n In {2:n_windings}
            Domain_Lin += Region[{Winding~{n}, StrandedWinding~{n}}];
        EndFor
    EndIf
    Domain_Lin_NoJs = Region[{Air, Iron}];
    Domain_NonLin   = Region[{}];
  EndIf


  Domain = Region[{DomainC, DomainCC}] ;
  //DomainCond1 = Region[{Winding1, StrandedWinding1}];
  //DomainCond2 = Region[{Winding2, StrandedWinding2}];
  //DomainCond3 = Region[{Winding3, StrandedWinding3}];
  /*
  For n In {1:n_windings} // added by Othman
      DomainCond~{n} = Region[{Winding~{n}, StrandedWinding~{n}}] ;
  EndFor
  */


  DomainCond~{1} = Region[{Winding~{1}, StrandedWinding~{1}}] ;
  If(Flag_Transformer)
    For n In {2:n_windings} // added by Othman
        DomainCond~{n} = Region[{Winding~{n}, StrandedWinding~{n}}] ;
    EndFor
  EndIf


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

  Inductance_Cir  = Region[ {} ];
  //Capacitance1_Cir = Region[ {} ] ;
  //Capacitance2_Cir = Region[ {} ] ;
  //Capacitance3_Cir = Region[ {} ] ;
  Capacitance_Cir~{1} = Region[{}] ;
  /*
  For n In {1:n_windings}   //added by Othman
        Capacitance_Cir~{n} = Region[{}] ;
  EndFor

  For n In {1:n_windings}
        Capacitance_Cir = Region[ {Capacitance_Cir~{n}}] ;
  EndFor
  */


  If (Flag_Transformer)
     For n In {2:n_windings}   //added by Othman
        Capacitance_Cir~{n} = Region[{}] ;
     EndFor
  EndIf
  Capacitance_Cir = Region[ {Capacitance_Cir~{1}}] ;
  If(Flag_Transformer)
    For n In {2:n_windings}
        Capacitance_Cir += Region[ {Capacitance_Cir~{n}}] ;
    EndFor
  EndIf


  //Capacitance_Cir = Region[ {Capacitance_Cir~{1}, Capacitance_Cir~{n}}] ;





  //Capacitance_Cir = Region[ {Capacitance1_Cir, Capacitance2_Cir, Capacitance3_Cir} ] ;

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
  /*
  For n In {1:n_windings}
      AreaCell[#{StrandedWinding~{n}}] = AreaCell~{n} ;
  EndFor
  */

  AreaCell[#{StrandedWinding~{1}}] = AreaCell~{1};
  If(Flag_Transformer)
    For n In {2:n_windings}
        AreaCell[#{StrandedWinding~{n}}] = AreaCell~{n} ;
    EndFor
  EndIf

  /*
  If(Flag_Transformer)
    AreaCell[#{StrandedWinding2}] = AreaCell2;
  EndIf
  If(Flag_Three_Transformer)
    AreaCell[#{StrandedWinding3}] = AreaCell3;
  EndIf
  */
  // in non-stranded domains, def. AreaCell to 1 (neutral element of mutiplication)
  /*
  For n In {1:n_windings}
      AreaCell[#{Winding~{n}}] = 1.;
  EndFor
  */


  AreaCell[#{Air, Iron, Winding~{1}}] = 1.;
  If(Flag_Transformer)
    For n In {2:n_windings}
        AreaCell[#{Winding~{n}}] = 1.;
    EndFor
  EndIf



  // Material Properties

  // sigma: conductivity (= imaginary part of complex permitivity)
  //rho[] = 1/sigma[];
  /*
  For n In {1:n_windings}
      sigma[#{Winding~{n}}] = sigma_winding~{n} ;
  EndFor
  */


  sigma[#{Winding~{1}}] = sigma_winding~{1} ;
  If (Flag_Transformer)
    For n In {2:n_windings}
        sigma[#{Winding~{n}}] = sigma_winding~{n} ;
    EndFor
  EndIf






  /*
  If(Flag_Three_Transformer)
    sigma[#{Winding3}] = sigma_winding_3 ;
  EndIf
  */
  If(Flag_Conducting_Core)
    sigma[#{Iron}] = sigma_core;
    sigma[#{Air}] = 0.;
  EndIf
  If(!Flag_Conducting_Core)
    sigma[#{Air, Iron}] = 0.;
  EndIf

  // nu: reluctivity
  // nu = 1/mu
  nu[#{Air}] = Complex[nu0, 0];
  mu[#{Air}] = Complex[mu0, 0];
  //nu[#{Winding1, Winding2, Winding3}] = Complex[nu0, 0];
  nu[#{Winding~{1}}] = Complex[nu0, 0];
  /*
  For n In {1:n_windings}
      nu[#{Winding~{n}}] = Complex[nu0, 0];
  EndFor
  */

  If(Flag_Transformer)
    For n In {2:n_windings}
        nu[#{Winding~{n}}] = Complex[nu0, 0];
    EndFor
  EndIf



  /*
  If(Flag_Transformer)
    For n In {n_windings} //added by Othman
        nu[#{Winding~{n}] = Complex[nu0, 0];
    EndFor
  EndIf
  */
  // Hysteresis Loss
  // Imaginary Part Of Permeability
  // Liste von Lukas hinterlegen
  //mu_imag[ #{Iron} ] = mu0 * f_mu_imag[$1, $2];

  If(!Flag_NL)
    If(Flag_Fixed_Loss_Angle)
        mu[#{Iron}]   = Complex[mu0*mur_real, -mu0*mur_imag] ;
        nu[#{Iron}]   = 1/mu[$1, $2] ;
    ElseIf(Flag_Permeability_From_Data)
        //mu[#{Iron}]   = Complex[mu0*(mur^2-f_mu_imag[$1, $2]^2)^(0.5), mu0*f_mu_imag[$1, $2]] ;  // TODO
        mu[#{Iron}]   = Complex[mu0*f_mu_real[$1], -mu0*f_mu_imag[$1]] ;
        nu[#{Iron}]   = 1/mu[$1, $2] ;
    Else
        mu[#{Iron}]   = mu0*mur ;
        nu[#{Iron}]   = 1/mu[$1, $2] ;
    EndIf

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

  /*
  For n In {1:n_windings}
      FSinusoidal~{n}[] = F_Cos_wt_p[]{2*Pi*Freq, Phase~{n}}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
      Fct_Src~{n}[] = FSinusoidal~{n}[];
      Signn~{n} = (Phase~{n}==Pi) ? -1 : 1;  //TODO: Inductance Calc
  EndFor
  */


  FSinusoidal~{1}[] = F_Cos_wt_p[]{2*Pi*Freq, Phase~{1}}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Fct_Src~{1}[] = FSinusoidal~{1}[];
  Signn~{1} = (Phase~{1}==Pi) ? -1 : 1;  //TODO: Inductance Calc
  If(Flag_Transformer)
    For n In {2:n_windings}
       FSinusoidal~{n}[] = F_Cos_wt_p[]{2*Pi*Freq, Phase~{n}}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
       Fct_Src~{n}[] = FSinusoidal~{n}[];
       Signn~{n} = (Phase~{n}==Pi) ? -1 : 1;  //TODO: Inductance Calc
    EndFor
  EndIf

  /*
  FSinusoidal1[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_1}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Fct_Src1[] = FSinusoidal1[];
  Sign1 = (Phase_1==Pi) ? -1 : 1;  //TODO: Inductance Calc

  If(Flag_Transformer)
    FSinusoidal2[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_2}; //Complex_MH[1, 0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
    Fct_Src2[] = FSinusoidal2[];
    Sign2 = (Phase_2==Pi) ? -1 : 1;  //TODO: Inductance Calc
  EndIf

  If(Flag_Three_Transformer)
    FSinusoidal3[] = F_Cos_wt_p[]{2*Pi*Freq, Phase_3}; //Complex_MH[1, 0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
    Fct_Src3[] = FSinusoidal3[];
    Sign3 = (Phase_3==Pi) ? -1 : 1;  //TODO: Inductance Calc
  EndIf
  */

  // Auxiliary functions for post-processing
  nuOm[#{Air}] = nu[]*Complex[0.,1.];
  nuOm[#{Iron}] = -nu[$1]*Complex[0.,1.];
  //nuOm[#{Winding1, Winding2, Winding3}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ];


  // Resistive/Skin Coefficient - will be multiplied with current "ir", which is zero except in windings
  // kkk[#{Iron, Air, Winding1, Winding2, Winding3}] =  0 ; // choose arbitrary value
  //kkk[#{Iron, Air}] =  0 ;
  /*
  For n In {1:n_windings}
      kkk[#{Winding~{n}}] =  0 ;
  EndFor
  */


  kkk[#{Iron, Air, Winding~{1}}] =  0 ; // choose arbitrary value
  If (Flag_Transformer)
    For n In {2:n_windings}
        kkk[#{Winding~{n}}] =  0 ;
    EndFor
  EndIf

  /*
  For n In {1:n_windings}
      If(Flag_HomogenisedModel~{n})
         // Secondary
         file_ZSkinRe~{n}  = Sprintf("Strands_Coefficients/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
         file_ZSkinIm~{n}  = Sprintf("Strands_Coefficients/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
         file_NuProxRe~{n} = Sprintf("Strands_Coefficients/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
         file_NuProxIm~{n} = Sprintf("Strands_Coefficients/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
         skin_rhor_list~{n}() = ListFromFile[ file_ZSkinRe~{n} ];
         skin_rhoi_list~{n}() = ListFromFile[ file_ZSkinIm~{n} ];
         prox_nur_list~{n}()  = ListFromFile[ file_NuProxRe~{n} ];
         prox_nui_list~{n}()  = ListFromFile[ file_NuProxIm~{n} ];
         skin_rhor~{n}[] = InterpolationLinear[$1]{ skin_rhor_list~{n}() };
         skin_rhoi~{n}[] = InterpolationLinear[$1]{ skin_rhoi_list~{n}() };
         prox_nur~{n}[]  = InterpolationLinear[$1]{ prox_nur_list~{n}() } ;
         prox_nui~{n}[]  = InterpolationLinear[$1]{ prox_nui_list~{n}() } ;
         // Formula from Paper:
         nu[#{StrandedWinding~{n}}] = nu0*Complex[prox_nur~{n}[Rr~{n}], prox_nui~{n}[Rr~{n}]*Fill~{n}*Rr~{n}^2/2];
         nuOm[#{StrandedWinding~{n}}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill
         kkk[#{StrandedWinding~{n}}] =  SymFactor * skin_rhor~{n}[Rr~{n}] / sigma_winding~{n} / Fill~{n} ;
         sigma[#{StrandedWinding~{n}}] = SymFactor * skin_rhor~{n}[Rr~{n}] / sigma_winding~{n} / Fill~{n} ;
      EndIf
  EndFor
  */



  If(Flag_HomogenisedModel~{1})
    // Homogenization coefficients
    // Primary (Inductor + Transformer)
    file_ZSkinRe~{1}  = Sprintf("Strands_Coefficients/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill~{1}, NbrLayers~{1});
    file_ZSkinIm~{1}  = Sprintf("Strands_Coefficients/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill~{1}, NbrLayers~{1});
    file_NuProxRe~{1}= Sprintf("Strands_Coefficients/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill~{1}, NbrLayers~{1});
    file_NuProxIm~{1} = Sprintf("Strands_Coefficients/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill~{1}, NbrLayers~{1});
    skin_rhor_list~{1}() = ListFromFile[ file_ZSkinRe~{1} ];
    skin_rhoi_list~{1}() = ListFromFile[ file_ZSkinIm~{1} ];
    prox_nur_list~{1}()  = ListFromFile[ file_NuProxRe~{1} ];
    prox_nui_list~{1}()  = ListFromFile[ file_NuProxIm~{1} ];
    skin_rhor~{1}[] = InterpolationLinear[$1]{ skin_rhor_list~{1}() };
    skin_rhoi~{1}[] = InterpolationLinear[$1]{ skin_rhoi_list~{1}() };
    prox_nur~{1}[]  = InterpolationLinear[$1]{ prox_nur_list~{1}() } ;
    prox_nui~{1}[]  = InterpolationLinear[$1]{ prox_nui_list~{1}() } ;
    nu[#{StrandedWinding~{1}}] = nu0*Complex[prox_nur~{1}[Rr~{1}], prox_nui~{1}[Rr~{1}]*Fill~{1}*Rr~{1}^2/2];
    nuOm[#{StrandedWinding~{1}}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill
    kkk[#{StrandedWinding~{1}}] =  SymFactor * skin_rhor~{1}[Rr~{1}] / sigma_winding~{1} / Fill~{1} ;  // TODO: Add qI (skin_rhoi_1[]); effects the reactive power (usually has minor effect)
    sigma[#{StrandedWinding~{1}}] = SymFactor * skin_rhor~{1}[Rr~{1}] / sigma_winding~{1} / Fill~{1} ;
  EndIf

  If(Flag_Transformer)
    For n In {2:n_windings}
        If(Flag_HomogenisedModel~{n})
            // Secondary
            file_ZSkinRe~{n}  = Sprintf("Strands_Coefficients/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
            file_ZSkinIm~{n}  = Sprintf("Strands_Coefficients/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
            file_NuProxRe~{n} = Sprintf("Strands_Coefficients/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
            file_NuProxIm~{n} = Sprintf("Strands_Coefficients/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill~{n}, NbrLayers~{n});
            skin_rhor_list~{n}() = ListFromFile[ file_ZSkinRe~{n} ];
            skin_rhoi_list~{n}() = ListFromFile[ file_ZSkinIm~{n} ];
            prox_nur_list~{n}()  = ListFromFile[ file_NuProxRe~{n} ];
            prox_nui_list~{n}()  = ListFromFile[ file_NuProxIm~{n} ];
            skin_rhor~{n}[] = InterpolationLinear[$1]{ skin_rhor_list~{n}() };
            skin_rhoi~{n}[] = InterpolationLinear[$1]{ skin_rhoi_list~{n}() };
            prox_nur~{n}[]  = InterpolationLinear[$1]{ prox_nur_list~{n}() } ;
            prox_nui~{n}[]  = InterpolationLinear[$1]{ prox_nui_list~{n}() } ;
            // Formula from Paper:
            nu[#{StrandedWinding~{n}}] = nu0*Complex[prox_nur~{n}[Rr~{n}], prox_nui~{n}[Rr~{n}]*Fill~{n}*Rr~{n}^2/2];
            nuOm[#{StrandedWinding~{n}}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill
            kkk[#{StrandedWinding~{n}}] =  SymFactor * skin_rhor~{n}[Rr~{n}] / sigma_winding~{n} / Fill~{n} ;
            sigma[#{StrandedWinding~{n}}] = SymFactor * skin_rhor~{n}[Rr~{n}] / sigma_winding~{n} / Fill~{n} ;
        EndIf
    EndFor
  EndIf



  /*
  If(Flag_Transformer)
    If(Flag_HomogenisedModel2)
      // Secondary
      file_ZSkinRe_2  = Sprintf("Strands_Coefficients/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_ZSkinIm_2  = Sprintf("Strands_Coefficients/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_NuProxRe_2= Sprintf("Strands_Coefficients/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
      file_NuProxIm_2 = Sprintf("Strands_Coefficients/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill2, NbrLayers2);
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
      kkk[#{StrandedWinding2}] =  SymFactor * skin_rhor_2[Rr2] / sigma_winding_2 / Fill2 ;
      sigma[#{StrandedWinding2}] = SymFactor * skin_rhor_2[Rr2] / sigma_winding_2 / Fill2 ;
    EndIf
  EndIf
  If(Flag_Three_Transformer)
    If(Flag_HomogenisedModel3)
      // Tertiary
      file_ZSkinRe_3  = Sprintf("Strands_Coefficients/coeff/pI_RS_la%.2g_%.2glayer.dat", Fill3, NbrLayers3);
      file_ZSkinIm_3  = Sprintf("Strands_Coefficients/coeff/qI_RS_la%.2g_%.2glayer.dat", Fill3, NbrLayers3);
      file_NuProxRe_3= Sprintf("Strands_Coefficients/coeff/qB_RS_la%.2g_%.2glayer.dat", Fill3, NbrLayers3);
      file_NuProxIm_3 = Sprintf("Strands_Coefficients/coeff/pB_RS_la%.2g_%.2glayer.dat", Fill3, NbrLayers3);
      skin_rhor_list_3() = ListFromFile[ file_ZSkinRe_3 ];
      skin_rhoi_list_3() = ListFromFile[ file_ZSkinIm_3 ];
      prox_nur_list_3()  = ListFromFile[ file_NuProxRe_3 ];
      prox_nui_list_3()  = ListFromFile[ file_NuProxIm_3 ];
      skin_rhor_3[] = InterpolationLinear[$1]{ skin_rhor_list_3() };
      skin_rhoi_3[] = InterpolationLinear[$1]{ skin_rhoi_list_3() };
      prox_nur_3[]  = InterpolationLinear[$1]{ prox_nur_list_3() } ;
      prox_nui_3[]  = InterpolationLinear[$1]{ prox_nui_list_3() } ;
      // Formula from Paper:
      nu[#{StrandedWinding3}] = nu0*Complex[prox_nur_3[Rr3], prox_nui_3[Rr3]*Fill3*Rr3^2/2];
      nuOm[#{StrandedWinding3}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ]; // sTill
      kkk[#{StrandedWinding3}] =  SymFactor * skin_rhor_3[Rr3] / sigma_winding_3 / Fill3 ;
      sigma[#{StrandedWinding3}] = SymFactor * skin_rhor_3[Rr3] / sigma_winding_3 / Fill3 ;
    EndIf
  EndIf
  */
  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  // List of nodes related to circuit
  // Inductor
  // Primary

  N~{1}~{1}() = {1:nbturns~{1}};   // Node 1 for each turn
  N~{1}~{2}() = {2:nbturns~{1}+1}; // Node 2 for each turn

  If(Flag_Transformer)
    For n In {2:n_windings}
        N~{n}~{1}() = {1:nbturns~{n}};
        N~{n}~{2}() = {2:nbturns~{n}+1};
    EndFor
  EndIf

  /*
  For n In {1:n_windings}
      N~{n}~{1}() = {1:nbturns~{n}};
      N~{n}~{2}() = {2:nbturns~{n}+1};
  EndFor
  */


  /*
  // Transformer
  If(Flag_Transformer)
    // Secondary
    N2_1() = {1:nbturns2};   // Node 1 for each turn
    N2_2() = {2:nbturns2+1}; // Node 2 for each turn
  EndIf

  If(Flag_Three_Transformer)
    // Tertiary
    N3_1() = {1:nbturns3};   // Node 1 for each turn
    N3_2() = {2:nbturns3+1}; // Node 2 for each turn
  EndIf
  */

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
          If(Flag_Circuit==0 && Flag_HomogenisedModel~{1}==0)
            { Region Winding~{1} ; Value Val_EE~{1}; TimeFunction Fct_Src~{1}[] ; }
          EndIf
          If(Flag_Circuit==0 && Flag_HomogenisedModel~{1}==1)
            { Region StrandedWinding~{1} ; Value Val_EE~{1}; TimeFunction Fct_Src~{1}[] ; }
          EndIf
      EndIf
      If(Flag_Transformer)
         For n In {2:n_windings}
             If(1)
                  If(Flag_Circuit==0 && Flag_HomogenisedModel~{n}==0)
                     { Region Winding~{n} ; Value Val_EE~{n}; TimeFunction Fct_Src~{n}[] ; }
                  EndIf
                  If(Flag_Circuit==0 && Flag_HomogenisedModel~{n}==1)
                     { Region StrandedWinding~{n} ; Value Val_EE~{n}; TimeFunction Fct_Src~{n}[] ; }
                  EndIf
             EndIf
         EndFor
      EndIf
      /*
      For n In {1:n_windings}
          If(1)
             If(Flag_Circuit==0 && Flag_HomogenisedModel~{n}==0)
               { Region Winding~{n} ; Value Val_EE~{n}; TimeFunction Fct_Src~{n}[] ; }
             EndIf
             If(Flag_Circuit==0 && Flag_HomogenisedModel~{n}==1)
               { Region StrandedWinding~{n} ; Value Val_EE~{n}; TimeFunction Fct_Src~{n}[] ; }
             EndIf
          EndIf
      EndFor
      */

      /*
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

      If(Flag_Three_Transformer)
        //If(Val_EE_3!=0)
        If(1)
            If(Flag_Circuit==0 && Flag_HomogenisedModel3==0)
              { Region Winding3 ; Value Val_EE_3; TimeFunction Fct_Src3[] ; }
            EndIf
            If(Flag_Circuit==0 && Flag_HomogenisedModel3==1)
              { Region StrandedWinding3 ; Value Val_EE_3; TimeFunction Fct_Src3[] ; }
            EndIf
        EndIf
      EndIf
      */
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
      /*
      CreateDir[DirResValsWinding~{1}];
      For n In {2:n_windings}
          CreateDir[DirResValsWinding~{n}];
      EndFor
      */

      CreateDir[DirResValsWinding~{1}];
      If(Flag_Transformer)
        For n In {2:n_windings}
            CreateDir[DirResValsWinding~{n}];
        EndFor
      EndIf


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
// definition of varoius post-processing quantities with the help of the mag.vec.pot. a
// ----------------------
PostProcessing {

  { Name MagDyn_a ; NameOfFormulation MagDyn_a ; NameOfSystem A;
    PostQuantity {



      // ------------------------------------------------------------------------------------------------
      // Potentials

      { Name a ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az_int ; Value { Integral { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ; Integration II ; } } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name ur ; Value { Term { [ {ur} ] ; In Domain  ; Jacobian Vol  ;} } }



      // ------------------------------------------------------------------------------------------------
      // Electrical Field

      { Name e ; Value { Term { [ -1*(Dt[{a}]+{ur}/CoefGeo) ] ; In Domain ; Jacobian Vol ; } } }
      { Name MagEz ; Value { Term { [  -1*(Dt[{a}]+{ur}/CoefGeo)  ] ; In Domain ; Jacobian Vol ; } } }



      // ------------------------------------------------------------------------------------------------
      // Magnetic Field

      { Name h ; Value { Term { [ nu[{d a}, Freq]*{d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name Magh ; Value { Term { [ Norm[ nu[{d a}, Freq]*{d a} ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name Mag_h_real ; Value { Term { [ Norm[ Re [ nu[{d a}, Freq]*{d a} ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name Mag_h_imag ; Value { Term { [ Norm[ Im [ nu[{d a}, Freq]*{d a} ] ] ] ; In Domain ; Jacobian Vol ; } } }



      // ------------------------------------------------------------------------------------------------
      // Magnetic Flux Density

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name b_pol ; Value { Term { [ Norm [ Re [ Cart2Pol[ {d a} ] ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name im_b_pol ; Value { Term { [ Norm [ Im [ Cart2Pol[ {d a} ] ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name Magb ; Value { Term { [ Norm[ {d a} ] ]; In Domain ; Jacobian Vol ; } } }
      { Name Mag_b_real ; Value { Term { [ Norm[ Re [ {d a} ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name Mag_b_imag ; Value { Term { [ Norm[ Im [ {d a} ] ] ] ; In Domain ; Jacobian Vol ; } } }



      // ------------------------------------------------------------------------------------------------
      // Permeability Plot

      { Name nur ; Value { Term { [ Norm[ nu[{d a}, Freq] / mu0 ] ] ; In Domain ; Jacobian Vol ; } } }
      //{ Name mur ; Value { Term { [ 1 / Norm[ nu[{d a}, Freq] / mu0 ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name mur ; Value { Term { [ 1 / Norm [Im[ nu[{d a}, Freq]] * mu0 ] ] ; In Iron ; Jacobian Vol ; } } }
      { Name mur_norm ; Value { Term { [ Norm [Im[ mu[{d a}, Freq]] / mu0 ] ] ; In Iron ; Jacobian Vol ; } } }
      { Name mur_re ; Value { Term { [ Re[ 1/nu[{d a}, Freq] / mu0 ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name mur_im ; Value { Term { [ Norm [ Im[ 1/nu[{d a}, Freq] / mu0 ] ] ] ; In Domain ; Jacobian Vol ; } } }
      { Name nur_re ; Value { Term { [ Re[ nu[{d a}, Freq] * mu0 ] ] ; In Domain ; Jacobian Vol ; } } }  // := mur_re / (mur_re^2 + mur_im^2)
      { Name nur_im ; Value { Term { [ Im[ nu[{d a}, Freq] * mu0 ] ] ; In Domain ; Jacobian Vol ; } } }  // := mur_im / (mur_re^2 + mur_im^2)


      // ------------------------------------------------------------------------------------------------
      // Current Density

      { Name j ; Value {
            Term { [ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ; In DomainC ; Jacobian Vol ; }
            Term { [ -1/AreaCell[]*{ir} ] ; In DomainS ; Jacobian Vol ; } } }

      { Name ir ; Value {
            Term { [ {ir} ] ; In DomainS ; Jacobian Vol ; } } }

      { Name ir_re ; Value {
            Term { [ Re[{ir}] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name ir_im ; Value {
            Term { [ Im[{ir}] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name ir_norm ; Value {
            Term { [ Norm[{ir}] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name jz ; Value {
            Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
            Term { [ CompZ[ -1/AreaCell[]*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name J_rms ; Value {
            Term { [ Norm[ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ] ; In DomainC ; Jacobian Vol ; }
            Term { [ Norm[ CompZ[ {ir}/AreaCell[] ] ] ] ; In DomainS ; Jacobian Vol ; }
            Term { [ Norm[ 0 ] ] ; In Domain_Lin_NoJs ; Jacobian Vol ; } } } //non_lin case is missing



      // ------------------------------------------------------------------------------------------------
      // Power Loss


      // DomainC (Solid Conductors)

      If(Freq==0.0)
          { Name j2F ; Value { Integral {
             [ CoefGeo*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;
             In DomainC ; Jacobian Vol ; Integration II ; } } }
          { Name j2F_density ; Value { Integral {
             [ CoefGeo/ElementVol[]*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2F ; Value { Integral {
             [ 0.5*CoefGeo*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }

           { Name j2F_density ; Value { Integral {
             [ 0.5*CoefGeo/ElementVol[]*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;  // 0.5* added by Till
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      EndIf


      // DomainS (Stranded Conductors)

      If(Freq==0.0)
           { Name j2H ; Value { Integral {
             [ CoefGeo*( Re[-{d a}*Conj[nuOm[]*{d a}]] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
             In DomainS ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2H ; Value { Integral {
             [ 0.5*CoefGeo*( Norm[ Re[{d a}*Conj[nuOm[]*{d a}]] ] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 added
             In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2H_density ; Value { Integral {
             [ 0.5*CoefGeo/ElementVol[]*( Norm[ Re[{d a}*Conj[nuOm[]*{d a}]] ] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 added
             In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2Hprox ; Value { Integral {
            [ 0.5*CoefGeo*Norm[ Re[{d a}*Conj[nuOm[]*{d a}]] ] ] ;// 0.5 added by Till
            In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2Hskin ; Value { Integral {
            [ 0.5*CoefGeo*kkk[]*SquNorm[-1/AreaCell[]*{ir}] ] ;// 0.5 added by Till
            In DomainS ; Jacobian Vol ; Integration II ; } } }
      EndIf



      // ------------------------------------------------------------------------------------------------
      // Steinmetz Core Loss

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



      // ------------------------------------------------------------------------------------------------
      // Hysteresis Losses (According To Complex Core Parameters)

      { Name p_hyst ; Value { Integral {
        // [ 0.5 * CoefGeo * 2*Pi*Freq * Im[mu[Norm[{d a}], Freq]] * SquNorm[nu[Norm[{d a}], Freq] * Norm[{d a}]] ] ;
        [ - 0.5 * CoefGeo * 2*Pi*Freq * Im[mu[{d a}, Freq]] * SquNorm[nu[{d a}, Freq] * {d a}] ] ;
        In Iron ; Jacobian Vol ; Integration II ;} } }          // TODO: mur 2350 | general mur; multiplication at simulation begin with loss angle

      { Name p_hyst_density ; Value { Integral {
        [ - 0.5 * CoefGeo/ElementVol[] * 2*Pi*Freq * Im[mu[{d a}, Freq]] * SquNorm[nu[Norm[{d a}], Freq] * {d a}] ] ;
        In Iron ; Jacobian Vol ; Integration II ;} } }



      // ------------------------------------------------------------------------------------------------
      // Energy

      { Name MagEnergy ; Value {
          Integral { [ 1/4*CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }

      { Name ElectEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }

      //{ Name b2F ; Value { Integral { // Magnetic Energy
      //      [ CoefGeo*nu[{d a}]*SquNorm[{d a}] ] ;
      //     In Domain ; Jacobian Vol ; Integration II ; } } }
      //{ Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell[] ] ;
      //      In Domain ; Jacobian Vol  ; Integration II ; } } }



      // ------------------------------------------------------------------------------------------------
      // Complex Power


      // DomainC (Solid Conductors)

      { Name SoF ; Value {
          Integral {
            [ CoefGeo * Complex[ sigma[] * SquNorm[(Dt[{a}]+{ur}/CoefGeo)], nu[{d a}]*SquNorm[{d a}] ] ] ; // added {d a} in nu[...]
            In DomainC ; Jacobian Vol ; Integration II ; }
        } }
        // together: Domain -> DomainC


      // DomainS (Stranded Conductors)

      { Name SoH ; Value { Integral { // Complex power = Active power +j * Reactive power => S = P+j*Q
            [ CoefGeo * ({d a}*Conj[nuOm[{d a}]*{d a}] + kkk[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
            In DomainC ; Jacobian Vol ; Integration II ; } } } //Complex power
            // xfmr changed Domain to DomainC
            // to prevent from div by zero error in "Air" and "Core" domains



      // ------------------------------------------------------------------------------------------------
      // Voltage (Voltage_i = dFlux_Linkage_i / dt)
      // Distinguish between litz wire case and solid case

      For n In {1:n_windings}
          If(Flag_HomogenisedModel~{n})
             { Name Voltage~{n} ; Value { Integral { [ CoefGeo / AreaCell~{n} * (CompZ[Dt[{a}]] + kkk[]*CompZ[{ir}] / AreaCell~{n}) ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
          Else
             { Name Voltage~{n} ; Value { Integral { [ CompZ[{ur}] / AreaCell~{n} ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
          EndIf
      EndFor
      /*
      If(Flag_HomogenisedModel~{1})
        { Name Voltage~{1} ; Value { Integral { [ CoefGeo / AreaCell~{1} * (CompZ[Dt[{a}]] + kkk[]*CompZ[{ir}] / AreaCell~{1}) ]; In DomainCond~{1}; Jacobian Vol; Integration II; } } }
      Else
        { Name Voltage~{1} ; Value { Integral { [ CompZ[{ur}] / AreaCell~{1} ]; In DomainCond~{1}; Jacobian Vol; Integration II; } } }
      EndIf
      If(Flag_Transformer)
         For n In {2:n_windings}
             If(Flag_HomogenisedModel~{n})
                { Name Voltage~{n} ; Value { Integral { [ CoefGeo / AreaCell~{n} * (CompZ[Dt[{a}]] + kkk[]*CompZ[{ir}] / AreaCell~{n}) ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
             Else
                { Name Voltage~{n} ; Value { Integral { [ CompZ[{ur}] / AreaCell~{n} ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
             EndIf
         EndFor
      EndIf
      */


      /*
      If(Flag_Transformer)
        If(Flag_HomogenisedModel2)
          { Name Voltage_2 ; Value { Integral { [ CoefGeo / AreaCell2 * (CompZ[Dt[{a}]] + kkk[]*CompZ[{ir}] / AreaCell2) ]; In DomainCond2; Jacobian Vol; Integration II; } } }
        Else
          { Name Voltage_2 ; Value { Integral { [ CompZ[{ur}] / AreaCell2 ]; In DomainCond2; Jacobian Vol; Integration II; } } }
        EndIf
      EndIf

      If(Flag_Three_Transformer)
        If(Flag_HomogenisedModel3)
          { Name Voltage_3 ; Value { Integral { [ CoefGeo / AreaCell3 * (CompZ[Dt[{a}]] + kkk[]*CompZ[{ir}] / AreaCell3) ]; In DomainCond3; Jacobian Vol; Integration II; } } }
        Else
          { Name Voltage_3 ; Value { Integral { [ CompZ[{ur}] / AreaCell3 ]; In DomainCond3; Jacobian Vol; Integration II; } } }
        EndIf
      EndIf
      */



      // ------------------------------------------------------------------------------------------------
      // Flux (Linkage)

      For n In {1:n_windings}
          { Name Flux_Linkage~{n} ; Value {
            Integral { [ CoefGeo / AreaCell~{n} * CompZ[{a}] ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
      EndFor

      /*
      { Name Flux_Linkage~{1} ; Value {
        Integral { [ CoefGeo / AreaCell~{1} * CompZ[{a}] ]; In DomainCond~{1}; Jacobian Vol; Integration II; } } }

      If(Flag_Transformer)
        For n In {2:n_windings}
            { Name Flux_Linkage~{n} ; Value {
              Integral { [ CoefGeo / AreaCell~{n} * CompZ[{a}] ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
        EndFor
      EndIf
      */

      /*
      If(Flag_Transformer)
        { Name Flux_Linkage_2 ; Value {
          Integral { [ CoefGeo / AreaCell2 * CompZ[{a}] ]; In DomainCond2; Jacobian Vol; Integration II; } } }
      EndIf

      If(Flag_Three_Transformer)
        { Name Flux_Linkage_3 ; Value {
          Integral { [ CoefGeo / AreaCell3 * CompZ[{a}] ]; In DomainCond3; Jacobian Vol; Integration II; } } }
      EndIf
      */


      // ------------------------------------------------------------------------------------------------
      // (Self) Inductances

      For n In {1:n_windings}
          If(Val_EE~{n}!=0)
            { Name L~{n}~{n} ; Value { Integral {
              [ Signn~{n} * CoefGeo / AreaCell~{n} * CompZ[{a}] / Val_EE~{n} ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
            { Name LFromMagEnergy~{n}~{n} ; Value { Integral {
              [ 2 * CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) / (Val_EE~{n}*Val_EE~{n}) ]; In Domain; Jacobian Vol; Integration II; } } }
          EndIf
      EndFor


      /*
      If(Val_EE~{1}!=0)
        { Name L_1_1 ; Value { Integral {
          [ Signn~{1} * CoefGeo / AreaCell~{1} * CompZ[{a}] / Val_EE~{1} ]; In DomainCond~{1}; Jacobian Vol; Integration II; } } }
        { Name LFromMagEnergy_1_1 ; Value { Integral {
          [ 2 * CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) / (Val_EE~{1}*Val_EE~{1}) ]; In Domain; Jacobian Vol; Integration II; } } }
      EndIf

      If(Flag_Transformer)
        For n In {2:n_windings}
           If(Val_EE~{n}!=0)
             { Name L~{n}~{n} ; Value { Integral {
               [ Signn~{n} * CoefGeo / AreaCell~{n} * CompZ[{a}] / Val_EE~{n} ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
             { Name LFromMagEnergy~{n}~{n} ; Value { Integral {
               [ 2 * CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) / (Val_EE~{n}*Val_EE~{n}) ]; In Domain; Jacobian Vol; Integration II; } } }
           EndIf
        EndFor
      EndIf
      */


      /*
      If(Flag_Transformer)
        For n In {2:n_windings}
            For j In {2:n_windings}
                If (n == j)
                    If(Val_EE~{n}!=0)
                      { Name L~{n}~{j} ; Value { Integral {
                        [ Signn~{n} * CoefGeo / AreaCell~{n} * CompZ[{a}] / Val_EE~{n} ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
                      { Name L~{n}~{j}_from_MagEnergy ; Value { Integral {
                        [ 2 * CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) / (Val_EE~{n}*Val_EE~{n}) ]; In Domain; Jacobian Vol; Integration II; } } }
                    EndIf
                EndIf
            EndFor
        EndFor
      EndIf
      */
      /*
      If(Flag_Transformer)
        If(Val_EE_2!=0)
          { Name L_22 ; Value { Integral {
            [ Sign2 * CoefGeo / AreaCell2 * CompZ[{a}] / Val_EE_2 ]; In DomainCond2; Jacobian Vol; Integration II; } } }
          { Name L_22_from_MagEnergy ; Value { Integral {
            [ 2 * CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) / (Val_EE_2*Val_EE_2) ]; In Domain; Jacobian Vol; Integration II; } } }
        EndIf
      EndIf

      If(Flag_Three_Transformer)
        If(Val_EE_3!=0)
          { Name L_33 ; Value { Integral {
            [ Sign3 * CoefGeo / AreaCell3 * CompZ[{a}] / Val_EE_3 ]; In DomainCond3; Jacobian Vol; Integration II; } } }
          { Name L_33_from_MagEnergy ; Value { Integral {
            [ 2 * CoefGeo*nu[{d a}, Freq]*({d a}*{d a}) / (Val_EE_3*Val_EE_3) ]; In Domain; Jacobian Vol; Integration II; } } }
        EndIf
      EndIf
      */


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


Include "fields.pro";
Include "values.pro";





