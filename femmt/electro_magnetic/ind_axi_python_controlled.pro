// ----------------------
// Files and Directories
// ----------------------
Include "Parameter.pro";
Include "postquantities.pro";
If(Flag_Permeability_From_Data)
  Include "core_materials_temp.pro";
EndIf
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
CORE_PN                = 120000;

//physical numbers of conductors in n transformer
For n In {1:n_windings}
       iCOND~{n} = 130000 + 1000*(n-1);
       istrandedCOND~{n} = 140000 + 1000*(n-1);
EndFor



// ----------------------
// Groups
// ----------------------
Group{
  // ----------------------
  // Physical Domains
  // ----------------------
  Air  = Region[{AIR, AIR_EXT}];

  Core = Region[{}];
  // Core Domain
  For n In {1:nCoreParts}
    CorePart~{n} = Region[{(CORE_PN+n-1)}];
    Core += Region[{(CORE_PN+n-1)}];
  EndFor

  // Non Conducting Domain:
  // Initialize the core-shell domain region to air
  DomainCC = Region[{Air}];

  // Add the Core region to the core-shell domain region
  If(!Flag_Conducting_Core)
    DomainCC += Region[{Core}];
  EndIf

  // Boundary Conditions
  // including symmetry
  OuterBoundary = Region[{OUTBND}];

  // Current Conducting Domains
  // Create a region for the winding
  For n In {1:n_windings} // loop over each winding //added by Othman
      Winding~{n} = Region[{}]; // create a region for the winding
      StrandedWinding~{n} = Region[{}]; // create a region for the stranded winding
  EndFor

  // Loop over each winding
  For n In {1:n_windings}
      nbturns~{n} = NbrCond~{n} / SymFactor;
       // Loop over each turn in this winding to create a region for turns and then adding it to the winding
      For winding_number In {1:nbturns~{n}}
        Turn~{n}~{winding_number} = Region[{(iCOND~{n}+winding_number-1)}];
        Winding~{n} += Region[{(iCOND~{n}+winding_number-1)}];
        TurnStrand~{n}~{winding_number} = Region[{(istrandedCOND~{n}+winding_number-1)}];
        StrandedWinding~{n} += Region[{(istrandedCOND~{n}+winding_number-1)}];
      EndFor
  EndFor
   // Add this winding to the core domain region
  For n In {1:n_windings}
      DomainC += Region[{Winding~{n}}] ;
  EndFor
   // Add the Core region to the core domain region
  If(Flag_Conducting_Core)
    DomainC += Region[{Core}];
  EndIf
   // Add this stranded winding to the shell domain region
  For n In {1:n_windings}
      DomainS += Region[{StrandedWinding~{n}}] ;
  EndFor
  // Add the shell domain to the core-shell domain region
  DomainCC          += Region[{DomainS}] ;
   //  the linear and non linear domains to air and all windings

  If(Flag_NL)
    Domain_Lin      = Region[{Air}];
    For n In {1:n_windings}
        Domain_Lin += Region[{Winding~{n}, StrandedWinding~{n}}];
    EndFor
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Core}];
  Else
    Domain_Lin      = Region[{Air, Core}];
    For n In {1:n_windings}
        Domain_Lin += Region[{Winding~{n}, StrandedWinding~{n}}];
    EndFor
    Domain_Lin_NoJs = Region[{Air, Core}];
    Domain_NonLin   = Region[{}];
  EndIf
  // Initialize the main domain to the core and core-shell domains
  Domain = Region[{DomainC, DomainCC}] ;

 // Loop over each winding and add its regions to the corresponding conductor domain
  For n In {1:n_windings} // added by Othman
      DomainCond~{n} += Region[{Winding~{n}, StrandedWinding~{n}}] ;
  EndFor


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

  For n In {1:n_windings}   //added by Othman
      Capacitance_Cir~{n} = Region[{}] ;
  EndFor
  Capacitance_Cir = {};
  For n In {1:n_windings}
      Capacitance_Cir += Region[ {Capacitance_Cir~{n}}] ;
  EndFor

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

  For n In {1:n_windings}
      AreaCell[#{StrandedWinding~{n}}] = AreaCell~{n} ;
  EndFor

  // in non-stranded domains, def. AreaCell to 1 (neutral element of mutiplication)

  AreaCell[#{Air, Core}] = 1.;
  For n In {1:n_windings}
      AreaCell[#{Winding~{n}}] = 1.;
  EndFor


  // Material Properties

  // sigma: conductivity (= imaginary part of complex permitivity)
  //rho[] = 1/sigma[];
  // Set the conductivity of the winding region
  For n In {1:n_windings}
      sigma[#{Winding~{n}}] = sigma_winding~{n} ;
  EndFor

  If(Flag_Conducting_Core)
    sigma[#{Core}] = sigma_core;
    sigma[#{Air}] = 0.;
  EndIf
  If(!Flag_Conducting_Core)
    sigma[#{Air, Core}] = 0.;
  EndIf

  // nu: reluctivity
  // nu = 1/mu
  nu[#{Air}] = Complex[nu0, 0];
  mu[#{Air}] = Complex[mu0, 0];

  For n In {1:n_windings}
      nu[#{Winding~{n}}] = Complex[nu0, 0];
  EndFor

  // Hysteresis Loss
  // Imaginary Part Of Permeability
  // Liste von Lukas hinterlegen
  //mu_imag[ #{Core} ] = mu0 * f_mu_imag[$1, $2];

  If(!Flag_NL)
    If(Flag_Fixed_Loss_Angle)
        mu[#{Core}]   = Complex[mu0*mur_real, -mu0*mur_imag] ;
        nu[#{Core}]   = 1/mu[$1, $2] ;
    ElseIf(Flag_Permeability_From_Data)
        //mu[#{Core}]   = Complex[mu0*(mur^2-f_mu_imag[$1, $2]^2)^(0.5), mu0*f_mu_imag[$1, $2]] ;  // TODO
        mu[#{Core}]   = Complex[mu0*f_mu_real[$1], -mu0*f_mu_imag[$1]] ;
        nu[#{Core}]   = 1/mu[$1, $2] ;
    Else
        mu[#{Core}]   = mu0*mur ;
        nu[#{Core}]   = 1/mu[$1, $2] ;
    EndIf

  Else
    //nu[ #{Core} ] = nu_3kW[$1] ;
    //h[ #{Core} ]  = h_3kW[$1];
    //dhdb_NL[ #{Core} ]= dhdb_3kW_NL[$1] ;
    //dhdb[ #{Core} ]   = dhdb_3kW[$1] ;

    //nu[ #{Core} ] = nu_N95[$1] ;
    nu[ #{Core} ] = nu~{Core_Material}[$1] ;
    h[ #{Core} ]  = h~{Core_Material}[$1];
    //dhdb_NL[ #{Core} ]= dhdb_95_NL[$1] ;
    dhdb_NL[ #{Core} ]= dhdb_95100_NL[$1] ;
    dhdb[ #{Core} ]   = dhdb~{Core_Material}[$1] ;
  EndIf

  // Excitation Current


  For n In {1:n_windings}
      FSinusoidal~{n}[] = F_Cos_wt_p[]{2*Pi*Freq, Phase~{n}}; //Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
      Fct_Src~{n}[] = FSinusoidal~{n}[];
      Signn~{n} = (Phase~{n}==Pi) ? -1 : 1;
  EndFor

  // Auxiliary functions for post-processing
  nuOm[#{Air}] = nu[]*Complex[0.,1.];
  nuOm[#{Core}] = -nu[$1]*Complex[0.,1.];
  //nuOm[#{Winding1, Winding2, Winding3}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ];


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
         nuOm[#{StrandedWinding~{n}}] = Complex[ 2 * Pi * Freq * Im[nu[]], -Re[nu[]] ];
         // sigma[#{StrandedWinding~{n}}] = SymFactor * Complex[ skin_rhor~{n}[Rr~{n}] / sigma_winding~{n} / Fill~{n},  2*Pi*Freq*skin_rhoi~{n}[Rr~{n}]*mu0/(8*Pi*Fill~{n})];
         sigma[#{StrandedWinding~{n}}] = SymFactor * Complex[ skin_rhor~{n}[Rr~{n}] / sigma_winding~{n} / Fill~{n}, 0];
      EndIf
  EndFor

  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  // List of nodes related to circuit

  For n In {1:n_windings}
      N~{n}~{1}() = {1:nbturns~{n}};
      N~{n}~{2}() = {2:nbturns~{n}+1};
  EndFor

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

    }
  }


  // ---- Only For Voltage Excitation ----
  { Name Voltage_2D ;
    Case{
      If(Flag_Conducting_Core)
        { Region Core ; Value 0; }
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

      CreateDir[DirResValsCore];
      For n In {1:n_windings}
          CreateDir[DirResValsWinding~{n}];
      EndFor

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
      { Name mur ; Value { Term { [ 1 / Norm [Im[ nu[{d a}, Freq]] * mu0 ] ] ; In Core ; Jacobian Vol ; } } }
      { Name mur_norm ; Value { Term { [ Norm [Im[ mu[{d a}, Freq]] / mu0 ] ] ; In Core ; Jacobian Vol ; } } }
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
             [ CoefGeo/ElementVol[]*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2F ; Value { Integral {
             [ 0.5*CoefGeo*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;// 0.5 for frequency domain
             In DomainC ; Jacobian Vol ; Integration II ; } } }

           { Name j2F_density ; Value { Integral {
             [ 0.5*CoefGeo/ElementVol[]*sigma[]*SquNorm[ {ur}/CoefGeo - Dt[{a}] ] ] ;// 0.5 for frequency domain
             In DomainC ; Jacobian Vol ; Integration II ; } } }
      EndIf


      // DomainS (Stranded Conductors)

      If(Freq==0.0)
           { Name j2H ; Value { Integral {
             [ CoefGeo*( Re[-{d a}*Conj[nuOm[]*{d a}]] + sigma[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
             In DomainS ; Jacobian Vol ; Integration II ; } } }
      Else
           { Name j2H ; Value { Integral {
             [ 0.5*CoefGeo*( Norm[ Re[{d a}*Conj[nuOm[]*{d a}]] ] + sigma[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 for frequency domain
             In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2H_density ; Value { Integral {
             [ 0.5*CoefGeo/ElementVol[]*( Norm[ Re[{d a}*Conj[nuOm[]*{d a}]] ] + sigma[]*SquNorm[-1/AreaCell[]*{ir}]) ] ; // 0.5 for frequency domain
             In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2Hprox ; Value { Integral {
            [ 0.5*CoefGeo*Norm[ Re[{d a}*Conj[nuOm[]*{d a}]] ] ] ;// 0.5 for frequency domain
            In DomainS ; Jacobian Vol ; Integration II ; } } }

           { Name j2Hskin ; Value { Integral {
            [ 0.5*CoefGeo*sigma[]*SquNorm[-1/AreaCell[]*{ir}] ] ;// 0.5 for frequency domain
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
                                   ) ] ; In Core ; Jacobian Vol ; Integration II ;} } }
      EndIf

      If(Flag_Steinmetz_loss)
        { Name pSE ; Value { Integral { [ CoefGeo * ki * Freq^alpha * (Norm[{d a}])^beta
                                     ] ; In Core ; Jacobian Vol ; Integration II ;} } }

        { Name pSE_density ; Value { Integral { [ CoefGeo* ki * Freq^alpha * (Norm[{d a}])^beta
                                     ] ; In Core ; Jacobian Vol ; Integration II ;} } }
      EndIf



      // ------------------------------------------------------------------------------------------------
      // Hysteresis Losses (According To Complex Core Parameters)

      { Name p_hyst ; Value { Integral {
        // [ 0.5 * CoefGeo * 2*Pi*Freq * Im[mu[Norm[{d a}], Freq]] * SquNorm[nu[Norm[{d a}], Freq] * Norm[{d a}]] ] ;
        [ - 0.5 * CoefGeo * 2*Pi*Freq * Im[mu[{d a}, Freq]] * SquNorm[nu[{d a}, Freq] * {d a}] ] ;
        In Core ; Jacobian Vol ; Integration II ;} } }          // TODO: mur 2350 | general mur; multiplication at simulation begin with loss angle

      { Name p_hyst_density ; Value { Integral {
        [ - 0.5 * CoefGeo/ElementVol[] * 2*Pi*Freq * Im[mu[{d a}, Freq]] * SquNorm[nu[Norm[{d a}], Freq] * {d a}] ] ;
        In Core ; Jacobian Vol ; Integration II ;} } }



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
            [ CoefGeo * ({d a}*Conj[nuOm[{d a}]*{d a}] + sigma[]*SquNorm[-1/AreaCell[]*{ir}]) ] ;
            In DomainC ; Jacobian Vol ; Integration II ; } } } //Complex power
            // xfmr changed Domain to DomainC
            // to prevent from div by zero error in "Air" and "Core" domains



      // ------------------------------------------------------------------------------------------------
      // Voltage (Voltage_i = dFlux_Linkage_i / dt)
      // Distinguish between litz wire case and solid case for n-windings

      For n In {1:n_windings}
          If(Flag_HomogenisedModel~{n})
             { Name Voltage~{n} ; Value { Integral { [ CoefGeo / AreaCell~{n} * (CompZ[Dt[{a}]] + sigma[]*CompZ[{ir}] / AreaCell~{n}) ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
          Else
             { Name Voltage~{n} ; Value { Integral { [ CompZ[{ur}] / AreaCell~{n} ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
          EndIf
      EndFor

      // ------------------------------------------------------------------------------------------------
      // Flux (Linkage)

      For n In {1:n_windings}
          { Name Flux_Linkage~{n} ; Value {
            Integral { [ CoefGeo / AreaCell~{n} * CompZ[{a}] ]; In DomainCond~{n}; Jacobian Vol; Integration II; } } }
      EndFor

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





