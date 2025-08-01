// ----------------------
Include "Parameter.pro";
Include "postquantities.pro";
// All Variables - remove or create in python
// ----------------------
Nb_max_iter             = 20;
relaxation_factor       = 1.;
stop_criterion          = 1e-8;

// ----------------------
// half inductor with axisymmetry
// 1 means full cylinder
SymFactor               = 1. ;
CoefGeo                 = 2*Pi*SymFactor ; // axisymmetry +/* symmetry factor */
n_windings              = Number_of_Windings; // Number of windings

// ----------------------
// Physical numbers
// ----------------------
OUTBND                  = 111111;
BOUND_LEFT              = 111112;
AIR                     = 110000;
AIR_EXT                 = 110001;
AIR_COND                = 1000000;
Bobbin_Insulation       = 110010;
Layer_Insulation        = 110050;
//Cond_Insulation         = 2000000;
CORE_PN                 = 120000;
ExtGmsh                 = ".pos";
po = "Output 2D/";


//physical numbers of conductors in n transformer
For n In {1:n_windings}
       iCOND~{n} = 130000 + 1000*(n-1);
       istrandedCOND~{n} = 150000 + 1000*(n-1);
EndFor

//physical numbers of conductor insulation in n winding
For n In {1:n_windings}
       Cond_Insulation~{n} = 2000000 + 1000*(n-1);
EndFor

// ----------------------
// Groups (Physical Domains Setup for Electrostatics)
// ----------------------
Group {

    Air  = Region[{AIR, AIR_EXT}];
    // Insulation = Region[{Insulation}];

    // Core
    Core = Region[{}];
    // Core Domain
    For n In {1:nCoreParts}
      CorePart~{n} = Region[{(CORE_PN+n-1)}];
      Core += Region[{(CORE_PN+n-1)}];
    EndFor
    // cond_insulation domain
    /*ConductorInsulation = Region[{}];
    For n In {1:n_windings}
      CondIsoPart~{n} = Region[{(Cond_Insulation+1000*(n-1))}];
      ConductorInsulation += Region[{CondIsoPart~{n}}];
    EndFor*/


    // Winding Domains Setup
    For n In {1:n_windings}  // Define regions for each winding
        Winding~{n} = Region[{}];  // Initialize a region for the winding
    EndFor
    // Loop to Define Regions for Each Turn within Windings
    For winding_number In {1:n_windings}
        nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

        // Loop over each turn in this winding to create a region for turns and add it to the winding
        For turn_number In {1:nbturns~{winding_number}}
            Turn~{winding_number}~{turn_number} = Region[{(iCOND~{winding_number} + turn_number - 1)}];  // Define a region for each turn
            Winding~{winding_number} += Region[{(iCOND~{winding_number} + turn_number - 1)}];  // Add each turn to the respective winding region
            Air += Region[{(AIR_COND + iCOND~{winding_number} + turn_number - 1)}];  // Include air surrounding each turn in the overall air domain
        EndFor
    EndFor

    // Domain Definitions
    // Non-Conducting Domain (Air and Insulation)
    DomainCC = Region[{Air}];
    DomainCC += Region[{Bobbin_Insulation}];
    DomainCC += Region[{Layer_Insulation}];
    // DomainCC += Region[{ConductorInsulation}];
    For n In {1:n_windings}
        DomainCC += Region[{Cond_Insulation~{n}}];  // Initialize a region for the winding
    EndFor
    //DomainCC += Region[{Cond_Insulation}];

    // Add the Core region to the appropriate domain region
    /*If (Flag_ground_OutBoundary)
        If (Flag_ground_core || Flag_excite_core)
            // If both the core and the outer boundary are grounded, add the core to DomainC
            DomainC += Region[{Core}];
        Else
            // If the outer boundary is grounded but not the core, the core should be floating
            DomainCC += Region[{Core}];
        EndIf
    Else
        // If the outer boundary is NOT grounded
        If (Flag_ground_core || Flag_excite_core)
            DomainC += Region[{Core}];
        Else
            DomainCC += Region[{Core}];
        EndIf
    EndIf*/

    // Add the Core region to the appropriate domain region
    If (Flag_ground_OutBoundary)
        If (Flag_excite_core)
            // If both the core and the outer boundary are grounded, add the core to DomainC
            DomainC += Region[{Core}];
        Else
            // If the outer boundary is grounded but not the core, the core should be floating
            DomainCC += Region[{Core}];
        EndIf
    Else
        // If the outer boundary is NOT grounded
        If (Flag_excite_core)
            DomainC += Region[{Core}];
        Else
            DomainCC += Region[{Core}];
        EndIf
    EndIf

    // Boundary Conditions (Outer boundary to define potential)
    // Neumann boundary conditions should be applied at the left du its symmetry
    // Neumann boundary conditions specify the derivative of the function ( here is the derivative of the voltage(electric flux)) normal to the boundary.
    // This means that there is no electric flux across a boundary
    // This implies insulation or symmetry (no flow of electric charge across the boundary).
    // Note that Since there are no non-homogeneous Neumann conditions in our code, Sur_Neu_Ele can be defined as empty.
    // Otherwise, surface with imposed non-homogeneous Neumann boundary conditions (on n.d = -n . (epsilon grad v)), so nd[] function should be defined in our
    // code as (epsilon[] * Norm[{d u0}]), look to the solver.
    // non-homogeneous indicates that the derivative of the voltage (electric flux) is not equal to zero
    // But as we do not have non-homogeneous boundary conditions, this will not affect our results here.
    Sur_Neu_Ele = Region[{BOUND_LEFT}];
    // Dirichlet boundary conditions specify the value of the function (here is the voltage) directly on the boundary.
    // For instance, setting a boundary to be at 0 volts (ground).
    // OUTBND does not include the BOUND_TOP_LEFT and BOUND_BOT_LEFT.
    OuterBoundary = Region[{OUTBND}];

    // Add Winding Regions to the Core Domain
    For n In {1:n_windings}
        DomainC += Region[{Winding~{n}}];
    EndFor

    // Initialize the Main Domain
    For n In {1:n_windings}
      DomainCond~{n} += Region[{Winding~{n}}] ;
    EndFor

    // Define the Conductor Domain for Each Winding (Turn-wise Capacitance Calculation)
    For winding_number In {1:n_windings}
        DomainCond~{winding_number} = Region[{Winding~{winding_number}}];  // Assign each winding to a conductor domain
        // Define the Conductor Domain for Each Turn within the Winding
        nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;  // Calculate the number of turns in each winding
        For turn_number In {1:nbturns~{winding_number}}
            DomainCond~{winding_number}~{turn_number} = Region[{Turn~{winding_number}~{turn_number}}];  // Assign each turn to a conductor domain
        EndFor
    EndFor

    Domain = Region[{DomainC, DomainCC}];
    // Dummy Region for post-processing with functions
    DomainDummy = Region[12345];
}
// Excitation
// ----------------------
Function {

  // Define relative permittivity (epsilon) for all regions
  AreaCell[#{Air, Bobbin_Insulation, Core}] = 1.;
  For n In {1:n_windings}
      AreaCell[#{Winding~{n}}] = 1.;
  EndFor
  // Core area: needed for finding the avg voltage in the core
  SurfCore[] = SurfaceArea[]{CORE_PN} ;
  // Materials
  er_air = 1;
  er_core = 100000;
  //er_cond_insulation = 3;
  For n In {1:n_windings}
      epsilon[#{Cond_Insulation~{n}}] = e0 * er_turns_insulation~{n};
  EndFor
  epsilon[#{Air}] = e0 * er_air;
  epsilon[#{Core}] = e0 * er_core;
  epsilon[#{Bobbin_Insulation}] = e0 * er_bobbin;
  epsilon[#{Layer_Insulation}] = e0 * er_layer_insulation;
  //epsilon[#{Cond_Insulation}] = e0 * er_cond_insulation;
  // The winding permittivity is set to 1, but it does not play any role in the simulation.
  For winding_number In {1:n_windings}
      er_winding~{winding_number} = 1;
      nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;
      For turn_number In {1:nbturns~{winding_number}}
          er_turn~{winding_number}~{turn_number} = er_winding~{winding_number};
          epsilon[#{Turn~{winding_number}~{turn_number}}] = e0 * er_turn~{winding_number}~{turn_number};
      EndFor
  EndFor

  If (Flag_voltage)
      // Assign a voltage to each turn of each turn
      For winding_number In {1:n_windings}
        nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;
        For turn_number In {1:nbturns~{winding_number}}
          //Val_Potential_Turn~{winding_number}~{turn_number} = Voltage~{winding_number} * (turn_number / nbturns~{winding_number});
          Val_Potential_Turn~{winding_number}~{turn_number} = Voltage~{winding_number}~{turn_number};
        EndFor
      EndFor
  EndIf

  If (Flag_charge)
      // Assign a charge to each turn of each turn
      For winding_number In {1:n_windings}
        nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;
        For turn_number In {1:nbturns~{winding_number}}
          Cha_Potential_Turn~{winding_number}~{turn_number} = Charge~{winding_number}~{turn_number};
        EndFor
      EndFor
  EndIf
}

// Definition of Constraints for Electrostatic Analysis
// ----------------------
Constraint {

    { Name Electrostatic_Potential ; Type Assign;
        Case {
            If (Flag_voltage)
                // Assigning a fixed potential to each turn within each winding
                For winding_number In {1:n_windings}
                    nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;
                    For turn_number In {1:nbturns~{winding_number}}
                        { Region Turn~{winding_number}~{turn_number} ; Type Assign ; Value Val_Potential_Turn~{winding_number}~{turn_number}; }
                    EndFor
                EndFor
            EndIf

            /*If(Flag_ground_core)
                { Region Core ; Type Assign ; Value 0; }
            ElseIf(Flag_excite_core)
                { Region Core ; Type Assign ; Value v_core; }
            EndIf*/
            If(Flag_excite_core)
                { Region Core ; Type Assign ; Value v_core; }
            EndIf
        }
    }
    { Name SetArmatureCharge; Type Assign;
        Case {
            If (Flag_charge)
                // Assigning a fixed charge to each turn within each winding
                For winding_number In {1:n_windings}
                    nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;
                    For turn_number In {1:nbturns~{winding_number}}
                        { Region Turn~{winding_number}~{turn_number} ; Type Assign ; Value Cha_Potential_Turn~{winding_number}~{turn_number}; }
                    EndFor
                EndFor
                If(Flag_excite_core)
                    { Region Core ; Type Assign ; Value v_core; }
                EndIf
            EndIf

        }
    }
    { Name Dirichlet_Ele; Type Assign;
        Case {
           If(Flag_ground_OutBoundary)
              { Region OuterBoundary; Type Assign ; Value 0; }
           EndIf
        }
    }

}

// ----------------------
// Include solver file(s)
// ----------------------
Include "solver_electrostatic.pro"
// Resolution
Resolution {
  { Name EleSta_v;
    System {
      { Name Sys_Ele; NameOfFormulation Electrostatic_Potential; }
    }
    Operation {
      CreateDir[DirResValsCore];
      CreateDir[DirResValsCapacitance];
      For n In {1:n_windings}
          CreateDir[DirResValsWinding~{n}];
          CreateDir[DirResValsCharge~{n}];
          CreateDir[DirResValsVoltage~{n}];
          CreateDir[DirResValsCapacitanceFromQV~{n}];
          nbturns~{n} = NbrCond~{n} / SymFactor;
          For turn In {1:nbturns~{n}}
             CreateDir[DirResValsTurn~{turn}];
          EndFor

      EndFor
      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
      PostOperation[Map_local] ;
      PostOperation[Get_global] ;
    }
  }
}

// Post processing
PostProcessing {
 { Name EleSta; NameOfFormulation Electrostatic_Potential; NameOfSystem Sys_Ele;
    PostQuantity {
      // Voltage
      { Name u0; Value {
          Term { [ {u0} ]; In Domain; Jacobian Vol; }
        }
      }
      // Electric field
      { Name e; Value {
          Term { [ -{d u0} ]; In Domain; Jacobian Vol; }
        }
      }
      // Electric field
      { Name Welocal; Value {
          Term { [ epsilon[] / 2. * SquNorm[{d u0}] ]; In Domain; Jacobian Vol; }
        }
      }
      // electric field density
      { Name d; Value {
          Term { [ -epsilon[] * {d u0} ]; In Domain; Jacobian Vol; }
        }
      }
      // Magnitude of Electric Field (|E|)
      { Name MagE; Value {
          Term { [ Norm[{d u0}] ]; In Domain; Jacobian Vol; }
        }
      }

      // Magnitude of Displacement Field (|D|)
      { Name MagD; Value {
          Term { [ Norm[epsilon[] * {d u0}] ]; In Domain; Jacobian Vol; }
        }
      }

      // Real and Imaginary Components of Electric Field
      { Name Real_E; Value {
          Term { [ Re[{d u0}] ]; In DomainCC; Jacobian Vol; }
        }
      }

      { Name Imag_E; Value {
          Term { [ Im[{d u0}] ]; In DomainCC; Jacobian Vol; }
        }
      }

      // Real and Imaginary Components of Displacement Field
      { Name Real_D; Value {
          Term { [ Re[epsilon[] * {d u0}] ]; In DomainCC; Jacobian Vol; }
        }
      }

      { Name Imag_D; Value {
          Term { [ Im[epsilon[] * {d u0}] ]; In DomainCC; Jacobian Vol; }
        }
      }
      // Global voltages
      { Name U; Value {
          Term { [ {U} ]; In DomainC; }
        }
      }
      // Global charges
      { Name Q; Value {
          Term { [ CoefGeo * {Q} ]; In DomainC; }
        }
      }
      // Stored Energy
      { Name energy; Value {
          Integral { Type Global;
            [ CoefGeo * epsilon[] / 2. * SquNorm[{d u0}] ];
            In Domain; Jacobian Vol; Integration II;
          }
	  }
      }

      { Name Charge; Value {
              Integral { [ CoefGeo * epsilon[] * Norm[{d u0}] ];
                In Domain; Jacobian Vol; Integration II; }
            }
       }
      // Average voltage of the core
      { Name Avg_Voltage_Core; Value {Integral { [ {u0} / SurfCore[] ]; In Region[{Core}]; Jacobian Sur; Integration II; }}}

      // Voltages of every turn to print them in different files
      For winding_number In {1:n_windings}
        nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

        // Print charge for each turn in the winding
        For turn_number In {1:nbturns~{winding_number}}
            { Name U~{winding_number}~{turn_number}; Value {
                 Term { [ {U} ]; In Turn~{winding_number}~{turn_number}; }
                 }
            }
        EndFor
      EndFor
      // Charges of every turn to print them in different files
      For winding_number In {1:n_windings}
        nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

        // Print charge for each turn in the winding
        For turn_number In {1:nbturns~{winding_number}}
            { Name Q~{winding_number}~{turn_number}; Value {
                 Term { [ CoefGeo * {Q} ]; In Turn~{winding_number}~{turn_number}; }
                 }
            }
        EndFor
      EndFor


      // Capacitance Calculation Between Turns for Each Winding ( from charges)
      For winding_number In {1:n_windings}
            nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

            // Loop through each turn to define its capacitance relative to other turns
            For turn_number1 In {1:nbturns~{winding_number}}
                    { Name Capacitance_Turn_Core~{winding_number}~{turn_number1}; Value {
                          Term { Type Global;
                                 [ $Charge_Core / ($U~{winding_number}~{turn_number1} - $Avg_Voltage_Core) ];
                                 In Turn~{winding_number}~{turn_number1}; }
                        }
                    }

                For turn_number2 In {1:nbturns~{winding_number}}
                    // Calculate capacitance for every pair of turns, where turn_number1 is the reference voltage turn
                    If (turn_number1 == turn_number2)
                        // Self-capacitance (Q / V, with V being the turn voltage)
                        { Name Capacitance_Turn_Self~{winding_number}~{turn_number1}; Value {
                              Term { Type Global;
                                     [ $Q~{winding_number}~{turn_number1}  / $U~{winding_number}~{turn_number1} ];
                                     In Turn~{winding_number}~{turn_number1}; }
                            }
                        }
                    Else
                        // Mutual capacitance between different turns (Q_turn2 / V_turn1)
                        { Name Capacitance_Turn~{winding_number}~{turn_number1}~{turn_number2}; Value {
                              Term { Type Global;
                                     [ $Q~{winding_number}~{turn_number2} / $U~{winding_number}~{turn_number1} ];
                                     In Turn~{winding_number}~{turn_number1}; }
                            }
                        }
                    EndIf
                EndFor
            EndFor
      EndFor

      // Capacitance Calculation Between Turns of Different Windings
      For winding_number1 In {1:n_windings}
        nbturns1~{winding_number1} = NbrCond~{winding_number1} / SymFactor;

        For winding_number2 In {1:n_windings}
            If (winding_number1 != winding_number2)
                nbturns2~{winding_number2} = NbrCond~{winding_number2} / SymFactor;

                For turn_number1 In {1:nbturns1~{winding_number1}}
                    For turn_number2 In {1:nbturns2~{winding_number2}}
                        // Mutual capacitance between different turns in different windings (Q_turn2 / V_turn1)
                        { Name Capacitance_Turn_Windings~{winding_number1}~{turn_number1}~{winding_number2}~{turn_number2}; Value {
                              Term { Type Global;
                                     [ $Q~{winding_number2}~{turn_number2} / $U~{winding_number1}~{turn_number1} ];
                                     In Turn~{winding_number1}~{turn_number1}; }
                            }
                        }
                    EndFor
                EndFor
            EndIf
        EndFor
      EndFor

      // Calculate capacitance through the stored energy
      If (Flag_voltage)
        For winding_number In {1:n_windings}
            nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

            // Loop through each turn as the reference turn
            For turn_number1 In {1:nbturns~{winding_number}}
                // Loop through every other turn, including the reference turn itself
                // Calculate the capacitance between the turn and the core
                { Name Capacitance_Between_Turns_Core~{winding_number}~{turn_number1}; Value {
                    Term { Type Global;
                      [ 2 * $energy_Component / ((Val_Potential_Turn~{winding_number}~{turn_number1} - $Avg_Voltage_Core) * (Val_Potential_Turn~{winding_number}~{turn_number1} - $Avg_Voltage_Core)) ];
                      In DomainCond~{winding_number}~{turn_number1}; }
                    }
                }
                For turn_number2 In {1:nbturns~{winding_number}}
                    // Calculate the voltage difference between the two turns
                    VoltageDifference_Pair~{winding_number}~{turn_number1}~{turn_number2} = Val_Potential_Turn~{winding_number}~{turn_number2} - Val_Potential_Turn~{winding_number}~{turn_number1};

                    // Calculate the capacitance between the two turns
                    /*{ Name Capacitance_Between_Turns~{winding_number}~{turn_number1}~{turn_number2}; Value {
                        Integral { [ CoefGeo * epsilon[]  * SquNorm[{d u0}] / (VoltageDifference_Pair~{winding_number}~{turn_number1}~{turn_number2} * VoltageDifference_Pair~{winding_number}~{turn_number1}~{turn_number2}) ];
                          In Domain; Jacobian Vol; Integration II; }
                      }
                    }*/
                    // using $... will increase the computation of our simulation
                    { Name Capacitance_Between_Turns~{winding_number}~{turn_number1}~{turn_number2}; Value {
                        Term { Type Global;
                          [ 2 * $energy_Component / (VoltageDifference_Pair~{winding_number}~{turn_number1}~{turn_number2} * VoltageDifference_Pair~{winding_number}~{turn_number1}~{turn_number2}) ];
                          In DomainCond~{winding_number}~{turn_number1}; }
                        }
                    }

                EndFor
                // Calculate capacitance between this turn and turns in other windings
                For other_winding_number In {1:n_windings}
                  If (other_winding_number != winding_number)
                    nbturns~{other_winding_number} = NbrCond~{other_winding_number} / SymFactor;

                    // Loop through each turn in the other winding
                    For turn_number2 In {1:nbturns~{other_winding_number}}
                      // Calculate the voltage difference between the turn in the current winding and a turn in the other winding
                      VoltageDifference_Cross~{winding_number}~{turn_number1}~{other_winding_number}~{turn_number2} =
                      Val_Potential_Turn~{other_winding_number}~{turn_number2} - Val_Potential_Turn~{winding_number}~{turn_number1};
                      { Name Capacitance_Cross~{winding_number}~{turn_number1}~{other_winding_number}~{turn_number2} ;
                       Value { Term { Type Global;
                        [ 2 * $energy_Component / (VoltageDifference_Cross~{winding_number}~{turn_number1}~{other_winding_number}~{turn_number2} *
                         VoltageDifference_Cross~{winding_number}~{turn_number1}~{other_winding_number}~{turn_number2}) ] ; In DomainCond~{winding_number}~{turn_number1}; } } }
                    EndFor
                  EndIf
                EndFor
            EndFor
        EndFor
      EndIf

    } // Quantity Section
 } // Name EleSta_v
} // PostProcessing


Include "fields_electrostatic.pro";
Include "values_electrostatic.pro";
