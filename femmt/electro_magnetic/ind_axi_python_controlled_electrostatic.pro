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
AIR                     = 110000;
AIR_EXT                 = 110001;
AIR_COND                = 1000000;
CORE_PN                 = 120000;
ExtGmsh                 = ".pos";

//physical numbers of conductors in n transformer
For n In {1:n_windings}
       iCOND~{n} = 130000 + 1000*(n-1);
       istrandedCOND~{n} = 150000 + 1000*(n-1);
EndFor

// ----------------------
// Groups (Physical Domains Setup for Electrostatics)
// ----------------------
Group {

    Air  = Region[{AIR, AIR_EXT}];

    // Core
    Core = Region[{}];
    // Core Domain
    For n In {1:nCoreParts}
      CorePart~{n} = Region[{(CORE_PN+n-1)}];
      Core += Region[{(CORE_PN+n-1)}];
    EndFor

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
    // Non-Conducting Domain (Air)
    DomainCC = Region[{Air}];

    // Add the Core region to the non-conducting domain region
    If (Flag_ground_OutBoundary)
        // If the outer boundary is grounded, the core should be floating
        DomainCC += Region[{Core}];
    Else
        // If the outer boundary is NOT grounded, handle the core according to the core flag
        If (Flag_ground_core)
            DomainC += Region[{Core}];
        Else
            DomainCC += Region[{Core}];
        EndIf
    EndIf

    // Boundary Conditions (Outer boundary to define potential)
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
  AreaCell[#{Air, Core}] = 1.;
  For n In {1:n_windings}
      AreaCell[#{Winding~{n}}] = 1.;
  EndFor
  // Core area: needed for finding the avg voltage in the core
  SurfCore[] = SurfaceArea[]{CORE_PN} ;
  // Materials
  er_air = 1;
  er_core = 5000;
  epsilon[#{Air}] = e0 * er_air;
  epsilon[#{Core}] = e0 * er_core;
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

            If(Flag_ground_core)
                { Region Core ; Type Assign ; Value 0; }
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
            EndIf
            If(Flag_ground_core)
                { Region Core ; Type Assign ; Value 0; }
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
          Term { [ {Q} ]; In DomainC; }
        }
      }
      // These C need to be reviewed
      { Name C; Value {
          Term { [ {Q}/{U} ]; In DomainC; }
        }
      }
      // Stored Energy in the whole domain
      { Name energy_Component; Value {
              Integral { [ CoefGeo * epsilon[] / 2. * SquNorm[{d u0}] ];
                In DomainCC; Jacobian Vol; Integration II; }
            }
          }
      // Stored Energy in the air
      { Name energy_Air; Value {
          Integral { Type Global;
            [ CoefGeo * epsilon[] / 2. * SquNorm[{d u0}] ];
            In Air; Jacobian Vol; Integration II;
          }
	  }
      }
      // Stored Energy in the whole domain
      { Name energy_Core; Value {
              Integral { [ CoefGeo * epsilon[] / 2. * SquNorm[{d u0}] ];
                In Core; Jacobian Vol; Integration II; }
            }
          }
      // Charges on the non conducting domain
      { Name Charge; Value {
              Integral { [ CoefGeo * epsilon[] * Norm[{d u0}] ];
                In DomainCC; Jacobian Vol; Integration II; }
            }
       }
      // Average voltage of the core
      { Name Avg_Voltage_Core; Value {Integral { [ {u0} / SurfCore[] ]; In Region[{Core}]; Jacobian Sur; Integration II; }}}

      If (Flag_voltage)
         // Calculate capacitance through the stored energy
        For winding_number In {1:n_windings}
            nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

            // Loop through each turn as the reference turn
            For turn_number1 In {1:nbturns~{winding_number}}
                // Loop through every other turn, including the reference turn itself
                // Calculate the capacitance between the turn and the core
                { Name Capacitance_Between_Turns_Core~{winding_number}~{turn_number1}; Value {
                    Term { Type Global;
                      [ 2 * $energy_Component / ((Val_Potential_Turn~{winding_number}~{turn_number1} - $Avg_Voltage_Core) * (Val_Potential_Turn~{winding_number}~{turn_number1} - $Avg_Voltage_Core)) ];
                      In Domain; }
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
                          In Domain; }
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
                         VoltageDifference_Cross~{winding_number}~{turn_number1}~{other_winding_number}~{turn_number2}) ] ; In Domain; } } }
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
