PostOperation Get_global UsingPost EleSta {
  // energy stored in air
  Print[ energy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Energy_Stored_Component.dat"], LastTimeStepOnly, StoreInVariable $energy_Component];
  Print[ energy[Air], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Energy_Stored_Air.dat"], LastTimeStepOnly, StoreInVariable $energy_Air];
  Print[ energy[Core], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Energy_Stored_Core.dat"], LastTimeStepOnly, StoreInVariable $energy_Core];
  // average voltage of the core
  Print[ Avg_Voltage_Core[Core], OnGlobal, Format TimeTable, File > StrCat[DirResCirc,"Avg_Core_voltage.dat"], LastTimeStepOnly, StoreInVariable $Avg_Voltage_Core];

  // Charges
  Print[ Charge[Air], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"charge_Air.dat"], LastTimeStepOnly, StoreInVariable $Charge_Air];
  Print[ Charge[Core], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"charge_Core.dat"], LastTimeStepOnly, StoreInVariable $Charge_Core];
  // voltage and charges on windings
  For n In {1:n_windings}
         Print[ U, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_%g.dat"], n] , LastTimeStepOnly];
  EndFor
  For n In {1:n_windings}
         Print[ Q, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"Q_%g.dat"], n] , LastTimeStepOnly];
  EndFor

  // print charge and voltage for each turn in separate file (Just another way)
  For winding_number In {1:n_windings}
    nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

    // Print charge for each turn in the winding
    For turn_number In {1:nbturns~{winding_number}}
        Print[ U, OnRegion Turn~{winding_number}~{turn_number}, Format TimeTable,
               File > Sprintf[StrCat[DirResCirc, "U_%g_%g.dat"], winding_number, turn_number], LastTimeStepOnly];
        Print[ Q, OnRegion Turn~{winding_number}~{turn_number}, Format TimeTable,
               File > Sprintf[StrCat[DirResCirc, "Q_%g_%g.dat"], winding_number, turn_number], LastTimeStepOnly];
    EndFor
 EndFor

  // print voltages for each turn
  For winding_number In {1:n_windings}
    nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

    // Print voltage for each turn in the winding
    For turn_number In {1:nbturns~{winding_number}}
        Print[ U~{winding_number}~{turn_number},
                       OnRegion Turn~{winding_number}~{turn_number}, Format Table, File > Sprintf[StrCat[DirResValsVoltage~{winding_number},
                        "voltage_%g_%g.dat"], winding_number, turn_number], LastTimeStepOnly, StoreInVariable $U~{winding_number}~{turn_number}];
    EndFor
  EndFor

  // print charges for each turn
  For winding_number In {1:n_windings}
    nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

    // Print charge for each turn in the winding
    For turn_number In {1:nbturns~{winding_number}}
        Print[ Q~{winding_number}~{turn_number},
                       OnRegion Turn~{winding_number}~{turn_number}, Format Table, File > Sprintf[StrCat[DirResValsCharge~{winding_number},
                        "Charge_%g_%g.dat"], winding_number, turn_number], LastTimeStepOnly, StoreInVariable $Q~{winding_number}~{turn_number}];
    EndFor
  EndFor


  // Calculate and print capacitance through the stored energy
  If (Flag_voltage)
      For winding_number In {1:n_windings}
          nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

          // Loop through each turn as the reference turn
          For turn_number1 In {1:nbturns~{winding_number}}
              // Loop through every other turn, including the reference turn itself
              Print[ Capacitance_Between_Turns_Core~{winding_number}~{turn_number1},
                       OnRegion DomainCond~{winding_number}~{turn_number1}, Format Table, File > Sprintf[StrCat[DirResValsTurn~{turn_number1},
                        "C_%g_%g_Core.dat"], winding_number, turn_number1], LastTimeStepOnly];
              For turn_number2 In {1:nbturns~{winding_number}}
                  // Print the calculated capacitance between each pair of turns
                  /*Print[ Capacitance_Between_Turns~{winding_number}~{turn_number1}~{turn_number2}[Domain],
                         OnGlobal, Format TimeTable,
                         File > Sprintf[StrCat[DirResValsTurn~{turn_number1}, "C_%g_%g_%g.dat"], winding_number, turn_number1, turn_number2] ];*/
                  // Using this for $...
                  Print[ Capacitance_Between_Turns~{winding_number}~{turn_number1}~{turn_number2},
                       OnRegion DomainCond~{winding_number}~{turn_number1}, Format Table, File > Sprintf[StrCat[DirResValsTurn~{turn_number1},
                        "C_%g_%g_%g.dat"], winding_number, turn_number1, turn_number2], LastTimeStepOnly];
              EndFor
              // Print capacitance between this turn and turns in other windings
               For other_winding_number In {1:n_windings}
                  If (other_winding_number != winding_number)
                    nbturns~{other_winding_number} = NbrCond~{other_winding_number} / SymFactor;

                    // Loop through each turn in the other winding
                    For turn_number2 In {1:nbturns~{other_winding_number}}
                      // Print the calculated capacitance between each pair of turns in different windings
                      Print[ Capacitance_Cross~{winding_number}~{turn_number1}~{other_winding_number}~{turn_number2},
                       OnRegion DomainCond~{winding_number}~{turn_number1}, Format Table, File > Sprintf[StrCat[DirResValsTurn~{turn_number1},
                        "C_Cross_%g_%g_%g_%g.dat"], winding_number, turn_number1, other_winding_number, turn_number2], LastTimeStepOnly];
                    EndFor
                  EndIf
               EndFor
          EndFor
      EndFor
  EndIf



  // Capacitances from QV relation on each turn of the same winding
  If (Flag_voltage)
      // Capacitance Calculation Between Turns for Each Winding (from charges)
        For winding_number In {1:n_windings}
            nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;

            // Loop through each turn to define its capacitance relative to other turns
            For turn_number1 In {1:nbturns~{winding_number}}
                // Print the self-capacitance value
                Print[Capacitance_Turn_Core~{winding_number}~{turn_number1},
                      OnRegion Turn~{winding_number}~{turn_number1},
                      Format Table,
                      File > Sprintf[StrCat[DirResValsCapacitanceFromQV~{winding_number}, "C_%g_%g_Core.dat"], winding_number, turn_number1],
                      LastTimeStepOnly];
                For turn_number2 In {1:nbturns~{winding_number}}
                    // Calculate capacitance for every pair of turns, where turn_number1 is the reference voltage turn
                    If (turn_number1 == turn_number2)
                        // Print the self-capacitance value
                        Print[Capacitance_Turn_Self~{winding_number}~{turn_number1},
                              OnRegion Turn~{winding_number}~{turn_number1},
                              Format Table,
                              File > Sprintf[StrCat[DirResValsCapacitanceFromQV~{winding_number}, "C_%g_%g_%g.dat"], winding_number, turn_number1, turn_number2],
                              LastTimeStepOnly];
                    Else
                        // Print the mutual capacitance value
                        Print[Capacitance_Turn~{winding_number}~{turn_number1}~{turn_number2},
                              OnRegion Turn~{winding_number}~{turn_number1},
                              Format Table,
                              File > Sprintf[StrCat[DirResValsCapacitanceFromQV~{winding_number}, "C_%g_%g_%g.dat"], winding_number, turn_number1, turn_number2],
                              LastTimeStepOnly];
                    EndIf
                EndFor
            EndFor
        EndFor
  EndIf

  // Print Capacitance Between Turns of Different Windings
  For winding_number1 In {1:n_windings}
    nbturns1~{winding_number1} = NbrCond~{winding_number1} / SymFactor;

    For winding_number2 In {1:n_windings}
        If (winding_number1 != winding_number2)
            nbturns2~{winding_number2} = NbrCond~{winding_number2} / SymFactor;

            For turn_number1 In {1:nbturns1~{winding_number1}}
                For turn_number2 In {1:nbturns2~{winding_number2}}
                    // Print the mutual capacitance value between different windings
                    Print[Capacitance_Turn_Windings~{winding_number1}~{turn_number1}~{winding_number2}~{turn_number2},
                          OnRegion Turn~{winding_number1}~{turn_number1},
                          Format Table,
                          File > Sprintf[StrCat[DirResValsCapacitanceFromQV~{winding_number1},
                          "C_%g_%g_%g_%g.dat"], winding_number1, turn_number1, winding_number2, turn_number2],
                          LastTimeStepOnly];
                EndFor
            EndFor
        EndIf
    EndFor
  EndFor



}