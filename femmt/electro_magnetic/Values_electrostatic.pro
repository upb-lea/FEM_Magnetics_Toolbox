PostOperation Get_global UsingPost EleSta {
  // energy stored in air
  Print[ energy_Component[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Energy_Stored_Component.dat"], LastTimeStepOnly, StoreInVariable $energy_Component];
  Print[ energy_Air[Air], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Energy_Stored_Air.dat"], LastTimeStepOnly, StoreInVariable $energy_Air];
  Print[ Avg_Voltage_Core[Core], OnGlobal, Format TimeTable, File > StrCat[DirResCirc,"Avg_Core_voltage.dat"], LastTimeStepOnly, StoreInVariable $Avg_Voltage_Core];

  // Charge
  Print[ Charge[Air], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"charge.dat"], LastTimeStepOnly, StoreInVariable $charge];
  // Capacitance
  If (Flag_voltage)
      // Calculate and print capacitance through the stored energy
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


  For n In {1:n_windings}
         Print[ U, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_%g.dat"], n] , LastTimeStepOnly];
  EndFor
  For n In {1:n_windings}
         Print[ Q, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"Q_%g.dat"], n] , LastTimeStepOnly];
  EndFor
  For n In {1:n_windings}
         Print[ C, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"C_%g.dat"], n] , LastTimeStepOnly];
  EndFor
  //Print[ u0_avg_core, OnRegion Core, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"Core_voltage.dat"], n] , LastTimeStepOnly];
  //Print[ u0_avg_core[ Core ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Core_voltage.dat"]];
  //Print[ u0_avg_core, OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Core_voltage.dat"]];

}