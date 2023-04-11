// === Global/Integrated Quantities ===

PostOperation Get_global UsingPost MagDyn_a {

  // Losses


  // Windings Total
  // Solid
  //Print[ SoF[ DomainC ], OnGlobal, Format TimeTable,  File > Sprintf("results/SF_iron.dat")] ; // TODO: Complex power
  Print[ j2F[ Winding~{1} ], OnGlobal, Format TimeTable, File > StrCat[DirResVals, "j2F_1.dat"]] ;

  If(Flag_Transformer)
    For n In {2:n_windings}
        Print[ j2F[ Winding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals, "j2F_%g.dat"], n]] ;
    EndFor
  EndIf
  /*
  For n In {1:n_windings}
      Print[ j2F[ Winding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals, "j2F_%g.dat"], n]] ;
  EndFor
  */


  //Print[ j2F[ Winding3 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2F_3.dat"]] ;
  // Stranded
  //Print[ SoH[ StrandedWinding1 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"SH_1.dat"] ] ;  // TODO: Complex power
  //Print[ SoH[ DomainS ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"SH.dat"] ] ;
  //Print[ j2H[ DomainS ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"LossesStrandedWindings.dat"] ] ;
  Print[ j2H[ StrandedWinding~{1} ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_1.dat"] ] ;
  //Print[ j2H[ StrandedWinding1 ], OnGlobal, Format Table];

  If(Flag_Transformer)
    For n In {2:n_windings}
        Print[ j2H[ StrandedWinding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"j2H_%g.dat"], n]] ;
    EndFor
  EndIf
  /*
  For n In {1:n_windings}
      Print[ j2H[ StrandedWinding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"j2H_%g.dat"], n]] ;
  EndFor
  */

  //Print[ j2H[ StrandedWinding2 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_2.dat"] ] ;
  //Print[ j2H[ StrandedWinding3 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_3.dat"] ] ;
  //Print[ j2H[ StrandedWinding2 ], OnGlobal, Format Table];
  //Print[ j2Hskin[StrandedWinding1],   OnGlobal , Format Table];
  //Print[ j2Hprox[StrandedWinding1],   OnGlobal , Format Table];
  //Print[ j2Hskin[StrandedWinding2],   OnGlobal , Format Table];
  //Print[ j2Hprox[StrandedWinding2],   OnGlobal , Format Table];

  // Single Turns

  If(Flag_HomogenisedModel~{1}) // Differentiate fine und hom
    For isF In {1:nbturns~{1}}
      Print[ j2H[ TurnStrand~{1}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{1},"Losses_turn_%g.dat"], isF] ] ;
    EndFor
  Else
    For isF In {1:nbturns~{1}}
      Print[ j2F[ Turn~{1}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{1},"Losses_turn_%g.dat"], isF] ] ;
      //Print[ az_int[ Turn~{1}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{1},"a_turn_%g.dat"], isF] ] ;
    EndFor
  EndIf

  If(Flag_Transformer)
    For n In {2:n_windings}
        If(Flag_HomogenisedModel~{n}) // Differentiate fine und hom
          For isF In {1:nbturns~{n}}
              Print[ j2H[ TurnStrand~{n}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{n},"Losses_turn_%g.dat"], isF] ] ;
          EndFor
        Else
          For isF In {1:nbturns~{n}}
              Print[ j2F[ Turn~{n}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{n},"Losses_turn_%g.dat"], isF] ] ;
          EndFor
        EndIf
    EndFor
  EndIf

   /*
   For n In {1:n_windings}
       If(Flag_HomogenisedModel~{n}) // Differentiate fine und hom
          For isF In {1:nbturns~{n}}
              Print[ j2H[ TurnStrand~{n}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals~{n},"Losses_turn_%g.dat"], isF] ] ;
          EndFor
       Else
          For isF In {1:nbturns~{n}}
              Print[ j2F[ Turn~{n}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals~{n},"Losses_turn_%g.dat"], isF] ] ;
          EndFor
       EndIf
   EndFor
   */

  /*
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


  If(Flag_Three_Transformer)
    If(Flag_HomogenisedModel3) // Differentiate fine und hom
      For isF In {1:nbturns3}
        Print[ j2H[ TurnStrand3~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsTertiary,"Losses_turn_%g.dat"], isF] ] ;
      EndFor
    Else
      For isF In {1:nbturns3}
        Print[ j2F[ Turn3~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsTertiary,"Losses_turn_%g.dat"], isF] ] ;
      EndFor
    EndIf
  EndIf
  */


  // Core

  // Eddy Current Losses according to sigma in core/iron
  Print[ j2F[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"CoreEddyCurrentLosses.dat"]] ;

  // Hysteresis Losses according to complex permeability in core/iron
  Print[ p_hyst[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"p_hyst.dat"]] ;// Core losses

  // Steinmetz Core Losses
  If(Flag_Generalized_Steinmetz_loss)
    Print[ piGSE[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"piGSE.dat"]] ;// Core losses
    Print[ piGSE[ Iron ], OnGlobal, Format Table];
  EndIf

  If(Flag_Steinmetz_loss)
    Print[ pSE[ Iron ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"pSE.dat"]] ;// Core losses
    Print[ pSE[ Iron ], OnGlobal, Format Table];
  EndIf


  // Stored Energy
  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"ME.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];
  // Print[ MagEnergy[Iron], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"ME_iron.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];
  // Print[ MagEnergy[Air], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"ME_air.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];
  // Print[ MagEnergy[Winding1], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"ME_winding1.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy];

  // Flux (Linkage)

  Print[ Flux_Linkage~{1}[DomainCond~{1}], OnGlobal, Format Table, File > StrCat[DirResVals,"Flux_Linkage_1.dat"]];
  // Print[ Flux_Linkage_1[DomainCond1], OnGlobal, Format Table];

  If(Flag_Transformer)
    For n In {2:n_windings}
        Print[ Flux_Linkage~{n}[DomainCond~{n}], OnGlobal, Format Table, File > Sprintf[StrCat[DirResVals,"Flux_Linkage_%g.dat"], n]];
        //Print[ Flux_Linkage~{n}[DomainCond~{n}], OnGlobal, Format Table];
    EndFor
  EndIf
  /*
  For n In {1:n_windings}
      Print[ Flux_Linkage~{n}[DomainCond~{n}], OnGlobal, Format Table, File > Sprintf[StrCat[DirResVals,"Flux_Linkage_%g.dat"], n]];
      //Print[ Flux_Linkage~{n}[DomainCond~{n}], OnGlobal, Format Table];
  EndFor
  */


  /*
  If(Flag_Transformer)
    Print[ Flux_Linkage_2[DomainCond2], OnGlobal, Format Table, File > StrCat[DirResVals,"Flux_Linkage_2.dat"]];
    // Print[ Flux_Linkage_2[DomainCond2], OnGlobal, Format Table];
  EndIf

  If(Flag_Three_Transformer)
    Print[ Flux_Linkage_3[DomainCond3], OnGlobal, Format Table, File > StrCat[DirResVals,"Flux_Linkage_3.dat"]];
    // Print[ Flux_Linkage_3[DomainCond3], OnGlobal, Format Table];
  EndIf
  */

  // Inductances

  If(Val_EE_1!=0)
    Print[ L_1_1[DomainCond~{1}], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_1_1.dat"]] ;
    Print[ LFromMagEnergy_1_1[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"LFromMagEnergy_1_1.dat"]] ;
  EndIf

  If(Flag_Transformer)
    For n In {2:n_windings}
       If(Val_EE~{n}!=0)
         Print[ L~{n}~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"L_%g_%g.dat"], n, n]] ;
         Print[ LFromMagEnergy~{n}~{n}[Domain], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"LFromMagEnergy_%g_%g.dat"], n, n]] ;
       EndIf
    EndFor
  EndIf
  /*
  For n In {1:n_windings}
      If(Val_EE~{n}!=0)
         Print[ L~{n}~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"L_%g_%g.dat"], n, n]] ;
         Print[ LFromMagEnergy~{n}~{n}[Domain], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"LFromMagEnergy_%g_%g.dat"], n, n]] ;
      EndIf
  EndFor
  */
  /*
  If(Flag_Transformer)
      If(Val_EE_2!=0)
        Print[ L_22[DomainCond2], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_22.dat"]] ;
        Print[ L_22_from_MagEnergy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_22_from_MagEnergy.dat"]] ;
      EndIf
  EndIf

  If(Flag_Three_Transformer)
      If(Val_EE_3!=0)
        Print[ L_33[DomainCond3], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_33.dat"]] ;
        Print[ L_33_from_MagEnergy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"L_33_from_MagEnergy.dat"]] ;
      EndIf
  EndIf
  */

  // Voltage

  Print[ Voltage~{1}[DomainCond~{1}], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Voltage_1.dat"]];
  If(Flag_Transformer)
    For n In {2:n_windings}
        Print[ Voltage~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"Voltage_%g.dat"], n]];
    EndFor
  EndIf
  /*
  For n In {1:n_windings}
      Print[ Voltage~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"Voltage_%g.dat"], n]];
  EndFor
  */

  /*
  If(Flag_Transformer)
    Print[ Voltage_2[DomainCond2], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Voltage_2.dat"]];
  EndIf

  If(Flag_Three_Transformer)
    Print[ Voltage_3[DomainCond3], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Voltage_3.dat"]];
  EndIf
  */

  // Circuit Quantities

  If(!Flag_Circuit)
    If(!Flag_HomogenisedModel~{1})
      Print[ I, OnRegion Winding~{1}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_1_f%g.dat"], Freq] , LastTimeStepOnly];
      Print[ U, OnRegion Winding~{1}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_1_f%g.dat"], Freq] , LastTimeStepOnly];
    Else
      Print[ I, OnRegion StrandedWinding~{1}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_1_f%g.dat"], Freq] , LastTimeStepOnly];
      Print[ U, OnRegion StrandedWinding~{1}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_1_f%g.dat"], Freq] , LastTimeStepOnly];
    EndIf

    If(Flag_Transformer)
       For n In {2:n_windings}
           If(!Flag_HomogenisedModel~{n})
             Print[ I, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_%g_f%g.dat"], n, Freq] , LastTimeStepOnly];
             Print[ U, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_%g_f%g.dat"], n, Freq] , LastTimeStepOnly];
           Else
             Print[ I, OnRegion StrandedWinding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_%g_f%g.dat"], n, Freq] , LastTimeStepOnly];
             Print[ U, OnRegion StrandedWinding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_%g_f%g.dat"], n, Freq] , LastTimeStepOnly];
           EndIf
       EndFor
    EndIf
  Else
     Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("results/I_f%g.dat", Freq), LastTimeStepOnly];
     Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("results/U_f%g.dat", Freq), LastTimeStepOnly];
  EndIf

}
