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
  Print[ SoH[ StrandedWinding1 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"SH_1.dat"] ] ;  // TODO

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
      Print[ az_int[ Turn1~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsPrimary,"a_turn_%g.dat"], isF] ] ;
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

  Print[ SoH[ DomainS ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"SH.dat"] ] ;

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

