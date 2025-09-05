// === Global/Integrated Quantities ===


PostOperation Get_global UsingPost MagDyn_a {

  // Losses


  // Windings Total
  // Solid
  //Print[ SoF[ DomainC ], OnGlobal, Format TimeTable,  File > Sprintf("results/SF_iron.dat")] ; // TODO: Complex power
  //Print[ j2F[ Winding~{1} ], OnGlobal, Format TimeTable, File > StrCat[DirResVals, "j2F_1.dat"]] ;

 If (Flag_Stream_Visualization)
    // initialize
    solid_exist = 0;
    litz_exist  = 0;

    For n In {1:n_windings}
        If(!Flag_HomogenisedModel~{n})
            solid_exist = 1; // found a solid conductor
        Else
            litz_exist  = 1; // found a litz conductor
        EndIf
    EndFor
  If (solid_exist)
      For n In {1:n_windings}
          Print[ j2F[ Winding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals, "j2F_%g.dat"], n], LastTimeStepOnly, StoreInVariable $j2F, SendToServer Sprintf[StrCat[po,"14j2F_%g [J]"], n],  Color "LightYellow" ];
      EndFor
  EndIf
  If(litz_exist)
      For n In {1:n_windings}
          Print[ j2H[ StrandedWinding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"j2H_%g.dat"], n], LastTimeStepOnly, StoreInVariable $j2H, SendToServer Sprintf[StrCat[po,"15j2H_%g [J]"], n],  Color "LightYellow" ] ;
      EndFor
  EndIf

Else
    For n In {1:n_windings}
          Print[ j2F[ Winding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals, "j2F_%g.dat"], n], LastTimeStepOnly, StoreInVariable $j2F, SendToServer Sprintf[StrCat[po,"14j2F_%g [J]"], n],  Color "LightYellow" ];
           Print[ j2H[ StrandedWinding~{n} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"j2H_%g.dat"], n], LastTimeStepOnly, StoreInVariable $j2H, SendToServer Sprintf[StrCat[po,"15j2H_%g [J]"], n],  Color "LightYellow" ] ;
      EndFor
EndIf

  //Print[ j2F[ Winding3 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2F_3.dat"]] ;
  // Stranded
  //Print[ SoH[ StrandedWinding1 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"SH_1.dat"] ] ;  // TODO: Complex power
  //Print[ SoH[ DomainS ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"SH.dat"] ] ;
  //Print[ j2H[ DomainS ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"LossesStrandedWindings.dat"] ] ;
  //Print[ j2H[ StrandedWinding~{1} ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_1.dat"] ] ;
  //Print[ j2H[ StrandedWinding1 ], OnGlobal, Format Table];


  //Print[ j2H[ StrandedWinding2 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_2.dat"] ] ;
  //Print[ j2H[ StrandedWinding3 ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"j2H_3.dat"] ] ;
  //Print[ j2H[ StrandedWinding2 ], OnGlobal, Format Table];
  //Print[ j2Hskin[StrandedWinding1],   OnGlobal , Format Table];
  //Print[ j2Hprox[StrandedWinding1],   OnGlobal , Format Table];
  //Print[ j2Hskin[StrandedWinding2],   OnGlobal , Format Table];
  //Print[ j2Hprox[StrandedWinding2],   OnGlobal , Format Table];

  // Single Turns

   For n In {1:n_windings}
       If(Flag_HomogenisedModel~{n}) // Differentiate fine und hom
          For isF In {1:nbturns~{n}}
              Print[ j2H[ TurnStrand~{n}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{n},"Losses_turn_%g.dat"], isF], LastTimeStepOnly, StoreInVariable $j2H ] ;
          EndFor
       Else
          For isF In {1:nbturns~{n}}
              Print[ j2F[ Turn~{n}~{isF} ], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResValsWinding~{n},"Losses_turn_%g.dat"], isF], LastTimeStepOnly, StoreInVariable $j2F] ;
          EndFor
       EndIf
   EndFor


  // Core

  // Eddy Current Losses according to sigma in core/iron
  Print[ j2F[ Core ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"CoreEddyCurrentLosses.dat"], LastTimeStepOnly, StoreInVariable $j2F] ;
 // Hysteresis Losses according to complex permeability in Core
 // Hysteresis Losses according to complex permeability in core/iron (not true in time domain)
 // Hysteresis Losses according to complex permeability in Core
//  Print[ p_hyst[ Core ], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"p_hyst.dat"]] ;// Core losses
//  For n In {1:nCoreParts}
//      Print[ p_hyst[ CorePart~{n} ], OnGlobal , Format TimeTable, File > Sprintf[StrCat[DirResValsCore, "p_hyst_%g.dat"], n]] ;
//      Print[ j2F[ CorePart~{n} ], OnGlobal , Format TimeTable, File > Sprintf[StrCat[DirResValsCore, "CoreEddyCurrentLosses_%g.dat"], n]] ;
//  EndFor


  // Stored Energy
  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"ME.dat"], LastTimeStepOnly, StoreInVariable $MagEnergy, SendToServer StrCat[po,"12ME"],  Color "LightYellow"];

  //EndFor
   For n In {1:n_windings}
      Print[ Flux_Linkage~{n}[DomainCond~{n}], OnGlobal, Format Table, File > Sprintf[StrCat[DirResVals,"Flux_Linkage_%g.dat"], n], LastTimeStepOnly, StoreInVariable $Flux,
    SendToServer Sprintf[StrCat[po,"13Phi_%g [Wb]"], n],  Color "LightYellow"];
      //Print[ Flux_Linkage~{n}[DomainCond~{n}], OnGlobal, Format Table];
  EndFor

  // Inductances
  For n In {1:n_windings}
      If(Val_EE~{n}!=0)
         Print[ L~{n}~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"L_%g_%g.dat"], n, n], LastTimeStepOnly] ;
         Print[ LFromMagEnergy~{n}~{n}[Domain], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"LFromMagEnergy_%g_%g.dat"], n, n], LastTimeStepOnly] ;
      EndIf
  EndFor


  // Voltage
  //Print[ Voltage~{1}[DomainCond~{1}], OnGlobal, Format TimeTable, File > StrCat[DirResVals,"Voltage_1.dat"]];

  For n In {1:n_windings}
      Print[ Voltage~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"Voltage_%g.dat"], n], LastTimeStepOnly];
  EndFor


//  For n In {1:n_windings}
//      Print[ Current~{n}[DomainCond~{n}], OnGlobal, Format TimeTable, File > Sprintf[StrCat[DirResVals,"Current_%g.dat"], n], LastTimeStepOnly];
//  EndFor

  // Circuit Quantities



//  If(!Flag_Circuit)
//     For n In {1:n_windings}
//         If(!Flag_HomogenisedModel~{n})
//           Print[ I, OnRegion Winding~{n}, Format TimeTable, File > StrCat[DirResCirc,"I_%g_f%g.dat"] , LastTimeStepOnly, SendToServer StrCat[po,"20Fct_Src~{n}[] [A]"], Color "LightYellow" ];
//           Print[ U, OnRegion Winding~{n}, Format TimeTable, File > StrCat[DirResCirc,"U_%g_f%g.dat"] , LastTimeStepOnly, SendToServer StrCat[po,"20Voltage~{n}[] [A]"], Color "LightYellow" ];
//         Else
//           Print[ I, OnRegion StrandedWinding~{n}, Format TimeTable, File > StrCat[DirResCirc,"I_%g_f%g.dat"] , LastTimeStepOnly , SendToServer StrCat[po,"20Fct_Src~{n}[] [A]"], Color "LightYellow"];
//           Print[ U, OnRegion StrandedWinding~{n}, Format TimeTable, File > StrCat[DirResCirc,"U_%g_f%g.dat"] , LastTimeStepOnly, , SendToServer StrCat[po,"20Voltage~{n}[] [A]"], Color "LightYellow"];
//         EndIf
//     EndFor
//  Else
//     Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("results/I_f%g.dat", Freq), LastTimeStepOnly];
//     Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("results/U_f%g.dat", Freq), LastTimeStepOnly];
//  EndIf
  If(!Flag_Circuit)
     For n In {1:n_windings}
         If(!Flag_HomogenisedModel~{n})
           Print[ I, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_%g_T%g.dat"], n, T] , LastTimeStepOnly, SendToServer Sprintf[StrCat[po,"14I_%g [A]"], n], Color "LightYellow" ];
           Print[ U, OnRegion Winding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_%g_T%g.dat"], n, T] , LastTimeStepOnly, SendToServer Sprintf[StrCat[po,"15V_%g [V]"], n], Color "LightYellow"];
         Else
           Print[ I, OnRegion StrandedWinding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"I_%g_T%g.dat"], n, T] , LastTimeStepOnly, SendToServer Sprintf[StrCat[po,"14I_%g [A]"], n], Color "LightYellow"];
           Print[ U, OnRegion StrandedWinding~{n}, Format TimeTable, File > Sprintf[StrCat[DirResCirc,"U_%g_T%g.dat"], n, T] , LastTimeStepOnly, SendToServer Sprintf[StrCat[po,"15V_%g [V]"], n], Color "LightYellow"];
         EndIf
     EndFor
  Else
     Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("results/I_f%g.dat", Freq), LastTimeStepOnly];
     Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("results/U_f%g.dat", Freq), LastTimeStepOnly];
  EndIf

}
//PostOperation {
//  { Name T_resampled; NameOfPostProcessing The; ResampleTime[time0, timemax, dtResample];
//    Operation {
//      Print[ j2F, OnElementsOf Winding~{n} , File "j2F.dat"];
//    }
//  }
//}


