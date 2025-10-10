PostOperation Map_local UsingPost EleSta {

  // Potentials for the entire domain
  ExtGmsh     = ".pos";
  Print[ u0, OnElementsOf Region[{Domain}], Name "Potential / V", File StrCat[DirResFields, "Potential", ExtGmsh], LastTimeStepOnly ] ;
  If (!Flag_Stream_Visualization)
    Print[ u0, OnElementsOf Region[{Core}], Name "Average Potential on Core / V", File StrCat[DirResFields, "Voltage_Core_Map", ExtGmsh], LastTimeStepOnly ];
  EndIf

  // Electric Field vector in the entire domain
  // Print[ e, OnElementsOf Region[{Domain}], Name "Electric Field", File StrCat[DirResFields, "Efield", ExtGmsh], LastTimeStepOnly ] ;
  If (Flag_Stream_Visualization)
    Print[ Welocal, OnElementsOf Region[{Domain}], Name "Stored Energy", File StrCat[DirResFields, "We", ExtGmsh], LastTimeStepOnly ] ;
  EndIf
  Print[ MagE, OnElementsOf Region[{Domain}], Name "Magnitude Electric Field / V/m", File StrCat[DirResFields, "MagE", ExtGmsh], LastTimeStepOnly ] ;

  // Displacement Field vector in the entire domain
  // Print[ d, OnElementsOf Region[{Domain}], Name "Electric Field Density", File StrCat[DirResFields, "Dfield", ExtGmsh], LastTimeStepOnly ] ;
  Print[ MagD, OnElementsOf Region[{Domain}], Name "Magnitude Electric Field Density / C/m^2", File StrCat[DirResFields, "MagD", ExtGmsh], LastTimeStepOnly ] ;


  // Settings for visualization output (optional)
  // The electric field and electric field density are shown in linear
  If (Flag_Stream_Visualization)
  Echo[ Str[
  "For k In {0:PostProcessing.NbViews-1}",
  "  View[k].RangeType  = 3;",  // per timestep
  "  View[k].NbIso = 40;",
  "  View[k].IntervalsType= 3;",
  "  View[k].AutoPosition = 3;",
  "  If (!StrCmp(View[k].Name, 'Magnitude Electric Field / V/m'))",
  "    View[k].ScaleType = 1;",
  "  EndIf",
  "  If (!StrCmp(View[k].Name, 'Magnitude Electric Field Density / C/m^2'))",
  "    View[k].ScaleType = 1;",
  "  EndIf",
  "EndFor"
 ], File "option.pos"];
 EndIf

}