PostOperation Map_local UsingPost EleSta {

  // Potentials for the entire domain
  ExtGmsh     = ".pos";
  //Print[ u0, OnElementsOf Domain, File StrCat[DirResFields, "Potential_Map.pos"], LastTimeStepOnly, Format Gmsh];
  //Print[ u0,  OnElementsOf Domain, File StrCat[DirResFields,"Potential",ExtGmsh], LastTimeStepOnly ] ;
  //Print[ u0,  OnElementsOf Domain, Name "Potential / V" , File StrCat[DirResFields, "Potential", ExtGmsh]];
  Print[ u0, OnElementsOf Region[{Domain}], Name "Potential / V", File StrCat[DirResFields, "Potential", ExtGmsh], LastTimeStepOnly ] ;
//  Print [ u0, OnLine {{0.0075,-0.01475,0}{0.0075,-0.00025,0}} {125}, Name "Potential on Core BottomLeft / V",File StrCat[DirResFields, "Potential_BotLeftCore", ExtGmsh], LastTimeStepOnly ] ;
//  Print [ u0, OnLine {{0.0075,0.00025,0}{0.0075,0.01475,0}} {125}, Name "Potential on Core TopLeft / V",File StrCat[DirResFields, "Potential_TopLeftCore", ExtGmsh], LastTimeStepOnly ] ;
//  Print [ u0, OnLine {{0.0075,0.01475,0}{0.0195,0.01475,0}} {125}, Name "Potential on Core Top / V",File StrCat[DirResFields, "TopCore", ExtGmsh], LastTimeStepOnly ] ;
//  Print [ u0, OnLine {{0.0195,0.01475,0}{0.0195,-0.01475,0}} {125}, Name "Potential on Core Right / V",File StrCat[DirResFields, "RightCore", ExtGmsh], LastTimeStepOnly ] ;
//  Print [ u0, OnLine {{0.0195,-0.01475,0}{0.0075,-0.01475,0}} {125}, Name "Potential on Core Bot / V",File StrCat[DirResFields, "BotCore", ExtGmsh], LastTimeStepOnly ] ;
//  Print [ u0, OnPoint {0.0195,-0.01475,0}, Name "Potential on Point CoreBotRight / V",File StrCat[DirResFields, "CorePoint", ExtGmsh], LastTimeStepOnly ];
//  Print [ u0, OnPoint {0.0154,-0.01245,0}, Name "Potential on PointLastConductorMid / V",File StrCat[DirResFields, "TurnPointMid", ExtGmsh], LastTimeStepOnly ];
//  Print [ u0, OnPoint {0.0167,-0.01245,0}, Name "Potential on PointLastConductorRight/ V",File StrCat[DirResFields, "TurnPointRight", ExtGmsh], LastTimeStepOnly ];
  Print[ u0, OnElementsOf Region[{Core}], Name "Voltage on Core / V", File StrCat[DirResFields, "Voltage_Core_Map", ExtGmsh], LastTimeStepOnly ];

  // Electric Field vector in the entire domain
  // Print[ e, OnElementsOf Domain, File StrCat[DirResFields, "Electric_Field_Map.pos"], LastTimeStepOnly, Format Gmsh ];
  //Print[ e,  OnElementsOf Domain, File StrCat[DirResFields,"Efield",ExtGmsh], LastTimeStepOnly ] ;
  Print[ e, OnElementsOf Region[{Domain}], Name "Electric Field", File StrCat[DirResFields, "Efield", ExtGmsh], LastTimeStepOnly ] ;
  Print[ Welocal, OnElementsOf Region[{Domain}], Name "Stored Energy", File StrCat[DirResFields, "We", ExtGmsh], LastTimeStepOnly ] ;
  Print[ MagE, OnElementsOf Region[{Domain}], Name "Magnitude Electric Field / V/m", File StrCat[DirResFields, "MagE", ExtGmsh], LastTimeStepOnly ] ;

  // Displacement Field vector in the entire domain

  //Print[ d, OnElementsOf Domain, File StrCat[DirResFields, "Displacement_Field_Map.pos"], LastTimeStepOnly, Format Gmsh];
  //Print[ d,  OnElementsOf Domain, File StrCat[DirResFields,"Dfield",ExtGmsh], LastTimeStepOnly ] ;
  Print[ d, OnElementsOf Region[{Domain}], Name "Electric Field Density", File StrCat[DirResFields, "Dfield", ExtGmsh], LastTimeStepOnly ] ;
  Print[ MagD, OnElementsOf Region[{Domain}], Name "Magnitude Electric Field Density / C/m^2", File StrCat[DirResFields, "MagD", ExtGmsh], LastTimeStepOnly ] ;
  Print[ Q, OnElementsOf Region[{DomainC}], Name "Charge/ C", File StrCat[DirResFields, "Q", ExtGmsh], LastTimeStepOnly ] ;


  // Settings for visualization output (optional)
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File OptionPos];
}