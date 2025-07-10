PostOperation Map_local UsingPost EleSta {

  // Potentials for the entire domain
  ExtGmsh     = ".pos";
  Print[ u0, OnElementsOf Region[{Domain}], Name "Potential / V", File StrCat[DirResFields, "Potential", ExtGmsh], LastTimeStepOnly ] ;
  Print[ u0, OnElementsOf Region[{Core}], Name "Average Potential on Core / V", File StrCat[DirResFields, "Voltage_Core_Map", ExtGmsh], LastTimeStepOnly ];

  // Electric Field vector in the entire domain
  Print[ e, OnElementsOf Region[{Domain}], Name "Electric Field", File StrCat[DirResFields, "Efield", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ Welocal, OnElementsOf Region[{Domain}], Name "Stored Energy", File StrCat[DirResFields, "We", ExtGmsh], LastTimeStepOnly ] ;
  Print[ MagE, OnElementsOf Region[{Domain}], Name "Magnitude Electric Field / V/m", File StrCat[DirResFields, "MagE", ExtGmsh], LastTimeStepOnly ] ;

  // Displacement Field vector in the entire domain
  Print[ d, OnElementsOf Region[{Domain}], Name "Electric Field Density", File StrCat[DirResFields, "Dfield", ExtGmsh], LastTimeStepOnly ] ;
  Print[ MagD, OnElementsOf Region[{Domain}], Name "Magnitude Electric Field Density / C/m^2", File StrCat[DirResFields, "MagD", ExtGmsh], LastTimeStepOnly ] ;


  // Settings for visualization output (optional)
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File OptionPos];
}