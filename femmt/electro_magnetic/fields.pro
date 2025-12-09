// === Field Quantities ===

PostOperation Map_local UsingPost MagDyn_a {

  // Potentials
  //Print[ raz,OnElementsOf Domain,  File StrCat[DirResFields, "raz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ az,OnElementsOf Domain,  File StrCat[DirResFields, "az", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ a,OnElementsOf Domain,  File StrCat[DirResFields, "a", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ur,OnElementsOf Domain,  File StrCat[DirResFields, "ur", ExtGmsh], LastTimeStepOnly ] ;

  // Electrical Field
  //Print[ e,  OnElementsOf Domain,  File StrCat[DirResFields, "e", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ MagEz,  OnElementsOf Domain,  File StrCat[DirResFields, "MagEz", ExtGmsh],  LastTimeStepOnly ] ;

  // Magnetic Field
  //Print[ h,  OnElementsOf Domain,  File StrCat[DirResFields, "h", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ Magh,  OnElementsOf Domain,  File StrCat[DirResFields, "Magh", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ Mag_h_real,  OnElementsOf Domain,  File StrCat[DirResFields, "Mag_h_real", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ Mag_h_imag,  OnElementsOf Domain,  File StrCat[DirResFields, "Mag_h_imag", ExtGmsh],  LastTimeStepOnly ] ;


  //Print[ mur,  OnElementsOf Domain,  File StrCat[DirResFields, "mur", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ mur_norm,  OnElementsOf Domain,  File StrCat[DirResFields, "mur_norm", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ nur_re,  OnElementsOf Domain,  File StrCat[DirResFields, "nur_re", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ nur_im,  OnElementsOf Domain,  File StrCat[DirResFields, "nur_im", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ mur_re,  OnElementsOf Domain,  File StrCat[DirResFields, "mur_re", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ mur_im,  OnElementsOf Domain,  File StrCat[DirResFields, "mur_im", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ sigma_re,  OnElementsOf DomainS,  File StrCat[DirResFields, "sigma_re", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ sigma_im,  OnElementsOf DomainS,  File StrCat[DirResFields, "sigma_im", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ p_hyst,  OnElementsOf Domain,  File StrCat[DirResFields, "p_hyst", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ p_hyst_density,  OnElementsOf Domain,  File StrCat[DirResFields, "p_hyst_density", ExtGmsh],  LastTimeStepOnly ] ;

  // Magnetic Flux (Density)
  //Print[ b,  OnElementsOf Domain,  File StrCat[DirResFields, "b", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ b_pol,  OnElementsOf Domain,  File StrCat[DirResFields, "b_pol", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ im_b_pol,  OnElementsOf Domain,  File StrCat[DirResFields, "im_b_pol", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ Mag_b_real,  OnElementsOf Domain,  File StrCat[DirResFields, "Mag_b_real", ExtGmsh],  LastTimeStepOnly] ;
  //Print[ Mag_b_imag,  OnElementsOf Domain,  File StrCat[DirResFields, "Mag_b_imag", ExtGmsh],  LastTimeStepOnly] ;
  If(Flag_show_standard_fields)
     Print[ Magb,  OnElementsOf Domain, Name "Magnitude B-Field / T" , File StrCat[DirResFields, "Magb", ExtGmsh]];
  EndIf
  //  , StoreInVariable $Magb maybe use this for Core Loss

  // Energy
  //Print[ MagEnergy,  OnElementsOf Domain,  File StrCat[DirResFields, "MagEnergy", ExtGmsh],  LastTimeStepOnly ] ;

  // Current Density
  //Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirResFields, "jz", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ j, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirResFields, "j", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ J_rms, OnElementsOf Region[{Domain}], File StrCat[DirResFields, "J_rms", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ir, OnElementsOf Region[{Domain}], File StrCat[DirResFields, "ir", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ir_re, OnElementsOf Region[{Domain}], File StrCat[DirResFields, "ir_re", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ir_im, OnElementsOf Region[{Domain}], File StrCat[DirResFields, "ir_im", ExtGmsh], LastTimeStepOnly ] ;
  //Print[ ir_norm, OnElementsOf Region[{Domain}], File StrCat[DirResFields, "ir_norm", ExtGmsh], LastTimeStepOnly ] ;

  // Ohmic Loss
  // This code is written to avoid the duplication of the printed losses, and also to avoid printing both losses where do we have just solid losses
  // This is just in stream visualization
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

    If(Flag_show_standard_fields)
        If(solid_exist)
            Print[ j2F_density, OnElementsOf Region[{DomainC}], Name "Solid wire and core eddy current loss density / W/m^3",
                   File StrCat[DirResFields, "j2F_density", ExtGmsh], LastTimeStepOnly] ;
        EndIf
        If(litz_exist)
            Print[ j2H_density, OnElementsOf DomainS, Name "Litz wire loss density / W/m^3", File StrCat[DirResFields,"j2H_density",ExtGmsh], LastTimeStepOnly] ;
        EndIf
    EndIf
 Else
    If(Flag_show_standard_fields)
         Print[ j2F_density, OnElementsOf Region[{DomainC}], Name "Solid wire and core eddy current loss density / W/m^3", File StrCat[DirResFields, "j2F_density", ExtGmsh]] ;
         Print[ j2H_density, OnElementsOf DomainS, Name "Litz wire loss density / W/m^3", File StrCat[DirResFields,"j2H_density",ExtGmsh]] ;
    EndIf
 EndIf

  // Settings
  /*
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].ScaleType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "Option.pos"]; */
 // Settings
 If (Flag_Stream_Visualization)
  Echo[ Str[
  "For k In {0:PostProcessing.NbViews-1}",
  "  View[k].RangeType  = 3;",  // per timestep
  "  View[k].NbIso = 25;",
  "  View[k].IntervalsType= 3;",
  "  View[k].AutoPosition = 3;",
  "  If (!StrCmp(View[k].Name, 'Solid wire and core eddy current loss density / W/m^3'))",
  "    View[k].ScaleType = 2;",
  "    View[k].SaturateValues = 1;",
  "  EndIf",
  "  If (!StrCmp(View[k].Name, 'Litz wire loss density / W/m^3'))",
  "    View[k].ScaleType = 2;",
  "  EndIf",
  "  If (!StrCmp(View[k].Name, 'Litz wire loss density / W/m^3'))",
  "    View[k].ScaleType = 2;",
  "  EndIf",
  "EndFor"
 ], File "option.pos"];
 EndIf
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

}