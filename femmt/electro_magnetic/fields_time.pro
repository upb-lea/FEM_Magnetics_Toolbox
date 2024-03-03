

// === Field Quantities ===

PostOperation Map_local UsingPost MagDyn_a {

  // Potentials
  //Print[ raz,OnElementsOf Domain,  File StrCat[DirResFields, "raz", ExtGmsh], LastTimeStepOnly ] ;
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

  // Core Loss Density
  If(Flag_Generalized_Steinmetz_loss)
    Print[ piGSE,  OnElementsOf Domain,  File StrCat[DirResFields, "piGSE", ExtGmsh],  LastTimeStepOnly ] ;
  EndIf

  If(Flag_Steinmetz_loss)
    Print[ pSE,  OnElementsOf Domain,  File StrCat[DirResFields, "pSE", ExtGmsh],  LastTimeStepOnly ] ;
    Print[ pSE_density,  OnElementsOf Domain,  File StrCat[DirResFields, "pSE_density", ExtGmsh],  LastTimeStepOnly ] ;
  EndIf

  //Print[ mur,  OnElementsOf Domain,  File StrCat[DirResFields, "mur", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ mur_norm,  OnElementsOf Domain,  File StrCat[DirResFields, "mur_norm", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ nur_re,  OnElementsOf Domain,  File StrCat[DirResFields, "nur_re", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ nur_im,  OnElementsOf Domain,  File StrCat[DirResFields, "nur_im", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ mur_re,  OnElementsOf Domain,  File StrCat[DirResFields, "mur_re", ExtGmsh],  LastTimeStepOnly ] ;
  //Print[ mur_im,  OnElementsOf Domain,  File StrCat[DirResFields, "mur_im", ExtGmsh],  LastTimeStepOnly ] ;
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
     //Print[ Magb,  OnElementsOf Domain, Name "Magnitude B-Field / T" , File StrCat[DirResFields, "Magb", ExtGmsh], LastTimeStepOnly];  // for  compute command -v2
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
  If(Flag_show_standard_fields)
    // to show the losses in every step
    Print[ j2F_density, OnElementsOf Region[{DomainC}], Name "Solid wire and core eddy current loss density / W/m^3", File StrCat[DirResFields, "j2F_density", ExtGmsh]] ;
    // it can be like this
    //Print[ j2F_density, OnElementsOf Region[{DomainC}], Name "Solid wire and core eddy current loss density / W/m^3", File StrCat[DirResFields, "j2F_density", ExtGmsh, LastTimeStepOnly ] ;
  EndIf
  If(Flag_show_standard_fields)
    //Print[ j2H, OnElementsOf DomainS, Name "Litz wire losses / W" , File StrCat[DirResFields,"jH",ExtGmsh] ] ;
    Print[ j2H_density, OnElementsOf DomainS, Name "Litz wire loss density / W/m^3", File StrCat[DirResFields,"j2H_density",ExtGmsh]] ;
  EndIf
  //Print[ j2Hprox,   OnElementsOf DomainS, File StrCat[DirResFields,"jHprox",ExtGmsh] ] ;
  //Print[ j2Hskin,   OnElementsOf DomainS, File StrCat[DirResFields,"jHskin",ExtGmsh] ] ;


  // Settings
  Echo[Str[ "For k In {0:PostProcessing.NbViews-1}",
      "View[k].RangeType = 3;" ,// per timestep
      "View[k].NbIso = 25;",
      "View[k].IntervalsType = 3;",
      "View[k].AutoPosition = 3",
      "EndFor"// iso values
    ], File "option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

}
