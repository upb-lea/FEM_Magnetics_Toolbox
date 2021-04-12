// GUI

mgeo = "{0Geometrical parameters/";
mfem = "{1Analysis parameters/";

col1 = "AliceBlue";
col2 = "Blue";
col3 = "Ivory";

DefineConstant[
  Flag_HomogenisedModel = {0, Choices{
      0="Fine reference", 1="Homogenized"}, Name StrCat[mfem,"001Model"], Highlight Str[col2] }

  Fine = !Flag_HomogenisedModel
  Homo = Flag_HomogenisedModel  // homogenized windings
];

DefineConstant[
  MD = {1,  Name StrCat[ mgeo,"01Mesh density"], Highlight Str[col1]}
  MDH = {2,  Name StrCat[ mgeo,"01Mesh density in winding window"], Highlight Str[col1]}
  Flag_HalfModel = {1,  Choices{0,1},  Name StrCat[ mgeo,"00Symmetry"], Highlight Str[col1], Closed 1}
  NbrLayersX = {6,  Min 1, Max 6, Name StrCat[ mgeo,"10Number of layers X"], Highlight Str[col3]}
  NbrLayersY = {20, Min 2, Max 20, Step 2, Name StrCat[ mgeo,"11Number of layers Y"],
    Highlight Str[col3]} // It has to be even
  NbrCond = {NbrLayersX * NbrLayersY, ReadOnly 1, Highlight "LightGrey",
    Name  StrCat[ mgeo,"12Number of turns (total)"]}

  Rc = { Sqrt[1/Pi]*1e-3, Name StrCat[ mgeo,"20Conductor radius"], Highlight Str[col3]}
  Dex = { 2.2*Rc, Name StrCat[ mgeo,"21Packing side"], Highlight Str[col3]}
  Dey = Dex
  AreaCell = Dex*Dey
  AreaCond = Pi*Rc^2
  Fill = {AreaCond/AreaCell, Name StrCat[ mgeo,"30Fill factor"], ReadOnly 1, Highlight "LightGrey"}

  sgap = {1., Min 1, Max 8,  Name StrCat[ mgeo,"40Gap size factor"], Highlight Str[col3]}
];

// Some dimensions
X1 = 0. ;
X2 = 10e-3;
X3 = 20e-3;
X4 = 25e-3;

Y2 = 19e-3;
Y3 = 13.5e-3;

Y4 = 1.5e-3*sgap ; //Half height of the gap

// coordinates of rectangular winding window
Xw1 = 11e-3 ;
Xw2 = Xw1 + 8.e-3 ;
Yw1 = -25.5e-3/2 ;
Yw2 =  25.5e-3/2 ;

If(Fmod[NbrLayersY,2])
  Printf("Warning: Number of layers along Y has to be even %g => %g", NbrLayersY, NbrLayersY+1);
  NbrLayersY = NbrLayersY+1;
EndIf

// Printf("Round conductor of radius %g mm",Rc*1e3);
// Printf("Fill factor = %g mm^2 / %g mm^2 = %g ", AreaCond*1e6, AreaCell*1e6, Fill);

SigmaCu = 6e7;
mu0 = 4.e-7 * Pi ;
nu0 = 1./mu0;

//DC Resistance R_dc = L*rho/A ; L = length (m); A = area (m^2)
//Inductance of a solenoid L = mu0*mur*NbrTurns*NbrTurns*Section/Length
NbrTurns = NbrCond ;
Len = (2*Pi*(Xw1+Xw2)/2)*NbrTurns ; // Length = reserved word?

gx = (Xw2-Xw1)-NbrLayersX * Dex ;
gy = (Yw2-Yw1)-NbrLayersY * Dey ;

// correcting of rectangular winding window coordinates (fine model=> a bit smaller)
Xw1_ = Xw1+gx/2 ; // radius of gap
Xw2_ = Xw2-gx/2 ;
Yw1_ = 0 ;
Yw2_ = Yw2-gy/2 ;

DefineConstant[
  // Some values computed directly from the geometry
  Rdc = { Len/SigmaCu/AreaCond, Name StrCat[ mgeo,"50Total DC resistance [Ω]"], ReadOnly 1, Highlight "LightGrey"}
  L   = { mu0*Pi*Xw1_^2*NbrTurns^2/(2*Y4), Name StrCat[ mgeo,"51Inductance [H m⁻¹]"], ReadOnly 1, Highlight "LightGrey"}
  tau = { L/Rdc, Name StrCat[ mgeo,"Time constant [s]"], ReadOnly 1, Highlight "LightGrey"}
  b_gap = { mu0*NbrTurns*1./(2*Y4), Name StrCat[ mgeo,"53Induction in airgap [T]"], ReadOnly 1, Highlight "LightGrey"}

  // Printf("DC resistance %g [Ohm]", Rdc);
  // Printf("Inductance  %g [H/m]", L);
  // Printf("Induction in airgap %g [T]", b_gap);
];

//----------------------------------
// Physical numbers
//----------------------------------
OUTBND  = 1111;

AIR    = 1000;
AIRGAP = 1100;
IRON   = 2000;
INSULATION = 3000;

iCOND = 4000;

ALLCOND = 5000;
