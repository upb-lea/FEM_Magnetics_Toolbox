// half inductor with axisymmetry
// 1 means full cylinder
SymFactor               = 1. ;
CoefGeo                 = 2*Pi*SymFactor ; // axisymmetry +/* symmetry factor */
n_windings = Number_of_Windings;

// ----------------------
// Physical numbers
// ----------------------
OUTBND              = 111111;
AIR                 = 110000;
AIR_EXT             = 110001;
AIR_COND            = 1000000;
CORE_PN             = 120000;


//physical numbers of conductors in n transformer
For n In {1:n_windings}
       iCOND~{n} = 130000 + 1000*(n-1);
       istrandedCOND~{n} = 150000 + 1000*(n-1);
EndFor

// ----------------------
// Groups
// ----------------------
Group{
  // ----------------------
  // Physical Domains
  // ----------------------
  Air  = Region[{AIR, AIR_EXT}];

  Core = Region[{}];
  // Core Domain
  For n In {1:nCoreParts}
    CorePart~{n} = Region[{(CORE_PN+n-1)}];
    Core += Region[{(CORE_PN+n-1)}];
  EndFor


  // Current Conducting Domains
  // Create a region for the winding
  For n In {1:n_windings} // loop over each winding
      Winding~{n} = Region[{}]; // create a region for the winding
      StrandedWinding~{n} = Region[{}]; // create a region for the stranded winding
  EndFor


  // Loop over each winding
  For winding_number In {1:n_windings}
      nbturns~{winding_number} = NbrCond~{winding_number} / SymFactor;
       // Loop over each turn in this winding to create a region for turns and then adding it to the winding
      For turn_number In {1:nbturns~{winding_number}}
        // Solid
        Turn~{winding_number}~{turn_number} = Region[{(iCOND~{winding_number}+turn_number-1)}];
        Winding~{winding_number} += Region[{(iCOND~{winding_number}+turn_number-1)}];
        Air += Region[{(AIR_COND+iCOND~{winding_number}+turn_number-1)}];
        // Stranded
        TurnStrand~{winding_number}~{turn_number} = Region[{(istrandedCOND~{winding_number}+turn_number-1)}];
        StrandedWinding~{winding_number} += Region[{(istrandedCOND~{winding_number}+turn_number-1)}];
        Air += Region[{(AIR_COND+istrandedCOND~{winding_number}+turn_number-1)}];
      EndFor
  EndFor


  // Non Conducting Domain:
  // Initialize the core-shell domain region to air
  DomainCC = Region[{Air}];

  // Add the Core region to the core-shell domain region
  If(!Flag_Conducting_Core)
    DomainCC += Region[{Core}];
  EndIf

  // Boundary Conditions
  // including symmetry
  OuterBoundary = Region[{OUTBND}];

   // Add the winding to the core domain region
  For n In {1:n_windings}
      DomainC += Region[{Winding~{n}}] ;
  EndFor
   // Add the Core region to the core domain region
  If(Flag_Conducting_Core)
    DomainC += Region[{Core}];
  EndIf
   // Add this stranded winding to the shell domain region
  For n In {1:n_windings}
      DomainS += Region[{StrandedWinding~{n}}] ;
  EndFor
  // Add the shell domain to the core-shell domain region
  DomainCC          += Region[{DomainS}] ;
   //  the linear and non linear domains to air and all windings

  If(Flag_NL)
    Domain_Lin      = Region[{Air}];
    For n In {1:n_windings}
        Domain_Lin += Region[{Winding~{n}, StrandedWinding~{n}}];
    EndFor
    Domain_Lin_NoJs = Region[{Air}];
    Domain_NonLin   = Region[{Core}];
  Else
    Domain_Lin      = Region[{Air, Core}];
    For n In {1:n_windings}
        Domain_Lin += Region[{Winding~{n}, StrandedWinding~{n}}];
    EndFor
    Domain_Lin_NoJs = Region[{Air, Core}];
    Domain_NonLin   = Region[{}];
  EndIf
  // Initialize the main domain to the core and core-shell domains
  Domain = Region[{DomainC, DomainCC}] ;

 // Loop over each winding and add its regions to the corresponding conductor domain
  For n In {1:n_windings}
      DomainCond~{n} += Region[{Winding~{n}, StrandedWinding~{n}}] ;
  EndFor

  // Dummy region number for postpro with functions
  DomainDummy = Region[ 12345 ] ;
}
