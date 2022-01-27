//getdp cell -setnumber Mode 1 -setnumber CellType 0 -sol MagDyn_a

Include "cell_dat.pro";

ExtGmsh    = ".pos" ;
ExtGnuplot = ".dat" ;
StringOut = Sprintf("_RS_la%.2g_%.2glayer.dat", Fill, NbrLayers); // Just two digits
ResDir = "coeff/";

xx_ = "2000"; // Active X (second pair corresponds to X')
yy_ = "0200"; // Active Y (second pair corresponds to Y')
null_="0000"; // Not active for that graph
plot_coefs = 0; // Choices{0,1}, Name "Output/0Plot skin and proximity coefficients ?"}


Function {
  //Sigma = 6e7;
  Sigma = 5.8e7;
  mu0 = 4.e-7 * Pi ;
  nu0 = 1./mu0;

  delta_cell = Rc/Rr_cell;
  Omega = 2/(delta_cell*delta_cell*mu0*Sigma);
  Freq = Omega/(2*Pi);  Period = 1./Freq;

  Printf("delta_cell      %g mm ", delta_cell*1000);
  Printf("Frequency  %g Hz   Pulsation  %g rad/s", Freq, Omega);
}

Group {
  TheCond = Region[{COND}];
  TheIsol = Region[{ISOL}];
  TheCell = Region[{TheCond, TheIsol}];

  Isol  = Region[{@ ISOL:ISOL+NbrCond-1 @}];
  Cond  = Region[{@ COND:COND+NbrCond-1 @}];
  Bound = Region[{BOUND}];

  DomainC  = Region[{Cond}];
  DomainCC = Region[{Isol}];
  Domain   = Region[{DomainC, DomainCC}] ;
  DomainDummy = Region[{1234}];
}

Function {
  nu[] = nu0;
  sigma[Cond] = Sigma;
  sigma[TheIsol] = 0;

  AreaTheCell[] = SurfaceArea[]{COND}+SurfaceArea[]{ISOL} ;
  AreaTheCond[] = SurfaceArea[]{COND} ;

  Rdc[] = 1. /AreaTheCond[] /Sigma ;
  Req_[] = Sqrt[AreaTheCond[]/Pi];

  SkinRef_i[] = mu0/(8*Pi) ;

  If(Omega==0)
    ProxRef_[] = AreaTheCell[]*Sigma*Fill*(Req_[]*1e-10)^2/4 ; // Added by Till to prevent from div by zero
  Else
    ProxRef_[] = AreaTheCell[]*Sigma*Fill*(Req_[]*Omega)^2/4 ;
  EndIf

  Constraint_a[] = (Mode==1) ? 0. :
    ((Mode==2) ? Vector[1,0,0] * Vector[$Y,-$X,0] : Vector[0,1,0] * Vector[$Y,-$X,0]) ;
  Constraint_I[] = (Mode==1) ? 1. :  0. ;
}

Constraint {
  { Name MagneticVectorPotential_2D ;
    Case {
      { Region Bound ; Value Constraint_a[] ; }
    }
  }

  { Name Current ;
    Case {
      { Region DomainC ; Value Constraint_I[] ; }
    }
  }

  { Name Voltage ;
    Case {
    }
  }

}

Jacobian {
  { Name Vol ; Case { { Region All ;  Jacobian Vol ; } } }
}
Integration {
  { Name CurlCurl ; Case { { Type Gauss ; Case { { GeoElement Triangle ; NumberOfPoints 1 ; } } } } }
}

FunctionSpace {
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef ae ; EntityType NodesOf ; NameOfConstraint MagneticVectorPotential_2D ; }
    }
  }

  { Name Hregion_u_2D ; Type Form1P ;
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC; Entity DomainC; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType Region ; NameOfConstraint Voltage ; }
      { NameOfCoef I ; EntityType Region ; NameOfConstraint Current ; }
    }
  }
}


Formulation {
 { Name MagDyn_a ; Type FemEquation ;
    Quantity {
       { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
       { Name ur ; Type Local  ; NameOfSpace Hregion_u_2D  ; }
       { Name I  ; Type Global ; NameOfSpace Hregion_u_2D  [I] ; }
       { Name U  ; Type Global ; NameOfSpace Hregion_u_2D  [U] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Domain ; Jacobian Vol ; Integration CurlCurl ; }
      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration CurlCurl ; }

      Galerkin { [ sigma[] * Dof{ur} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration CurlCurl ; }
      Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration CurlCurl ; }
      Galerkin { [ sigma[] * Dof{ur} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration CurlCurl ; }
      GlobalTerm { [ Dof{I} , {U} ] ; In DomainC ; }
    }
  }
}


Resolution {
  { Name MagDyn_a ;
    System {
      { Name A ; NameOfFormulation MagDyn_a ; Frequency Freq ; }
    }
    Operation {
      CreateDir[ResDir];

      SetTime[Rr_cell];
      Generate[A]; Solve[A]; SaveSolution[A];

      PostOperation[Map_local] ;
      PostOperation[Get_coeffs] ;
    }
  }
}


PostProcessing {
  { Name MagDyn_a ; NameOfFormulation MagDyn_a ;
    PostQuantity {
      { Name a ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol ;} } }
      { Name jz ; Value { Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}) ] ] ; In DomainC ; Jacobian Vol ;} } }	
      { Name jrms ; Value { Term { [ SquNorm[(Dt[{a}]+{ur})] ] ; In DomainC ; Jacobian Vol ;} } }	
      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name bav ; Value { Integral { [{d a}/AreaTheCell[]]  ;
            In Domain ; Jacobian Vol ; Integration CurlCurl ; } } }

      // Skin-effect coefficients
      { Name pI ; Value { Integral { [ SquNorm[(Sigma*(Dt[{a}]+{ur}))]*AreaTheCond[] ] ;
            In TheCell ; Jacobian Vol ; Integration CurlCurl ; } } }
      { Name qI ; Value { Integral { [ nu0 * SquNorm[{d a}]/SkinRef_i[] ] ;
            In TheCell ; Jacobian Vol ; Integration CurlCurl ; } } }

      // Proximity-effect coefficients
      { Name qB ; Value { Integral { [ SquNorm[{d a}]/AreaTheCell[] ] ;
            In TheCell ; Jacobian Vol  ; Integration CurlCurl ; } } }

      // -------------------------------------------------------
      // changed by Till to prevent from a zero at 0 Hz / needs a review
      // doesnt work at the moment
      // see FEMMT.py - pre_simulate()

      If(Omega==0)
           { Name pB ; Value { Term { [ Complex[1, 0]] ; In TheCond ; } } }

      Else
           { Name pB ; Value { Integral { [ SquNorm[sigma[]*(Dt[{a}]+{ur})]/Sigma/ProxRef_[] ] ;
             In TheCond ; Jacobian Vol  ; Integration CurlCurl ; } } }
      EndIf

      //{ Name pB ; Value { Integral { [ SquNorm[sigma[]*(Dt[{a}]+{ur})]/Sigma/ProxRef_[] ] ;
      //      In TheCond ; Jacobian Vol  ; Integration CurlCurl ; } } }
      // -------------------------------------------------------

      { Name nuRe ; Value { Term { Type Global; [ nu0*$qB ] ; In DomainDummy ; } } }
      { Name nuIm ; Value { Term { Type Global; [ nu0*$pB*Fill*Rr_cell^2/2 ] ; In DomainDummy ; } } }

      // Normalized by nu0 or mu0
      { Name muRe ; Value { Term { Type Global; [ Re[1/Complex[$qB, $pB*Fill*Rr_cell^2/2] ]] ; In DomainDummy ; } } }
  { Name muIm ; Value { Term { Type Global; [ Im[1/Complex[$qB, $pB*Fill*Rr_cell^2/2] ]] ; In DomainDummy ; } } }

    }
  }
}

PostOperation Map_local UsingPost MagDyn_a {
  //Print[ b, OnElementsOf Domain, File StrCat[ResDir, "b",ExtGmsh]] ;
  //Print[ a, OnElementsOf Domain, File StrCat[ResDir, "a",ExtGmsh]] ;
  //Print[ jz, OnElementsOf Domain, File StrCat[ResDir, "jz",ExtGmsh]] ;
  //Print[ jrms, OnElementsOf Domain, File StrCat[ResDir, "jrms",ExtGmsh]] ;
  Echo[ Str["i=PostProcessing.NbViews-1;
             View[i].Light=0;
             View[i].LineWidth = 2;
             View[i].RangeType=3;
             View[i].IntervalsType=3;
             View[i].NbIso = 25;"], File StrCat[ResDir,"option_cell.pos"] ];
}

PostOperation Get_coeffs UsingPost  MagDyn_a {
  If(Mode==1)
    Print[ pI[TheCond], OnGlobal, Format TimeTable, File > StrCat[ResDir, "pI", StringOut],
      SendToServer "Output/0Skin/pI", StoreInVariable $skin_pI ] ;
    Print[ qI[TheCell], OnGlobal, Format TimeTable, File > StrCat[ResDir, "qI", StringOut],
      SendToServer "Output/0Skin/qI", StoreInVariable $skin_qI ] ;
  EndIf
  If(Mode==2)
    Print[ qB[TheCell], OnGlobal, Format TimeTable, File > StrCat[ResDir, "qB", StringOut],
      SendToServer "Output/1Prox/qB", StoreInVariable $qB] ;
    Print[ pB[TheCond], OnGlobal, Format TimeTable, File > StrCat[ResDir, "pB", StringOut],
      SendToServer "Output/1Prox/pB", StoreInVariable $pB ] ;
    Print[ muRe, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
      SendToServer "Output/2Prox/Re(mu)", File > StrCat[ResDir, "muRe", StringOut]];
    Print[ muIm, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
      SendToServer "Output/2Prox/Im(mu)", File > StrCat[ResDir, "muIm", StringOut]];
  EndIf
  If(Mode==3)
    Print[ qB[TheCell], OnGlobal, Format TimeTable, File > StrCat[ResDir, "qBy", StringOut],
      SendToServer "Output/Prox/qBy", StoreInVariable $qBy ] ;
    Print[ pB[TheCond], OnGlobal, Format TimeTable, File > StrCat[ResDir, "pBy", StringOut],
      SendToServer "Output/Prox/pBy", StoreInVariable $pBy ] ;
  EndIf
}


DefineConstant[
  R_ = {"MagDyn_a", Name "GetDP/1ResolutionChoices", Visible 1, Closed 1},
  C_ = {"-solve -v 4 -v2 -bin", Name "GetDP/9ComputeCommand", Visible 1, Closed 1}
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 1,Closed 1}
];


If(plot_coefs)
  DefineConstant[
    pI_= {1, Name "Output/0Skin/pI", Graph StrCat[yy_,null_,null_,null_], Visible (Mode==1)}
    qI_= {1, Name "Output/0Skin/qI", Graph StrCat[null_,yy_,null_,null_], Visible (Mode==1)}
    qB_= {1, Name "Output/1Prox/qB", Graph StrCat[null_,null_,yy_,null_], Visible (Mode==2)}
    pB_= {1, Name "Output/1Prox/pB", Graph StrCat[null_,null_,null_,yy_], Visible (Mode==2)}
  ];
EndIf
