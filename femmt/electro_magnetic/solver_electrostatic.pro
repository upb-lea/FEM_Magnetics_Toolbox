// Template solver file for FEMMT framework for Electrostatic Analysis
// Most parts are based on GetDP tutorials
// ----------------------
// ----------------------
// Functions

Jacobian {
  { Name Vol ; Case { { Region All ; Jacobian VolAxiSqu ; } } }
  { Name Sur ; Case { { Region All ; Jacobian SurAxi ; } } }
}

Integration {
  { Name II ; Case {
      { Type Gauss ; Case {
          { GeoElement Triangle ;    NumberOfPoints 4 ; }
          { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
          { GeoElement Line       ; NumberOfPoints  13 ; }
        } }
  } }
}

// ----------------------
// ----------------------
// FunctionSpace for Electrostatics
FunctionSpace {

  // Electric Scalar Potential
  { Name Hgrad_v_Ele ; Type Form0 ;
    BasisFunction {
//      { Name s0 ; NameOfCoef u0 ; Function BF_Node ;
//        Support Region[{Domain}] ; Entity NodesOf [ All ] ; }
      //For n In {1:n_windings}
            { Name sn; NameOfCoef vn; Function BF_Node;
                Support Region[{Domain}]; Entity NodesOf[ All, Not DomainC ]; }
      //EndFor
      //For n In {1:n_windings}
          { Name s0; NameOfCoef u0; Function BF_GroupOfNodes;
            Support Region[{Domain}]; Entity GroupsOfNodesOf[ DomainC ]; }
      //EndFor
          // Include additional basis functions for nodes within windings to show internal voltage
          /*For n In {1:n_windings}
            { Name sw; NameOfCoef vw; Function BF_Node;
              Support Region[{Winding~{n}}]; Entity NodesOf [ Winding~{n} ]; }
          EndFor*/
    }
    GlobalQuantity {
      { Name GlobalPotential; Type AliasOf       ; NameOfCoef u0; }
      { Name ArmatureCharge ; Type AssociatedWith; NameOfCoef u0; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_Ele; }
      { NameOfCoef GlobalPotential ; EntityType GroupsOfNodesOf ; NameOfConstraint Electrostatic_Potential ; }
      { NameOfCoef ArmatureCharge ; EntityType GroupsOfNodesOf ; NameOfConstraint SetArmatureCharge ; }
      // Dirichlet boundary condition for nodes inside the winding to reflect excitation potential
      /*For n In {1:n_windings}
        { NameOfCoef vw; EntityType NodesOf; NameOfConstraint Electrostatic_Potential; }
      EndFor*/

    }
  }
}
Formulation {
  { Name Electrostatic_Potential ; Type FemEquation ;
    Quantity {
       { Name u0  ; Type Local  ; NameOfSpace Hgrad_v_Ele ; }  // Electric potential
       { Name U   ; Type Global ; NameOfSpace Hgrad_v_Ele [GlobalPotential]; }
       { Name Q   ; Type Global; NameOfSpace Hgrad_v_Ele [ArmatureCharge]; }
    }
    Equation {
      // Poisson's Equation for Electrostatics: -div(epsilon * grad(u)) = 0
      /*Galerkin { [ epsilon[] * Dof{d u0}, {d u0} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }*/
      // Laplace's equation for electrostatics
//	  Galerkin { [ epsilon[{d u0}] * Dof{d u0}, {d u0} ] ;
//		In DomainCC ; Jacobian Vol ; Integration II ; }
      Integral { [ epsilon[] * Dof{d u0} , {d u0} ];
        In Domain; Jacobian Vol; Integration II; }
      //Integral { [ epsilon[] * Norm[{d u0}] , {u0} ]; In Sur_Neu_Ele; Jacobian Sur; Integration II; }
      //Integral { [ 0 , {u0} ]; In Sur_Neu_Ele; Jacobian Sur; Integration II; }
	  For n In {1:n_windings}
	    GlobalTerm { [ -Dof{Q} , {U} ]; In Winding~{n}; }
	  EndFor

    }
  }
}





