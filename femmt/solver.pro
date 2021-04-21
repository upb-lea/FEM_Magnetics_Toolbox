// Template solver files for FEMMT framework
// Most parts are taken from GetDP tuorials
// ----------------------
// ----------------------
// Functions

Jacobian {
  { Name Vol ; Case { { Region All ; Jacobian VolAxiSqu ; } } }
  { Name Sur ; Case { { Region All ; Jacobian SurAxi ; } } }
}

// ... so oder so ähnlich für den nicht-axialsymmetrischen Fall
//Jacobian {
//  { Name Vol;
//    Case {
//      { Region DomainInf ; Jacobian VolSphShell{Val_Rint, Val_Rext} ; }
//      { Region All ; Jacobian Vol; }
//    }
//  }
//}

Integration {
  { Name II ; Case {
      { Type Gauss ; Case {
          { GeoElement Triangle ;    NumberOfPoints 4 ; }
          { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
        } }
    } }
}

// ----------------------
// ----------------------

FunctionSpace {

  { Name Hcurl_a_2D ; Type Form1P ; // Split for TD + homog: subspace isolating unknowns
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ All, Not Winding ] ; }
      { Name seh ; NameOfCoef aeh ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ Winding ] ; }
    }
    SubSpace {
      { Name aH ; NameOfBasisFunction {seh}; } // Subspace only used in TD homog
    }

    Constraint {
      { NameOfCoef ae  ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
     }
  }

  { Name Hregion_u_2D ; Type Form1P ; // Gradient of Electric scalar potential (2D)
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC ; Entity DomainC ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef I ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS ; Entity DomainS ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


  // For circuit equations
  { Name Hregion_Z ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainZt_Cir ; Entity DomainZt_Cir ; }
    }
    GlobalQuantity {
      { Name Iz ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Uz ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }

}

// ----------------------
// ----------------------
// formulation of the magnetodynmic problem in terms of a magnetic vector potential
// frequency domain
Formulation {
  { Name MagDyn_a ; Type FemEquation ;
    Quantity {
       { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

       { Name ur ; Type Local  ; NameOfSpace Hregion_u_2D  ; }  // = nabla*el.potential = grad*phi
       { Name I  ; Type Global ; NameOfSpace Hregion_u_2D[I] ; }
       { Name U  ; Type Global ; NameOfSpace Hregion_u_2D[U] ; }

       { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
       { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
       { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

       { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
       { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Domain_Lin ; Jacobian Vol ; Integration II ; }
      If(Flag_NL)
        Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
        Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{ur}/CoefGeo , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{ur}/CoefGeo , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{I}, {U} ] ; In DomainC ; }

      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ CoefGeo*Ns[]/Sc[] * Dof{a}, {ir} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }

      Galerkin { [ Ns[]/Sc[] / sigma[] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
        In DomainS ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{Us}/CoefGeo , {Is} ] ; In DomainS ; }

      If(Flag_Circuit)
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }

        GlobalTerm { [ Dof{Uz}                      , {Iz} ] ; In Inductance_Cir ; }
        GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ] ; In Inductance_Cir ; }

        GlobalTerm { [ Dof{Iz}        , {Iz} ] ; In Capacitance1_Cir ; }
        GlobalTerm { NeverDt[ Dof{Iz} , {Iz} ] ; In Capacitance2_Cir ; }
        GlobalTerm { DtDof [ Capacitance[] * Dof{Uz} , {Iz} ] ; In Capacitance_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {I};  Loop {U};  Equation {I};  In DomainC ; }
          { Node {Is}; Loop {Us}; Equation {Us}; In DomainS ; }
          { Node {Iz}; Loop {Uz}; Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }
}