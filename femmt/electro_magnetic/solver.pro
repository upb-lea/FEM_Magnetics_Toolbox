// Template solver files for FEMMT framework
// Most parts are taken from GetDP tuorials
// ----------------------
// ----------------------
// Functions

Jacobian {
  { Name Vol ; Case { { Region All ; Jacobian VolAxiSqu ; } } }
  { Name Sur ; Case { { Region All ; Jacobian SurAxi ; } } }
}

// ... that or almost that way for the non radial symmetrc case (z-symmetry)
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
          { GeoElement Triangle ;    NumberOfPoints 4 ; } //changed
          { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
          { GeoElement Line       ; NumberOfPoints  13 ; } //added
        } }
  } }
}

// ----------------------
// ----------------------

FunctionSpace {

  // Magnetic Vector Potential
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se1 ; NameOfCoef ae1 ; Function BF_PerpendicularEdge ;
        Support Region[{Domain}] ; Entity NodesOf [ All ] ; }
   }
    Constraint {
      { NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }

  // Gradient of Electric scalar potential (2D)
  { Name Hregion_u_2D ; Type Form1P ;
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

  // Imprinted Current Density
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
       //{ Name tmp  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
       //{ Name jc  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

       { Name ur ; Type Local  ; NameOfSpace Hregion_u_2D  ; }  // = nabla*el.potential = grad(phi)
       { Name I  ; Type Global ; NameOfSpace Hregion_u_2D[I] ; }
       { Name U  ; Type Global ; NameOfSpace Hregion_u_2D[U] ; }

       { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
       { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
       { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

       { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
       { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    /*
    Equation {
      Galerkin { [  Dof{tmp} , {tmp} ]; In Domain; Integration II; Jacobian Vol;  }
      Galerkin { [ nu[{d a}] * {d a} , {tmp} ]; In Domain; Integration II; Jacobian Vol;  }
    }

    Equation {
      Galerkin { [  Dof{jc} , {jc} ]; In Domain; Integration II; Jacobian Vol;  }
      Galerkin { [ sigma[]*(Dt[{a}]+{ur}/CoefGeo) , {jc} ]; In Domain; Integration II; Jacobian Vol;  }
    }
    */

    Equation {

      // Nabla x ( 1/mu Nabla x A)
      Galerkin { [ nu[Norm[{d a}], Freq] * Dof{d a} , {d a} ]  ;
        In Domain_Lin ; Jacobian Vol ; Integration II ; }
      If(Flag_NL)
        Galerkin { [ nu[{d a}, Freq] * Dof{d a} , {d a} ]  ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
        Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      // sigma d/dt A
      // technical current direction can be applied with neg. sigma (right hand rule)
      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      // sigma Nabla PHI
      Galerkin { [ -sigma[] * Dof{ur}/CoefGeo , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ -sigma[] * Dof{ur}/CoefGeo , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      // -Je (imprinted current density)
      Galerkin { [ -1/AreaCell[] *  Dof{ir}, {a} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ 1/AreaCell[] * Dof{a}, {ir} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      //GlobalTerm { [ Dof{Us}/CoefGeo, {Is} ] ; In DomainS ; }

      /*
      If(Flag_Conducting_Core)
        GlobalTerm { [ Dof{I}, {U} ] ; In Iron ; }
      EndIf
      */


      // GlobalTerms are used for the voltage - current relation
      For n In {1:n_windings}
          If(!Flag_HomogenisedModel~{n})
            If(Val_EE~{n}!=0)
              GlobalTerm { [ Dof{I}, {U} ] ; In Winding~{n} ; }
            EndIf
//          Else
//            Galerkin { [ NbrCond~{n}/CoefGeo/AreaCell[] / sigma[] * NbrCond~{n}/CoefGeo/AreaCell[]* Dof{ir} , {ir} ] ;
//                        In StrandedWinding~{n} ; Jacobian Vol ; Integration II ; }
//            GlobalTerm { [ Dof{Us}/CoefGeo, {Is} ] ; In StrandedWinding~{n} ; }
          EndIf
      EndFor











      //Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
      //  In DomainS ; Jacobian Vol ; Integration II ; }
      //Galerkin { DtDof [ CoefGeo*Ns[]/Sc[] * Dof{a}, {ir} ] ;
      //  In DomainS ; Jacobian Vol ; Integration II ; }

      //Galerkin { [ Ns[]/Sc[] / sigma[] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
      //  In DomainS ; Jacobian Vol ; Integration II ; }
      //GlobalTerm { [ Dof{Us}/CoefGeo , {Is} ] ; In DomainS ; }

      /* together
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
      */
    }
  }


}