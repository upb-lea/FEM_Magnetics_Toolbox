Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian VolAxi ; }
    }
  }
  { Name JSur ;
    Case {
      { Region All ; Jacobian SurAxi ; }
    }
  }
}

Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
	  { GeoElement Point       ; NumberOfPoints  1 ; }
	  { GeoElement Line        ; NumberOfPoints  3 ; }
	  { GeoElement Triangle    ; NumberOfPoints  4 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	  { GeoElement Tetrahedron ; NumberOfPoints  4 ; }
	  { GeoElement Hexahedron  ; NumberOfPoints  6 ; }
	  { GeoElement Prism       ; NumberOfPoints  6 ; }
	}
      }
    }
  }
}

FunctionSpace {
  { Name Hgrad_T; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef Tn; Function BF_Node; Support Total;
        Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef Tn; EntityType NodesOf ; NameOfConstraint Temperature; }
    }
  }
	
  { Name Hgrad_T2; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef Tn; Function BF_Node; Support Warm;
        Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef Tn; EntityType NodesOf ; NameOfConstraint Temperature; }
    }
  }

}

Formulation {

  { Name The_T ; Type FemEquation;
    Quantity {
      { Name T;  Type Local; NameOfSpace Hgrad_T; }
	  { Name influx; Type Local; NameOfSpace Hgrad_T2; }
    }
    Equation {
      Galerkin { [ k[] * Dof{d T}, {d T} ];
                 In Total; Integration I1; Jacobian JVol;  }

      Galerkin { [ -qVol[] , {T} ];
                 In Total; Integration I1; Jacobian JVol;  }

      Integral { [ -qVol[] , {influx} ];
                 In Warm; Integration I1; Jacobian JVol;  }
				 
      Integral { [ Dof{influx} , {influx} ];
                 In Warm; Integration I1; Jacobian JVol;  }
    }
  }
}
