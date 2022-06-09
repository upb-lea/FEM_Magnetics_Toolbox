Include "Parameters.pro";
Include "Group.pro";
Include "Function.pro";
Include "Constraint.pro";
Include "Solver.pro";

Resolution {
  { Name analysis;
    System {
      { Name T; NameOfFormulation The_T; }
    }
    Operation {
      Generate[T] ; Solve T ; SaveSolution T; PostOperation[map] ;
	}
  }
}

PostProcessing {
  { Name The; NameOfFormulation The_T; NameOfSystem T;
    PostQuantity {
      { Name T; Value{ Local{ [ {T} ] ; In Total; Jacobian JVol; } } }
      { Name influx; Value{ Local{ [ {influx} ] ; In Warm; Jacobian JVol; } } }
	  { Name material; Value{ Local{ [ k[] ] ; In Total; Jacobian JVol; } } }
    }
  }
}

Include "PostOperation.pro";
