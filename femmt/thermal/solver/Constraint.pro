Constraint {
  { Name Temperature ;
	Case {
		If(flag_boundary_top==1)
			{ Region region_boundary_top ; Type Assign; Value value_boundary_top ; }
		EndIf
		If(flag_boundary_top_right==1)
			{ Region region_boundary_top_right ; Type Assign; Value value_boundary_top_right ; }
		EndIf
		If(flag_boundary_right_top==1)
			{ Region region_boundary_right_top ; Type Assign; Value value_boundary_right_top ; }
		EndIf
		If(flag_boundary_right==1)
			{ Region region_boundary_right ; Type Assign; Value value_boundary_right ; }
		EndIf
		If(flag_boundary_right_bottom==1)
			{ Region region_boundary_right_bottom ; Type Assign; Value value_boundary_right_bottom ; }
		EndIf
		If(flag_boundary_bottom_right==1)
			{ Region region_boundary_bottom_right ; Type Assign; Value value_boundary_bottom_right ; }
		EndIf
		If(flag_boundary_bottom==1)
			{ Region region_boundary_bottom ; Type Assign; Value value_boundary_bottom ; }
		EndIf
		}
	}
}