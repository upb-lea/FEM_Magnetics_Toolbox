from femmt import *
import copy


def set_center_tapped_windings(core, primary_turns, primary_radius, secondary_parallel_turns, secondary_thickness_foil):
    # following quantities may be refractored later
    winding_isolations = define_center_tapped_insulation(primary_to_primary=2e-4,
                                                         secondary_to_secondary=2e-4,
                                                         primary_to_secondary=5e-4)

    # TODO: sum following code in a function Center Tapped

    insulation = model.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([winding_isolations.primary_to_primary,
                                        winding_isolations.secondary_to_secondary,
                                        winding_isolations.primary_to_secondary], 0.0005)

    winding1 = Conductor(0, Conductivity.Copper)
    winding1.set_litz_round_conductor(primary_radius, 50, 0.00011, None, conductor_arrangement=ConductorArrangement.SquareFullWidth)

    winding2 = Conductor(1, Conductivity.Copper)
    winding2.set_rectangular_conductor(thickness=secondary_thickness_foil)
    winding2.parallel = True

    winding3 = Conductor(2, Conductivity.Copper)
    winding3.set_rectangular_conductor(thickness=secondary_thickness_foil)
    winding3.parallel = True

    primary_row = single_row(number_of_conds_per_winding=primary_turns,
                                 window_width=core.window_w - insulation.core_cond[2] - insulation.core_cond[3],
                                 winding_tag=WindingTag.Primary,
                                 conductor_type=ConductorType.RoundLitz,
                                 radius=primary_radius,
                                 cond_cond_isolation=winding_isolations.primary_to_primary)

    secondary_row = single_row(number_of_conds_per_winding=secondary_parallel_turns,
                                   window_width=core.window_w - insulation.core_cond[2] - insulation.core_cond[3],
                                   winding_tag=WindingTag.Secondary,
                                   conductor_type=ConductorType.RectangularSolid,
                                   thickness=secondary_thickness_foil)

    tertiary_row = copy.deepcopy(secondary_row)
    tertiary_row.winding_tag = WindingTag.Tertiary

    stack = stack_center_tapped_transformer(primary_row, secondary_row, tertiary_row,
                                                window_height=core.window_h, isolations=winding_isolations)
    # print(stack)
    winding_window = WindingWindow(core, insulation)

    vwws, winding_scheme_type = winding_window.split_with_stack(stack)
    # print(winding_scheme_type)

    # . add conductor to vww and add winding window to MagneticComponent
    set_winding_counter = 0
    primary_turns_in_groups, secondary_turns_in_groups = get_number_of_turns_in_groups(stack)
    # print(primary_turns_in_groups, secondary_turns_in_groups)
    primary_conductors_to_be_placed = primary_turns - primary_turns_in_groups
    for row_element in stack.order:
        if type(row_element) == StackIsolation:
            pass
        else:
            if type(row_element) == ConductorRow:
                # TODO: kann man sicher viel eleganter lösen ...
                if row_element.winding_tag == WindingTag.Primary:
                    primary_conductors_to_be_placed -= row_element.number_of_conds_per_row
                    if primary_conductors_to_be_placed >= 0:
                        vwws[set_winding_counter].set_winding(winding1,
                                                              row_element.number_of_conds_per_row,
                                                              winding_scheme_type[set_winding_counter])
                    elif primary_conductors_to_be_placed < 0:
                        # In the last row,only th rest shall be placed
                        vwws[set_winding_counter].set_winding(winding1,
                                                              row_element.number_of_conds_per_row + primary_conductors_to_be_placed,
                                                              winding_scheme_type[set_winding_counter])
                        primary_conductors_to_be_placed = 0
                elif row_element.winding_tag == WindingTag.Secondary:
                    vwws[set_winding_counter].set_winding(winding2,
                                                          row_element.number_of_conds_per_row,
                                                          winding_scheme_type[set_winding_counter])
                elif row_element.winding_tag == WindingTag.Tertiary:
                    vwws[set_winding_counter].set_winding(winding3,
                                                          row_element.number_of_conds_per_row,
                                                          winding_scheme_type[set_winding_counter])

            elif type(row_element) == CenterTappedGroup:
                if primary_conductors_to_be_placed < 0:
                    turns1 = primary_conductors_to_be_placed
                    primary_conductors_to_be_placed = 0
                else:
                    turns1 = 0
                turns2 = 0

                for row in row_element.stack:
                    if type(row) == StackIsolation:
                        pass
                    elif type(row) == ConductorRow:
                        if row.winding_tag == WindingTag.Primary:
                            turns1 += row.number_of_conds_per_row
                        elif row.winding_tag == WindingTag.Secondary:
                            turns2 += row.number_of_conds_per_row

                vwws[set_winding_counter].set_center_tapped_winding(conductor1=winding1, turns1=turns1,
                                                                    conductor2=winding2, turns2=turns2,
                                                                    conductor3=winding3, turns3=turns2,
                                                                    isolation_primary_to_primary=winding_isolations.primary_to_primary,
                                                                    isolation_secondary_to_secondary=winding_isolations.secondary_to_secondary,
                                                                    isolation_primary_to_secondary=winding_isolations.primary_to_secondary)

            set_winding_counter += 1

    # vwws = set_center_tapped_winding ...
    #
    # top_left.set_winding(winding1, 8, WindingType.Single)
    # top_right.set_winding(winding2, 6, WindingType.Single)
    # bot_right.set_winding(winding3, 12, WindingType.Single)
    # print(f"{insulation = }")
    return insulation, winding_window
