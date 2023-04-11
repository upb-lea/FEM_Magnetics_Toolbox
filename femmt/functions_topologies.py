from femmt import *
import copy


def create_stacked_winding_windows(core, insulation):
    """

    :param core:
    :param insulation:
    :return:
    """
    winding_window_top = WindingWindow(core, insulation)
    winding_window_bot = WindingWindow(core, insulation)

    max_pos = core.window_h_bot / 2 + core.core_inner_diameter / 4  # TODO: could also be done arbitrarily
    min_pos = core.window_h_bot / 2
    distance = max_pos - min_pos
    horizontal_split = min_pos + distance / 2
    insulation.vww_insulations = distance + 2 * min(insulation.core_cond)  # TODO: enhance the insulations situation!!!

    winding_window_top.max_bot_bound = horizontal_split + insulation.vww_insulations / 2
    winding_window_top.max_top_bound = winding_window_top.max_top_bound
    winding_window_top.max_left_bound = winding_window_top.max_left_bound
    winding_window_top.max_right_bound = winding_window_top.max_right_bound

    winding_window_bot.max_bot_bound = winding_window_bot.max_bot_bound
    winding_window_bot.max_top_bound = horizontal_split - insulation.vww_insulations / 2
    winding_window_bot.max_left_bound = winding_window_bot.max_left_bound
    winding_window_bot.max_right_bound = winding_window_bot.max_right_bound

    return winding_window_top, winding_window_bot


def place_windings_A(vwws, winding_scheme_type, transformer_stack, primary_turns,
                     winding1, winding2, winding3, winding_isolations):
    # . add conductor to vww and add winding window to MagneticComponent
    set_winding_counter = 0
    primary_turns_in_groups, secondary_turns_in_groups = get_number_of_turns_in_groups(transformer_stack)
    # print(primary_turns_in_groups, secondary_turns_in_groups)
    primary_conductors_to_be_placed = primary_turns - primary_turns_in_groups
    for row_element in transformer_stack.order:
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

    return vwws

def place_windings_B(vwws, winding_scheme_type, transformer_stack, primary_turns,
                     winding1, winding2, winding3, winding_isolations):


    # . add conductor to vww and add winding window to MagneticComponent
    set_winding_counter = 0
    primary_turns_in_groups, secondary_turns_in_groups = get_number_of_turns_in_groups(transformer_stack)
    # print(primary_turns_in_groups, secondary_turns_in_groups)
    primary_conductors_to_be_placed = primary_turns - primary_turns_in_groups
    for row_element in transformer_stack.order:
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

    return vwws


def set_center_tapped_windings(core, primary_turns, primary_radius, secondary_parallel_turns, secondary_thickness_foil, primary_coil_turns=None):
    """

    :param core:
    :param primary_turns:
    :param primary_radius:
    :param secondary_parallel_turns:
    :param secondary_thickness_foil:
    :param primary_coil_turns:
    :return:
    """
    # Define the insulation
    insulation = model.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    winding_isolations = define_center_tapped_insulation(primary_to_primary=2e-4,
                                                         secondary_to_secondary=2e-4,
                                                         primary_to_secondary=5e-4)
    insulation.add_winding_insulations([winding_isolations.primary_to_primary,
                                        winding_isolations.secondary_to_secondary,
                                        winding_isolations.primary_to_secondary], 0.0005)

    # Define windings
    winding1 = Conductor(0, Conductivity.Copper)
    winding1.set_litz_round_conductor(primary_radius, 50, 0.00011, None, conductor_arrangement=ConductorArrangement.SquareFullWidth)

    winding2 = Conductor(1, Conductivity.Copper)
    winding2.set_rectangular_conductor(thickness=secondary_thickness_foil)
    winding2.parallel = True

    winding3 = Conductor(2, Conductivity.Copper)
    winding3.set_rectangular_conductor(thickness=secondary_thickness_foil)
    winding3.parallel = True

    # Create single winding rows
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

    # Depending on core geometry, define the winding window
    if core.core_type == CoreType.Single:
        ww_bot = WindingWindow(core, insulation)
        ww_bot_height = core.window_h
    elif core.core_type == CoreType.Stacked:
        ww_top, ww_bot = fmt.create_stacked_winding_windows(core, insulation)
        vww_top = ww_top.split_window(fmt.WindingWindowSplit.NoSplit)
        ww_bot_height = core.window_h_bot

    # Define the transformer winding stack and split the transformer winding window
    interleaving_type = CenterTappedInterleavingType.TypeA
    transformer_stack = stack_center_tapped_transformer(primary_row, secondary_row, tertiary_row, window_height=ww_bot_height, isolations=winding_isolations,
                                                        interleaving_type=CenterTappedInterleavingType.TypeA)

    vwws_bot, winding_scheme_type = ww_bot.split_with_stack(transformer_stack)

    if interleaving_type == CenterTappedInterleavingType.TypeA:
        vwws_bot = place_windings_A(vwws_bot, winding_scheme_type, transformer_stack, primary_turns,
                                    winding1, winding2, winding3, winding_isolations)
    elif interleaving_type == CenterTappedInterleavingType.TypeB:
        vwws_bot = place_windings_B(vwws_bot, winding_scheme_type, transformer_stack, primary_turns,
                                    winding1, winding2, winding3, winding_isolations)


    # If "stacked-core", also set primary coil turns
    if core.core_type == CoreType.Stacked:
        vww_top.set_winding(winding1, primary_coil_turns, None)

    # Depending on core geometry, return virtual winding window(s) and insulation
    if core.core_type == CoreType.Single:
        return insulation, ww_bot
    elif core.core_type == CoreType.Stacked:
        return insulation, ww_top, ww_bot
