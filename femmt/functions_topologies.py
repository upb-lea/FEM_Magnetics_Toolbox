# python libraries
import copy
from typing import Optional

# 3rd party libraries

# femmt libraries
from femmt.dtos import *
from femmt.model import WindingWindow, Conductor, Insulation
from femmt.functions_drawing import *
from femmt.functions_model import define_center_tapped_insulation


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


def check_if_primary_conductor_row_fits_in_vww(vww, row_element: ConductorRow, winding_element, winding_insulations):
    """

    :param vww:
    :param row_element:
    :return:
    """
    if row_element.row_height > np.round(vww.top_bound-vww.bot_bound, 6):
        raise Exception(f"Row is too high!")
    elif row_element.number_of_conds_per_row * winding_element.conductor_radius*2 + (row_element.number_of_conds_per_row-1)*winding_insulations.primary_to_primary >= vww.right_bound-vww.left_bound:
        # remark: additional bobbin is already included into the vww!
        raise Exception(f"Row does not fit into virtual winding window!")


def place_center_tapped_conductor_row(vwws, row_element, row_winding_scheme_type, no_vww, primary_conductors_to_be_placed,
                                      winding1, winding2, winding3, winding_insulations):
    """

    :param vwws: list of virtual winding windows
    :param row_element: row element to be placed in vwws
    :param row_winding_scheme_type:
    :param no_vww:
    :param primary_conductors_to_be_placed:
    :param winding1:
    :param winding2:
    :param winding3:
    :return:
    """
    if row_element.winding_tag == WindingTag.Primary:
        check_if_primary_conductor_row_fits_in_vww(vww=vwws[no_vww], row_element=row_element, winding_element=winding1, winding_insulations=winding_insulations)
        primary_conductors_to_be_placed -= row_element.number_of_conds_per_row
        if primary_conductors_to_be_placed >= 0:
            vwws[no_vww].set_winding(winding1, row_element.number_of_conds_per_row, row_winding_scheme_type)
        elif primary_conductors_to_be_placed < 0:
            # In the last row,only th rest shall be placed
            vwws[no_vww].set_winding(winding1, row_element.number_of_conds_per_row + primary_conductors_to_be_placed, row_winding_scheme_type)
            primary_conductors_to_be_placed = 0

    elif row_element.winding_tag == WindingTag.Secondary:
        vwws[no_vww].set_winding(winding2, row_element.number_of_conds_per_row, row_winding_scheme_type)

    elif row_element.winding_tag == WindingTag.Tertiary:
        vwws[no_vww].set_winding(winding3, row_element.number_of_conds_per_row, row_winding_scheme_type)

    return primary_conductors_to_be_placed


def place_windings_in_vwws(vwws, winding_scheme_type, transformer_stack, primary_turns,
                           winding1, winding2, winding3, winding_insulations):
    """

    :param vwws:
    :param winding_scheme_type: list with the winding schemes according to the transformer stack
                                !explicitly does not contain the insulations!
    :param transformer_stack:
    :param primary_turns:
    :param winding1:
    :param winding2:
    :param winding3:
    :param winding_insulations:
    :return:
    """
    # Count how many virtual winding windows were set
    set_vwws = 0

    # The number of turns in groups is needed to stop adding conductors, when all conductors are placed
    # Some conductor placing algorithms will overdefine the number, some will directly make it adequately
    primary_turns_in_groups, secondary_turns_in_groups = get_number_of_turns_in_groups(transformer_stack)
    primary_conductors_to_be_placed = primary_turns - primary_turns_in_groups

    # Iterate over the rows and place them
    for row_element in transformer_stack.order:
        if type(row_element) == StackIsolation:
            pass

        elif type(row_element) == ConductorRow:
            primary_conductors_to_be_placed = place_center_tapped_conductor_row(vwws=vwws, row_element=row_element, row_winding_scheme_type=winding_scheme_type[set_vwws],
                                                                                no_vww=set_vwws, primary_conductors_to_be_placed=primary_conductors_to_be_placed,
                                                                                winding1=winding1, winding2=winding2, winding3=winding3, winding_insulations=winding_insulations)
            set_vwws += 1

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

            vwws[set_vwws].set_center_tapped_winding(conductor1=winding1, turns1=turns1,
                                                     conductor2=winding2, turns2=turns2,
                                                     conductor3=winding3, turns3=turns2,
                                                     isolation_primary_to_primary=winding_insulations.primary_to_primary,
                                                     isolation_secondary_to_secondary=winding_insulations.secondary_to_secondary,
                                                     isolation_primary_to_secondary=winding_insulations.primary_to_secondary)

            set_vwws += 1

    return vwws


def set_center_tapped_windings(core,
                               primary_turns, primary_radius, primary_number_strands, primary_strand_radius, primary_additional_bobbin,
                               secondary_parallel_turns, secondary_thickness_foil, center_foil_additional_bobbin,
                               iso_top_core, iso_bot_core, iso_left_core, iso_right_core,
                               iso_primary_to_primary, iso_secondary_to_secondary, iso_primary_to_secondary,
                               interleaving_type: CenterTappedInterleavingType, interleaving_scheme: InterleavingSchemesFoilLitz,
                               bobbin_coil_top=None, bobbin_coil_bot=None, bobbin_coil_left=None, bobbin_coil_right=None,
                               primary_coil_turns=None, winding_temperature: Optional[float] = None):
    """
    Set center tapped windings

    :param interleaving_scheme:
    :param center_foil_additional_bobbin:
    :param primary_additional_bobbin:
    :param interleaving_type:
    :param primary_strand_radius:
    :param primary_number_strands:
    :param iso_primary_to_secondary:
    :param iso_secondary_to_secondary:
    :param iso_primary_to_primary:
    :param iso_right_core:
    :param iso_left_core:
    :param iso_bot_core:
    :param iso_top_core:
    :param core:
    :param primary_turns:
    :param primary_radius:
    :param secondary_parallel_turns:
    :param secondary_thickness_foil:
    :param bobbin_coil_right: 
    :param bobbin_coil_left: 
    :param bobbin_coil_bot: 
    :param bobbin_coil_top: 
    :param primary_coil_turns:
    :param winding_temperature: winding temperature in °C
    :type winding_temperature: Optional[float]
    :return:
    """
    def define_insulations():
        insulation = Insulation(flag_insulation=False)
        insulation.add_core_insulations(iso_top_core, iso_bot_core, iso_left_core, iso_right_core)
        insulation.add_winding_insulations([[iso_primary_to_primary, iso_primary_to_secondary, iso_primary_to_secondary],
                                            [iso_primary_to_secondary, iso_secondary_to_secondary, iso_primary_to_secondary],
                                            [iso_primary_to_secondary, iso_primary_to_secondary, iso_secondary_to_secondary]])
        return insulation
    insulation = define_insulations()
    # TODO: the following statement does not provide any new information to the model at the moment -> MERGE both insulation concepts together (FEMMT globally)
    winding_insulations = define_center_tapped_insulation(primary_to_primary=iso_primary_to_primary,
                                                          secondary_to_secondary=iso_secondary_to_secondary,
                                                          primary_to_secondary=iso_primary_to_secondary)

    def define_windings(winding_temperature: float):
        winding1 = Conductor(0, Conductivity.Copper, winding_material_temperature=winding_temperature)
        winding1.set_litz_round_conductor(primary_radius, primary_number_strands, primary_strand_radius, None, conductor_arrangement=ConductorArrangement.SquareFullWidth)
        # winding1.set_solid_round_conductor(primary_radius, conductor_arrangement=ConductorArrangement.SquareFullWidth)

        winding2 = Conductor(1, Conductivity.Copper, winding_material_temperature=winding_temperature)
        winding2.set_rectangular_conductor(thickness=secondary_thickness_foil)
        winding2.parallel = True

        winding3 = Conductor(2, Conductivity.Copper, winding_material_temperature=winding_temperature)
        winding3.set_rectangular_conductor(thickness=secondary_thickness_foil)
        winding3.parallel = True
        return winding1, winding2, winding3
    winding1, winding2, winding3 = define_windings(winding_temperature)

    def define_rows():
        primary_row = single_row(number_of_conds_per_winding=primary_turns,
                                 window_width=core.window_w - insulation.core_cond[2] - insulation.core_cond[3],
                                 winding_tag=WindingTag.Primary,
                                 conductor_type=ConductorType.RoundLitz,
                                 radius=primary_radius,
                                 cond_cond_isolation=insulation.cond_cond[0][0])

        secondary_row = single_row(number_of_conds_per_winding=secondary_parallel_turns,
                                   window_width=core.window_w - insulation.core_cond[2] - insulation.core_cond[3],
                                   winding_tag=WindingTag.Secondary,
                                   conductor_type=ConductorType.RectangularSolid,
                                   thickness=secondary_thickness_foil)

        tertiary_row = copy.deepcopy(secondary_row)
        tertiary_row.winding_tag = WindingTag.Tertiary
        return primary_row, secondary_row, tertiary_row
    primary_row, secondary_row, tertiary_row = define_rows()

    # Depending on core geometry, define the winding window
    if core.core_type == CoreType.Single:
        ww_bot = WindingWindow(core, insulation)
        available_height = core.window_h - iso_top_core - iso_bot_core
    elif core.core_type == CoreType.Stacked:
        ww_top, ww_bot = create_stacked_winding_windows(core, insulation)
        vww_top = ww_top.split_window(WindingWindowSplit.NoSplitWithBobbin, top_bobbin=bobbin_coil_top, bot_bobbin=bobbin_coil_bot, left_bobbin=bobbin_coil_left, right_bobbin=bobbin_coil_right)
        available_height = core.window_h_bot - iso_top_core - iso_bot_core
    else:
        raise Exception(f"Unknown core type {core.core_type}")

    # Define the transformer winding stack
    transformer_stack = stack_center_tapped_transformer(primary_row, secondary_row, tertiary_row,
                                                        isolations=winding_insulations, available_height=available_height,
                                                        interleaving_type=interleaving_type, interleaving_scheme=interleaving_scheme,
                                                        primary_additional_bobbin=primary_additional_bobbin, center_foil_additional_bobbin=center_foil_additional_bobbin)

    # Split the transformer winding window (ww_bot) in n virtual winding windows (vwws)
    vwws_bot, winding_scheme_type = ww_bot.split_with_stack(transformer_stack)

    # Place the windings in the virtual winding windows
    vwws_bot = place_windings_in_vwws(vwws_bot, winding_scheme_type, transformer_stack, primary_turns, winding1, winding2, winding3, winding_insulations)

    # If "stacked-core", also set primary coil turns
    if core.core_type == CoreType.Stacked:
        vww_top.set_winding(winding1, primary_coil_turns, None)

    # Depending on core geometry, return winding window(s) and insulation
    if core.core_type == CoreType.Single:
        return insulation, ww_bot
    elif core.core_type == CoreType.Stacked:
        return insulation, ww_top, ww_bot
