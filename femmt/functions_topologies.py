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


def place_windings(vwws, winding_scheme_type, transformer_stack, primary_turns,
                   winding1, winding2, winding3, winding_isolations):
    """

    :param vwws:
    :param winding_scheme_type:
    :param transformer_stack:
    :param primary_turns:
    :param winding1:
    :param winding2:
    :param winding3:
    :param winding_isolations:
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
            # TODO: kann man sicher viel eleganter lösen ...
            if row_element.winding_tag == WindingTag.Primary:
                primary_conductors_to_be_placed -= row_element.number_of_conds_per_row
                if primary_conductors_to_be_placed >= 0:
                    vwws[set_vwws].set_winding(winding1, row_element.number_of_conds_per_row, winding_scheme_type[set_vwws])
                elif primary_conductors_to_be_placed < 0:
                    # In the last row,only th rest shall be placed
                    vwws[set_vwws].set_winding(winding1, row_element.number_of_conds_per_row + primary_conductors_to_be_placed, winding_scheme_type[set_vwws])
                    primary_conductors_to_be_placed = 0
            elif row_element.winding_tag == WindingTag.Secondary:
                vwws[set_vwws].set_winding(winding2, row_element.number_of_conds_per_row, winding_scheme_type[set_vwws])
            elif row_element.winding_tag == WindingTag.Tertiary:
                vwws[set_vwws].set_winding(winding3, row_element.number_of_conds_per_row, winding_scheme_type[set_vwws])
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
                                                                isolation_primary_to_primary=winding_isolations.primary_to_primary,
                                                                isolation_secondary_to_secondary=winding_isolations.secondary_to_secondary,
                                                                isolation_primary_to_secondary=winding_isolations.primary_to_secondary)

            set_vwws += 1

    return vwws


def set_center_tapped_windings(core, primary_additional_bobbin,
                               primary_turns, primary_radius, primary_number_strands, primary_strand_radius,
                               secondary_parallel_turns, secondary_thickness_foil,
                               iso_top_core, iso_bot_core, iso_left_core, iso_right_core,
                               iso_primary_to_primary, iso_secondary_to_secondary, iso_primary_to_secondary,
                               interleaving_type: CenterTappedInterleavingType,
                               bobbin_coil_top=None, bobbin_coil_bot=None, bobbin_coil_left=None, bobbin_coil_right=None,
                               primary_coil_turns=None, winding_temperature: Optional[float] = None):
    """
    Set center tapped windings

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
        insulation = Insulation()
        insulation.add_core_insulations(iso_top_core, iso_bot_core, iso_left_core, iso_right_core)
        insulation.add_winding_insulations([[iso_primary_to_primary, iso_primary_to_secondary, iso_primary_to_secondary],
                                            [iso_primary_to_secondary, iso_secondary_to_secondary, iso_primary_to_secondary],
                                            [iso_primary_to_secondary, iso_primary_to_secondary, iso_secondary_to_secondary]])
        return insulation
    insulation = define_insulations()
    # TODO: the following statement does not provide any new information to the model at the moment -> MERGE both insulation concepts together (FEMMT globally)
    winding_isolations = define_center_tapped_insulation(primary_to_primary=iso_primary_to_primary,
                                                         secondary_to_secondary=iso_secondary_to_secondary,
                                                         primary_to_secondary=iso_primary_to_secondary)

    def define_windings(winding_temperature: float):
        winding1 = Conductor(0, Conductivity.Copper, winding_material_temperature=winding_temperature)
        winding1.set_litz_round_conductor(primary_radius, primary_number_strands, primary_strand_radius, None, conductor_arrangement=ConductorArrangement.SquareFullWidth)

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
    transformer_stack = stack_center_tapped_transformer(primary_row, secondary_row, tertiary_row, available_height=available_height, isolations=winding_isolations,
                                                        interleaving_type=interleaving_type, primary_additional_bobbin=primary_additional_bobbin)

    # Split the transformer winding window (ww_bot) in n virtual winding windows (vwws)
    vwws_bot, winding_scheme_type = ww_bot.split_with_stack(transformer_stack)

    # Place the windings in the virtual winding windows
    vwws_bot = place_windings(vwws_bot, winding_scheme_type, transformer_stack, primary_turns, winding1, winding2, winding3, winding_isolations)

    # If "stacked-core", also set primary coil turns
    if core.core_type == CoreType.Stacked:
        vww_top.set_winding(winding1, primary_coil_turns, None)

    # Depending on core geometry, return winding window(s) and insulation
    if core.core_type == CoreType.Single:
        return insulation, ww_bot
    elif core.core_type == CoreType.Stacked:
        return insulation, ww_top, ww_bot
