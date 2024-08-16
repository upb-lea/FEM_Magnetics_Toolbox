"""Includes functions to generate different winding topologies."""
# python libraries
import copy
from typing import Optional

# 3rd party libraries

# femmt libraries
from femmt.dtos import *
from femmt.model import WindingWindow, Conductor, Insulation, Core, VirtualWindingWindow
from femmt.functions_drawing import *
from femmt.functions_model import define_center_tapped_insulation


def create_stacked_winding_windows(core: Core, insulation: Insulation) -> (WindingWindow, WindingWindow):
    """
    Create stacked winding windows.

    :param core: Core class
    :type core: femmt.Core
    :param insulation: Insulation class
    :type insulation: Insulation
    :return: winding_window_top, winding_window_bot
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


def check_if_primary_conductor_row_fits_in_vww(vww: VirtualWindingWindow, row_element: ConductorRow,
                                               winding_element: Conductor,
                                               winding_insulations: ThreeWindingIsolation):
    """
    Check if the primary conductor row fits into the virtual winding window.

    :param vww: Virtual winding window
    :type vww: VirtualWindingWindow
    :param row_element: row element to be placed in vwws
    :type row_element: ConductorRow
    :param winding_element: Winding element as a conductor
    :type winding_element: Conductor
    :param winding_insulations: Winding insulations
    :type winding_insulations: ThreeWindingIsolation
    """
    if row_element.row_height > np.round(vww.top_bound-vww.bot_bound, 6):
        raise Exception("Row is too high!")
    elif row_element.number_of_conds_per_row * winding_element.conductor_radius*2 + (row_element.number_of_conds_per_row-1) * \
            winding_insulations.primary_to_primary >= vww.right_bound-vww.left_bound:
        # remark: additional bobbin is already included into the vww!
        raise Exception("Row does not fit into virtual winding window!")


def place_center_tapped_conductor_row(vwws: list, row_element: ConductorRow, row_winding_scheme_type: WindingScheme, no_vww: int,
                                      primary_conductors_to_be_placed: int, winding1: Conductor, winding2: Conductor, winding3: Conductor,
                                      winding_insulations: ThreeWindingIsolation, wrap_para_type: WrapParaType,
                                      foil_horizontal_placing_strategy: FoilHorizontalDistribution = None) -> int:
    """
    Place center-tapped conductor row.

    :param vwws: list of virtual winding windows
    :type vwws: list
    :param row_element: row element to be placed in vwws
    :type row_element: ConductorRow
    :param row_winding_scheme_type: Type of the winding scheme of the row
    :type row_winding_scheme_type: WindingScheme
    :param no_vww: Number of virtual winding window
    :type no_vww: int
    :param primary_conductors_to_be_placed: Number of primary conductors to be placed
    :type primary_conductors_to_be_placed: int
    :param winding1: Winding 1 as Conductor
    :type winding1: Conductor
    :param winding2: Winding 2 as Conductor
    :type winding2: Conductor
    :param winding3: Winding 3 as Conductor
    :type winding3: Conductor
    :param winding_insulations: Winding insulations
    :type winding_insulations: ThreeWindingIsolation
    :param wrap_para_type: wrap parameter type
    :type wrap_para_type: WrapParaType
    :param foil_horizontal_placing_strategy: strategy for placing foil horizontal windings
    :type foil_horizontal_placing_strategy: FoilHorizontalDistribution
    :return: Number of primary conductors to be placed
    :rtype: int
    """
    if row_element.winding_tag == WindingTag.Primary:
        check_if_primary_conductor_row_fits_in_vww(vww=vwws[no_vww], row_element=row_element, winding_element=winding1, winding_insulations=winding_insulations)
        primary_conductors_to_be_placed -= row_element.number_of_conds_per_row
        if primary_conductors_to_be_placed >= 0:
            vwws[no_vww].set_winding(winding1, row_element.number_of_conds_per_row, row_winding_scheme_type, wrap_para_type=wrap_para_type)
        elif primary_conductors_to_be_placed < 0:
            # In the last row, only the rest shall be placed
            vwws[no_vww].set_winding(winding1, row_element.number_of_conds_per_row + primary_conductors_to_be_placed, row_winding_scheme_type,
                                     wrap_para_type=wrap_para_type)
            primary_conductors_to_be_placed = 0

    elif row_element.winding_tag == WindingTag.Secondary:
        vwws[no_vww].set_winding(winding2, row_element.number_of_conds_per_row, row_winding_scheme_type,
                                 wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=foil_horizontal_placing_strategy)

    elif row_element.winding_tag == WindingTag.Tertiary:
        vwws[no_vww].set_winding(winding3, row_element.number_of_conds_per_row, row_winding_scheme_type,
                                 wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=foil_horizontal_placing_strategy)

    return primary_conductors_to_be_placed


def place_windings_in_vwws(vwws: list, winding_scheme_type: list, transformer_stack, primary_turns: int,
                           winding1: Conductor, winding2: Conductor, winding3: Conductor,
                           winding_insulations: ThreeWindingIsolation, wrap_para_type: WrapParaType,
                           foil_horizontal_placing_strategy: FoilHorizontalDistribution = None):
    """
    Place windings in virtual winding windows.

    :param vwws: list of virtual winding windows
    :type vwws: list
    :param winding_scheme_type: list with the winding schemes according to the transformer stack
                                !explicitly does not contain the insulations!
    :type winding_scheme_type: list
    :param transformer_stack:
    :type transformer_stack:
    :param primary_turns: Number of primary turns
    :type primary_turns: int
    :param winding1: Winding 1 as Conductor
    :type winding1: Conductor
    :param winding2: Winding 2 as Conductor
    :type winding2: Conductor
    :param winding3: Winding 3 as Conductor
    :type winding3: Conductor
    :param wrap_para_type: wrap parameter type
    :type wrap_para_type: WrapParaType
    :param winding_insulations: Winding insulations
    :type winding_insulations: ThreeWindingIsolation
    :param foil_horizontal_placing_strategy: strategy for placing foil horizontal windings
    :type foil_horizontal_placing_strategy: FoilHorizontalDistribution
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
        if isinstance(row_element, StackIsolation):
            pass

        elif isinstance(row_element, ConductorRow):
            primary_conductors_to_be_placed = place_center_tapped_conductor_row(vwws=vwws, row_element=row_element,
                                                                                row_winding_scheme_type=winding_scheme_type[set_vwws],
                                                                                no_vww=set_vwws,
                                                                                primary_conductors_to_be_placed=primary_conductors_to_be_placed,
                                                                                winding1=winding1, winding2=winding2, winding3=winding3,
                                                                                winding_insulations=winding_insulations,
                                                                                wrap_para_type=wrap_para_type,
                                                                                foil_horizontal_placing_strategy=foil_horizontal_placing_strategy)
            set_vwws += 1

        elif isinstance(row_element, CenterTappedGroup):
            if primary_conductors_to_be_placed < 0:
                turns1 = primary_conductors_to_be_placed
                primary_conductors_to_be_placed = 0
            else:
                turns1 = 0
            turns2 = 0

            for row in row_element.stack:
                if isinstance(row, StackIsolation):
                    pass
                elif isinstance(row, ConductorRow):
                    if row.winding_tag == WindingTag.Primary:
                        turns1 += row.number_of_conds_per_row
                    elif row.winding_tag == WindingTag.Secondary:
                        turns2 += row.number_of_conds_per_row

            vwws[set_vwws].set_center_tapped_winding(conductor1=winding1, turns1=turns1,
                                                     conductor2=winding2, turns2=turns2,
                                                     conductor3=winding3, turns3=turns2,
                                                     insulation_primary_to_primary=winding_insulations.primary_to_primary,
                                                     insulation_secondary_to_secondary=winding_insulations.secondary_to_secondary,
                                                     insulation_primary_to_secondary=winding_insulations.primary_to_secondary)

            set_vwws += 1

    return vwws


def set_center_tapped_windings(core: Core,
                               primary_turns: int, primary_radius: float, primary_number_strands: int, primary_strand_radius: float,
                               primary_additional_bobbin: float,
                               secondary_parallel_turns: int, secondary_thickness_foil: float, center_foil_additional_bobbin: float,
                               iso_top_core: float, iso_bot_core: float, iso_left_core: float, iso_right_core: float,
                               iso_primary_to_primary: float, iso_secondary_to_secondary: float, iso_primary_to_secondary: float,
                               interleaving_type: CenterTappedInterleavingType, interleaving_scheme: InterleavingSchemesFoilLitz,
                               bobbin_coil_top: float = None, bobbin_coil_bot: float = None, bobbin_coil_left: float = None, bobbin_coil_right: float = None,
                               primary_coil_turns: int = None, winding_temperature: Optional[float] = None,
                               wrap_para_type: WrapParaType = WrapParaType.FixedThickness,
                               foil_horizontal_placing_strategy: FoilHorizontalDistribution = None):
    """
    Set center tapped windings.

    :param interleaving_scheme: Interleaving scheme
    :type interleaving_scheme: InterleavingSchemesFoilLitz
    :param center_foil_additional_bobbin: Additional bobbin
    :type center_foil_additional_bobbin: float
    :param primary_additional_bobbin:
    :type primary_additional_bobbin: float
    :param interleaving_type: interleaving type
    :type interleaving_type: CenterTappedInterleavingType
    :param primary_strand_radius: Primary strand radius in m
    :type primary_strand_radius: float
    :param primary_number_strands: Primary number of strands
    :type primary_number_strands: in
    :param iso_primary_to_secondary: Insulation primary to secondary winding in m
    :type iso_primary_to_secondary: float
    :param iso_secondary_to_secondary: Insulation secondary to secondary winding in m
    :type iso_secondary_to_secondary: float
    :param iso_primary_to_primary: Insulation primary to primary winding in m
    :type iso_primary_to_primary: float
    :param iso_right_core: Insulation right core side in m
    :type iso_right_core: float
    :param iso_left_core: Insulation left core side in m
    :type iso_left_core: float
    :param iso_bot_core: Insulation bottom core side in m
    :type iso_bot_core: float
    :param iso_top_core: Insulation top core side in m
    :type iso_top_core: float
    :param core: Core class
    :type core: Core
    :param primary_turns: Number of primary turns
    :type primary_turns: int
    :param primary_radius: Radius of primary winding conductors in m
    :type primary_radius: float
    :param secondary_parallel_turns: Number of secondary parallel turns
    :type secondary_parallel_turns: int
    :param secondary_thickness_foil: Foil thickness of secondary winding in m
    :type secondary_thickness_foil: float
    :param bobbin_coil_right: Right bobbin thickness of the coil part in m
    :type bobbin_coil_right: float
    :param bobbin_coil_left: Left bobbin thickness of the coil part in m
    :type bobbin_coil_left: float
    :param bobbin_coil_bot: Bottom bobbin thickness of the coil part in m
    :type bobbin_coil_bot: float
    :param bobbin_coil_top: Top bobbin thickness of the coil part in m
    :type bobbin_coil_top: float
    :param primary_coil_turns: Number of primary coil turns
    :type primary_coil_turns: int
    :param winding_temperature: winding temperature in °C
    :type winding_temperature: Optional[float]
    :param wrap_para_type: wrap parameter type
    :type wrap_para_type: WrapParaType
    :param foil_horizontal_placing_strategy: strategy for placing foil horizontal windings
    :type foil_horizontal_placing_strategy: FoilHorizontalDistribution
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
        winding1.set_litz_round_conductor(primary_radius, primary_number_strands, primary_strand_radius, None,
                                          conductor_arrangement=ConductorArrangement.SquareFullWidth)
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
        vww_top = ww_top.split_window(WindingWindowSplit.NoSplitWithBobbin, top_bobbin=bobbin_coil_top, bot_bobbin=bobbin_coil_bot,
                                      left_bobbin=bobbin_coil_left, right_bobbin=bobbin_coil_right)
        available_height = core.window_h_bot - iso_top_core - iso_bot_core
    else:
        raise Exception(f"Unknown core type {core.core_type}")

    # Define the transformer winding stack
    transformer_stack = stack_center_tapped_transformer(primary_row, secondary_row, tertiary_row,
                                                        isolations=winding_insulations, available_height=available_height,
                                                        interleaving_type=interleaving_type, interleaving_scheme=interleaving_scheme,
                                                        primary_additional_bobbin=primary_additional_bobbin,
                                                        center_foil_additional_bobbin=center_foil_additional_bobbin)

    # Split the transformer winding window (ww_bot) in n virtual winding windows (vwws)
    vwws_bot, winding_scheme_type = ww_bot.split_with_stack(transformer_stack)

    # Place the windings in the virtual winding windows
    vwws_bot = place_windings_in_vwws(vwws_bot, winding_scheme_type, transformer_stack, primary_turns, winding1, winding2, winding3, winding_insulations,
                                      wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=foil_horizontal_placing_strategy)

    # If "stacked-core", also set primary coil turns
    if core.core_type == CoreType.Stacked:
        vww_top.set_winding(winding1, primary_coil_turns, None)

    # Depending on core geometry, return winding window(s) and insulation
    if core.core_type == CoreType.Single:
        return insulation, ww_bot
    elif core.core_type == CoreType.Stacked:
        return insulation, ww_top, ww_bot
