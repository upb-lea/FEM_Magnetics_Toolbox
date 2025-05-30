"""Functions to draw different conductor schemes for the FEM simulation."""
# python libraries
import copy

# 3rd party libraries
import numpy as np

# femmt libraries
from femmt import *
from femmt.dtos import *


def number_of_rows(row: ConductorRow) -> int:
    """
    Get the number of rows needed to place some conductors in a given winding window.

    :param row: ConductorRow DTO
    :type row: ConductorRow
    :return: number of rows
    :rtype: int
    """
    number_rows = int(row.number_of_conds_per_winding / row.number_of_conds_per_row) + \
                     (row.number_of_conds_per_winding % row.number_of_conds_per_row > 0)
    return number_rows


def single_row(number_of_conds_per_winding: int, window_width: float, winding_tag: WindingTag, conductor_type: ConductorType,
               thickness: float | None = None, radius: float | None = None, cond_cond_isolation=None,
               additional_bobbin: float = 0) -> ConductorRow:
    """
    Define a full row of the defined conductor in the specified window width.

    It is assumed, that the row is full, which means, as many conductors are
    place horizontally next to each other, as possible.
    It is assumed, that the row fits into the winding window. A later control
    is needed.

    :param number_of_conds_per_winding: number of conductors per winding
    :type number_of_conds_per_winding: int
    :param window_width: window width in meter
    :type window_width: float
    :param winding_tag: winding tag (Primary / Secondary / Tertiary)
    :type winding_tag: WindingTag
    :param conductor_type: conductor type (RoundSolid / RoundLitz / RectangularSolid)
    :type conductor_type: ConductorType
    :param thickness: Optional[float]
    :type thickness: Optional[float]
    :param radius: radius in meter
    :type radius: Optional[float]
    :param cond_cond_isolation: distance between two conductors in meter
    :type cond_cond_isolation: float
    :param additional_bobbin: extra bobbin inner distance in meter
    :type additional_bobbin: float
    :return:
    """
    # Create ConductorRow dataclass
    conductor_row = ConductorRow(number_of_conds_per_winding=number_of_conds_per_winding,
                                 number_of_conds_per_row=None,
                                 row_height=None,
                                 winding_tag=winding_tag,
                                 number_of_rows=None,
                                 additional_bobbin=additional_bobbin)

    if conductor_type == ConductorType.RoundLitz or conductor_type == ConductorType.RoundSolid:
        conductor_row.row_height = radius * 2
        # How many full conductor cross-sections fit into one row of the winding window?
        conductor_row.number_of_conds_per_row = 0
        while window_width > 2 * (conductor_row.number_of_conds_per_row + 1) * radius + cond_cond_isolation * \
                conductor_row.number_of_conds_per_row:
            conductor_row.number_of_conds_per_row += 1
    elif conductor_type == ConductorType.RectangularSolid:
        # Rectangular Solid Conductors are meant to be as possible, so their cross-section is
        # simply adjusted to the winding window width
        conductor_row.row_height = thickness
        conductor_row.number_of_conds_per_row = 1

    conductor_row.number_of_rows = number_of_rows(conductor_row)
    return conductor_row


def check_secondary_and_tertiary_are_the_same(secondary_row: ConductorRow, tertiary_row: ConductorRow):
    """
    Define the single rows relate to a center-tapped transformer.

    :param secondary_row: row of secondary conductors
    :type secondary_row: ConductorRow
    :param tertiary_row: row of tertiary conductors
    :type tertiary_row: ConductorRow
    :return:
    """
    if secondary_row.row_height == tertiary_row.row_height and \
            secondary_row.number_of_conds_per_row == tertiary_row.number_of_conds_per_row and \
            secondary_row.number_of_rows == tertiary_row.number_of_rows and \
            secondary_row.number_of_conds_per_winding == tertiary_row.number_of_conds_per_winding:
        return True
    else:
        return False


# def stack_group_of_rows():
#     # How many complete groups have to be stacked in the window?
#     if center_tapped_group.secondary_number_of_rows == 1 and center_tapped_group.secondary_rest == 0:
#         number_of_groups = secondary_row.number_of_rows
#     elif center_tapped_group.primary_number_of_rows == 1 and center_tapped_group.primary_rest == 0:
#         number_of_groups = primary_row.number_of_rows
#
#     # Initialize the winding objects "to be placed"
#     single_primary_rows_to_be_placed = [primary_row] * (primary_row.number_of_rows - number_of_groups * \
#         center_tapped_group.primary_number_of_rows)
#     single_secondary_rows_to_be_placed = [secondary_row] * (secondary_row.number_of_rows - number_of_groups * \
#         center_tapped_group.secondary_number_of_rows)
#     groups_to_be_placed = [center_tapped_group] * number_of_groups
#
#     # Depending on the rows that could not be covered by the groups, integrate them among the groups
#     if len(single_secondary_rows_to_be_placed) == 0:
#         stack_order = mix_x_and_I(single_primary_rows_to_be_placed, groups_to_be_placed)
#     elif len(single_primary_rows_to_be_placed) == 0:
#         stack_order = mix_x_and_I(single_secondary_rows_to_be_placed, groups_to_be_placed)
#
#     # Add the tertiary winding to the stack_order
#     add_tertiary_winding_to_stack(stack_order, tertiary_row)
#
#     # Insert insulations into the stack_order
#     def insert_insulations_to_stack(stack_order, isolations):
#         # here: zero index means bot_row
#         for i, bot_row in enumerate(stack_order[0:-1]):
#             top_row = stack_order[i + 1]
#             # Which insulation is needed depends on the bot and top row neighbours
#             if type(top_row


def group_center_tapped(primary_number_of_rows: int, secondary_number_of_rows: int,
                        primary_row, secondary_row, tertiary_row,
                        isolations: ThreeWindingIsolation):
    """
    Group center tapped windings.

    :param primary_number_of_rows: number of primary rows
    :type primary_number_of_rows: int
    :param secondary_number_of_rows: number of secondary rows
    :type secondary_number_of_rows: int
    :param primary_row: row of primary conductors
    :type primary_row: ConductorRow
    :param secondary_row: row of secondary conductors
    :type secondary_row: ConductorRow
    :param tertiary_row: row of tertiary conductors
    :type tertiary_row: ConductorRow
    :param isolations: isolation DTO
    :type isolations: ThreeWindingIsolation
    """
    group = CenterTappedGroup(primary_number_of_rows=None, secondary_number_of_rows=None,
                              primary_rest=None, secondary_rest=None, stack=None)
    # TODO: rests do not really belong into each group but only once
    if primary_number_of_rows > secondary_number_of_rows:
        group.primary_number_of_rows = int(primary_number_of_rows / secondary_number_of_rows)
        group.secondary_number_of_rows = 1
        group.primary_rest = primary_number_of_rows % secondary_number_of_rows
        group.secondary_rest = 0
    elif secondary_number_of_rows > primary_number_of_rows:
        group.primary_number_of_rows = 1
        group.secondary_number_of_rows = int(secondary_number_of_rows / primary_number_of_rows)
        group.primary_rest = 0
        group.secondary_rest = secondary_number_of_rows % primary_number_of_rows
    elif primary_number_of_rows == secondary_number_of_rows:
        group.primary_number_of_rows = 1
        group.secondary_number_of_rows = 1
        group.primary_rest = 0
        group.secondary_rest = 0

    # Currently the group is simply stacked with all primaries and after that all secondaries & tertiaries
    # TODO: make positioning of conductor rows more flexible
    # TODO: make insertion of inner insulations also flexible and in a dedicated function
    stack = []
    for i in range(0, group.primary_number_of_rows):
        stack.append(primary_row)
        if i < group.primary_number_of_rows - 1:
            stack.append(StackIsolation(thickness=isolations.primary_to_primary))
        elif i == group.primary_number_of_rows - 1:
            stack.append(StackIsolation(thickness=isolations.primary_to_secondary))

    for i in range(0, group.secondary_number_of_rows):
        stack.append(secondary_row)
        stack.append(StackIsolation(thickness=isolations.secondary_to_tertiary))
        stack.append(tertiary_row)
        if i < int(secondary_number_of_rows/2) - 1:
            stack.append(StackIsolation(thickness=isolations.tertiary_to_secondary))
        elif i == int(secondary_number_of_rows/2) - 1:
            pass

    group.stack = stack

    return group


def get_height_of_group(group: CenterTappedGroup) -> float:
    """
    Return the total height of the conductors and insulation.

    :param group: center tapped group
    :type group: CenterTappedGroup
    :return: height of the group
    :rtype: float
    """
    total_height = 0
    for row_element in group.stack:
        if isinstance(row_element, ConductorRow):
            total_height += row_element.row_height
        elif isinstance(row_element, StackIsolation):
            total_height += row_element.thickness
    return total_height


def add_tertiary_winding_to_stack(stack_order_without_tertiary: list[ConductorRow],
                                  tertiary_row_to_be_added: list[ConductorRow]):
    """
    Add tertiary winding to stack.

    :param stack_order_without_tertiary: list of ConductorRows
    :type stack_order_without_tertiary: list[ConductorRow]
    :param tertiary_row_to_be_added: list of ConductorRows
    :type tertiary_row_to_be_added: list[ConductorRow]
    """
    secondary_tags = []
    number_of_added_tertiaries = 0
    for i, obj in enumerate(stack_order_without_tertiary):
        if isinstance(obj, ConductorRow):
            if obj.winding_tag == WindingTag.Secondary:
                secondary_tags.append(i + 1 + number_of_added_tertiaries)
                number_of_added_tertiaries += 1
    for tag in secondary_tags:
        stack_order_without_tertiary.insert(tag, tertiary_row_to_be_added)


# Insert insulations into the stack_order
def insert_insulations_to_stack(stack_order: list, isolations: ThreeWindingIsolation):
    """
    Insert insulations to a stack.

    :param stack_order: order of the stack
    :type stack_order: list
    :param isolations: isolations according to the ThreeWindingIsolation class
    :type isolations: ThreeWindingIsolation
    """
    # Which insulation is needed depends on the bot and top row neighbours
    # TODO: insert insulations into the stack depending on isolation matrix

    # TODO: remove following fix insert
    insulation_positions = []
    insulation_tags = []
    insulation_stack = []
    number_of_added_insulations = 0
    # here: zero index means bot_row
    for i, _ in enumerate(stack_order[0:-1]):
        insulation_string = ""
        # print(f"{bot_row = }")
        # print(f"{top_row = }\n")
        insulation_positions.append(i + 1 + number_of_added_insulations)
        number_of_added_insulations += 1

        # TODO: dynamic isolations (vertical)
        for n, row in enumerate([stack_order[i], stack_order[i + 1]]):
            if isinstance(row, ConductorRow):
                if row.winding_tag == WindingTag.Primary:
                    insulation_string += "primary"
                elif row.winding_tag == WindingTag.Secondary:
                    insulation_string += "secondary"
                elif row.winding_tag == WindingTag.Tertiary:
                    insulation_string += "tertiary"
            elif isinstance(row, CenterTappedGroup):
                position_in_group = n-1  # 0 for the bot and -1 for the top element
                if row.stack[position_in_group].winding_tag == WindingTag.Primary:
                    insulation_string += "primary"
                elif row.stack[position_in_group].winding_tag == WindingTag.Secondary:
                    insulation_string += "secondary"
                elif row.stack[position_in_group].winding_tag == WindingTag.Tertiary:
                    insulation_string += "tertiary"

            if n == 0:
                insulation_string += "_to_"

        insulation_tags.append(insulation_string)
        insulation_stack.append(StackIsolation(thickness=getattr(isolations, insulation_string)))

    # Fill in the isolations
    for i, position in enumerate(insulation_positions):
        stack_order.insert(position, insulation_stack[i])

    return insulation_tags
    # print(f"{insulation_positions = }")
    # print(f"{insulation_stack = }")
    # print(f"{insulation_tags = }")


def get_set_of_integers_from_string_list(string_list: list[str]):
    """Get the list of the set of integers in a list of strings.

    Used by winding_scheme key.

    :param string_list:
    :type string_list: list[str]
    """
    integer_list = []
    for single_string in string_list:
        try:
            if int(single_string) not in integer_list:
                integer_list.append(int(single_string))
        except:
            pass
    return integer_list


def stack_order_from_interleaving_scheme(interleaving_scheme: InterleavingSchemesFoilLitz, stack_order: list,
                                         primary_additional_bobbin: float, center_foil_additional_bobbin: float,
                                         primary_row: ConductorRow, secondary_row: ConductorRow,
                                         tertiary_row: ConductorRow, insulations: ThreeWindingIsolation):
    """
    Get the stack order from a given interleaving scheme.

    :param interleaving_scheme: interleaving scheme
    :type interleaving_scheme: InterleavingSchemesFoilLitz
    :param stack_order: order of the stack
    :type stack_order: list
    :param primary_additional_bobbin: primary additional bobbin thickness in m in the center
    :type primary_additional_bobbin: float
    :param center_foil_additional_bobbin: additional bobbin thickness in m for the center foil
    :type center_foil_additional_bobbin: float
    :param primary_row: row of primary conductors
    :type primary_row: ConductorRow
    :param secondary_row: row of secondary conductors
    :type secondary_row: ConductorRow
    :param tertiary_row: row of tertiary conductors
    :type tertiary_row: ConductorRow
    :param insulations: insulation distances from the ThreeWindingIsolation class
    :type insulations: ThreeWindingIsolation
    """
    # Init the winding counters (needed for vertical insulation adjustments)
    number_of_primary_rows, number_of_secondary_rows, number_of_tertiary_rows = 0, 0, 0

    # Split the interleaving scheme directive to obtain the single stacked winding parts
    string_stack_order = interleaving_scheme.split(sep="_")

    # Needed to shift all other primary rows half a position to the right (quasi hexagonal)
    max_primary_conductors_per_row = max(get_set_of_integers_from_string_list(string_stack_order))

    # Iterate over the stack
    for no_element, element in enumerate(string_stack_order):

        # Foil windings are handled first: Secondary and Tertiary (usually symmetric)
        if element == "sec" or element == "ter":
            if element == "sec":
                temp_row = secondary_row
                number_of_secondary_rows += 1
            elif element == "ter":
                temp_row = tertiary_row
                number_of_tertiary_rows += 1

            # Stacking of the foils does not depend on secondary and tertiary
            # center_foil_additional_bobbin affects only, if the foil is in the inner 50 % of all rows
            # (close to the air gap)
            center_interval = 0.5
            if no_element < len(string_stack_order)*(center_interval/2) or no_element >= len(string_stack_order) * \
                    (1/2+center_interval/2):
                stack_order.append(temp_row)
            else:
                center_temp_row = copy.deepcopy(temp_row)
                center_temp_row.additional_bobbin += center_foil_additional_bobbin
                stack_order.append(center_temp_row)

        # Primary litz wire is handled secondly
        else:
            temp_row = copy.deepcopy(primary_row)
            temp_row.number_of_conds_per_row = int(element)
            temp_row.additional_bobbin = primary_additional_bobbin

            if temp_row.number_of_conds_per_row < max_primary_conductors_per_row:
                temp_row.additional_bobbin += temp_row.row_height / 2 + insulations.primary_to_primary / 2
            stack_order.append(temp_row)
            number_of_primary_rows += 1

    return number_of_primary_rows, number_of_secondary_rows, number_of_tertiary_rows


def adjust_vertical_insulation_center_tapped_stack(interleaving_scheme: InterleavingSchemesFoilLitz, primary_row: ConductorRow,
                                                   secondary_row: ConductorRow, tertiary_row: ConductorRow,
                                                   primary_additional_bobbin: float, center_foil_additional_bobbin: float,
                                                   isolations: ThreeWindingIsolation, available_height: float):
    """Adjust the vertical insulation of a center tapped stack.

    The insulation is equally distributed.

    :param primary_row: primary row according to the ConductorRow class
    :type primary_row: ConductorRow
    :param secondary_row: secondary row according to the ConductorRow class
    :type secondary_row: ConductorRow
    :param tertiary_row: tertiary row according to the ConductorRow class
    :type tertiary_row: ConductorRow
    :param isolations: isolations according to the ThreeWindingIsolation class
    :type isolations: ThreeWindingIsolation
    :param available_height: available height in meter
    :type available_height: float
    :param primary_additional_bobbin: additional bobbin thickness in meter
    :type primary_additional_bobbin: float
    :param center_foil_additional_bobbin: additional bobbin in meter for the center foil
    :type center_foil_additional_bobbin: float
    :param interleaving_scheme: interleaving scheme (e.g. ter_3_4_ter_sec_4_3_sec / ter_4_3_ter_sec_3_4_sec / ...)
    :type interleaving_scheme: InterleavingSchemesFoilLitz
    """
    # Initial stacking with predefined minimum insulations
    initial_stack_order = []

    number_of_primary_rows, number_of_secondary_rows, number_of_tertiary_rows = \
        stack_order_from_interleaving_scheme(interleaving_scheme, initial_stack_order, primary_additional_bobbin,
                                             center_foil_additional_bobbin,
                                             primary_row, secondary_row, tertiary_row, isolations)

    insulation_tags = insert_insulations_to_stack(initial_stack_order, isolations)
    number_of_insulations_primary_to_primary = insulation_tags.count("primary_to_primary")
    number_of_insulations_primary_to_secondary = insulation_tags.count("primary_to_secondary") + insulation_tags.count("secondary_to_primary") + \
        insulation_tags.count("primary_to_tertiary") + insulation_tags.count("tertiary_to_primary")
    number_of_insulations_secondary_to_secondary = insulation_tags.count("secondary_to_secondary") + insulation_tags.count("tertiary_to_tertiary") + \
        insulation_tags.count("secondary_to_tertiary") + insulation_tags.count("tertiary_to_secondary")

    # 1
    minimum_needed_height = number_of_primary_rows * primary_row.row_height + number_of_secondary_rows * secondary_row.row_height + \
        number_of_tertiary_rows * secondary_row.row_height + number_of_insulations_primary_to_primary * isolations.primary_to_primary + \
        number_of_insulations_secondary_to_secondary * isolations.secondary_to_tertiary + \
        number_of_insulations_primary_to_secondary * isolations.primary_to_secondary

    # 2
    new_iso_parameter = np.round((available_height - minimum_needed_height + number_of_insulations_primary_to_secondary * \
                                  isolations.primary_to_secondary) / number_of_insulations_primary_to_secondary, 9) - 1e-9
    if new_iso_parameter > isolations.primary_to_secondary:
        isolations.primary_to_secondary = new_iso_parameter
        isolations.primary_to_tertiary = new_iso_parameter
        isolations.tertiary_to_primary = new_iso_parameter
        isolations.secondary_to_primary = new_iso_parameter


def stack_center_tapped_transformer(primary_row: ConductorRow, secondary_row: ConductorRow, tertiary_row: ConductorRow,
                                    isolations: ThreeWindingIsolation, available_height: float,
                                    interleaving_type: CenterTappedInterleavingType,
                                    interleaving_scheme: InterleavingSchemesFoilLitz,
                                    primary_additional_bobbin: float, center_foil_additional_bobbin: float):
    """
    Define the vertical stacking of previously defined ConductorRows.

    IMPORTANT DEFINITION: the rows are calculated without taking into account
    any vertical insulation. Only the horizontal insulation from conductors
    of the same WindingTag is taken into account in ConductorRow definition.

    :param primary_row: primary row according to the ConductorRow class
    :type primary_row: ConductorRow
    :param secondary_row: secondary row according to the ConductorRow class
    :type secondary_row: ConductorRow
    :param tertiary_row: tertiary row according to the ConductorRow class
    :type tertiary_row: ConductorRow
    :param isolations: isolations according to the ThreeWindingIsolation class
    :type isolations: ThreeWindingIsolation
    :param available_height: available height in meter
    :type available_height: float
    :param interleaving_type: custom / TypeA / TypeB / TypeC / TypeD
    :type interleaving_type: CenterTappedInterleavingType
    :param interleaving_scheme: interleaving scheme (e.g. ter_3_4_ter_sec_4_3_sec / ter_4_3_ter_sec_3_4_sec / ...)
    :type interleaving_scheme: InterleavingSchemesFoilLitz
    :param primary_additional_bobbin: additional bobbin thickness in meter
    :type primary_additional_bobbin: float
    :param center_foil_additional_bobbin: additional bobbin in meter for the center foil
    :type center_foil_additional_bobbin: float
    """
    if not check_secondary_and_tertiary_are_the_same(secondary_row, tertiary_row):
        logger.warning("Secondary and tertiary winding are not defined similar. "
                       "That is not a nice center-tapped transformer :(")

    elif interleaving_type == CenterTappedInterleavingType.TypeA:
        # Usually for center-tapped it is the secondary/=tertiary winding number
        # But for parallel windings this can also differ -> a general approach is implemented
        # what is the smallest number of rows?
        center_tapped_group = group_center_tapped(primary_row.number_of_rows, secondary_row.number_of_rows,
                                                  primary_row, secondary_row, tertiary_row,
                                                  isolations)

        # How many complete groups have to be stacked in the window?
        if center_tapped_group.secondary_number_of_rows == 1 and center_tapped_group.secondary_rest == 0:
            number_of_groups = secondary_row.number_of_rows
        elif center_tapped_group.primary_number_of_rows == 1 and center_tapped_group.primary_rest == 0:
            number_of_groups = primary_row.number_of_rows

        # Initialize the winding objects "to be placed"
        # TODO: is it possible handle the rest number of primary conductors for the last row? not definitely needed here
        single_primary_rows_to_be_placed = [primary_row] * (primary_row.number_of_rows - number_of_groups * center_tapped_group.primary_number_of_rows)
        single_secondary_rows_to_be_placed = [secondary_row] * (secondary_row.number_of_rows - number_of_groups * center_tapped_group.secondary_number_of_rows)
        groups_to_be_placed = [center_tapped_group] * number_of_groups

        # Depending on the rows that could not be covered by the groups, integrate them among the groups
        number_of_single_rows = 0
        if len(single_primary_rows_to_be_placed) != 0:
            stack_order = mix_x_and_i(single_primary_rows_to_be_placed, groups_to_be_placed)
            number_of_single_rows = len(single_primary_rows_to_be_placed)
        elif len(single_secondary_rows_to_be_placed) != 0:
            stack_order = mix_x_and_i(single_secondary_rows_to_be_placed, groups_to_be_placed)
            number_of_single_rows = len(single_secondary_rows_to_be_placed)
        else:
            stack_order = mix_x_and_i(single_secondary_rows_to_be_placed, groups_to_be_placed)

        # Add the tertiary winding to the stack_order
        add_tertiary_winding_to_stack(stack_order, tertiary_row)

        insert_insulations_to_stack(stack_order, isolations)

        # Create the complete ConductorStack from the stack_order
        return ConductorStack(number_of_groups=number_of_groups, number_of_single_rows=number_of_single_rows,
                              order=stack_order)

    elif interleaving_type == CenterTappedInterleavingType.TypeB:
        number_of_single_rows = None
        stack_order = []

        stack_order.append(tertiary_row)
        stack_order.append(primary_row)
        stack_order.append(secondary_row)
        stack_order.append(tertiary_row)
        stack_order.append(primary_row)
        stack_order.append(secondary_row)
        stack_order.append(tertiary_row)
        stack_order.append(primary_row)
        stack_order.append(secondary_row)

        insert_insulations_to_stack(stack_order, isolations)

        # Create the complete ConductorStack from the stack_order
        return ConductorStack(number_of_groups=0, number_of_single_rows=number_of_single_rows, order=stack_order)

    elif interleaving_type == CenterTappedInterleavingType.TypeC:

        # Implicit usage of the available height to ensure symmetrical placement of secondary and tertiary winding
        adjust_vertical_insulation_center_tapped_stack(interleaving_scheme, primary_row, secondary_row, tertiary_row,
                                                       primary_additional_bobbin, center_foil_additional_bobbin,
                                                       isolations, available_height)

        # Placing with adjusted isolation
        stack_order = []
        stack_order_from_interleaving_scheme(interleaving_scheme, stack_order,
                                             primary_additional_bobbin, center_foil_additional_bobbin,
                                             primary_row, secondary_row, tertiary_row,
                                             isolations)

        # Add insulations afterward
        insert_insulations_to_stack(stack_order, isolations)

        # Create the complete ConductorStack from the stack_order
        return ConductorStack(number_of_groups=0, number_of_single_rows=None, order=stack_order)

    elif interleaving_type == CenterTappedInterleavingType.TypeD:
        number_of_single_rows = None
        stack_order = []
        for _ in range(0, primary_row.number_of_rows):
            stack_order.append(primary_row)
        for _ in range(0, secondary_row.number_of_rows):
            stack_order.append(secondary_row)
        for _ in range(0, tertiary_row.number_of_rows):
            stack_order.append(tertiary_row)

        insert_insulations_to_stack(stack_order, isolations)

        # Create the complete ConductorStack from the stack_order
        return ConductorStack(number_of_groups=0, number_of_single_rows=number_of_single_rows, order=stack_order)


def get_number_of_turns_in_groups(stack: ConductorStack):
    """
    Get the number of turns in one center tapped stack for Type A winding scheme. Experimental.

    :param stack: conductor stack
    :type stack: ConductorStack
    """
    turns1 = 0
    turns2 = 0

    for row_element in stack.order:
        if isinstance(row_element, StackIsolation):
            pass
        else:
            if isinstance(row_element, ConductorRow):
                pass
            elif isinstance(row_element, CenterTappedGroup):
                for row in row_element.stack:
                    if isinstance(row, StackIsolation):
                        pass
                    elif isinstance(row, ConductorRow):
                        if row.winding_tag == WindingTag.Primary:
                            turns1 += row.number_of_conds_per_row
                        elif row.winding_tag == WindingTag.Secondary:
                            turns2 += row.number_of_conds_per_row
    return turns1, turns2


def is_even(x: int):
    """Check if given number is even or odd.

    :param x: input number to check
    :type x: int
    """
    if x % 2:
        return False
    else:
        return True


def center(l_list: list):
    """
    Return the center index of a list. Rounds off.

    :param l_list: list to return the center index
    :type l_list: list
    """
    return int(len(l_list) / 2)


def mix_x_and_i(input_x: list, input_i: list):
    """General usage to interleave windings. One winding could be input_x and the other input_i. Experimental.

    Example: 16 primary windings (input_x), 3 secondary windings (input_i).
    Tries to fit these windings symmetric into the winding window.

    :param input_x: List of e.g. primary windings
    :type input_x: list
    :param input_i: List of e.g. secondary windings
    :type input_i: list
    """
    len_x = len(input_x)
    len_i = len(input_i)
    if len_x > len_i:
        logger.warning("x must be smaller or equal I")
    else:
        if len_x == 0:
            return input_i
        else:
            x = [0] * len_x
            current = [1] * len_i

            if is_even(len_x):
                # If an even number shall be distributed, always start from the outsides
                count = 1
                mix_distance = 1 + int(len_i / len_x)  # TODO:distance adaptive not always 1
                while len(x) > 0:
                    if count > 0:
                        current.insert(count, x[0])
                        count = - count
                    elif count < 0:
                        current.insert(count, x[0])
                        count = -count + mix_distance
                    x.pop(0)
            elif not is_even(len_x):
                temp_x = x[0]
                temp_mixed = mix_x_and_i(x[0:-1], current)
                temp_mixed.insert(center(temp_mixed), temp_x)
                current = temp_mixed
            input_list = [input_x[0], input_i[0]]
            return [input_list[i] for i in current]
