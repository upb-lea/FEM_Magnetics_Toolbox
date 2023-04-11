from femmt import *
from femmt.enumerations import *
from femmt.dtos import *


def number_of_rows(row: ConductorRow):
    """
    :return:
    """
    number_rows = int(row.number_of_conds_per_winding / row.number_of_conds_per_row) + (row.number_of_conds_per_winding % row.number_of_conds_per_row > 0)
    return number_rows


def single_row(number_of_conds_per_winding, window_width, winding_tag: WindingTag, conductor_type: ConductorType, thickness=None, radius=None, cond_cond_isolation=None):
    """
    Defines a full row of the defined conductor in the specified window width.
    It is assumed, that the row is full, which means, as many conductors are
    place horizontally next to each other, as possible.
    It is assumed, that the row fits into the winding window. A later control
    is needed.
    :return:
    """
    # Create ConductorRow dataclass
    conductor_row = ConductorRow(number_of_conds_per_winding=number_of_conds_per_winding,
                                 number_of_conds_per_row=None,
                                 row_height=None,
                                 winding_tag=winding_tag,
                                 number_of_rows=None)

    if conductor_type == ConductorType.RoundLitz or conductor_type == ConductorType.RoundSolid:
        conductor_row.row_height = radius * 2
        # How many full conductor cross sections fit into one row of the winding window?
        conductor_row.number_of_conds_per_row = 0
        while window_width > 2 * (conductor_row.number_of_conds_per_row + 1) * radius + cond_cond_isolation * conductor_row.number_of_conds_per_row:
            conductor_row.number_of_conds_per_row += 1
    elif conductor_type == ConductorType.RectangularSolid:
        # Rectangular Solid Conductors are meant to be as possible, so their cross section is
        # simply adjusted to the winding window width
        conductor_row.row_height = thickness
        conductor_row.number_of_conds_per_row = 1

    conductor_row.number_of_rows = number_of_rows(conductor_row)
    return conductor_row


def check_secondary_and_tertiary_are_the_same(secondary_row: ConductorRow, tertiary_row: ConductorRow):
    """
    Does the definition of the single rows relate to a center-tapped transformer?
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
#     single_primary_rows_to_be_placed = [primary_row] * (primary_row.number_of_rows - number_of_groups * center_tapped_group.primary_number_of_rows)
#     single_secondary_rows_to_be_placed = [secondary_row] * (secondary_row.number_of_rows - number_of_groups * center_tapped_group.secondary_number_of_rows)
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


def group_center_tapped(primary_number_of_rows, secondary_number_of_rows,
                        primary_row, secondary_row, tertiary_row,
                        isolations: ThreeWindingIsolation):
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
        if i < secondary_number_of_rows - 1:
            stack.append(StackIsolation(thickness=isolations.tertiary_to_secondary))
        elif i == secondary_number_of_rows - 1:
            pass

    group.stack = stack

    return group


def get_height_of_group(group: CenterTappedGroup):
    """returns the total height of thr conductors and insulation
    """
    total_height = 0
    for row_element in group.stack:
        if type(row_element) == ConductorRow:
            total_height += row_element.row_height
        elif type(row_element) == StackIsolation:
            total_height += row_element.thickness
    return total_height


def add_tertiary_winding_to_stack(stack_order_without_tertiary, tertiary_row_to_be_added):
    secondary_tags = []
    number_of_added_tertiaries = 0
    for i, obj in enumerate(stack_order_without_tertiary):
        if type(obj) == ConductorRow:
            if obj.winding_tag == WindingTag.Secondary:
                secondary_tags.append(i + 1 + number_of_added_tertiaries)
                number_of_added_tertiaries += 1
    for tag in secondary_tags:
        stack_order_without_tertiary.insert(tag, tertiary_row_to_be_added)


# Insert insulations into the stack_order
def insert_insulations_to_stack(stack_order, isolations: ThreeWindingIsolation):
    # Which insulation is needed depends on the bot and top row neighbours
    # TODO: insert insulations into the stack depending on isolation matrix
    # if type(top_row

    # TODO: remove following fix insertation
    insulation_positions = []
    insulation_tags = []
    number_of_added_insulations = 0
    # here: zero index means bot_row
    for i, bot_row in enumerate(stack_order[0:-1]):
        top_row = stack_order[i + 1]
        # print(f"{bot_row = }")
        # print(f"{top_row = }\n")
        insulation_positions.append(i + 1 + number_of_added_insulations)
        number_of_added_insulations += 1

        # if type(bot_row) == ConductorRow:  # TODO: so ähnlich nach group oder row fragen und dann die infos über die isolation rausszierehn bei bot row und top row
        #         if obj.winding_tag == WindingTag.Secondary:
        insulation_tags.append(StackIsolation(thickness=isolations.secondary_to_secondary))

    # Fill in the isolations
    for i, position in enumerate(insulation_positions):
        stack_order.insert(position, insulation_tags[i])


def stack_center_tapped_transformer(primary_row: ConductorRow, secondary_row: ConductorRow, tertiary_row: ConductorRow,
                                    window_height, isolations: ThreeWindingIsolation, interleaving_type: CenterTappedInterleavingType):
    """Defines the vertical stacking of previously defined ConductorRows.
    IMPORTANT DEFINITION: the rows are calculated without taking into account
    any vertical insulation. Only the horizontal insulation from conductors
    of the same WindingTag is taken into account in ConductorRow definition.
    """

    if not check_secondary_and_tertiary_are_the_same(secondary_row, tertiary_row):
        print("Secondary and tertiary winding are not defined similar. "
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

        # Initialize the winding objects "to be placed"  # TODO: is it possible handle the rest number of primary conductors for the last row? not definitely needed here...
        single_primary_rows_to_be_placed = [primary_row] * (primary_row.number_of_rows - number_of_groups * center_tapped_group.primary_number_of_rows)
        single_secondary_rows_to_be_placed = [secondary_row] * (secondary_row.number_of_rows - number_of_groups * center_tapped_group.secondary_number_of_rows)
        groups_to_be_placed = [center_tapped_group] * number_of_groups

        # Depending on the rows that could not be covered by the groups, integrate them among the groups
        number_of_single_rows = 0
        if len(single_primary_rows_to_be_placed) != 0:
            stack_order = mix_x_and_I(single_primary_rows_to_be_placed, groups_to_be_placed)
            number_of_single_rows = len(single_primary_rows_to_be_placed)
        elif len(single_secondary_rows_to_be_placed) != 0:
            stack_order = mix_x_and_I(single_secondary_rows_to_be_placed, groups_to_be_placed)
            number_of_single_rows = len(single_secondary_rows_to_be_placed)
        else:
            stack_order = mix_x_and_I(single_secondary_rows_to_be_placed, groups_to_be_placed)

        # Add the tertiary winding to the stack_order
        add_tertiary_winding_to_stack(stack_order, tertiary_row)


        insert_insulations_to_stack(stack_order, isolations)

        # Create the complete ConductorStack from the stack_order
        return ConductorStack(number_of_groups=number_of_groups, number_of_single_rows=number_of_single_rows, order=stack_order)

    elif interleaving_type == CenterTappedInterleavingType.TypeB:

        number_of_single_rows = None
        stack_order = None
        # Create the complete ConductorStack from the stack_order
        return ConductorStack(number_of_groups=0, number_of_single_rows=number_of_single_rows, order=stack_order)


def get_number_of_turns_in_groups(stack):
    turns1 = 0
    turns2 = 0

    for row_element in stack.order:
        if type(row_element) == StackIsolation:
            pass
        else:
            if type(row_element) == ConductorRow:
                pass
            elif type(row_element) == CenterTappedGroup:
               for row in row_element.stack:
                    if type(row) == StackIsolation:
                        pass
                    elif type(row) == ConductorRow:
                        if row.winding_tag == WindingTag.Primary:
                            turns1 += row.number_of_conds_per_row
                        elif row.winding_tag == WindingTag.Secondary:
                            turns2 += row.number_of_conds_per_row
    return turns1, turns2


def is_even(x: int):
    if x % 2:
        return False
    else:
        return True


def center(l: List):
    return int(len(l) / 2)


def mix_x_and_I(input_x, input_I):
    len_x = len(input_x)
    len_I = len(input_I)
    if len_x > len_I:
        print("x must be smaller or equal I")
    else:
        if len_x == 0:
            return input_I
        else:
            x = [0] * len_x
            I = [1] * len_I

            if is_even(len_x):
                # If an even number shall be distributed, always start from the outsides
                count = 1
                mix_distance = 1 + int(len_I / len_x)  # TODO:distance adaptive not always 1
                while len(x) > 0:
                    if count > 0:
                        I.insert(count, x[0])
                        count = - count
                    elif count < 0:
                        I.insert(count, x[0])
                        count = -count + mix_distance
                    x.pop(0)
            elif not is_even(len_x):
                temp_x = x[0]
                temp_mixed = mix_x_and_I(x[0:-1], I)
                temp_mixed.insert(center(temp_mixed), temp_x)
                I = temp_mixed
            # print(f"{I = }")
            input_list = [input_x[0], input_I[0]]
            return [input_list[i] for i in I]

