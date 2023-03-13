from femmt import *
from femmt.dtos import *

def define_center_tapped_insulation(primary_to_primary, secondary_to_secondary, primary_to_secondary):
    """
    This function defines the ThreeWindingIsolation dto for the special case of a center-tapped
    transformer. It's assumed, that the secondary windings are isolated symmetrically
    :param primary_to_primary:
    :param secondary_to_secondary:
    :param primary_to_secondary:
    :return:
    """
    return ThreeWindingIsolation(primary_to_primary=primary_to_primary,
                                 primary_to_secondary=primary_to_secondary,
                                 primary_to_tertiary=primary_to_secondary,
                                 secondary_to_primary=primary_to_secondary,
                                 secondary_to_secondary=secondary_to_secondary,
                                 secondary_to_tertiary=secondary_to_secondary,
                                 tertiary_to_primary=primary_to_secondary,
                                 tertiary_to_secondary=secondary_to_secondary,
                                 tertiary_to_tertiary=secondary_to_secondary)
