"""Different functions to describe a model."""
from femmt.dtos import *


def define_center_tapped_insulation(primary_to_primary: float, secondary_to_secondary: float, primary_to_secondary: float):
    """
    Define the ThreeWindingIsolation dto for the special case of a center-tapped transformer.

    It's assumed, that the secondary windings are isolated symmetrically.
    :param primary_to_primary: Primary winding to primary winding insulation in m
    :type primary_to_primary: float
    :param secondary_to_secondary: Secondary winding to secondary winding insulation in m
    :type secondary_to_secondary: float
    :param primary_to_secondary: Primary winding to secondary winding insulation in m
    :type primary_to_secondary: float
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
