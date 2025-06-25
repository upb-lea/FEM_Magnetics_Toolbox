"""Different functions to describe a model."""
from femmt.dtos import *
import numpy as np
import os

def write_permeability_pro_file(parent_directory,
                                b_ref_vec: np.ndarray | list = None,
                                mu_r_real_vec: np.ndarray | list = None,
                                mu_r_imag_vec: np.ndarray | list = None,
                                silent: bool = False):
    """
    Export data from the material database in a certain file format.

    :param parent_directory: path to parent directory
    :type parent_directory: str
    :param b_ref_vec: reference vector for mu_r_real and mu_r_imag
    :type b_ref_vec: ndarray or list
    :param mu_r_real_vec: real part of mu_r_abs as a vector
    :type mu_r_real_vec: ndarray or list
    :param mu_r_imag_vec: imaginary part of mu_r_abs as a vector
    :type mu_r_imag_vec: ndarray or list
    :param silent: enables/disables print
    :type silent: bool
    """
    # fix numpy array inside normal python list problem
    # converts everything from scratch to a list, unified file format.
    b_ref_vec = np.array(b_ref_vec).tolist()
    mu_r_real_vec = np.array(mu_r_real_vec).tolist()
    mu_r_imag_vec = np.array(mu_r_imag_vec).tolist()
    with open(os.path.join(parent_directory, "core_materials_temp.pro"), "w") as file:
        file.write('Include "Parameter.pro";\n')
        file.write(
            f"Function{{\n  b = {str(b_ref_vec).replace('[', '{').replace(']', '}')} ;\n  "
            f"mu_real = {str(mu_r_real_vec).replace('[', '{').replace(']', '}')} ;"
            f"\n  mu_imag = {str(mu_r_imag_vec).replace('[', '{').replace(']', '}')} ;\n  "
            f"mu_imag_couples = ListAlt[b(), mu_imag()] ;\n  "
            f"mu_real_couples = ListAlt[b(), mu_real()] ;\n  "
            f"f_mu_imag_d[] = InterpolationLinear[Norm[$1]]{{List[mu_imag_couples]}};\n  "
            f"f_mu_real_d[] = InterpolationLinear[Norm[$1]]{{List[mu_real_couples]}};\n  "
            f"f_mu_imag[] = f_mu_imag_d[$1];\n  "
            f"f_mu_real[] = f_mu_real_d[$1];\n }}  ")

    if not silent:
        print(f"Data is exported to {parent_directory} in a .pro-file.")


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
