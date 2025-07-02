"""Contains different functions, used by the whole FEMMT functions."""
# Python standard libraries
# Python standard libraries
import json
import pkg_resources
import subprocess
import sys
import os
import warnings
from scipy.integrate import quad
import logging

# Third parry libraries
import gmsh
from matplotlib import pyplot as plt
import numpy.typing as npt
import numpy as np

# Local libraries
from femmt.constants import *
from femmt.enumerations import ConductorType
from femmt.dtos import *
import femmt.functions_reluctance as fr

logger = logging.getLogger(__name__)

colors_femmt_default = {"blue": (28, 113, 216),
                        'red': (192, 28, 40),
                        "green": (46, 194, 126),
                        "orange": (230, 97, 0),
                        "purple": (129, 61, 156),
                        "brown": (134, 94, 60),
                        "gray": (119, 118, 123),
                        "yellow": (245, 194, 17),
                        "black": (0, 0, 0),
                        "white": (255, 255, 255)
                        }

colors_geometry_femmt_default = {
    "core": "black",
    "air_gap": "yellow",
    "winding": ["orange", "brown", "yellow", "green", "red", "black", "gray", "blue", "orange", "purple", "gray",
                "blue", "orange", "purple"],
    "insulation": "blue",
    "potting_inner": "gray",
    "potting_outer": "gray",
}

colors_ba_jonas = {"blue": (28, 113, 216),
                   'red': (213, 6, 6),
                   "green": (6, 213, 6),
                   "orange": (230, 97, 0),
                   "purple": (129, 61, 156),
                   "brown": (134, 94, 60),
                   "gray": (193, 193, 193),
                   "yellow": (255, 171, 6),
                   "black": (58, 58, 58),
                   "white": (255, 255, 255),
                   "gray_dark": (109, 109, 109),
                   "gray_dark_dark": (50, 50, 50)
                   }

colors_geometry_ba_jonas = {
    "core": "black",
    "air_gap": "yellow",
    "winding": ["green", "red", "yellow"],
    "insulation": "gray_dark",
    "potting_inner": "gray",
    "potting_outer": "gray_dark_dark",
}

colors_geometry_draw_only_lines = {
    "core": "gray_dark",
    "air_gap": "gray_dark",
    "winding": ["green", "green", "green"],
    "insulation": "gray_dark",
    "potting_inner": "gray_dark",
    "potting_outer": "gray_dark",
}


def core_database() -> dict:
    """
    Return a core geometry for defined core structure.

    All dimensions are nominal dimensions without consideration of tolerances.

    For PQ core sizes (e.g. PQ 40/40), it has been found out that
    core_dimension_x / core_dimension_y = 1.45, the error over all available shapes is maximum 7%
    (compared to datasheet value)
    Derivation:
    core_list: ['PQ 20/20', 'PQ 26/25', 'PQ 32/30', 'PQ 35/35', 'PQ 40/40', 'PQ 50/50']
    factor core_dimension_x / core_dimension_y = [1.46, 1.39, 1.45, 1.35, 1.44, 1.56]
    mean over the factors = 1.45
    max derivation / mean = 1.07 (< 7% accuracy)

    :return: Dict including core_h, core_inner_diameter, window_h, window_w
    :rtype: dict
    """
    core_dict = {}

    # -----------------------
    # PQ Cores
    # -----------------------

    core_dict["PQ 16/11.6"] = {
        "core_h": 11.6e-3,
        "core_inner_diameter": 7e-3,
        "window_h": 6.7e-3,
        "window_w": 3.7e-3,
        "core_dimension_x": 16.4e-3,
        "core_dimension_y": 11.2e-3,
    }
    core_dict["PQ 20/16"] = {
        "core_h": 16.2e-3,
        "core_inner_diameter": 8.8e-3,
        "window_h": 10.3e-3,
        "window_w": (18 - 8.8) / 2 * 1e-3,
        "core_dimension_x": 20.5e-3,
        "core_dimension_y": 14.0e-3,
    }
    core_dict["PQ 20/20"] = {
        "core_h": 20.2e-3,
        "core_inner_diameter": 8.8e-3,
        "window_h": 14.3e-3,
        "window_w": 4.6e-3,
        "core_dimension_x": 20.5e-3,
        "core_dimension_y": 14.0e-3,
    }
    core_dict["PQ 26/20"] = {
        "core_h": 20.16e-3,
        "core_inner_diameter": 12e-3,
        "window_h": 11.5e-3,
        "window_w": 5.25e-3,
        "core_dimension_x": 26.5e-3,
        "core_dimension_y": 19.0e-3,
    }
    core_dict["PQ 26/25"] = {
        "core_h": 24.76e-3,
        "core_inner_diameter": 12e-3,
        "window_h": 16.1e-3,
        "window_w": (22.5 - 12) / 2 * 1e-3,
        "core_dimension_x": 26.5e-3,
        "core_dimension_y": 19.0e-3,
    }
    core_dict["PQ 32/20"] = {
        "core_h": 20.5e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 11.5e-3,
        "window_w": (27.5 - 13.45) / 2 * 1e-3,
        "core_dimension_x": 32.0e-3,
        "core_dimension_y": 22.0e-3,
    }
    core_dict["PQ 32/30"] = {
        "core_h": 30.35e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 21.3e-3,
        "window_w": (27.5 - 13.45) / 2 * 1e-3,
        "core_dimension_x": 32.0e-3,
        "core_dimension_y": 22.0e-3,
    }
    core_dict["PQ 35/35"] = {
        "core_h": 34.8e-3,
        "core_inner_diameter": 14.35e-3,
        "window_h": 25e-3,
        "window_w": (32 - 14.35) / 2 * 1e-3,
        "core_dimension_x": 35.1e-3,
        "core_dimension_y": 26.0e-3,
    }
    core_dict["PQ 40/30"] = {
        "core_h": 30.3e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 20e-3,
        "window_w": (37 - 14.9) / 2 * 1e-3,
        "core_dimension_x": 40.5e-3,
        "core_dimension_y": 28.0e-3,
    }
    core_dict["PQ 40/40"] = {
        "core_h": 39.8e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 29.5e-3,
        "window_w": (37 - 14.9) / 2 * 1e-3,
        "core_dimension_x": 40.5e-3,
        "core_dimension_y": 28.0e-3,
    }
    core_dict["PQ 50/40"] = {
        "core_h": 40e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 26.1e-3,
        "window_w": (44 - 20) / 2 * 1e-3,
        "core_dimension_x": 50.0e-3,
        "core_dimension_y": 32.0e-3,
    }
    core_dict["PQ 50/50"] = {
        "core_h": 50e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 36.1e-3,
        "window_w": (44 - 20) / 2 * 1e-3,
        "core_dimension_x": 50.0e-3,
        "core_dimension_y": 32.0e-3,
    }
    core_dict["PQ 65/60"] = {
        "core_h": 60e-3,
        "core_inner_diameter": 26e-3,
        "window_h": 42e-3,
        "window_w": (65 - 26) / 2 * 1e-3,
        "core_dimension_x": 65.0e-3,
        "core_dimension_y": 45.0e-3,
    }
    # -----------------------
    # PM Cores
    # -----------------------

    core_dict["PM 114/93"] = {
        "core_h": 93 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(43e-3, 5.4e-3),
        "window_h": 63 * 1e-3,
        "window_w": (88 - 43) / 2 * 1e-3,
    }
    core_dict["PM 50/39"] = {
        "core_h": 39 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(20e-3, 5.4e-3),
        "window_h": 26.4 * 1e-3,
        "window_w": (39 - 20) / 2 * 1e-3,
    }
    core_dict["PM 62/49"] = {
        "core_h": 49 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(25.5e-3, 5.4e-3),
        "window_h": 33.4 * 1e-3,
        "window_w": (48.8 - 25.5) / 2 * 1e-3,
    }
    core_dict["PM 74/59"] = {
        "core_h": 59 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(29.5e-3, 5.4e-3),
        "window_h": 40.7e-3,
        "window_w": (57.5 - 29.5) / 2 * 1e-3,
    }
    core_dict["PM 87/70"] = {
        "core_h": 70 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(31.7e-3, 8.5e-3),
        "window_h": 48 * 1e-3,
        "window_w": (67.1 - 31.7) / 2 * 1e-3,
    }
    return core_dict

def bobbin_database() -> dict:
    """
    Return a dictionary containing bobbin dimensions for various core structures.

    :return: Dictionary containing bobbin dimensions from Datasheet.
    :rtype: Dict
    """
    bobbin_dict = {}

    # -----------------------
    # PQ Bobbins
    # -----------------------

    # bobbin_dict["PQ 16/11.6"] = {
    #     "bobbin_h": None,
    #     "bobbin_inner_diameter": None,
    #     "bobbin_window_h": None,
    #     "bobbin_window_w": None
    # }
    bobbin_dict["PQ 20/16"] = {
        "bobbin_h": 18.3e-3,
        "bobbin_inner_diameter": 10.95e-3,
        "bobbin_window_h": 9.8e-3,
        "bobbin_window_w": (23.3 - 9.15) / 2 * 1e-3
    }
    bobbin_dict["PQ 20/20"] = {
        "bobbin_h": 22.3e-3,
        "bobbin_inner_diameter": 10.9e-3,
        "bobbin_window_h": 13.9e-3,
        "bobbin_window_w": (23.3 - 9.2) / 2 * 1e-3
    }
    bobbin_dict["PQ 26/25"] = {
        "bobbin_h": 29.3e-3,
        "bobbin_inner_diameter": 14.2e-3,
        "bobbin_window_h": 15.5e-3,
        "bobbin_window_w": (26.5 - 14.2) / 2 * 1e-3
    }
    bobbin_dict["PQ 32/20"] = {
        "bobbin_h": 18.8e-3,
        "bobbin_inner_diameter": 16.2e-3,
        "bobbin_window_h": 10.7e-3,
        "bobbin_window_w": (32.3 - 16.2) / 2 * 1e-3,
    }
    bobbin_dict["PQ 32/30"] = {
        "bobbin_h": 32.8e-3,
        "bobbin_inner_diameter": 15.9e-3,
        "bobbin_window_h": 20.7e-3,
        "bobbin_window_w": (26.4 -15.9) / 2 * 1e-3
    }
    bobbin_dict["PQ 40/40"] = {
        "bobbin_h": 45.3e-3,
        "bobbin_inner_diameter": 15.55e-3,
        "bobbin_window_h": 28.75e-3,
        "bobbin_window_w": (40.3 - 15.55) / 2 * 1e-3
    }
    bobbin_dict["PQ 50/50"] = {
        "bobbin_h": 51.5e-3,
        "bobbin_inner_diameter": 23.2e-3,
        "bobbin_window_h": 35.2e-3,
        "bobbin_window_w": (51.3 - 23.2) / 2 * 1e-3
    }
    bobbin_dict["PQ 65/60"] = {
        "bobbin_h": 65.5e-3,
        "bobbin_inner_diameter": 27.3e-3,
        "bobbin_window_h": 40.7e-3,
        "bobbin_window_w": (66.5 - 27.3) / 2 * 1e-3
    }

    # -----------------------
    # PM Bobbins
    # -----------------------

    bobbin_dict["PM 114/93"] = {
        "bobbin_h": 91 * 1e-3,
        "bobbin_inner_diameter": 44e-3,
        "bobbin_window_h": 62.3e-3,
        "bobbin_window_w": (91 - 44) / 2 * 1e-3,
     }
    bobbin_dict["PM 50/39"] = {
        "bobbin_h": 33.8 * 1e-3,
        "bobbin_inner_diameter": 20.4e-3,
        "bobbin_window_h": 25.9e-3,
        "bobbin_window_w": (38.5 - 20.4) / 2 * 1e-3,
     }
    bobbin_dict["PM 62/49"] = {
        "bobbin_h": 48 * 1e-3,
        "bobbin_inner_diameter": 25.7e-3,
        "bobbin_window_h": 33.0e-3,
        "bobbin_window_w": (48 - 25.7) / 2 * 1e-3,
    }
    bobbin_dict["PM 74/59"] = {
        "bobbin_h": 57.8 * 1e-3,
        "bobbin_inner_diameter": 30e-3,
        "bobbin_window_h": 40.3e-3,
        "bobbin_window_w": (57.3 - 30) / 2 * 1e-3,
    }
    bobbin_dict["PM 87/70"] = {
        "bobbin_h": 68.2 * 1e-3,
        "bobbin_inner_diameter": 32.5e-3,
        "bobbin_window_h": 47.2e-3,
        "bobbin_window_w": (66.2 - 32.5) / 2 * 1e-3,
    }
    return bobbin_dict

def insulation_materials_database() -> dict:
    """
    Return insulation properties for different type of materials

    :return: Dict including insulation parameters
    :rtype: dict
    """
    # To see the shape and the properties of these materials of the wire insulation, review this website:
    # https://www.awcwire.com/customersupport/techinfo/insulation-materials?srsltid=AfmBOoqXkVrB6ITF-R9nL_UlgGVpzX2xB2ENjMNQQHmczLRcf0-y6YwG

    insulation_materials = {
        # wire_insulation materials
        "wire_insulation": {
            # Plastic materials
            "plastic_insulation": {
                # PVC
                # 1.a PVC (pure)
                "Polyvinyl Chloride (PVC)": {"dielectric_constant": 4.0,
                                                       "thermal_conductivity": None,
                                                       "max_temperature": None},
                # """Reference: Huang, J., Zhang, X., Liu, R., Ding, Y., & Guo, D. (2023). Polyvinyl chloride-based dielectric elastomer with high permittivity and
                #     # low viscoelasticity for actuation and sensing. Nature communications, 14(1), 1483."""

                # 1.b Semi-Rigid PVC (SR-PVC). it is from 2.7 to 6.5
                "Semi-Rigid PVC (SR-PVC)": {"dielectric_constant": 3.6,
                                                       "thermal_conductivity": None,
                                                       "max_temperature": None},
                # 1.c Plenum Polyvinyl Chloride (Plenum PVC)
                "Plenum Polyvinyl Chloride (Plenum PVC)": {"dielectric_constant": 3.5,
                                                     "thermal_conductivity": None,
                                                     "max_temperature": None},
                #     """Reference: https://www.anixter.com/content/dam/Anixter/Guide/7H0011X0_W&C_Tech_Handbook_Sec_03.pdf"""

                # 2. Polyethylene (PE): it ranges from 2.7 to 2.8 (30% glass fiber)
                "Polyethylene (PE)": {"dielectric_constant": 2.7,
                                                    "thermal_conductivity": None,
                                                    "max_temperature": None},
                # """Reference: https://passive-components.eu/what-is-dielectric-constant-of-plastic-materials/"""
                #     """Reference:
                #      https://www.awcwire.com/customersupport/techinfo/insulation-materials#:~:text=Semi%2DRigid%20PVC%20(SR%2DPVC)%20is%20mainly%20used,
                #      degrees%20Celsius%2C%20300%20volts)."""

                # 3.Polypropylene (PP) 10-20 glass fiber 2.2 - 2.3
                "Polypropylene (PP)": {"dielectric_constant": 2.3,
                                               "thermal_conductivity": None,
                                               "max_temperature": None},
                # """Reference:
                #          https://www.awcwire.com/customersupport/techinfo/insulation-materials#:~:text=Semi%2DRigid%20PVC%20(SR%2DPVC)%20is%20mainly%20used,
                #          degrees%20Celsius%2C%20300%20volts)."""

                # 4.Polyurethane (PUR): the permittivity differs from 1.065 to 3.35 based on the density (kg/m^3). (need to be reviewed)
                "Polyurethane (PUR)": {"dielectric_constant": 3.35,
                                               "thermal_conductivity": None,
                                               "max_temperature": None},
                # """Reference: Beverte, I. (2025). Investigation of the Partial Permittivity of Rigid Polyurethane Foams by a Circular One-Side-Access Capacitive Sensor.
                #      Polymers, 17(5), 602."""
                #     """https://www.anixter.com/content/dam/Anixter/Guide/7H0011X0_W&C_Tech_Handbook_Sec_03.pdf"""

                # 5. Chlorinated Polyethylene (CPE)
                "Chlorinated Polyethylene (CPE)": {"dielectric_constant": 2.3,
                                               "thermal_conductivity": None,
                                               "max_temperature": None},
                # """https://www.anixter.com/content/dam/Anixter/Guide/7H0011X0_W&C_Tech_Handbook_Sec_03.pdf"""

                # 6. Nylon : 3.2 - 5
                "Nylon": {"dielectric_constant": 5,
                                       "thermal_conductivity": None,
                                       "max_temperature": None},
                #     """https://www.anixter.com/content/dam/Anixter/Guide/7H0011X0_W&C_Tech_Handbook_Sec_03.pdf"""
            },
            # Rubber Materials
            "rubber_insulation": {
                # 1. Thermoplastic Rubber (TPR) 3.30 - 5.10. It can be called Thermoplastic Elastomer (TPE)
                "Thermoplastic Rubber(TPR)": {"dielectric_constant": 5.10,
                                               "thermal_conductivity": None,
                                               "max_temperature": None},
                #     """Reference: https://www.matweb.com/search/datasheet.aspx?matguid=0619837e5f584a1f8c5e6f692952898a"""

                # 2. Neoprene (Polychloroprene): 4-6.7
                "Neoprene (Polychloroprene)": {"dielectric_constant": 6.7,
                                               "thermal_conductivity": None,
                                               "max_temperature": None},
                #     """Reference: https://hep.physics.illinois.edu/home/serrede/p435/lecture_notes/dielectric_constants.pdf"""

                # 3. Styrene-Butadiene Rubber (SBR): 2.5 - 3
                "Styrene-Butadiene Rubber (SBR)": {"dielectric_constant": 3.0,
                                                   "thermal_conductivity": None,
                                                   "max_temperature": None},
                # """Reference: https://www.azom.com/properties.aspx?ArticleID=1844"""

                # 4. Silicone: 2.9 - 4
                "Silicone": {"dielectric_constant": 4.0,
                               "loss_tangent": None,
                               "thermal_conductivity": None,
                               "max_temperature": None},
                # """Reference: https://en.wikipedia.org/wiki/Relative_permittivity"""

                # 5. Fiberglass 3.0 - 4.0
                "Fiberglass": {"dielectric_constant": 4.0,
                                 "thermal_conductivity": None,
                                 "max_temperature": None},
                #     """Reference: https://passive-components.eu/what-is-dielectric-constant-of-plastic-materials/"""

                # 6. Ethylene Propylene Rubber (EPR): 2.4 and  can reach 4
                "Ethylene Propylene Rubber (EPR)": {"dielectric_constant": 2.4,
                                                   "thermal_conductivity": None,
                                                   "max_temperature": None},
                # """Reference: https://passive-components.eu/what-is-dielectric-constant-of-plastic-materials/"""

                # 7. Rubber ( refers to natural rubber and SBR compounds.)
                # 7.a: natural rubber: 2.7 – 4.0 (low freq); 2.4–2.7 (GHz range)
                "Natural Rubber": {"dielectric_constant": 2.7,
                                        "loss_tangent": None,
                                        "thermal_conductivity": None,
                                        "max_temperature": None},
                # """Reference: Al-Hartomy, O. A., Al-Ghamdi, A., Dishovsky, N., Shtarkova, R., Iliev, V., Mutlay, I., & El-Tantawy, F. (2012).
                #    Dielectric and microwave properties of natural rubber based nanocomposites containing graphene."""

                # 7.b: SBR: 2.5 to 3 (low freq); up to 6.6 (Ghz freq)
                "Rubber (SBR)": {"dielectric_constant": 3,
                                           "thermal_conductivity": None,
                                           "max_temperature": None},
                # """Reference: Gunasekaran,S.,Natarajan, R. K., Kala, A., & Jagannathan, R. (2008).
                # Dielectric studies of some rubber materials at microwave frequencies."""

                # 8. Chlorosulfonated Polyethylene (CSPE): Measured dielectric constant: 8-10
                "Chlorosulfonated Polyethylene (CSPE)": {"dielectric_constant": 8.5,
                                                            "thermal_conductivity": None,
                                                            "max_temperature": None},
                # """Reference: Ganguly, S., & Das, N. C. (2015). Chlorosulphonated polyethylene and its composites for electronic applications.
                #  In Flexible and stretchable electronic composites (pp. 229-259). Cham: Springer International Publishing."""

                # 9. Chlorosulfonated Polyethylene (CSPE): Measured dielectric constant: 8-10
                "Chlorosulfonated Polyethylene (CSPE)": {"dielectric_constant": 8.5,
                                                                 "thermal_conductivity": None,
                                                                 "max_temperature": None},
                # """Reference: Ganguly, S., & Das, N. C. (2015). Chlorosulphonated polyethylene and its composites for electronic applications.
                #  In Flexible and stretchable electronic composites (pp. 229-259). Cham: Springer International Publishing."""
            },
            # Fluoropolymer Insulation Types
            "fluoropolymer_insulation": {
                # 1. Perfluoroalkoxy (PFA):  dielectric constant: 2.06 to 2.10:
                "Perfluoroalkoxy (PFA)": {"dielectric_constant": 2.06,
                                         "thermal_conductivity": None,
                                         "max_temperature": None},
                # """Reference: https://adtech.co.uk/application/files/1816/0500/0871/Adtech_PFA_General_Properties_2020.pdf"""
                # """Reference: https://www.fluorotherm.com/technical-information/materials-overview/pfa-properties/"""

                # 2. Polytetrafluoroethylene (PTFE):  dielectric constant: 2.12 to 2.01
                "Polytetrafluoroethylene (PTFE)": {"dielectric_constant": 2.12,
                                                         "thermal_conductivity": None,
                                                         "max_temperature": None},
                # """Reference: Li, L., Bowler, N., Kessler, M. R., & Yoon, S. H. (2010). Dielectric response of PTFE and ETFE wiring insulation to thermal exposure.
                #  IEEE Transactions on Dielectrics and Electrical Insulation, 17(4), 1234-1241."""

                # 3. Fluorinated Ethylene Propylene (FEP):  dielectric constant is about 2.2
                "Fluorinated Ethylene Propylene (FEP)": {"dielectric_constant": 2.2,
                                                          "thermal_conductivity": None,
                                                          "max_temperature": None},
                # """Reference: Lv, X., Yv, J., Wang, X., & Huang, P. (2022). Flexible low dielectric polyimide/fluorinated ethylene propylene composite films for
                #  flexible integrated circuits. Polymer Science, Series B, 64(2), 219-228.."""

                # 4. Ethylene Tetrafluoroethylene (ETFE) :  dielectric constant is about 2.2 to 2.6
                "Ethylene Tetrafluoroethylene (ETFE)": {"dielectric_constant": 2.6,
                                                                        "thermal_conductivity": None,
                                                                        "max_temperature": None},
                # """Reference: Wang, M., He, Y., Yang, X., Hou, X., Li, W., Tan, S., & Zhang, Z. (2024). Optimizing thermal and dielectric properties of
                #  ethylene-tetrafluoroethylene (ETFE)/h-BN composites via interface engineering: activation of C–F bonds on ETFE for surface grafting.
                #   Journal of Materials Chemistry A, 12(45), 31424-31431."""

                # 5. Ethylenechlorotrifluoroethylene (ECTFE) :  dielectric constant is 2.5
                # Note 4. and 5. have approximately the same properties
                "Ethylenechlorotrifluoroethylene (ECTFE)": {"dielectric_constant": 2.5,
                                                           "thermal_conductivity": None,
                                                           "max_temperature": None},
                # """Reference: https://www.polyfluor.nl/assets/files/datasheet-ectfe-uk.pdf"""

                # 6. Polyvinylidene Fluoride (PVDF) :  dielectric constant is 8.5 at 1MHz
                "Polyvinylidene Fluoride (PVDF)": {"dielectric_constant": 8.5,
                                                   "thermal_conductivity": None,
                                                   "max_temperature": None},
                # """Reference: https://www.ipolymer.com/pdf/PVDF.pdf"""

                # 7.Thermoplastic Elastomers (TPE) :  dielectric constant is 3.3 to 5.1
                "Thermoplastic Elastomers (TPE)": {"dielectric_constant": 4.5,
                                                  "thermal_conductivity": None,
                                                  "max_temperature": None},
                # """https://www.matweb.com/search/datasheet.aspx?matguid=0619837e5f584a1f8c5e6f692952898a&"""

            }
        },

        # Bobbin materials is based on these references:
        # 1. https://www.cosmocorp.com/docs/en/cosmo-bcat-en-mdres.pdf
        # 2. https://pearl-hifi.com/06_Lit_Archive/06_Mchy_Methods_Catalogues/Coil_Winding/Transformer_Bobbin_and_Core_Selection.pdf

        "core_insulation": {
            "bobbins": {
                "Thermoplastic": {
                    "Polyamide (Nylon 66)": {
                        "dielectric_constant": 3.8,
                        "thermal_conductivity": 0.25,
                        "max_temperature": 130
                    },
                    "Polybutylene Terephthalate (PBT)": {
                        "dielectric_constant": 3.7,
                        "thermal_conductivity": 0.25,
                        "max_temperature": 155
                    },
                    "Polyphenylene Sulfide (PPS)": {
                        "dielectric_constant": 3.8,
                        "thermal_conductivity": 0.30,
                        "max_temperature": 200
                    },
                    "Liquid Crystal Polymer (LCP)": {
                        "dielectric_constant": 3.6,
                        "thermal_conductivity": 0.50,
                        "max_temperature": 240
                    },
                    "Polyethylene Terephthalate (PET)": {
                        "dielectric_constant": 3.6,
                        "thermal_conductivity": 0.15,
                        "max_temperature": 150
                    }
                },
                "Thermoset": {
                    "Diallyl Phthalate (DAP)": {
                        "dielectric_constant": 4.4,
                        "thermal_conductivity": 0.20,
                        "max_temperature": 200
                    },
                    "Phenolic": {
                        "dielectric_constant": 4.5,
                        "thermal_conductivity": 0.20,
                        "max_temperature": 220
                    }
                }
            }
        },

        # The kapton has different types, but the dielectric constant is approximately around 3.4-3.5
        # Reference: https://www.dupont.com/content/dam/electronics/amer/us/en/electronics/public/documents/en/EI-10167-Kapton-General-Specifications.pdf
        "film_insulation": {
            "Kapton": {
                "dielectric_constant": 3.5,
                "thermal_conductivity": None,
                "max_temperature": None,
            },
            "Nomex 410": {
                "dielectric_constant": 1.6,
                "thermal_conductivity": None,
                "max_temperature": None,
            },
            "Mylar": {
                "dielectric_constant": 3.2,
                "thermal_conductivity": None,
                "max_temperature": None,
            },
            "PET / bOPET": {
                "dielectric_constant": 3.3,
                "thermal_conductivity": None,
                "max_temperature": None,
            },
            "PVC Tape": {
                "dielectric_constant": 4.0,
                "max_temperature": None,
                "thermal_conductivity": None
            }
        }
    }

    return insulation_materials


def litz_database() -> dict:
    """
    Return litz parameters for defined litz wires.

    :return: Dict including litz parameters like strand_numbers, strand_radii and conductor_radii
    :rtype: dict
    """
    litz_dict = {}

    litz_dict["1.5x105x0.1"] = {"strands_numbers": 105,
                                "strand_radii": 0.1e-3 / 2,
                                "conductor_radii": 1.5e-3 / 2,
                                "ff": None,
                                "manufacturer": "PACK",
                                "material_number": "",
                                "litz": "RUPALIT V155",
                                "insulation": "textile"}
    litz_dict["1.4x200x0.071"] = {"strands_numbers": 200,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 1.4e-3 / 2,
                                  "ff": None,
                                  "manufacturer": "PACK",
                                  "material_number": "",
                                  "litz": "RUPALIT V155",
                                  "insulation": "textile"}
    litz_dict["2.0x405x0.071"] = {"strands_numbers": 405,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 2.0e-3 / 2,
                                  "ff": None,
                                  "manufacturer": "",
                                  "material_number": "",
                                  "litz": "",
                                  "insulation": "unknown blue plastic"}
    litz_dict["2.0x800x0.05"] = {"strands_numbers": 800,
                                 "strand_radii": 0.05e-3 / 2,
                                 "conductor_radii": 2e-3 / 2,
                                 "ff": None,
                                 "manufacturer": "Elektrisola",
                                 "material_number": "12104184",
                                 "litz": "",
                                 "insulation": ""
                                 }
    litz_dict["1.1x60x0.1"] = {"strands_numbers": 60,
                               "strand_radii": 0.1e-3 / 2,
                               "conductor_radii": 1.1e-3 / 2,
                               "ff": None,
                               "manufacturer": "PACK",
                               "material_number": "",
                               "litz": "RUPALIT V155",
                               "insulation": "textile"
                               }
    litz_dict["1.35x200x0.071"] = {"strands_numbers": 200,
                                   "strand_radii": 0.071e-3 / 2,
                                   "conductor_radii": 1.35e-3 / 2,
                                   "ff": None,
                                   "manufacturer": "PACK",
                                   "material_number": "",
                                   "litz": "RUPALIT V155",
                                   "insulation": "textile"}

    litz_dict["3.2x2100x0.05"] = {"strands_numbers": 2100,
                                  "strand_radii": 0.05e-3 / 2,
                                  "conductor_radii": 3.2e-3 / 2,
                                  "ff": None,
                                  "manufacturer": "PACK",
                                  "material_number": "AB21220373",
                                  "litz": "RUPALIT V155",
                                  "insulation": "textile"
                                  }

    litz_dict["4.6x2160x0.071"] = {"strands_numbers": 2160,
                                   "strand_radii": 0.071e-3 / 2,
                                   "conductor_radii": 4.6e-3 / 2,
                                   "ff": None,
                                   "manufacturer": "PACK",
                                   "material_number": "AB21225497",
                                   "litz": "RUPALIT V155",
                                   "insulation": "textile"
                                   }

    litz_dict["2.9x1200x0.06"] = {"strands_numbers": 1200,
                                  "strand_radii": 0.06e-3 / 2,
                                  "conductor_radii": 2.9e-3 / 2,
                                  "ff": None,
                                  "manufacturer": "Elektrisola",
                                  "material_number": "",
                                  "litz": "",
                                  "insulation": "unknown plastic"}

    litz_dict["2.6x1000x0.06"] = {"strands_numbers": 1000,
                                  "strand_radii": 0.06e-3 / 2,
                                  "conductor_radii": 2.6e-3 / 2,
                                  "ff": None,
                                  "manufacturer": "Elektrisola",
                                  "material_number": "",
                                  "litz": "",
                                  "insulation": "unknown plastic"}

    litz_dict["1.8x512x0.05"] = {"strands_numbers": 512,
                                 "strand_radii": 0.05e-3 / 2,
                                 "conductor_radii": 1.8e-3 / 2,
                                 "ff": None,
                                 "manufacturer": "PACK",
                                 "material_number": "AB21217207",
                                 "litz": "RUPALIT Safety VB155",
                                 "insulation": "3 layers Mylar"}

    litz_dict["2.3x600x0.071"] = {"strands_numbers": 600,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 2.3e-3 / 2,
                                  "ff": None,
                                  "manufacturer": "PACK",
                                  "material_number": "AB21220522",
                                  "litz": "RUPALIT Safety Profil V155",
                                  "insulation": "3 layers Mylar"}

    litz_dict["2.8x400x0.1"] = {"strands_numbers": 400,
                                "strand_radii": 0.1e-3 / 2,
                                "conductor_radii": 2.8e-3 / 2,
                                "ff": None,
                                "manufacturer": "PACK",
                                "material_number": "AB21222210",
                                "litz": "RUPALIT Safety V155",
                                "insulation": "3 layers Mylar"}

    litz_dict["1.71x140x0.1"] = {"strands_numbers": 140,
                                 "strand_radii": 0.1e-3 / 2,
                                 "conductor_radii": 1.71e-3 / 2,
                                 "ff": None,
                                 "manufacturer": "",
                                 "material_number": "",
                                 "litz": "",
                                 "insulation": ""}

    litz_dict["1.7x500x0.06"] = {"strands_numbers": 500,
                                 "strand_radii": 0.06e-3 / 2,
                                 "conductor_radii": 1.7e-3 / 2,
                                 "ff": None,
                                 "manufacturer": "",
                                 "material_number": "",
                                 "litz": "",
                                 "insulation": ""}

    return litz_dict


def wire_material_database() -> dict[str, WireMaterial]:
    """
    Return wire materials e.g. copper, aluminum in a dictionary.

    :return: Dict with materials and conductivity
    :rtype: dict
    """
    wire_material = {}

    wire_material["Copper"] = WireMaterial(
        name="copper",
        sigma=5.8e7,
        temperature=25,
        temperature_coefficient=3.9e-3,
        thermal_conductivity=400,
        volumetric_mass_density=8920,
    )

    wire_material["Aluminum"] = WireMaterial(
        name="aluminum",
        sigma=3.7e7,
        temperature=25,
        temperature_coefficient=3.9e-3,
        thermal_conductivity=235,
        volumetric_mass_density=2699,
    )

    return wire_material


def conductivity_temperature(material: str, temperature: float) -> float:
    """
    Calculate the conductivity for a certain temperature of the material.

    :param material: material name, e.g. "copper"
    :type material: str
    :param temperature: temperature in °C
    :type temperature: float
    :return: conductivity of material at given temperature
    :rtype: float
    """
    material_from_database = wire_material_database()[material]

    sigma_database = material_from_database.sigma
    temperature_database = material_from_database.temperature
    temperature_coefficient_database = material_from_database.temperature_coefficient

    resistance_temperature = 1 / sigma_database * (
        1 + temperature_coefficient_database * (temperature - temperature_database))
    sigma_temperature = 1 / resistance_temperature

    return sigma_temperature


def create_folders(*args: str) -> None:
    """
    Create folders for every given folder path (if it does not exist).

    :param args: Folder names
    :type args: str
    """
    for folder in list(args):
        if not os.path.exists(folder):
            os.mkdir(folder)


def cost_material_database() -> dict:
    """
    Return costs for core and winding. This is split in material and fabrication costs.

    Both, material and fabrication costs have a euro_per_kilogram and a euro_per_unit (fix costs) price.

    Source: R. Burkart and J. Kolar 'Component Cost Models for Multi-Objective Optimizations of #
    Switched-Mode Power Converter' 2013.

    These are outdated prices (year 2013). Update needed in the future.
    """
    cost_database = {}
    cost_database["core_euro_per_kilogram"] = {"ferrite": 5.5,
                                               "amorphous": 16,
                                               "nanocristalline": 23,
                                               "high_si_steel": 12,
                                               "goes": 2.5}
    cost_database["winding_material_euro_per_kilogram"] = {ConductorType.RoundSolid.name: 10,
                                                           "flat": 10,
                                                           ConductorType.RectangularSolid.name: 20,
                                                           ConductorType.RoundLitz.name: -1}

    cost_database["winding_material_euro_per_unit"] = {ConductorType.RoundSolid.name: 1,
                                                       "flat": 2,
                                                       ConductorType.RectangularSolid.name: 2,
                                                       ConductorType.RoundLitz.name: 1}

    cost_database["winding_fabrication_euro_per_kilogram"] = {ConductorType.RoundSolid.name: 7,
                                                              "flat": 21,
                                                              ConductorType.RectangularSolid.name: 14,
                                                              ConductorType.RoundLitz.name: 7}
    cost_database["winding_fabrication_euro_per_unit"] = {ConductorType.RoundSolid.name: 2,
                                                          "flat": 4,
                                                          ConductorType.RectangularSolid.name: 2.5,
                                                          ConductorType.RoundLitz.name: 2}

    cost_database["winding_material_euro_per_kilogram_for_litz"] = {"sigma_numerator": 15,
                                                                    "sigma_denominator": 0.45}

    cost_database["gross_margin"] = 0.25

    return cost_database


def pm_core_inner_diameter_calculator(inner_core_diameter: float, hole_diameter: float) -> np.array:
    """
    Calculate the effective inner core diameter without the hole often used in PM-cores.

    :param inner_core_diameter: inner core diameter
    :type inner_core_diameter: float
    :param hole_diameter: hole diameter
    :type hole_diameter: float
    :return: effective inner core diameter without hole
    :rtype: np.array
    """
    area_inner_core_inner_diameterithout_hole = (inner_core_diameter / 2) ** 2 * np.pi
    area_hole = (hole_diameter / 2) ** 2 * np.pi
    area_total = area_inner_core_inner_diameterithout_hole - area_hole

    return np.around(2 * np.sqrt(area_total / np.pi), decimals=4)


def install_pyfemm_if_missing() -> None:
    """Installs femm-software pip package in case of running on Windows machine. Windows users only."""
    required = {'pyfemm'}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed

    if missing:
        logger.info("Missing 'pyfemm' installation.")
        logger.info("Installing 'pyfemm' ...")
        python = sys.executable
        subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
        logger.info("'pyfemm' is now installed!")

def litz_calculate_number_strands(n_layers: int) -> int:
    """
    Return the number of strands in a hexagonal litz winding with a specified number of layers (n_layers).

    CAUTION: Zero number of layers corresponds to a single strand.

    :param n_layers: number of litz_layers
    :type n_layers: int

    :return: number of strands in a litz wire
    :rtype: int

    """
    return 3 * (n_layers + 1) ** 2 - 3 * (n_layers + 1) + 1


def litz_calculate_number_layers(n_strands: int) -> int:
    """
    Return the number of layers in a hexagonal litz winding with a specified number of strands (n_strands).

    .. note:: Zero number of layers corresponds to a single strand.

    :param n_strands: Number of strands in a litz
    :type n_strands: int

    :return: number of layers for a litz
    :rtype: int
    """
    return np.sqrt(0.25 + (n_strands - 1) / 3) - 0.5


def fft(period_vector_t_i: npt.ArrayLike, sample_factor: int = 1000, plot: str = 'no', mode: str = 'rad',
        f0: float | None = None, title: str = 'ffT', filter_type: str = 'factor',
        filter_value_factor: float = 0.01, filter_value_harmonic: int = 100,
        figure_size: tuple = None, figure_directory: str = None) -> npt.NDArray[list]:
    """
    Calculate the FFT for a given input signal. Input signal is in vector format and should include one period.

    Output vector includes only frequencies with amplitudes > 1% of input signal

    :Minimal Example:

    >>> import femmt as fmt
    >>> import numpy as np
    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> out = fmt.fft(example_waveform, plot='yes', mode='rad', f0=25000, title='ffT input current')

    :param period_vector_t_i: numpy-array [[time-vector[,[current-vector]]. One period only
    :type period_vector_t_i: np.array
    :param sample_factor: f_sampling/f_period, defaults to 1000
    :type sample_factor: int
    :param plot: insert anything else than "no" or 'False' to show a plot to visualize input and output
    :type plot: str
    :param mode: 'rad'[default]: full period is 2*pi, 'deg': full period is 360°, 'time': time domain.
    :type mode: str
    :param f0: fundamental frequency. Needs to be set in 'rad'- or 'deg'-mode
    :type f0: float
    :param title: plot window title, defaults to 'ffT'
    :type title: str
    :param filter_type: 'factor'[default] or 'harmonic' or 'disabled'.
    :type filter_type: str
    :param filter_value_factor: filters out amplitude-values below a certain factor of max. input amplitude.
        Should be 0...1, default to 0.01 (1%)
    :type filter_value_factor: float
    :param filter_value_harmonic: filters out harmonics up to a certain number. Default value is 100.
        Note: count 1 is DC component, count 2 is the fundamental frequency
    :type filter_value_harmonic: int
    :param figure_directory: full path with file extension
    :type figure_directory: tuple
    :param figure_size: None for auto-fit; fig_size for matplotlib (width, length)
    :type figure_size: tuple

    :return: numpy-array [[frequency-vector],[amplitude-vector],[phase-vector]]
    :rtype: npt.NDArray[list]
    """
    # check for correct input parameter
    if (mode == 'rad' or mode == 'deg') and f0 is None:
        raise ValueError("if mode is 'rad' or 'deg', a fundamental frequency f0 must be set")
    # check for input is list. Convert to numpy-array
    if isinstance(period_vector_t_i, list):
        if plot != 'no' and plot is not False:
            logger.warning("Input is list, convert to np.array()")
        period_vector_t_i = np.array(period_vector_t_i)

    # first value of time vector must be zero
    if period_vector_t_i[0][0] != 0:
        raise ValueError("Period vector must start with 0 seconds!")

    # mode pre-calculation
    if mode == 'rad':
        period_vector_t_i[0] = period_vector_t_i[0] / (2 * np.pi * f0)
    elif mode == 'deg':
        period_vector_t_i[0] = period_vector_t_i[0] / (360 * f0)
    elif mode != 'time':
        raise ValueError("Mode not available. Choose: 'rad', 'deg', 'time'")

    t = period_vector_t_i[0]
    i = period_vector_t_i[1]

    # fft-function works per default in time domain
    t_interp = np.linspace(0, t[-1], sample_factor)
    i_interp = np.interp(t_interp, t, i)

    f0 = round(1 / t[-1])
    Fs = f0 * sample_factor

    # frequency domain
    f = np.linspace(0, (sample_factor - 1) * f0, sample_factor)
    x = np.fft.fft(i_interp)
    x_mag = np.abs(x) / sample_factor
    phi_rad = np.angle(x)

    f_corrected = f[0:int(sample_factor / 2 + 1)]
    x_mag_corrected = 2 * x_mag[0:int(sample_factor / 2 + 1)]
    x_mag_corrected[0] = x_mag_corrected[0] / 2
    phi_rad_corrected = phi_rad[0:int(sample_factor / 2 + 1)]

    f_out = []
    x_out = []
    phi_rad_out = []
    if filter_type.lower() == 'factor':
        for count, _ in enumerate(x_mag_corrected):
            if x_mag_corrected[count] > filter_value_factor * max(abs(i)):
                f_out.append(f_corrected[count])
                x_out.append(x_mag_corrected[count])
                phi_rad_out.append(phi_rad_corrected[count])
    elif filter_type.lower() == 'harmonic':
        for count, _ in enumerate(x_mag_corrected):
            if count < filter_value_harmonic:
                f_out.append(f_corrected[count])
                x_out.append(x_mag_corrected[count])
                phi_rad_out.append(phi_rad_corrected[count])
    elif filter_type.lower() == 'disabled':
        f_out = f_corrected
        x_out = x_mag_corrected
        phi_rad_out = phi_rad_corrected
    else:
        raise ValueError(
            f"filter_type '{filter_value_harmonic}' not available: Must be 'factor','harmonic' or 'disabled ")

    if plot != 'no' and plot is not False:
        logger.info(f"{title=}")
        logger.info(f"{t[-1]=}")
        logger.info(f"{f0=}")
        logger.info(f"{Fs=}")
        logger.info(f"{sample_factor=}")
        logger.info(f"f_out = {np.around(f_out, decimals=0)}")
        logger.info(f"x_out = {np.around(x_out, decimals=3)}")
        logger.info(f"phi_rad_out = {np.around(phi_rad_out, decimals=3)}")

        reconstructed_signal = 0
        for i_range in range(len(f_out)):
            reconstructed_signal += x_out[i_range] * np.cos(
                2 * np.pi * f_out[i_range] * t_interp + phi_rad_out[i_range])

        fig, [ax1, ax2, ax3] = plt.subplots(num=title, nrows=3, ncols=1, figsize=figure_size)
        ax1.plot(t, i, label='original signal')
        ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')
        ax1.grid()
        ax1.set_title('Signal')
        ax1.set_xlabel('time in s')
        ax1.set_ylabel('Amplitude')
        ax1.legend()
        ax2.stem(f_out, x_out)
        ax2.grid()
        ax2.set_title('ffT')
        ax2.set_xlabel('Frequency in Hz')
        ax2.set_ylabel('Amplitude')

        ax3.stem(f_out, phi_rad_out)
        ax3.grid()
        ax3.set_title('ffT')
        ax3.set_xlabel('Frequency in Hz')
        ax3.set_ylabel('Phase in rad')

        plt.tight_layout()
        if figure_directory is not None:
            plt.savefig(figure_directory, bbox_inches="tight")
        plt.show()

    return np.array([f_out, x_out, phi_rad_out])


def plot_fourier_coefficients(frequency_list: list, amplitude_list: list, phi_rad_list: list,
                              sample_factor: int = 1000, figure_directory: str = None):
    """
    Plot fourier coefficients in a visual figure.

    :param frequency_list: List of frequencies in Hz
    :type frequency_list: list
    :param amplitude_list: List of amplitudes in A
    :type amplitude_list: list
    :param phi_rad_list: List of angles in rad
    :type phi_rad_list: list
    :param sample_factor: sample factor
    :type sample_factor: int
    :param figure_directory: directory of figure to save
    :type figure_directory: str
    """
    # dc and ac handling
    nonzero_frequencies = [f for f in frequency_list if f != 0]
    if nonzero_frequencies:
        time_period = 1 / min(nonzero_frequencies)
    else:
        time_period = 1
    # time_period = 1 / min(frequency_list)

    t_interp = np.linspace(0, time_period, sample_factor)
    reconstructed_signal = 0
    for i_range, _ in enumerate(frequency_list):
        reconstructed_signal += amplitude_list[i_range] * np.cos(
            2 * np.pi * frequency_list[i_range] * t_interp + phi_rad_list[i_range])

    fig, [ax1, ax2, ax3] = plt.subplots(nrows=3, ncols=1)
    ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')
    ax1.grid()
    ax1.set_title('Signal')
    ax1.set_xlabel('time in s')
    ax1.set_ylabel('Amplitude')
    ax1.legend()
    ax2.stem(frequency_list, amplitude_list)
    ax2.grid()
    ax2.set_title('ffT')
    ax2.set_xlabel('Frequency in Hz')
    ax2.set_ylabel('Amplitude')

    ax3.stem(frequency_list, phi_rad_list)
    ax3.grid()
    ax3.set_title('ffT')
    ax3.set_xlabel('Frequency in Hz')
    ax3.set_ylabel('Phase in rad')

    plt.tight_layout()
    if figure_directory is not None:
        plt.savefig(figure_directory, bbox_inches="tight")
    plt.close('all')  # close the figures to remove the warning when you run many figures

    # plt.show()


def compare_fft_list(input_data_list: list, sample_factor: int = 1000, mode: str = 'rad',
                     f0: float | None = None) -> None:
    """
    Generate fft curves from input curves and compare them to each other.

    :Minimal Example:

    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> example_waveform_2 = np.array([[0, 0.55, 3.14, 3.69, 6.28],[-138.37, 257.58, 138.37, -257.58, -138.37]])
    >>> compare_fft_list([example_waveform, example_waveform_2], mode='rad', f0=25000)

    :param input_data_list: list of fft-compatible numpy-arrays [element, element, ... ], each element format
        like [[time-vector[,[current-vector]]. One period only
    :param mode: 'rad'[default]: full period is 2*pi, 'deg': full period is 360°, 'time': time domain.
    :type mode: str
    :param f0: fundamental frequency. Needs to be set in 'rad'- or 'deg'-mode
    :type f0: float
    :param sample_factor: sample factor, defaults to 1000
    :type sample_factor: int
    """
    out = []
    for count, _ in enumerate(input_data_list):
        out.append([fft(input_data_list[count], sample_factor=sample_factor, plot='no', mode=mode, f0=f0)])

    fig, axs = plt.subplots(2, len(input_data_list), sharey=True)
    for count, _ in enumerate(input_data_list):
        axs[0, count].plot(input_data_list[count][0], input_data_list[count][1], label='original signal')
        axs[0, count].grid()
        axs[0, count].set_xlabel('time in s')
        axs[0, count].set_ylabel('Amplitude')
        axs[1, count].stem(out[count][0][0], out[count][0][1])
        axs[1, count].grid()
        axs[1, count].set_xlabel('frequency in Hz')
        axs[1, count].set_ylabel('Amplitude')

        # ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')

    plt.tight_layout()
    plt.show()


def store_as_npy_in_directory(dir_path: str, file_name: str, numpy_data) -> None:
    """
    Store a numpy array in a given directory.

    :param dir_path: directory path
    :type dir_path: str
    :param file_name: file name
    :type file_name: str
    :param numpy_data: numpy array
    :type numpy_data:
    """
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    np.save(dir_path + "/" + file_name, numpy_data)


def get_dicts_with_keys_and_values(data, **kwargs) -> dict:
    """
    Return a list of dictionaries out of a list of dictionaries which contains pairs of the given key(s) and value(s).

    :param data: list of dicts
    :type data: list
    :param kwargs: keys and values in dicts
    """
    invalid_index = []
    for n, dictionary in enumerate(data):
        for key, value in kwargs.items():
            if not (key in dictionary and value == dictionary[key]):
                invalid_index.append(n)
                break
    valid_data = np.delete(data, invalid_index)
    return valid_data


def get_dict_with_unique_keys(data: list[dict], *keys) -> dict:
    """
    Return a dictionary out of a list of dictionaries which contains the given key(s).

    :param data: list of dicts
    :type data: list[dict]
    :param keys: keys in dicts
    :return:
    """
    invalid_index = []
    for n, dictionary in enumerate(data):
        for key in keys:
            if key not in dictionary:
                invalid_index.append(n)
                break
    valid_data = np.delete(data, invalid_index)
    if len(valid_data) != 1:
        warnings.warn("Keyword(s) not unique!", stacklevel=2)
    # Only one dictionary shall survive --> choose element 0
    return valid_data[0]


def find_common_frequencies(frequency_list_1: list, amplitude_list_1: list, phase_list_1_rad_or_deg: list,
                            frequency_list_2: list, amplitude_list_2: list, phase_list_2_rad_or_deg: list) -> list:
    """
    Find common frequencies and returns a list of intersections.

    :param amplitude_list_1: Amplitudes signal 1
    :type amplitude_list_1: list
    :param phase_list_1_rad_or_deg: Phases signal 1, can be degree or rad. return is same as input.
    :type phase_list_1_rad_or_deg: list
    :param frequency_list_1: Frequencies signal 1
    :type frequency_list_1: list
    :param amplitude_list_2: Amplitudes signal 2
    :type amplitude_list_2: list
    :param phase_list_2_rad_or_deg: Phases signal 2, can be degree or rad. return is same as input
    :type phase_list_2_rad_or_deg: list
    :param frequency_list_2: Frequencies signal 2
    :type frequency_list_2: list
    :return: [current_pair_list, phase_pair_list, common_frequency_list]
    :rtype: tuple

    :Example:

    >>> import femmt as fmt
    >>> frequency_1 = [50, 100, 150, 200]
    >>> frequency_2 = [50, 100, 150, 170, 200]
    >>> amplitude_1 = [1, 2, 3, 4]
    >>> amplitude_2 = [5, 6, 7, 8, 9]
    >>> phase_1 = [10, 20, 30, 40]
    >>> phase_2 = [50, 60, 70, 80, 90]
    >>> common_f, common_a, common_phase = fmt.find_common_frequencies(frequency_1, amplitude_1, phase_1,
    >>>     frequency_2, amplitude_2, phase_2)
    :Returns:
    >>> common_f = [200, 50, 100, 150]
    >>> common_a = [[4, 9], [1, 5], [2, 6], [3, 7]]
    >>> common_phase = [[40, 90], [10, 50], [20, 60], [30, 70]]
    """
    common_amplitude_list_1 = []
    common_phase_list_1 = []
    common_amplitude_list_2 = []
    common_phases_list_2 = []

    common_frequency_list = list(set(frequency_list_1).intersection(frequency_list_2))
    logger.info(f"{common_frequency_list=}")
    common_frequency_list.sort()
    logger.info(f"{common_frequency_list=}")

    # Delete the corresponding phases and amplitudes
    if isinstance(amplitude_list_1, list):
        for frequency in common_frequency_list:
            common_amplitude_list_1.append(amplitude_list_1[frequency_list_1.index(frequency)])
            common_phase_list_1.append(phase_list_1_rad_or_deg[frequency_list_1.index(frequency)])
            common_amplitude_list_2.append(amplitude_list_2[frequency_list_2.index(frequency)])
            common_phases_list_2.append(phase_list_2_rad_or_deg[frequency_list_2.index(frequency)])
    elif isinstance(amplitude_list_1, np.ndarray):
        for frequency in common_frequency_list:
            common_amplitude_list_1.append(amplitude_list_1[np.where(frequency_list_1 == frequency)][0])
            common_phase_list_1.append(phase_list_1_rad_or_deg[np.where(frequency_list_1 == frequency)][0])
            common_amplitude_list_2.append(amplitude_list_2[np.where(frequency_list_2 == frequency)][0])
            common_phases_list_2.append(phase_list_2_rad_or_deg[np.where(frequency_list_2 == frequency)][0])
    else:
        warnings.warn("Either a list or a np.ndarray must be provided!", stacklevel=2)

    current_pair_list = list(map(list, zip(common_amplitude_list_1, common_amplitude_list_2)))
    phase_pair_list = list(map(list, zip(common_phase_list_1, common_phases_list_2)))

    return [common_frequency_list, current_pair_list, phase_pair_list]


def sort_out_small_harmonics(frequency_list: list, amplitude_pair_list: list,
                             phase_pair_list_rad_or_deg: list, sort_out_factor: float) -> list:
    """
    Sort out small harmonics from a given fft-output of a signal.

    :param frequency_list: list of input frequencies
    :type frequency_list: list
    :param amplitude_pair_list: list of amplitude pairs
    :type amplitude_pair_list: list
    :param phase_pair_list_rad_or_deg: list of phase pairs (can be rad or degree)
    :type phase_pair_list_rad_or_deg: list
    :param sort_out_factor: sort out factor [0...1]
    :type sort_out_factor: float
    :return: [frequency_list, amplitude_pair_list, phase_pair_list_rad_or_deg]
    :rtype: list
    """
    # Calculate geometric lengths
    amp_tot = np.sqrt(np.sum(np.array(amplitude_pair_list) ** 2, axis=0))
    # amp_tot = np.max(amplitude_pairs, axis=0)

    invalid_index = []
    for n, amplitude_pair in enumerate(amplitude_pair_list):
        if all(amplitude / amp_tot[i] < sort_out_factor for i, amplitude in enumerate(amplitude_pair)):
            invalid_index.append(n)

    phase_pair_list_rad_or_deg = np.delete(phase_pair_list_rad_or_deg, invalid_index, axis=0)
    amplitude_pair_list = np.delete(amplitude_pair_list, invalid_index, axis=0)
    frequency_list = np.delete(frequency_list, invalid_index)

    return [frequency_list, amplitude_pair_list, phase_pair_list_rad_or_deg]


# Reluctance Model [with calculation]

def calculate_cylinder_volume(cylinder_diameter: float, cylinder_height: float):
    """
    Calculate the volume of an ideal cylinder.

    This function is uses e.g. to calculate the volume of the inner core part.
    :param cylinder_height: height of cylinder
    :type cylinder_height: float
    :param cylinder_diameter: diameter of cylinder
    :type cylinder_diameter: float
    :returns: volume
    :rtype: float
    """
    return (cylinder_diameter / 2) ** 2 * np.pi * cylinder_height


def create_physical_group(dim: int, entities: int, name: str):
    """
    Create a physical group, what is used inside ONELAB.

    :param dim: dim inside onelab
    :type dim: int
    :param entities: entity inside onelab
    :type entities: int
    :param name: name
    :type name: str
    """
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)

    return tag


def visualize_simulation_results(simulation_result_file_path: str, store_figure_file_path: str, show_plot: bool = True) -> None:
    """
    Visualize the simulation results by a figure.

    :param simulation_result_file_path: file path for the simulation results
    :type simulation_result_file_path: str
    :param store_figure_file_path: file path for the figure to store
    :type store_figure_file_path: str
    :param show_plot: True to show the plot
    :type show_plot: bool
    """
    with open(simulation_result_file_path, "r") as fd:
        loaded_results_dict = json.loads(fd.read())

    # Initialize accumulators for cumulative losses and inductances
    cumulative_core_hysteresis = 0
    cumulative_core_eddy = 0
    cumulative_losses = []
    windings_labels = []
    # Determine if this is a single simulation or a sweep
    is_single_simulation = len(loaded_results_dict["single_sweeps"]) == 1

    for index, sweep in enumerate(loaded_results_dict["single_sweeps"]):
        freq = sweep['f']
        loss_core_eddy_current = sweep.get("core_eddy_losses", 0)
        loss_core_hysteresis = sweep.get("core_hyst_losses", 0)

        # Accumulate core losses
        cumulative_core_hysteresis += loss_core_hysteresis
        cumulative_core_eddy += loss_core_eddy_current

        # Plotting for each frequency
        fig, ax = plt.subplots()
        ax.bar(0, loss_core_hysteresis, width=0.35, label='Core Hysteresis Loss')
        ax.bar(0, loss_core_eddy_current, bottom=loss_core_hysteresis, width=0.35, label='Core Eddy Current Loss')

        for i in range(1, 4):
            winding_key = f"winding{i}"
            if winding_key in sweep:
                inductance = sweep[winding_key].get("flux_over_current", [0])[0]
                loss = sweep[winding_key].get("winding_losses", 0)

                if len(cumulative_losses) < i:
                    cumulative_losses.append(loss)
                    windings_labels.append(f"Winding {i}")
                else:
                    cumulative_losses[i - 1] += loss

                # Plot for current frequency
                ax.bar(i, loss, width=0.35, label=f'{windings_labels[i - 1]} Loss at {freq} Hz')

        ax.set_ylabel('Losses in W')
        ax.set_title(f'Loss Distribution at {freq} Hz')
        ax.legend()
        plt.grid(True)

        # Save plot for the current frequency
        base_path, ext = os.path.splitext(store_figure_file_path)
        filename = f"{base_path}_{index}{ext}"
        plt.savefig(filename, bbox_inches="tight")
        plt.close(fig)

    # Plot cumulative results for core and windings
    fig, ax = plt.subplots()
    if is_single_simulation:
        ax.bar(0, cumulative_core_hysteresis, width=0.35, label='Cumulative Core Hysteresis Loss')
        ax.bar(0, cumulative_core_eddy, bottom=cumulative_core_hysteresis, width=0.35, label='Cumulative Core Eddy Current Loss')
    else:
        ax.bar(0, cumulative_core_eddy, width=0.35, label='Cumulative Core Eddy Current Loss')

    for index, loss in enumerate(cumulative_losses):
        ax.bar(index + 1, loss, width=0.35, label=f'{windings_labels[index]} Cumulative Loss')

    ax.set_ylabel('total Losses in W')
    ax.set_title('Loss Distribution in Magnetic Components for all frequencies')
    ax.set_xticks(range(len(windings_labels) + 1))
    ax.set_xticklabels(['Core'] + windings_labels)
    ax.legend()
    plt.grid(True)
    base_path, ext = os.path.splitext(store_figure_file_path)
    cumulative_filename = f"{base_path}_total_freq{ext}"
    plt.savefig(cumulative_filename, bbox_inches="tight")

    if show_plot:
        plt.show()

    return loaded_results_dict

def point_is_in_rect(x: float, y: float, rect: list):
    """
    Check if a given x-y point is inside a rectangular field (e.g. inside a conductor).

    :param x: x coordinate of the point to check
    :type x: float
    :param y: y coordinate of the point to check
    :type y: float
    :param rect: rectangular
    :type rect: list
    """
    # x, y of the point
    # List of 4 points given as tuples with (x, y) in the order top-right, top-left, bottom-right, bottom-left

    # Return true if point is in rect
    if y < rect[0][1] and y > rect[3][1] and x > rect[0][0] and x < rect[1][0]:
        return True
    return False


def get_number_of_turns_of_winding(winding_windows: list, windings: list, winding_number: int):
    """
    Get the number of turns of a winding.

    :param winding_windows: list of winding windows
    :type winding_windows: list
    :param windings: list of windings
    :type windings: list
    :param winding_number: number of winding
    :type winding_number: int
    """
    turns = 0
    for ww in winding_windows:
        for vww in ww.virtual_winding_windows:
            for index, winding in enumerate(windings):
                if winding.winding_number == winding_number:
                    # TODO: change index_turns right no. of winding numbers, right position in list and
                    #  length of list is needed
                    try:
                        turns += vww.turns[index]
                    except:
                        pass
    return turns


def cost_function_core(core_weight: float, core_type: str = "ferrite") -> float:
    """
    Calculate core material costs depending on material and weight.

    :param core_weight: core weight in kg
    :type core_weight: float
    :param core_type: core type. Can be "ferrite", "amorphous", "nanocristalline", "high_si_steel", "goes"
    :type core_type: str
    :return: costs of core in euro
    :rtype: float
    """
    cost_database = cost_material_database()
    sigma_core = cost_database["core_euro_per_kilogram"][core_type]

    return sigma_core * core_weight


def cost_function_winding(wire_weight_list: list[float], wire_type_list: list[str],
                          single_strand_cross_section_list: list[float] | None = None):
    """
    Calculate single winding material and fabrication costs depending on winding-type and weight.

    Reference: Ralph Burkart and Johann W. Kolar: "Component Cost Models for Multi-Objective Optimizations
    of Switched-Mode Power Converters"

    :param wire_weight_list: winding weight in kg in list-form
    :type wire_weight_list: list[float]
    :param wire_type_list: winding type. Must fit to enum-names in ConductorType-Enum
    :type wire_type_list: list[str]
    :param single_strand_cross_section_list: single strand cross-section in list-form
    :type single_strand_cross_section_list: list[float]
    :return: winding cost of single winding
    :rtype: float
    """
    if single_strand_cross_section_list is None:
        single_strand_cross_section_list = []

    cost_database = cost_material_database()
    winding_cost_list = []

    for winding_count, winding_weight in enumerate(wire_weight_list):
        # material cost (per kilogram and per unit)
        sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram"][
            wire_type_list[winding_count]]
        if sigma_material_winding_euro_per_kilogram == -1:
            # case for special litz wire calculation. Additional data is loaded from cost_database.
            sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram_for_litz"][
                "sigma_numerator"] / (single_strand_cross_section_list[winding_count] * 1e6 + cost_database[
                    "winding_material_euro_per_kilogram_for_litz"]["sigma_denominator"])

        winding_material_euro_per_unit = cost_database["winding_material_euro_per_unit"][wire_type_list[winding_count]]

        winding_material_cost = sigma_material_winding_euro_per_kilogram * winding_weight + \
            winding_material_euro_per_unit

        # fabrication cost (per kilogram and per unit)
        sigma_fabrication_euro_per_kilogram = cost_database["winding_fabrication_euro_per_kilogram"][
            wire_type_list[winding_count]]
        fabrication_material_euro_per_unit = cost_database["winding_fabrication_euro_per_unit"][
            wire_type_list[winding_count]]

        winding_fabrication_cost = sigma_fabrication_euro_per_kilogram * winding_weight + fabrication_material_euro_per_unit

        winding_cost_list.append(winding_material_cost + winding_fabrication_cost)

    return winding_cost_list


def cost_function_total(core_weight: float, core_type: str, wire_weight_list: list[float], wire_type_list: list[str],
                        single_strand_cross_section_list: None | list[float] = None) -> float:
    """
    Calculate the total costs for an inductive element.

    This includes material costs for core and winding, fabrication costs for core and winding and manufacturer margin

    Reference: Ralph Burkart and Johann W. Kolar: "Component Cost Models for Multi-Objective Optimizations of
    Switched-Mode Power Converters"

    :param core_weight: core weight in kg
    :type core_weight: float
    :param core_type: core type. Can be "ferrite", "amorphous", "nanocristalline", "high_si_steel", "goes"
    :type core_type: str
    :param wire_weight_list: winding weight in kg
    :type wire_weight_list: float
    :param wire_type_list: winding type in list-form. Must fit to enum-names in ConductorType-Enum
    :type wire_type_list: list[str]
    :param single_strand_cross_section_list: single strand cross-section in list-form
    :type single_strand_cross_section_list: list[float]
    :return: total costs for inductive element
    :rtype: float
    """
    if single_strand_cross_section_list is None:
        single_strand_cross_section_list = []

    cost_database = cost_material_database()

    cost_core = cost_function_core(core_weight, core_type)

    cost_winding_list = cost_function_winding(wire_weight_list, wire_type_list, single_strand_cross_section_list)
    cost_winding = sum(cost_winding_list)

    total_cost_including_margin = 1 / (1 - cost_database["gross_margin"]) * (cost_core + cost_winding)

    return total_cost_including_margin


def find_result_log_file(result_log_folder: str, keyword_list: list, value_min_max: list):
    """
    Find a result log-file in a folder with many result-log files.

    Check a dictionary keyword list for matching a certain value (equal, greater equal, smaller equal).

    :param result_log_folder: filepath to result-log folder
    :type result_log_folder: str
    :param keyword_list: list with hierarchical keywords for dictionary structure, e.g. ["simulation_settings", "core", "core_inner_diameter"]
    :type keyword_list: list
    :param value_min_max: value to check for
    :type value_min_max: list

    :Example:

    Check for files with a core inner diameter smaller equal than 0.02 m.
    >>> import femmt as fmt
    >>> fmt.find_result_log_file("/home/filepath/fem_simulation_data", ["simulation_settings", "core",
    >>>     "core_inner_diameter"],[0.015, 0.02])
    """
    files_list = os.listdir(result_log_folder)

    value_min = value_min_max[0]
    value_max = value_min_max[1]

    for file in files_list:
        file_path = os.path.join(result_log_folder, file)
        with open(file_path, "r") as fd:
            full_data = json.loads(fd.read())

        if len(keyword_list) == 2:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]]
        elif len(keyword_list) == 3:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]]
        elif len(keyword_list) == 4:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]][keyword_list[3]]
        elif len(keyword_list) == 5:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]][keyword_list[3]][
                keyword_list[4]]

        if value_min <= data_to_compare <= value_max:
            logger.info(f"{value_min} <= {data_to_compare} <= {value_max} for file named {file}")


def wave_vector(f: float, complex_permeability: complex, complex_permittivity: complex, conductivity: float):
    """
    Calculate the wave-vector of a signal inside the core material with its material parameters.

    :param f: frequency
    :type f: float
    :param complex_permeability: complex permeability
    :type complex_permeability: complex
    :param complex_permittivity: complex permittivity
    :type complex_permittivity: complex
    :param conductivity: conductivity of the core material
    :type conductivity: float
    """
    omega = 2 * np.pi * f
    j = complex(0, 1)
    complex_equivalent_permittivity = complex_permittivity - j * conductivity / omega
    return omega * np.sqrt(complex_permeability * complex_equivalent_permittivity)


def axial_wavelength(f: float, complex_permeability: float, complex_permittivity: float, conductivity: float):
    """
    Calculate the axial wavelength for a given frequency.

    :param f: Frequency in Hz
    :type f: float
    :param complex_permeability: complex permeability
    :type complex_permeability: float
    :param complex_permittivity: complex permittivity
    :type complex_permittivity: float
    :param conductivity: electrical conductivity
    :type conductivity: float
    """
    k = wave_vector(f, complex_permeability, complex_permittivity, conductivity)
    return 2 * np.pi / k.real


def check_mqs_condition(radius: float, frequency: float, complex_permeability: float, complex_permittivity: float,
                        conductivity: float, relative_margin_to_first_resonance: float = 0.5):
    """
    Check if the condition for a magnetoquasistatic simulation is fulfilled.

    Calculates the ratio (core-diameter / wavelength) and includes a safety margin factor of 0.5.
    In case of ratio > 1, the simulated frequency is too high. A magnetoquasistatic simulation will not lead to good
    results. It is recommended to reduce the frequency or use a full-wave solver (not supported by FEMMT).

    :param radius: core radius
    :type radius: float
    :param frequency: frequency in Hz
    :type frequency: float
    :param complex_permeability: complex permeability
    :type complex_permeability: float
    :param complex_permittivity: complex permittivity
    :type complex_permittivity: float
    :param conductivity: core conductivity
    :type conductivity: float
    :param relative_margin_to_first_resonance: relative margin to the first resonance. Defaults to 0.5.
    :type relative_margin_to_first_resonance: float
    """
    if frequency == 0:
        raise ValueError("check_mqs_condition() only works for frequencies != 0")

    axial_lambda = axial_wavelength(frequency, complex_permeability, complex_permittivity, conductivity)
    diameter_to_wavelength_ratio_of_first_resonance = 0.7655
    diameter_to_wavelength_ratio = 2 * radius / axial_lambda
    if diameter_to_wavelength_ratio > diameter_to_wavelength_ratio_of_first_resonance * relative_margin_to_first_resonance:
        # raise Warning(f"Resonance Ratio: {diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance} - "
        #               f"1 means 1st resonance - should be kept well below 1 to ensure MQS approach to be correct! ")
        logger.info(f"Resonance Ratio: {diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance}")


def create_open_circuit_excitation_sweep(I0: float, n: float, frequency: float) -> list[list[float]]:
    """
    Create a circuit excitation sweep with the other windings unloaded.

    :param I0: current in A
    :type I0: float
    :param n: turns ratio n
    :type n: float
    :param frequency: Frequency in Hz
    :type frequency: float
    """
    frequencies = [frequency] * n
    currents = [[0] * n for _ in range(n)]
    phases = [[180] * n for _ in range(n)]

    for x in range(0, n):
        for y in range(0, n):
            if x == y:
                currents[x][y] = I0
                phases[x][y] = 0

    return frequencies, currents, phases


def list_to_complex(complex_list: list):
    """
    Brings a list of two numbers (where first is real part, second is imaginary part) into a python specific complex number.

    :param complex_list:
    :type complex_list: list
    :return: complex number
    :rtype: complex
    """
    return complex(complex_list[0], complex_list[1])


def get_self_inductances_from_log(log: dict) -> list:
    """
    Read the self-inductances from the result log file (dictionary).

    :param log: Result log dictionary
    :type log: dict
    :return: self-inductances in a list
    :rtype: list
    """
    self_inductances = []
    for ol_index, open_loop_result in enumerate(log["single_sweeps"]):
        active_winding_name = f"winding{ol_index + 1}"
        self_inductances.append(list_to_complex(open_loop_result[active_winding_name]["flux_over_current"]))
    return self_inductances


def get_flux_linkages_from_log(log: dict) -> list:
    """
    Read the flux-linkages from the result log file (dictionary).

    :param log: Result log dictionary
    :type log: dict
    :return: flux-linkages in a list
    :rtype: list
    """
    flux_linkages = []
    for ol_index, open_loop_result in enumerate(log["single_sweeps"]):
        flux_linkages.append([])
        for winding_index in range(0, len(log["single_sweeps"])):
            flux_linkages[ol_index].append(list_to_complex(open_loop_result[f"winding{winding_index + 1}"]["flux"]))
    return flux_linkages


def get_coupling_matrix(flux_linkages: list) -> np.array:
    """
    Calculate the coupling factors from the given flux linkages.

    :param flux_linkages: flux-linkages
    :type flux_linkages: list
    :return: coupling-matrix in a matrix (np.array)
    :rtype: np.array
    """
    coupling_matrix = [[None] * len(flux_linkages) for _ in range(len(flux_linkages))]
    for self_index in range(0, len(flux_linkages)):
        for cross_index in range(0, len(flux_linkages)):
            coupling_matrix[cross_index][self_index] = flux_linkages[cross_index][self_index].real / \
                flux_linkages[self_index][self_index].real
    return coupling_matrix


def get_mean_coupling_factors(coupling_matrix: np.array):
    """
    Calculate the mean coupling factors from the coupling matrix.

    :param coupling_matrix: matrix with coupling factors between windings
    :type coupling_matrix: np.array
    """
    mean_coupling_factors = [[None] * len(coupling_matrix) for _ in range(len(coupling_matrix))]
    for self_index in range(0, len(coupling_matrix)):
        for cross_index in range(0, len(coupling_matrix)):
            mean_coupling_factors[cross_index][self_index] = (coupling_matrix[cross_index][self_index] * \
                                                              coupling_matrix[self_index][cross_index]) ** 0.5
    return mean_coupling_factors


def get_inductance_matrix(self_inductances: np.array, mean_coupling_factors: np.array, coupling_matrix: np.array):
    """
    Get the inductance matrix from self_inductances, mean_coupling_factors and the coupling_matrix.

    :param self_inductances: matrix with self inductances in H
    :type self_inductances: np.array
    :param mean_coupling_factors: mean coupling factors
    :type mean_coupling_factors: np.array
    :param coupling_matrix: matrix with coupling factors
    :type coupling_matrix: np.array
    """
    inductance_matrix = [[None] * len(mean_coupling_factors) for _ in range(len(mean_coupling_factors))]
    for x in range(0, len(coupling_matrix)):
        for y in range(0, len(coupling_matrix)):
            inductance_matrix[x][y] = mean_coupling_factors[x][y] * (self_inductances[x] * self_inductances[y]) ** 0.5
    return inductance_matrix


def visualize_flux_linkages(flux_linkages: list) -> None:
    """
    Print the flux linkages to the terminal (or file-) output.

    :param flux_linkages: flux-linkages in a list
    :type flux_linkages: list
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        for y in range(0, len(flux_linkages)):
            string_to_print += f"Phi_{x+1}{y+1} = {flux_linkages[x][y]}     Induced by I_{y+1} in Winding{x+1}\n"

    logger.info("\nFluxes: ")
    logger.info(string_to_print)


def visualize_self_inductances(self_inductances: list | np.ndarray, flux_linkages: list | np.ndarray) -> None:
    """
    Print the self-inductances to the terminal (or file-) output.

    :param self_inductances: self-inductances in H in a list or numpy array
    :type self_inductances: list | np.ndarray
    :param flux_linkages: flux linkages
    :type flux_linkages: st | np.ndarray
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        string_to_print += f"L_{x+1}_{x+1} = {self_inductances[x]}\n"
    logger.info("\n"
                "Self Inductances: ")
    logger.info(string_to_print)


def visualize_self_resistances(self_inductances: list, flux_linkages: list, frequency: float) -> None:
    """
    Calculate and print the self resistances to the terminal (or file-) output.

    :param self_inductances: self-inductances in a list
    :type self_inductances: list
    :param flux_linkages: flux-linkage
    :type flux_linkages: list
    :param frequency: Frequency in Hz
    :type frequency: float
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        string_to_print += f"Z_{x+1}_{x+1} = {self_inductances[x].imag*2*np.pi*frequency}\n"
    logger.info("\n"
                "Self Resistances: ")
    logger.info(string_to_print)


def visualize_coupling_factors(coupling_matrix: np.array, flux_linkages: list):
    """
    Print the coupling factors to the terminal (or file-) output.

    :param coupling_matrix: matrix with coupling factors between the windings
    :type coupling_matrix: np.array
    :param flux_linkages: flux-linkages in a list
    :type flux_linkages: list
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        for y in range(0, len(coupling_matrix)):
            string_to_print += f"K_{x + 1}{y + 1} = Phi_{x + 1}{y + 1} / Phi_{y + 1}{y + 1} = {coupling_matrix[x][y]}\n"
    logger.info("\n"
                "Coupling Factors: ")

    logger.info(string_to_print)


def visualize_mean_coupling_factors(mean_coupling_factors: list):
    """
    Print the mean coupling factors to the terminal (or file-) output.

    :param mean_coupling_factors: mean_coupling_factors in a list
    :type mean_coupling_factors: list
    """
    string_to_print = ""
    for x in range(0, len(mean_coupling_factors)):
        for y in range(0, len(mean_coupling_factors)):
            string_to_print += f"k_{x + 1}{y + 1} = Sqrt(K_{x + 1}{y + 1} * K_{y + 1}{x + 1}) = M_{x + 1}{y + 1} " \
                               f"/ Sqrt(L_{x + 1}_{x + 1} * L_{y + 1}_{y + 1}) = {mean_coupling_factors[x][y]}\n"
    logger.info("\nMean Coupling Factors: ")
    logger.info(string_to_print)


def visualize_mean_mutual_inductances(inductance_matrix: np.array):
    """
    Print the mean mutual inductances to the terminal (or file-) output.

    :param inductance_matrix: inductance matrix
    :type inductance_matrix: np.array

    e.g.  M_12 = M_21 = k_12 * (L_11 * L_22) ** 0.5
    """
    string_to_print = ""
    for x in range(0, len(inductance_matrix)):
        for y in range(0, len(inductance_matrix)):
            if x == y:
                pass
            else:
                string_to_print += f"M_{x + 1}{y + 1} = {inductance_matrix[x][y].real}\n"
    logger.info("\nMean Mutual Inductances: ")
    logger.info(string_to_print)


def visualize_mutual_inductances(self_inductances: list, coupling_factors: list):
    """
    Print the mutual inductances to the terminal (or file-) output.

    :param self_inductances: Matrix with self inductances
    :type self_inductances: list
    :param coupling_factors: Matrix with coupling factors
    :type coupling_factors: list

    e.g. M_12 = L_11 * K_21  !=   M_21 = L_22 * K_12   (ideally, they are the same)
    """
    string_to_print = ""
    for x in range(0, len(coupling_factors)):
        for y in range(0, len(coupling_factors)):
            if x == y:
                pass
            else:
                string_to_print += f"M_{x + 1}{y + 1} = {self_inductances[y].real * coupling_factors[x][y]}\n"
    logger.info("\nMutual Inductances: ")
    logger.info(string_to_print)


def visualize_inductance_matrix_coefficients(inductance_matrix: np.array):
    """Visualize the inductance matrix coefficients in the terminal.

    e.g. M_12 = L_11 * K_21  !=   M_21 = L_22 * K_12   (ideally, they are the same)

    :param inductance_matrix: inductance matrix of transformer
    :type inductance_matrix: np.array
    """
    string_to_print = ""
    for x in range(0, len(inductance_matrix)):
        for y in range(0, len(inductance_matrix)):
            if x == y:
                string_to_print += f"L_{x + 1}{y + 1} = {inductance_matrix[x][y].real}\n"
            else:
                string_to_print += f"M_{x + 1}{y + 1} = {inductance_matrix[x][y].real}\n"
    logger.info("Inductance Matrix Coefficients: ")
    logger.info(string_to_print)


def visualize_inductance_matrix(inductance_matrix: np.array) -> None:
    """
    Visualize the inductance matrix in the terminal.

    :param inductance_matrix: inductance matrix in H
    :type inductance_matrix: np.array

    e.g. M_12 = L_11 * K_21  !=   M_21 = L_22 * K_12   (ideally, they are the same)
    """
    string_to_print = ""
    for x in range(0, len(inductance_matrix)):
        for y in range(0, len(inductance_matrix)):
            string_to_print += f"{np.round(inductance_matrix[x][y].real, 12)} "
        string_to_print += "\n"

    logger.info("\nInductance Matrix: ")
    logger.info(string_to_print)

def calculate_quadrature_integral(time_steps: list[float], data: list[float]) -> float:
    """
    Calculate the integral of given data over specific time steps using the quad method.

    :param time_steps: List of time steps.
    :type time_steps: list[float]
    :param data: List of data corresponding to each timestep.
    :type data: list[float]
    :return: The calculated integral.
    :rtype: float
    """
    func = lambda x: np.interp(x, time_steps, data)
    return quad(func, time_steps[0], time_steps[-1])[0]

def calculate_squared_quadrature_integral(time_steps: list[float], data: list[float]) -> float:
    """
    Calculate the integral of squared given data over specific time steps using the quad method.

    :param time_steps: List of time steps.
    :type time_steps: list[float]
    :param data: List of data corresponding to each timestep.
    :type data: list[float]
    :return: The calculated integral.
    :rtype: float
    """
    func = lambda x: np.interp(x, time_steps, data) ** 2
    return quad(func, time_steps[0], time_steps[-1])[0]

def calculate_average(integral: float, time_steps: list[float]) -> float:
    """
    Compute the average in general.

    :param integral: The integral value.
    :type integral: float
    :param time_steps: List of time steps.
    :type time_steps: list[float]

    Returns:
    :return: The calculated average.
    :rtype: float.
    """
    total_time = time_steps[-1] - time_steps[0]
    if total_time == 0:
        raise ValueError("Total time cannot be zero.")
    return integral / total_time

def calculate_rms(squared_integral: float, time_steps: list[float]) -> float:
    """
    Compute the RMS.

    :param squared_integral: The integral value.
    :type squared_integral: float
    :param time_steps: List of time steps.
    :type time_steps: list[float]

    Returns:
    :return: The calculated average.
    :rtype: float.
    """
    total_time = time_steps[-1] - time_steps[0]
    if total_time == 0:
        raise ValueError("Total time cannot be zero.")
    mean_square = squared_integral / total_time  # Calculate the mean of the square of the data
    return np.sqrt(mean_square)  # Take the square root to get RMS value


def convert_air_gap_corner_points_to_center_and_distance(corner_points: list) -> list:
    """
    Convert the list-defined air_gap_corner_points from a "two_d_axi" object to center points and lengths as to separate lists.

    :param corner_points: in usage of magnetic component -> "self.two_d_axi.p_air_gaps.tolist()"
    :type corner_points: list
    :returns: centers and heights of the air gaps
    :rtype: list
    """
    centers = []
    heights = []
    # 4 corner points make up one air gap
    for n_air_gap in range(0, int(len(corner_points)/4)):
        width = corner_points[n_air_gap*4 + 1][0] - corner_points[n_air_gap*4 + 0][0]
        height = corner_points[n_air_gap * 4 + 2][1] - corner_points[n_air_gap * 4 + 0][1]
        heights.append(height)
        centers.append(
            [
                corner_points[n_air_gap*4 + 0][0] + width/2,  # x-coordinate of air gap's center
                corner_points[n_air_gap*4 + 2][1] - height/2  # y-coordinate of air gap's center
            ]
        )
    return centers, heights


def time_current_vector_to_fft_excitation(time_current_vectors: list[list[list[float]]], fft_filter_value_factor: float = 0.01):
    """
    Perform FFT to get the primary and secondary currents e.g. to calculate the wire losses.

    For further calculations e.g. calculating wire losses, the single frequencies can be 'linear added' to get the total winding losses.

    :param time_current_vectors: primary and secondary current waveforms over time
    :type time_current_vectors: list[list[list[float]]]
    :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default]
        removes all amplitudes below 1 % of the maximum amplitude from the result-frequency list
    :type fft_filter_value_factor: float
    """
    # winding losses
    frequency_current_phase_deg_list = []
    # collect winding losses simulation input parameters
    for time_current_vector in time_current_vectors:
        [frequency_list, amplitude, phi_rad] = fft(time_current_vector, mode='time', filter_value_factor=fft_filter_value_factor)
        phi_deg = np.rad2deg(phi_rad)
        frequency_current_phase_deg_list.append([frequency_list, amplitude, phi_deg])

    # check if all frequency vectors include the same frequencies
    # WORKAROUND: if any frequency is not included in one of the vectors it is
    # added with amplitude  = 0 and phase = 0
    # TODO: recalculate the fft at the "missing frequencies and add their values...
    all_frequencies = set()
    for count in range(len(frequency_current_phase_deg_list) - 1):
        if not np.array_equal(frequency_current_phase_deg_list[count][0], frequency_current_phase_deg_list[count + 1][0]):
            all_frequencies = all_frequencies | set(frequency_current_phase_deg_list[count][0]) | set(
                frequency_current_phase_deg_list[count + 1][0])

    for frequency in list(all_frequencies):
        for count in range(0, len(frequency_current_phase_deg_list)):
            if frequency not in frequency_current_phase_deg_list[count][0]:
                ii = np.searchsorted(frequency_current_phase_deg_list[count][0], frequency)
                frequency_current_phase_deg_list[count][0] = np.insert(
                    frequency_current_phase_deg_list[count][0], ii, frequency)
                frequency_current_phase_deg_list[count][1] = np.insert(
                    frequency_current_phase_deg_list[count][1], ii, 0)
                frequency_current_phase_deg_list[count][2] = np.insert(
                    frequency_current_phase_deg_list[count][2], ii, 0)

    # transfer format from fft()-output to excitation_sweep()-input
    current_list_list = []
    phi_deg_list_list = []
    for count_frequency, _ in enumerate(frequency_list):
        currents_single_frequency = []
        phi_deg_single_frequency = []
        for count_current, _ in enumerate(time_current_vectors):
            currents_single_frequency.append(frequency_current_phase_deg_list[count_current][1][count_frequency])
            phi_deg_single_frequency.append(frequency_current_phase_deg_list[count_current][2][count_frequency])
        current_list_list.append(currents_single_frequency)
        phi_deg_list_list.append(phi_deg_single_frequency)
    return frequency_list, current_list_list, phi_deg_list_list


def hysteresis_current_excitation(input_time_current_vectors: list[list[list[float]]]):
    """
    Collect the peak current and the corresponding phase shift for the fundamental frequency for all windings.

    Results are used for calculating the hysteresis losses by another function.
    In case of a center-tapped transformer, halving the amplitudes will be done by split_hysteresis_loss_excitation_center_tapped.

    :param input_time_current_vectors: e.g. [[time_vec, i_primary_vec], [time_vec, i_secondary_vec]]
    :type input_time_current_vectors: list[list[list[float]]]
    :raises ValueError: if time vector does not start at zero seconds.
    :return: hyst_frequency, hyst_current_amplitudes, hyst_phases_deg, e.g. 200400.80170764355 [6.13, 26.65] [49.13, 229.49]
    :rtype: list[list[float]]
    """
    if input_time_current_vectors[0][0][0] != 0:
        raise ValueError("time must start at 0 seconds!")

    # collect simulation input parameters from time_current_vectors
    hyst_current_amplitudes = []
    hyst_phases_deg = []
    hyst_frequency = 1 / (input_time_current_vectors[0][0][-1])
    for time_current_vector in input_time_current_vectors:
        # collect hysteresis loss simulation input parameters
        hyst_current_amplitudes.append(fr.max_value_from_value_vec(time_current_vector[1])[0])
        hyst_phases_deg.append(
            fr.phases_deg_from_time_current(time_current_vector[0], time_current_vector[1])[0])
    return hyst_frequency, hyst_current_amplitudes, hyst_phases_deg

def close_excel_file_if_open(filepath):
    """
    Close the specified Excel file if it is currently open.

    :param filepath: The path to the Excel file to close.
    :type filepath: str
    """
    # Get the absolute path
    filepath = os.path.abspath(filepath)

    try:
        excel = win32com.client.Dispatch("Excel.Application")
        for workbook in excel.Workbooks:
            if workbook.FullName.lower() == filepath.lower():
                workbook.Close(SaveChanges=False)
                return
        excel.Quit()
    except Exception as e:
        print(f"Unable to close Excel. Error: {e}")

def json_to_excel(json_file_path: str, output_excel_path: str) -> None:
    """
    Extract data from the electrostatic simulation and write it into a log file.

    :param json_file_path: Path to the JSON input file containing simulation results.
    :type json_file_path: str
    :param output_excel_path: Path where the Excel (.xlsx) file will be saved.
    :type output_excel_path: str
    :rtype: None
    """
    # Trying to close the Excel file if it's open
    close_excel_file_if_open(output_excel_path)
    # Load the JSON data from the file
    with open(json_file_path, 'r') as json_file:
        data = json.load(json_file)

    # Prepare the different data sections
    charges_data = []
    energy_data = []
    average_voltages_data = []
    capacitances_within_data = []
    capacitances_between_data = []
    capacitances_between_turns_core_data = []

    # Extract charges
    charge_value = data.get("charges", None)
    if charge_value is not None:
        charges_data.append({"Charge Type": "Total Charge", "Value (Coulombs)": charge_value})

    # Extract energy
    for key, value in data.get("energy", {}).items():
        energy_data.append({"Energy Type": key, "Value (Joules)": value})

    # Extract average voltages
    for region, voltage in data.get("average_voltages", {}).items():
        average_voltages_data.append({"Region": region, "Average Voltage (V)": voltage})

    # Extract capacitance within windings
    for winding, turns in data.get("capacitances", {}).get("within_winding", {}).items():
        for turn, connections in turns.items():
            for target_turn, capacitance_value in connections.items():
                capacitances_within_data.append({
                    "Winding": winding,
                    "Turn": turn,
                    "To Turn": target_turn,
                    "Capacitance (F)": capacitance_value
                })

    # Extract capacitance between windings
    for winding1, windings in data.get("capacitances", {}).get("between_windings", {}).items():
        for winding2, turns in windings.items():
            for turn1, connections in turns.items():
                for turn2, capacitance_value in connections.items():
                    capacitances_between_data.append({
                        "Winding 1": winding1,
                        "Turn 1": turn1,
                        "Winding 2": winding2,
                        "Turn 2": turn2,
                        "Capacitance (F)": capacitance_value
                    })
    # Extract capacitance between turns and core
    for winding, turns in data.get("capacitances", {}).get("between_turns_core", {}).items():
        for turn, capacitance_value in turns.items():
            capacitances_between_turns_core_data.append({
                "Winding": winding,
                "Turn": turn,
                "Capacitance to Core (F)": capacitance_value
            })

    # Create DataFrames for each section
    charges_df = pd.DataFrame(charges_data)
    energy_df = pd.DataFrame(energy_data)
    average_voltages_df = pd.DataFrame(average_voltages_data)
    capacitances_within_df = pd.DataFrame(capacitances_within_data)
    capacitances_between_df = pd.DataFrame(capacitances_between_data)
    capacitances_between_turns_core_df = pd.DataFrame(capacitances_between_turns_core_data)

    # Write to Excel file with multiple sheets
    with pd.ExcelWriter(output_excel_path) as writer:
        if not charges_df.empty:
            charges_df.to_excel(writer, sheet_name='Charges', index=False)
            worksheet = writer.sheets['Charges']
            worksheet.set_column('A:B', 30)
        if not energy_df.empty:
            energy_df.to_excel(writer, sheet_name='Energy', index=False)
            worksheet = writer.sheets['Energy']
            worksheet.set_column('A:B', 30)
        if not average_voltages_df.empty:
            average_voltages_df.to_excel(writer, sheet_name='Average_Voltages', index=False)
            worksheet = writer.sheets['Average_Voltages']
            worksheet.set_column('A:E', 30)
        if not capacitances_within_df.empty:
            capacitances_within_df.to_excel(writer, sheet_name='Capacitances_Within', index=False)
            worksheet = writer.sheets['Capacitances_Within']
            worksheet.set_column('A:D', 30)
        if not capacitances_between_df.empty:
            capacitances_between_df.to_excel(writer, sheet_name='Capacitances_Between', index=False)
            worksheet = writer.sheets['Capacitances_Between']
            worksheet.set_column('A:E', 30)
        if not capacitances_between_turns_core_df.empty:
            capacitances_between_turns_core_df.to_excel(writer, sheet_name='Turns_Core', index=False)
            worksheet = writer.sheets['Turns_Core']
            worksheet.set_column('A:C', 30)

def compare_excel_files(femmt_excel_path: str, femm_excel_path: str, comparison_output_path: str) -> None:
    """
    Compare two Excel files (FEMMT and FEMM) and generate a new Excel file with comparison results.

    This function loads two Excel files, one generated by FEMMT and the other by FEMM, compares the data across
    all common sheets, and calculates the differences between corresponding values. The results include:
    - Absolute Difference
    - Relative Error
    - Relative Error Percentage

    The comparison is saved into a new Excel file, with each comparison in a separate sheet named after
    the original sheet with the "_Comparison" suffix.

    :param femmt_excel_path: Path to the Excel file generated by FEMMT.
    :type femmt_excel_path: str
    :param femm_excel_path: Path to the Excel file generated by FEMM.
    :type femm_excel_path: str
    :param comparison_output_path: Path to save the resulting comparison Excel file.
    :type comparison_output_path: str
    """
    # Trying to close the Excel file if it's open
    close_excel_file_if_open(comparison_output_path)
    # Load both Excel files, get all sheets
    femmt_sheets = pd.read_excel(femmt_excel_path, sheet_name=None)
    femm_sheets = pd.read_excel(femm_excel_path, sheet_name=None)

    # Define sheets to compare
    sheets_to_compare = femmt_sheets.keys()

    # Create an Excel writer for the output
    with pd.ExcelWriter(comparison_output_path, engine='xlsxwriter') as writer:
        # Iterate through each sheet to compare
        for sheet_name in sheets_to_compare:
            if sheet_name in femm_sheets:
                # Load DataFrames for the current sheet
                femmt_df = femmt_sheets[sheet_name]
                femm_df = femm_sheets[sheet_name]
                # Rename columns (FEMMT and FEMM)
                femmt_df.columns = [f"{col}_FEMMT" for col in femmt_df.columns]
                femm_df.columns = [f"{col}_FEMM" for col in femm_df.columns]

                # Concatenate both DataFrames side by side
                comparison_df = pd.concat([femmt_df, femm_df], axis=1)

                # Calculating difference, relative error, and relative error in percentage for columns
                for femmt_col, femm_col in zip(femmt_df.columns, femm_df.columns):
                    col_name = femmt_col.replace("_FEMMT", "")
                    if np.issubdtype(comparison_df[femmt_col].dtype, np.number) and np.issubdtype(comparison_df[femm_col].dtype, np.number):
                        comparison_df[f"{col_name}_Difference"] = comparison_df[femmt_col] - comparison_df[femm_col]
                        comparison_df[f"{col_name}_Relative_Error"] = comparison_df[f"{col_name}_Difference"] / comparison_df[femmt_col].replace(0, np.nan)
                        comparison_df[f"{col_name}_Error_Percent"] = comparison_df[f"{col_name}_Relative_Error"] * 100

                # Writing to the Excel output file
                comparison_df.to_excel(writer, sheet_name=f"{sheet_name}_Comparison", index=False)
                worksheet = writer.sheets[f"{sheet_name}_Comparison"]
                worksheet.set_column('A:Z', 35)

if __name__ == '__main__':
    pass
