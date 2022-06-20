from dataclasses import dataclass
from typing import List
from .femmt_enumerations import *
from .femmt_functions import NbrLayers, wire_material_database
import numpy as np

class Winding:
    """
    A winding defines a conductor which is wound around a magnetic component such as transformer or inductance.
    The winding is defined by its conductor and the way it is placed in the magnetic component. To allow different
    arrangements of the conductors in several winding windows (hexagonal or square packing, interleaved, ...) in
    this class only the conductor parameters are specified. 

    TODO More documentation
    """

    turns: List[int] #0: primary, 1: secondary
    conductor_type: ConductorType = None
    ff: float = None
    strand_radius: float = None
    n_strands: int = 0
    n_layers: int
    conductor_radius: float = None
    a_cell: float
    thickness: float = None
    wrap_para: WrapParaType = None
    cond_sigma: float
    parallel: int = 1 # TODO What is this parameter?
    winding_type: WindingType
    winding_scheme: WindingScheme

    def __init__(self, turns_primary: int, turns_secondary: int, conductivity: Conductivity, winding_type: WindingType, winding_scheme: WindingScheme):
        if turns_primary < 1 and turns_secondary < 1:
            raise Exception("Either number of primary or number of secondary turns need to be at least 1.")

        self.winding_type = winding_type

        if winding_scheme in [WindingScheme.Hexagonal, WindingScheme.Square, WindingScheme.Square_Full_Width]:
            # winding type needs to be primary or secondary
            if self.winding_type != WindingType.Primary and self.winding_type != WindingType.Secondary:
                raise Exception(f"For {winding_scheme} winding type needs to be primary or secondary")
        elif winding_scheme in [WindingScheme.Horizontal, WindingScheme.Vertical, WindingScheme.Bifilar, WindingScheme.Blockwise]:
            # winding type needs to be interleaved
            if self.winding_type != WindingType.Interleaved:
                raise Exception(f"For {winding_scheme} winding type needs to be interleaved")
        else:
            raise Exception(f"{winding_scheme} does not fit into any possible winding schemes")
        
        self.winding_scheme = winding_scheme

        self.turns = [turns_primary, turns_secondary]

        dict_material_database = wire_material_database()
        if conductivity.value in dict_material_database:
            self.cond_sigma = dict_material_database[conductivity.value]["sigma"]
        else:
            raise Exception(f"Material {conductivity.value} not found in database")

    def set_stacked_conductor(self, thickness: float, wrap_para: WrapParaType):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Stacked
        self.thickness = thickness
        self.wrap_para = wrap_para
        self.a_cell = 1 # TODO Surface size needed?
        self.conductor_radius = 1 # Revisit

    def set_full_conductor(self, thickness: float, wrap_para: WrapParaType):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Full
        self.thickness = thickness
        self.wrap_para = wrap_para
        self.a_cell = 1 # TODO Surface size needed?
        self.conductor_radius = 1 # Revisit

    def set_foil_conductor(self, thickness: float, wrap_para: WrapParaType):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Foil
        self.thickness = thickness
        self.wrap_para = wrap_para
        self.a_cell = 1 # TODO Surface size needed?
        self.conductor_radius = 1 # Revisit

    def set_solid_conductor(self, conductor_radius: float):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Solid
        self.conductor_radius = conductor_radius
        self.a_cell = np.pi * conductor_radius ** 2

    def set_litz_conductor(self, conductor_radius: float, number_strands: int, strand_radius: float, fill_factor: float):
        """
        Only 3 of the 4 parameters are needed. The other one needs to be none
        """
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Litz
        self.conductor_radius = conductor_radius
        self.n_strands = number_strands
        self.strand_radius = strand_radius
        self.ff = fill_factor

        if number_strands is None:
            self.n_strands = conductor_radius ** 2 / strand_radius ** 2 * fill_factor
        elif conductor_radius is None:
            self.conductor_radius = np.sqrt(number_strands * strand_radius ** 2 / fill_factor)
        elif fill_factor is None:
            ff_exact = number_strands * strand_radius ** 2 / conductor_radius ** 2
            self.ff = np.around(ff_exact, decimals=2)
        else:
            raise Exception("Wrong inputs for litz conductor")

        self.n_layers = NbrLayers(number_strands)
        self.a_cell = number_strands * strand_radius ** 2 * np.pi / fill_factor

        print(f"Updated Litz Configuration: \n"
              f" ff: {fill_factor} \n"
              f" Number of layers/strands: {self.n_layers}/{number_strands} \n"
              f" Strand radius: {strand_radius} \n"
              f" Conductor radius: {conductor_radius}\n"
              f"---")

class Core:
    """
    This creates the core base for the model.

    frequency = 0: mu_rel only used if non_linear == False
    frequency > 0: mu_rel is used
    
    TODO More Documentation
    """
    type: str

    # Standard material data
    material: str  # "95_100" := TDK-N95 | Currently only works with Numbers corresponding to BH.pro

    # Permeability
    # TDK N95 as standard material:
    permeability_type: PermeabilityType
    mu_rel: float           # Relative Permeability [if complex: mu_complex = re_mu_rel + j*im_mu_rel with mu_rel=|mu_complex|]
    phi_mu_deg: float       # mu_complex = mu_rel * exp(j*phi_mu_deg)
    # re_mu_rel: float      # Real part of relative Core Permeability  [B-Field and frequency-dependent]
    # im_mu_rel: float      # Imaginary part of relative Core Permeability

    # Permitivity - [Conductivity in a magneto-quasistatic sense]
    sigma: float            # Imaginary part of complex equivalent permittivity [frequency-dependent]

    # Dimensions
    core_w: float           # Axi symmetric case | core_w := core radius
    core_h: float
    window_w: float         # Winding window width
    window_h: float         # Winding window height
    core_type: str = "EI"   # Basic shape of magnetic conductor
    
    steinmetz_loss: int = 0
    generalized_steinmetz_loss: int = 0

    def __init__(self, core_w: float, window_w: float, window_h: float, material: str = "custom",  # "95_100" 
                   loss_approach: LossApproach = LossApproach.LossAngle, mu_rel: float = 3000,
                   phi_mu_deg: float = None, sigma: float = None, non_linear: bool = False, **kwargs):
        # Set parameters
        self.core_w = core_w
        self.core_h = None # TODO Set core_h to not none
        self.window_w = window_w
        self.window_h = window_h
        self.type = "axi_symmetric"
        self.material = material
        self.non_linear = non_linear
        self.mu_rel = mu_rel
        self.phi_mu_deg = phi_mu_deg

        # Check loss approach
        if loss_approach == LossApproach.Steinmetz:
            self.sigma = 0
            if self.material != "custom":
                self.permeability_type = PermeabilityType.FromData
                self.sigma = f"sigma_from_{self.material}"
            else:
                raise Exception(f"When steinmetz losses are set a material needs to be set as well.")
        elif loss_approach == LossApproach.LossAngle:
            if self.material == "custom":
                self.sigma = sigma
            else:
                self.sigma = f"sigma_from_{self.material}"

            if phi_mu_deg is not None and phi_mu_deg != 0:
                self.permeability_type = PermeabilityType.FixedLossAngle
            else:
                self.permeability_type = PermeabilityType.RealValue
        else:
            raise Exception("Loss approach {loss_approach.value} is not implemented")

        # Set attributes of core with given keywords
        # TODO Should we allow this? Technically this is not how an user interface should be designed
        for key, value in kwargs.items():
            setattr(self, key, value)

class AirGaps:
    """
    Contains methods and arguments to describe the air gaps in a magnetic component

    An air gap can be added with the add_air_gap function. It is possible to set different positions and heights.
    """

    core: Core
    midpoints: List[List[float]]  #: list: [position_tag, air_gap_position, air_gap_h]
    number: int

    def __init__(self, method: AirGapMethod, core: Core):
        self.method = method
        self.core = core
        self.midpoints = []
        self.number = 0

    def add_air_gap(self, leg_position: AirGapLegPosition, position_value: float, height: float):
        for index, midpoint in enumerate(self.midpoints):
            if midpoint[0] == leg_position and midpoint[1] + midpoint[2] < position_value - height \
                    and midpoint[1] - midpoint[2] > position_value + height:
                raise Exception(f"Air gaps {index} and {len(self.midpoints)} are overlapping")

        if self.method == AirGapMethod.Center:
            if self.number >= 1:
                raise Exception("The 'center' position for air gaps can only have 1 air gap maximum")
            else:
                self.midpoints.append([0, 0, height])
                self.number += 1
                
        elif self.method == AirGapMethod.Manually:
            self.midpoints.append([leg_position.value, position_value, height])
            self.number += 1
        elif self.method == AirGapMethod.Percent:
            position = position_value / 100 * self.core.window_h - self.core.window_h / 2
            self.midpoints.append([leg_position.value, position, height])
            self.number += 1
        else:
            raise Exception(f"Method {self.method} is not supported.")


class Isolation:
    """
    This class defines isolations for the model.
    An isolation between the winding window and the core can always be set.
    When having a inductor only the primary2primary isolation is necessary.
    When having a (integrated) transformer secondary2secondary and primary2secondary isolations can be set as well.
    """

    cond_cond: List[float] = []
    core_cond: List[float] = []

    def add_winding_isolations(self, primary2primary, secondary2secondary = None, primary2secondary = None):
        self.cond_cond = [primary2primary, secondary2secondary, primary2secondary]

    def add_core_isolations(self, top_core, bot_core, left_core, right_core):
        self.core_cond = [top_core, bot_core, left_core, right_core]

@dataclass
class StrayPath:
    """
    This class is needed when an integrated transformer shall be created.

    TODO: Thickness of the stray path must be fitted for the real Tablet (effective area of the
    "stray air gap" is different in axi-symmetric approximation
    """

    start_index: int        # lower air gap that characterizes the stray path
    radius: float
    width: float
    midpoint: List[List[float]]

class VirtualWindingWindow:
    """
    A virtual winding window is the area, where either some kind of interleaved conductors or a one winding
    (primary, secondary,...) is placed in a certain way.

    An instance of this class will be automatically created when the Winding is added to the MagneticComponent
    """

    # Rectangular frame:
    bot_bound: float
    top_bound: float
    left_bound: float
    right_bound: float

    # Arrangement of the Conductors in the virtual winding window
    # Obviously depends on the chosen conductor type
    winding: List[WindingType]
    scheme: List[WindingScheme]

    def __init__(self, winding: WindingType, scheme: WindingScheme):
        self.winding = winding
        self.scheme = scheme