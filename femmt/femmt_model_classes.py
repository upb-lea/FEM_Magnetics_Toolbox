from dataclasses import dataclass
from typing import List
from femmt_enumerations import *
from femmt_functions import NbrLayers, wire_material_database
import numpy as np

class VirtualWindingWindow:
    """
    A virtual winding window is the area, where either some kind of interleaved conductors or a one winding
    (primary, secondary,...) is placed in a certain way.

    TODO More documentation
    """

    # Rectangular frame:
    bot_bound: float
    top_bound: float
    left_bound: float
    right_bound: float

    # Arrangement of the Conductors in the virtual winding window
    # Obviously depends on the chosen conductor type
    winding: List[WindingType]  # "interleaved" | "primary"/"secondary"
    scheme: List[WindingScheme]  # "bifilar", "vertical", "horizontal", ["hexa", "square"] | "hexa", "square"

    def __init__(self, winding: WindingType, scheme: WindingScheme):
        self.winding = winding
        self.scheme = scheme

class Winding:
    """
    A winding defines a conductor which is wound around a magnetic component such as transformer or inductance.
    The winding is defined by its conductor and the way it is placed in the magnetic component. To allow different
    arrangements of the conductors in several winding windows (hexagonal or square packing, interleaved, ...) in
    this class only the conductor parameters are specified. Then, by calling class:Winding in
    class:VirtualWindingWindow the arrangement of the conductors is specified.

    TODO More documentation
    """

    turns: List[int] #0: primary, 1: secondary
    conductor_type: ConductorType = None
    ff: float = None
    strand_radius: float = None
    n_strands: int = None
    n_layers: int
    conductor_radius: float = None
    a_cell: float
    thickness: float = None
    wrap_para: WrapParaType = WrapParaType.default
    cond_sigma: float
    parallel: int = None # TODO What is this parameter?
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

        self.turns_primary = turns_primary
        self.turns_secondary = turns_secondary

        dict_material_database = wire_material_database()
        if conductivity.value in dict_material_database:
            self.cond_sigma = dict_material_database[conductivity]["sigma"]
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
    frequency = 0: mu_rel only used if non_linear == False
    frequency > 0: mu_rel is used
    TODO Documentation
    """

    # Standard material data
    material: str  # "95_100" := TDK-N95 | Currently only works with Numbers corresponding to BH.pro

    # Permeability
    # TDK N95 as standard material:
    permeability_type: str
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
    
    steinmetz_loss: float = 0

    def __init__(self, core_w: float, window_w: float, window_h: float, material: str = "custom",  # "95_100" 
                   loss_approach: str = None, loss_data_source: str = "custom", mu_rel: float = 3000,
                   phi_mu_deg: float = None, sigma: float = None, non_linear: bool = False, **kwargs):
        # TODO This still needs to be reworked to the new way of handling input parameters
        self.core_w = core_w
        self.core_h = None # TODO Set core_h to not none
        self.window_w = window_w
        self.window_h = window_h

        print(f"Update the magnetic Core to {self.type}-type with following parameters: {kwargs}\n"
                f"---")

        # Material Properties
        self.core.material = material
        self.core.non_linear = non_linear

        # Conductivity
        if self.core.material == "custom":  # user defines the conductivity
            self.core.sigma = sigma
        if loss_approach == "Steinmetz":  # conductivity must be set to 0 for Steinmetz approach
            self.core.sigma = 0
        else:
            self.core.sigma = f"sigma_from_{self.core.material}"

        # Permeability
        self.mu_rel = mu_rel
        self.phi_mu_deg = phi_mu_deg

        # Check for which kind of permeability definition is used
        if loss_approach == "loss_angle":
            if self.core.phi_mu_deg is not None and self.core.phi_mu_deg != 0:
                self.core.permeability_type = "fixed_loss_angle"
            else:
                self.core.permeability_type = "real_value"
        else:
            if self.core.material != "custom":
                self.core.permeability_type = "from_data"
            # else:
            #     raise Exception(f"Permeability must be specified with real_value, fixed_loss_angle or from_data")


        # Set attributes of core with given keywords
        for key, value in kwargs.items():
            setattr(self, key, value)

@dataclass
class StrayPath:
    """
    TODO: Thickness of the stray path must be fitted for the real Tablet (effective area of the
    "stray air gap" is different in axi-symmetric approximation
    """

    start_index: int        # lower air gap that characterizes the stray path
    radius: float
    width: float
    midpoint: List[float]   # TODO correct Datatype?

class AirGaps:
    """
    Contains methods and arguments to describe the air gaps in a magnetic component

    TODO Add documentation
    """

    core: Core
    midpoints: List[float]  #: list: [position_tag, air_gap_position, air_gap_h, c_air_gap]

    def __init__(self, method: AirGapMethod, core: Core):
        if method != AirGapMethod.Center:
            raise Exception(f"The method {method} is currently not supported")
        self.method = method
        self.core = core

    def add_air_gap(self, leg_position: AirGapLegPosition, position_value: float, height: float):
        for index, midpoint in enumerate(self.midpoints):
            if midpoint[0] == leg_position and midpoint[1] + midpoint[2] < position_value - height \
                    and midpoint[1] - midpoint[2] > position_value + height:
                raise Exception(f"Air gaps {index} and {len(self.midpoints)} are overlapping")

        if self.method == AirGapMethod.Center:
            if len(self.number) >= 1:
                raise Exception("The 'center' position for air gaps can only have 1 air gap maximum")
            else:
                self.midpoints.append([0, 0, height])
                
        if self.method == AirGapMethod.Manually:
            self.midpoints.append([leg_position.value, position_value, height])
        if self.method == AirGapMethod.Percent:
            position = position_value / 100 * self.core.window_h - self.core.window_h / 2
            self.midpoints.append([leg_position.value, position, height])

class Isolation:
    """
    Isolation
        - Between two turns of common conductors: first n_conductor arguments of cond_cond
        - Between two neighboured conductors: last n_conductor-1 arguments

    :param cond_cond: list of floats to describe the isolations between conductors
    :type cond_cond: List
    :param core_cond: list of floats to describe the isolations between conductors and the core
    :type core_cond: List
    :return: None
    :rtype: None

    :Inductor Example:

    core_cond_isolation=[windings2top_core, windings2bot_core, windings2left_core, windings2right_core],
    cond_cond_isolation=[winding2primary]

    :Transformer Example:

    core_cond_isolation=[windings2top_core, windings2bot_core, windings2left_core, windings2right_core],
    cond_cond_isolation=[primary2primary, secondary2secondary, primary2secondary]
    """

    cond_cond: List[float] = []
    core_cond: List[float] = []

    def add_winding_isolations(self, primary2primary, secondary2secondary, primary2secondary):
        self.cond_cond = [primary2primary, secondary2secondary, primary2secondary]

    def add_core_isolations(self, top_core, bot_core, left_core, right_core):
        self.core_cond = [top_core, bot_core, left_core, right_core]
