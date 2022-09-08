# Python standard libraries
import numpy as np
from dataclasses import dataclass
from typing import List, Optional

# Local libraries
import Functions as ff
from Enumerations import *
from Data import MeshData

class Conductor:
    """
    A winding defines a conductor which is wound around a magnetic component such as transformer or inductance.
    The winding is defined by its conductor and the way it is placed in the magnetic component. To allow different
    arrangements of the conductors in several winding windows (hexagonal or square packing, interleaved, ...) in
    this class only the conductor parameters are specified. 

    TODO More documentation
    """

    conductor_type: ConductorType
    conductor_arrangement: ConductorArrangement
    wrap_para: WrapParaType = None
    conductor_radius: float = None
    winding_number: int
    thickness: float = None
    ff: float = None
    strand_radius: float = None
    n_strands: int = 0
    n_layers: int
    a_cell: float
    cond_sigma: float
    parallel: int = 1 # TODO What is this parameter?

    conductor_is_set: bool

    # Not used in femmt_classes. Only needed for to_dict()
    conductivity: Conductivity = None

    def __init__(self, winding_number: int, conductivity: float):
        if winding_number < 0:
            raise Exception("Winding index cannot be negative.")

        self.winding_number = winding_number
        self.conductivity = conductivity
        self.conductor_is_set = False

        dict_material_database = ff.wire_material_database()
        if conductivity.name in dict_material_database:
            self.cond_sigma = dict_material_database[conductivity.name]["sigma"]
        else:
            raise Exception(f"Material {conductivity.name} not found in database")

    def set_rectangular_conductor(self, thickness: float):
        if self.conductor_is_set:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_is_set = True
        self.conductor_type = ConductorType.RectangularSolid
        self.thickness = thickness
        self.a_cell = 1 # TODO Surface size needed?
        self.conductor_radius = 1 # Revisit

    def set_solid_round_conductor(self, conductor_radius: float, conductor_arrangement: ConductorArrangement):
        if self.conductor_is_set:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_is_set = True
        self.conductor_type = ConductorType.RoundSolid
        self.conductor_arrangement = conductor_arrangement
        self.conductor_radius = conductor_radius
        self.a_cell = np.pi * conductor_radius ** 2

    def set_litz_round_conductor(self, conductor_radius: float, number_strands: int, strand_radius: float, fill_factor: float, conductor_arrangement: ConductorArrangement):
        """
        Only 3 of the 4 parameters are needed. The other one needs to be none
        """
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_is_set = True
        self.conductor_type = ConductorType.RoundLitz
        self.conductor_arrangement = conductor_arrangement
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
        elif strand_radius is None:
            self.strand_radius = np.sqrt(conductor_radius**2*fill_factor/number_strands)
        else:
            raise Exception("1 of the 4 parameters need to be None.")

        self.n_layers = ff.NbrLayers(number_strands)
        self.a_cell = self.n_strands * self.strand_radius ** 2 * np.pi / self.ff

        print(f"Updated Litz Configuration: \n"
              f" ff: {self.ff} \n"
              f" Number of layers/strands: {self.n_layers}/{self.n_strands} \n"
              f" Strand radius: {self.strand_radius} \n"
              f" Conductor radius: {self.conductor_radius}\n"
              f"---")

    """
    def to_dict(self):
        conductor_settings = {
            "conductor_type": self.conductor_type.name
        }
        if self.conductor_type in [ConductorType.Foil, ConductorType.Full, ConductorType.Stacked]:
            conductor_settings["thickness"] = self.thickness
            conductor_settings["wrap_para"] = self.wrap_para
        elif self.conductor_type == ConductorType.Litz:
            conductor_settings["conductor_radius"] = self.conductor_radius
            conductor_settings["n_strands"] = self.n_strands
            conductor_settings["strand_radius"] = self.strand_radius
            conductor_settings["ff"] = self.ff
        elif self.conductor_type == ConductorType.Solid:
            conductor_settings["conductor_radius"] = self.conductor_radius
        else:
            raise Exception(f"Unknown conductor type {self.conductor_type}")

        contents = {
            "turns_primary": self.turns_primary,
            "turns_secondary": self.turns_secondary,
            "conductivity": self.conductivity.name,
            "winding_type": self.winding_type.name,
            "winding_scheme": self.winding_scheme.name,
            "conductor_settings": conductor_settings
        }

        return contents
    """

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
    core_w: float           # Axi symmetric case | core_w := 2x core radius
    core_h: float
    window_w: float         # Winding window width
    window_h: float         # Winding window height
    core_type: str = "EI"   # Basic shape of magnetic conductor
    
    steinmetz_loss: int = 0
    generalized_steinmetz_loss: int = 0
    
    # TODO Does this represent the number of windows the EI core has?
    number_core_windows: int

    # Needed for to_dict
    loss_approach: LossApproach = None

    # TODO explanation
    r_inner: float
    r_outer: float

    correct_outer_leg: bool

    def __init__(self, core_w: float, window_w: float, window_h: float, material: str = "custom",  # "95_100" 
                   loss_approach: LossApproach = LossApproach.LossAngle, mu_rel: float = 3000,
                   phi_mu_deg: float = None, sigma: float = None, non_linear: bool = False, correct_outer_leg: bool = False, **kwargs):
        # Set parameters
        self.core_w = core_w
        self.core_h = None  # TODO Set core_h to not none
        self.window_w = window_w
        self.window_h = window_h
        self.type = "axi_symmetric"
        self.material = material
        self.non_linear = non_linear
        self.mu_rel = mu_rel
        self.phi_mu_deg = phi_mu_deg

        self.loss_approach = loss_approach

        self.number_core_windows = 2
        self.correct_outer_leg = correct_outer_leg

        self.r_inner = window_w + core_w / 2
        if correct_outer_leg:
            A_out = 200 * 10 ** -6
            self.r_outer = np.sqrt(A_out / np.pi + self.r_inner ** 2)  # Hardcode for PQ 40/40
        else:
            self.r_outer = np.sqrt((core_w / 2) ** 2 + self.r_inner ** 2)

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

        # Needed because of to_dict
        self.kwargs = kwargs

    def to_dict(self):
        return {
            "core_w": self.core_w,
            "window_w": self.window_w, 
            "window_h": self.window_h, 
            "material": self.material, 
            "loss_approach": self.loss_approach.name,
            "mu_rel": self.mu_rel,
            "phi_mu_deg": self.phi_mu_deg, 
            "sigma": self.sigma, 
            "non_linear": self.non_linear, 
            "kwargs": self.kwargs
        }

class AirGaps:
    """
    Contains methods and arguments to describe the air gaps in a magnetic component

    An air gap can be added with the add_air_gap function. It is possible to set different positions and heights.
    """

    core: Core
    midpoints: List[List[float]]  #: list: [position_tag, air_gap_position, air_gap_h]
    number: int

    # Needed for to_dict
    air_gap_settings: List

    def __init__(self, method: AirGapMethod, core: Core):
        self.method = method
        self.core = core
        self.midpoints = []
        self.number = 0
        self.air_gap_settings = []

    def add_air_gap(self, leg_position: AirGapLegPosition, position_value: Optional[float], height: float):
        """
        Brings a single air gap to the core.

        :param leg_posistion: CenterLeg, OuterLeg
        :type leg_position: AirGapLegPosition
        :param position_value: if AirGapMethod == Percent: 0...100, elif AirGapMethod == Manually: position hight in [m]
        :type position_value: float
        :param height: Air gap height in [m]
        :type height: float
        """
        self.air_gap_settings.append({
            "leg_position": leg_position.name, 
            "position_value": position_value,
            "height": height})

        for index, midpoint in enumerate(self.midpoints):
            if midpoint[0] == leg_position and midpoint[1] + midpoint[2] < position_value - height \
                    and midpoint[1] - midpoint[2] > position_value + height:
                raise Exception(f"Air gaps {index} and {len(self.midpoints)} are overlapping")

        if leg_position == AirGapLegPosition.LeftLeg or leg_position == AirGapLegPosition.RightLeg:
            raise Exception("Currently the legpositions LeftLeg and RightLeg are not supported")

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
            if position_value > 100 or position_value < 0:
                raise Exception("AirGap position values for the percent method need to be between 0 and 100.")
            position = position_value / 100 * self.core.window_h - self.core.window_h / 2

            # When the position is above the winding window it needs to be adjusted
            if position + height / 2 > self.core.window_h / 2:
                position -= (position + height / 2) - self.core.window_h / 2
            elif position - height / 2 < -self.core.window_h / 2:
                position += -self.core.window_h / 2 - (position - height / 2) 

            self.midpoints.append([leg_position.value, position, height])
            self.number += 1
        else:
            raise Exception(f"Method {self.method} is not supported.")

    def to_dict(self):
        content = {
            "method": self.method.name,
            "air_gap_number": len(self.air_gap_settings)
        }

        if self.number > 0:
            content["air_gaps"] = self.air_gap_settings

        return content

class Isolation:
    """
    This class defines isolations for the model.
    An isolation between the winding window and the core can always be set.
    When having a inductor only the primary2primary isolation is necessary.
    When having a (integrated) transformer secondary2secondary and primary2secondary isolations can be set as well.
    """

    cond_cond: List[float] = []
    core_cond: List[float] = []

    def add_winding_isolations(self, primary2primary, secondary2secondary = 0, primary2secondary = 0):
        if primary2primary is None:
            primary2primary = 0
        if secondary2secondary is None:
            secondary2secondary = 0
        if primary2secondary is None:
            primary2secondary = 0

        self.cond_cond = [primary2primary, secondary2secondary, primary2secondary]

    def add_core_isolations(self, top_core, bot_core, left_core, right_core):
        if top_core is None:
            top_core = 0
        if bot_core is None:
            bot_core = 0
        if left_core is None:
            left_core = 0
        if right_core is None:
            right_core = 0

        self.core_cond = [top_core, bot_core, left_core, right_core]

    def to_dict(self):
        return {
            "winding_isolations": self.cond_cond,
            "core_isolations": self.core_cond
        }

@dataclass
class StrayPath:
    """
    This class is needed when an integrated transformer shall be created.

    TODO: Thickness of the stray path must be fitted for the real Tablet (effective area of the
    "stray air gap" is different in axi-symmetric approximation
    """

    start_index: int        # Air gaps are sorted from lowest to highest. This index refers to the air_gap index bottom up
    length: float           # Resembles the length of the whole tablet starting from the y-axis

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

    winding_type: WindingType
    winding_scheme: WindingScheme # Or InterleavedWindingScheme in case of an interleaved winding

    windings: List[Conductor]
    turns: List[int]

    winding_is_set: bool

    def __init__(self, bot_bound: float, top_bound: float, left_bound: float, right_bound: float):
        self.bot_bound = bot_bound
        self.top_bound = top_bound
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.winding_is_set = False

    def set_winding(self, conductor: Conductor, turns: int, winding_scheme: WindingScheme):
        self.winding_type = WindingType.Single
        self.winding_scheme = winding_scheme
        self.windings = [conductor]
        self.turns = [turns]
        self.winding_is_set = True

    def set_interleaved_winding(self, conductor1: Conductor, turns1: int, conductor2: Conductor, turns2: int, winding_scheme: InterleavedWindingScheme):
        self.winding_type = WindingType.Interleaved
        self.winding_scheme = winding_scheme
        self.windings = [conductor1, conductor2]
        self.turns = [turns1, turns2]
        self.winding_is_set = True

    def __repr__(self):
        return f"WindingType: {self.winding_type}, WindingScheme: {self.winding_scheme}| Bounds: bot: {self.bot_bound}, top: {self.top_bound}, left: {self.left_bound}, right: {self.right_bound}"

class WindingWindow:
    max_bot_bound: float
    max_top_bound: float
    max_left_bound: float
    max_right_bound: float
    
    isolation_vww: float

    virtual_winding_windows: List[VirtualWindingWindow]

    def __init__(self, core, isolation):
        
        self.max_bot_bound = -core.window_h / 2 + isolation.core_cond[0]
        self.max_top_bound = core.window_h / 2 - isolation.core_cond[1]
        self.max_left_bound = core.core_w / 2 + isolation.core_cond[2]
        self.max_right_bound = core.r_inner - isolation.core_cond[3]

        # Isolation between vwws
        self.isolation_vww = isolation.cond_cond[2]

    def split_window(self, split_type: WindingWindowSplit, horizontal_split_factor: float = 0, vertical_split_factor: float = 0):
        # Calculate split lengths
        horizontal_split = self.max_left_bound + (self.max_right_bound - self.max_left_bound) * horizontal_split_factor
        vertical_split = self.max_top_bound - abs(self.max_bot_bound - self.max_top_bound) * vertical_split_factor

        # Check for every possible split type and return the corresponding VirtualWindingWindows
        if split_type == WindingWindowSplit.NoSplit:
            complete = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                                top_bound = self.max_top_bound,
                                                left_bound = self.max_left_bound,
                                                right_bound = self.max_right_bound)
            self.virtual_winding_windows = [complete]
            return complete
        elif split_type == WindingWindowSplit.HorizontalSplit:
            right = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                                top_bound = self.max_top_bound,
                                                left_bound = vertical_split + self.isolation_vww / 2,
                                                right_bound = self.max_right_bound)

            left = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                                top_bound = self.max_top_bound,
                                                left_bound = self.max_left_bound,
                                                right_bound = vertical_split - self.isolation_vww / 2)
            self.virtual_winding_windows = [left, right]
            return left, right
        elif split_type == WindingWindowSplit.VerticalSplit:
            top = VirtualWindingWindow(bot_bound = horizontal_split + self.isolation_vww / 2,
                                                top_bound = self.max_top_bound,
                                                left_bound = self.max_left_bound,
                                                right_bound = self.max_right_bound)

            bot = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                                top_bound = horizontal_split - self.isolation_vww / 2,
                                                left_bound = self.max_left_bound,
                                                right_bound = self.max_right_bound)
            self.virtual_winding_windows = [top, bot]
            return top, bot
        elif split_type == WindingWindowSplit.HorizontalAndVerticalSplit:
            top_left = VirtualWindingWindow(bot_bound = horizontal_split + self.isolation_vww / 2,
                                                top_bound = self.max_top_bound,
                                                left_bound = self.max_left_bound,
                                                right_bound = vertical_split - self.isolation_vww / 2)

            top_right = VirtualWindingWindow(bot_bound = horizontal_split + self.isolation_vww / 2,
                                                top_bound = self.max_top_bound,
                                                left_bound = vertical_split + self.isolation_vww / 2,
                                                right_bound = self.max_right_bound)

            bot_left = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                                top_bound = horizontal_split - self.isolation_vww / 2,
                                                left_bound = self.max_left_bound,
                                                right_bound = vertical_split - self.isolation_vww / 2)

            bot_right = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                                top_bound = horizontal_split - self.isolation_vww / 2,
                                                left_bound = vertical_split + self.isolation_vww / 2,
                                                right_bound = self.max_right_bound)
            self.virtual_winding_windows = [top_left, top_right, bot_left, bot_right]
            return top_left, top_right, bot_left, bot_right
        else:
            raise Exception(f"Winding window split type {split_type} not found")

    def combine_vww(self, vww1, vww2):
        index1 = self.virtual_winding_windows.index(vww1)
        index2 = self.virtual_winding_windows.index(vww2)

        if abs(index2-index1) == 3:
            raise Exception("Cannot combine top left and bottom right.")

        self.virtual_winding_windows.remove(vww1)
        self.virtual_winding_windows.remove(vww2)

        new_vww = VirtualWindingWindow(bot_bound = min(vww1.bot_bound, vww2.bot_bound), 
                                    top_bound = max(vww1.top_bound, vww2.top_bound), 
                                    left_bound = min(vww1.left_bound, vww2.left_bound), 
                                    right_bound = max(vww1.right_bound, vww2.right_bound))

        self.virtual_winding_windows.append(new_vww)

        return new_vww

class TwoDaxiSymmetric:
    """
    This class creates the model of the magnetic component. This is done by creating lists of points which
    will be later added to the gmsh file.

    :return:

    """

    core: Core
    virtual_winding_windows: List[VirtualWindingWindow]
    air_gaps: AirGaps
    stray_path: StrayPath
    isolation: Isolation
    component_type: ComponentType
    mesh_data: MeshData
    number_of_windings: int

    # List of points which represent the model
    # Every List is a List of 4 Points: x, y, z, mesh_factor
    p_outer: List[List[float]]
    p_region_bound: List[List[float]] 
    p_window: List[List[float]] 
    p_air_gaps: List[List[float]]
    p_conductor: List[List[float]]
    p_iso_core: List[List[float]]
    p_iso_pri_sec: List[List[float]]

    def __init__(self, core: Core, mesh_data: MeshData, air_gaps: AirGaps, virtual_winding_windows: List[VirtualWindingWindow], 
                stray_path: StrayPath, isolation: Isolation, component_type: ComponentType, number_of_windings: int):
        self.core = core
        self.mesh_data = mesh_data
        self.virtual_winding_windows = virtual_winding_windows
        self.air_gaps = air_gaps
        self.component_type = component_type
        self.stray_path = stray_path
        self.isolation = isolation
        self.number_of_windings = number_of_windings

        # -- Arrays for geometry data -- 
        self.p_outer = np.zeros((4, 4))
        self.p_region_bound = np.zeros((4, 4))
        self.p_window = np.zeros((4 * core.number_core_windows, 4))
        self.p_air_gaps = np.zeros((4 * air_gaps.number, 4)) 
        self.p_conductor = []
        self.p_iso_core = []
        self.p_iso_pri_sec = []

        for i in range(number_of_windings):
            self.p_conductor.insert(i, [])

        self.r_inner = core.r_inner
        self.r_outer = core.r_outer

    def draw_outer(self):
        """
        Draws the outer points

        :return:
        """
        # Outer Core
        # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
        # TODO Case core_h is not None
        if self.core.core_h is None:
            self.p_outer[0][:] = [-self.r_outer,
                                -(self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]

            self.p_outer[1][:] = [self.r_outer,
                                -(self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]

            self.p_outer[2][:] = [-self.r_outer,
                                (self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]

            self.p_outer[3][:] = [self.r_outer,
                                (self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]
        else:
            self.p_outer[0][:] = [-self.r_outer,
                                -self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

            self.p_outer[1][:] = [self.r_outer,
                                -self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

            self.p_outer[2][:] = [-self.r_outer,
                                self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

            self.p_outer[3][:] = [self.r_outer,
                                self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

    def draw_window(self):
        # Window
        # At this point both windows (in a cut) are modeled
        # print(f"win: c_window: {self.component.mesh.c_window}")
        self.p_window[0] = [-self.r_inner,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[1] = [-self.core.core_w / 2,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[2] = [-self.r_inner,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[3] = [-self.core.core_w / 2,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[4] = [self.core.core_w / 2,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[5] = [self.r_inner,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[6] = [self.core.core_w / 2,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[7] = [self.r_inner,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

    def draw_air_gaps(self):
        # Air gaps
        # "air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
        #   - position_tag: specifies the gapped "leg"
        #   - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
        #   - air_gap_h: height/length of the air gap
        #   - c_air_gap: mesh accuracy factor
        # at this point the 4 corner points of each air gap are generated out of "air_gaps"

        mesh_accuracy = self.core.window_w / 20 * self.mesh_data.global_accuracy

        for i in range(0, self.air_gaps.number):

            # # Left leg (-1)
            # if self.component.air_gaps.midpoints[i][0] == -1:
            #     self.p_air_gaps[i * 4] = [-(self.component.core.core_w + self.component.core.window_w),
            #                               self.component.air_gaps.midpoints[i][1] -
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #     self.p_air_gaps[i * 4 + 1] = [-(self.component.core.core_w / 2 + self.component.core.window_w),
            #                                   self.component.air_gaps.midpoints[i][1] -
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #     self.p_air_gaps[i * 4 + 2] = [-(self.component.core.core_w + self.component.core.window_w),
            #                                   self.component.air_gaps.midpoints[i][1] +
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #     self.p_air_gaps[i * 4 + 3] = [-(self.component.core.core_w / 2 + self.component.core.window_w),
            #                                   self.component.air_gaps.midpoints[i][1] +
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #
            # # Right leg (+1)
            # if self.component.air_gaps.midpoints[i][0] == 1:
            #     self.p_air_gaps[i * 4] = [self.component.core.core_w / 2 + self.component.core.window_w,
            #                               self.component.air_gaps.midpoints[i][1] -
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #     self.p_air_gaps[i * 4 + 1] = [self.component.core.core_w + self.component.core.window_w,
            #                                   self.component.air_gaps.midpoints[i][1] -
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #     self.p_air_gaps[i * 4 + 2] = [self.component.core.core_w / 2 + self.component.core.window_w,
            #                                   self.component.air_gaps.midpoints[i][1] +
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
            #     self.p_air_gaps[i * 4 + 3] = [self.component.core.core_w + self.component.core.window_w,
            #                                   self.component.air_gaps.midpoints[i][1] +
            #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]

            # Center leg (0)
            if self.air_gaps.midpoints[i][0] == 0:
                # The center points are transformed each into 4 corner points

                air_gap_y_position = self.air_gaps.midpoints[i][1]
                air_gap_height = self.air_gaps.midpoints[i][2]
                air_gap_length_top = self.core.core_w / 2
                air_gap_length_bot = self.core.core_w / 2

                # Check for stray_paths in integrated transformers
                if self.component_type == ComponentType.IntegratedTransformer:
                    if self.stray_path.start_index == i:
                        # Stray path is above current air_gap
                        air_gap_length_top = self.stray_path.length
                    elif self.stray_path.start_index + 1 == i:
                        # Stray path is below current air_gap
                        air_gap_length_bot = self.stray_path.length

                # Bottom left
                self.p_air_gaps[i * 4 + 0] = [0,
                                                air_gap_y_position -
                                                air_gap_height / 2,
                                                0,
                                                self.mesh_data.c_core]

                # Bottom right
                self.p_air_gaps[i * 4 + 1] = [air_gap_length_bot,
                                                air_gap_y_position -
                                                air_gap_height / 2,
                                                0,
                                                mesh_accuracy]

                # Top left
                self.p_air_gaps[i * 4 + 2] = [0,
                                                air_gap_y_position +
                                                air_gap_height / 2,
                                                0,
                                                self.mesh_data.c_core]

                # Top right
                self.p_air_gaps[i * 4 + 3] = [air_gap_length_top,
                                                air_gap_y_position +
                                                air_gap_height / 2,
                                                0,
                                                mesh_accuracy]

        # In order to close the air gap when a stray_path is added, additional points need to be added
        if self.component_type == ComponentType.IntegratedTransformer:
            top_point = [self.core.core_w / 2,
                            self.air_gaps.midpoints[self.stray_path.start_index+1][1] -
                            air_gap_height / 2,
                            0,
                            mesh_accuracy]
            bot_point = [self.core.core_w / 2,
                            self.air_gaps.midpoints[self.stray_path.start_index][1] +
                            air_gap_height / 2,
                            0,
                            mesh_accuracy]
            self.p_close_air_gaps = [top_point, bot_point]

    def draw_conductors(self, draw_top_down = True):
        # Draw every conductor type based on the virtual_winding_window bounds

        for virtual_winding_window in self.virtual_winding_windows:
            # Get bounds from virtual winding window
            bot_bound = virtual_winding_window.bot_bound
            top_bound = virtual_winding_window.top_bound
            left_bound = virtual_winding_window.left_bound
            right_bound = virtual_winding_window.right_bound

            # Check the possible WindingTypes and draw accordingly
            if virtual_winding_window.winding_type == WindingType.Interleaved:
                # Two windings in the virtual winding window
                windings = virtual_winding_window.windings
                winding0 = windings[0]
                winding1 = windings[1]

                turns = virtual_winding_window.turns
                turns0 = turns[0]
                turns1 = turns[1]

                winding_numbers = [winding0.winding_number, winding1.winding_number]

                # Now check for every winding scheme which is possible for interleaved winding type
                if virtual_winding_window.winding_scheme == InterleavedWindingScheme.Bifilar:
                    """
                    - Bifilar interleaving means a uniform winding scheme of two conductors (prim. and sec.)
                    - Can only be used for conductors of identical radius (in terms of litz radius for 
                        stranded wires)
                    - Excess windings are placed below the bifilar ones
    
                    """
                    if winding0.conductor_radius != winding1.conductor_radius:
                        print("For bifilar winding scheme both conductors must be of the same radius!")
                    else:
                        print("Bifilar winding scheme is applied")

                    raise Exception("Bifilar winding scheme is not implemented yet.")

                if virtual_winding_window.winding_scheme == InterleavedWindingScheme.Vertical:
                    """
                    - Vertical interleaving means a winding scheme where the two conductors are alternating 
                        in vertical (y-)direction
                    - This is practically uncommon
                    - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                        conductor
                    """

                    raise Exception("Vertical winding scheme is not implemented yet.")

                if virtual_winding_window.winding_scheme == InterleavedWindingScheme.Horizontal:
                    """
                    - Horizontal interleaving means a winding scheme where the two conductors are alternating in 
                    horizontal  (x-)direction (Tonnenwicklung)
                    - This is practically most common
                    - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                        conductor
                        
                    """

                    # assume 2 winding transformer and dedicated stray path:
                    if draw_top_down:
                        # This will draw from top bound down

                        # Initialize the list, that counts the already placed conductors
                        N_completed = [0, 0]

                        # Initialize the starting conductor
                        if turns0 >= turns1:
                            # Primary starts first
                            col_cond = 0
                        else:
                            # Secondary starts fist
                            col_cond = 1

                        # Get winding number of current winding
                        winding_number = winding_numbers[col_cond]

                        # Initialize the x and y coordinate
                        x = left_bound + windings[col_cond].conductor_radius
                        y = top_bound - windings[col_cond].conductor_radius
                        top_window_iso_counter = 0
                        # Continue placing as long as not all conductors have been placed
                        while (turns0 - N_completed[0] != 0) or \
                                (turns1 - N_completed[1] != 0):
                            if turns[col_cond] - N_completed[col_cond] != 0:
                                # is this winding not already finished?
                                if x < right_bound - windings[col_cond].conductor_radius:
                                    while y > bot_bound + windings[col_cond].conductor_radius and \
                                            N_completed[col_cond] < turns[col_cond]:
                                        self.p_conductor[winding_number].append([
                                            x, 
                                            y, 
                                            0, 
                                            self.mesh_data.c_center_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x - windings[col_cond].conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x,
                                            y + windings[col_cond].conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x + windings[col_cond].conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([x,
                                                                            y - windings[col_cond].conductor_radius,
                                                                            0,
                                                                            self.mesh_data.c_conductor[winding_number]])

                                        N_completed[col_cond] += 1

                                        y -= windings[col_cond].conductor_radius * 2 + self.isolation.cond_cond[col_cond]  # one from bot to top

                                    x += windings[col_cond].conductor_radius + \
                                            windings[(col_cond + 1) % 2].conductor_radius + \
                                            self.isolation.cond_cond[2]  # from left to right

                                    # Reset y
                                    col_cond = (col_cond + 1) % 2
                                    winding_number = winding_numbers[col_cond]
                                    y = top_bound - windings[col_cond].conductor_radius
                                    top_window_iso_counter += 1
                                else:
                                    break

                            else:
                                # is this winding already finished? - continue with the other one
                                col_cond = (col_cond + 1) % 2

                                # Correct the reset of y and correct x displacement
                                x += windings[col_cond].conductor_radius - \
                                        windings[(col_cond + 1) % 2].conductor_radius \
                                        - self.isolation.cond_cond[2] + self.isolation.cond_cond[
                                            col_cond]

                                y = top_bound - windings[col_cond].conductor_radius
                                top_window_iso_counter -= 1
                    else:
                        # This will draw from bottom bound up.

                        # Initialize the list, that counts the already placed conductors
                        N_completed = [0, 0]

                        # Initialize the starting conductor
                        if turns0 >= turns1:
                            col_cond = 0
                        else:
                            col_cond = 1

                        # Get winding number of current winding
                        winding_number = winding_numbers[col_cond]

                        # Initialize the x and y coordinate
                        x = left_bound + windings[col_cond].conductor_radius
                        y = bot_bound + windings[col_cond].conductor_radius

                        # Continue placing as long as not all conductors have been placed
                        while (turns0 - N_completed[0] != 0) or \
                                (turns1 - N_completed[1] != 0):
                            if turns[col_cond] - N_completed[col_cond] != 0:
                                # is this winding not already finished?
                                if x < right_bound - windings[col_cond].conductor_radius:
                                    while y < top_bound - windings[col_cond].conductor_radius and \
                                            N_completed[col_cond] < turns[col_cond]:
                                        self.p_conductor[winding_number].append([
                                            x, 
                                            y, 
                                            0, 
                                            self.mesh_data.c_center_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x - windings[col_cond].conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x,
                                            y + windings[col_cond].conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x + windings[col_cond].conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        self.p_conductor[winding_number].append([
                                            x,
                                            y - windings[col_cond].conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number]])

                                        N_completed[col_cond] += 1

                                        y += windings[col_cond].conductor_radius * 2 + \
                                                self.isolation.cond_cond[col_cond]  # one from bot to top

                                    x += windings[col_cond].conductor_radius + \
                                            windings[(col_cond + 1) % 2].conductor_radius + \
                                            self.isolation.cond_cond[2]  # from left to right

                                    # Reset y
                                    col_cond = (col_cond + 1) % 2
                                    y = bot_bound + windings[col_cond].conductor_radius

                                else:
                                    break

                            else:
                                # is this winding already finished? - continue with the other one
                                col_cond = (col_cond + 1) % 2
                                winding_number = winding_numbers[col_cond]

                                # Correct the reset of y and correct x displacement
                                x += windings[col_cond].conductor_radius - \
                                        windings[(col_cond + 1) % 2].conductor_radius \
                                        - self.isolation.cond_cond[2] + \
                                        self.isolation.cond_cond[
                                            col_cond]

                                y = bot_bound + windings[col_cond].conductor_radius
                
                if virtual_winding_window.winding_scheme == InterleavedWindingScheme.VerticalStacked:
                    winding_number0, winding_number1 = winding_numbers

                    # First winding from bottom to top
                    if winding0.conductor_arrangement == ConductorArrangement.Square:
                        while y < top_bound - winding0.conductor_radius and \
                                i < turns0:
                            while x < right_bound - winding0.conductor_radius and \
                                    i < turns0:
                                self.p_conductor[winding_number0].append([
                                    x, 
                                    y, 
                                    0, 
                                    self.mesh_data.c_center_conductor[winding_number0]])
                                self.p_conductor[winding_number0].append([
                                    x - winding0.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number0]])
                                self.p_conductor[winding_number0].append([
                                    x,
                                    y + winding0.conductor_radius,
                                    0,
                                    self.mesh_data.c_conductor[winding_number0]])
                                self.p_conductor[winding_number0].append([
                                    x + winding0.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number0]])
                                self.p_conductor[num].append([
                                    x, 
                                    y - winding0.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                x += winding0.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from left to right
                            y += winding0.conductor_radius * 2 + \
                                    self.isolation.cond_cond[
                                        num]  # one step from bot to top
                            x = left_bound + winding0.conductor_radius  # always the same
                    elif winding0.conductor_arrangement == ConductorArrangement.Hexagonal:
                            base_line = True

                            while y < top_bound - winding0.conductor_radius and \
                                    i < turns0:
                                while x < right_bound - winding0.conductor_radius and \
                                        i < turns0:
                                    self.p_conductor[winding_number0].append([
                                        x, 
                                        y, 
                                        0, 
                                        self.mesh_data.c_center_conductor[winding_number0]])
                                    self.p_conductor[winding_number0].append([
                                        x - winding0.conductor_radius, 
                                        y, 
                                        0,
                                        self.mesh_data.c_conductor[winding_number0]])
                                    self.p_conductor[winding_number0].append([
                                        x, 
                                        y + winding0.conductor_radius, 
                                        0,
                                        self.mesh_data.c_conductor[winding_number0]])
                                    self.p_conductor[winding_number0].append([
                                        x + winding0.conductor_radius, 
                                        y, 
                                        0,
                                        self.mesh_data.c_conductor[winding_number0]])
                                    self.p_conductor[winding_number0].append([
                                        x, 
                                        y - winding0.conductor_radius, 
                                        0,
                                        self.mesh_data.c_conductor[winding_number0]])
                                    i += 1

                                    x += 2 * np.cos(np.pi / 6) * (
                                            winding0.conductor_radius +
                                            self.isolation.cond_cond[num] / 2)

                                    # depending on what line, hexa scheme starts shifted
                                    # reset y to "new" bottom
                                    base_line = (not base_line)
                                    if base_line:
                                        y -= (winding0.conductor_radius +
                                                self.isolation.cond_cond[num])
                                    else:
                                        y += (winding0.conductor_radius +
                                                self.isolation.cond_cond[num])

                                # Undo last base_line reset
                                if base_line:
                                    y += (winding0.conductor_radius +
                                            self.isolation.cond_cond[num])
                                else:
                                    y -= (winding0.conductor_radius +
                                            self.isolation.cond_cond[num])

                                base_line = True
                                x = left_bound + winding0.conductor_radius
                                y += winding0.conductor_radius + \
                                        self.component.isolation.cond_cond[num]
                    else:
                        raise Exception(f"Unknown conductor_arrangement {winding0.conductor_arrangement}")

                    # Second winding from top to bottom
                    if winding1.conductor_arrangement == ConductorArrangement.Square:
                        while y > bot_bound + winding1.conductor_radius and i < \
                                turns1:
                            while x < right_bound - winding1.conductor_radius and i < \
                                    turns1:
                                self.p_conductor[winding_number1].append([
                                    x, 
                                    y, 
                                    0, 
                                    self.mesh_data.c_center_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x - winding1.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x, 
                                    y + winding1.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x + winding1.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x, 
                                    y - winding1.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                i += 1

                                x += winding1.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from left to right
                            y += -(winding1.conductor_radius * 2) - \
                                    self.isolation.cond_cond[
                                        num]  # one step from bot to top
                            x = left_bound + winding1.conductor_radius  # always the same
                    elif winding1.conductor_arrangement == ConductorArrangement.Hexagonal:
                        base_line = True

                        while y > bot_bound + winding1.conductor_radius and \
                                i < turns1:
                            while x < right_bound - winding1.conductor_radius and \
                                    i < turns1:
                                print(f"i: {i} "
                                        f"x: {x} "
                                        f"y: {y} ")

                                self.p_conductor[winding_number1].append([
                                    x, 
                                    y, 
                                    0, 
                                    self.mesh_data.c_center_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x - winding1.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x, 
                                    y + winding1.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x + winding1.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])
                                self.p_conductor[winding_number1].append([
                                    x, 
                                    y - winding1.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[winding_number1]])

                                i += 1
                                x += 2 * np.cos(np.pi / 6) * (
                                        winding1.conductor_radius +
                                        self.isolation.cond_cond[num] / 2)

                                # depending on what line, hexa scheme starts shifted
                                # reset y to "new" bottom
                                base_line = (not base_line)
                                if base_line:
                                    y += (winding1.conductor_radius +
                                            self.isolation.cond_cond[num])
                                else:
                                    y -= (winding1.conductor_radius +
                                            self.isolation.cond_cond[num])

                            # Undo last base_line reset
                            if base_line:
                                y -= (winding1.conductor_radius +
                                        self.isolation.cond_cond[num])
                            else:
                                y += (winding1.conductor_radius +
                                        self.isolation.cond_cond[num])

                            base_line = True
                            x = left_bound + winding1.conductor_radius
                            # from top to bottom
                            y -= (winding1.conductor_radius +
                                    self.isolation.cond_cond[num])
                    else:
                        raise Exception(f"Unknown conductor_arrangement {winding1.conductor_arrangement}")
            elif virtual_winding_window.winding_type == WindingType.Single:
                # One winding in the virtual winding window
                winding = virtual_winding_window.windings[0]
                turns = virtual_winding_window.turns[0]
                conductor_type = winding.conductor_type
                conductor_arrangement = winding.conductor_arrangement
                winding_scheme = virtual_winding_window.winding_scheme

                num = winding.winding_number

                # Check if the coil is round or rectangular
                if conductor_type == ConductorType.RectangularSolid:    
                    # Now check for each possible winding scheme 
                    if winding_scheme == WindingScheme.Full:
                        # Full window conductor
                        self.p_conductor[num].append([
                            left_bound, 
                            bot_bound, 
                            0, 
                            self.mesh_data.c_conductor[num]])
                        self.p_conductor[num].append([
                            right_bound, 
                            bot_bound, 
                            0, 
                            self.mesh_data.c_conductor[num]])
                        self.p_conductor[num].append([
                            left_bound, 
                            top_bound, 
                            0, 
                            self.mesh_data.c_conductor[num]])
                        self.p_conductor[num].append([
                            right_bound, 
                            top_bound, 
                            0, 
                            self.mesh_data.c_conductor[num]])
                    elif winding_scheme == WindingScheme.FoilVertical:
                        # Foil conductors where each conductor is very high and the conductors are expanding in the x-direction
                        if winding.wrap_para == WrapParaType.Fixed_Thickness:
                            # Wrap defined number of turns and chosen thickness
                            for i in range(turns):
                                # CHECK if right bound is reached
                                if (left_bound + (i + 1) * winding.thickness +
                                    i * self.isolation.cond_cond[num]) <= right_bound:
                                    # Foils
                                    self.p_conductor[num].append([
                                        left_bound + i * winding.thickness + i * self.isolation.cond_cond[num],
                                        bot_bound, 
                                        0, 
                                        self.mesh_data.c_conductor[num]])
                                    self.p_conductor[num].append([
                                        left_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num], 
                                        bot_bound, 
                                        0,
                                        self.mesh_data.c_conductor[num]])
                                    self.p_conductor[num].append([
                                        left_bound + i * winding.thickness + i * self.isolation.cond_cond[num],
                                        top_bound, 
                                        0, 
                                        self.mesh_data.c_conductor[num]])
                                    self.p_conductor[num].append([
                                        left_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num], 
                                        top_bound, 
                                        0,
                                        self.mesh_data.c_conductor[num]])
                        elif winding.wrap_para == WrapParaType.Interpolate:
                            # Fill the allowed space in the Winding Window with a chosen number of turns
                            x_interpol = np.linspace(left_bound, right_bound + self.isolation.cond_cond[num],
                                                        winding.turns + 1)
                            for i in range(turns):
                                # Foils
                                self.p_conductor[num].append([
                                    x_interpol[i], 
                                    bot_bound, 
                                    0, 
                                    self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x_interpol[i + 1] - self.component.isolation.cond_cond[num], 
                                    bot_bound, 
                                    0,
                                    self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x_interpol[i],
                                    top_bound, 
                                    0, 
                                    self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x_interpol[i + 1] - self.component.isolation.cond_cond[num], 
                                    top_bound, 
                                    0,
                                    self.component.mesh.c_conductor[num]])
                        else:
                            raise Exception(f"Unknown wrap para type {winding.wrap_para}")
                    elif winding_scheme == WindingScheme.FoilHorizontal:
                        # Foil conductors where each conductor is very long and the conductors are expanding the y-direction
                        # Stack defined number of turns and chosen thickness
                        for i in range(turns):
                            # CHECK if top bound is reached
                            if (bot_bound + (i + 1) * winding.thickness +
                                i * self.isolation.cond_cond[num]) <= top_bound:
                                # stacking from the ground
                                self.p_conductor[num].append([
                                    left_bound, 
                                    bot_bound + i * winding.thickness + i * self.isolation.cond_cond[num], 
                                    0, 
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    bot_bound + i * winding.thickness + i * self.isolation.cond_cond[num],
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    left_bound,
                                    bot_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num], 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    bot_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num],
                                    0,
                                    self.mesh_data.c_conductor[num]])      
                    else:
                        raise Exception(f"Winding scheme {winding_scheme} is not implemented.")
                elif conductor_type == ConductorType.RoundSolid or conductor_type == ConductorType.RoundLitz:
                    # Since round conductors have no winding scheme check for each conductor_arrangement
                    if conductor_arrangement == ConductorArrangement.Square:
                        y = bot_bound + winding.conductor_radius
                        x = left_bound + winding.conductor_radius
                        i = 0
                        # Case n_conductors higher that "allowed" is missing
                        while x < right_bound - winding.conductor_radius and i < turns:
                            while y < top_bound - winding.conductor_radius and i < turns:
                                self.p_conductor[num].append([
                                    x, 
                                    y, 
                                    0, 
                                    self.mesh_data.c_center_conductor[num]])
                                self.p_conductor[num].append([
                                    x - winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x, 
                                    y + winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append(
                                    [x + winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x, 
                                    y - winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                y += winding.conductor_radius * 2 + self.isolation.cond_cond[num]  # one step from left to right
                            x += winding.conductor_radius * 2 + self.isolation.cond_cond[num]  # from left to top
                            y = bot_bound + winding.conductor_radius
                    elif conductor_arrangement == ConductorArrangement.Hexagonal:
                        y = bot_bound + winding.conductor_radius
                        x = left_bound + winding.conductor_radius
                        i = 0
                        base_line = True
                        # Case n_conductors higher that "allowed" is missing
                        while x < right_bound - winding.conductor_radius \
                                and i < turns:
                            while y < top_bound - winding.conductor_radius and \
                                    i < turns:
                                self.p_conductor[num].append([
                                    x, 
                                    y, 
                                    0, 
                                    self.mesh_data.c_center_conductor[num]])
                                self.p_conductor[num].append([
                                    x - winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x, 
                                    y + winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x + winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x, 
                                    y - winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                y += winding.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from bottom to top
                            x += 2 * np.cos(np.pi / 6) * (
                                    winding.conductor_radius +
                                    self.isolation.cond_cond[num] / 2)
                            # * np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from left to right
                            # depending on what line, hexa scheme starts shifted
                            # reset y to "new" bottom
                            base_line = (not base_line)
                            if base_line:
                                y = bot_bound + winding.conductor_radius
                            else:
                                y = bot_bound + 2 * winding.conductor_radius + \
                                    self.isolation.cond_cond[
                                        num] / 2
                    elif conductor_arrangement == ConductorArrangement.SquareFullWidth:
                        y = bot_bound + winding.conductor_radius
                        x = left_bound + winding.conductor_radius
                        i = 0
                        # Case n_conductors higher that "allowed" is missing
                        while y < top_bound - winding.conductor_radius \
                                and i < turns:
                            while x < right_bound - winding.conductor_radius \
                                    and i < turns:
                                self.p_conductor[num].append([
                                    x, 
                                    y, 
                                    0, 
                                    self.mesh_data.c_center_conductor[num]])
                                self.p_conductor[num].append([
                                    x - winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x,
                                    y + winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x + winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x, 
                                    y - winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                x += winding.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from left to top
                            y += winding.conductor_radius * 2 + \
                                    self.isolation.cond_cond[
                                        num]  # one step from left to right
                            x = left_bound + winding.conductor_radius  # always the same
                    else:
                        raise Exception(f"Conductor arrangement {conductor_arrangement} is not implemented.")
                else:
                    raise Exception(f"Conductor shape {winding.conductor_shape} is not implemented.")
            else:
                raise Exception(f"Unknown winding type {virtual_winding_window.winding_type}")

            # Checking the Conductors
            for winding in virtual_winding_window.windings:
                num = winding.winding_number

                # Convert to numpy
                # Check if all Conductors could be resolved
                self.p_conductor[num] = np.asarray(self.p_conductor[num])

                # TODO:CHECKS for rect. conductors
                """ CHECK: rectangle conductors with 4 points
                if self.component.windings[num].conductor_type == "full" or 
                        self.component.windings[num].conductor_type == "stacked" or \
                        self.component.windings[num].conductor_type == "foil":
                    if int(self.p_conductor[num].shape[0]/4) < self.component.windings[num].turns:
                        warnings.warn("Too many turns that do not fit in the winding window.")
                        # self.component.windings[num].turns = int(self.p_conductor[num].shape[0]/4)
                        self.component.valid = None
                """

                # CHECK: round conductors with 5 points
                if winding.conductor_type in [ConductorType.RoundSolid, ConductorType.RoundLitz]:
                    if int(self.p_conductor[num].shape[0] / 5) < virtual_winding_window.turns[num]:
                        # Warning: warnings.warn("Too many turns that do not fit in the winding window.")
                        # Correct: self.component.windings[num].turns = int(self.p_conductor[num].shape[0]/5)
                        # TODO: break, but remove warning. valid bit should be set to False
                        #  Code must go to the next parameter-iteration step for geometric sweep
                        self.component.valid = False
                        # TODO Tell the user which winding window
                        raise Exception(f"Too many turns that do not fit in the winding window {str(virtual_winding_window)}")

            # Region for Boundary Condition
            self.p_region_bound[0][:] = [-self.r_outer * self.mesh_data.padding,
                                        -(self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]
            self.p_region_bound[1][:] = [self.r_outer * self.mesh_data.padding,
                                        -(self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]
            self.p_region_bound[2][:] = [-self.r_outer * self.mesh_data.padding,
                                        (self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]
            self.p_region_bound[3][:] = [self.r_outer * self.mesh_data.padding,
                                        (self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]

    def draw_isolations(self, isolation_deltas):
        """
        DISCLAIMER
        Because the distance from the core to the winding is set by
        iso.core_cond, a delta, which is used in order for no overlapping lines will cause
        the "real" isolation to be slightly smaller than set by the user.
        """

        window_h = self.core.window_h
        iso = self.isolation
        mesh = self.mesh_data

        # Since there are many cases in which alternating conductors would lead to slightly different
        # mesh densities a simplification is made: Just use the lowest mesh density to be safe all the time.
        mesh_density_to_winding = min(mesh.c_conductor)

        mesh_density_to_core = mesh.c_window
            

        # Using the delta the lines and points from the isolation and the core/windings are not overlapping
        # which makes creating the mesh more simpler
        # Isolation between winding and core
        iso_core_delta_left = isolation_deltas["core_left"]
        iso_core_delta_top = isolation_deltas["core_top"]
        iso_core_delta_right = isolation_deltas["core_right"]
        iso_core_delta_bot = isolation_deltas["core_bot"]
        iso_iso_delta = isolation_deltas["iso_iso"]
        iso_winding_delta_left = isolation_deltas["winding_left"]
        iso_winding_delta_top = isolation_deltas["winding_top"]
        iso_winding_delta_right = isolation_deltas["winding_right"]
        iso_winding_delta_bot = isolation_deltas["winding_bot"]

        self.p_iso_core = [] # Order: Left, Top, Right, Bot
        self.p_iso_pri_sec = []

        if self.component_type == ComponentType.IntegratedTransformer:
            # TODO implement for integrated_transformers
            # TODO Change back to warnings?
            print("Isolations are not set because they are not implemented for integrated transformers.")
        else:
            # Useful points for isolation creation
            left_x = self.core.core_w / 2 + iso_core_delta_left
            top_y = window_h / 2 - iso_core_delta_top
            right_x = self.r_inner - iso_core_delta_right
            bot_y = -window_h / 2 + iso_core_delta_bot

            # Useful lengths
            left_iso_width = iso.core_cond[2] - iso_core_delta_left - iso_winding_delta_left
            top_iso_height = iso.core_cond[0] - iso_core_delta_top - iso_winding_delta_top
            right_iso_width = iso.core_cond[3] - iso_core_delta_right - iso_winding_delta_right
            bot_iso_height = iso.core_cond[1] - iso_core_delta_bot - iso_winding_delta_bot

            # Core to Pri isolation
            iso_core_left = [
                [
                    left_x,
                    top_y - top_iso_height - iso_iso_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x + left_iso_width,
                    top_y - top_iso_height - iso_iso_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x + left_iso_width,
                    bot_y + bot_iso_height + iso_iso_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    bot_y + bot_iso_height + iso_iso_delta,
                    0,
                    mesh_density_to_core
                ]
            ]
            iso_core_top = [
                [
                    left_x,
                    top_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    top_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    top_y - top_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    top_y - top_iso_height,
                    0,
                    mesh_density_to_winding
                ]
            ]
            iso_core_right = [
                [
                    right_x - right_iso_width,
                    top_y - top_iso_height - iso_iso_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    top_y - top_iso_height - iso_iso_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    bot_y + bot_iso_height + iso_iso_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x - right_iso_width,
                    bot_y + bot_iso_height + iso_iso_delta,
                    0,
                    mesh_density_to_winding
                ]
            ]
            iso_core_bot = [
                [
                    left_x,
                    bot_y + bot_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    bot_y + bot_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    bot_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x,
                    bot_y,
                    0,
                    mesh_density_to_core
                ]
            ]

            self.p_iso_core = [iso_core_left, iso_core_top, iso_core_right, iso_core_bot]

            """
            Since the way virtual winding windows are handled is changed this will be changed too

            # Isolation between virtual winding windows
            if self.component.vw_type == VirtualWindingType.FullWindow:
                # Only one vww -> The winding can be interleaved 
                # -> still an isolation between pri and sec necessary
                if len(self.component.virtual_winding_windows) == 1 and self.component.virtual_winding_windows[0].winding == WindingType.Interleaved:
                    # vertical isolations needed between the layers
                    # bifilar and vertical do not exist yet
                    vww = self.component.virtual_winding_windows[0]
                    winding_0 = self.component.windings[0]
                    winding_1 = self.component.windings[1]
                    current_x = vww.left_bound + 2 * winding_0.conductor_radius

                    mesh_density_left_side = 0
                    mesh_density_right_side = 0
                    if winding_0.turns[0] > winding_1.turns[0]:
                        mesh_density_left_side = mesh.c_conductor[0]
                        mesh_density_right_side = mesh.c_conductor[1]
                    else:
                        mesh_density_left_side = mesh.c_conductor[1]
                        mesh_density_right_side = mesh.c_conductor[0]

                    if vww.scheme == WindingScheme.Horizontal:
                        self.p_iso_pri_sec = []
                        for index in range(self.top_window_iso_counter-1):
                            self.p_iso_pri_sec.append([
                                [
                                    current_x + iso_winding_delta_left,
                                    top_y - top_iso_height - iso_iso_delta,
                                    0,
                                    mesh_density_left_side
                                ],
                                [
                                    current_x + iso.cond_cond[2] - iso_winding_delta_right,
                                    top_y - top_iso_height - iso_iso_delta,
                                    0,
                                    mesh_density_right_side
                                ],
                                [
                                    current_x + iso.cond_cond[2] - iso_winding_delta_right,
                                    bot_y + bot_iso_height + iso_iso_delta,
                                    0,
                                    mesh_density_right_side
                                ],
                                [
                                    current_x + iso_winding_delta_left,
                                    bot_y + bot_iso_height + iso_iso_delta,
                                    0,
                                    mesh_density_left_side
                                ]
                            ])
                            # The sec winding can start first when it has a higher number of turns
                            # col_cond
                            if index % 2 == self.col_cond_start:
                                current_x += iso.cond_cond[2] + 2 * winding_1.conductor_radius
                            else:
                                current_x += iso.cond_cond[2] + 2 * winding_0.conductor_radius
                    elif vww.scheme == WindingScheme.Vertical:
                        raise Exception("Vertical scheme not implemented yet!")
                    elif vww.scheme == WindingScheme.Bifilar:
                        raise Exception("Bifilar scheme not implemented yet!")
                    else:
                        raise Exception(f"The winding scheme {vww.scheme} is unknown.")
            elif self.component.vw_type == VirtualWindingType.Split2:
                # Two vwws -> a horizontal isolation is needed
                vww_bot = self.component.virtual_winding_windows[0]
                vww_top = self.component.virtual_winding_windows[1]
                self.p_iso_pri_sec.append([
                    [
                        vww_top.left_bound,
                        vww_top.bot_bound - iso_winding_delta_top,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        vww_top.right_bound,
                        vww_top.bot_bound - iso_winding_delta_top,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        vww_top.right_bound,
                        vww_bot.top_bound + iso_winding_delta_bot,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        vww_top.left_bound,
                        vww_bot.top_bound + iso_winding_delta_bot,
                        0,
                        mesh_density_to_winding
                    ]
                ])
            else:
                warnings.warn(f"Isolations are not implemented for components with type {self.component.vw_type}")
            """

    def draw_model(self, isolation_deltas = None):

        self.draw_outer()

        self.draw_window()

        self.draw_air_gaps()

        self.draw_conductors()

        if isolation_deltas is None:
            isolation_deltas = {
                "core_left": 0.00001,
                "core_top": 0.00001,
                "core_bot": 0.00001,
                "core_right": 0.00001,
                "iso_iso" : 0.00001,
                "winding_left": 0.00001,
                "winding_top": 0.00001,
                "winding_right": 0.00001,
                "winding_bot": 0.00001
            }

        self.draw_isolations(isolation_deltas)