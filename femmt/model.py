# Python standard libraries
from dataclasses import dataclass
from typing import List, Tuple, Optional, Union

# 3rd party libraries
import numpy as np

# Local libraries
import femmt.functions as ff
from femmt.functions_model import *
import femmt.functions_reluctance as fr
from femmt.enumerations import *
from femmt.constants import *
from femmt.functions_drawing import *
import materialdatabase as mdb


class Conductor:
    """
    A winding defines a conductor which is wound around a magnetic component such as transformer or inductance.
    The winding is defined by its conductor and the way it is placed in the magnetic component. To allow different
    arrangements of the conductors in several winding windows (hexagonal or square packing, interleaved, ...) in
    this class only the conductor parameters are specified. 
    """
    # TODO More documentation
    conductor_type: ConductorType
    conductor_arrangement: ConductorArrangement = None
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

    conductor_is_set: bool

    # Not used in femmt_classes. Only needed for to_dict()
    conductivity: Conductivity = None

    def __init__(self, winding_number: int, conductivity: Conductivity, parallel: bool = False,
                 winding_material_temperature: float = 100):
        """Creates a conductor object.
        The winding_number sets the order of the conductors. Every conductor needs to have a unique winding number.
        The conductor with the lowest winding number (starting from 0) will be treated as primary, second-lowest number as secondary and so on.

        :param winding_number: Unique number for the winding
        :type winding_number: int
        :param conductivity: Sets the conductivity for the conductor
        :type conductivity: float
        :param winding_material_temperature: temperature of winding material, default set to 100 °C
        :type winding_material_temperature: float
        :param parallel: Set to True to introduce parallel conductors. Default set to False
        :type parallel: bool
        """
        if winding_number < 0:
            raise Exception("Winding index cannot be negative.")

        self.winding_number = winding_number
        self.conductivity = conductivity
        self.conductor_is_set = False
        self.parallel = parallel

        dict_material_database = ff.wire_material_database()
        if conductivity.name in dict_material_database:
            self.cond_sigma = ff.conductivity_temperature(conductivity.name, winding_material_temperature)
        else:
            raise Exception(f"Material {conductivity.name} not found in database")

    def set_rectangular_conductor(self, thickness: float):
        if self.conductor_is_set:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_is_set = True
        self.conductor_type = ConductorType.RectangularSolid
        self.thickness = thickness
        self.a_cell = None  # can only be set after the width is determined
        self.conductor_radius = 1  # Revisit

    def set_solid_round_conductor(self, conductor_radius: float, conductor_arrangement: ConductorArrangement):
        if self.conductor_is_set:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_is_set = True
        self.conductor_type = ConductorType.RoundSolid
        self.conductor_arrangement = conductor_arrangement
        self.conductor_radius = conductor_radius
        self.a_cell = np.pi * conductor_radius ** 2

    def set_litz_round_conductor(self, conductor_radius: Optional[float], number_strands: Optional[int],
                                 strand_radius: Optional[float],
                                 fill_factor: Optional[float], conductor_arrangement: ConductorArrangement):
        """
        Only 3 of the 4 parameters are needed. The other one needs to be none
        """
        if self.conductor_is_set:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_is_set = True
        self.conductor_type = ConductorType.RoundLitz
        self.conductor_arrangement = conductor_arrangement
        self.conductor_radius = conductor_radius
        self.n_strands = number_strands
        self.strand_radius = strand_radius
        self.ff = fill_factor

        if number_strands is None:
            self.n_strands = int(conductor_radius ** 2 / strand_radius ** 2 * fill_factor)
        elif conductor_radius is None:
            self.conductor_radius = np.sqrt(number_strands * strand_radius ** 2 / fill_factor)
        elif fill_factor is None:
            ff_exact = number_strands * strand_radius ** 2 / conductor_radius ** 2
            self.ff = np.around(ff_exact, decimals=2)
            if self.ff > 0.99:
                raise Exception(f"A fill factor of {self.ff} is unrealistic!")
        elif strand_radius is None:
            self.strand_radius = np.sqrt(conductor_radius ** 2 * fill_factor / number_strands)
        else:
            raise Exception("1 of the 4 parameters need to be None.")

        self.n_layers = ff.litz_calculate_number_layers(number_strands)
        self.a_cell = self.n_strands * self.strand_radius ** 2 * np.pi / self.ff

        ff.femmt_print(f"Updated Litz Configuration: \n"
                       f" ff: {self.ff} \n"
                       f" Number of layers/strands: {self.n_layers}/{self.n_strands} \n"
                       f" Strand radius: {self.strand_radius} \n"
                       f" Conductor radius: {self.conductor_radius}\n"
                       f"---")

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return self.__dict__ != other.__dict__

    def to_dict(self):
        return {
            "winding_number": self.winding_number,
            "conductivity": self.conductivity.name,
            "conductor_type": self.conductor_type.name,
            "thickness": self.thickness,
            "conductor_radius": self.conductor_radius,
            "conductor_arrangement": self.conductor_arrangement.name if self.conductor_arrangement is not None else None,
            "number_strands": self.n_strands,
            "strand_radius": self.strand_radius,
            "fill_factor": self.ff
        }


class Core:
    """
    This creates the core base for the model.
    # TODO More documentation and get rid of double initializations
    frequency = 0: mu_r_abs only used if non_linear == False
    frequency > 0: mu_r_abs is used
    """


    # Standard material data
    material: str

    # Permeability
    # TDK N95 as standard material:
    permeability_type: PermeabilityType
    mu_r_abs: float  # Relative Permeability [if complex: mu_complex = re_mu_rel + j*im_mu_rel with mu_rel=|mu_complex|]
    phi_mu_deg: float  # mu_complex = mu_r_abs * exp(j*phi_mu_deg)

    # Permitivity - [Conductivity in a magneto-quasistatic sense]
    sigma: float  # Imaginary part of complex equivalent permittivity [frequency-dependent]

    steinmetz_loss: int = 0
    generalized_steinmetz_loss: int = 0

    # Needed for to_dict
    loss_approach: LossApproach = None

    # Database
    # material_database is variable to load in material_database
    temperature: float  # temperature at which data is required
    material: str  # material to be accessed from data base
    datasource: str  # type of data to be accessed ( datasheet or measurement)

    file_path_to_solver_folder: str  # location to create temporary pro file

    def __init__(self,
                 # dimensions
                 core_type: CoreType = CoreType.Single,
                 core_dimensions=None,
                 correct_outer_leg: bool = False,

                 # material data
                 material: str = "custom",
                 temperature: float = None,
                 loss_approach: LossApproach = LossApproach.LossAngle,
                 mu_r_abs: float = 3000,

                 # permeability
                 permeability_datasource: MaterialDataSource = None,
                 permeability_datatype: MeasurementDataType = None,
                 permeability_measurement_setup: str = None,
                 phi_mu_deg: float = None,
                 non_linear: bool = False,

                 # permittivity
                 permittivity_datasource: str = None,
                 permittivity_datatype: str = None,
                 permittivity_measurement_setup: str = None,
                 sigma: float = None,
                 steinmetz_parameter: list = None,
                 generalized_steinmetz_parameter: list = None,
                 **kwargs):
        """TODO Doc

        :param core_inner_diameter: diameter of the inner core
        :type core_inner_diameter: float
        :param window_w: width of the winding window
        :type window_w: float
        :param window_h: height of the winding window
        :type window_h: float
        :param material: Material name, e.g. 'N95', defaults to "custom"mu_rel
        :type material: str, optional
        :param mu_r_abs: absolute value of the permeability, defaults to 3000
        :type mu_r_abs: float, optional
        :param phi_mu_deg: loss angle of the material in degree, defaults to None
        :type phi_mu_deg: float, optional
        :param sigma: core conductivity, defaults to None
        :type sigma: float, optional
        :param non_linear: _description_, defaults to False
        :type non_linear: bool, optional
        :param correct_outer_leg: Manual correction so cross-section of inner leg is not same as outer leg (PQ 40/40 only!!!), defaults to False (recommended!)
        :type correct_outer_leg: bool, optional
        """
        # Set parameters
        self.core_type = core_type  # Basic shape of magnetic conductor

        # Geometric Parameters
        self.core_inner_diameter = core_dimensions.core_inner_diameter
        self.window_w = core_dimensions.window_w
        self.correct_outer_leg = correct_outer_leg
        self.r_inner = self.window_w + self.core_inner_diameter / 2

        if self.core_type == CoreType.Single:
            self.window_h = core_dimensions.window_h
            self.number_core_windows = 2
            self.core_h = self.window_h + self.core_inner_diameter / 2  # TODO: could also be done arbitrarily
        if self.core_type == CoreType.Stacked:
            self.window_h_bot = core_dimensions.window_h_bot
            self.window_h_top = core_dimensions.window_h_top
            # self.core_h = self.window_h_bot + self.window_h_top + self.core_inner_diameter * 3 / 4  # TODO: could also be done arbitrarily
            self.core_h = self.window_h_bot + self.core_inner_diameter / 2  # TODO: could also be done arbitrarily
            self.number_core_windows = 4

        if correct_outer_leg:
            # hard-coded case for PQ 40/40-cores:
            # the outer cross-section differs from inner cross-section and is corrected here.
            # Note: for PQ 40/40 cores only!!
            A_out = 200 * 10 ** -6
            self.r_outer = np.sqrt(A_out / np.pi + self.r_inner ** 2)  # Hardcode for PQ 40/40
        else:
            # set r_outer, so cross-section of outer leg has same cross-section as inner leg
            # this is the recommended default-case
            self.r_outer = fr.calculate_r_outer(self.core_inner_diameter, self.window_w)

        # Material Parameters
        # General
        # Initialize database
        self.material_database = mdb.MaterialDatabase(ff.silent)
        self.material = material
        self.file_path_to_solver_folder = None
        self.temperature = temperature
        # Permeability
        self.permeability = {"datasource": permeability_datasource,
                             "measurement_setup": permeability_measurement_setup,
                             "datatype": permeability_datatype}
        self.non_linear = non_linear
        self.mu_r_abs = mu_r_abs
        self.phi_mu_deg = phi_mu_deg
        self.loss_approach = loss_approach
        # Permittivity
        self.complex_permittivity = None
        self.permittivity = {"datasource": permittivity_datasource,
                             "measurement_setup": permittivity_measurement_setup,
                             "datatype": permittivity_datatype}
        self.sigma = sigma

        # Check loss approach
        if loss_approach == LossApproach.Steinmetz:
            self.sigma = 0
            if self.material != "custom":
                self.permeability_type = PermeabilityType.FromData
                self.mu_r_abs = self.material_database.get_material_attribute(material_name=self.material,
                                                                             attribute="initial_permeability")

                steinmetz_data = self.material_database.get_steinmetz_data(material_name=self.material, type="Steinmetz",
                                                              datasource="measurements")
                self.ki = steinmetz_data['ki']
                self.alpha = steinmetz_data['alpha']
                self.beta = steinmetz_data['beta']

            if self.material == "custom":  # ----steinmetz_parameter consist of list of ki, alpha , beta from the user
                self.ki = steinmetz_parameter[0]
                self.alpha = steinmetz_parameter[1]
                self.beta = steinmetz_parameter[2]
                ff.femmt_print(f"{self.ki, self.alpha, self.beta = }")
            else:
                raise Exception(f"When steinmetz losses are set a material needs to be set as well.")
        # if loss_approach == LossApproach.Generalized_Steinmetz:
        #     raise NotImplemented
            # self.sigma = 0
            # if self.material != "custom":
            #     self.permeability_type = PermeabilityType.FromData
            #     self.mu_rel = Database.get_initial_permeability(material_name=self.material)
            #     self.t_rise = Database.get_steinmetz_data(material_name=self.material, type="Generalized_Steinmetz")[
            #         't_rise']
            #     self.t_fall = Database.get_steinmetz_data(material_name=self.material, type="Generalized_Steinmetz")[
            #         't_fall']
            # elif self.material == "custom":  # ----generalized_steinmetz_parameter consist of list of ki, alpha , beta from the user
            #     self.t_rise = generalized_steinmetz_parameter[0]
            #     self.t_fall = generalized_steinmetz_parameter[1]

        if loss_approach == LossApproach.LossAngle:

            if self.permittivity["datasource"] == MaterialDataSource.Custom:
                self.sigma = sigma  # from user
            else:
                self.sigma = 1 / self.material_database.get_material_attribute(material_name=self.material, attribute="resistivity")

            if self.permeability["datasource"] == MaterialDataSource.Custom:
                 # this is a service for the user:
                # In case of not giving a fixed loss angle phi_mu_deg,
                # the initial permeability from datasheet is used (RealValue)
                if phi_mu_deg is not None and phi_mu_deg != 0:
                    self.permeability_type = PermeabilityType.FixedLossAngle
                else:
                    self.permeability_type = PermeabilityType.RealValue
            else:
                self.permeability_type = PermeabilityType.FromData
                self.mu_r_abs = self.material_database.get_material_attribute(material_name=self.material, attribute="initial_permeability")

        else:
            raise Exception("Loss approach {loss_approach.value} is not implemented")

        # Set attributes of core with given keywords
        # TODO Should we allow this? Technically this is not how an user interface should be designed
        for key, value in kwargs.items():
            setattr(self, key, value)

        # Needed because of to_dict
        self.kwargs = kwargs

    def update_sigma(self, frequency: bool) -> None:
        """
        Updates the core conductivity.
        The core conductivity is used to calculate eddy current losses inside the FEM simulation.

        In case of datasheet parameters, the DC-resistance is loaded from the material database.
        In case of measurement parameters, the AC-resistance is loaded from the material database.
         * The AC-resistance is interpolated by the material database for the certain frequency
         * AC-resistance = imaginary part of complex permittivity

        The AC-resistance decreases when increasing the frequency.
        Using AC-resistance (interpolated from measurements) is more accurate than the static DC-datasheet resistance.

        The conductivity 'sigma' = (1 / AC-resistance) is used for the simulation.

        :param frequency: operating frequency in Hz
        :type frequency: float
        """
        if self.permittivity["datasource"] == MaterialDataSource.Measurement:
            epsilon_r, phi_epsilon_deg = self.material_database.get_permittivity(temperature=self.temperature, frequency=frequency, material_name=self.material,
                                                                                 datasource=self.permittivity["datasource"],
                                                                                 measurement_setup=self.permittivity["measurement_setup"],
                                                                                 datatype=self.permittivity["datatype"])
            self.complex_permittivity = epsilon_0 * epsilon_r * complex(np.cos(np.deg2rad(phi_epsilon_deg)), np.sin(np.deg2rad(phi_epsilon_deg)))
            self.sigma = complex(2 * np.pi * frequency * self.complex_permittivity.imag, 2 * np.pi * frequency * self.complex_permittivity.real)

        if self.permittivity["datasource"] == MaterialDataSource.ManufacturerDatasheet:
            self.sigma = 1 / self.material_database.get_material_attribute(material_name=self.material, attribute="resistivity")

    def update_core_material_pro_file(self, frequency, electro_magnetic_folder, plot_interpolation: bool = False):
        # This function is needed to update the pro file for the solver depending on the frequency of the
        # upcoming simulation
        ff.femmt_print(f"{self.permeability['datasource'] = }")
        self.material_database.permeability_data_to_pro_file(temperature=self.temperature, frequency=frequency,
                                                             material_name=self.material,
                                                             datasource=self.permeability["datasource"],
                                                             datatype=self.permeability["datatype"],
                                                             measurement_setup=self.permeability["measurement_setup"],
                                                             parent_directory=electro_magnetic_folder,
                                                             plot_interpolation=plot_interpolation)

    def to_dict(self):
        # TODO: mdb
        if self.core_type == CoreType.Single:
            return {
                "core_inner_diameter": self.core_inner_diameter,
                "window_w": self.window_w,
                "window_h": self.window_h,
                "material": self.material,
                "loss_approach": self.loss_approach.name,
                "mu_r_abs": self.mu_r_abs,
                "phi_mu_deg": self.phi_mu_deg,
                "sigma": self.sigma,
                "non_linear": self.non_linear,
                "correct_outer_leg": self.correct_outer_leg,
                "temperature": self.temperature,
                "permeability_datasource": self.permeability["datasource"],
                "permeability_measurement_setup": self.permeability["measurement_setup"],
                "permeability_datatype": self.permeability["datatype"],
                "permittivity_datasource": self.permittivity["datasource"],
                "permittivity_measurement_setup": self.permittivity["measurement_setup"],
                "permittivity_datatype": self.permittivity["datatype"]
            }

        elif self.core_type == CoreType.Stacked:
            return {
                "core_inner_diameter": self.core_inner_diameter,
                "window_w": self.window_w,
                "window_h_bot": self.window_h_bot,
                "window_h_top": self.window_h_top,
                "material": self.material,
                "loss_approach": self.loss_approach.name,
                "mu_r_abs": self.mu_r_abs,
                "phi_mu_deg": self.phi_mu_deg,
                "sigma": self.sigma,
                "non_linear": self.non_linear,
                "correct_outer_leg": self.correct_outer_leg,
                "temperature": self.temperature,
                "permeability_datasource": self.permeability["datasource"],
                "permeability_measurement_setup": self.permeability["measurement_setup"],
                "permeability_datatype": self.permeability["datatype"],
                "permittivity_datasource": self.permittivity["datasource"],
                "permittivity_measurement_setup": self.permittivity["measurement_setup"],
                "permittivity_datatype": self.permittivity["datatype"]
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
        """Creates an AirGaps object. An AirGapMethod needs to be set. This determines the way the air gap will be added to the model.
        In order to calculate the air gap positions the core object needs to be given.

        :param method: The method determines the way the air gap position is set.
        :type method: AirGapMethod
        :param core: The core object
        :type core: Core
        """
        self.method = method
        self.core = core
        self.midpoints = []
        self.number = 0
        self.air_gap_settings = []

    def add_air_gap(self, leg_position: AirGapLegPosition, height: float, position_value: Optional[float] = 0, stacked_position: StackedPosition = None):
        """
        Brings a single air gap to the core.

        :param leg_position: CenterLeg, OuterLeg
        :type leg_position: AirGapLegPosition
        :param position_value: if AirGapMethod == Percent: 0...100, elif AirGapMethod == Manually: position height in [m]
        :type position_value: float
        :param height: Air gap height in [m]
        :type height: float
        :param stacked_position: Top, Bot
        :type stacked_position: StackedPosition
        """
        self.air_gap_settings.append({
            "leg_position": leg_position.name,
            "position_value": position_value,
            "height": height,
            "stacked_position": stacked_position})

        for index, midpoint in enumerate(self.midpoints):
            if midpoint[0] == leg_position and midpoint[1] + midpoint[2] < position_value - height \
                    and midpoint[1] - midpoint[2] > position_value + height:
                raise Exception(f"Air gaps {index} and {len(self.midpoints)} are overlapping")

        if leg_position == AirGapLegPosition.LeftLeg or leg_position == AirGapLegPosition.RightLeg:
            raise Exception("Currently the leg positions LeftLeg and RightLeg are not supported")

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

        elif self.method == AirGapMethod.Stacked:
            # set midpoints
            # TODO: handle top and bot
            if stacked_position == StackedPosition.Bot:
                self.midpoints.append([0, 0, height])
                self.number += 1
            if stacked_position == StackedPosition.Top:
                self.midpoints.append([0, self.core.window_h_bot/2 + self.core.core_inner_diameter / 4 + height / 2, height])  # TODO: could also be done arbitrarily
                self.number += 1

        else:
            raise Exception(f"Method {self.method} is not supported.")

    def to_dict(self):
        if self.number == 0:
            return {}

        content = {
            "method": self.method.name,
            "air_gap_number": len(self.air_gap_settings)
        }

        if self.number > 0:
            content["air_gaps"] = self.air_gap_settings

        return content


class Insulation:
    """
    This class defines insulation for the model.
    An insulation between the winding window and the core can always be set.
    When having an inductor only the primary2primary insulation is necessary.
    When having a (integrated) transformer secondary2secondary and primary2secondary insulations can be set as well.

    Only the isolation between winding window and core is drawn as a "physical" isolation (4 rectangles). All other isolations
    are only describing a set distance between the object.

    In general, it is not necessary to add an insulation object at all when no insulation is needed.
    """
    cond_cond: List[List[float]]  # two-dimensional list with size NxN, where N is the number of windings (symmetrical isolation matrix)
    core_cond: List[float]  # list with size 4x1, with respectively isolation of cond_n -> [top_core, bot_core, left_core, right_core]

    insulation_delta: float

    def __init__(self):
        """Creates an insulation object.

        Sets an insulation_delta value. In order to simplify the drawing of the isolations between core and winding window the isolation rectangles
        are not exactly drawn at the specified position. They are slightly smaller and the offset can be changed with the insulation_delta variable.
        In general, it is not recommended to change this value.
        """
        # Default value for all insulations
        # If the gaps between insulations and core (or windings) are to big/small just change this value
        self.insulation_delta = 0.00001

    def add_winding_insulations(self, inner_winding_insulation: List[List[float]]):
        """Adds insulations between turns of one winding and insulation between virtual winding windows.
        Insulation between virtual winding windows is not always needed.
        :param inner_winding_insulation: List of floats which represent the insulations between turns of the same winding. This does not correspond to the order conductors are added to the winding! Instead, the winding number is important. The conductors are sorted by ascending winding number. The lowest winding number therefore is combined with index 0. The second lowest with index 1 and so on.
        :type inner_winding_insulation: List[List[float]]
        """
        if inner_winding_insulation is [[]]:
            raise Exception("Inner winding insulations list cannot be empty.")

        self.cond_cond = inner_winding_insulation

    def add_core_insulations(self, top_core: float, bot_core: float, left_core: float, right_core: float):
        """Adds insulations between the core and the winding window. Creating those will draw real rectangles in the model.

        :param top_core: Insulation between winding window and top core
        :type top_core: float
        :param bot_core: Insulation between winding window and bottom core
        :type bot_core: float
        :param left_core: Insulation between winding window and left core
        :type left_core: float
        :param right_core: Insulation between winding window and right core
        :type right_core: float
        """
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
        if len(self.cond_cond) == 0:
            return {}

        return {
            "inner_winding_insulations": self.cond_cond,
            "core_insulations": self.core_cond
        }


@dataclass
class StrayPath:
    """
    This class is needed when an integrated transformer shall be created.
    A start_index and a length can be given. The start_index sets the position of the tablet.
    start_index=0 will create the tablet between the lowest and second-lowest air gaps. start_index=1 will create the tablet
    between the second lowest and third-lowest air gap. Therefore, it is necessary for the user to make sure that enough air gaps exist!
    The length parameter sets the length of the tablet starting at the y-axis (not the right side of the center core). It therefore
    determines the air gap between the tablet and the outer core leg.
    """
    # TODO: Thickness of the stray path must be fitted for the real Tablet (effective area of the "stray air gap" is different in axi-symmetric approximation)
    start_index: int  # Air gaps are sorted from lowest to highest. This index refers to the air_gap index bottom up
    length: float  # Resembles the length of the whole tablet starting from the y-axis


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
    winding_scheme: Union[WindingScheme, InterleavedWindingScheme]  # Either WindingScheme or InterleavedWindingScheme (Depending on the winding)
    wrap_para: WrapParaType

    windings: List[Conductor]
    turns: List[int]

    winding_is_set: bool
    winding_insulation: float

    def __init__(self, bot_bound: float, top_bound: float, left_bound: float, right_bound: float):
        """Creates a virtual winding window with given bounds. By default, a virtual winding window is created by the WindingWindow class.
        The parameter values are given in metres and depend on the axisymmetric coordinate system.

        :param bot_bound: Bottom bound
        :type bot_bound: float
        :param top_bound: Top bound
        :type top_bound: float
        :param left_bound: Left bound
        :type left_bound: float
        :param right_bound: Right bound
        :type right_bound: float
        """
        self.bot_bound = bot_bound
        self.top_bound = top_bound
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.winding_is_set = False

    def set_winding(self, conductor: Conductor, turns: int, winding_scheme: WindingScheme,
                    wrap_para_type: WrapParaType = None):
        """Sets a single winding to the current virtual winding window. A single winding always contains one conductor.

        :param conductor: Conductor which will be set to the vww.
        :type conductor: Conductor
        :param turns: Number of turns of the conductor
        :type turns: int
        :param winding_scheme: Winding scheme defines the way the conductor is wrapped. Can be set to None.
        :type winding_scheme: WindingScheme
        :param wrap_para_type: Additional wrap parameter. Not always needed, defaults to None
        :type wrap_para_type: WrapParaType, optional
        """
        self.winding_type = WindingType.Single
        self.winding_scheme = winding_scheme
        self.windings = [conductor]
        self.turns = [0] * (conductor.winding_number + 1)  # TODO: find another soultion for this (is needed in mesh.py for air_stacked)
        # self.turns = [0] * (3)  # TODO: find another soultion for this (is needed in mesh.py for air_stacked)
        self.turns.insert(conductor.winding_number, turns)
        self.winding_is_set = True
        self.wrap_para = wrap_para_type

        if winding_scheme is WindingScheme.FoilVertical and wrap_para_type is None:
            raise Exception("When winding scheme is FoilVertical a wrap para type must be set.")

    def set_interleaved_winding(self, conductor1: Conductor, turns1: int, conductor2: Conductor, turns2: int, winding_scheme: InterleavedWindingScheme):
        """Sets an interleaved winding to the current virtual winding window. An interleaved winding always contains two conductors.
        If a conductor is primary or secondary is determined by the value of the winding number of the conductor. The order of the function parameters
        is irrelevant.

        :param conductor1: Conductor 1 which will be added to the vww. Not equal to primary winding.
        :type conductor1: Conductor
        :param turns1: Turns of conductor 1
        :type turns1: int
        :param conductor2: Conductor 2 which will be added to the vww. Not equal to secondary winding.
        :type conductor2: Conductor
        :param turns2: Turns of conductor 2
        :type turns2: int
        :param winding_scheme: Interleaved winding scheme defines the way the conductors are drawn
        :type winding_scheme: InterleavedWindingScheme
        """
        self.winding_type = WindingType.TwoInterleaved
        self.winding_scheme = winding_scheme
        self.windings = [conductor1, conductor2]
        self.turns = [turns1, turns2]
        self.winding_is_set = True
        self.wrap_para = None

    def set_center_tapped_winding(self,
                                  conductor1: Conductor, turns1: int,
                                  conductor2: Conductor, turns2: int,
                                  conductor3: Conductor, turns3: int,
                                  isolation_primary_to_primary: float,
                                  isolation_secondary_to_secondary: float,
                                  isolation_primary_to_secondary: float):
        # TODO: centertapped is following line allowed to set winding insulation this way?
        # self.winding_insulation = define_center_tapped_insulation(primary_to_primary=2e-4,
        #                                                           secondary_to_secondary=2e-4,
        #                                                           primary_to_secondary=5e-4)
        self.winding_insulation = [isolation_primary_to_primary, isolation_secondary_to_secondary, isolation_primary_to_secondary]
        self.winding_type = WindingType.CenterTappedGroup
        self.winding_scheme = None  # TODO: centertapped maybe add vertical or sth. like this
        self.windings = [conductor1, conductor2, conductor3]
        self.turns = [turns1, turns2, turns3]
        # TODO: centertapped is also a deepcopy needed somewhere: tertiary...?
        self.winding_is_set = True
        self.wrap_para = None

    def __repr__(self):
        return f"WindingType: {self.winding_type}, WindingScheme: {self.winding_scheme}, Bounds: bot: {self.bot_bound}, top: {self.top_bound}, left: {self.left_bound}, right: {self.right_bound}"

    def to_dict(self):
        if hasattr(self, 'winding_insulation'):
            return {
                "bot_bound": self.bot_bound,
                "top_bound": self.top_bound,
                "left_bound": self.left_bound,
                "right_bound": self.right_bound,
                "winding_type": self.winding_type.name,
                "winding_scheme": self.winding_scheme.name if self.winding_scheme is not None else None,
                "wrap_para": self.wrap_para.name if self.wrap_para is not None else None,
                "windings": [winding.to_dict() for winding in self.windings],
                "turns": self.turns,
                "winding_insulation": self.winding_insulation
            }

        else:
            return {
                "bot_bound": self.bot_bound,
                "top_bound": self.top_bound,
                "left_bound": self.left_bound,
                "right_bound": self.right_bound,
                "winding_type": self.winding_type.name,
                "winding_scheme": self.winding_scheme.name if self.winding_scheme is not None else None,
                "wrap_para": self.wrap_para.name if self.wrap_para is not None else None,
                "windings": [winding.to_dict() for winding in self.windings],
                "turns": self.turns,
            }

    # TODO Since in combine_vww it is necessary to compare vwws maybe a __eq__ and __ne__ 
    # function should be implemented.


class WindingWindow:
    """Represents the winding window which is necessary for every model in FEMMT.
    Depending on the type different virtual winding windows are created by this class which then contain the different conductors.
    """
    max_bot_bound: float
    max_top_bound: float
    max_left_bound: float
    max_right_bound: float

    # 4 different insulations which can be Null if there is a vww overlapping
    # The lists contain 4 values x1, y1, x2, y2 where (x1, y1) is the lower left and (x2, y2) the upper right point 
    top_iso: List[float]
    left_iso: List[float]
    bot_iso: List[float]
    right_iso: List[float]

    insulations: Insulation
    split_type: WindingWindowSplit
    stray_path: StrayPath
    air_gaps: AirGaps

    horizontal_split_factor: float
    vertical_split_factor: float

    virtual_winding_windows: List[VirtualWindingWindow]

    def __init__(self, core: Core, insulations: Insulation, stray_path: StrayPath = None, air_gaps: AirGaps = None):
        """Creates a winding window which then creates up to 4 virtual winding windows. In order to correctly calculate the
        virtual winding windows the core, isolations, stray_path and air_gaps objects are needed.

        The stray_path and air_gaps objects are only needed when having an integrated transformer.

        :param core: Core object
        :type core: Core
        :param insulations: Insulation object
        :type insulations: Insulation
        :param stray_path: Stray path object. Only needed for integrated transformer, defaults to None
        :type stray_path: StrayPath, optional
        :param air_gaps: Air gaps path object. Only needed for integrated transformer, defaults to None
        :type air_gaps: AirGaps, optional
        """
        self.core = core
        self.stray_path = stray_path
        self.air_gaps = air_gaps

        if self.core.core_type == CoreType.Single:
            self.max_bot_bound = -core.window_h / 2 + insulations.core_cond[1]
            self.max_top_bound = core.window_h / 2 - insulations.core_cond[0]
            self.max_left_bound = core.core_inner_diameter / 2 + insulations.core_cond[2]
            self.max_right_bound = core.r_inner - insulations.core_cond[3]
        elif self.core.core_type == CoreType.Stacked:  # top, bot, left, right
            self.max_bot_bound = -core.window_h_bot / 2 + insulations.core_cond[1]
            self.max_top_bound = core.window_h_bot / 2 + core.window_h_top + core.core_inner_diameter / 4 - insulations.core_cond[0]
            self.max_left_bound = core.core_inner_diameter / 2 + insulations.core_cond[2]
            self.max_right_bound = core.r_inner - insulations.core_cond[3]

        # Insulations
        self.insulations = insulations

    def to_dict(self):
        return {
            "max_bot_bound": self.max_bot_bound,
            "max_top_bound": self.max_top_bound,
            "max_left_bound": self.max_left_bound,
            "max_right_bound": self.max_right_bound,
            "virtual_winding_windows": [virtual_winding_window.to_dict() for virtual_winding_window in self.virtual_winding_windows],
        }


    def split_window(self, split_type: WindingWindowSplit, split_distance: float = 0,
                     horizontal_split_factor: float = 0.5, vertical_split_factor: float = 0.5,
                     top_bobbin: float = None, bot_bobbin: float = None, left_bobbin: float = None, right_bobbin: float = None) -> Tuple[VirtualWindingWindow]:
        """Creates up to 4 virtual winding windows depending on the split type and the horizontal and vertical split factors.
        The split factors are values between 0 and 1 and determine a horizontal and vertical line at which the window is split.
        Not every value is needed for every split type:
        - NoSplit: No factor is needed
        - HorizontalSplit: Horizontal split factor needed
        - VerticalSplit: Vertical split factor needed
        - HorizontalAndVerticalSplit: Both split factors needed
        
        Up to 4 virtual winding windows are returned:
        - NoSplit: complete
        - HorizontalSplit: left, right
        - VerticalSplit: top, bottom
        - HorizontalAndVerticalSplit: top_left, top_right, bot_left, bot_right

        :param split_type: Determines the arrangement in which virtual winding windows are created
        :type split_type: WindingWindowSplit
        :param horizontal_split_factor: Horizontal split factor, defaults to 0.5
        :type horizontal_split_factor: float, optional
        :param vertical_split_factor: Vertical split factor, defaults to 0.5
        :type vertical_split_factor: float, optional
        :return: Tuple containing the virtual winding windows
        :rtype: Tuple[VirtualWindingWindow]
        """

        self.split_type = split_type

        self.horizontal_split_factor = horizontal_split_factor
        self.vertical_split_factor = vertical_split_factor

        # Calculate split lengths
        if self.stray_path is not None and self.air_gaps is not None and self.air_gaps.number > self.stray_path.start_index:
            air_gap_1_position = self.air_gaps.midpoints[self.stray_path.start_index][1]
            air_gap_2_position = self.air_gaps.midpoints[self.stray_path.start_index + 1][1]
            max_pos = max(air_gap_2_position, air_gap_1_position)
            min_pos = min(air_gap_2_position, air_gap_1_position)
            distance = max_pos - min_pos    # TODO: this is set in accordance to the midpoint of the air gap:
                                            # TODO: should be changed to the core-cond isolation
            horizontal_split = min_pos + distance / 2
            vertical_split = self.max_left_bound + (self.max_right_bound - self.max_left_bound) * vertical_split_factor
            split_distance = distance  # here, the distance between the two vwws is set automatically
        else:
            horizontal_split = self.max_top_bound - abs(self.max_bot_bound - self.max_top_bound) * horizontal_split_factor
            vertical_split = self.max_left_bound + (self.max_right_bound - self.max_left_bound) * vertical_split_factor

        # Check for every possible split type and return the corresponding VirtualWindingWindows
        if split_type == WindingWindowSplit.NoSplit:

            complete = VirtualWindingWindow(
                bot_bound=self.max_bot_bound,
                top_bound=self.max_top_bound,
                left_bound=self.max_left_bound,
                right_bound=self.max_right_bound)

            self.virtual_winding_windows = [complete]
            return complete
        elif split_type == WindingWindowSplit.NoSplitWithBobbin:
            bobbin_def = [top_bobbin, bot_bobbin, left_bobbin, right_bobbin]
            for index, element in enumerate(bobbin_def):
                if element is not None and element > self.insulations.core_cond[index]:
                    bobbin_def[index] = self.insulations.core_cond[index] - element
                else:
                    bobbin_def[index] = 0

            complete = VirtualWindingWindow(
                top_bound=self.max_top_bound + bobbin_def[0],
                bot_bound=self.max_bot_bound - bobbin_def[1],
                left_bound=self.max_left_bound - bobbin_def[2],
                right_bound=self.max_right_bound + bobbin_def[3])
            self.virtual_winding_windows = [complete]
            return complete

        elif split_type == WindingWindowSplit.VerticalSplit:
            right = VirtualWindingWindow(
                bot_bound=self.max_bot_bound,
                top_bound=self.max_top_bound,
                left_bound=vertical_split + split_distance / 2,
                right_bound=self.max_right_bound)

            left = VirtualWindingWindow(
                bot_bound=self.max_bot_bound,
                top_bound=self.max_top_bound,
                left_bound=self.max_left_bound,
                right_bound=vertical_split - split_distance / 2)

            self.virtual_winding_windows = [left, right]
            return left, right
        elif split_type == WindingWindowSplit.HorizontalSplit:
            top = VirtualWindingWindow(
                bot_bound=horizontal_split + split_distance / 2,
                top_bound=self.max_top_bound,
                left_bound=self.max_left_bound,
                right_bound=self.max_right_bound)

            bot = VirtualWindingWindow(
                bot_bound=self.max_bot_bound,
                top_bound=horizontal_split - split_distance / 2,
                left_bound=self.max_left_bound,
                right_bound=self.max_right_bound)

            self.virtual_winding_windows = [top, bot]
            return top, bot
        elif split_type == WindingWindowSplit.HorizontalAndVerticalSplit:
            top_left = VirtualWindingWindow(
                bot_bound=horizontal_split + split_distance / 2,
                top_bound=self.max_top_bound,
                left_bound=self.max_left_bound,
                right_bound=vertical_split - split_distance / 2)

            top_right = VirtualWindingWindow(
                bot_bound=horizontal_split + split_distance / 2,
                top_bound=self.max_top_bound,
                left_bound=vertical_split + split_distance / 2,
                right_bound=self.max_right_bound)

            bot_left = VirtualWindingWindow(
                bot_bound=self.max_bot_bound,
                top_bound=horizontal_split - split_distance / 2,
                left_bound=self.max_left_bound,
                right_bound=vertical_split - split_distance / 2)

            bot_right = VirtualWindingWindow(
                bot_bound=self.max_bot_bound,
                top_bound=horizontal_split - split_distance / 2,
                left_bound=vertical_split + split_distance / 2,
                right_bound=self.max_right_bound)

            self.virtual_winding_windows = [top_left, top_right, bot_left, bot_right]
            return top_left, top_right, bot_left, bot_right
        else:
            raise Exception(f"Winding window split type {split_type} not found")

    def split_with_stack(self, stack: ConductorStack):
        """
        The winding window is splitted according to a ConductorStack dataclass.
        :param stack:
        :return:
        """
        bottom_up = True
        vertical_space_used = 0  # initialize the counter with zero
        self.virtual_winding_windows = []
        winding_scheme_type = []
        if bottom_up:
            for i, row_element in enumerate(stack.order):
                vertical_space_used_last = vertical_space_used
                if type(row_element) == StackIsolation:
                    vertical_space_used += row_element.thickness
                else:
                    additional_bobbin = 0
                    if type(row_element) == ConductorRow:
                        vertical_space_used += row_element.row_height
                        additional_bobbin = row_element.additional_bobbin
                        if row_element.winding_tag == WindingTag.Primary:
                            winding_scheme_type.append(WindingType.Single)
                        elif row_element.winding_tag == WindingTag.Secondary or row_element.winding_tag == WindingTag.Tertiary:
                            winding_scheme_type.append(WindingScheme.FoilHorizontal)

                    elif type(row_element) == CenterTappedGroup:
                        vertical_space_used += get_height_of_group(group=row_element)
                        winding_scheme_type.append(WindingType.CenterTappedGroup)

                    self.virtual_winding_windows.append(
                        VirtualWindingWindow(
                            bot_bound=self.max_bot_bound + vertical_space_used_last,
                            top_bound=self.max_bot_bound + vertical_space_used,
                            left_bound=self.max_left_bound + additional_bobbin,
                            right_bound=self.max_right_bound))

        return self.virtual_winding_windows, winding_scheme_type

    def combine_vww(self, vww1: VirtualWindingWindow, vww2: VirtualWindingWindow) -> VirtualWindingWindow:
        """Combines the borders of two virtual winding windows.

        :param vww1: Virtual winding window 1
        :type vww1: VirtualWindingWindow
        :param vww2: Virtual winding window 2
        :type vww2: VirtualWindingWindow
        :return: Virtual winding window with new bounds
        :rtype: VirtualWindingWindow
        """
        index1 = self.virtual_winding_windows.index(vww1)
        index2 = self.virtual_winding_windows.index(vww2)

        if abs(index2 - index1) == 3:
            raise Exception("Cannot combine top left and bottom right.")
        # TODO add check for top right and bottom left

        self.virtual_winding_windows.remove(vww1)
        self.virtual_winding_windows.remove(vww2)

        new_vww = VirtualWindingWindow(bot_bound=min(vww1.bot_bound, vww2.bot_bound),
                                       top_bound=max(vww1.top_bound, vww2.top_bound),
                                       left_bound=min(vww1.left_bound, vww2.left_bound),
                                       right_bound=max(vww1.right_bound, vww2.right_bound))

        self.virtual_winding_windows.append(new_vww)

        return new_vww

    def NHorizontalSplit(self, horizontal_split_factors: list[float] = None, vertical_split_factors: list[float] = None):
        if vertical_split_factors == None:
            vertical_split_factors = [None] * (len(horizontal_split_factors)+1)

        # Convert horizontal_split_factors to a numpy array
        horizontal_split_factors = np.array(horizontal_split_factors)
        if self.stray_path is not None and self.air_gaps is not None and self.air_gaps.number > self.stray_path.start_index:
            air_gap_1_position = self.air_gaps.midpoints[self.stray_path.start_index][1]
            air_gap_2_position = self.air_gaps.midpoints[self.stray_path.start_index + 1][1]
            max_pos = max(air_gap_2_position, air_gap_1_position)
            min_pos = min(air_gap_2_position, air_gap_1_position)
            distance = max_pos - min_pos  # TODO: this is set in accordance to the midpoint of the air gap:
            # TODO: should be changed to the core-cond isolation
            horizontal_splits = min_pos + distance / 2
            vertical_split = self.max_left_bound + (
                    self.max_right_bound - self.max_left_bound) * vertical_split_factor
            split_distance = distance  # here, the distance between the two vwws is set automatically
        else:
            horizontal_splits = [self.max_left_bound]
            horizontal_splits = horizontal_splits + list(self.max_left_bound + (self.max_right_bound-self.max_left_bound) * horizontal_split_factors)
            horizontal_splits.append(self.max_right_bound)
            print(f"{horizontal_split_factors = }")
            print(f"{horizontal_splits = }")
            print(f"{self.max_left_bound = }")

        # Initialize lists for cells and virtual_winding_windows
        self.cells = []
        self.virtual_winding_windows = []

        # Create the remaining pairs of virtual winding windows
        for i in range(0, len(horizontal_splits)-1):
            if vertical_split_factors[i] != None:
                self.cells.append(VirtualWindingWindow(
                    bot_bound=self.max_bot_bound,
                    top_bound=self.max_top_bound - (self.max_top_bound - self.max_bot_bound) * vertical_split_factors[i],
                    left_bound=horizontal_splits[i],
                    right_bound=horizontal_splits[i+1]))
                self.cells.append(VirtualWindingWindow(
                    bot_bound=self.max_top_bound - (self.max_top_bound - self.max_bot_bound) * vertical_split_factors[i],
                    top_bound=self.max_top_bound,
                    left_bound=horizontal_splits[i],
                    right_bound=horizontal_splits[i+1]))
            else:
                self.cells.append(VirtualWindingWindow(
                    bot_bound=self.max_bot_bound,
                    top_bound=self.max_top_bound,
                    left_bound=horizontal_splits[i],
                    right_bound=horizontal_splits[i+1]))

        # Loop through the number of horizontal splits (plus one to account for the last cells)
        # Append each pair of left and right cells to the virtual_winding_windows list.
        # This creates a list of all virtual winding windows in the sequence: Left, Right, Left, Right, ...
        for i in range(0, len(self.cells)):
            self.virtual_winding_windows.append(self.cells[i])
        # Instead of returning the separate lists of Left_cells and Right_cells,
        # we return the combined list of all virtual winding windows.
        return self.virtual_winding_windows