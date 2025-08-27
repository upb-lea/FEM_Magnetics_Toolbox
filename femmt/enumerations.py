"""Enumeration for FEMMT."""
from enum import IntEnum, Enum

class Verbosity(IntEnum):
    """State of verbosity."""

    # TODO Currently in ToFile the FEMMT console outputs are just suppressed not written to a file
    # TODO Add Verbosity for materialdatabase
    Silent = 1  # No outputs
    ToConsole = 2  # Outputs to console
    ToFile = 3  # Outputs to file

class VisualizationMode(IntEnum):
    """How simulation results are shown in Gmsh."""

    Final = 1   # Only final static result
    Post = 2    # Animate after simulation
    Stream = 3  # Live results during simulation

class WindingTag(IntEnum):
    """Names of windings."""

    Primary = 1
    Secondary = 2
    Tertiary = 3


class WindingWindowSplit(IntEnum):
    """Determines how many virtual winding windows are created by the winding window. Used in Winding window class."""

    NoSplit = 1
    """Virtual winding window same size as winding window
    """
    HorizontalSplit = 2
    """Splits winding window in two virtual winding windows which are separated by a horizontal line
    """
    VerticalSplit = 3
    """Splits winding window in two virtual winding windows which are separated by a vertical line
    """
    HorizontalAndVerticalSplit = 4
    """Splits winding window in four virtual winding windows separated by a horizontal and vertical line
    """
    VerticalStack = 5
    """Based on a vertical stack, the winding window is split into several virtual winding windows
    """
    TenCells_Split = 6
    """The winding window is split into 2x5 virtual winding windows
    """
    NCells_Split = 7
    """The winding window is split into 2xN/2 virtual winding windows
    """
    NoSplitWithBobbin = 8
    """
    Acts like "NoSplit", but takes bobbin geometry instead of core-cond insulation to define the virtual winding window.
    """


class ComponentType(IntEnum):
    """Sets the component type for the whole simulation. Needs to be given to the MagneticComponent on creation."""

    Inductor = 1
    Transformer = 2
    IntegratedTransformer = 3

class SimulationType(IntEnum):
    """Sets the simulation type. The static is just to show the fields."""

    FreqDomain = 1
    TimeDomain = 2
    ElectroStatic = 3

class CoreType(IntEnum):
    """Sets the core type for the whole simulation. Needs to be given to the MagneticComponent on creation."""

    Single = 1  # one axisymmetric core
    Stacked = 2  # one and a half cores

class AirGapMethod(IntEnum):
    """Sets the method how the air gap position (vertical) is set.

    Used in AirGaps class.
    """

    Center = 1
    """Only valid for one air gap. This air gap will always be placed in the middle (vertically)
    """
    Percent = 2
    """A value between 0 and 100 will determine the vertical position.
    """
    Manually = 3
    """The vertical position needs to be given manually. In metres.
    """
    Stacked = 4
    """
    Two air gaps can be defined with "bot" (center of lower core) and "top" (between backside of upper core and
    stacked core).
    """


class AirGapLegPosition(IntEnum):
    """Sets the core at which the air gap will be added. Currently only CenterLeg is supported.

    Used when adding an air gap to the model.
    """

    LeftLeg = -1
    """Air gap in left leg.
    """
    CenterLeg = 0
    """Air gap in center leg.
    """
    RightLeg = 1
    """Air gap in right leg.
    """


class StackedPosition(IntEnum):
    """For stacked cores: options to place air gaps.

    1 to place air gap in the top winding window
    2 to place air gap in the bot winding window
    """

    Top = 1
    Bot = 2


class WindingType(IntEnum):
    """Internally used in VirtualWindingWindow class."""

    TwoInterleaved = 1
    """only two winding (transformer) interleaving
    """
    Single = 2
    """only one winding (with n turns) in the virtual winding window
    """
    CenterTappedGroup = 3
    """special 3 winding topology with typical center gapped winding schemes
    """


class WindingScheme(IntEnum):
    """Used when adding a single winding to the virtual winding window.

    Only used with a rectangular solid conductor.
    """

    Full = 1
    """The whole virtual winding window is filled with one conductor.
    """
    SquareFullWidth = 2
    """Not implemented. Foils are drawn along x-axis first and then along y-axis.
    """
    FoilVertical = 3
    """Foils are very tall (y-axis) and drawn along x-axis.
    """
    FoilHorizontal = 4
    """Foils are very wide (x-axis) and drawn along y-axis.
    """


class InterleavedWindingScheme(IntEnum):
    """Used when adding an interleaved winding to the virtual winding window."""

    Bifilar = 1
    """Not implemented.
    """
    VerticalAlternating = 2
    """Not implemented. First and second winding are interleaved vertically (rows)
    """
    HorizontalAlternating = 3
    """First and second winding are interleaved horizontally (cols)
    """
    VerticalStacked = 4
    """First winding is drawn bottom to top. Second winding is drawn top to bottom.
    """


class ConductorArrangement(IntEnum):
    """Set for round conductors when having a single conductor in the virtual winding window."""

    Square = 1
    """Turns are drawn in a grid (perfectly aligned). First drawn in y-direction then x-direction.
    """
    SquareFullWidth = 2
    """Turns are drawn in a grid. First drawn in x-direction then in y-direction .
    """
    Hexagonal = 3
    """Turns are drawn more compact. The turn of the next line slides in the empty space between two
    turns of the previous line. First drawn in y-direction then x-direction.
    """

class Align(IntEnum):
    """Specifies the distribution direction for conductors when starting from the center. This can be done for having single windings in vww."""

    ToEdges = 1
    """Conductors are placed according to the specified peripheral placing strategy without adjusting for central alignment."""

    CenterOnVerticalAxis = 2
    """Conductors are centered across the middle line of the vertical axis."""

    CenterOnHorizontalAxis = 3
    """Conductors are centered across the middle line of the horizontal axis."""

class ConductorDistribution(IntEnum):
    """Defines specific strategies for placing conductors starting from the peripheral (edges) of the virtual winding window."""

    VerticalUpward_HorizontalRightward = 1
    """Places conductors vertically upwards first, then moves horizontally rightward for the next set with consistent direction."""

    VerticalUpward_HorizontalLeftward = 2
    """Places conductors vertically upwards first, then moves horizontally leftward for the next set with consistent direction."""

    VerticalDownward_HorizontalRightward = 3
    """Places conductors vertically downwards first, then moves horizontally rightward for the next set with consistent direction."""

    VerticalDownward_HorizontalLeftward = 4
    """Places conductors vertically downwards first, then moves horizontally leftward for the next set with consistent direction."""

    HorizontalRightward_VerticalUpward = 5
    """Places conductors horizontally rightward first, then moves vertically upward for the next set with consistent direction."""

    HorizontalRightward_VerticalDownward = 6
    """Places conductors horizontally rightward first, then moves vertically downward for the next set with consistent direction."""

    HorizontalLeftward_VerticalUpward = 7
    """Places conductors horizontally leftward first, then moves vertically upward for the next set with consistent direction."""

    HorizontalLeftward_VerticalDownward = 8
    """Places conductors horizontally leftward first, then moves vertically downward for the next set with consistent direction."""

class FoilVerticalDistribution(IntEnum):
    """Defines specific strategies for placing rectangular foil vertical conductors starting from the peripheral (edges) of the virtual winding window."""

    HorizontalRightward = 1
    """Moves horizontally rightward for the next set with consistent direction."""

    HorizontalLeftward = 2
    """Moves horizontally leftward for the next set with consistent direction."""

class FoilHorizontalDistribution(IntEnum):
    """Defines specific strategies for placing rectangular foil horizontal conductors starting from the peripheral (edges) of the virtual winding window."""

    VerticalUpward = 1
    """Moves vertically upward for the next set with consistent direction."""

    VerticalDownward = 2
    """Moves vertically downward for the next set with consistent direction."""

class CenterTappedInterleavingType(IntEnum):
    """Contains different interleaving types for the center tapped transformer."""

    custom = 1
    TypeA = 2
    TypeB = 3
    TypeC = 4
    TypeD = 5


class InterleavingSchemesFoilLitz(str, Enum):
    """
    Contains interleaving schemes for mixed foil/litz windings.

    ----sec---
    ooo-primary-ooo
    ooo-primary-ooo
    ---ter---
    ---sec---
    ooo-primary-ooo
    ooo-primary-ooo
    ---ter---
    """

    ter_3_4_ter_sec_4_3_sec = "ter_3_4_ter_sec_4_3_sec"
    ter_4_3_ter_sec_3_4_sec = "ter_4_3_ter_sec_3_4_sec"
    ter_3_4_sec_ter_4_3_sec = "ter_3_4_sec_ter_4_3_sec"
    ter_4_3_sec_ter_3_4_sec = "ter_4_3_sec_ter_3_4_sec"
    ter_sec_3_4_4_3_sec_ter = "ter_sec_3_4_4_3_ter_sec"
    ter_sec_4_3_3_4_sec_ter = "ter_sec_4_3_3_4_ter_sec"
    _4_3_ter_sec__sec_ter_3_4 = "4_3_ter_sec_sec_ter_3_4"
    ter_sec_5_ter_sec_5_ter_sec_5_ter_sec = "ter_sec_5_ter_sec_5_ter_sec_5_ter_sec"


class ConductorType(IntEnum):
    """Sets the type of the conductor."""

    RoundSolid = 1
    RoundLitz = 2
    RectangularSolid = 3


class WrapParaType(IntEnum):
    """
    Sets the wrap para type.

    Only necessary for a single conductor in a virtual winding window and a FoilVertical winding scheme.
    """

    FixedThickness = 1
    """The foils have a fixed thickness given when creating the conductor. The virtual winding window
    may not be fully occupied.
    """
    Interpolate = 2
    """
    The foils will have a dynamic thickness. The thickness is chosen in such way that the virtual winding window is
    fully occupied. The thickness parameter when creating the conductor is irrelevant.
    """
    CustomDimensions = 3
    """
    The foils will have a given width and thickness. The virtual winding window may not be fully occupied..
    """

class ConductorMaterial(IntEnum):
    """Sets the conductivity of the conductor."""

    Copper = 1
    Aluminium = 2


class LossApproach(IntEnum):
    """Sets the way how losses will be calculated."""

    Steinmetz = 1
    LossAngle = 2


class PermeabilityType(IntEnum):
    """Sets the way how permeability data is received."""

    FixedLossAngle = 1
    RealValue = 2
    FromData = 3


class ExcitationMeshingType(IntEnum):
    """When running an excitation it is possible to not mesh at every frequency."""

    MeshOnlyLowestFrequency = 1
    MeshOnlyHighestFrequency = 2
    MeshEachFrequency = 3


# Following Enums must always be consistent with the materialdatabase
class MaterialDataSource(str, Enum):
    """Sets the source from where data is taken."""

    Custom = "custom"
    Measurement = "measurements"
    ManufacturerDatasheet = "manufacturer_datasheet"


class MeasurementDataType(str, Enum):
    """Sets the type of measurement data."""

    ComplexPermeability = "complex_permeability"
    ComplexPermittivity = "complex_permittivity"
    Steinmetz = "Steinmetz"

class CoreMaterialType(str, Enum):
    """Sets the core material type for the whole simulation. Needs to be given to the MagneticComponent on creation."""

    Linear = "linear"
    Imported = "imported"
