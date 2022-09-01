from enum import Enum

class ComponentType(Enum):
    Inductor = "inductor"
    Transformer = "transformer"
    IntegratedTransformer = "integrated_transformer"

class AirGapMethod(Enum):
    Center = "center"
    Percent = "percent"
    Manually = "manually"

class AirGapLegPosition(Enum):
    LeftLeg = -1
    CenterLeg = 0
    RightLeg = 1

class WindingType(Enum):
    Interleaved = "interleaved"
    Single = "single"

class WindingScheme(Enum):
    Full = "full"
    Stacked = "stacked"
    SquareFullWidth = "square_full_width"
    FoilVertical = "foil_vertical"
    FoilHorizontal = "foil_horizontal"

class InterleavedWindingScheme(Enum):
    Bifilar = "bifilar"
    VerticalAlternating = "vertical_alternating"
    HorizontalAlternating = "horizontal_alternating"
    VerticalStacked = "vertical_stacked"

class ConductorArrangement(Enum):
    Square = "square"
    SquareFullWidth = "square_full_width"
    Hexagonal = "hexagonal"

class ConductorShape(Enum):
    Round = "round"
    Rectangular = "rect"

class WrapParaType(Enum):
    FixedThickness = "fixed_thickness"
    Interpolate = "interpolate"

class Conductivity(Enum):
    Copper = "copper"
    Aluminium = "aluminium"

class LossApproach(Enum):
    Steinmetz = "steinmetz"
    LossAngle = "loss_angle"

class PermeabilityType(Enum):
    FixedLossAngle = "fixed_loss_angle"
    RealValue = "real_value"
    FromData = "from_data"

class ExcitationMeshingType(Enum):
    MeshOnce = "mesh_once"
    MeshOnlyHighestFrequency = "mesh_only_highest_frequency"
    MeshEachFrequency = "mesh_each_frequency"