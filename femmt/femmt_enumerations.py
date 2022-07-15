from enum import Enum

class ComponentType(Enum):
    Inductor = "inductor"
    Transformer = "transformer"
    IntegratedTransformer = "integrated_transformer"

class VirtualWindingType(Enum):
    FullWindow = "full_window"
    Split2 = "center"

class AirGapMethod(Enum):
    Center = "center"
    Percent = "percent"
    Manually = "manually"

class AirGapLegPosition(Enum):
    LeftLeg = -1
    CenterLeg = 0
    RightLeg = 1

class ConductorType(Enum):
    Stacked = "stacked"
    Full = "full"
    Foil = "foil"
    Solid = "solid"
    Litz = "litz"

class WindingType(Enum):
    Interleaved = "interleaved"
    Primary = "primary"
    Secondary = "secondary"

class WindingScheme(Enum):
    # If winding is primary or secondary
    Hexagonal = "hexa"
    Square = "square"
    Square_Full_Width = "square_full_width"
    
    # If winding is interleaved
    Horizontal = "horizontal"
    Vertical = "vertical"
    Bifilar = "bifilar"
    Blockwise = "blockwise"

class WrapParaType(Enum):
    Fixed_Thickness = "fixed_thickness"
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