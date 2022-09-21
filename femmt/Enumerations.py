from enum import IntEnum

class WindingWindowSplit(IntEnum):
    NoSplit = 1                         # Whole window as one
    HorizontalSplit = 2                 # Horizontal line which separates top and bot
    VerticalSplit = 3                   # Vertical line which separates left and right
    HorizontalAndVerticalSplit = 4      # Split into a grid with 4 separate cells

class ComponentType(IntEnum):
    Inductor = 1
    Transformer = 2
    IntegratedTransformer = 3

class AirGapMethod(IntEnum):
    Center = 1
    Percent = 2
    Manually = 3

class AirGapLegPosition(IntEnum):
    LeftLeg = -1
    CenterLeg = 0
    RightLeg = 1

class WindingType(IntEnum):
    Interleaved = 1
    Single = 2

class WindingScheme(IntEnum):
    Full = 1
    Stacked = 2
    SquareFullWidth = 3
    FoilVertical = 4
    FoilHorizontal = 5

class InterleavedWindingScheme(IntEnum):
    Bifilar = 1
    VerticalAlternating = 2
    HorizontalAlternating = 3
    VerticalStacked = 4

class ConductorArrangement(IntEnum):
    Square = 1
    SquareFullWidth = 2
    Hexagonal = 3

class ConductorType(IntEnum):
    RoundSolid = 1
    RoundLitz = 2
    RectangularSolid = 3

class WrapParaType(IntEnum):
    FixedThickness = 1
    Interpolate = 2

class Conductivity(IntEnum):
    Copper = 1
    Aluminium = 2

class LossApproach(IntEnum):
    Steinmetz = 1
    LossAngle = 2

class PermeabilityType(IntEnum):
    FixedLossAngle = 1
    RealValue = 2
    FromData = 3

class ExcitationMeshingType(IntEnum):
    MeshOnlyLowestFrequency = 1
    MeshOnlyHighestFrequency = 2
    MeshEachFrequency = 3