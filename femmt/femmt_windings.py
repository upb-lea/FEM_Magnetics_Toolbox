from enum import Enum
from typing import List

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

class LitzImplicitType(Enum):
    default = ""
    Implicit_litz_radius = "implicit_litz_radius"
    Implicit_ff = "implicit_ff"
    Implicit_strands_number = "implicit_strands_number"

class WrapParaType(Enum):
    default = ""
    Fixed_Thickness = "fixed_thickness"
    Interpolare = "interpolate"

class Conductivity(Enum):
    default = ""
    Copper = "copper"

class Winding:
    turns_primary: int
    turns_secondary: int
    conductor_type: ConductorType = None
    ff: float = None
    strand_radius: float = None
    number_strands: int = None
    conductor_radius: float = None
    thickness: float = None
    wrap_para: WrapParaType = WrapParaType.default
    conductivity: Conductivity = Conductivity.default
    implicit_litz_type: LitzImplicitType = LitzImplicitType.default
    parallel: int = None # TODO What is this parameter?
    winding_type: WindingType
    winding_scheme: WindingScheme

    def __init__(self, turns_primary: int, turns_secondary:int, conductivity: Conductivity, winding_type: WindingType, winding_scheme: WindingScheme):
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
        self.conductivity = conductivity

    def set_stacked_conductor(self, thickness: float, wrap_para: WrapParaType):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Stacked
        self.thickness = thickness
        self.wrap_para = wrap_para

    def set_full_conductor(self, thickness: float, wrap_para: WrapParaType):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Full
        self.thickness = thickness
        self.wrap_para = wrap_para

    def set_foil_conductor(self, thickness: float, wrap_para: WrapParaType):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Foil
        self.thickness = thickness
        self.wrap_para = wrap_para

    def set_solid_conductor(self, conductor_radius: float):
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Solid
        self.conductor_radius = conductor_radius

    def set_litz_conductor(self, conductor_radius: float, number_strands: int, strand_radius: float, fill_factor: float):
        """
        Only 3 of the 4 parameters are needed. The other one needs to be none
        """
        if self.conductor_type is not None:
            raise Exception("Only one conductor can be set for each winding!")

        self.conductor_type = ConductorType.Litz
        self.conductor_radius = conductor_radius
        self.number_strands = number_strands
        self.strand_radius = strand_radius
        self.ff = fill_factor 

        if number_strands is None:
            self.implicit_litz_type = LitzImplicitType.Implicit_strands_number
        elif conductor_radius is None:
            self.implicit_litz_type = LitzImplicitType.Implicit_litz_radius
        elif fill_factor is None:
            self.implicit_litz_type = LitzImplicitType.Implicit_ff
        else:
            raise Exception("Wrong inputs for litz conductor")

class Windings:
    windings: List[Winding] = []
    
    isolation_primary2primary: float = 0
    isolation_secondary2secondary: float = 0
    isolation_primary2secondary: float = 0

    isolation_winding_to_top_core: float = 0
    isolation_winding_to_right_core: float = 0
    isolation_winding_to_bot_core: float = 0
    isolation_winding_to_left_core: float = 0

    def add_winding_isolations(self, primary2primary, secondary2secondary, primary2secondary):
        self.isolation_primary2primary = primary2primary
        self.isolation_secondary2secondary = secondary2secondary
        self.isolation_primary2secondary = primary2secondary

    def add_core_isolations(self, top_core, bot_core, left_core, right_core):
        self.isolation_winding_to_top_core = top_core
        self.isolation_winding_to_bot_core = bot_core
        self.isolation_winding_to_left_core = left_core
        self.isolation_winding_to_right_core = right_core
        
    def add_winding(self, winding: Winding):
        self.windings.append(winding)

    def add_windings(self, windings: Winding):
        self.windings.extend(windings)

    def get_all_parameters(self):
        n_turns = []
        conductor_type = []
        litz_para_type = []
        ff = []
        strands_numbers = []
        strand_radii = []
        conductor_radii = []
        winding = []
        scheme = []
        conductivity_sigma = []
        wrap_para = []
        parallel = []
        thickness = []

        for index, current_winding in enumerate(self.windings):
            if current_winding.conductor_type is None:
                raise Exception(f"Conductor type for winding {index} must be set")
            n_turns.append([current_winding.turns_primary, current_winding.turns_secondary])
            conductor_type.append(current_winding.conductor_type.value)
            litz_para_type.append(current_winding.implicit_litz_type.value)
            ff.append(current_winding.ff)
            strands_numbers.append(current_winding.number_strands)
            strand_radii.append(current_winding.strand_radius)
            conductor_radii.append(current_winding.conductor_radius)
            winding.append(current_winding.winding_type.value)
            scheme.append(current_winding.winding_scheme.value)
            conductivity_sigma.append(current_winding.conductivity.value)
            wrap_para.append(current_winding.wrap_para.value)
            parallel.append(current_winding.parallel)
            thickness.append(current_winding.thickness)

        return {
            "cond_cond_isolation": [self.isolation_primary2primary, self.isolation_primary2secondary, self.isolation_primary2secondary],
            "core_cond_isolation": [self.isolation_winding_to_top_core, self.isolation_winding_to_bot_core, 
                self.isolation_winding_to_left_core, self.isolation_winding_to_right_core],
            "n_turns": n_turns,
            "conductor_type": conductor_type,
            "litz_para_type": litz_para_type,
            "ff": ff,
            "strands_numbers": strands_numbers,
            "strand_radii": strand_radii,
            "conductor_radii": conductor_radii,
            "winding": winding,
            "scheme": scheme,
            "conductivity_sigma": conductivity_sigma,
            "wrap_para": wrap_para,
            "parallel": parallel,
            "thickness": thickness
        }