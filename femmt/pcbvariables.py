from enum import IntEnum, Enum
from dataclasses import dataclass
"""Nachher Teilung von ENUM und dataclasses"""

""" werte siehe https://www.multi-circuit-boards.eu/leiterplatten-design-hilfe/lagenaufbau/leiterplatten-materialien.html"""
class PcbMaterial(IntEnum):
    """Defines Values of Material PCB is made of (Prepeq)"""
    fr4_standard = 1
    fr4_alternative = 2
    fr4_tracking_resistant = 3
    fr4_halogenfree = 4
    fr4_midTg = 5
    fr4_midTg_alternative = 6
    fr4_HTg = 7
    fr4_HTg_alternative = 8
    """1-8 standard board, 9-18 HF boards"""
    Rogers4350B = 9
    Rogers4003C = 10
    Megtron6 = 11
    RogersRO3003 = 12
    RogersRO3006 = 13
    RogersRO3010 = 14
    TaconicRF35 = 15
    TaconicTLX = 16
    RogersRO3001 = 17
    TaconicTLC = 18


@dataclass
class Layers:
    """defines Values of Layers in PCB"""
    LayerThickness : float
    """Thickness of Layers of Cupper"""
    LayerNumber : int
    """Number of Layers in PCB"""
    WindingPerLayer : int
    """Number of Windings in each Layer"""


