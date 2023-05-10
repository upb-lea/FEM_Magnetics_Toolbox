# python libraries

# 3rd party libraries
import optuna

# femmt libraries
from femmt.optimization.ito_dtos import *


class OptunaFemmtParser:
    @staticmethod
    def parse(frozen_trial: optuna.trial.FrozenTrial) -> ItoSingleResultFile:

         return ItoSingleResultFile(
            case=frozen_trial.number,
            # geometry parameters
            air_gap_top=frozen_trial.user_attrs["air_gap_top"],
            air_gap_bot=frozen_trial.user_attrs["air_gap_bot"],
            air_gap_middle=frozen_trial.user_attrs["air_gap_middle"],
            n_p_top=frozen_trial.params["n_p_top"],
            n_p_bot=frozen_trial.params["n_p_bot"],
            n_s_top=frozen_trial.params["n_s_top"],
            n_s_bot=frozen_trial.params["n_s_bot"],
            window_h_top=frozen_trial.params["window_h_top"],
            window_h_bot = frozen_trial.params["window_h_bot"],
            window_w=frozen_trial.params["window_w"],
            core_material=frozen_trial.params["material"],
            core_inner_diameter=frozen_trial.params["core_inner_diameter"],
            primary_litz_wire=frozen_trial.params["primary_litz_wire"],
            secondary_litz_wire=frozen_trial.params["secondary_litz_wire"],

            # reluctance model results
            flux_top_max=frozen_trial.user_attrs["flux_top_max"],
            flux_bot_max=frozen_trial.user_attrs["flux_bot_max"],
            flux_stray_max=frozen_trial.user_attrs["flux_stray_max"],
            flux_density_top_max=frozen_trial.user_attrs["flux_density_top_max"],
            flux_density_bot_max=frozen_trial.user_attrs["flux_density_bot_max"],
            flux_density_stray_max=frozen_trial.user_attrs["flux_density_stray_max"],
            p_hyst=frozen_trial.user_attrs["p_hyst"],
            primary_litz_wire_loss=frozen_trial.user_attrs["primary_litz_wire_loss"],
            secondary_litz_wire_loss=frozen_trial.user_attrs["secondary_litz_wire_loss"],
            core_2daxi_total_volume=frozen_trial.values[0],
            total_loss=frozen_trial.values[1],
        )
