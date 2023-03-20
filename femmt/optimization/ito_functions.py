# python libraries
import shutil
import os

# 3rd party libraries

# femmt libraries
from femmt.optimization.integrated_transformer_dtos import *
import femmt.functions_reluctance as fr
import femmt.functions as ff
import femmt as fmt

def integrated_transformer_fem_simulations_from_result_dtos(config_dto: ItoSingleInputConfig,
                                                            simulation_dto_list: List[ItoSingleResultFile],
                                                            visualize: bool = False,
                                                            ):
    ito_target_and_fixed_parameters_dto = fmt.optimization.ito_optuna.ItoOptuna.calculate_fix_parameters(config_dto)

    time_extracted, current_extracted_1_vec = fr.time_vec_current_vec_from_time_current_vec \
        (config_dto.time_current_1_vec)
    time_extracted, current_extracted_2_vec = fr.time_vec_current_vec_from_time_current_vec \
        (config_dto.time_current_2_vec)
    fundamental_frequency = int(1 / time_extracted[-1])

    phase_deg_1, phase_deg_2 = fr.phases_deg_from_time_current(time_extracted, current_extracted_1_vec, current_extracted_2_vec)
    i_peak_1, i_peak_2 = fr.max_value_from_value_vec(current_extracted_1_vec, current_extracted_2_vec)

    for dto in simulation_dto_list:
        try:
            # 1. chose simulation type
            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                        working_directory=ito_target_and_fixed_parameters_dto.fem_working_directory,
                                        silent=True)

            window_h = dto.window_h_bot + dto.window_h_top + dto.core_inner_diameter /4

            core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=dto.core_inner_diameter,
                                                            window_w=dto.window_w,
                                                            window_h=window_h)

            core = fmt.Core(core_type=fmt.CoreType.Single,
                            core_dimensions=core_dimensions,
                            material = dto.core_material,
                            temperature = config_dto.temperature,
                            frequency = fundamental_frequency,
                            permeability_datasource = fmt.MaterialDataSource.ManufacturerDatasheet,
                            permittivity_datasource = fmt.MaterialDataSource.ManufacturerDatasheet)

            geo.set_core(core)

            table_length = dto.core_inner_diameter / 2 + dto.window_w - dto.air_gap_middle
            # 2.1 set stray path parameters
            stray_path = fmt.StrayPath(start_index=0,
                                       length=table_length)
            geo.set_stray_path(stray_path)

            # Note: bot air gap needs to be subtracted, top air gap needs to be added
            air_gap_top_position_percent = \
                                                       (dto.window_h_bot + dto.core_inner_diameter / 4 + dto.air_gap_top / 2) / window_h * 100
            air_gap_bot_position_percent = (dto.window_h_bot - dto.air_gap_bot / 2) / window_h * 100

            # 3. set air gap parameters
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, dto.air_gap_bot, air_gap_bot_position_percent)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, dto.air_gap_top, air_gap_top_position_percent)
            geo.set_air_gaps(air_gaps)

            # 4. set insulations
            insulation = fmt.Insulation()
            insulation.add_core_insulations(0.001, 0.001, 0.001, 0.001)
            insulation.add_winding_insulations([0.0002, 0.0002], 0.0001)
            geo.set_insulation(insulation)

            # 5. create winding window and virtual winding windows (vww)
            # For an integrated transformer it is not necessary to set horizontal and vertical split factors
            # since this is determined by the stray_path
            winding_window = fmt.WindingWindow(core, insulation, stray_path, air_gaps)
            top, bot = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit)

            # 6. set conductor parameters
            primary_litz = ff.litz_database()[dto.primary_litz_wire]
            secondary_litz = ff.litz_database()[dto.secondary_litz_wire]

            winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
            winding1.set_litz_round_conductor(primary_litz["conductor_radii"], primary_litz["strands_numbers"],
                                              primary_litz["strand_radii"], None, fmt.ConductorArrangement.Square)

            winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
            winding2.set_litz_round_conductor(secondary_litz["conductor_radii"], secondary_litz["strands_numbers"],
                                              secondary_litz["strand_radii"], None, fmt.ConductorArrangement.Square)

            # 7. add conductor to vww and add winding window to MagneticComponent
            top.set_interleaved_winding(winding1, dto.n_p_top, winding2, dto.n_s_top,
                                        fmt.InterleavedWindingScheme.HorizontalAlternating,
                                        0.0005)
            bot.set_interleaved_winding(winding1, dto.n_p_bot, winding2, dto.n_s_bot,
                                        fmt.InterleavedWindingScheme.HorizontalAlternating,
                                        0.0005)
            geo.set_winding_window(winding_window)

            # 8. start simulation with given frequency, currents and phases
            geo.create_model(freq=fundamental_frequency, visualize_before=visualize)
            geo.single_simulation(freq=fundamental_frequency,
                                  current=[i_peak_1, i_peak_2],
                                  phi_deg=[phase_deg_1, phase_deg_2],
                                  show_results=visualize)

            source_json_file = os.path.join(ito_target_and_fixed_parameters_dto.fem_working_directory, "results", "log_electro_magnetic.json")
            destination_json_file = os.path.join(ito_target_and_fixed_parameters_dto.fem_simulation_results_directory,
                                                f'case_{dto.case}.json')

            shutil.copy(source_json_file, destination_json_file)

            del geo

        except Exception as e:
            print(f"Exception: {e}")
