import pycodestyle


def test_files():

    files = ['../../femmt/constants.py',
             '../../femmt/logparser.py',
             '../../femmt/enumerations.py',
             '../../femmt/dtos.py',
             '../../femmt/hpc.py',
             '../../femmt/functions_reluctance.py',
             '../../femmt/functions.py',
             '../../femmt/functions_drawing.py',
             '../../femmt/functions_model.py',
             '../../femmt/mesh.py',
             '../../femmt/component.py',
             '../../femmt/drawing.py',
             '../../femmt/data.py',
             '../../femmt/model.py',
             '../../femmt/reluctance.py',

             # examples
             '../../femmt/examples/advanced_inductor_sweep.py',
             '../../femmt/examples/advanced_sto.py',
             '../../femmt/examples/basic_inductor.py',
             '../../femmt/examples/basic_inductor_foil_vertical.py',
             '../../femmt/examples/basic_transformer.py',
             '../../femmt/examples/basic_transformer_5_windings.py',
             '../../femmt/examples/basic_transformer_center_tapped.py',
             '../../femmt/examples/basic_transformer_integrated.py',
             '../../femmt/examples/basic_transformer_interleaved.py',
             '../../femmt/examples/basic_transformer_n_winding.py',
             '../../femmt/examples/basic_transformer_stacked.py',
             '../../femmt/examples/basic_transformer_stacked_center_tapped.py',
             '../../femmt/examples/basic_transformer_three_winding.py',
             '../../femmt/examples/excitation_sweep_example.py',
             '../../femmt/examples/hpc_examples.py',
             '../../femmt/examples/inductor_optimization.py',
             '../../femmt/examples/ito_brute_force_example.py',
             '../../femmt/examples/ito_optuna_example.py',
             '../../femmt/examples/parameter_sweep_example.py',
             '../../femmt/examples/reluctance_example.py',
             '../../femmt/examples/thermal_examples.py',
             '../../femmt/examples/Three_Winding_Transformer.py',
             '../../femmt/examples/winding_test.py',

             # tests
             'test_femmt.py',
             'lint.py',

             # optimization
             '../../femmt/optimization/sto_dtos.py',
             '../../femmt/optimization/sto.py',
             '../../femmt/optimization/optuna_femmt_parser.py',
             '../../femmt/optimization/to.py',
             '../../femmt/optimization/to_dtos.py',
             '../../femmt/optimization/ito_dtos.py',
             '../../femmt/optimization/ito_functions.py',
             '../../femmt/optimization/ito.py',
             '../../femmt/optimization/functions_optimization.py',

             # thermal files
             '../../femmt/thermal/thermal_functions.py',
             '../../femmt/thermal/thermal_simulation.py',
             '../../femmt/thermal/thermal_classes.py']

    style = pycodestyle.StyleGuide(config_file='../../tox.ini')
    result = style.check_files(files)
    print(result.total_errors)
    assert result.total_errors == 0
