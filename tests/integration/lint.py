import pycodestyle


def test_files():

    files = ['../../femmt/constants.py',
             '../../femmt/logparser.py',
             '../../femmt/enumerations.py',
             '../../femmt/dtos.py',
             '../../femmt/hpc.py',
             '../../femmt/functions_reluctance.py',
             '../../femmt/functions.py',
             '../../femmt/mesh.py',
             '../../femmt/component.py',

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

             # tests
             'test_femmt.py',
             'lint.py',

             # optimization
             '../../femmt/optimization/sto_dtos.py',
             '../../femmt/optimization/optuna_femmt_parser.py',
             '../../femmt/optimization/to.py',
             '../../femmt/optimization/to_dtos.py',
             '../../femmt/optimization/ito_dtos.py',
             '../../femmt/optimization/ito_functions.py',

             # thermal files
             '../../femmt/thermal/thermal_functions.py']

    style = pycodestyle.StyleGuide(config_file='../../tox.ini')
    result = style.check_files(files)
    print(result.total_errors)
    assert result.total_errors == 0
