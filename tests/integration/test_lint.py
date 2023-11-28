import pycodestyle


def test_files():

    files = ['../../femmt/examples/advanced_inductor_sweep.py',
             '../../femmt/constants.py',
             '../../femmt/logparser.py',
             '../../femmt/enumerations.py',
             '../../femmt/dtos.py',
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
             '../../tests/integration/test_femmt.py',
             'test_lint.py'
             ]
    style = pycodestyle.StyleGuide(config_file='../../tox.ini')
    result = style.check_files(files)
    print(result.total_errors)
    assert result.total_errors == 0
