************
Installation
************

=======================
ONELAB installation
=======================
* Go to https://onelab.info/
* Download the Desktop Version for your OS (Windows, Linux or macOS)
* Unpack the software and remember the file path. This will be needed later when installing FEMMT.

=========================================
FEMM installation [for Windows User only]
=========================================
* Go to https://www.femm.info/wiki/Download
* Install FEMM as described
* FEMM can be used as an alternative FEM solver for 2D simulations

=============
Install FEMMT
=============
Chose to install the development version of FEMMT or the release version.

FEMMT release version (recommended)
***********************************
This is the stable release version.
.. code-block::

	pip install femmt

FEMMT development version (for developers only)
***********************************************
This is the latest development version with the latest features.
Note: You may need to install [git](https://git-scm.com/downloads).

.. code-block::

    cd /Documents/Folder/of/Interest   
    git clone git@github.com:upb-lea/FEM_Magnetics_Toolbox.git
    pip install -e .

=============================
Minimal example and first run
=============================
Run the example from here: [basic_example.py](/femmt/examples/basic_example.py).
FEMMT will ask you for the installation path of ONELAB during first use.

