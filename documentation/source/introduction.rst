FEM Magnetics Toolbox (FEMMT)
=============================

Python toolbox to generate preconfigured figures for FEM simulation
tools in power electronics.

The toolbox contains two parts, a reluctance module and a FEM module.

* The reluctance module is for pre-calculations 
* The FEM module is for detailed calculations

Installation
---------------

To run FEMMT, python (version 3.8 or above) and onelab is needed.

ONELAB installation
~~~~~~~~~~~~~~~~~~~~~~~

-  Go to https://onelab.info/
-  Download the Desktop Version for your OS (Windows, Linux or macOS)
-  Unpack the software and remember the file path. This will be needed
   later when installing FEMMT.

Install FEMMT
~~~~~~~~~~~~~~~~~

FEMMT can be installed using the python pip package manager.
Either a release version can be installed using pip or a development version by downloading this repository.

FEMMT release version (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This installs the stable release version.

::

   pip install femmt

FEMMT development version (for developers only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the latest development version with the latest features. Note:
You may need to install `git <https://git-scm.com/downloads>`__.
Also have a look at the `developers notes </developers_notes.md>`__.


::

   cd /Documents/Folder/of/Interest/FEMMT   
   git clone git@github.com:upb-lea/FEM_Magnetics_Toolbox.git
   pip install -e .

FEMMT is using the material database. To use the latest version for developing, also install the material database in developer mode.

::

   cd /Documents/Folder/of/Interest/materialdatabase   
   git clone git@github.com:upb-lea/materialdatabase.git
   pip install -e .



Examples
-----------

This toolbox is able to build a complete FEM simulation from simple
Python code. The following figure shows the Python code on the left and
the corresponding FEM simulation on the right. |image_femmt_screenshot|

Code examples can be found in this `example file </femmt/examples/basic_example.py>`__. This file is updated
regulary.

Basics
~~~~~~~~~~

The magnetic component can be an inductor, a transformer, or a
transformer with integrated stray path. The parameterization process is
divided into the following steps: 

1. Choose the simulation type, whether it's frequency domain or time domain simulation,
2. set core parameters (geometry, material), 
3. set air gap parameters (position, height),
4. set insulation distances
5. set conductor parameters (litz/solid wire),
6. start the frequency domain simulation by specifying the given frequencies, currents, phases, and the time domain simulation by specifying the given currents, and time parameters.

Please have a look at the `basic_example </femmt/examples/basic_inductor.py>`__ for frequency domain simulation, and
at the `basic_example </femmt/examples/basic_inductor_time_domain.py>`__ for time domain simulation.

The examples contain among other things: 

* Geometries: Coil, transformer, transformer with integrated stray path, 
* wire and stranded wire definition, 
* air gaps definition, 
* excitation with different frequencies, amplitudes and phases for frequency domain simulation, and with amplitudes and time steps for time domain simulation.

The simulation results can be found in `working_directory/results/result_log_electro_magnetic.json`. Working directory can either be set by the user otherwise it located at `/python-side-packages-path/femmt`.
In it you can find 

* power loss in the core: hysteresis losses and eddy current losses, 
* losses per winding and for each individual winding,
* self- and mutual inductances.

For a more detailed guide on how to create a model, please have a look :ref:`here <user_guide_model_creation>`.

Counting arrow system
~~~~~~~~~~~~~~~~~~~~~~~~~

Defined as depicted here:

|image_counting_arrow_system|

GUI (Experimental)
-------------------

There is a first preview for a GUI. Installing this is a bit cumbersome
at first, but will be simplified in the future: 

* Download the complete repository via ``Code`` -> ``Download ZIP`` and unpack it. 
* install the development version of femmt as described above 
* run python ``downloads/path-to_femmt/femmt/gui/femmt_gui.py``

Please note, the GUI is experimental.

|image_gui_definition|

FEMM Validation (for developers only)
--------------------------------------

For verification purposes a FEMM model can be created in FEMMT. To do this
FEMM needs to be installed as well as the FEMM python package.
The functionality is limited, e.g. the air gaps are limited to the 'center'-type. Other types, like 'percent' are not implemented.


.. |image_femmt_screenshot| image:: ../images/FEMMT_Screenshot.png
.. |image_counting_arrow_system| image:: ../images/counting_arrow_system.png
.. |image_gui_definition| image:: ../images/femmt_gui_definition.png
