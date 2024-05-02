FEM Magnetics Toolbox (FEMMT)
=============================

This README only provides a short overview. For more information please have a look at the detailed documentation `here <https://upb-lea.github.io/FEM_Magnetics_Toolbox/main/intro.html>`__.

Python toolbox to generate preconfigured figures for FEM simulation
tools in power electronics.

The toolbox contains two parts, a reluctance module and a FEM module. 

* The reluctance module is for pre-calculations 
* The FEM module is for detailed calculations

The toolbox is accessible via python code or a graphical user interface
(GUI), which current development status is experimental. 

|image0|

Functionality examples 

* work with pre-defined standard core structures
* work with pre-defined litz wires 
* use python to perform parametersweeps, e.g. perform several automated simulations of different air gap sizes 
* read the results automated with python from the FEM simulation tool
* run a thermal simulation to see the temperatures, once the magnetoquasistatic simulation has finished

**Note: Alpha Version!** 

* GUI is experimental, 
* reluctance module is currently working for a single optimization example and not fully implemented yet.

Documentation
-------------------
Please have a look at the `documentation <https://upb-lea.github.io/FEM_Magnetics_Toolbox/intro.html>`__. You will find tutorials and a function description.

Installation
---------------

To run FEMMT python (version 3.8 or above) and onelab is needed.

ONELAB installation
~~~~~~~~~~~~~~~~~~~~~~~

-  Go to https://onelab.info/
-  Download the Desktop Version for your OS (Windows, Linux or macOS)
-  Unpack the software and remember the file path. This will be needed
   later when installing FEMMT.

Install FEMMT
~~~~~~~~~~~~~~~~~

FEMMT can be installed using the python pip package manager.
This is the stable release version (recommended). 

::

   pip install femmt

For working with the latest version, refer to the `documentation <https://upb-lea.github.io/FEM_Magnetics_Toolbox/main/intro.html>`__.

Minimal example
------------------

This toolbox is able to build a complete FEM simulation from simple
Python code. The following figure shows the Python code on the left and
the corresponding FEM simulation on the right. |image1|

To run a minimal example please have a look at the `basic_example.py </femmt/examples/basic_example.py>`__.

GUI (Experimental)
-------------------

There is a first preview for a GUI. Installing this is a bit cumbersome
at first, but will be simplified in the future: 

* Download the complete repository via ``Code`` -> ``Download ZIP`` and unpack it. 
* install the development version of femmt as described above 
* run python ``downloads/path-to_femmt/femmt/gui/femmt_gui.py``

Please note, the GUI is experimental.

|image2|

Roadmap
----------

Planned features in 2022: 

* Software stability and general improvements, 
* add more Functionality to the GUI, 

Bug Reports
--------------

Please use the issues report button within github to report bugs.

Contributing
---------------

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change. For contributing, please refer
to this `section <Contributing.rst>`__.

Changelog
------------

Find the changelog `here <CHANGELOG.md>`__

License
----------

`GPLv3 <https://choosealicense.com/licenses/gpl-3.0/>`__

History and project status
------------------------------

This project was initially written in matlab using FEMM simulation tool.
It became clear that the project was no longer a small project. The
project should be completely rewritten, because many new complex levels
have been added. To place the project in the open source world, the
programming language python is used.

.. |image0| image:: documentation/images/femmt.png
.. |image1| image:: documentation/images/FEMMT_Screenshot.png
.. |image2| image:: documentation/images/femmt_gui_definition.png
.. |image3| image:: documentation/images/counting_arrow_system.png
