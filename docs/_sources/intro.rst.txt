.. sectnum::

************
Introduction
************
==============================
 FEM Magnetics Toolbox (FEMMT)
==============================
Python toolbox to generate preconfigured figures for FEM simulation tools in power electronics.

The toolbox contains two parts, a reluctance module and a FEM module.
 * The reluctance module is for pre-calculations
 * The FEM module is for detailed calculations

The toolbox is accessible via python code or a graphical user interface (GUI), which current development status is experimental.

.. image:: ../../documentation/femmt.png
    :width: 300px
    :align: center
    :alt: femmt overview

Functionality examples
    * work with pre-defined standard core structures
    * work with pre-defined litz wires
    * use python to perform parametersweeps, e.g. perform several automated simulations of different air gap sizes
    * read the results automated with python from the FEM simulation tool

.. note::
      Development status: Alpha
          * GUI is experimental,
          * reluctance module is currently working for a single optimization example and not fully implemented yet.

.. include:: install.rst

.. include:: usage.rst

.. include:: others.rst
