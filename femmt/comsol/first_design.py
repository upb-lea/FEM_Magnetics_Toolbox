"""
Creation of a .mph file from scratch for magnetics modelling
"""
import mph

# The start() function returns a client object,
# i.e. an instance of the Client class.
client = mph.start()

pymodel = client.create('my_first_design')
model = pymodel.java

#######################################################################################################################
# Create component 1
model.component().create("comp1", True)

# Setting dimension to 2D-axis-symmetric
model.component("comp1").geom().create("geom1", 2)
model.component("comp1").geom("geom1").axisymmetric(True)

# Creation of mesh grid
model.component("comp1").mesh().create("mesh1")

# Setting Physics: magnetic and electric fields
model.component("comp1").physics().create("mef", "ElectricInductionCurrents", "geom1")

# Setting study: frequency domain
model.study().create("std1")
model.study("std1").create("freq", "Frequency")
model.study("std1").feature("freq").setSolveFor("/physics/mef", True)

#######################################################################################################################
# Set global definition parameters
model.param().set("core_w", "0.015 [m]", "core inner diameter")

# Remove global definition parameters
# model.param().remove("core_w")

model.save('my_first_design.mph')
