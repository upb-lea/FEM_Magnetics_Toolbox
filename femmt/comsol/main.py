import mph

# The start() function returns a client object,
# i.e. an instance of the Client class.
client = mph.start()

# Creates a new model and returns it as a [`Model`](#Model) instance.
model = client.load('C:\\Users\\tpiepe\\sciebo\\01_Makoleme-DFG\\11_Paper\\PCIM25\\comsol\\2024_11_11-comsol.mph')

print(client.names())
# model.component("comp1").geom().create("geom1", 2);
# model.component("comp1").geom("geom1").axisymmetric('true');
print(model.parameters())
for (name, value) in model.parameters().items():

    description = model.description(name)

    print(f'{description:20} {name} = {value}')

# print(mph.tree(model))

# Modify a parameter
model.parameter('lam_rad', '5e-4')
print(f"lam_rad= {model.parameter('lam_rad')}")

print(f"components: {model.components()}")

model.component("comp1").func().create("int1", "Interpolation")
model.component("comp1").func("int1").setIndex("table", 0.1, 0, 0)
# model.component("comp1").func("int1").setIndex("table", 0.05, 0, 0);
# model.component("comp1").func("int1").setIndex("table", 0, 0, 0);
# model.component("comp1").func("int1").setIndex("table", 0.05, 1, 0);
# model.component("comp1").func("int1").setIndex("table", 0.1, 2, 0);
# model.component("comp1").func("int1").setIndex("table", 0.15, 3, 0);
# model.component("comp1").func("int1").setIndex("table", 0.2, 4, 0);
# model.component("comp1").func("int1").setIndex("table", 1500, 0, 1);
# model.component("comp1").func("int1").setIndex("table", 1500, 1, 1);
# model.component("comp1").func("int1").setIndex("table", 1500, 2, 1);
# model.component("comp1").func("int1").setIndex("table", 1500, 3, 1);
# model.component("comp1").func("int1").setIndex("table", 1500, 4, 1);
# model.component("comp1").func().duplicate("int2", "int1");
# model.component("comp1").func("int1").label("mu_real");
# model.component("comp1").func("int2").label("mu_imag");
# model.component("comp1").func("int2").setIndex("table", 0, 0, 1);
# model.component("comp1").func("int2").setIndex("table", 100, 0, 1);
# model.component("comp1").func("int2").setIndex("table", 200, 1, 1);
# model.component("comp1").func("int2").setIndex("table", 400, 2, 1);
# model.component("comp1").func("int2").setIndex("table", 800, 3, 1);
# model.component("comp1").func("int2").setIndex("table", 150, 1, 1);
# model.component("comp1").func("int2").setIndex("table", 50, 0, 1);
# model.component("comp1").func("int2").setIndex("table", 100, 1, 1);
# model.component("comp1").func("int2").setIndex("table", 200, 2, 1);
# model.component("comp1").func("int2").setIndex("table", 400, 3, 1);
# model.component("comp1").func("int2").setIndex("table", 800, 4, 1);
# model.component("comp1").func("int2").createPlot("pg4");

# model.mesh()
# model.solve()
# model.export('Image', 'magB.png')
# model.save('C:\\Users\\tpiepe\\sciebo\\01_Makoleme-DFG\\11_Paper\\PCIM25\\comsol\\2024_11_11-comsol_modified.mph')
