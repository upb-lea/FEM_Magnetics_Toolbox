from femmt import MagneticComponent

geo = MagneticComponent(component_type="inductor")

print(geo.load_result("ME", 1, "real"))



