import femmt as fmt

# 1. chose simulation type
geo = fmt.MagneticComponent(component_type="inductor")

# 2. set core parameters
core = fmt.core_database()["PQ 40/40"]
# geo.core.update(window_h=0.04, window_w=0.00745,
#                mu_rel=3100, phi_mu_deg=12,
#                sigma=0.6)
geo.core.update(core_w=core["core_w"], window_w=core["window_w"], window_h=core["window_h"],
                mu_rel=3100, phi_mu_deg=12,
                sigma=0.6)

# 3. set air gap parameters
geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.0005], position_tag=[0])

geo.update_conductors(n_turns=[[5]], conductor_type=["foil"], thickness=[1e-3], wrap_para=["fixed_thickness"],
                      winding=["primary"], scheme=[None],
                      core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0005],
                      conductivity_sigma=["copper"])

geo.create_model(freq=100000, visualize_before=True, do_meshing=True, save_png=False)
geo.single_simulation(freq=100000, current=[3], show_results=True)

