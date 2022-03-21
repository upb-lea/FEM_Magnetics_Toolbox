from femmt import *
import json
import matplotlib.pyplot as plt

def extract_losses(losses_file):
    with open(losses_file, "r") as fd:
        losses = json.loads(fd.read())

        return losses["Losses"]["Core_Eddy_Current"], losses["Losses"]["Core_Hysteresis"]

if __name__ == "__main__":

    frequencies = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000]
    losses_hyst = []
    losses_eddy = []

    geo = MagneticComponent(component_type="inductor")

    geo.core.update(window_h=0.03, window_w=0.011,
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)

    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.0005], position_tag=[0])

    geo.update_conductors(n_turns=[[18]], conductor_type=["solid"], conductor_radii=[0.0015],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0001],
                          conductivity_sigma=["copper"])

    for frequency in frequencies:
        geo.create_model(freq=frequency, visualize_before=False, do_meshing=True, save_png=False)

        geo.single_simulation(freq=frequency, current=[10], show_results=False)

        eddy, hyst = extract_losses(geo.e_m_results_log_path)
        losses_eddy.append(eddy)
        losses_hyst.append(hyst)

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.grid()
    plt.plot(frequencies, losses_eddy, "r--", label=r"$P_eddy$")
    plt.plot(frequencies, losses_hyst, "b-", label=r"$P_hyst$")
    plt.title(fr"Hysterese und Wirbelstromverluste")
    plt.xlabel(r"Frequenz in Hz")
    plt.ylabel(r"Verluste in W")
    #plt.xscale("log")
    #plt.yscale("log")
    plt.legend(loc="upper right")
    plt.savefig("verluste.png")
    plt.show()
