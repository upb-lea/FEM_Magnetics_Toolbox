from DESSCA import dessca_model
import numpy as np

my_dessca_instance = dessca_model(box_constraints=[[0, 1], [0, 1], [0, 1], [0, 1],
                                                       [0, 1], [0, 1], [0, 1]],
                                       state_names=["Strom", "Kerndurchmesser",
                                                    "FensterHöhe", "Fenster-Breite", "Luftspalthöhe", "Kernisolation-Links", "Leiterradius"],
                                       bandwidth=0.1,
                                       render_online= False,
                                       )

next_sample_suggest = my_dessca_instance.update_and_sample()
for i in range(10000):
    next_sample_suggest = my_dessca_instance.update_and_sample(np.transpose([next_sample_suggest]))
    print (i)