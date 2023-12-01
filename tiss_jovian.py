import numpy as np
from numpy import sqrt, cos, sin, arccos, arcsin
import matplotlib.pyplot as plt
import core_tisserand as tiss

 

# Jupiter
muJupiter = 126686534.921
 
# Io
muIo = 5959.916
RIo = 421.8e3
rpIo = 1821.6+150

# Europa
muEuropa = 3202.739
REuropa = 671.1e3
rpEuropa = 1560.8+150

# Ganymede
muGanymede = 9887.834
RGanymede = 1070e3
rpGanymede = 2634.1+150

# Callisto
muCallisto = 7179.289
RCallisto = 1882e3
rpCallisto = 2410.3+150

 


# # Figura 1
 
# TP = tiss.tisserand_plot(muJupiter, RIo)
# TP.add_secondary(muIo, RIo, rpIo, 'Io', 'orange')
# TP.add_secondary(muEuropa, REuropa, rpEuropa, 'Europa', 'green')
# TP.add_secondary(muGanymede, RGanymede, rpGanymede, 'Ganymede', 'blue')
# TP.add_secondary(muCallisto, RCallisto, rpCallisto, 'Callisto', 'magenta')


# vinf_list_IO = np.linspace(1, 8, 8)/TP.vconv
# vinf_list_EU = np.linspace(0.5, 5, 11)/TP.vconv
# vinf_list_GA = np.linspace(0.5, 5, 11)/TP.vconv
# vinf_list_CA = np.linspace(0.5, 5, 11)/TP.vconv

# TP.add_box(0, 5, 0, 10)
# TP.show_body_names()
# TP.add_vinf_contour(vinf_list_IO, 0)
# TP.add_vinf_contour(vinf_list_EU, 1)
# TP.add_vinf_contour(vinf_list_GA, 2)
# TP.add_vinf_contour(vinf_list_CA, 3)

# TP.save_to_file('tisserand_jovian.png')


# Figura 2

TP = tiss.tisserand_plot(muJupiter, RIo)
TP.add_secondary(muIo, RIo, rpIo, 'Io', 'grey') #, 'orange')
TP.add_secondary(muEuropa, REuropa, rpEuropa, 'Europa', 'grey') #'green')
TP.add_secondary(muGanymede, RGanymede, rpGanymede, 'Ganymede', 'grey') #blue
TP.add_secondary(muCallisto, RCallisto, rpCallisto, 'Callisto', 'grey')
TP.add_box(0, 5, 0, 30)
TP.show_body_names()


vinf_v = np.arange(0.5, 10, 0.5)/TP.vconv
# vinf_v = np.linspace(0.5, 5, 11)/TP.vconv

ipl = 2
# TP.add_vinf_contour(np.linspace(0.5, 5, 11)/TP.vconv, ipl)
TP.add_vinf_contour(vinf_v, ipl)
list_resonance = [  
    (1, 1, 'blue'),
    (6, 7, 'magenta'),
    (2, 1, 'red'),
    (3, 1, 'green'),
    (4, 1, 'orange'),
]
TP.add_periodo_contour(list_resonance, ipl, 'lr' ) #,'rl')



TP.save_to_file('tisserand_jovian_resonance.png')


