import numpy as np
from numpy import sqrt, cos, sin, arccos, arcsin
import matplotlib.pyplot as plt
import core_tisserand as tiss

# Sun
muSun = 1.3271e11       # [km^3/s^2]

# Earth
muEarth = 398600.4418
REarth = 149.5e6
rpEarth = 6378.137+150

# Venus
muVenus = 324858.63
RVenus = 108.2e6
rpVenus = 6051.8+150

# Mars
muMars = 42828.314258067
RMars = 227.9e6
rpMars = 3389.5+150

# Jupiter
muJupiter = 126686534.921
RJupiter = 778.6e6
rpJupiter = 69911+150


 


        
 
TP = tiss.tisserand_plot(muSun, REarth)

TP.add_secondary(muEarth, REarth, rpEarth, 'Earth', 'green')
TP.add_secondary(muVenus, RVenus, rpVenus, 'Venus', 'red')
TP.add_secondary(muMars, RMars, rpMars, 'Mars', 'orange')   
# TP.add_secondary(muJupiter, RJupiter, rpJupiter, 'Jupiter', 'blue')

vinf_list = np.array([1., 2., 3, 4, 5, 6, 7, 8])/TP.vconv


TP.add_box(0, 2, 0, 6)
TP.show_body_names()
TP.add_vinf_contour(vinf_list, 0)
TP.add_vinf_contour(vinf_list, 1)
TP.add_vinf_contour(vinf_list, 2)
# TP.add_vinf_contour(vinf_list, 3)

TP.save_to_file('tisserand_solar.png')


