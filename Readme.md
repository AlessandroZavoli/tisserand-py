# Tissenrand Plot


This is a simple python script to plot the Tissenrand graph (RA, RP). 

## Class: TisserandPlot

A class representing a Tisserand plot.

Attributes:
- mu_p_dim: Gravitational parameter of the primary body [adim]
- rconv_dim: Conversion factor for distance [adim]
- rconv: Orbital radius of the primary body [adim]
- vconv: Orbital velocity of the primary body [adim]
- tconv: Orbital period of the primary body [adim]
- mu_p_dim: Gravitational parameter of the primary body [adim]
- secondarys: List of secondary bodies added to the Tisserand plot
- fig: Figure object for the plot
- ax: Axes object for the plot

### Methods

#### `__init__(self, mu_p_dim, rconv_dim=1)`
Initializes the Tisserand plot with the given parameters

#### `add_secondary(self, mu_s_dim, a_s_dim, rp_min_dim, name, color='black')`
Adds a secondary body to the Tisserand plot

#### `show_body_names(self)`
Displays the names of the secondary bodies on the plot

#### `add_alfa_contour(self, alfa, body_id)`
Adds an alpha contour to the Tisserand plot

#### add_periodo_contour(self, list_resonance, body_id, opz=''): 
Adds a period contour to the Tisserand plot

#### add_box(self, R_P_min, R_P_max, R_A_min, R_A_max): 
Adds a box to the Tisserand plot

#### add_vinf_contour(self, vinf, body_id): 
Adds a vinf contour to the Tisserand plot

#### save_to_file(self, filename): 
Saves the Tisserand plot to a file


<!-- ## Class: TisserandPlot

The `TisserandPlot` class is used to create a Tisserand graph, which is a plot of the Tisserand parameter against the semi-major axis of a celestial body.

### Methods

#### `__init__(self, ra, rp)`

This is the constructor method for the `TisserandPlot` class. It initializes the class with the right ascension (ra) and the radius of perigee (rp).

#### `plot(self)`

This method generates the Tisserand plot using the provided right ascension and radius of perigee.

### Usage

Here is a basic example of how to use the `TisserandPlot` class:

```python
from tisserand_plot import TisserandPlot

# Initialize the class with right ascension and radius of perigee
tp = TisserandPlot(ra=1.0, rp=2.0)

# Generate the plot
tp.plot() -->