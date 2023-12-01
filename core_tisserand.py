import numpy as np
from numpy import sqrt, cos, sin, arccos, arcsin
import matplotlib.pyplot as plt


def direct_transform(vinf, alfa, a_s):  
    """ 
    R_A, R_P = direct_tisserand(vinf, alfa, a_s)
    where:
    - vinf [adim]
    - alfa \in [0, \pi]
    - a_s orbital radius of the secondary body [adim]
    
    """

    V_s = sqrt(1/a_s)

    V_sc = sqrt(vinf**2 + V_s**2 + 2*vinf*V_s*np.cos(alfa))

    V_sc_t = V_s + vinf*cos(alfa)

    a = -0.5/(V_sc**2/2 - 1/a_s)
    h = a_s*V_sc_t
    e = sqrt(1-h**2/(a))

    R_P = a*(1-e)
    R_A = a*(1+e)

    R_A = np.where(R_A > 0, R_A, np.nan)
    R_P = np.where(R_A > 0, R_P, np.nan)
    return R_A, R_P


def inverse_transform(R_A, R_P, a_s):
    """
    vinf, alfa = inverse_tisserand(R_A, R_P, a_s)
    """
    V_s = sqrt(1/a_s)
    a_sc = 0.5*(R_P + R_A)
    e_sc = 1-R_P/a_sc
    p_sc = a_sc*(1-e_sc**2)
    cos_nu =  (p_sc-a_s)/(a_s*e_sc) 
    nu = np.arccos(cos_nu)

    V_sc_t = sqrt(1/p_sc)*(1+e_sc*cos_nu)
    V_sc_r = sqrt(1/p_sc)*e_sc*sin(nu)

    vinf_r = V_sc_r
    vinf_t = V_sc_t - V_s
    vinf = sqrt(vinf_r**2 + vinf_t**2)

    alfa = np.arccos(vinf_t/vinf)

    return vinf, alfa


def deflection_angle(vinf, R_s, mu_s, rp_min):
    """
    vinif [adim]
    R_s [adim]
    rp_min [adim]
    """
    sin_deltaM = (mu_s/rp_min)/(vinf**2 + mu_s/rp_min)
    delta = 2*np.arcsin(sin_deltaM)
    return delta




class tisserand_plot:
    def __init__(self, mu_p_dim, rconv_dim=1):

        self.rconv = rconv_dim
        self.vconv = sqrt(mu_p_dim/rconv_dim)
        self.tconv = self.rconv/self.vconv

        self.mu_p_dim = mu_p_dim
        self.secondarys = []

        self.fig = plt.figure(figsize=(12, 6), dpi=300)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel('$R_A$ [adim]')
        self.ax.set_ylabel('$R_P$ [adim]')
        

   
    def add_secondary(self, mu_s_dim, a_s_dim, rp_min_dim, name, color='black'):
        """
        Add a secondary to the Tisserand plot.
        a secondary body is saved as dictionary with the following keys:
        - name: name of the secondary
        - mu_s: gravitational parameter of the secondary [adim]
        - a_s: orbital radius of the secondary [adim]
        - rp_min: minimum equatorial radius of the hyperbola [adim]
        """
        self.secondarys.append(
            {
                'name': name,
                'color': color,
                'mu_s': mu_s_dim/self.mu_p_dim,
                'a_s': a_s_dim/self.rconv,
                'rp_min': rp_min_dim/self.rconv})
           
    def show_body_names(self):        
        for body in self.secondarys:
            self.ax.text(body['a_s']-0.01, body['a_s']+0.01, body['name'], 
                         horizontalalignment='right',
                         color=body['color'], fontsize=12)
        
 
    def add_alfa_contour(self, alfa, body_id):
        pass

    def add_periodo_contour(self, list_resonance, body_id, opz=''):
    
        a_s = self.secondarys[body_id]['a_s']
        mu_s = self.secondarys[body_id]['mu_s']
        rp_min = self.secondarys[body_id]['rp_min']
        color = self.secondarys[body_id]['color']

        for res in list_resonance:
            n, m, color = res

            a_sc = a_s * (n/m)**(2./3.)

            R_P_res = np.linspace(a_s*.999, a_s*0.01, 1000)
            R_A_res = 2*a_sc - R_P_res
    
            # # R_P_res = np.where(R_A_res >=a_s and R_P_res>=a_s, R_P_res, np.nan)
            # R_P_res.all((R_A_res >=a_s & R_P_res>=a_s).any())
            # R_A_res.all((R_A_res >=a_s) and (R_P_res>=a_s))
            mask = (R_A_res >=a_s) & (R_P_res<=a_s)
            R_P_res, R_A_res = R_P_res[mask], R_A_res[mask]

            self.ax.plot(R_A_res, R_P_res, '--', label=str(n) + ':' + str(m), color=color)

            self.ax.text(R_A_res[0], R_P_res[0]+0.1, str(n) + ':' + str(m), 
                        #  horizontalalignment='left',
                         color=color, fontsize=12)

            # per ogni coppia R_A, R_P lungo una risonanza, calcola anche l'angolo alfa corrispondente alla deflessione massima           
            vinf, alfa = inverse_transform(R_A_res, R_P_res, a_s)
            delta_max = deflection_angle(vinf, a_s, mu_s, rp_min)

            if 'l' in opz:
                alfa_left = alfa + delta_max
                alfa_left = np.where(alfa_left<=np.pi, alfa_left, np.pi)
                R_A_left, R_P_left = direct_transform(vinf, alfa_left, a_s)
                self.ax.plot(R_A_left, R_P_left, '-',color=color)

            if 'r' in opz:
                # alfa_right = np.max(np.zeros_like(alfa), alfa - delta_max)
                alfa_right = alfa - delta_max
                alfa_right = np.where(alfa_right>=0, alfa_right, 0)
                R_A_right, R_P_right = direct_transform(vinf, alfa_right, a_s)
                self.ax.plot(R_A_right, R_P_right, '-',color=color)




    def add_box(self, R_P_min, R_P_max, R_A_min, R_A_max):
        """
        Add a box to the Tisserand plot.
        """
        self.ax.set_xlim(R_A_min, R_A_max)
        self.ax.set_ylim(R_P_min, R_P_max)
        
        for body in self.secondarys:
        # body = self.secondarys[0]
        # print("body['a_s'] = ", body['a_s'])
            self.ax.vlines(x=body['a_s'], ymin=0., ymax=body['a_s'], color='black', linestyle='--')
            self.ax.hlines(body['a_s'], xmin=body['a_s'], xmax=100., color='black', linestyle='--')

    
    def add_vinf_contour(self, vinf, body_id):   

        if body_id > len(self.secondarys):
            raise ValueError('The body_id is not in the list of secondarys')

        a_s = self.secondarys[body_id]['a_s']
        mu_s = self.secondarys[body_id]['mu_s']
        rp_min = self.secondarys[body_id]['rp_min']
        color = self.secondarys[body_id]['color']

        n_alfa = 100
        alfa = np.linspace(0, np.pi, n_alfa)

        for vinf_i in vinf:
            R_A, R_P = direct_transform(vinf_i, alfa, a_s)

            self.ax.plot(R_A, R_P, label='vinf = ' + str(vinf_i), color=color, alpha=0.5)
                         
            

    def save_to_file(self, filename):
        self.fig.savefig(filename)