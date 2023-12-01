import numpy as np
from numpy import sqrt, cos, sin, arccos, arcsin


class tisserand_plot:
    def __init__(self, R_P, R_A):
        self.R_P = R_P
        self.R_A = R_A


mu_p = 1.327e11 # km^3/s^2



def direct_tisserand(vinf, alfa, R_s):

    V_s = sqrt(mu_p/R_s)

    V_sc = sqrt(vinf**2 + V_s**2 + 2*vinf*V_s*np.cos(alfa))

    V_sc_t = V_s + vinf*cos(alfa)

    a = -mu_p/(V_sc**2 - 2*mu_p/R_s)
    h = R_s*V_sc_t
    e = sqrt(1-h**2/(a*mu_p))

    R_P = a*(1-e)
    R_A = a*(1+e)

    return R_P, R_A

def inverse_tisserand(R_A, R_P, a_s):

    V_s = sqrt(mu_p/a_s)
    a_sc = 0.5*(R_P + R_A)
    e_sc = 1-R_P/a_sc
    p_sc = a_sc*(1-e_sc**2)
    cos_nu =  (p_sc-a_s)/(a_s*e_sc) 
    nu = np.arccos(cos_nu)

    V_sc_t = sqrt(mu_p/p_sc)*(1+e_sc*cos_nu)
    V_sc_r = sqrt(mu_p/p_sc)*e_sc*sin(nu)

    vinf_r = V_sc_r
    vinf_t = V_sc_t - V_s
    vinf = sqrt(vinf_r**2 + vinf_t**2)

    alfa = np.arccos(vinf_t/vinf)

    return vinf, alfa




a_s = 149.5e6 # km
rconv = a_s
vconv = sqrt(mu_p/rconv)
mu_s = 398600 # km^3/s^2



import matplotlib.pyplot as plt

# plt.figure()

# # Boundaries of the plot
# plt.plot([R_s/rconv, R_s/rconv], [0.4, 1], '-k')
# plt.plot([R_s/rconv, 3], [R_s/rconv, R_s/rconv], '-k')

# # iso-contours of vinf
# list_alfa = np.linspace(0, np.pi, 100)
# list_vinf = np.linspace(1, 8, 8)

# for vinf in list_vinf:   
#     list_R_P, list_R_A = [], []
#     for alfa in list_alfa:
#         R_P, R_A = direct_tisserand(vinf, alfa, R_s)
#         list_R_P.append(R_P)
#         list_R_A.append(R_A)
#     arr_R_P = np.array(list_R_P)/R_s
#     arr_R_A = np.array(list_R_A)/R_s
#     plt.plot(arr_R_A, arr_R_P, '-g')

# # iso-contours of alpha


# list_alfa = np.linspace(0, np.pi, 13)
# list_vinf = np.linspace(0.01, 8, 100)

# for alfa in list_alfa:
#     list_R_P, list_R_A = [], []
#     for vinf in list_vinf:   
#         R_P, R_A = direct_tisserand(vinf, alfa, R_s)
#         list_R_P.append(R_P)
#         list_R_A.append(R_A)
#     arr_R_P = np.array(list_R_P)/R_s
#     arr_R_A = np.array(list_R_A)/R_s
#     plt.plot(arr_R_A, arr_R_P, '--k')

# plt.xlabel('$R_A$ [adim]')
# plt.ylabel('$R_P$ [adim]')


# -------------------------------------------------------------------
# Tisserand plot for the Earth, showing iso-resonalances contours


def deflection_angle(vinf, R_s):

    rp_flyby = 6500 # km
    sin_deltaM = (mu_s/rp_flyby)/(vinf**2 + mu_s/rp_flyby)
    delta = 2*np.arcsin(sin_deltaM)
    return delta





def resonant(a_s, n, m):
    a_sc = a_s * (n/m)**(2/3)
    return a_sc


def iso_period_contour(a_s, n, m):
    
    R_P_min = 0.01*rconv
    
    # a_sc = resonant(a_s, n, m)
    a_sc = a_s * (n/m)**(2./3.)

    R_P = np.linspace(a_s, R_P_min, 10000)
    R_A = 2*a_sc - R_P


    return R_P, R_A
    

plt.figure()



# list_resonance = [(1,1), (2,1), (3,1)]
list_resonance = [(2,1,'m'), (3,1,'g'), (5,1,'b')]
for resonance in list_resonance:
    n, m, color = resonance
 
    # Plot the resonance contour
    arr_R_P_res, arr_R_A_res = iso_period_contour(a_s, n, m)
    # plt.plot(arr_R_A/a_s, arr_R_P/a_s, '--k')

    # per ogni coppia R_A, R_P lungo una risonanza, calcola anche l'angolo alfa corrispondente alla deflessione massima
    arr_R_P_left, arr_R_A_left = [], []
    arr_R_P_right, arr_R_A_right = [], []
    for R_A_res, R_P_res in zip(arr_R_A_res, arr_R_P_res):
        
        vinf, alfa = inverse_tisserand(R_A_res, R_P_res, a_s)
        delta_max = deflection_angle(vinf, a_s)

        alfa_right = max(0., alfa - delta_max)
        alfa_left = min(np.pi, alfa + delta_max)
        
        
        R_P_left, R_A_left = direct_tisserand(vinf, alfa_left, a_s)
        R_P_right, R_A_right = direct_tisserand(vinf, alfa_right, a_s)

        arr_R_P_left.append(R_P_left)
        arr_R_A_left.append(R_A_left)
        arr_R_P_right.append(R_P_right)
        arr_R_A_right.append(R_A_right)

    arr_R_A_res, arr_R_P_res = np.array(arr_R_A_res)/rconv, np.array(arr_R_P_res)/rconv
    arr_R_A_left, arr_R_P_left = np.array(arr_R_A_left)/rconv, np.array(arr_R_P_left)/rconv
    arr_R_A_right, arr_R_P_right = np.array(arr_R_A_right)/rconv, np.array(arr_R_P_right)/rconv

    plt.plot(arr_R_A_res, arr_R_P_res, '--',color=color)
    plt.plot(arr_R_A_left, arr_R_P_left, '-',color=color)
    plt.plot(arr_R_A_right, arr_R_P_right, '-',color=color)

# Boundaries of the plot
plt.plot([a_s/rconv, a_s/rconv], [0, 1], '-k')
plt.plot([a_s/rconv, 20], [a_s/rconv, a_s/rconv], '-k')
plt.xlabel('$R_A$ [adim]')
plt.ylabel('$R_P$ [adim]')
plt.xlim(0,20)
plt.ylim(0,1.1)
plt.show()