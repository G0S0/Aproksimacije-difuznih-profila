#Najvjv greške: oxy i deoxy sigma, zr i zv, zbrajanje oduzimanje termonva, zbunio dermis i epidermins

import numpy as np
import matplotlib.pyplot as plt
#Parameters fit for Caucassian skin:
Cm = 0.005
Bm = 0.7
Ch = 0.005
gama = 0.75

sigma_a_oxyHemoRed = 0.75
sigma_a_oxyHemoGreen = 30
sigma_a_oxyHemoBlue = 50

sigma_a_deoxyHemoRed = 2
sigma_a_deoxyHemoGreen = 30
sigma_a_deoxyHemoBlue = 20

thickness_epidermis = 0.25

n = 1.4
F_dr = (-1.440 / n**2 + 0.710 / n + 0.668 + 0.0636 * n)
A = (1+F_dr)/(1-F_dr)

def D(mu_s, mu_a):
    return 1.0 / (3.0 * (mu_s + mu_a))


#sigma_a - absorption coefficient
#sigma_s - scattering coefficient
#epi - epidermis, derm - dermis, apos - apostrof
#lamb - light wavelength, R=680 G=530 B=470
lambRed = 680
lambGreen = 530
lambBlue = 470

def sigma_a_em(lamb):
    return 6.6*pow(10, 10)*pow(lamb, -3.33)

def sigma_a_pm(lamb):
    return 2.9*pow(10, 14)*pow(lamb, -4.75)

def sigma_a_baseline(lamb):
    return 0.0244 + 8.53*np.exp(-1*(lamb-154)/66.2)

def sigma_a_oxyHemo(lamb):
    if(lamb==680):
        return sigma_a_oxyHemoRed
    elif(lamb ==530):
        return sigma_a_oxyHemoGreen
    elif(lamb == 470):
        return sigma_a_oxyHemoBlue
    else:
        return 0
    
def sigma_a_deoxyHemo(lamb):
    if(lamb==680):
        return sigma_a_deoxyHemoRed
    elif(lamb ==530):
        return sigma_a_deoxyHemoGreen
    elif(lamb == 470):
        return sigma_a_deoxyHemoBlue
    else:
        return 0

def sigma_a_epi(lamb):
    return Cm*(Bm*sigma_a_em(lamb)+(1-Bm)*sigma_a_pm(lamb)) + (1-Cm)*sigma_a_baseline(lamb)
    #return 0.1

def sigma_a_derm(lamb):
    return Ch*(gama*sigma_a_oxyHemo(lamb)+(1-gama)*sigma_a_deoxyHemo(lamb))+(1-Ch)*sigma_a_baseline(lamb)
    #return 0.2

def sigma_s_apos(lamb):
    #return 14.74*pow(lamb, -0.22)+2.2*pow(10, 11)*pow(lamb, -4)
    return 1.0

def sigma_s_apos_derm(lamb):
    #return 14.74*pow(lamb, -0.22)+2.2*pow(10, 11)*pow(lamb, -4)*1.5
    return 1.5

def sigma_t_apos_derm(lamb):
    return sigma_a_derm(lamb) + sigma_s_apos_derm(lamb)

def sigma_t_apos_epi(lamb):
    return sigma_a_epi(lamb) + sigma_s_apos(lamb)

def alpha_apos_epi(lamb):
    return sigma_s_apos(lamb)/sigma_t_apos_epi(lamb)

def alpha_apos_derm(lamb):
    return sigma_s_apos_derm(lamb)/sigma_t_apos_derm(lamb)

def sigma_tr_epi(lamb):
    return np.sqrt(3*sigma_a_epi(lamb)*sigma_t_apos_epi(lamb))

def sigma_tr_derm(lamb):
    return np.sqrt(3*sigma_a_derm(lamb)*sigma_t_apos_derm(lamb))

def z_r_epi(lamb):
    return 1/sigma_t_apos_epi(lamb) 

def z_r_derm(lamb):
    return 1/sigma_t_apos_derm(lamb)

def z_v_epi(lamb):
    return (1+4*A/3)/sigma_t_apos_epi(lamb)

def z_v_derm(lamb):
    return (1+4*A/3)/sigma_t_apos_derm(lamb)

def d_r_epi(r, z_r_epi):
    return np.sqrt(pow(r,2) + pow(z_r_epi,2))

def d_r_derm(r, z_r_derm):
    return np.sqrt(pow(r,2) + pow(z_r_derm,2))

def d_v_epi(r, z_v_epi):
    return np.sqrt(pow(r,2) + pow(z_v_epi,2))

def d_v_derm(r, z_v_derm):
    return np.sqrt(pow(r,2) + pow(z_v_derm,2))

def dermis_dipole(r, lamb, z_r_derm, z_v_derm):
    term1 = (alpha_apos_derm(lamb)*z_r_derm*(1+sigma_tr_derm(lamb)*d_r_derm(r, z_r_derm))*np.exp(-1*sigma_tr_derm(lamb)*d_r_derm(r, z_r_derm)))/(4*np.pi*pow(d_r_derm(r, z_r_derm),3))
    term2 = (alpha_apos_derm(lamb)*z_v_derm*(1+sigma_tr_derm(lamb)*d_v_derm(r, z_v_derm))*np.exp(-1*sigma_tr_derm(lamb)*d_v_derm(r, z_v_derm)))/(4*np.pi*pow(d_v_derm(r, z_v_derm),3))
    return term1 + term2

def epidermis_dipole(r, lamb, z_r_epi, z_v_epi):
    term1 = (alpha_apos_epi(lamb)*z_r_epi*(1+sigma_tr_epi(lamb)*d_r_epi(r, z_r_epi))*np.exp(-1*sigma_tr_epi(lamb)*d_r_epi(r, z_r_epi)))/(4*np.pi*pow(d_r_epi(r, z_r_epi),3))
    term2 = (alpha_apos_epi(lamb)*z_v_epi*(1+sigma_tr_epi(lamb)*d_v_epi(r, z_v_epi))*np.exp(-1*sigma_tr_epi(lamb)*d_v_epi(r, z_v_epi)))/(4*np.pi*pow(d_v_epi(r, z_v_epi),3))
    return term1 + term2

def multipole(r, lamb):
    D1 = D(sigma_s_apos(lamb), sigma_a_epi(lamb))
    D2 = D(sigma_s_apos_derm(lamb), sigma_a_derm(lamb))

    zb1 = 2 * D1 * A
    zb2 = 2 * D2 * A

    # Source depths
    z0 = 1 / (sigma_s_apos(lamb) + sigma_s_apos_derm(lamb))
    z_r_epi = z0
    z_v_epi = z_r_epi - 2 * zb1
    # Shifted to dermis layer
    z_r_derm = thickness_epidermis + z0
    z_v_derm = -z_r_derm - 2 * zb2

    print(sigma_a_derm(lamb))
    print(sigma_a_epi(lamb))
    print(sigma_s_apos(lamb))
    print(sigma_s_apos_derm(lamb))
    return epidermis_dipole(r, lamb, z_r_epi, z_v_epi) + dermis_dipole(r, lamb, z_r_derm, z_v_derm) +epidermis_dipole(r, lamb, z_r_epi, z_v_epi) + dermis_dipole(r, lamb, z_r_derm, z_v_derm) 


r_values = np.linspace(0.01, 5, 1000)
R_values = multipole(r_values, lambRed)

plt.figure(figsize=(8, 5))
plt.plot(r_values, R_values, label="Red profile", color='red')

r_values = np.linspace(0.01, 2.5, 1000)
R_values = multipole(r_values, lambGreen)

plt.plot(r_values, R_values, label="Green profile", color='green')

r_values = np.linspace(0.01, 2.5, 1000)
R_values = multipole(r_values, lambBlue)
plt.plot(r_values, R_values, label="Blue profile", color='blue')


plt.title("Layered Multipole Diffuse Reflectance Profile R(r)")
plt.xlabel("Radius r (mm)")
plt.ylabel("R(r)")
#plt.yscale("log")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig("multipole.png")
print("Plot saved as 'layered_multipole_R_r_plot.png'")