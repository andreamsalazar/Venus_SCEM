from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from scipy.integrate import quad
from scipy.special import sph_harm
from scipy.special import gamma
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import time 
import logging 
from scipy.interpolate import interp1d

import sys

# Get the parameter value passed as a command-line argument

N_iterations = int(sys.argv[1]) # how many iterations for this batch
rheo_model = str(sys.argv[2]) # rheological model
tau_e = float(sys.argv[3]) # elastic relaxation number for Andrade in seconds
alph = float(sys.argv[4]) # Andrade scaling law
ps_min = float(sys.argv[5]) # minimum surface pressure (bar)
ps_max = float(sys.argv[6]) # maximum surface pressure (bar)
delta_t_s_min = float(sys.argv[7]) # minimum time in steam (Gyr)
delta_t_s_max = float(sys.argv[8]) # maximum time in steam (Gyr)
day_length_min = float(sys.argv[9])  # minimum initial rotation period (s)
day_length_max = float(sys.argv[10]) # maximum initial rotation period (s)
obl_o_min = np.radians(float(sys.argv[11])) # minimum initial obliquity (degrees)
obl_o_max = np.radians(float(sys.argv[12])) # maximum initial obliquity (degrees)
Q_min = float(sys.argv[13]) # minimum tidal quality factor, used for Constant-Q or Ocean_Andrade
Q_max = float(sys.argv[14]) # maximum tidal quality factor, used for Constant-Q or Ocean_Andrade
tau_sw_steam_min = float(sys.argv[15]) # minimum tau_sw in steam atmosphere
tau_sw_steam_max = float(sys.argv[16]) # maximum tau_sw in steam atmosphere
kappa_sw = float(sys.argv[17]) # atmospheric shortwave absorption coefficient
kappa_lw = float(sys.argv[18]) # atmospheric shortwave absorption coefficient
S_RG = float(sys.argv[19]) # runaway greenhouse limit [Wm-2]
visc = float(sys.argv[20]) # effective viscosity of core [m2s-1]
ramp_min = float(sys.argv[21]) # minimum ramping coefficient viscosity of mantle [m2s-1]
ramp_max = float(sys.argv[22]) # maximum ramping coefficient viscosity of mantle [m2s-1]

## interpolate albedo from Yang 2014

n_Yang = np.pi*2/(365*60*60*24)
rotation_periods = np.array([1,16,64,128,256,365])
sigs_Yang = 2*np.pi/(rotation_periods[::-1]*60*60*24)-n_Yang
all_solar_constants = np.array([1379.4943747591249, 1414.104535777566, 1466.9267604350466, 1509.9349498795164, 1558.9332174259953, 1611.1126981332104, 1806.115794878488, 2004.5074574305263, 2106.921933698868, 2393.6433608367, 2695.81674655818, 2993.4288665912286, 3282.220865044703, 3579.9226244639185])
from scipy.interpolate import RegularGridInterpolator

albedo_matrix = np.load('../input_data/albedo_Yang_matrix.npy')
interp_albedo = RegularGridInterpolator((sigs_Yang, all_solar_constants), albedo_matrix,bounds_error=False, fill_value=0) 


log_filename = f"Venus_evolution.log"
logging.basicConfig(
    filename=log_filename,     # Log file name
    level=logging.INFO,        # Log level (INFO, DEBUG, ERROR, etc.)
    format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
)    

G = 6.67e-11 # gravitational constant
AU = 1.49e11 # [m], astronomical unit
a = 0.723*AU # [m], Venus' semi-major axis
R_V = 6.0518e6 # [m], radius of Venus
M_V = 4.8685e24 # [kg], mass of Venus
M_sun = 1.989e30 # [kg], mass of Sun
Ifactor = 0.336 # moment of inertia factor for Venus
C = Ifactor*M_V*R_V**2 # moment of inertia Venus on spin axis
deltEd = 1.3e-5 # departure from dynamical ellipticity, CL2003 pg. 3

## for gravitational tide ##

k2 = 0.295 # 2nd order love number 
n = 2*np.pi/(224.701*60*60*24)# mean motion, 2pi/T Venus

sec_per_year = 60*60*24*365

Kg = (G*M_sun**2*R_V**5)/(C*a**6) # gravitational tide constant


## for CMF ##

Rc = 3.2e6 # m, suggested core radius of Venus, Correia 2003
Cc = 0.084*C # moment of inertia of core
deltE_c = 4*deltEd # core dynamical ellipticity anomaly 

gamma_el = 0.75 # C2003, I, pg. 8, eq. 50
ed_ec = 1/4 # non-hydrostatic, 4/3 for hydrostatic 


## for thermal tide

g = 8.87 # m/s2
sigma_SB = 5.67e-8 # Stefan-Boltzmann constant [W m-2 K-4]


## constants used in Salazar & Wordsworth, 2024
Cs_land = 1000000 # heat capacity of surface, [J kg-1]
cp = 1000 # specific heat of air, [J kg-1 K-1]
rho_mean = 5204 # density of planet, [kg m-3]
D = 1.66 # diffusivity constant 
R = 188 # specific gas constant, [J kg-1 K-1]
Cd = 0.0034 # surface drag coefficient
chi = 0.17 # reference thermal coupling factor
bar = 1e5

## lumiosity vs time
L_sun = 3.828e26 # W, luminosity of Sun
to = 4.6e9 # yr, age of solar system

Ka = (3*M_sun*R_V**3)/(5*C*rho_mean*a**3) # atmospheric thermal tide torque constant 

# Modern Venus atmosphere parameters

Omega_Venus = -2*np.pi/(244*60*60*24) # [rad/s], rotation rate
alpha_Venus = 0.77 # albeo
ps_Venus = 92e5 # [Pa] surface pressure
S_Venus = 2624.3 # [Wm-2], solar constant

# for Andrade rheology
ke = 0.25 # elastic love number
kf = 0.928 # fluid love number
tau_tot = kf*tau_e/ke

## load ocean tidal dissipation from Green 2019:

Q_data = np.load('../input_data/Q_withIT_330.npy')
sig_data = np.load('../input_data/sig_330_Green.npy')
f_interp_Qocean = interp1d(sig_data, Q_data, kind = 'linear',fill_value="extrapolate")

## albedo parameterization 
def alpha_fit(S, Omega):
    sigma = Omega-n
    sigma_crit = 0.744e-6* (S/1366)**(3/4)
    alpha_slow = 0.79*(1-np.exp(-(0.86*S/1366)))+0.08 
    alpha_rapid = 0.34
    alpha = alpha_rapid + (alpha_slow-alpha_rapid)/np.sqrt(1+(sigma/sigma_crit)**2) # based on dT 
    return(alpha)


def Tbar(S, alpha, tau_sw, tau_lw):
    k = tau_sw/tau_lw 
    F_bar = S*(1-alpha)*np.exp(-tau_sw)/np.pi # average incident stellar radiation at surface (Wm-2)
    SLW = S*(1-alpha)/8*(1+D/k - (1+D/k)*np.exp(-k*tau_lw)) # surface downwelling longwave (Wm-2), Equation (15)
    T_bar = np.power((F_bar+SLW)/sigma_SB, 1/4) # average surface temperature (K)
    return(T_bar)

def wind_speed(S, alpha, tau_sw, tau_lw, ps): # Us [m/s], Equation (23)
    T_bar = Tbar(S, alpha, tau_sw, tau_lw)
    T_eq = (S*(1-alpha)/(4*sigma_SB))**(1/4) # equilibrium surface temperature (K)
    Us = np.power(R/Cd * np.maximum((T_bar - T_eq),0) * (S/2)*(1-alpha)*np.exp(-tau_sw)*(1-np.exp(-tau_lw))/ps,1/3) 
    return(Us)

def qo_wo(S, alpha, tau_sw, tau_lw, ps, outgassing=False): # output q_o (defined in Equation 29) and w_o
    k = tau_sw/tau_lw
    F_bar = S*(1-alpha)*np.exp(-tau_sw)/np.pi # average incident stellar radiation at surface (Wm-2)
    SLW = S*(1-alpha)/8*(1+D/k - (1+D/k)*np.exp(-k*tau_lw)) # surface downwelling longwave (Wm-2), Equation (15)
    T_bar = np.power((F_bar+SLW)/sigma_SB, 1/4) # average surface temperature (K)
    delta_F = 1/6 * np.sqrt(15/(2*np.pi)) # sph. harm decomp coefficient of insolation field
    T_eq = (S*(1-alpha)/(4*sigma_SB))**(1/4) # equilibrium surface temperature (K)
    Us = wind_speed(S, alpha, tau_sw, tau_lw, ps)
    circ_strength = Us/Uso # circulation strength, Equation (22)
    delt_p = chi*ps*circ_strength # Equation (23)
    Cs = cp*delt_p/g + Cs_land # heat capacity of surface, including atmosphere
    if outgassing:
        w_o = 3.77e-7 # from Leconte 2015
    else:
        w_o = 4*sigma_SB*T_bar**3/Cs # thermal equilibrium frequency
    qo = -(1/4)*delta_F*S*(1-alpha)*np.exp(-tau_sw)*delt_p/(F_bar + SLW)*np.sqrt(10/(3*np.pi))
    return(qo, w_o)

def torque_analytic_forcingfreq(S, alpha, tau_sw, tau_lw, ps, Omega, n, m, l, outgassing=False): # q_tilde [Pa], Equation (29)
    qo, w_o = qo_wo(S, alpha, tau_sw, tau_lw, ps, outgassing)
    sigma = (m*Omega-l*n)
    torque = -qo*sigma/w_o/(1+(sigma/w_o)**2) # Equation (29)
    return(torque)


def thermal_tide_torque(M, R, a,S, alpha, tau_sw, tau_lw, ps, Omega_list, n, m, l, outgassing=False): # T_a, Equation (32)
    Ka = -3*M*R**3/((5*rho_mean*a**3))
    q_tilde = torque_analytic_forcingfreq(S, alpha, tau_sw, tau_lw, ps, Omega_list, n, m,l,outgassing)
    return(-Ka*q_tilde)

def gravitational_tide_torque(M, R, a,Q,Qn, k2, Omega_list, n, rheo_model='Constant-Q'): # T_g, Equation (4)
    Kg = -G*M**2*R**5/a**6
    if rheo_model == 'Constant-Q':
        bg = k2/Q * np.sign(Omega_list-n) 
    elif rheo_model == 'CL2003':
        bg = k2/Q * np.sign(Omega_list-n)  * (1-(1-Q/Qn)**(np.absolute(2*(Omega_list-n))/n))
    elif rheo_model == 'Andrade':
        bg = Andrade(2*(Omega_list-n))
    else:
        print('Error: please enter valid rheology model')
        bg = np.nan
    return(Kg*bg)

def Andrade(sig):
    B_sig = 1 + np.absolute(sig*tau_tot)**(1-alph) * (tau_e/tau_tot)**(1-alph)*np.sin(alph*np.pi/2)*gamma(1+alph)
    A_sig = (sig*tau_tot)*(1 + np.maximum(np.absolute(sig*tau_tot),1e-5)**(-alph) * (tau_e/tau_tot)**(1-alph)*np.cos(alph*np.pi/2)*gamma(1+alph) )
    k2Q = (kf-ke) * B_sig*sig*tau_tot/(A_sig**2 + B_sig**2) 
    return(k2Q)


Uso = wind_speed(1137, 0.2, 0.00001, 1, 1*bar) # U_so [m/s]    

## tune modern Venus atmosphere tide so that modern Venus rotation rate is equilibrium regardless of rheology
# def modern_Venus_tuning(x, rheo_model, Q, Qn):
#     tau_lw = x[0]
#     tau_sw = x[1]
#     k = tau_sw/tau_lw
#     S = S_Venus
#     alpha = 0.77
#     D = 1.66
#     F_bar = S*(1-alpha)*np.exp(-tau_sw)/np.pi # average incident stellar radiation at surface (Wm-2)
#     SLW = S*(1-alpha)/8*(1+D/k - (1+D/k)*np.exp(-k*tau_lw)) # surface downwelling longwave (Wm-2), Equation (15)
#     T_bar = np.power((F_bar+SLW)/sigma_SB, 1/4) # average surface temperature (K)
#     Tbar = T_bar-737 # choose solution with T_bar=737 K
#     rot_eq = gravitational_tide_torque(M_sun, R_V, a,Q,Qn, k2, Omega_Venus, n, rheo_model) + thermal_tide_torque(M_sun, R_V, a,S, alpha_Venus, tau_sw, tau_lw, ps_Venus, Omega_Venus, n,2,2)
#     return(Tbar, rot_eq)

# if rheo_model in ['Andrade', 'Ocean_Andrade']:
#     tau_lw_Venus, tau_sw_Venus = fsolve(modern_Venus_tuning, [200,0.5], args = ('Andrade', np.nan, np.nan))
#     logging.info("[tau_lw, tau_sw] for " + rheo_model+ ": {}, {}".format(tau_lw_Venus, tau_sw_Venus))

def modern_Venus_tuning(x, rheo_model, Q, Qn):
    A = x
    rot_eq = gravitational_tide_torque(M_sun, R_V, a,Q,Qn, k2, Omega_Venus, n, rheo_model) + thermal_tide_torque(M_sun, R_V, a,S_Venus, alpha_Venus, 0.54+A*(1+np.tanh((ps_Venus/bar - 50)/15))*0.5, 1 * (ps_Venus/bar)**(0.8), ps_Venus, Omega_Venus, n,2,2,True)
    return(rot_eq)
A_sw = float(fsolve(modern_Venus_tuning, [0.5], args = ('Andrade', np.nan, np.nan)))
tau_lw_Venus =  1 * (ps_Venus/bar)**(0.8)
tau_sw_Venus = 0.54+A_sw*(1+np.tanh((ps_Venus/bar - 50)/15))*0.5
logging.info("[tau_lw, tau_sw] for " + rheo_model+ ": {}, {}".format(tau_lw_Venus, tau_sw_Venus))



def spin_evolution_steam(t,x, rheo_model, tau_sw_steam):
    L = x[0] # L = C*omega
    obl = x[1] # obliquity [rad]
    omega = L/C 
    X = L*np.cos(obl)
    Lum = L_sun*(1+(2/5)*(1-t/to))**(-1) # Gough, 1981, solar evolution in time
    S = Lum/(4*np.pi*a**2) # solar constant, Wm-2

    pad_omega = np.sign(omega)*max(abs(omega), 1e-7) # for stable passage through omega = 0
        
    #########################################
    ########## Gravitational Tide ###########
    #########################################

    if rheo_model == 'Ocean_Andrade':
        def b_g_sig(sigma):
            k2Q_ocean = 0.2/f_interp_Qocean(sigma) * np.tanh((sigma)/(0.1*n)) * 0.5*(1+np.tanh((hab-0.15)/0.05)) # smooth change of sign and hab
            k2_Andrade = Andrade(sigma)
            return(k2Q_ocean + k2_Andrade)

    elif rheo_model == 'Andrade':
        def b_g_sig(sigma):
            return(Andrade(sigma))
    elif rheo_model == 'Constant-Q':
        def b_g_sig(sigma):
            if sigma >= 0:
                sign = 1
            else: 
                sign = -1
            return(k2*sign/Q * (1-(1-Q/Qn)**(abs(sigma)/n)))

    
        
    domegadt_grav = -Kg*(b_g_sig(omega)*(3/4 * X**2/L**2 * (1-X**2/L**2)) + b_g_sig(omega-2*n)*(3/16 * (1+X/L)**2 * (1-X**2/L**2))+ \
                b_g_sig(omega+2*n)*(3/16 * (1-X/L)**2 * (1-X**2/L**2) ) + b_g_sig(2*omega)*(3/8 * (1-X**2/L**2)**2) + \
                b_g_sig(2*omega-2*n)*(3/32 * (1+X/L)**4) + b_g_sig(2*omega+2*n)*(3/32 * (1-X/L)**4))

    
    
    dobldt_grav = -Kg  * np.sin(obl)/pad_omega * (b_g_sig(2*n)*(9/16 * np.sin(obl)**2) + \
                   b_g_sig(omega)*3/4 * np.cos(obl)**3 - b_g_sig(omega-2*n)*3/16*(1+np.cos(obl))**2 * (2-np.cos(obl)) +\
                   b_g_sig(omega+2*n)*3/16*(1-np.cos(obl))**2 * (2+np.cos(obl)) + b_g_sig(2*omega)*3/8*np.sin(obl)**2*np.cos(obl) +\
                   -b_g_sig(2*omega-2*n)*3/32*(1+np.cos(obl))**3 +  b_g_sig(2*omega+2*n)*3/32*(1-np.cos(obl))**3)
    
        
    ###########################################
    ########## Core Mantle Friction ###########
    ###########################################

    Ed = kf*R_V**5/(3*G*C) * pad_omega**2 + deltEd # elipticity factor 
    Ec = Ed/ed_ec # core elipticity factor
    alpha = 3*n**2/(2*pad_omega) * Ed # precession constant (rad/s)

    E = visc/(abs(pad_omega)*Rc**2) # Ekman number 
    
    kappa = 2.62*Cc*abs(omega)*np.sqrt(E) # viscous friction, CL003 Eq. 53
    
    Kf = kappa/(gamma_el*C) * (n/pad_omega)**4 *(3/2 * ed_ec)**2 # CL2003, Eq. 72, CMF torque constant

    domegadt_CMF = -omega*Kf*np.cos(obl)**2*np.sin(obl)**2 # CL2003, Eq. 71 
    dobldt_CMF = -Kf*np.cos(obl)**3*np.sin(obl) # CL2003, Eq. 74

    domegadt = domegadt_grav  + domegadt_CMF 
    dobldt = dobldt_grav   + dobldt_CMF 
    return(C*domegadt*sec_per_year, dobldt*sec_per_year)

def spin_evolution_hab(t,x, rheo_model,tau_lw_habitable, tau_sw_habitable, p_s_habitable, delta_t_s, Q, Qn):
    L = x[0] # L = C*omega
    obl = x[1] # obliquity [rad]
    omega = L/C 
    X = L*np.cos(obl)
    Lum = L_sun*(1+(2/5)*(1-t/to))**(-1) # Gough, 1981, solar evolution in time
    S = Lum/(4*np.pi*a**2) # solar constant, Wm-2

    pad_omega = np.sign(omega)*np.maximum(np.abs(omega), 1e-7) # for stable passage through omega = 0
    
    ## check if habitability is possible 
    #albedo = alpha_fit(S, pad_omega*np.sign(np.cos(obl))) # get albedo from cloud-rotation feedback
    sig_albedo = np.abs(pad_omega*np.sign(np.cos(obl)) - n) # |sigma|, assume symmetric in LOD
    albedo = float(interp_albedo([sig_albedo, S])) # interpolate Yang 2014 50-m data
    S_avg = S*(1-albedo)/4 # ASR
    ps = p_s_habitable
    tau_lw = tau_lw_habitable
    tau_sw = tau_sw_habitable
    albedo_in = albedo*1

    #########################################
    ########## Gravitational Tide ###########
    #########################################

    if rheo_model == 'Ocean_Andrade':
        def b_g_sig(sigma):
            k2Q_ocean = 0.2/f_interp_Qocean(sigma) * np.tanh((sigma)/(0.1*n)) # smooth change of sign and hab
            k2_Andrade = Andrade(sigma)
            return(k2Q_ocean + k2_Andrade)

    elif rheo_model == 'Andrade':
        def b_g_sig(sigma):
            return(Andrade(sigma))
    elif rheo_model == 'Constant-Q':
        def b_g_sig(sigma):
            if sigma >= 0:
                sign = 1
            else: 
                sign = -1
            return(k2*sign/Q * (1-(1-Q/Qn)**(abs(sigma)/n)))

    
        
    domegadt_grav = -Kg*(b_g_sig(omega)*(3/4 * X**2/L**2 * (1-X**2/L**2)) + b_g_sig(omega-2*n)*(3/16 * (1+X/L)**2 * (1-X**2/L**2))+ \
                b_g_sig(omega+2*n)*(3/16 * (1-X/L)**2 * (1-X**2/L**2) ) + b_g_sig(2*omega)*(3/8 * (1-X**2/L**2)**2) + \
                b_g_sig(2*omega-2*n)*(3/32 * (1+X/L)**4) + b_g_sig(2*omega+2*n)*(3/32 * (1-X/L)**4))

    
    
    dobldt_grav = -Kg  * np.sin(obl)/pad_omega * (b_g_sig(2*n)*(9/16 * np.sin(obl)**2) + \
                   b_g_sig(omega)*3/4 * np.cos(obl)**3 - b_g_sig(omega-2*n)*3/16*(1+np.cos(obl))**2 * (2-np.cos(obl)) +\
                   b_g_sig(omega+2*n)*3/16*(1-np.cos(obl))**2 * (2+np.cos(obl)) + b_g_sig(2*omega)*3/8*np.sin(obl)**2*np.cos(obl) +\
                   -b_g_sig(2*omega-2*n)*3/32*(1+np.cos(obl))**3 +  b_g_sig(2*omega+2*n)*3/32*(1-np.cos(obl))**3)
    
        
    ###########################################
    ########## Core Mantle Friction ###########
    ###########################################

    Ed = kf*R_V**5/(3*G*C) * pad_omega**2 + deltEd # elipticity factor 
    Ec = Ed/ed_ec # core elipticity factor
    alpha = 3*n**2/(2*pad_omega) * Ed # precession constant (rad/s)

    E = visc/(abs(pad_omega)*Rc**2) # Ekman number 
    
    kappa = 2.62*Cc*abs(omega)*np.sqrt(E) # viscous friction, CL003 Eq. 53
    
    Kf = kappa/(gamma_el*C) * (n/pad_omega)**4 *(3/2 * ed_ec)**2 # CL2003, Eq. 72, CMF torque constant

    domegadt_CMF = -omega*Kf*np.cos(obl)**2*np.sin(obl)**2 # CL2003, Eq. 71 
    dobldt_CMF = -Kf*np.cos(obl)**3*np.sin(obl) # CL2003, Eq. 74
    
    
###########################
###### Thermal Tide #######
###########################

    ## ExoTides (Salazar & Wordsworth, 2024)

    
    def b_a_sig(omega, m, l):
        q_tilde = torque_analytic_forcingfreq(S, albedo_in, tau_sw, tau_lw, ps, omega, n, m, l)
        return(-q_tilde)
    


    
    domegadt_atm = -Ka*(b_a_sig(omega,1,0)*(3/4 * X**2/L**2 * (1-X**2/L**2)) + b_a_sig(omega,1,2)*(3/16 * (1+X/L)**2 * (1-X**2/L**2))+ \
                b_a_sig(omega,1,-2)*(3/16 * (1-X/L)**2 * (1-X**2/L**2) ) + b_a_sig(omega,2,0)*(3/8 * (1-X**2/L**2)**2) + \
                b_a_sig(omega,2,2)*(3/32 * (1+X/L)**4) + b_a_sig(omega,2,-2)*(3/32 * (1-X/L)**4))
    

    
    dobldt_atm = -Ka  * np.sin(obl)/pad_omega * (b_a_sig(omega, 0, -2)*(9/16 * np.sin(obl)**2) + \
                    b_a_sig(omega,1,0)*3/4 * np.cos(obl)**3 - b_a_sig(omega,1,2)*3/16*(1+np.cos(obl))**2 * (2-np.cos(obl)) +\
                    b_a_sig(omega,1,-2)*3/16*(1-np.cos(obl))**2 * (2+np.cos(obl)) + b_a_sig(omega,2,0)*3/8*np.sin(obl)**2*np.cos(obl) +\
                    -b_a_sig(omega,2,2)*3/32*(1+np.cos(obl))**3 +  b_a_sig(omega,2,-2)*3/32*(1-np.cos(obl))**3)

    ## sum torques together ## 

    domegadt = domegadt_grav  + domegadt_CMF + domegadt_atm
    dobldt = dobldt_grav   + dobldt_CMF + dobldt_atm
    return(C*domegadt*sec_per_year, dobldt*sec_per_year)

def spinout(t, x, rheo_model,tau_lw_habitable, tau_sw_habitable, p_s_habitable, delta_t_s, Q, Qn):
        L, obl  = x
        omega = (L/C)
        pad_omega = np.sign(omega)*max(abs(omega), 1e-7)
        Lum = L_sun*(1+(2/5)*(1-t/to))**(-1)
        S = Lum/(4*np.pi*a**2)
        sig_albedo = abs(pad_omega*np.sign(np.cos(obl)) - n) # |sigma|, assume symmetric in LOD
        albedo = float(interp_albedo([sig_albedo, S])) # interpolate Yang 2014 50-m data
        #albedo = alpha_fit(S, pad_omega*np.sign(np.cos(obl)))
        S_avg = S*(1-albedo)/4
        return(S_avg - S_RG)
spinout.terminal = True   # stop integration at first crossing
spinout.direction = 1.0   # prefer upward crossing 


def spin_evolution_Venus(t,x, rheo_model,p_s_habitable, t_end_hab, Q, Qn,Q_ocean, ramp):
    L = x[0] # L = C*omega
    obl = x[1] # obliquity [rad]
    omega = L/C 
    X = L*np.cos(obl)
    Lum = L_sun*(1+(2/5)*(1-t/to))**(-1) # Gough, 1981, solar evolution in time
    S = Lum/(4*np.pi*a**2) # solar constant, Wm-2

    pad_omega = np.sign(omega)*np.maximum(np.abs(omega), 1e-7) # for stable passage through omega = 0

    ps = (ps_Venus-p_s_habitable) * ((t-t_end_hab)/(to-t_end_hab))**ramp + p_s_habitable # outgassing rate set with ramp 
    tau_lw = 1 * (ps/bar)**(0.8)
    tau_sw = 0.54+A_sw*(1+np.tanh((ps/bar - 50)/15))*0.5
    albedo_in = alpha_Venus
            
    #########################################
    ########## Gravitational Tide ###########
    #########################################

    if rheo_model in ['Andrade','Ocean_Andrade']:
        def b_g_sig(sigma):
            return(Andrade(sigma))
    elif rheo_model == 'Constant-Q':
        def b_g_sig(sigma):
            if sigma >= 0:
                sign = 1
            else: 
                sign = -1
            return(k2*sign/Q * (1-(1-Q/Qn)**(abs(sigma)/n)))

    
        
    domegadt_grav = -Kg*(b_g_sig(omega)*(3/4 * X**2/L**2 * (1-X**2/L**2)) + b_g_sig(omega-2*n)*(3/16 * (1+X/L)**2 * (1-X**2/L**2))+ \
                b_g_sig(omega+2*n)*(3/16 * (1-X/L)**2 * (1-X**2/L**2) ) + b_g_sig(2*omega)*(3/8 * (1-X**2/L**2)**2) + \
                b_g_sig(2*omega-2*n)*(3/32 * (1+X/L)**4) + b_g_sig(2*omega+2*n)*(3/32 * (1-X/L)**4))

    
    
    dobldt_grav = -Kg  * np.sin(obl)/pad_omega * (b_g_sig(2*n)*(9/16 * np.sin(obl)**2) + \
                   b_g_sig(omega)*3/4 * np.cos(obl)**3 - b_g_sig(omega-2*n)*3/16*(1+np.cos(obl))**2 * (2-np.cos(obl)) +\
                   b_g_sig(omega+2*n)*3/16*(1-np.cos(obl))**2 * (2+np.cos(obl)) + b_g_sig(2*omega)*3/8*np.sin(obl)**2*np.cos(obl) +\
                   -b_g_sig(2*omega-2*n)*3/32*(1+np.cos(obl))**3 +  b_g_sig(2*omega+2*n)*3/32*(1-np.cos(obl))**3)
    
        
    ###########################################
    ########## Core Mantle Friction ###########
    ###########################################

    Ed = kf*R_V**5/(3*G*C) * pad_omega**2 + deltEd # elipticity factor 
    Ec = Ed/ed_ec # core elipticity factor
    alpha = 3*n**2/(2*pad_omega) * Ed # precession constant (rad/s)

    E = visc/(abs(pad_omega)*Rc**2) # Ekman number 
    
    kappa = 2.62*Cc*abs(omega)*np.sqrt(E) # viscous friction, CL003 Eq. 53
    
    Kf = kappa/(gamma_el*C) * (n/pad_omega)**4 *(3/2 * ed_ec)**2 # CL2003, Eq. 72, CMF torque constant

    domegadt_CMF = -omega*Kf*np.cos(obl)**2*np.sin(obl)**2 # CL2003, Eq. 71 
    dobldt_CMF = -Kf*np.cos(obl)**3*np.sin(obl) # CL2003, Eq. 74
    
    
###########################
###### Thermal Tide #######
###########################

    ## ExoTides (Salazar & Wordsworth, 2024)

    
    def b_a_sig(omega, m, l):
        q_tilde = torque_analytic_forcingfreq(S, albedo_in, tau_sw, tau_lw, ps, omega, n, m, l,True)
        return(-q_tilde)
    


    
    domegadt_atm = -Ka*(b_a_sig(omega,1,0)*(3/4 * X**2/L**2 * (1-X**2/L**2)) + b_a_sig(omega,1,2)*(3/16 * (1+X/L)**2 * (1-X**2/L**2))+ \
                b_a_sig(omega,1,-2)*(3/16 * (1-X/L)**2 * (1-X**2/L**2) ) + b_a_sig(omega,2,0)*(3/8 * (1-X**2/L**2)**2) + \
                b_a_sig(omega,2,2)*(3/32 * (1+X/L)**4) + b_a_sig(omega,2,-2)*(3/32 * (1-X/L)**4))
    

    
    dobldt_atm = -Ka  * np.sin(obl)/pad_omega * (b_a_sig(omega, 0, -2)*(9/16 * np.sin(obl)**2) + \
                    b_a_sig(omega,1,0)*3/4 * np.cos(obl)**3 - b_a_sig(omega,1,2)*3/16*(1+np.cos(obl))**2 * (2-np.cos(obl)) +\
                    b_a_sig(omega,1,-2)*3/16*(1-np.cos(obl))**2 * (2+np.cos(obl)) + b_a_sig(omega,2,0)*3/8*np.sin(obl)**2*np.cos(obl) +\
                    -b_a_sig(omega,2,2)*3/32*(1+np.cos(obl))**3 +  b_a_sig(omega,2,-2)*3/32*(1-np.cos(obl))**3)

    ## sum torques together ## 

    domegadt = domegadt_grav  + domegadt_CMF + domegadt_atm
    dobldt = dobldt_grav   + dobldt_CMF + dobldt_atm
    return(np.array([C*domegadt*sec_per_year, dobldt*sec_per_year], dtype=float))
    

logging.info("Starting integration")

t_eval = np.arange(0, 4.6, 0.001)*1e9 # years
t_resurface = 4e9 # resurfacing event 
def run_Venus_evolution(delta_t_s, L_o, obl_o, tau_lw, tau_sw, ps, Q, Qn,Q_ocean,tau_lw_Venus, tau_sw_Venus,tau_sw_steam, ramp_rate): 
    ## steam state
    t_steam_end = delta_t_s             # when steam ends (already in seconds)
    t_resurface = 4.0e9                 # resurfacing time
    # t_hab_end will be determined by your event; initialize to t_resurface for now

    # Subset for steam/hab BEFORE you know t_hab_end:
    mask_steam = (t_eval <= t_steam_end)
    t_eval_steam = t_eval[mask_steam]

    mask_hab_pre = (t_eval > t_steam_end) & (t_eval <= t_resurface)
    t_eval_hab = t_eval[mask_hab_pre]

    sol_steam = solve_ivp(spin_evolution_steam, [0,delta_t_s],[L_o, obl_o], method='BDF',max_step=1e8,t_eval=t_eval_steam,args = (rheo_model, tau_sw_steam))
    obl_steam = 180/np.pi * sol_steam.y[1] # rad --> degrees
    omega_steam = sol_steam.y[0]/C
    hab_steam = np.zeros(len(obl_steam)) 
    pst_steam = np.zeros(len(obl_steam)) + 10e5
    
    ## hab state
    Lum_endsteam = L_sun*(1+(2/5)*(1-delta_t_s/to))**(-1)
    S_endsteam = Lum_endsteam/(4*np.pi*a**2)
    sig_albedo = abs(omega_steam[-1]*np.sign(np.cos(np.radians(obl_steam[-1]))) - n) # |sigma|, assume symmetric in LOD
    albedo_endsteam = interp_albedo([sig_albedo, S_endsteam]) # interpolate Yang 2014 50-m data
    #albedo_endsteam = alpha_fit(S_endsteam, omega_steam[-1]*np.sign(np.cos(np.radians(obl_steam[-1]))))

    S_avg_endsteam = S_endsteam*(1-albedo_endsteam)/4
    
    if S_avg_endsteam >= S_RG: # if already not hab
        t_hab_end = delta_t_s
        obl_hab = [] # rad --> degrees
        omega_hab = []
        hab_hab = []
        pst_hab = []
        initial_rot = sol_steam.y[0][-1]/C
        initial_obl = 180/np.pi * sol_steam.y[1][-1]
    else:
        sol_hab = solve_ivp(spin_evolution_hab, [t_eval_hab[0],t_resurface],[C*float(omega_steam[-1]), np.radians(float(obl_steam[-1]))], method='BDF',max_step=1e8,t_eval=t_eval_hab,args = (rheo_model,tau_lw, tau_sw, ps, delta_t_s, Q, Qn), events=[spinout])
        t_hab_end = float(sol_hab.t_events[0][0]) if sol_hab.t_events[0].size > 0 else t_resurface
        obl_hab = 180/np.pi * sol_hab.y[1] # rad --> degrees
        omega_hab = sol_hab.y[0]/C
        initial_rot = omega_hab[-1]
        initial_obl = obl_hab[-1]
        hab_hab = np.zeros(len(obl_hab)) + 1
        pst_hab = np.zeros(len(obl_hab)) + ps
    mask_steam = (t_eval <= t_steam_end)
    mask_hab   = (t_eval > t_steam_end) & (t_eval <= t_hab_end)
    mask_venus = (t_eval > t_hab_end)
    
    mask_steam = (t_eval <= t_steam_end)
    mask_hab   = (t_eval > t_steam_end) & (t_eval <= t_hab_end)
    mask_venus = (t_eval > t_hab_end)

    # Subsets to pass to the last phase
    t_eval_Venus = t_eval[mask_venus]
    ## venus state 
    y0_Venus = np.array([float(C*initial_rot), float(np.radians(initial_obl))], dtype=float)

    sol_Venus = solve_ivp(spin_evolution_Venus, [t_eval_Venus[0],t_eval_Venus[-1]],y0_Venus, method='BDF',max_step=1e8,t_eval=t_eval_Venus,args = (rheo_model,ps,t_hab_end, Q, Qn,Q_ocean, ramp_rate))
    obl_Venus = 180/np.pi * sol_Venus.y[1] # rad --> degrees
    omega_Venus = sol_Venus.y[0]/C
    hab_Venus = np.zeros(len(obl_Venus))-1
    pst_Venus = (ps_Venus-ps) * ((t_eval_Venus-t_hab_end)/(4.6e9-t_hab_end))**ramp_rate + ps

    eq_rotation_rate = np.hstack([omega_steam,omega_hab, omega_Venus])
    eq_obl = np.hstack([obl_steam,obl_hab, obl_Venus]) 
    habitable = np.hstack([hab_steam,hab_hab, hab_Venus])
    ps_t = np.hstack([pst_steam,pst_hab, pst_Venus])

    return(eq_rotation_rate, eq_obl, habitable, ps_t)


all_solutions = pd.DataFrame(columns=['final state','pathway','hab end','Q','Q_ocean','$t_s$','$p_s$', 'ramp','$T_o$', '$\epsilon_o$', '$\tau^{s}_{SW}$','$\epsilon_f$','ever hab?', 't_start', 't_end', '$\Omega_f$', '$\Omega(t)$', '$\epsilon(t)$', '$h(t)$', '$ps(t)$'])
for m in range(N_iterations):
    ## randomly sample ## 
    tau_sw_steam = np.random.uniform(low = tau_sw_steam_min, high = tau_sw_steam_max)
    delta_t_s = 10**(np.random.uniform(low = np.log10(delta_t_s_min), high = np.log10(delta_t_s_max)))*1e9
    ps = 10**(np.random.uniform(low = np.log10(ps_min), high = np.log10(ps_max)))*bar
    # day_length = np.random.uniform(low=day_length_min, high=day_length_max) # uniform in period space
    # L_o = C*2*np.pi/day_length
    L_o = C*(np.random.uniform(low = 2*np.pi/(day_length_min), high = 2*np.pi/(day_length_max))) # uniform in rate space
    day_length_o = 2*np.pi/(L_o/C)
    obl_o = np.random.uniform(low = obl_o_min, high = obl_o_max)
    tau_lw = ps*kappa_lw/g
    tau_sw = ps*kappa_sw/g
    ramp_rate = 10**np.random.uniform(low=np.log10(ramp_min), high=np.log10(ramp_max))
    Q = np.random.uniform(low = Q_min, high = Q_max) # only used for Constant-Q rheology
    Qn = Q+0.1
    Q_ocean = Q*1
    eq_rotation_rate_all, eq_obl_all, habitability_all,ps_t_all = run_Venus_evolution(delta_t_s, L_o, obl_o,tau_lw, tau_sw, ps, Q, Qn,Q_ocean,tau_lw_Venus, tau_sw_Venus,tau_sw_steam, ramp_rate)
    if len(eq_rotation_rate_all) < len(t_eval): # mark crashed solutions 
        final_state = 'crashed'
        pathway = 'crashed'
        endhab = 'crashed'
        new_row = {'final state': final_state,
        'pathway': pathway, 'hab end' : endhab, 'Q':Q, 'Q_ocean':Q_ocean,'$t_s$': delta_t_s, '$p_s$': ps/bar, 'ramp':ramp_rate,
        '$T_o$': day_length_o/60/60/24, '$\epsilon_o$': np.degrees(obl_o), '$\tau^{s}_{SW}$': tau_sw_steam,
        '$\epsilon_f$': eq_obl_all[-1], 'ever hab?': np.nan,
        't_start': np.nan, 't_end': np.nan, 
        '$\Omega_f$': eq_rotation_rate_all[-1], '$\Omega(t)$' : eq_rotation_rate_all[::10],
        '$\epsilon(t)$' : eq_obl_all[::10], '$h(t)$' : habitability_all[::10], '$ps(t)$' : ps_t_all[::10]
        }
    else:
        rotation_period = np.round(abs(2*np.pi/(eq_rotation_rate_all[-1]*60*60*24)),0)            
        if np.any(habitability_all> 0.1): # if ever habitable
            ever_hab = 1
            start_hab = t_eval[np.where(habitability_all > 0.1)[0][0]]/1e9 # start of habitability
            end_hab = t_eval[np.where(habitability_all > 0.1)[-1][-1]]/1e9 # end of habitability
        else: # never habitable
            ever_hab = 0
            start_hab = 0
            end_hab = 0
        T_t = 2*np.pi/(eq_rotation_rate_all[:]*60*60*24) # rotation period
        T_f = T_t[-1] # final rotation period, for sorting
        eps_f = eq_obl_all[-1] # final obliquity [deg]
        hab = habitability_all[:] # dummy habitability variable 
        t_hab = t_eval[hab > 0.1] # all habitable times
        t_end = end_hab*1e9 
        ## specify why habitable state ended, if applicable 
        if len(t_hab) <= 2: # very short-lived
            endhab = 'never hab'
        elif 3.9e9 <= t_end <= 4.1e9: # near resurfacing event
            endhab = 'resurface'
        elif np.mean(np.gradient(abs(T_t[hab > 0.1])[-10:])) < 0: # rotation rate was increasing when habitability ended
            endhab = 'spinout'
        else: # was stable spin rate but Sun got too bright
            endhab = 'solar evolution' 
        ### TIDALLY LOCKED SOLUTIONS ###
        if (abs(abs(np.round(T_f, 0)) - 225)) <= 1 and (abs(abs(np.cos(np.radians(eps_f))) - 1) < 0.01):
            final_state = 'TL'
            if delta_t_s > 0 and (abs(abs(T_t[t_eval <= delta_t_s][-1]) - 225)) <= 1: # became tidally locked during steam state
                pathway = 'steam too long'
            else:
                pathway = 'too close' # started to close to TL equilibria and fell in 

        ### MODERN VENUS SOLUTIONS ###
        elif (abs(abs(np.round(T_f, 0)) - 245)) <= 6 and (abs(abs(np.cos(np.radians(eps_f))) - 1) < 0.01):
            final_state = 'Venus'
            if np.cos(np.radians(eps_f))*np.cos(obl_o) < 0: # final orientation different from initial 
                pathway = 'flipped obliquity'
            elif T_f < 0 :
                pathway = 'flipped rotation'
            elif np.cos(obl_o) < 0:
                pathway = 'started retrograde'
            else:
                pathway = 'other'

        ## PROGRADE VENUS SOLUTIONS ###
        elif (abs(abs(np.round(T_f, 0)) - 77)) <= 6 and (abs(abs(np.cos(np.radians(eps_f))) - 1) < 0.01):
            final_state = 'Prograde Venus'
            if np.cos(np.radians(eps_f))*np.cos(obl_o) < 0:
                pathway = 'flipped obliquity'
            elif T_f < 0 :
                pathway = 'flipped rotation'
            elif np.cos(obl_o) > 0:
                pathway = 'started prograde'
            else:
                pathway = 'other'

        ## NOT EQUILIBRATED SOLUTIONS ##            
        else:
            final_state = 'not eq'
            if len(t_hab) <= 2:
                pathway = 'started too fast'
            elif abs(eps_f-90) < 5: # very rare
                pathway = 'unstable equilibria'
            else:
                pathway = 'other'
        new_row = {'final state': final_state,
        'pathway': pathway, 'hab end' : endhab, 'Q':Q, 'Q_ocean':Q_ocean,'$t_s$': delta_t_s, '$p_s$': ps/bar, 'ramp':ramp_rate,
        '$T_o$': day_length_o/60/60/24, '$\epsilon_o$': np.degrees(obl_o), '$\tau^{s}_{SW}$': tau_sw_steam,
        '$\epsilon_f$': eq_obl_all[-1], 'ever hab?': ever_hab,
        't_start': start_hab, 't_end': end_hab, 
        '$\Omega_f$': eq_rotation_rate_all[-1], '$\Omega(t)$' : eq_rotation_rate_all[::10],
        '$\epsilon(t)$' : eq_obl_all[::10], '$h(t)$' : habitability_all[::10], '$ps(t)$':ps_t_all[::10]
        }
    all_solutions.loc[len(all_solutions)] = new_row
    N_tot = len(all_solutions)
    total_TL = len(all_solutions.loc[all_solutions['final state'] == 'TL'])
    total_Venus = len(all_solutions.loc[all_solutions['final state'] == 'Venus'])
    total_otherVenus = len(all_solutions.loc[all_solutions['final state'] == 'Prograde Venus'])
    total_noteq = len(all_solutions.loc[all_solutions['final state'] == 'not eq'])
    total_crashed = len(all_solutions.loc[all_solutions['final state'] == 'crashed'])
    long_hab_Venus = len(all_solutions.loc[(all_solutions['final state'] == 'Venus') & (all_solutions['t_end']-all_solutions['t_start'] > 1)])

    logging.info(f"Saving .pkl, breakdown: crashed: {np.round(total_crashed/N_tot,2)*100}% | TL: {np.round(total_TL/N_tot,2)*100}% | Not Eq: {np.round(total_noteq/N_tot,2)*100}% | 77 day: {np.round(total_otherVenus/N_tot,2)*100}% | 244 day: {np.round(total_Venus/N_tot,2)*100}% | Long hab: {long_hab_Venus}")
        
all_solutions.to_pickle('./all_solutions.pkl')



      