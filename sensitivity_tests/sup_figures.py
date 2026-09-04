import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as mlines
from constants_functions import *
from scipy.optimize import fsolve
#from scipy.integrate import odeint
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

import sys

path = str(sys.argv[1]) 
S_RG = float(sys.argv[2]) 

bar = 1e5
sec_per_year = 60*60*24*365

plt.rcParams.update({
    # FIGURE & AXES SIZING
    "figure.figsize": (6, 4),  # Set figure size (width, height) in inches
    "figure.dpi": 300,  # High resolution for publication quality
    "axes.titlesize": 14,  # Title font size
    "axes.labelsize": 10,  # Axis label font size
    "axes.labelpad": 3,  # Padding for axis labels

    # TICKS & GRIDLINES
    "xtick.labelsize": 12,  # X-axis tick label font size
    "ytick.labelsize": 12,  # Y-axis tick label font size
    "xtick.major.size": 6,  # Major tick size
    "ytick.major.size": 6,
    "xtick.minor.size": 5,  # Minor tick size
    "ytick.minor.size": 5,
    "xtick.direction": "in",  # Ticks inside the plot
    "ytick.direction": "in",

    # GRID SETTINGS
    "axes.grid": True,  # Enable grid by default
    "grid.color": "gray",  # Grid line color
    "grid.linestyle": "--",  # Dashed grid lines for major ticks
    "grid.linewidth": 0.6,  # Thin grid lines
    "grid.alpha": 0.5,  # Slight transparency for better readability

    # MINOR GRID SETTINGS
    "grid.linestyle": ":",  # Dotted style for minor grid lines
    "grid.linewidth": 0.4,  # Thinner minor grid lines

    # LINES & MARKERS
    "lines.linewidth": 2,  # Line thickness
    "lines.markersize": 6,  # Marker size

    # FONT SETTINGS
    "font.family": "sans-serif",  # Use serif font (e.g., Times New Roman)
    "font.size": 12,  # General font size

    # LEGEND SETTINGS
    "legend.frameon": False,  # No legend border
    "legend.fontsize": 10,
    "legend.loc": "best",

    # SAVEFIG SETTINGS
    "savefig.dpi": 300,  # High resolution for saving figures
    "savefig.transparent": True,  # Transparent background for easy overlay in publications
    "savefig.bbox": "tight",  # Prevents cutting off labels when saving

    # AXES SPINES
    "axes.spines.top": True,  # Show top spine
    "axes.spines.right": True,  # Show right spine
    "axes.spines.left": True,
    "axes.spines.bottom": True,
})

def solar_evolution(time):
    L_sun = 3.828e26 # W
    to = 4.57e9
    t_years = time*1
    Lum = L_sun*(1+(2/5)*(1-t_years/to))**(-1)
    S = Lum/(4*np.pi*a**2)
    return(S)

def alpha_fit(S, Omega):
    sigma = Omega-n_orb
    sigma_crit = 0.744e-6
    alpha_slow = 0.075*(S-1243.3)**(1/4) + 0.18 + 0.08  # based on Tbar
    alpha_rapid = 0.34
    alpha = alpha_rapid + (alpha_slow-alpha_rapid)/np.sqrt(1+(sigma/sigma_crit)**2) # based on dT 
    return(alpha)
def Omega_crit_habitable(time, sign, S_RG):
    S = solar_evolution(time)
    sigma_crit = 0.744e-6 * (S/1366)**(3/4)
    alpha_slow = 0.075*(S-1243.3)**(1/4) + 0.18 + 0.08  # based on Tbar
    alpha_rapid = 0.34
    Omega_end_hab = n_orb + sign*sigma_crit*np.sqrt(((alpha_slow-alpha_rapid)/np.maximum(1-alpha_rapid-(S_RG*4)/S,0))**2 -1)
    return(2*np.pi/(60*60*24*Omega_end_hab))

n_orb=n*1
t_eval = np.linspace(0, 4.6, 5000)*1e9 # years
t_save = t_eval[::10]
all_solutions = pd.read_pickle('./' + path + '/all_solutions.pkl')

## 

TL_all = all_solutions.loc[all_solutions['final state']=='TL']
TL_everhab = TL_all.loc[TL_all['ever hab?'] == 1]
TL_neverhab = TL_all.loc[TL_all['ever hab?'] == 0]

not_eq_all = all_solutions.loc[all_solutions['final state']=='not eq']
not_eq_neverhab = not_eq_all.loc[not_eq_all['ever hab?'] == 0]
not_eq_everhab = not_eq_all.loc[not_eq_all['ever hab?'] == 1]

Venus_all = all_solutions.loc[all_solutions['final state']=='Venus']
Venus_everhab = Venus_all.loc[Venus_all['ever hab?']==1]
Venus_neverhab = Venus_all.loc[Venus_all['ever hab?']==0]

altVenus_all = all_solutions.loc[all_solutions['final state']=='Prograde Venus']
altVenus_everhab = altVenus_all.loc[altVenus_all['ever hab?']==1]
altVenus_neverhab = altVenus_all.loc[altVenus_all['ever hab?']==0]

Venuscolor = 'darkorange'
altVenuscolor = 'skyblue'
TLcolor = 'silver'
noteqcolor = '#1b9e77'

### hab length
from matplotlib.ticker import ScalarFormatter
To_list = np.linspace(4, 400, 1000)
sigma_crit = -2*np.pi/(10*60*60*24) - n
sigma_o = -2*np.pi/(To_list*60*60*24)-n 
delta_ts = 5e7*60*60*24*365
sigma_steam = sigma_o #+ (3/2)*Kg*C*k2/(50*C) * delta_ts
wo = 1.49e-6
qo = -1711
from matplotlib.ticker import ScalarFormatter
fig, ax = plt.subplots(2,1,dpi=150, figsize=(6,7))
delta_t_hab = np.maximum(C*wo/(C*Ka*qo) * (np.log(sigma_steam/sigma_crit) + (1/(2*wo**2))*(sigma_steam**2 - sigma_crit**2)),0)
ax[0].scatter(Venus_everhab['$T_o$'], Venus_everhab['t_end']-Venus_everhab['t_start'], color=Venuscolor, edgecolor='k',lw=0.3)
ax[0].scatter(Venus_neverhab['$T_o$'], Venus_neverhab['t_end']-Venus_neverhab['t_start'], color=Venuscolor, edgecolor='r',lw=0.3)

ax[0].set_ylim(-0.05,2.2)
ax[0].set_ylabel('Duration of habitable state (Gyr)')
ax[0].set_xlabel('Initial rotation period (days)')
ax[0].set_xscale('log',base=2)

ax[0].xaxis.set_major_formatter(ScalarFormatter())
ax[0].xaxis.get_major_formatter().set_scientific(False)
ax[0].set_xticks([1,2,4,8,16,32,64,128,256])
t_end_list = np.linspace(-0.03,2.2,100)

ax[1].scatter(4.6-Venus_everhab['t_end'], Venus_everhab['t_end']-Venus_everhab['t_start'], c=Venuscolor, edgecolor='k',lw=0.3)
ax[1].plot(4.6-t_end_list+0.05, t_end_list, 'k-')
ax[1].set_ylim(-0.03,2.2)
ax[1].set_xlim(5,0)

#ax[1].set_xscale('log')
ax[1].set_ylabel('Duration of habitable state (Gyr)')
ax[1].set_xlabel('Last habitable time (Gya)')

ax[1].axvline(x=0.5, color='r')
ax[1].fill_betweenx(t_end_list, 0.5, 1.4, color='r', alpha=0.1)

plt.savefig('./Figures/' + path +'_hab_scaling.pdf', bbox_inches='tight')


# ## pathways figure 


# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import numpy as np
# import matplotlib.patches as mpatches
# import matplotlib.ticker as ticker

# fig, ax = plt.subplots(2,1,figsize=(10, 9),tight_layout=True)

# def rate_to_period(Omega, _):
#     if Omega == 0:
#         return("∞")
#     period= 2*np.pi/Omega/(60*60*24)
#     return(f"{period:.0f}")

# period_ticks = np.array([-1,-10,-100,-400,np.inf,400,100,10,1])
# period_ticks_top = np.array([400,100,10,1])
# period_ticks_bottom = np.array([1,10,100,400])

# omega_ticks = 2*np.pi/(period_ticks*60*60*24)
# omega_ticks_top = 2*np.pi/(period_ticks_top*60*60*24)
# omega_ticks_bottom = 2*np.pi/(period_ticks_bottom*60*60*24)
# Venuscolor = 'darkorange'
# altVenuscolor = 'skyblue'
# TLcolor = 'silver'
# noteqcolor = '#1b9e77'

# varx ='$t_s$'
# vary ='$T_o$'
# ## tidally locked
# n=1
# size=[15, 15, 15,50]
# everhab = [TL_everhab, not_eq_everhab, altVenus_everhab, Venus_everhab]
# neverhab = [TL_neverhab, not_eq_neverhab, altVenus_neverhab, Venus_neverhab]
# prograde_zorder = [0, 1, 2, 3]
# retrograde_zorder = [0,2,3,2]
# everhab_alpha = 1
# neverhab_alpha = 1
# lw=0.3
# div=2
# everhab_line = 'k'
# neverhab_line='r'
# n_orb = 2*np.pi/(225*60*60*24)
# colors = [TLcolor, noteqcolor, altVenuscolor, Venuscolor]

# # Central plots--scatters
# ax_center_top = ax[0]


# ax_center_bottom = ax[1]



# for i in range(len(everhab)):
#     prograde_everhab = everhab[i].loc[everhab[i]['$\\epsilon_o$'] < 90]
#     retrograde_everhab = everhab[i].loc[everhab[i]['$\\epsilon_o$'] >= 90]
#     prograde_neverhab = neverhab[i].loc[neverhab[i]['$\\epsilon_o$'] < 90]
#     retrograde_neverhab = neverhab[i].loc[neverhab[i]['$\\epsilon_o$'] >= 90]

#     ax_center_top.scatter(prograde_everhab[varx][::n],2*np.pi/(60*60*24*prograde_everhab[vary][::n]), facecolor = colors[i],edgecolors=everhab_line,lw=lw,alpha=everhab_alpha,s=size[i], zorder=prograde_zorder[i])
#     ax_center_top.scatter(prograde_neverhab[varx][::n], 2*np.pi/(60*60*24*prograde_neverhab[vary][::n]), facecolor = colors[i],edgecolors=neverhab_line,lw=lw,alpha=neverhab_alpha,s=size[i], zorder=prograde_zorder[i])
#     ax_center_bottom.scatter(retrograde_everhab[varx][::n], 2*np.pi/(60*60*24*retrograde_everhab[vary][::n]), facecolor = colors[i],edgecolors=everhab_line,lw=lw,alpha=everhab_alpha,s=size[i], zorder=retrograde_zorder[i])
#     ax_center_bottom.scatter(retrograde_neverhab[varx][::n], 2*np.pi/(60*60*24*retrograde_neverhab[vary][::n]), facecolor = colors[i],edgecolors=neverhab_line,lw=lw,alpha=neverhab_alpha,s=size[i], zorder=retrograde_zorder[i])
#     ax_center_top.get_xaxis().set_visible(False)
#     ax_center_top.yaxis.set_major_formatter(FuncFormatter(rate_to_period))
#     ax_center_bottom.yaxis.set_major_formatter(FuncFormatter(rate_to_period))

#     ax_center_top.set_xscale('log')
#     ax_center_top.set_xticks([])
#     ax_center_bottom.set_xscale('log')
#     ax_center_bottom.set_xlim(9e5,4e9)
#     ax_center_top.set_xlim(9e5,4e9)
#     #plt.ylabel(labely)
#     ax_center_bottom.set_xlabel('$\\Delta t_s$ (yr)', fontsize=20)
#     ax_center_top.grid(False)
#     ax_center_bottom.grid(False)

#     fig.text(-0.05, 0.5, '$T_o$ (days)', rotation = -270,va="center", fontsize=20)
    
#     ax_center_top.tick_params(axis='both', which='major', labelsize=12)
#     ax_center_bottom.tick_params(axis='both', which='major', labelsize=12)
    
#     handles = [
#         mlines.Line2D([], [], color=TLcolor, marker='o', linestyle='None', markersize=10, label='Tidally Locked'),
#         mlines.Line2D([], [], color=noteqcolor, marker='o', linestyle='None', markersize=10, label='Not Equilibrated'),
#         mlines.Line2D([], [], color=altVenuscolor, marker='o', linestyle='None', markersize=10, label='Prograde Venus'),
#         mlines.Line2D([], [], color=Venuscolor, marker='o', linestyle='None', markersize=10, label='Venus')
#     ]

#     # Add legend to top central plot (you could also add it to fig or bottom axis if you prefer)
#     ax_center_top.legend(handles=handles, loc=(0,1.02), frameon=True, ncol=4, fontsize=14)


        
#     ax_center_top.set_ylim(1e-6,1e-4)
#     ax_center_top.set_yscale('log')#, linthresh=0.1e-10)
#     ax_center_top.set_yticks(omega_ticks_top)
#     ax_center_top.yaxis.set_major_formatter(FuncFormatter(rate_to_period))


#     ax_center_bottom.set_ylim(1e-4,1e-6)
#     ax_center_bottom.set_yscale('log')#, linthresh=0.1e-10)
#     ax_center_bottom.set_yticks(omega_ticks_bottom)
#     ax_center_bottom.yaxis.set_major_formatter(FuncFormatter(rate_to_period))


    
#     ax_center_top.annotate('$\\epsilon_o < 90^{\\circ}$', xy = (1.25e9,0.25e-6), bbox=dict(boxstyle="round", facecolor="w", alpha=1),fontsize=20)
#     ax_center_bottom.annotate('$\\epsilon_o > 90^{\\circ}$', xy = (1.25e9,0.35e-6), bbox=dict(boxstyle="round", facecolor="w", alpha=1),fontsize=20)

# plt.savefig('./Figures/' + path + '_pathways_fig.pdf', bbox_inches='tight')
# plt.show()
