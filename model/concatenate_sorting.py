import pandas as pd
import numpy as np

import sys

num = int(sys.argv[1])

all_solutions = pd.DataFrame(columns=['final state','pathway','hab end','Q','Q_ocean','$t_s$','$p_s$', '$T_o$', '$\epsilon_o$', '$\tau^{s}_{SW}$','$\epsilon_f$','ever hab?', 't_start', 't_end', '$\Omega_f$', '$\Omega(t)$', '$\epsilon(t)$', '$h(t)$'])

for i in range(1,num+1):
    new_sol = pd.read_pickle('./num_' + str(i) + '/' + 'all_solutions.pkl')
    all_solutions = pd.concat([all_solutions, new_sol])
all_solutions.to_pickle('./all_solutions.pkl')

TL_all = all_solutions.loc[all_solutions['final state']=='TL']
TL_everhab = TL_all.loc[TL_all['ever hab?'] == 1]

not_eq_all = all_solutions.loc[all_solutions['final state']=='not eq']
noteq_everhab = not_eq_all.loc[not_eq_all['ever hab?'] == 1]

altVenus_all = all_solutions.loc[all_solutions['final state']=='Prograde Venus']
altVenus_everhab = altVenus_all.loc[altVenus_all['ever hab?'] == 1]

crashed_all = all_solutions.loc[all_solutions['final state']=='crashed']

Venus_all = all_solutions.loc[all_solutions['final state']=='Venus']
Venus_everhab = Venus_all.loc[Venus_all['ever hab?'] == 1]
Venus_neverhab = Venus_all.loc[Venus_all['ever hab?'] == 0]
Venus_longhab = Venus_everhab.loc[Venus_everhab['t_end']-Venus_everhab['t_start']>1]

## summary 

print('crashed: ' + str(len(crashed_all)/len(all_solutions) * 100))

print('TL: ' + str(len(TL_all)/len(all_solutions) * 100) + ' % (' + str(len(TL_everhab)/len(TL_all) * 100) + ' %)')
print('Not Eq.: ' + str(len(not_eq_all)/len(all_solutions) * 100)+ ' % (' + str(len(noteq_everhab)/len(not_eq_all) * 100) + ' %)')
print('Prograde Venus: ' + str(len(altVenus_all)/len(all_solutions) * 100) + ' % (' + str(len(altVenus_everhab)/len(altVenus_all) * 100) + ' %)')
print('Venus: ' + str(len(Venus_all)/len(all_solutions) * 100) + ' %')
print('Venus ever hab?: ' + str(len(Venus_everhab)/len(Venus_all) * 100) + ' %')
print('Venus never hab?: ' + str(len(Venus_neverhab)/len(Venus_all) * 100) + ' %')
print('Venus hab > 1 Gyr: ' + str(len(Venus_longhab)/len(Venus_all) * 100)) 
print('Avg. t_hab Venus: ' + str(np.mean(Venus_everhab['t_end']-Venus_everhab['t_start'])))
print('Max. t_hab Venus: ' + str(np.max(Venus_everhab['t_end']-Venus_everhab['t_start'])))

