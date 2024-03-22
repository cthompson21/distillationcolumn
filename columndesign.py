# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 10:44:39 2024

@author: chris
"""

import pandas as pd
import numpy as np
from matplotlib.pyplot import *
import distillation.amundson_1958.main as am
import itertools
import seaborn as sns
import os
import pdb




close('all')
# model=solve_for_model(16, _plot=True)  
# suptitle('array')         
df = pd.DataFrame( columns = ['Qcond', 'Qreboil', 'hfeed',
                              'TotalEnergyInput_kJ_mol', 'MinEnergy', 'recovery',
                              'purity', 'reflux_ratio', 'num_stages', 'feedstage1', 'feedstage2'])

product = 'butane'
os.chdir('H:/My Drive/PyScripts/ViaPy')
num_stages = 24
for repeat in [0]:#np.arange(1,2.5,0.5):
    F =np.zeros((num_stages,))
    # feedstage1=feedstage2=6
    # if feedstage1==feedstage2:
    #     F[feedstage1] = 1
    # else:    
    #     F[feedstage1] = 0.5
    # F[feedstage2] += (1-F[feedstage1])
    
    # z = [np.zeros((num_stages,)),np.zeros((num_stages,))]
    # if feedstage1==feedstage2:
    #     z[1][feedstage1] = 0.5
    #     z[0][feedstage1] = 0.5
    # else:
    #     z[1][feedstage1] = 0.95
    #     z[1][feedstage2] = 0.05
    #     z[0][feedstage1] += (1-z[1][feedstage1])
    #     z[0][feedstage2] += (1-z[1][feedstage2])
    
    # model = am.Model(
    # components=['butane','pentane'],
    # F=F, # kmol/h
    # P=4e5, # Pa
    # z_feed = z,
    # RR=15,
    # D=0.4,
    # )
    
    # pur = 0.98
    # model.solve_for_min_purity(pur,'butane')
        
    
    # if model.get_purity('butane')>pur:
    #     print("got it")
    # model.plot_molefractions()
    
    
    # for num in np.arange(5,12,1):
    
    for feedstage1, feedstage2 in itertools.product(np.arange(6,18,1), np.arange(6,18,1)):
        if feedstage1<feedstage2:
            continue
        
       
        # model.num_stage = num_stages
        # model.RR= reflux_ratio
        F =np.zeros((num_stages,))
        
        if feedstage1==feedstage2:
            F[feedstage1] = 1
        else:    
            F[feedstage1] = 0.5
        F[feedstage2] += (1-F[feedstage1])
        
        z = [np.zeros((num_stages,)),np.zeros((num_stages,))]
        if feedstage1==feedstage2:
            z[1][feedstage1] = 0.5
            z[0][feedstage1] = 0.5
        else:
            z[1][feedstage1] = 0.85
            z[1][feedstage2] = 0.15
            z[0][feedstage1] += (1-z[1][feedstage1])
            z[0][feedstage2] += (1-z[1][feedstage2])
        if np.any(z[0]>1):
            print('z0')
            break
        if np.any(z[0]>1):
            print('z1')
            break
        if np.any(F>1):
            print('F')
            break
        
        
        model = am.Model(
        components=['acetone','2-propanol'],
        F=F, # kmol/h
        P=4e5, # Pa
        z_feed = z,
        RR=1,
        D=0.4,
        )
        product='acetone'
        
        model.solve_self()
        
        pur=0.98
        model.solve_for_min_purity(pur, product)
            
        if model.get_purity(product)>pur:
            print("got it")
        
        model.plot_molefractions()
        
        
        
        

        TotalEnergyInput = (model.Q_reboiler_rule()-model.Q_condenser_rule())/1e6
        MinimumSeparationEnergy=0#-8.314*298*(sum((model.z_feed[component]*np.log(model.z_feed[component]) for component in model.components)))
        
        print("TotalEnergyInput: {0:.3f} kJ/mol".format(TotalEnergyInput))
       
        recovery = model.y[product][1]*model.D/np.sum(model.z_feed[product]*model.F_feed)
        print("Recovery: {0:.4f} ".format(recovery))
        purity = model.y[product][1]
        b=[model.Q_condenser_rule()/1000,
                      model.Q_reboiler_rule()/1000,
                            model.h_feed_rule(model.feed_stage)/1000 , 
                            TotalEnergyInput,
                            MinimumSeparationEnergy,
                            recovery,
                            purity,
                            model.RR,
                            num_stages,
                            feedstage1,
                            feedstage2,
                           ]
        
        a = pd.DataFrame(columns=df.columns, data=[b],index = [df.shape[0]])
        df = pd.concat((df, a),axis=0)
    
df.index = np.arange(df.shape[0])
fig, ((ax1,ax2), (ax3,ax4)) = subplots(2,2)



df.to_excel('C:/Users/cthompson/Desktop/DistillationModel1.xlsx')
sns.lineplot(data=df, x='feedstage1', y='TotalEnergyInput_kJ_mol', hue='feedstage2', ax=ax1)
sns.lineplot(data=df, x='feedstage1', y='purity', hue='feedstage2', ax=ax2)
sns.lineplot(data=df, x='feedstage1', y='Qcond', hue='feedstage2', ax=ax3)
sns.lineplot(data=df, x='feedstage1', y='Qreboil', hue='feedstage2',  ax=ax4)
ax1.set_ylabel('energy balance')
ax2.set_ylabel('purity')
ax3.set_ylabel('Qcond')
ax4.set_ylabel('Qreboil')
ax2.set_ylim(0,1)

fig, ((ax1,ax2), (ax3,ax4)) = subplots(2,2)

sns.lineplot(data=df, x='feedstage2', y='TotalEnergyInput_kJ_mol', hue='feedstage1', ax=ax1)
sns.lineplot(data=df, x='feedstage2', y='purity', hue='feedstage1', ax=ax2)
sns.lineplot(data=df, x='feedstage2', y='Qcond', hue='feedstage1', ax=ax3)
sns.lineplot(data=df, x='feedstage2', y='Qreboil', hue='feedstage1', ax=ax4)
ax1.set_ylabel('energy balance')
ax2.set_ylabel('purity')
ax3.set_ylabel('Qcond')
ax4.set_ylabel('Qreboil')
ax2.set_ylim(0,1)



