# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 10:44:39 2024

@author: chris
"""

import pandas as pd
import numpy as np
from matplotlib.pyplot import *
import distillation.amundson_1958.main_old as am
import itertools
import seaborn as sns






def solve_for_model(num,R,feedstage1,feedstage2,  _plot=True):
    print('------------',num)
    global x, y
    model=None

    

    
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
        z[1][feedstage1] = 0.25
        z[1][feedstage2] = 0.75
        z[0][feedstage1] += (1-z[1][feedstage1])
        z[0][feedstage2] += (1-z[1][feedstage2])
 
    

    model = am.Model(
    components=['n-Butane','n-Pentane'],
    F=F, # kmol/h
    P=2*1e6, # Pa
    z_feed = z,
    RR=R,
    D=0.4,
    )
    
    
    
    model.solve_self()

    return model



close('all')

# model=solve_for_model(16, _plot=True)  
# suptitle('array')         
df = pd.DataFrame( columns = ['Qcond', 'Qreboil', 'hfeed',
                              'TotalEnergyInput', 'MinEnergy', 'recovery',
                              'purity', 'reflux_ratio', 'num_stages', 'feedstage1', 'feedstage2'])

product = 'n-Butane'

num_stages = 18
for reflux_ratio in [1]:
    model=solve_for_model(num_stages,reflux_ratio
                          ,6,6, _plot=True)
    
    
    
    # purity=0.95
    # while model.get_purity('n-Butane')<purity:
    #     print(model.get_purity('n-Butane'), model.RR)
    #     model.RR+=0.1
    #     model.solve_self(reinit=False)
        
        
    
    # if model.get_purity('n-Butane')>purity:
    #     print("got it")
    model.plot_molefractions()
    
    
    
    
    # for num in np.arange(5,12,1):
    # for feedstage1, feedstage2 in itertools.product(np.arange(2,11,1),[10,], ):
    for feedstage1, feedstage2 in itertools.product(np.arange(5,14,1),[13,], ):
        model.RR = 1.1
        
        
        # model.num_stage = num_stages
        # model.RR= reflux_ratio
        F =np.zeros((num_stages,))
        
        if feedstage1==feedstage2:
            F[feedstage1] = 1
        else:    
            F[feedstage1] = 0.5
            F[feedstage2] = 0.5
        # F[feedstage2] += (1-F[feedstage1])
        
        z = [np.zeros((num_stages,)),np.zeros((num_stages,))]
        if feedstage1==feedstage2:
            z[1][feedstage1] = 0.5
            z[0][feedstage1] = 0.5
        else:
            z[1][feedstage1] = 0.25  ### the higher feed stage should be poor in the heavy component
            z[1][feedstage2] = 0.75 ### the lower  feed stage should be richer in the heavy component
            z[0][feedstage1] += (1-z[1][feedstage1])
            z[0][feedstage2] += (1-z[1][feedstage2])
            
        print("the sum of feed is ",np.sum(z[1]*F), np.sum(z[0]*F))
        if np.any(z[0]>1):
            print('z0')
            break
        if np.any(z[0]>1):
            print('z1')
            break
        if np.any(F>1):
            print('F')
            break
        print(sum(F), max(z[0]), max(z[1]))
        model.set_feed(F,z)
        model.solve_self()
        
        purity=0.95
        delta=1
        while delta>0.001:
            delta= purity -model.get_purity('n-Butane')
            print(round(delta,4), rount(model.RR,3))
            if np.abs(delta)>0.2:
                model.RR+=0.1
            else: 
                model.RR+=0.02
            model.solve_self(reinit=False)
            
        if model.get_purity('n-Butane')>purity:
            print("got it")
        
        model.plot_molefractions()
        
        # suptitle(str(num_stages)+' stages')
        # EnergyBalance = model.Q_reboiler_rule()+model.Q_condenser_rule()-model.h_feed_rule(model.feed_stage)
        TotalEnergyInput = model.Q_reboiler_rule()-model.Q_condenser_rule()
        MinimumSeparationEnergy=0#-8.314*298*(sum((model.z_feed[component]*np.log(model.z_feed[component]) for component in model.components)))
        
        print("TotalEnergyInput (J/mol)", TotalEnergyInput/sum(model.F))
       
        recovery = model.y[product][1]*model.D/np.sum(model.z_feed[product]*np.array(model.F_feed))
        print('recovery is', recovery)
        purity = model.y[product][1]
        b=[model.Q_condenser_rule()/1000,
                      model.Q_reboiler_rule()/1000,
                            model.h_feed_rule(model.feed_stage)/1000 , 
                            TotalEnergyInput,
                            MinimumSeparationEnergy,
                            recovery,
                            purity,
                            reflux_ratio,
                            num_stages,
                            feedstage1,
                            feedstage2,
                           ]
        
        a = pd.DataFrame(columns=df.columns, data=[b],index = [df.shape[0]])
        df = pd.concat((df, a),axis=0)
    
df.index = np.arange(df.shape[0])
fig, ((ax1,ax2), (ax3,ax4)) = subplots(2,2)

df.to_excel('C:/Users/cthompson/Desktop/DistillationModel1.xlsx')
sns.lineplot(data=df, x='feedstage1', y='TotalEnergyInput', hue='reflux_ratio', ax=ax1)
sns.lineplot(data=df, x='feedstage1', y='purity', hue='reflux_ratio', ax=ax2)
sns.lineplot(data=df, x='feedstage1', y='Qcond', hue='reflux_ratio', ax=ax3)
# sns.lineplot(data=df, x='num_stages', y='Qreboil', hue='reflux_ratio', ax=ax4)
sns.lineplot(data=df, x='feedstage1', y='TotalEnergyInput', hue='reflux_ratio', ax=ax4)
ax1.set_ylabel('energy balance')
ax2.set_ylabel('purity')
ax3.set_ylabel('Qcond')
ax4.set_ylabel('Qreboil')
ax2.set_ylim(0,1)

fig, ((ax1,ax2), (ax3,ax4)) = subplots(2,2)

sns.lineplot(data=df, x='feedstage2', y='TotalEnergyInput', hue='reflux_ratio', ax=ax1)
sns.lineplot(data=df, x='feedstage2', y='purity', hue='reflux_ratio', ax=ax2)
sns.lineplot(data=df, x='feedstage2', y='Qcond', hue='reflux_ratio', ax=ax3)
# sns.lineplot(data=df, x='num_stages', y='Qreboil', hue='reflux_ratio', ax=ax4)
sns.lineplot(data=df, x='feedstage2', y='TotalEnergyInput', hue='reflux_ratio', ax=ax4)
ax1.set_ylabel('energy balance')
ax2.set_ylabel('purity')
ax3.set_ylabel('Qcond')
ax4.set_ylabel('Qreboil')
ax2.set_ylim(0,1)



