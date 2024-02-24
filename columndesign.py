# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 10:44:39 2024

@author: chris
"""

import pandas as pd
import numpy as np
from matplotlib.pyplot import *
import distillation.amundson_1958.main as am








def solve_for_model(num, _plot=True):
    print('------------',num)
    global x, y
    model=None

    
#     model = am.Model(
#     components=['n-Butane', 'n-Pentane', ],
#     F=[0,0,0,0,0,1,0,0,0,0], # kmol/h
#     P=3*1e5, # Pa
#     z_feed = [[0,0,0,0,0,0.2,0,0,0,0], [0,0,0,0,0,0.8,0,0,0,0]],
#     RR=0.001,
#     D=0.4,
#     N=num,
#     feed_stage=int(num/2),
# )
    
    F =np.zeros((num,))
    feedstage1 = int(num/3)
    feedstage2 = int(2*num/3)
    F[feedstage1] = 0.5
    F[feedstage2] = 0.5
    
    z = [np.zeros((num,)),np.zeros((num,))]
    z[1][feedstage1] = 0.25
    z[1][feedstage2] = 0.75
    z[0][feedstage1] = 0.75
    z[0][feedstage2] = 0.25
    
    print(z)

    
    
## The initial tutorial
    model = am.Model(
    components=['n-Butane','n-Pentane'],
    F=F, # kmol/h
    P=2*1e5, # Pa
    z_feed = z,
    RR=1.,
    D=0.4,
    

)
    

    
    
    product = 'n-Butane'

    model.add_parameters(verbose=True)
 
    model.T_feed = model.bubble_T_feed(stage=model.feed_stage)   ### XXX need to update this to take in the array
 
    for i in model.stages:
        model.T[i] = model.T_feed

    model.initialize_flow_rates()
   
    model.update_K_values()
    
    for i in model.components:
        model.solve_component_mass_bal(i)

    for stage in model.stages:
        model.T[stage] = model.bubble_T(stage)

    model.solve_energy_balances()
   
    
    iter = 0
   
    while not model.T_is_converged():
        model.update_K_values()
        for i in model.components:
            model.solve_component_mass_bal(i)
            
        for stage in model.stages:
            model.T[stage] = model.bubble_T(stage)
        iter += 1
        
    
    model.solve_energy_balances()
   
    
    outer_loop = 0
    inner_loop = 0
    while not model.flow_rates_converged():
        outer_loop += 1
        for i in model.components:
            model.solve_component_mass_bal(i)
            
        
        for stage in model.stages:
            model.T[stage] = model.bubble_T(stage)
       
        while not model.T_is_converged():
            inner_loop += 1
            model.update_K_values()
            for i in model.components:
                model.solve_component_mass_bal(i)
            
        
            for stage in model.stages:
                model.T[stage] = model.bubble_T(stage)
     
        model.solve_energy_balances()
        
        if outer_loop>16:
            break
        
    x = {}
    y= {}
    for i in model.components:
        x[i] = model.l[i][:]/model.L[:]
        y[i] = x[i]*model.K[i]
    model.y=y
    model.x=x
    
    if _plot:
        
        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(model.stages, model.T, 'o')
        ax.set_xlabel('Stage Number')
        ax.set_ylabel('Temperature [K]')
        
        # plot liquid-phase mole fractions
        fig2 = figure()
        ax2 = fig2.add_subplot(111)
        # calculate mole fractions
        for i in model.components:
            ax2.plot(model.stages, x[i], label=i)
        ax2.set_xticks(range(model.N+1))
        ax2.set_ylabel('Liquid phase mole fraction')
        ax2.set_xlabel('Stage Number')
    
        ax2.legend()


    
    
    return model

close('all')
# model=solve_for_model(16, _plot=True)  
# suptitle('array')         
df = pd.DataFrame( columns = ['Qcond', 'Qreboil', 'hfeed', 'EnergyBalance', 'MinEnergy', 'recovery', 'purity'])
product = 'n-Butane'
for num in np.arange(5,30,2):
    
    model=solve_for_model(num, _plot=True)
    suptitle(str(num)+' stages')
    EnergyBalance = model.Q_reboiler_rule()+model.Q_condenser_rule()-model.h_feed_rule(model.feed_stage)
    MinimumSeparationEnergy=0#-8.314*298*(sum((model.z_feed[component]*np.log(model.z_feed[component]) for component in model.components)))
    
    print("minimum energy (J/mol)", MinimumSeparationEnergy)
    print('Energy input (J/mol)', EnergyBalance/1000)
    
    recovery = model.y[product][1]*model.D/np.sum(model.z_feed[product]*np.array(model.F_feed))
    print('recovery is', recovery)
    purity = model.y[product][1]
    b=[model.Q_condenser_rule()/1000,
                  model.Q_reboiler_rule()/1000,
                        model.h_feed_rule(model.feed_stage)/1000 , 
                        EnergyBalance,
                        MinimumSeparationEnergy,
                        recovery,
                        purity,]
    print(b)
    a = pd.DataFrame(columns=df.columns, data=[b],index = [num],)
    df = pd.concat((df, a),axis=0)
fig, ((ax1,ax2), (ax3,ax4)) = subplots(2,2)

df['EnergyBalance'].plot(ax=ax1)
df['purity'].plot(ax=ax2)    
df['Qcond'].plot(ax=ax3)    
df['Qreboil'].plot(ax=ax4) 

ax1.set_ylabel('energy balance')
ax2.set_ylabel('purity')
ax3.set_ylabel('Qcond')
ax4.set_ylabel('Qreboil')
ax2.set_ylim(0,1)
