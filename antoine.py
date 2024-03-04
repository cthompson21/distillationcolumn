import pandas as pd
import os
import nistchempy as nist
import pdb

def read_chart(f_name):
    data = {}
    with open(f_name, 'r') as f:
        header = next(f).split(',')
        for line in f:
            compound, *vals = line.split(',')
            data[compound] = {
                key: val for key, val in zip(header[1:], vals)
            }

    return data


class Antoine:
    def __init__(self, compound, verbose=False):
        from distillation import ROOT_DIR, os
        
        
        x=pd.read_excel('distillation/equilibrium_data/BinaryDistillationSimulator.xlsx',
                        sheet_name='Antoine', index_col=0)
        self.name = compound
        if compound in list(x.index):
            self.A = x['A'][compound]
            self.B = x['B'][compound]
            self.C = x['C'][compound]
            
            self.delta_H_vap,self.delta_H_vap_A,self.delta_H_vap_beta, self.delta_H_vap_Tc, = \
                x.loc[compound,['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc',]]
                
            self.boiling_point=x.loc[compound, 'BP (K)'] 
            self.CpL=x.loc[compound, 'CpL']
            return 
        else:
            print('Compound not found in database!', compound)
            try:
                
                X=nist.Compound(compound, search_option='NAME')  ##hexane
                print("Searching NIST Webbook for", X.name)
                
                tables=X.get_thermo_condensed()
                self.A, self.B, self.C = X.antoine_parameters
                self.delta_H_vap,self.delta_H_vap_A,self.delta_H_vap_beta, self.delta_H_vap_Tc,= X.delta_vap_H_constants.loc[0,['Enthalpy Vaporization (kJ/mol)',
                                                                                                                                'delta_H_vap_A','delta_H_vap_beta',
                                                                                                                    'delta_H_vap_Tc',]]
                
                
                self.delta_H_vap = float(self.delta_H_vap)
                self.delta_H_vap*=1e6  ### convert from kJ/mol to J/kmol
                self.delta_H_vap_A*=1e6  ### convert from kJ/mol to J/kmol
                
                
                X.get_heat_capacities()
                X.get_boiling_point()
                self.CpL = X.Cp_l*1000  ## from J/mol-K to J/kmol-K
                self.boiling_point = X.boiling_point
            except:
                print("Could not find the compound in NIST webbook.  Try a different name")
                raise
                
        
        newcol=pd.DataFrame(columns=x.columns, index = [compound])
        
        newcol.loc[compound, ['A','B','C',]] = self.A, self.B, self.C
        
        newcol.loc[compound, ['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc',]] = \
            X.delta_vap_H_constants.loc[0, ['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc', ]]
        
        newcol.loc[compound, 'BP (K)'] = X.boiling_point
        
        newcol.loc[compound, 'CpL'] = X.Cp_l
        
        x = pd.concat((x,newcol), axis = 0)
        
        x.to_excel('distillation/equilibrium_data/BinaryDistillationSimulator.xlsx',
                        sheet_name='Antoine')

        
        
        if verbose:
            print('Setting Antoine parameters for %s:' % compound)
            print('             A:', self.A)
            print('             B:', self.B)
            print('             C:', self.C)
            
            print('     K value at 300 K, 1 bar= ', self.K_func(300., 1e5))
        return
        


    def K_func(self, T, p):
        """

        :param T: temperature in K
        :param p: pressure in Pa
        :return: K-value for component at specific *T* and *p*
        """
        return 10**(self.A-(self.B/(T+self.C)))
    
    def liquid_enthalpy(self,T):
        return (T-self.boiling_point)*self.CpL
    
    def vapor_enthalpy(self,T):
        return (T-self.boiling_point)*self.CpL + self.delta_H_vap
        
