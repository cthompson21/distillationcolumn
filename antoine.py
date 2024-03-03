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
            
            self.delta_H_vap,self.delta_H_vap_A,self.delta_H_vap_beta, self.delta_H_vap_Tc = \
                x.loc[compound,['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc']]
            return 
        else:
            print('compound not found in antoine!', compound)
            try:
                
                X=nist.Compound(compound, search_option='NAME')  ##hexane
                print("Searching NIST Webbook for", X.name)
                
                tables=X.get_thermo_condensed()
                self.A, self.B, self.C = X.antoine_parameters
                self.delta_H_vap,self.delta_H_vap_A,self.delta_H_vap_beta, self.delta_H_vap_Tc = X.delta_vap_H_constants.loc[0,['Enthalpy Vaporization (kJ/mol)',
                                                                                                                                 'delta_H_vap_A','delta_H_vap_beta',
                                                                                                                                 'delta_H_vap_Tc']]
            except:
                print("Could not find the compound in NIST webbook.  Try a different name")
                raise
                
        
        newcol=pd.DataFrame(columns=x.columns, index = [compound])
        
        newcol.loc[compound, ['A','B','C',]] = self.A, self.B, self.C
        
        newcol.loc[compound, ['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc',]] = \
            X.delta_vap_H_constants.loc[0, ['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc']]
        
        x = pd.concat((x,newcol), axis = 0)
        
        x.to_excel('distillation/equilibrium_data/BinaryDistillationSimulator.xlsx',
                        sheet_name='Antoine')

        
        
        if verbose:
            print('Setting Antoine parameters for %s:' % compound)
            print('             A:', self.A)
            print('             B:', self.B)
            print('             C:', self.C)
            
            print('     K value at 300 K, 1 bar= ', self.eval_SI(300., 1e5))
        return
        


    def eval_SI(self, T, p):
        """

        :param T: temperature in K
        :param p: pressure in Pa
        :return: K-value for component at specific *T* and *p*
        """
        return 10**(self.A-(self.B/(T+self.C)))
        
