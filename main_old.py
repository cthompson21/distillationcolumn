from distillation.bubble_point.calculation import bubble_point
from distillation.solvers import solve_diagonal

from scipy.sparse import linalg, diags
import numpy as np
import pdb
import pandas as pd
from matplotlib.pyplot import *




class Model:
    def __init__(self, components: list = None, F: list = None, P: float = 101325.,
                 z_feed: list = None, RR: float = 1, D: float = 0, N: int = 1, feed_stage: int = 0,
                 T_feed_guess: float = 300.):
        """Distillation column with partial reboiler and total condenser.
        Feed is saturated liquid.

        .. todo::
            implement component liquid flow rates (:attr:`Model.l`) throughout instead of mole fractions

        :param components: list of component names
        :param F: feed molar flow rate
        :param P: pressure (constant throughout column), [Pa]
        :param z_feed: mole fractions of each component, ordered
        :param RR: reflux ratio (L/D)
        :param D: distillate molar flow rate
        :param N: number of equilibrium contacts
        :param feed_stage: stage where feed is input
        :param T_feed_guess: guess temperature for feed stage
        """
        self.flow_rate_tol = 1.e-4
        self.temperature_tol = 1.e-2
        self.components = components
        self.F_feed = F
        self.P_feed = P
        self.p_tot = self.P_feed
        self.z_feed = {key: val for key, val in zip(components, z_feed)}
        self.RR = RR
        self.D = D
        self.B = sum(self.F_feed) - D
        self.N = len(F)-1  ### the index of the last stage (the reboiler)
        
        
        self.K_func = {}
        self.CpL_func = {}
        self.CpV_func = {}
        self.dH_func = {}
        self.T_ref = {}

        # create matrices for variables
        self.num_stages = self.N + 1
        self.stages = range(self.num_stages)
        self.L = np.zeros(self.num_stages)
        self.V = np.zeros(self.num_stages)
        self.L_old = np.zeros(self.num_stages)
        self.V_old = np.zeros(self.num_stages)
        self.T_feed_guess = np.ones((self.num_stages,))*T_feed_guess
        self.T_feed = self.T_feed_guess.copy()
        
        self.F = self.F_feed.copy()
        
        self.feed_stage = np.min(np.where(np.array(self.F_feed)>0)[0])##### XXXX this is a stand in for now
     
        self.x = {
            key: np.ones(self.num_stages)/len(components) for key in components
        }
        self.y = {
            key: np.zeros(self.num_stages)/len(components) for key in components
        }
        self.z = {
            key: np.zeros(self.num_stages) for key in components
        }
        self.l = {
            key: np.zeros(self.num_stages) for key in components
        }
        for component in self.components:
            self.z[component][:] = self.z_feed[component]
            


        
        self.T = np.zeros(self.num_stages)
        self.T_old = np.zeros(self.num_stages)
        self.K = {key: np.zeros(self.num_stages) for key in self.components}

        # solver parameters
        self.df = 1 # Dampening factor to prevent excessive oscillation of temperatures
        self.add_parameters(verbose=True)
        
        
    def get_purity(self, product):
        return self.x[product][0]
    def check_inputs(self):
        """ Checks the feed values and feed composition to avoid errors"""
        
        
        
        if not (np.sum(self.F)-self.D-self.B==0):
            print('error in F,D,B, sum')
            return -1
        
    
        total_mass_comp = np.sum(np.array(list((self.z[c] for c in self.components))),axis = 0)
        non_zero_feed = np.where(self.F>0)[0]
        if np.any(np.abs(total_mass_comp[non_zero_feed]-1)>0.01):
            print('input error in total_mass_comp of feed')
            print(total_mass_comp)
            return -1
        
        
        
        return 0
    
        
    
    def set_feed(self,F,z_feed):
        """allows setting of the feed to re-run the model. """
        
        self.F_feed = F
        
        self.z_feed = {key: val for key, val in zip(self.components, z_feed)}
        self.z = {
            key: np.zeros(self.num_stages) for key in self.components
        }
        self.l = {
            key: np.zeros(self.num_stages) for key in self.components
        }
        for component in self.components:
            self.z[component][:] = self.z_feed[component]
        self.L = np.zeros(self.num_stages)
        self.V = np.zeros(self.num_stages)
        self.L_old = np.zeros(self.num_stages)
        self.V_old = np.zeros(self.num_stages)
        
        self.F = self.F_feed.copy()
        
        self.feed_stage = np.min(np.where(np.array(self.F_feed)>0)[0])
        self.B = sum(self.F_feed) - self.D
        self.T_feed = self.T_feed_guess
        self.T = np.zeros(self.num_stages)
        self.T_old = np.zeros(self.num_stages)
        self.K = {key: np.zeros(self.num_stages) for key in self.components}
        
        
    def add_parameters(self, verbose=False):
        """Add thermodynamic parameters for calculation

        .. note::

            K values from DePriester charts
            CpL from Perrys
            CpV assumes ideal gas
            dH_vap from NIST Webbook

        """
        from distillation.equilibrium_data.depriester_charts import DePriester
        from distillation.equilibrium_data.heat_capacity_liquid import CpL
        from distillation.equilibrium_data.heat_capacity_vapor import CpV
        from distillation.equilibrium_data.heats_of_vaporization import dH_vap
        from distillation.equilibrium_data.antoine import Antoine
        
        
            
        self.K_func = {
            key: DePriester(key, verbose) for key in self.components
        }
        
        
        
        
        # self.K_func =             {
        #         key: Antoine(key, verbose) for key in self.components
        #     }
        
        
        self.CpL_func = {
            key: CpL(key, verbose) for key in self.components
        }
        self.CpV_func = {
            key: CpV(key, verbose) for key in self.components
        }
        self.dH_func = {
            key: dH_vap(key, verbose) for key in self.components
        }
        self.T_ref = {
            key: val.T_ref for key, val in self.dH_func.items()
        }

    def h_pure_rule(self, c, T):
        """rule for liquid enthalpy of pure component"""
        return self.CpL_func[c].integral_dT(self.T_ref[c], T)

    def h_j_rule(self, stage):
        """Enthalpy of liquid on stage *j*.
        Calculated for ideal mixture

        .. math::

            h_j = \\sum_i x_{ij}h^*_i(T_j)

        where the asterisk indicates the pure component enthalpy

        :return: :math:`h_j` [J/kmol]
        """
        return sum(
            self.x_ij_expr(c, stage) * self.h_pure_rule(c, self.T[stage]) for c in self.components
        )

    def x_ij_expr(self, i, j):
        """

        :param i: component name
        :param j: stage number
        :return: mole fraction on stage
        """
        return self.l[i][j] / self.L[j]

    def h_feed_rule(self, stage):
        """Enthalpy of liquid in feed mixture
        Calculated for ideal mixture

        .. math::

            h = \\sum_i x_{ij}h^*_i(T_j)

        where the asterisk indicates the pure component enthalpy

        :return: :math:`h` [J/kmol]
        """
        return sum(
            # self.z[c][stage] * self.h_pure_rule(c, self.T_feed) for c in self.components
            self.z[c][stage] * self.h_pure_rule(c, self.T[stage]) for c in self.components
        )

    def H_pure_rule(self, c, T):
        """Rule for vapor enthalpy of pure component"""
        return self.CpV_func[c].integral_dT(self.T_ref[c], T) + self.dH_func[c].eval()

    def H_j_rule(self, stage):
        """Enthalpy of vapor on stage *j*.
        Calculated for ideal mixture

        .. math::
            H_j = \\sum_i y_{ij}H^*_i(T_j)

        where the asterisk indicates the pure component enthalpy

        .. todo::
            convert y mole fractions to dynamic expression

        :return: :math:`H_j` [J/kmol]
        """
        return sum(
            self.y_ij_expr(c, stage) * self.H_pure_rule(c, self.T[stage]) for c in self.components
        )

    def y_ij_expr(self, i, j):
        """

        :param i: component name
        :param j: stage number
        :return: gas-phase mole fraction on stage
        """
        
        # return self.K[i][j] * self.x_ij_expr(i, j)
        return self.K_func[i].eval_SI(self.T[j], self.P_feed) * self.x_ij_expr(i, j)

    def Q_condenser_rule(self):
        """Condenser requirement can be determined from balances around total condenser"""
        # pdb.set_trace()
        return self.D * (1. + self.RR) * (self.h_j_rule(0) - self.H_j_rule(1))

    def Q_reboiler_rule(self):
        
        """ calculates the reboiler load using energy conservation in the column""" 
        
        return self.D * self.h_j_rule(0) + self.B * self.h_j_rule(self.N) \
               - np.sum(np.array(self.F_feed) * np.array(list(self.h_feed_rule(s) for s in range(self.num_stages))))  \
               - self.Q_condenser_rule()
               


    def step_3_to_step_6(self):
        num_iter = 0
        while not self.T_is_converged():
            self.update_K_values()
            for i in self.components:
                self.solve_component_mass_bal(i)
            self.update_T_values()
            num_iter += 1
        print('while loop exits with %i iterations' % num_iter)

    def run(self):
        self.generate_initial_guess()
        self.step_3_to_step_6()
        self.solve_energy_balances()
        main_loop = 0
        while not self.flow_rates_converged():
            for i in self.components:
                self.solve_component_mass_bal(i)
            self.update_T_values()
            self.step_3_to_step_6()
            self.solve_energy_balances()
            main_loop += 1
           

    def update_K_values(self):
        """
        .. include:: step3.rst

        """
        
        # pdb.set_trace()
        for c in self.components:
            self.K[c][:] = self.K_func[c].eval_SI(self.T[:], self.P_feed)
            # self.K[c][:] = self.K_func[c].eval_SI(self.T[:], self.p_tot)
            

        self.T_old = self.T.copy()

    def update_T_values(self):
        """Update temperatures in all stages
        by performing bubble point calculation

        .. todo::
            vectorize with matrix multiplication

        """
        
        
        
        # update from old calculations
        for i in range(self.num_stages):
            # calculate stage temperature now that all liquid-phase mole fractions are known
            self.T[i] = self.T_old[i] + self.df * (
                    self.bubble_T(i) - self.T_old[i]
            )

    def bubble_T(self, stage):
        l_total = sum(self.l[c][stage] for c in self.components)
        K_vals = [self.K_func[c].eval_SI for c in self.components]
        x_vals = [self.l[c][stage]/l_total for c in self.components]
        
        return bubble_point(x_vals, K_vals, self.P_feed, self.T_old[stage])

    def calculate_T_feed(self):
        self.T_feed = self.bubble_T_feed()
        self.initialize_stage_temperatures()

    def initialize_stage_temperatures(self):
        self.T[:] = self.T_feed

    def bubble_T_feed(self, stage=None):
        """
        .. include: calculate_feed_temperature.rst

        """
        if stage==None:
            print('USING THE HIGHEST FEED STAGE TO SET THE BUBBLE POINT FOR THE FEED.  THIS MAY CAUSE INCORRECT RESULTS')
            stage = self.feed_stage
        
        t_out = np.zeros(self.T_feed.shape)
        for stage in np.where(self.F>0)[0]:
            t_out[stage]=bubble_point(
            [self.z_feed[i][stage] for i in self.components],
            [self.K_func[i].eval_SI for i in self.components], self.P_feed, self.T_feed_guess[stage]
            )
        return t_out

    def initialize_flow_rates(self):
        """
        .. include:: step2.rst

        """
        # initialize L, V with CMO
        
       
        self.L[:] = self.RR * self.D + np.cumsum(self.F_feed)
        
        self.L[self.N] = self.B
        self.L = smooth(self.L)
        # self.V[1:] = self.RR * self.D + self.D
        self.V[1:] = self.L[0:self.N] + self.D -np.cumsum(self.F_feed)[0:self.N]
        self.V[0]=0
        
        # l = np.linspace(self.RR*self.D,self.B, self.N+1)
        # v=np.append([0],self.L[0:self.N]+self.D+np.cumsum(self.F[:self.N]))
        
    
        
        # self.L[:]=l+0.001
        # for comp in self.components:
        #     self.l[comp]=self.L*self.z_feed[comp]
        

    def T_is_converged(self):
        """
        .. include:: temp-converge.rst


        :return: True if T is converged, else False
        """
        eps = np.abs(self.T - self.T_old)
        return eps.max() < self.temperature_tol

    def solve_component_mass_bal(self, component):
        """Solve component mass balances

        .. todo:
            dont need to calculate D and A here as they are constant

        """
        A, B, C, D = make_ABC(
            self.V, self.L, self.K[component], self.F, self.z[component], self.D, self.B, self.N, None
        )
        
        result = solve_diagonal(A, B, C, D)
        ### i have made the result back into the x values, rather than the L values
        self.x[component][:]=result
        self.y[component][:]=self.x[component][:]*self.K[component]
        self.l[component][:] = result*self.L
        
        
        
        if np.any(np.isnan(result)):
            print("Nan result in solvecomponentmassbal")
            return -1
        
        self.L = sum(list((self.l[c] for c in self.components)))
        
    def mass_fraction_sum_check(self):### mass fractions add up to 1 check
    
        total_mass_frac = sum(list((self.x[c] for c in self.components)))
        if np.any(np.abs(total_mass_frac[:-1]-1)>0.03):
            # print('x')
            # print(self.x)
            # print('y')
            # print(self.y)
            return -1
        total_mass_frac = sum(list((self.y[c] for c in self.components)))
        if np.any(np.abs(total_mass_frac[:-1]-1)>0.03):
            # print('x')
            # print(self.x)
            # print('y')
            # print(self.y)
            return -1
        
        return 0
        
        
    def update_flow_rates(self):
        for i in self.stages:
            self.L[i] = sum(self.l[c][i] for c in self.components)

        self.V[0] = 0.  # total condenser
        self.V[1] = (self.RR + 1.) * self.D
        for i in range(2, self.num_stages):
            self.V[i] = self.L[i - 1] + self.D - sum(self.F[k] for k in range(i))

    def solve_energy_balances(self):
        """Solve energy balances"""

        self.L_old[:] = self.L.copy()  #### chris changed this to the copy
        self.V_old[:] = self.V.copy() ## chris changed this to the copy


        #### lower case in paper is liquid
        
        ### "mat" is the H matrix in the paper
        
        mat=np.ndarray((self.num_stages, self.num_stages))
        for p in range(1,self.N+1):
          
            mat[p,p]=self.h_j_rule(p-1) - self.H_j_rule(p)  
        for p in range(1,self.N):
            mat[p, p+1]=self.H_j_rule(p+1) - self.h_j_rule(p)

        
        G_array = np.ndarray((self.N+1))
        for p in range(1,self.N):
       
            G_array[p] = (-np.cumsum(self.F)[p-1]+self.D)*self.h_j_rule(p-1) + (np.cumsum(self.F)[p]-self.D)*self.h_j_rule(p) - self.F[p]*self.h_feed_rule(p)
        for p in range(self.N,self.N+1):
            G_array[p] = self.B*(self.h_j_rule(p) - self.h_j_rule(p-1))-self.Q_reboiler_rule()
        

        diagonal = np.array(list(mat[p,p] for p in range(0,self.N+1)))
        upperdiagonal = np.array(list(mat[p-1,p] for p in range(1,self.N+1)))
        lowerdiagonal = np.array(list(mat[p,p-1] for p in range(1,self.N+1)))
        
        
       
        A = diags(
            diagonals=[diagonal[1:], upperdiagonal[1:]],
            offsets=[0, 1],
            shape=(self.N, self.N),
            format='csr'
        )
        self.V[1:] = linalg.spsolve(A, G_array[1:])
        ###############end chris mod
        self.L[0] = self.RR * self.D
       
         
        for i in range(1, self.N):
            self.L[i] = self.V[i + 1] - self.D + np.sum(self.F[:i+1])   ### chris modified the sum to not include a loop
        for i in range(self.N, self.N+1): ### added on 3/3
            self.L[i]=self.B-self.V[i]


    def flow_rates_converged(self):
        """Determine if flow rates are converged

        Use the mathematical criterion in :meth:`Model.is_below_relative_error`

        """
        return self.is_below_relative_error(self.L, self.L_old) and self.is_below_relative_error(self.V[1:],
                                                                                                 self.V_old[1:])

    def is_below_relative_error(self, new, old):
        """Determine relative error between two vectors

        .. math::

            \\sqrt{\\left(\\frac{X_{\\mathrm{new}} - X_{\\mathrm{old}}}{X_{\\mathrm{new}}}\\right)^2} < \\epsilon

        The flow rate tolerance, :math:`\epsilon`,
        is found in the attribute :attr:`distillation.amundson_1958.main.Model.flow_rate_tol`

        :param new:
        :param old:
        :rtype: bool

        """
        return np.abs((new - old) / new).max() < self.flow_rate_tol
    
    def solve_self(self,reinit =True):
        if self.check_inputs():
            return -1
        
     
        if reinit:
            self.T_feed = self.bubble_T_feed(stage=self.feed_stage)   ### XXX need to update this to take in the array
     
            for i in self.stages:
                if self.T_feed[i] ==0:
                    self.T[i]=np.mean(self.T_feed[self.T_feed!=0])
                else:
                    self.T[i] = self.T_feed[i]
            
    
            self.initialize_flow_rates()
            
           
            self.update_K_values()
        
        for i in self.components:
            self.solve_component_mass_bal(i)
       

        for stage in self.stages:
            self.T[stage]=(self.bubble_T(stage)+self.T[stage])/2

        self.solve_energy_balances()
       
        
        iter = 0
       
        while not self.T_is_converged():
            self.update_K_values()
            for i in self.components:
                self.solve_component_mass_bal(i)
            for stage in self.stages:
                self.T[stage]=(self.bubble_T(stage)+self.T[stage])/2
            iter += 1
            
        
        self.solve_energy_balances()
       
        
        outer_loop = 0
        inner_loop = 0
        while not self.flow_rates_converged():
            outer_loop += 1
            while not self.T_is_converged():
                inner_loop += 1
                self.update_K_values()
                for i in self.components:
                    self.solve_component_mass_bal(i)
                for stage in self.stages:
                    self.T[stage] = (self.bubble_T(stage)+self.T[stage])/2
         
            self.solve_energy_balances()
            
            if outer_loop>256:
                break
            
        
            
        
        
    def plot_molefractions(self):
        

        # plot liquid-phase mole fractions
        fig2, (ax1,ax2,ax3,ax4) = subplots(1,4)
        
        # calculate mole fractions
        for i in self.components:
            ax3.plot(self.stages, self.x[i], label=i)
        ax3.legend()
            
        
        ax2.plot(self.stages,self.L, label='L')
        ax2.plot(self.stages, self.V, label='V')
        ax2.set_xticks(range(self.N+1))
        ax2.set_ylabel('Liquid phase mole fraction')
        ax2.set_xlabel('Stage Number')
        
        ax1.plot(self.stages, self.T, 'o')
        ax1.set_xlabel('Stage Number')
        ax1.set_ylabel('Temperature [K]')
        
        for idx in np.where(self.F>0)[0]:
            
            ax3.plot([idx,], [self.z[self.components[0]][idx]], 'or')
            
        
        
    
        ax2.legend()
        return fig2, (ax1,ax2,ax3,ax4) 



def make_ABC(V: np.array, L: np.array, K: np.array, F: np.array, z: np.array,
             Distillate: float, Bottoms: float, N: int, feed_stage: int):
    """
    Distillation column with partial reboiler and total condenser

    .. note::
        K_j is assumed to depend on *T* and *p*, but not composition

    :param V: vapor molar flow rate out of stage 0 to *N*
    :param L: liquid molar flow rate out of stage 0 to *N*
    :param K: equilibrium expressions for stage 0 to *N*
    :param F: feed flow rate into stage for stage 0 to *N*
    :param z: feed composition into stage for stage 0 to *N*
    :param Distillate: distillate flow rate
    :param Bottoms: bottoms flow rate
    :param N: number of equilibrium stages

    :return: A, B, C, D
    """

    
    W=Bottoms
    D=Distillate
    
    assert abs(V[0]) < 1e-8, 'Vapor flow rate out of total condenser is non-zero!'
    
    
    mat = np.ndarray((N+1, N+1))
    
    mat[0,0]=   -V[1] - (K[0]-1)*D
    
    for p in range(1,N):
        mat[p,p]= -(V[p+1]-D+V[p]*K[p])-np.cumsum(F)[p]
    
    for p in range(N, N+1):
        mat[p,p]= -(W+V[p]*K[p])
        
    for p in range(1,N+1):
        mat[p-1, p]=V[p]*K[p]
    for p in range(1,N+1):
        mat[p, p-1]=V[p]-D+np.cumsum(F)[p-1]
    
    diagonal = np.array(list(mat[p,p] for p in range(0,N+1)))
    upperdiagonal = np.array(list(mat[p-1,p] for p in range(1,N+1)))
    lowerdiagonal = np.array(list(mat[p,p-1] for p in range(1,N+1)))
    feedarray = -np.array(F[:]) * z[:]
    

    
    A=lowerdiagonal
    B=diagonal
    C=upperdiagonal
    D=feedarray
    
    return A, B, C, D


def solve_component_mass_balances(*args):
    """
    Distillation column with partial reboiler and total condenser

    .. note::
        K_j is assumed to depend on *T* and *p*, but not composition

    :param V: vapor molar flow rate out of stage 0 to *N*
    :param L: liquid molar flow rate out of stage 0 to *N*
    :param K: equilibrium expressions for stage 0 to *N*
    :param F: feed flow rate into stage for stage 0 to *N*
    :param z: feed composition into stage for stage 0 to *N*
    :param Distillate: distillate flow rate
    :param Bottoms: bottoms flow rate
    :param N: number of equilibrium stages
    :return: l
    """
    A, B, C, D = make_ABC(*args)
    
    return solve_diagonal(A, B, C, D)

def smooth(spectrum,window_len=5,window='flat'):
    
    if spectrum.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return spectrum
    
    
     
    #moving average
    w=np.ones(window_len,'d')
    s=np.r_[spectrum[window_len-1:0:-1],spectrum[:],spectrum[-1:-window_len:-1]]        
    out=np.convolve(w/w.sum(),s,mode='valid')[int((window_len-1)/2):int(-(window_len-1)/2)]
    
    return out
if __name__ == '__main__':
    import os
    os.chdir('C:/Users/chris/My Drive/PyScripts/ViaPy')
    feedstage1=6
    feedstage2=6
    num_stages=12
    R=1
    
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
        z[1][feedstage1] = 0.75
        z[1][feedstage2] = 0.25
        z[0][feedstage1] += (1-z[1][feedstage1])
        z[0][feedstage2] += (1-z[1][feedstage2])
 
    

    model = Model(
    components=['n-Butane','n-Pentane'],
    F=F, # kmol/h
    P=2*1e6, # Pa
    z_feed = z,
    RR=R,
    D=0.4,
    )
    
    print(model.h_j_rule(5))
    model.solve_self()
    model.plot_molefractions()
    