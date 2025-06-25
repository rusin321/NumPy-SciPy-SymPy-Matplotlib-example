import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import sympy as smp
from sympy import Symbol
from sympy import evaluate
from IPython.display import display
import sys
from sympy import init_printing


###Functions
#Defining differential equation
def dS_dr(r, S):
    phi, dphi_dr = S
    
    return [dphi_dr, dV_dphi_num(phi) - 2/r*dphi_dr]

def position_event(t, y):
    return y[0] + 0.1*T_init

def velocity_event(t, y):
    return y[1] - 0.00001

def t_max_find(t_min, t_max, N, phi_max):
    
    r = np.linspace(t_min, t_max, N)
    S = (phi_max, v_0)
    
    # Allow the event to stop the solver
    #event.terminal = True
    # Detect when y reaches given value from either direction
    #event.direction = 0
    
    position_event.terminal = True
    position_event.direction = 0
    
    velocity_event.terminal = True
    velocity_event.direction = 0
    
    sol = solve_ivp(dS_dr, t_span = (t_min, t_max), y0 = S, t_eval = r, events = [position_event, velocity_event])
    
    #plt.plot(np.linspace(t_min, len(sol.y[0])/N*t_max, len(sol.y[0])), sol.y[0])
    #plt.show()
    
    new_t_max = len(sol.y[0])
    
    #print("#")
    #print(new_t_max)
    
    return new_t_max/N
    

def bisection(phi_start, phi_end, t_min, t_max, r, phi_0, iteration = 0):
    #print(iteration)
    
    middle = (phi_end + phi_start)/2
    
    position_event.terminal = True
    position_event.direction = 0
    
    velocity_event.terminal = True
    velocity_event.direction = 0
    
    S = (middle, v_0)
    sol = solve_ivp(dS_dr, t_span = (t_min, t_max), y0 = S, t_eval = r, events = [position_event, velocity_event])
    
    phi_0.append([phi_start, phi_end])
    
    #print([phi_start, phi_end])
      
    if abs(middle - phi_start)/middle > 0.000001:
        if iteration < 50:
            if all(i >= 0 for i in sol.y[0]) == True:
                bisection(middle, phi_end, t_min, t_max, r, phi_0, iteration = iteration + 1)
                
                #phi_0.append([new_phi_start, new_phi_end])
                #print(1)

            elif any(i < 0 for i in sol.y[0]) == True:
                bisection(phi_start, middle, t_min, t_max, r, phi_0, iteration = iteration + 1)
                
                #phi_0.append([phi_start, new_phi_start])
                #print(2)
    else:
        pass
        #print("break")

def fun(T_init, e_init):
    limit_to_0 = smp.limit(V.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), phi, 0).evalf()
    V_norm = V - limit_to_0

    dV_dphi = smp.diff(V_norm, phi)

    V_norm_num = smp.utilities.lambdify(phi, V_norm.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), "numpy")
    
    #Plots
    #smp.plotting.plot((V.subs({e:0.3, a_b:a_b_2, pi:np.pi, T:T_init}), (phi, 0, 0.01)))
    #smp.plotting.plot((V_norm.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), (phi, 0, 3)))
    #smp.plotting.plot((dV_dphi.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), (phi, 0, 3)))

    global dV_dphi_num
    dV_dphi_num = smp.utilities.lambdify(phi, dV_dphi.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), "numpy")
    
    #Bisection
    phi_0 = []
    
    print("#")
    #print(np.pi*T_init/e_init)

    phi_max = np.pi*T_init/e_init
    #print(phi_max)

    end = 1000
    t_max = t_max_find(0.0001, end, 10000000, phi_max)*end*10
    #print(t_max)

    r = np.linspace(t_min, t_max, N)

    bisection(0.001, phi_max, t_min, t_max, r, phi_0)

    #Solving differential equation with calcuated bisection limits

    print(len(phi_0))
    #print(phi_0)

    print(phi_0[-1][1])

    r = np.linspace(t_min, t_max, N)
    
    S = (phi_0[-1][1], v_0)
    sol_ivp = solve_ivp(dS_dr, t_span = (t_min, t_max), y0 = S, t_eval = r)

    #plt.plot(r, sol_ivp.y[0], label = "sol_ivp")
    #plt.show()

    #print()
    return(phi_0[-1][1])

#####

#Definig potential
#lamda from 10^-5 to 10^-1 -> e from 0.46 to 0.02 (around)
#T from 0.01 to 100 (or T critical), it is best to "move" logarithmically

T_init = 1
e_init = 0.3

a_b_2 = 16*np.pi**2*np.exp(3/2-2*np.euler_gamma)

#Sympy
x, M, e, g, m_V, a_b, pi_V, T, lam, phi, pi = smp.symbols('x M e g m_V a_b pi_V T lambda phi pi', real = True)

lam = -e**4
M = e*phi
mu = pi*T
#a_b = 16*pi**2*smp.exp(3/2-2*smp.S.EulerGamma)
g = e
m_V = g**2*phi**2
pi_V = g**2*T**2/6
a = -20

V_tree = lam/4*phi**4
V_cw = 3/(64*pi**2)*M**4 * ( smp.log((M**2)/mu**2) - 5/6 )
V_T = 3*T**4/(2*pi**2)*( - pi**4/45 + pi**2/12*M**2/T**2 - pi/6*(M/T)**3 - 1/32*M**4/T**4*smp.log((M**2)/(a_b*T**2)) )
V_daisy = 3*T/(12*pi)*(m_V**3 - (m_V**2 + pi_V)**(3/2))

V = V_tree + V_cw + V_T + V_daisy

#smp.plotting.plot((V_tree.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), (phi, 0, 0.1)), (V_cw.subs({e:e_init, a_b:a_b_2, pi:np.pi, T:T_init}), (phi, 0, 0.1)))

#Initial conditions
v_0 = 0
N = 1000
t_min = 0.0001


escape_point = []
range_T = np.logspace(-2, 2.5, base = 10, num = 100)
range_e = np.linspace(0.1, 0.9, 9)
i = 1

#Calling functions for given T and e
for e_init in range_e:
    escape_point_temp = []
    
    for T_init in range_T:
        escape_point_temp.append(fun(T_init, e_init))
        
        print(i)
        i = i+1
    
    escape_point.append(escape_point_temp)

#print(escape_point)

for i in range(len(range_e)):
    plt.plot(range_T, escape_point[i], linestyle = '--', marker = 'o', ms = 3, label = "e = " + str(round(range_e[i], 2)))

plt.xlabel("$T$")
plt.ylabel("$\phi_0$ (Escape point)")
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize = "7")
plt.show()
