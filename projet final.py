import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


#Number of point and interval of time on what we base our opperation
time_s = np.linspace(0, 60, 1000)

#Effective delayed fraction of precursor group i, (no dimension)
beta_i = np.array([0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00027 
])

#Sum of all 6 precursor group's effective delayed fraction
beta_total = sum(beta_i)

#Decay constant of precursor group i, s-1
lambda_i = np.array([0.01244, 0.03051, 0.11144, 0.30137, 1.13631, 3.01368
])

#Mean genration time, seconde
LAMBDA = 4.3*10**-5

#Reactivity : 0.5$ step insertion (Ganapol benchmark case), (no dimension)
rhô = 0.5 * beta_total

#Power initial value, standardiezed at 1 (arbitrary value)
P_0 = 1

#Initial precursor concentration of precursor group i
Ci_0 = (beta_i/LAMBDA) * P_0/lambda_i



def integration(time_s, y) :
    
    (C1, C2, C3, C4, C5, C6, P
     ) = y
    
    Ci = y[:6]
    
    dydt = [beta_i[0]/LAMBDA * P - lambda_i[0] * C1, beta_i[1]/LAMBDA * P - lambda_i[1] * C2,
     beta_i[2]/LAMBDA * P - lambda_i[2] * C3, beta_i[3]/LAMBDA * P - lambda_i[3] * C4,
     beta_i[4]/LAMBDA * P - lambda_i[4] * C5, beta_i[5]/LAMBDA * P - lambda_i[5] * C6,
     (rhô - beta_total) / LAMBDA * P + sum(lambda_i * Ci)

]
    return dydt

#We regroup all the initial values to use them for solve_ivp
initial_values = np.append(Ci_0, P_0)

solution = solve_ivp(integration, 
                     (time_s.min(), time_s.max()),
                     initial_values,
                     method='Radau',
                     t_eval=time_s,
                     
)

#plt.style.use('dark_background')
#plt.plot(time_s, solution.y[6, 0:])
#plt.title('Power')
#plt.yscale('log')
#plt.show()

#Production yield of Iodine-135 by fission, (no dimension)
gamma_Iodine = 0.063

#Direct production yield of Xenon-135 by fission, (no dimension)
gamma_Xenon = 0.0024

#Probability of desintegration of Iodine-135, s-1
lambda_Iodine = 2.90 * 10 ** -5

#Probability of desintegration of Xenon-135, s-1
lambda_Xenon = 2.10 * 10 ** -5

#Probability of Xenon-135 to absorb a neutron, cm-2
sigma_a_Xenon = 2.6 * 10 ** -18

#Initial concentration of Iodine-135, atoms by cm3
#I_0 = (gamma_Iodine * P_0) / lambda_Iodine
I_0 = 0

#Initial concentration of Xenon-135, atoms by cm3
#Xe_0 = (gamma_Xenon * P_0 + lambda_Iodine * I_0) / (lambda_Xenon + sigma_a_Xenon * P_0)
Xe_0 = 0

time_s2 = np.linspace(0, 200000, 1000)

# Fission rate density, fissions cm-3 s-1 (~100 kW/L, typical PWR core average)
R_fission = 3.0e12
# Thermal neutron flux, n cm-2 s-1
phi = 3.0e13

def integration_bloc2(t, z):
    I, Xe = z
    dI = gamma_Iodine * R_fission - lambda_Iodine * I
    dXe = (gamma_Xenon * R_fission + lambda_Iodine * I
           - lambda_Xenon * Xe - sigma_a_Xenon * phi * Xe)
    return [dI, dXe]

initials_values = np.append(I_0, Xe_0)

solutions = solve_ivp(integration_bloc2,
                      (time_s2.min(), time_s2.max()),
                      initials_values,
                      method='Radau',
                      t_eval=time_s2
                      
 )

plt.style.use('dark_background') 
plt.plot(time_s2 / 3600, solutions.y[0], label="I-135")
plt.plot(time_s2 / 3600, solutions.y[1], label="Xe-135")
plt.title("Iodine-135 and xenon-135 buildup at constant power")
plt.xlabel("Time (h)")
plt.ylabel("Concentration (atoms/cm³)")
plt.legend()      
plt.show() 
print(solutions.y[0, -1], solutions.y[1, -1])
                      
                      
                      
                      
                      
                      
                      
                      
                     











