
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from typing import List

"""
We extract the data from four .txt file, create five vectors for each file, and fill them with the values of the temperature, ϵ, |m|, C_v and χ respectively for different lattice sizes
"""
#For a lattice with L=40
with  open("Ising_L40_tsp(4000000).txt", "r") as  infile:

    lines = infile.readlines()

    T: List[float] = []
    eps_40: List[float] = []
    m_40: List[float] = []
    C_v_40: List[float] = []
    Chi_40: List[float] = []

    for  line in  lines:
        vals = line.split()
        T.append(float(vals[0]))
        eps_40.append(float(vals[1]))
        m_40.append(float(vals[2]))
        C_v_40.append(float(vals[3]))
        Chi_40.append(float(vals[4]))

#For a lattice with L=60
with  open("Ising_L60_tsp(4000000).txt", "r") as  infile:

    lines = infile.readlines()


    eps_60: List[float] = []
    m_60: List[float] = []
    C_v_60: List[float] = []
    Chi_60: List[float] = []

    for  line in  lines:
        vals = line.split()
        eps_60.append(float(vals[1]))
        m_60.append(float(vals[2]))
        C_v_60.append(float(vals[3]))
        Chi_60.append(float(vals[4]))

#For a lattice with L=80
with  open("Ising_L80_tsp(4000000).txt", "r") as  infile:

    lines = infile.readlines()

    T_: List[float] = []    
    eps_80: List[float] = []
    m_80: List[float] = []
    C_v_80: List[float] = []
    Chi_80: List[float] = []

    for  line in  lines:
        vals = line.split()
        T_.append(float(vals[0]))
        eps_80.append(float(vals[1]))
        m_80.append(float(vals[2]))
        C_v_80.append(float(vals[3]))
        Chi_80.append(float(vals[4]))

#For a lattice with L=100
with  open("Ising_L100_tsp(4000000).txt", "r") as  infile:

    lines = infile.readlines()

    eps_100: List[float] = []
    m_100: List[float] = []
    C_v_100: List[float] = []
    Chi_100: List[float] = []

    for  line in  lines:
        vals = line.split()
        eps_100.append(float(vals[1]))
        m_100.append(float(vals[2]))
        C_v_100.append(float(vals[3]))
        Chi_100.append(float(vals[4]))

# Temperatures where the Cv is maximum

Tc1 = T[np.argmax(C_v_40)]
Tc2 = T[np.argmax(C_v_60)]
Tc3 = T[np.argmax(C_v_80)]
Tc4 = T[np.argmax(C_v_100)]

# Temperatures where the X is maximum

Tc1_b = T[np.argmax(Chi_40)]
Tc2_b = T[np.argmax(Chi_60)]
Tc3_b = T[np.argmax(Chi_80)]
Tc4_b = T[np.argmax(Chi_100)]

# Size of the matrix

L = np.array([40,60,80,100])

Linv = 1/L

# Resullts for the critical temperature for different sizes

TcL = np.array([Tc1, Tc2 , Tc3 , Tc4])

TcL_b = np.array([Tc1_b, Tc2_b , Tc3_b , Tc4_b])

#TcL = np.random.random(4)+2

# Linear regression with Cv

r1 = stats.linregress(Linv,TcL)

plt.plot( Linv , TcL , 'o', color = 'blue')

Llinspace = np.linspace(0.01 , 0.025 , 1000)

plt.plot( Llinspace , r1.slope*Llinspace + r1.intercept, color = 'red')

plt.ylabel(r'$T_c(L) / \ J \ k_B^{-1}$')

plt.xlabel(r'$L^{-1}$')

plt.grid()

plt.text( 0.023 ,2.28, r'$R^2$ = ' + str(round(-r1.rvalue,2)) )

print("Linear regression (Cv results)")

print( 'a = '+ str(round(r1.slope,1)) +' ± ' + str(round(r1.stderr,1)) +'  J/k_B')

print( 'T_c(L = inf) = ' + str(round(r1.intercept,3))+ ' ± ' + str(round(r1.intercept_stderr,3)) + '  J/k_B')

plt.title( r'Linear regression of $T_c(L)$ and $L^{-1}$ ($C_V$ results) ')

plt.savefig("Ising_Tcinf_linear_regression_Cv.pdf")

plt.show()


# Linear regression with X

r2 = stats.linregress(Linv,TcL_b)

plt.plot( Linv , TcL_b , 'o', color = 'blue')

Llinspace = np.linspace(0.01 , 0.025 , 1000)

plt.plot( Llinspace , r2.slope*Llinspace + r2.intercept, color = 'red')

plt.ylabel(r'$T_c(L) / \ J \ k_B^{-1}$')

plt.xlabel(r'$L^{-1}$')

plt.grid()

plt.text( 0.023 ,2.36, r'$R^2$ = ' + str(round(-r2.rvalue,2)) )

print("Linear regression (X results)")

print( 'a = '+ str(round(r2.slope,1)) +' ± ' + str(round(r2.stderr,1)) +'  J/k_B')

print( 'T_c(L = inf) = ' + str(round(r2.intercept,3))+ ' ± ' + str(round(r2.intercept_stderr,3)) + '  J/k_B')

plt.title( r'Linear regression of $T_c(L)$ and $L^{-1}$ ($\chi$ results) ')

plt.savefig("Ising_Tcinf_linear_regression_X.pdf")