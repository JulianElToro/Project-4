
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

# Size of the matrix

L = np.array([40,60,80,100])

Linv = 1/L

# Resullts for the critical temperature for different sizes

#TcL = np.array[("cositas")]

TcL = np.random.random(4)+2

# Linear regression

slope, intercept, r_value, p_value, std_err = stats.linregress(Linv,TcL)

plt.plot( Linv , TcL , 'o', color = 'blue')

Llinspace = np.linspace(0.01 , 0.025 , 1000)

plt.plot( Llinspace , slope*Llinspace + intercept, color = 'red')

plt.ylabel(r'$T_c(L) / \ J \ k_B^{-1}$')

plt.xlabel(r'$L^-1$')

plt.grid()

plt.text(0.019, 2.7, '$T_c(L) - T_c(L=\infty) = aL^{-1}$')

plt.text( 0.019, 2.67, 'a = '+ str(round(slope,2)) + '  $J/k_B$')

plt.text(0.019, 2.64, r'$T_c(L = \infty) \ = $' + str(round(intercept,2)) + '  $J/k_B$')

plt.text( 0.019 , 2.61, r'$R^2$ = ' + str(round(r_value,2)) )

plt.title( r'Linear regression of $T_c(L)$ and $L^{-1}$ ')