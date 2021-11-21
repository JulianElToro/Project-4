import numpy as np
import matplotlib.pyplot as plt
from typing import List

#The data gets extracted from a .txt file, in this case, the values of ϵ for each spin state

with  open("Ising_L20_T1_histogram.txt", "r") as  infile:
    
    lines = infile.readlines ()

    average_energy: List[float] = []
    
    for  line in  lines:
        vals = line.split()
        average_energy.append(float(vals [0]))
        
#Plot        

n=1000 

mu = np.mean(average_energy)

sigma = np.std(average_energy)

#Range of possible values of ϵ centered in the mean
e_values = np.linspace( mu - 3*sigma , mu + 3*sigma , n ) 

#The estimated probability function of ϵ
p_e = (1/(sigma * np.sqrt(2*np.pi))) * (np.exp(-(e_values - mu)**2 / (2 * sigma**2) ))
    
plt.figure()

#The histogram with the values centered in each bin
plt.hist(average_energy , bins = [-2.03, -2.01, -1.99 , -1.97 , -1.95 ], density = True , rwidth = 0.9 ,color = "b" , label = "Histogram of ⟨ϵ⟩")
#Plot of the probability distribution of ϵ 
plt.plot(e_values , p_e , color = "red", label = r'$p_ϵ(ϵ;T)$')
plt.title(r'Histogram of ⟨ϵ⟩ for $T=1 \ J/k_B$', fontsize=10)
plt.xlabel("⟨ϵ⟩/J" )
plt.ylabel("Probability density/J$^{-1}$")
plt.legend() #Legend gets painted
plt.grid(axis='y' , color = "black" , linewidth = 0.7 ) #Grid of horizontal lines
#plt.show()

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T1_histogram.pdf")