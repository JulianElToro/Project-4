import numpy as np
import matplotlib.pyplot as plt
from typing import List

"""
The code opens two data files, creates three vectors for each file, and fills them with the three 
components of the position and velocity of the particle respectively
"""
with  open("Ising_L20_T1_ordered.txt", "r") as  infile:
    
    lines = infile.readlines ()

    average_energy: List[float] = []
    average_magnetization: List[float] = []
    
    for  line in  lines:
        vals = line.split()
        average_energy.append(float(vals [0]))
        average_magnetization.append(float(vals [1]))
        
#Plot
n=100        
n_MMC = np.linspace(1 , n, 100)

plt.figure()
plt.plot(n_MMC, average_energy)
plt.title("Evolution of ⟨ϵ⟩ with the number of Monte Carlo cycles", fontsize=10)
plt.ylabel("⟨ϵ⟩" )
plt.xlabel("Monte Carlo cycles")
plt.grid(True) #Grids get painted
plt.show()

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T1_average_energy_ordered.pdf")

plt.figure()
plt.plot(n_MMC, average_magnetization)
plt.title("Evolution of ⟨|m|⟩ with the number of Monte Carlo cycles", fontsize=10)
plt.ylabel("⟨|m|⟩")
plt.xlabel("Monte Carlo cycles")
plt.grid(True) #Grids get painted
plt.show()

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T1_average_magnetization_ordered.pdf")