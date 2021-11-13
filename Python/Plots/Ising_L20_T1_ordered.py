import numpy as np
import matplotlib.pyplot as plt
from typing import List

with  open("Ising_L20_T1_ordered.txt", "r") as  infile:
    
    lines = infile.readlines ()

    average_energy: List[float] = []
    average_magnetization: List[float] = []
    
    for  line in  lines:
        vals = line.split()
        average_energy.append(float(vals [0]))
        average_magnetization.append(float(vals [1]))
        
#Plot
n=1000        
n_MMC = np.linspace(1 , n, 1000)

plt.figure()
plt.plot(n_MMC, average_energy, "ob", mfc="b", mec = "b", ms=2)
plt.title("Evolution of ⟨ϵ⟩ with the number of Monte Carlo cycles", fontsize=10)
plt.ylabel("⟨ϵ⟩" )
plt.xlabel("Monte Carlo cycles")
plt.grid(True) #Grids get painted
#plt.show()

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T1_average_energy_ordered.pdf")

plt.figure()
plt.plot(n_MMC, average_magnetization, "ob", mfc="b", mec = "b", ms=2)
plt.title("Evolution of ⟨|m|⟩ with the number of Monte Carlo cycles", fontsize=10)
plt.ylabel("⟨|m|⟩")
plt.xlabel("Monte Carlo cycles")
plt.grid(True) #Grids get painted
#plt.show()

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T1_average_magnetization_ordered.pdf")
