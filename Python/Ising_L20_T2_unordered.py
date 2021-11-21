import numpy as np
import matplotlib.pyplot as plt
from typing import List

with  open("Ising_L20_T2_unordered.txt", "r") as  infile:
    
    lines = infile.readlines ()

    average_energy: List[float] = []
    average_magnetization: List[float] = []
    
    for  line in  lines:
        vals = line.split()
        average_energy.append(float(vals [0]))
        average_magnetization.append(float(vals [1]))
        
#Plot
n=np.size(average_energy)      
n_MMC = np.linspace(1 , n ,n )

plt.figure()
plt.plot(n_MMC, average_energy, "ob", mfc="b", mec = "b", ms=2)
plt.title(r'Evolution of ⟨ϵ⟩ with the number of Monte Carlo cycles for T=2.4J/$k_B$ (unordered)', fontsize=10)
plt.ylabel("⟨ϵ⟩/$J$" )
plt.xlabel("Monte Carlo cycles")
plt.grid(True) #Grids get painted

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T2_e_unordered.pdf")

plt.figure()
plt.plot(n_MMC, average_magnetization, "ob", mfc="b", mec = "b", ms=2)
plt.title(r'Evolution of ⟨|m|⟩ with the number of Monte Carlo cycles for T=2.4J/$k_B$ (unordered)', fontsize=10)
plt.ylabel("⟨|m|⟩")
plt.xlabel("Monte Carlo cycles")
plt.grid(True) #Grids get painted

#The graph is saved in a PDF file
plt.savefig("Ising_L20_T2_m_unordered.pdf")

