from typing import List
import matplotlib.pyplot as plt

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

#Plot
plt.figure()

plt.plot(T, eps_40, 'b', label='L=40')
plt.plot(T, eps_60, 'r', label='L=60')
plt.plot(T_, eps_80, 'y', label='L=80')
plt.plot(T, eps_100, 'g', label='L=100')
plt.title('Energy per spin at different temperatures for diffferent lattice sizes')
plt.xlabel(r'Temperature [J/$k_b$]')
plt.ylabel(r'Energy per spin ($\epsilon$) [J]')
plt.legend()
plt.grid(True)

#We save the graph in a PDF file
plt.savefig('energy_per_spin_for_diff_temperatures.pdf')

plt.figure()

plt.plot(T, m_40, 'b', label='L=40')
plt.plot(T, m_60, 'r', label='L=60')
plt.plot(T_, m_80, 'y', label='L=80')
plt.plot(T, m_100, 'g', label='L=100')
plt.title('Magnetization per spin at different temperatures for diffferent lattice sizes')
plt.xlabel(r'Temperature [J/$k_b$]')
plt.ylabel(r'Magnetization per spin (|m|)')
plt.legend()
plt.grid(True)

#We save the graph in a PDF file
plt.savefig('magnetization_per_spin_for_diff_temperatures.pdf')

plt.figure()

plt.plot(T, C_v_40, 'b', label='L=40')
plt.plot(T, C_v_60, 'r', label='L=60')
plt.plot(T_, C_v_80, 'y', label='L=80')
plt.plot(T, C_v_100, 'g', label='L=100')
plt.title('Heat capacity at different temperatures for diffferent lattice sizes')
plt.xlabel(r'Temperature [J/$k_b$]')
plt.ylabel(r'Heat capacity ($C_v$) [$k_b$]')
plt.legend()
plt.grid(True)

plt.savefig('Cv_for_diff_temperatures.pdf')


plt.figure()

plt.plot(T, Chi_40, 'b', label='L=40')
plt.plot(T, Chi_60, 'r', label='L=60')
plt.plot(T_, Chi_80, 'y', label='L=80')
plt.plot(T, Chi_100, 'g', label='L=100')
plt.title('Susceptability at different temperatures for diffferent lattice sizes')
plt.xlabel(r'Temperature [J/$k_b$]')
plt.ylabel(r'Susceptability ($\chi$) [J$^{-1}$]')
plt.legend()
plt.grid(True)

plt.savefig('chi_for_diff_temperatures.pdf')
