# Project-4

We've used this repository to generate plots and data for supporting our report about the study of the critical temperature of a ferromagnetic system. For that we have implemented the Markov chain Monte Carlo(MCMC) method to the Ising model. 

In the folder named "Include" you will find the .hpp file that you must download for executing every .cpp, which you will find in the folder "Source". In folder "Phyton" you will find every .py file used to make the plots.

Now we are going to detail a little bit each .cpp file.

First, we have the file "Ising.cpp" which define the main class that we will use later in our study. It stores the main characteristics of our system(its size and temperature, mainly) and it also contains some functions that allow as to calculate the energy and the magnetization per spin of the system, and to perform the MCMC. This file is not meant to be directly executed.

In the file "Ising_L2.cp" we use the MCMC method to approximate the expected values for the energy and magnetization per spin(and their squares) and the heat capacity and susceptibiilty of the system.

In the files "One_particle_Euler.cpp" and "One_particle_RK4.cpp" we evolve only one particle through time using the Foward Euler and the Runge-Kutta method (respectively). With them, we generate the .txt files stored in the folder One particle (in Python/txt files), changing the step-size in the ones used to study the relative error. With this .txt we generate the plots "x vs y RK" and "z vs t RK" (using the file "One particle RK4.py"); "x vs y FE" and "z vs t FE" (using the file "One particle FE.py"); "absolute_error"(using "absolute_error.py"), and "relative_error" (using "relative_error.py").

In the files "Resonance_1" and "Resonance_2" we study the case of 100 particles trapped in our device. With them we produce the plots "Resonance_1" and "Resonance_2", using the python codes "Resonance_1.py" and "Resonance_2.py" (respectively) and the .txt files stored in the folder Resonance (in Python/txt files).

And, finally using the file "Two_particles_RK4.cpp" we simulate the case of having two particles in our trap. Changing different desired parameters we generate different .txt files store in the folder Two particles (in Python/txt files) and with them, we produce the rest of the plots of the folder Plots (in Python) using the rest of the pyhton files.

For building and execute each of the cpp mentioned before (except for the Particle and the PenningTrap one, as we already said), as, for example "Two_particles_RK4.cpp", you need to have all the files from the folders "Source" and "Include" in the same folder of your computer and write:

Build: g++ Two_particles_RK4.cpp PenningTrap.cpp Particle.cpp -o Two_particles_Rk4.exe -larmadillo Run: ./Two_particles_RK4.exe
