# Project-4

We've used this repository to generate plots and data for supporting our report about the study of the critical temperature of a ferromagnetic system. For that we have implemented the Markov chain Monte Carlo(MCMC) method to the Ising model. 

In the folder named "Include" you will find the .hpp file that you must download for executing every .cpp, which you will find in the folder "Source". In folder "Phyton" you will find every .py file used to make the plots.

Now we are going to detail a little bit each .cpp file.

First, we have the file "Ising.cpp" which define the main class that we will use later in our study. It stores the main characteristics of our system(its size and temperature, mainly) and it also contains some functions that allow as to calculate the energy and the magnetization per spin of the system, and to perform the MCMC. This file is not meant to be directly executed.

In the file "Ising_L2.cp" we use the MCMC method to approximate the expected values for the energy and magnetization per spin(and their squares) and the heat capacity and susceptibilty of the system. The generated .txt file is called "Ising_L2_means_Cv_X.txt" and you can find it in the folder Python/txt (as all the rest .txt files)

In the file "Ising_L20.cpp" we also approximate the expected values for the energy and magnetization per spin of our system, but, this time, in every step. With it, we have generated the .txt files that end in "ordered" or "unordered". With these txt we generate the plots also ended in "ordered" or "unordered" (using the python files ended in the same things).

In the file "Ising_L20_histogram.cpp" we produce the files ended in "histogram" that contain the value of the energy in each step of the MCMC (starting from the burn-in time). With these .txt we generate the plots ended in "histogram" using the python files ended in the same thing.

For building and execute each of the last three cpp mentioned before, as, for example "Ising_L2.cpp", you need to have all the files from the folders "Source" and "Include" in the same folder of your computer and write:

Build: g++ Ising_L2.cpp Ising.cpp -o Ising_L2.exe -larmadillo 
Run: ./Ising_L2.exe

And, finally using the files "Ising_temperature_study_noparallel.cpp" or "Ising_temperature_study_parallel.cpp" we approximate the expected values for the energy and magnetization per spin and the heat capacity and susceptibilty of the system for a range of temperatures, producing the txt ended in "ts" and the plots ended with "for_diff_temperatues" using the python file ended in "tsp_diffL". 

The first one of those is compiled an executed as the other thre mentioned before, but, if you want to execute the "parallel" one, you have to do like this:

Build: g++ Ising_temperature_study_parallel.cpp Ising.cpp -fopenmp -o Ising_temperature_study_parallel.exe -larmadillo
Then(before running): export OMP_NUM_THREADS=30 (you can put as many threads as you want until 100, and it will execute faster)
Run: ./Ising_study_parallel.exe
