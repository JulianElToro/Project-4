#include "Ising.hpp"

#include "Ising.hpp"

int main() {

	ofstream ofile1;
	ofile1.open("Ising_L20_unordered.txt");
	ofile1 << scientific;

	int L = 20;
	double T = 1;
	random_device rd;
	mt19937 mt(rd());

	Ising my_system(T, L, mt);

	mat S(L, L);
	int k, l;

	my_system.create_matrix(S, true);

	/*
	S.print("S=");
	my_system.flip_spin(S, k, l);
	S.print("New S=");
	cout << "Position of the fliped spin =" << k << l << endl;
	double e = my_system.energy_spin(S);
	cout << "e =" << e << endl;
	cout << my_system.magnetization_spin(S) << endl;
	*/

	int MC_cycles = 100;

	double e, m, e_mean, m_mean;

	for (int i = 0; i < MC_cycles; i++) {

		my_system.MCMC(S, k, l);

		e =+ my_system.energy_spin(S);
		e_mean = e / (1.0*i+1);

		m =+ abs(my_system.magnetization_spin(S));
		m_mean = m / (1.0*i+1);

		ofile1 << e_mean << m_mean << endl;

	}

	ofile1.close();

	return 0;

}
