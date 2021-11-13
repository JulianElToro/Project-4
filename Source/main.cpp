#include "Ising.hpp"

int main() {

	ofstream ofile1;
	ofile1.open("MC_test_1000_MC_cycles.txt");
	ofile1 << scientific;

	ofstream ofile2;
	ofile2.open("MC_test_1000_MC_cycles2.txt");
	ofile2 << scientific;

	int L = 2;
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

	int MC_cycles = 1000;

	double e_mean, m_mean, e2_mean, m2_mean, Cv_value, X_value;

	for (int i = 0; i < MC_cycles; i++) {

		my_system.MCMC(S, k, l);

		//S.print();

		//cout << "next" << endl;

		/*

		vector <double>  e_values, m_values;

		e_values.push_back(my_system.energy_spin(S));
		
		

		m_values.push_back(my_system.magnetization_spin(S));
		*/
		
		/*
		e_values.insert_rows( i, my_system.energy_spin(S) );

		m_values.insert_rows( i, my_system.magnetization_spin(S) );
		


		e_values.insert_rows(i, 1);

		m_values.insert_rows(i, 1);
		*/

		//double e_values = my_system.energy_spin(S);

		//double m_values = my_system.magnetization_spin(S);

		e_mean =+ my_system.energy_spin(S);

		m_mean =+ my_system.magnetization_spin(S);

		e2_mean =+ (my_system.energy_spin(S))*(my_system.energy_spin(S));

		m2_mean =+ (my_system.magnetization_spin(S))*(my_system.magnetization_spin(S));
			
		ofile1 << my_system.energy_spin(S) << my_system.magnetization_spin(S) << endl;
		
	}

	e_mean = e_mean / MC_cycles;

	m_mean = m_mean / MC_cycles;

	e2_mean = e2_mean / MC_cycles;

	m2_mean = m2_mean / MC_cycles;

	Cv_value = my_system.Cv(e_mean, e2_mean);

	X_value = my_system.X(m_mean, m2_mean);

	ofile2 << e_mean << "    " << m_mean << "    " << e2_mean << "    " << m2_mean << endl;

	ofile1.close();

	ofile2.close();

	
	return 0;

}