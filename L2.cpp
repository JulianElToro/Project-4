#include "Ising.hpp"

int main(){

	int L = 2;
	double T = 1;
	random_device rd;
	mt19937 mt(rd());

	Ising my_system(T, L, mt);

	mat S(L,L);
	int k,l;

	my_system.create_matrix(S, true);

	S.print("S=");

	my_system.flip_spin(S, k, l);

	S.print("New S=");

	cout << "Position of the fliped spin =" << k << l << endl;

	double e = my_system.energy_spin(S);

	cout << "e =" << e << endl;

	return 0;

}
