#include "Ising.hpp"


Ising::Ising(double T_in, int L_in, double N_in, mt19937 mt_in) {

    	//We assign the introduced values to the member variables

    	T_ = T_in;
    	L_ = L_in;
	N_ = N_in;
    	mt = mt_in;

}




//Method that creates a matrix (random or not) of spins filling it with 1 and -1

void Ising::create_matrix(mat& S, bool random) {

    uniform_int_distribution<int> dist_int1(0, 1);

    if (random) {

        for (int i = 0; i < L_; i++) {

            for (int j = 0; j < L_; j++) {

                S(i, j) = static_cast<double>(dist_int1(mt));

            }

        }


        S.replace(0, -1);
    }

    else {

        S = ones(L_, L_);
    }

}




//Method that flips one random spin, i.e it changes its sign

/*void Ising::flip_spin(mat& S, int& k, int& l) {

    uniform_int_distribution<int> dist_int(0, L_ - 1);

    k = (dist_int(mt));
    l = (dist_int(mt));

    S(k, l) = -S(k, l);

}*/




//Method that calculates the energy per spin of the system

double Ising::energy_spin(mat S) {

    double E = 0.0;

    for (int j = 0; j < L_; j++) {

        for (int i = 0; i < L_; i++) {

            E += -(S((L_ + (i - 1)) % L_, j) * S((L_ + i) % L_, j));

        }

    }

    for (int j = 0; j < L_; j++) {

        for (int i = 0; i < L_; i++) {

            E += -(S(i, (L_ + (j - 1)) % L_) * S(i, (L_ + j) % L_));

        }

    }


    double e = E / N_;

    return e;

}




//Method that calculates the magnetization per spin

double Ising::magnetization_spin(mat S) {

	double m = accu(S) /  N_;

    	return m;

}




//Method that calculates the Cv

double Ising::Cv(double mean_e, double mean_e2) {

    double Cv = ( N_ / (T_ * T_) ) *  ((mean_e2) - (mean_e * mean_e));

    return Cv;


}




//Method that calculates the X (susceptibility)

double Ising::X(double mean_m, double mean_m2) {

    double X = ( N_ / T_ ) * ((mean_m2) - (mean_m * mean_m));

    return X;


}




//Method that calculates the acceptance probability

double Ising::acceptance(mat S, int k, int l, double& dE) {

    double p;

    double s = S((L_ + (k - 1)) % L_, l) + S((L_ + (k + 1)) % L_, l) + S(k, (L_ + (l - 1)) % L_) + S(k, (L_ + (l + 1)) % L_);

    dE = 2.0 * S(k, l) * s;

    if(dE <= 0){

	p = 1.0;

    }

    else{

	p = exp(-dE/T_);

    }

    return p;

}




//Method that performs one loop of the Markov Chain Monte Carlo method

void Ising::MCMC(mat& S, int& k, int& l, double& q, double& dE, double& dM) {

	uniform_int_distribution<int> dist_int2(0, L_ - 1);

    	k = (dist_int2(mt));
    	l = (dist_int2(mt));

	//flip_spin(S, k, l);

	double a = acceptance(S, k, l, dE);

	uniform_real_distribution<double> dist_unif(0.0, 1.0);

	double r = dist_unif(mt);

    	if (r <= a) {

        	S(k, l) = -S(k, l);

		dM = S(k,l) * 2.0;

		q = 1;

	}


	else{

		q = 0;

	}


	//cout << "a=" << a << " " << "r=" << r << endl;

}
