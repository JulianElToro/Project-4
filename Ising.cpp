#include "Ising.hpp"


Ising::Ising(double T_in, int L_in) {

        //We assign the introduced values to the member variables

        T_ = T_in;
        L_ = L_in;

}




//Method that creates a matrix (random or not) of spins filling it with 1 and -1

void Ising::create_matrix(imat& S,bool random) {

        if (random) {

                S = randi(L_, L_, distr_param(0, 1));

		S.replace(0, -1);
        }

        else{

                S = ones(L_, L_);
        }

}




//Method that flips one random spin, i.e it changes its sign

void Ising::flip_spin(imat S, int& k, int& l){

	k = rand() % L_;
	l = rand() % L_;

	S(k, l) = - S(k, l);

}




//Method that adds the periodic boundary conditions

void Ising::boundary_conditions(imat& S_, imat S) {

        //First we extract the vectors in the borders

        rowvec up = S.row(0);
        rowvec down = S.row(L_ - 1);
        vec left = S.col(0);
        vec right = S.col(L_ - 1);


        //Then we modify some of them for matching the size of the matrix

        up.insert_rows(0, 1);
        up.insert_rows(L_ + 2, 1);
        down.insert_rows(0, 1);
        down.insert_rows(L_2 + 2, 1);


        //Now we create a new matrix S_ inserting these vectors in the exterior of S

        S_ = S;

        S_.insert.cols(0, left);
        S_.insert.cols(L_ + 2, right);
        S_.insert_rows(0, up);
        S_.insert_rows(L_ + 2, down);

}




//Method that calculates the energy per spin of the system

double Ising::energy_spin(imat S_) {

        double E = 0.0;

        for (int j = 1; j < L_ + 2; j++){

                for (int i = 0; i < L_ + 2; i++) {

                        E += -( S_(i, j) * S_(i + 1, j) );

                }

        }

        for (int j = 0; j < L_ + 2; j++) {

                for (int i = 1; i < L_ + 2; i++) {

                        E += -(S_(i, j) * S_(i, j + 1));

                }

        }


        double e = E / (L_ * L_);

        return e;

}




//Method that calculates the Cv

//Creo que esto deberíamos calcularlo después de sacar la media de e y e² que supongo que se hará usando Monte Carlo???

double Ising::Cv(imat S_, double mean_e, double mean_e2) {

        /*double E2 = 0.0;

        for (int j = 1; j < L_ + 2; j++) {

                for (int i = 0; i < L_ + 2; i++) {

                        E2 += (S_(i, j) * S_(i + 1, j))* (S_(i, j) * S_(i + 1, j));

                }

        }

        for (int j = 0; j < L_ + 2; j++) {

                for (int i = 1; i < L_ + 2; i++) {

                        E2 += (S_(i, j) * S_(i, j + 1))*(S_(i, j) * S_(i, j + 1));

                }

        }


        double e2 = E2 / (L_ * L_);

	double e2 = E2 / (L_ * L_);


        double Cv = ( e2 - energy_spin(S)*energy_spin(S) ) / (T_ * T_);*/


        double Cv = (mean_e2 - (mean_e)**2) / (T_ * T_)

        return Cv


}




//Method that calculates the X (susceptibility)

double Ising::X(double mean_m, double mean_m2){

	double X = (mean_m2 - (mean_m)**2) / T_

        return X


}




//Method that calculates the acceptance probability


double Ising::acceptance(imat S_, int k, int l){

        //double a;
        double p;

        vec exp_val(2/*5*/);
        exp_val(0) = exp(8/T_);
        exp_val(1) = exp(4/T_);
        /*exp_val(2) = 1;
        exp_val(3) = exp(-4/T_);
        exp_val(4) = exp(-8/T_);*/

        double s = S_(k+1, l) + S_(k+1, l+2) + S_(k, l+1) + S_(k+2, l+1);

	if (S_(k+1, l+1) == 1){

                /*if (s == 4){

                        p = exp_val(4);

                }


                if (s == 3){

                        p = exp_val(3);

                }


                if (s == 0){

                        p = exp_val(2);

                }*/


                if (s == 4 || s == 3 || s==0){

                        p = 1;

                }


                if (s == -3){

			p = exp_val(1);

                }


                if (s == -4){

                        p = exp_val(0);

                }

        }



        if (S(k+1, l+1) == -1){

                if (s == 4){

                        p = exp_val(0);

                }


                if (s == 3){

                        p = exp_val(1);

                }


                if (s == 0 || s == -3 || s == -4 ){

			p = 1;

                }


                /*if (s == 0){

                        p = exp_val(2);

                }


                if (s == -3){

                        p = exp_val(3);

                }


                if (s == -4){

                        p = exp_val(4);

                }*/

        }

        return p;

}




//Method that performs one loop of the Markov Chain Monte Carlo method

void Ising::MCMC(imat S, imat& S_,int k,int l){

        mat S0 = S

        flip_spin(S);

        boundary_conditions(S_,S);

        double r = randu();

        double a = acceptance(S_, k, l);

        if (r > a){

                S = S0;

        }

}
