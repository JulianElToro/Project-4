#include "Ising.hpp"


Ising::Ising(double T_in, int L_in) {

    //We assign the introduced values to the member variables

    T_ = T_in;
    L_ = L_in;

}




//Method that creates a matrix (random or not) of spins filling it with 1 and -1

void Ising::create_matrix(mat& S, bool random) {

    if (random) {

        for (int i = 0; i < L_   ; i++){

            for (int j = 0; j < L_ ; j++) {

                S(i, j) = static_cast<double>(rand() % 2);

            }

        }


        S.replace(0, -1);
    }

    else {

        S = ones(L_, L_);
    }

}




//Method that flips one random spin, i.e it changes its sign

void Ising::flip_spin(mat S, int& k, int& l) {

    k = ( rand() % L_);
    l = ( rand() % L_ );

    S(k, l) = -S(k, l);

}



/*
//Method that adds the periodic boundary conditions

void Ising::boundary_conditions(mat& S_, mat S) {

    //First we extract the vectors in the borders

    rowvec up = S.row(0);
    rowvec down = S.row(L_ - 1);
    vec left = S.col(0);
    vec right = S.col(L_ - 1);


    //Then we modify some of them for matching the size of the matrix

    up.insert_cols(0, 1);
    up.insert_cols(L_ + 1, 1);
    down.insert_cols(0, 1);
    down.insert_cols(L_ + 1, 1);


    //Now we create a new matrix S_ inserting these vectors in the exterior of S

    S_ = S;

    S_.insert_cols(L_, left);
    S_.insert_cols(0, right);
    S_.insert_rows(L_, up);
    S_.insert_rows(0, down);


}

*/


//Method that calculates the energy per spin of the system

double Ising::energy_spin(mat S) {

    double E = 0.0;

    for (int j = 0; j <= L_; j++) {

        for (int i = 0; i <= L_; i++) {

            E += -(S( ( L_ + (i-1) ) % L_, j) * S(( L_ + i ) % L_, j));

        }

    }

    for (int j = 0; j <= L_; j++) {

        for (int i = 0; i <= L_; i++) {

            E += -( S(i, (L_ + j - 1) % j)* S(i, (L_ + j) % j));

        }

    }


    double e = E / (L_ * L_);

    return e;

}




//Method that calculates the Cv

//Creo que esto deberíamos calcularlo después de sacar la media de e y e² que supongo que se hará usando Monte Carlo???

double Ising::Cv( double mean_e, double mean_e2) {

   double Cv = (mean_e2 - (mean_e) * (mean_e)) / (T_ * T_);

        return Cv;


}




//Method that calculates the X (susceptibility)

double Ising::X(double mean_m, double mean_m2) {

    double X = (mean_m2 - (mean_m) * (mean_m)) / T_;

        return X;


}




//Method that calculates the acceptance probability


double Ising::acceptance(mat S, int k, int l) {

    double p;

    vec exp_val(2);
    exp_val(0) = exp(8.0 / T_);
    exp_val(1) = exp(4.0 / T_);
   

    //double s = S(k , l) + S(k + 1, l + 2) + S(k, l + 1) + S(k + 2, l + 1);

    double s = S((L_ + (k-1))%L_, l) + S((L_+(k+1)) % L_, l ) + S(k, (L_ + (l-1)) % L_) + S(k, (L_ + (l+1)) % L_);

    if (S(k, l) == 1) {

       if (s == 4.0 || s == 3.0 || s == 0.0) {

            p = 1.0;

        }


        if (s == -3.0) {

            p = exp_val(1);

        }


        if (s == -4.0) {

            p = exp_val(0);

        }

    }



    if (S(k, l) == -1.0) {

        if (s == 4.0) {

            p = exp_val(0);

        }


        if (s == 3.0) {

            p = exp_val(1);

        }


        if (s == 0.0 || s == -3.0 || s == -4.0) {

            p = 1.0;

        }

    }

    return p;

}




//Method that performs one loop of the Markov Chain Monte Carlo method

void Ising::MCMC(mat S, int k, int l) {

    mat S0 = S;

    flip_spin(S, k , l);

    //boundary_conditions(S_, S);

    double r = randu();

    double a = acceptance(S, k, l);

    if (r > a) {

        S = S0;

    }

}