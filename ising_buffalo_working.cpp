// Ising Model in two dimensions

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

inline double std_rand()
{return rand() / (RAND_MAX + 1.0);}

double J = +1;                  // ferromagnetic coupling
int Lx, Ly;                     // number of spins in x and y
int N;                          // number of spins
int **s;                        // the spins
double T;                       // temperature
double H;                       // magnetic field

double w[17][3];                // Boltzmann factors

int steps = 0;                  // steps so far

void MonteCarlo_Metropolis ( ) {
     for (int ii = 0; ii < N; ii++){

     int i = int(Lx*std_rand());
     int j = int(Ly*std_rand());

     // find its neighbors using periodic boundary conditions
     int iPrev = i == 0 ? Lx-1 : i-1;
     int iNext = i == Lx-1 ? 0 : i+1;
     int jPrev = j == 0 ? Ly-1 : j-1;
     int jNext = j == Ly-1 ? 0 : j+1;

     // find sum of neighbors
     int sumNeighbors = s[iPrev][j] + s[iNext][j] + s[i][jPrev] + s[i][jNext];
     int delta_ss = 2*s[i][j]*sumNeighbors;

     // ratio of Boltzmann factors
     double ratio = w[delta_ss+8][1+s[i][j]];

     if (std_rand() < ratio) {
          s[i][j] = -s[i][j];}

     ++steps;}
}

double magnetizationPerSpin ( ) {
     int sSum = 0;
     for (int i = 0; i < Lx; i++)
     for (int j = 0; j < Ly; j++) {
          sSum += s[i][j];
     }
     return sSum / double(N);
}

double energyPerSpin ( ) {
     int sSum = 0, ssSum = 0;
     for (int i = 0; i < Lx; i++)
     for (int j = 0; j < Ly; j++) {
          sSum += s[i][j];
          int iNext = i == Lx-1 ? 0 : i+1;
          int jNext = j == Ly-1 ? 0 : j+1;
          ssSum += s[i][j]*(s[iNext][j] + s[i][jNext]);
     }
     return -(J*ssSum + H*sSum)/N;
}

int main (int argc, char *argv[]) {
double dT=0.02;
ofstream file("ising.data");
for(T=dT;T<=4;T+=dT){
     Lx=50; Ly=50; N = Lx*Ly;
     H=0;

     int MCSteps=2000;

     s = new int* [Lx];
     for (int i = 0; i < Lx; i++){s[i] = new int [Ly];}
     for (int i = 0; i < Lx; i++){
	for (int j = 0; j < Ly; j++){
//             s[i][j] = std_rand() < 0.5 ? +1 : -1;   // термализированный газ
               s[i][j] = 1;   // полностью намагниченная решетка
	}}
     for (int i = -8; i <= 8; i += 4) {
          w[i + 8][0] = exp( - (i * J + 2 * H) / T);
          w[i + 8][2] = exp( - (i * J - 2 * H) / T);
     }    steps = 0;



     int thermSteps = int(0.2 * MCSteps);
     cout << " Performing " << thermSteps << " steps to thermalize the system ..." << flush;
     for (int s = 0; s < thermSteps; s++)
          MonteCarlo_Metropolis();
//
cout<<"*"<<magnetizationPerSpin()<<"*";
//
     cout << " Done\n Performing production steps ..." << flush;
     double mAv = 0, m2Av = 0, eAv = 0, e2Av = 0;

     for (int s = 0; s < MCSteps; s++) {
          MonteCarlo_Metropolis();
          double m = magnetizationPerSpin();
          double e = energyPerSpin();
          mAv += m; m2Av += m * m;
          eAv += e; e2Av += e * e;
     }

     mAv /= MCSteps; m2Av /= MCSteps;
     eAv /= MCSteps; e2Av /= MCSteps;

     cout << " <m> = " << mAv << " +/- " << sqrt(m2Av - mAv*mAv) << endl;
     file << T <<"\t"<<mAv << endl;     
     cout << " <e> = " << eAv << " +/- " << sqrt(e2Av - eAv*eAv) << endl;
 }
     file.close();     
}
