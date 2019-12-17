#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const double pi = 3.14159;
const int qq=10;
const int m=100, n=100;
int nn;

// calculate energy for spin grid with extra edges
double e_nergy(int *S[], double* cos_theta[], int m, int n, double J){
    double s_um=0;
    for(int i=1; i<m-1; i++){
        for(int j=1; j<n-1; j++){
	s_um=s_um+(cos_theta[S[i][j]][S[i+1][j]]+
			   cos_theta[S[i][j]][S[i-1][j]]+
			   cos_theta[S[i][j]][S[i][j+1]]+
			   cos_theta[S[i][j]][S[i][j-1]]);}}
    s_um=s_um/2*J; return s_um;}

double dE(int *S[], double* cos_theta[], int i, int j, int nn_old, double J){

	double DE =(cos_theta[S[i][j]][S[i+1][j]]
			   +cos_theta[S[i][j]][S[i-1][j]]
			   +cos_theta[S[i][j]][S[i][j+1]]
			   +cos_theta[S[i][j]][S[i][j-1]])-
			   (cos_theta[nn_old][S[i+1][j]]
			   +cos_theta[nn_old][S[i-1][j]]
			   +cos_theta[nn_old][S[i][j+1]]
			   +cos_theta[nn_old][S[i][j-1]]);
   return DE/2*J;}
    

int main(){
	//cout<<setprecision(3);
	double r_atio=0.1;
	double J=1,E1t,E2t;
	long int total_steps=90000000;
	
	double f=2.71828;
	// Introduce arrays to keep log_gE and H_E data
	int n_odes=m*n*4+1;
	double* log_gE = new double[n_odes];
	double* H_E    = new double[n_odes];
	double* EpN = new double[n_odes];
	for(int i=0;i<n_odes;i++){
		log_gE[i]=0;H_E[i]=0;EpN[i]=double(i-2*m*n);}
		
	// create the main, basic 2D array of spins without extra edges
	int* s_data = new int[m*n];
	int** s = new int*[m];
	// initialize it
	for(int i = 0; i < m; i++){
		s[i] = s_data + n*i;}
	
	// see WikiPedia for Potts model description
	// set element values such that cos(2*pi*s[i][j]/qq) = -1 for qq=2 ==> cos(pi), for qq>2 see by yourself..
	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++){
			s[i][j] = qq-1;}}
			
	//		s[0][0]=1;s[0][1]=0;
	//		s[1][0]=0;s[1][1]=1;
			
	//  This array will be holding all data including the boundaries
	int* S_data = new int[(m+2)*(n+2)];
	int** S = new int*[(m+2)];
	for(int i = 0; i < (m+2); ++i){
		S[i] = S_data + (n+2)*i;}
	for(int i=1;i<(m+1);i++){
		for(int j=1;j<(n+1);j++){
			S[i][j] = s[i-1][j-1];}}
	
	// do I really need the next line to initialize the corner elements?
	S[0][0]=0;S[m+1][n+1]=0;S[0][n+1]=0;S[m+1][0]=0;
	// global top and bottom boundaries
	for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
	// global left and right boundaries
	for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}
	// NEED TO VISUALLY CHECK THE PREVIOUS STEPS
	
	// 2D array of the precomputed cos(theta_j-theta_i) values
	// makes program a little bit faster
	double* cost_heta = new double[qq*qq];
	double** cos_theta = new double*[qq];
	// initialize it
	for(int i = 0; i < qq; ++i){
		cos_theta[i] = cost_heta + qq*i;}
	// precomputing the array
	for(int t = 0; t < qq; t++){
		for (int q = 0; q < qq; q++){
			cos_theta[q][t] = cos(2*pi*double(q-t)/qq);}} // do I need two doubles?
			
	// first configuration of spins on the grid => E1t
	E1t=e_nergy(S, cos_theta, m+2, n+2, J);
	cout<<E1t<<"\n"<<"------"<<"\n"<<"\n";
	
	int c_ount=0;
	double d_E;
	while(c_ount<=total_steps){c_ount++;
	 //while(c_ount<=0){c_ount++;E2t=0;
		for (int ii=0;ii<m*n;ii++){ // we assume that we should randomly access pretty much all spins in the array

			int i = rand()%m+1; int j = rand()%n+1; // the indices should fall within the boundaries
			int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
			int nn_old=S[i][j];S[i][j]=nn;
			// update boundaries in order to use functions e_nergy and d_E!!
			if((i==1)||(i==m)||(j==1)||(j==n)){// update the big matrix boundaries
			 // top and bottom boundaries update
			 for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
			 // left and right boundaries update
			 for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}
			//d_E=dE(S,cos_theta,i,j,nn_old,J);
			//E2t=E1t+d_E;			
			//double E2_t=e_nergy(S, cos_theta, m+2, n+2, J);
	E2t=e_nergy(S, cos_theta, m+2, n+2, J);
			//cout<<nn_old<<"\t"<<nn<<"\t"<<d_E<<"\t"<<E2t<<"\t"<<i<<"-"<<j<<"\t"<<E2_t<<"\n";

			if (exp(log_gE[int(E1t)+2*m*n]-log_gE[int(E2t)+2*m*n])>=rand()/(RAND_MAX + 1.0)){
				E1t=E2t;
				log_gE[int(E2t)+2*m*n]+=log(f);
				H_E[int(E2t)+2*m*n]+=1;
				s[i-1][j-1] = nn; // update original small matrix
				//cout<<"*"<<E2t<<"*"<<"\n";
			}else{
				log_gE[int(E1t)+2*m*n]+=log(f);
				H_E[int(E1t)+2*m*n]+=1;
				// no spin configuration update here. Reverse S[][] to s[][]
				S[i][j]=nn_old;
				// update boundaries in order to use functions e_nergy and d_E!!
				if((i==1)||(i==m)||(j==1)||(j==n)){// update the big matrix boundaries
				 // top and bottom boundaries update
				 for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
				 // left and right boundaries update
				 for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}				
				}
		}

		//ofstream gfile("log_gE_HE.txt");
		//	for(int i=0;i<n_odes;i++){
		//		gfile<<i<<"\t"<<EpN[i]/(m*n)<<"\t"<<log_gE[i]<<"\t"<<H_E[i]<<"\n";}
		//	gfile.close();break;

		int max_H=1;
		for(int i=0;i<n_odes;i++){
			if(H_E[i]>=max_H){max_H=H_E[i];}}			
		int min_H=max_H;			
		for(int i=0;i<n_odes;i++){
			if(H_E[i]!=0 && H_E[i]<=min_H){min_H=H_E[i];}}

		r_atio=double(min_H)/double(max_H);  // do I need two doubles?
	cout<<c_ount<<"\t"<<H_E[1]<<"\t"<<H_E[4*m*n]<<"\t"<<r_atio<<"\t"<<f<<"\n";		
		if ((H_E[0]+H_E[1]+H_E[2]+H_E[64]>=5)&&(H_E[4*m*n-2]+H_E[4*m*n-1]+H_E[4*m*n]>=2)&&(c_ount>=10000)&&(r_atio>=0.90)){
			cout<<c_ount<<"\t"<<E1t/(m*n)<<"\t"<<min_H<<"\t"<<max_H<<"\t"<<r_atio<<"\t"<<f<<"\n";
			ofstream gfile("log_gE_HE.txt");
			for(int i=0;i<n_odes;i++){
				gfile<<EpN[i]/(m*n)<<"\t"<<log_gE[i]<<"\t"<<H_E[i]<<"\n";}
			gfile.close();
			for(int i=0;i<n_odes;i++){H_E[i]=0;}
			f=sqrt(f);
			c_ount=0;}
				
	}
	
	
	// dynamic arrays clean up
	//************************
	delete[] EpN;
	
	delete[] s_data;
	delete[] s;
	
	delete[] S_data;
	delete[] S;
	
	delete[] cost_heta;
	delete[] cos_theta;
	
	delete[] log_gE;
	delete[] H_E;
	
	return 0;
}
