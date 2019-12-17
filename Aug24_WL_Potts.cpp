// July 23, 2016 -- Code clean up
// clear all; close all; load log_gE_HE.txt;d_ata=log_gE_HE(:,2);m_in=min(d_ata(d_ata~=0)); plot(log_gE_HE(:,1),log_gE_HE(:,2)-m_in);axis([-2 2 0 max(log_gE_HE(:,2)-m_in)])

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const double pi = 3.14159;
const int qq=4;
const int m=4, n=4;
int nn;

// find index on current energy value in the spectrum
int i_ndex(double* EpN, int m, int n, double E){
	double m_in=abs(EpN[0]-E);
	int inde_x;
	for(int i=0; i<4*m*n+1;i++){
		if(abs(EpN[i]-E)<=m_in){
			m_in=abs(EpN[i]-E);inde_x=i;}}
			return inde_x;}

// calculate energy for spin grid with extra edges
double e_nergy(int *S[], double* cos_theta[], int m, int n, double J){
    double s_um=0;
    for(int i=1; i<m-1; i++){
        for(int j=1; j<n-1; j++){
	s_um=s_um+(cos_theta[S[i][j]][S[i+1][j]]+cos_theta[S[i][j]][S[i-1][j]]
                    +cos_theta[S[i][j]][S[i][j+1]]+cos_theta[S[i][j]][S[i][j-1]]);}}
    s_um=s_um/2*J; return s_um;}

int main(){
	//cout<<setprecision(4);
	//int flag=1;
	double r_atio=0.1;
	double J=1,E1t,E2t;
	long int total_steps=9000000000;
	
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
	
	int c_ount=0;
	while(c_ount<=total_steps){c_ount++;
		for (int ii=0;ii<m*n;ii++){ // we assume that we should randomly access pretty much all spins in the array
			int i = rand()%m+1; int j = rand()%n+1; // the indices should fall within the boundaries
			int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
			
			S[i][j]=nn;
			E2t=e_nergy(S, cos_theta, m+2, n+2, J);
			
			if (exp(log_gE[int(E1t)+2*m*n]-log_gE[int(E2t)+2*m*n])>=rand()/(RAND_MAX + 1.0)){
				E1t=E2t;
				log_gE[int(E2t)+2*m*n]+=log(f);
				H_E[int(E2t)+2*m*n]+=1;
				S[i][j]=nn;s[i-1][j-1] = nn;
				
				if((i==1)||(i==m)||(j==1)||(j==n)){
					// top and bottom boundaries update
					for(int j=0;j<(n+2);j++){
						S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
					// left and right boundaries update
					for(int i=0;i<(m+2);i++){
						S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}
						
			}else{
				log_gE[int(E1t)+2*m*n]+=log(f);
				H_E[int(E1t)+2*m*n]+=1;S[i][j]=s[i-1][j-1];} // no spin configuration update here. Reverse S[][] to s[][]
		}

		int max_H=1;
		for(int i=0;i<n_odes;i++){
			if(H_E[i]>=max_H){max_H=H_E[i];}}			
		int min_H=max_H;			
		for(int i=0;i<n_odes;i++){
			if(H_E[i]!=0 && H_E[i]<=min_H){min_H=H_E[i];}}

		r_atio=double(min_H)/double(max_H);  // do I need two doubles?
			
		if ((H_E[0]>2)&&(H_E[4*m*n]>2)&&(c_ount>=2000)&&(r_atio>=0.90)){
			cout<<H_E[0]<<"*******"<<"\n";
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
