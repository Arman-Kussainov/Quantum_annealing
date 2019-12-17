#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <iomanip>

using namespace std;

const double pi = 3.14159;
const int qq=7;
const int m=16, n=16;
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
    return s_um/2*J;}

int main(int argc, char *argv[]){
	int id, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	cout<<setprecision(5);

	double r_atio=0.1;
	double J=-1,E1t=-3*(m*n),E2t=-3*(m*n);
	double l_ow=-2,u_pper=0;
	double E_min=l_ow+double(id)*(u_pper-l_ow)/ntasks, E_max=l_ow+double(id+1)*(u_pper-l_ow)/ntasks;

	long int total_steps=90000000;
	
	double f=2.71828;
	// Introduce arrays to keep log_gE and H_E data
	int n_odes=m*n*4+1;
	double* log_gE = new double[n_odes];
	double* loggE_dis = new double[n_odes];

	double* log_gE_total = new double[n_odes*ntasks];
	double* H_E    = new double[n_odes];
	double* H_E_total = new double[n_odes*ntasks];

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

if(id>=6){
	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++){
//			s[i][j] = qq-1;
			s[i][j] = rand()%qq;
			}}}else{

	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++){
			s[i][j] = qq-1;
//			s[i][j] = rand()%qq;
			}}


}


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
			cos_theta[q][t] = cos(2*pi*double(q-t)/qq);
			(q==t)?(cos_theta[q][t]=1):(cos_theta[q][t]=0);
	}}



	while((E1t/(m*n)<E_min)||(E1t/(m*n)>E_max)){
		int i = rand()%m+1; int j = rand()%n+1; // the indices should fall within the boundaries
		int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
		int nn_old=S[i][j];S[i][j]=nn; s[i-1][j-1] = nn;
		// update boundaries in order to use functions e_nergy and d_E!!
		if((i==1)||(i==m)||(j==1)||(j==n)){// update the big matrix boundaries
		// top and bottom boundaries update
		for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
		// left and right boundaries update
		for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}
		//cout<<E_min<<"\t"<<E1t/(m*n)<<"\t"<<E_max<<"\n";
		E1t=e_nergy(S, cos_theta, m+2, n+2, J);}
		log_gE[int(E1t)+2*m*n]+=log(f);
		H_E[int(E1t)+2*m*n]+=1;


	int c_ount=0;

	while(c_ount<=total_steps){c_ount++;
		for (int ii=0;ii<3*m*n;ii++){ // we assume that we should randomly access pretty much all spins in the array

			int i = rand()%m+1; int j = rand()%n+1; // the indices should fall within the boundaries
			int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
			int nn_old=S[i][j];S[i][j]=nn;
			// update boundaries in order to use functions e_nergy and d_E!!
			if((i==1)||(i==m)||(j==1)||(j==n)){// update the big matrix boundaries
				// top and bottom boundaries update
				for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
				// left and right boundaries update
				for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}
									
			E2t=e_nergy(S, cos_theta, m+2, n+2, J);

			if((E_min<=E2t/(m*n)) && (E2t/(m*n)<E_max)){
//-----------------------------------------------------------------------------							
				if (exp(log_gE[int(E1t)+2*m*n]-log_gE[int(E2t)+2*m*n])>=rand()/(RAND_MAX + 1.0)){
					E1t=E2t;
//-----------------------------------------------------------------------------
					log_gE[int(E2t)+2*m*n]+=log(f);
					H_E[int(E2t)+2*m*n]+=1;
					s[i-1][j-1] = nn; // update original small matrix
							// boundaries have been updated earlier
//-----------------------------------------------------------------------------					
				}else{// walker stays. 
					log_gE[int(E1t)+2*m*n]+=log(f);
					H_E[int(E1t)+2*m*n]+=1;//}
					// no spin configuration update here. Reverse S[][] to s[][]
					S[i][j]=nn_old;
					// reverse boundaries back in order to use functions e_nergy and d_E!!
					if((i==1)||(i==m)||(j==1)||(j==n)){// update the big matrix boundaries
						// top and bottom boundaries update
						for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
						// left and right boundaries update
						for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}}
//-----------------------------------------------------------------------------											
			}else{// walker stays and do not leave the energy range
					log_gE[int(E1t)+2*m*n]+=log(f);
					H_E[int(E1t)+2*m*n]+=1;//}
					// no spin configuration update here. Reverse S[][] to s[][]
					S[i][j]=nn_old;
					// reverse boundaries back in order to use functions e_nergy and d_E!!
					if((i==1)||(i==m)||(j==1)||(j==n)){// update the big matrix boundaries
						// top and bottom boundaries update
						for(int j=0;j<(n+2);j++){S[0][j]=S[m][j];S[m+1][j]=S[1][j];}
						// left and right boundaries update
						for(int i=0;i<(m+2);i++){S[i][0]=S[i][n];S[i][n+1]=S[i][1];}}}
	 }
			int max_H=1;
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>=max_H){max_H=H_E[i];}}
			
			int min_H=max_H;			
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>1 && H_E[i]<min_H){min_H=H_E[i];}}
			(min_H!=max_H)?(r_atio=double(min_H)/double(max_H)):(r_atio=0);

//*********************************************************************************************
			if(r_atio>=0.70){
				cout<<id<<"\t"<<max_H<<"\t"<<min_H<<"\t"<<E_min<<"\t"<<E_max<<"\t"<<c_ount<<"\t"<<f<<"\n";				
				// normalize data
				// search backward for the first non zero element and subtruct it
				// may cause the appearance of the NEGATIVE values
				int non_zero=log_gE[n_odes-1];
				for(int i=n_odes-1;i>=0;i--){
					if(log_gE[i]!=0){non_zero=log_gE[i];}}	
				for(int i=0;i<n_odes;i++){
					if(log_gE[i]!=0){log_gE[i]=log_gE[i]-non_zero;}}	


				MPI_Gather(H_E,n_odes,MPI_DOUBLE,H_E_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Gather(log_gE,n_odes,MPI_DOUBLE,log_gE_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);

				// stitching pieces of energy range together
//---------------------------------------------------------------------------------------------				
				if(id==0){
					
					for(int i=0;i<n_odes;i++){loggE_dis[i]=0;} // to keep the stitched data

					double s_hift=0;
					for(int j=0;j<ntasks;j++){
						double task_shift=0;
						for(int i=0;i<n_odes;i++){
							if(log_gE_total[i+j*n_odes]!=0)
								{task_shift=log_gE_total[i+j*n_odes];}

							H_E[i]+=   H_E_total[i+j*n_odes];
							if(log_gE_total[i+j*n_odes]!=0){
								loggE_dis[i]+=(log_gE_total[i+j*n_odes]+s_hift);}}
							s_hift=s_hift+task_shift;}

					ofstream gfile("log_gE_HE.txt");
					for(int i=0;i<n_odes;i++){
						gfile<<EpN[i]/(m*n)<<"\t"<<loggE_dis[i]<<"\t"<<H_E[i]<<"\n";}
					gfile.close();}
//---------------------------------------------------------------------------------------------
				f=sqrt(f);
				for(int i=0;i<n_odes;i++){H_E[i]=0;}
				c_ount=0;}
//*********************************************************************************************				

			if (f<=1.0001){break;}


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
	delete[] loggE_dis;
	delete[] H_E;
	delete[] H_E_total;
	delete[] log_gE_total;

	MPI_Finalize();		
	return 0;
}
