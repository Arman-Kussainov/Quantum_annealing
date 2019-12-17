// December 01, 2016 : Need to figure out. what is going on with double counting in the E_nergy function
// December 05, 2016 : Does not go beyond 0
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <iomanip>

using namespace std;

const double pi = 3.14159;
const int qq=13;
const int m=8, n=8;
int nn;

// fixed to conventional treatment of the boundaries
double E_nergy(int *s[], double* cos_theta[], int m, int n, double J){
    double s_um=0;
    int up,down,left,right;

  for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){

	i==0   ? up=m-1  : up=i-1;
	i==m-1 ? down=0  : down=i+1;
	j==0   ? left=n-1: left=j-1;
	j==n-1 ? right=0 : right=j+1;

	s_um=s_um+(cos_theta[s[i][j]][s[down][j]]+
			   cos_theta[s[i][j]][s[up][j]]+
			   cos_theta[s[i][j]][s[i][right]]+
			   cos_theta[s[i][j]][s[i][left]]);}}

    return s_um*J;} // need to figure out. what is going on with double counting

int main(int argc, char *argv[]){
	int id, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	cout<<setprecision(5);

	double r_atio=0.1;
	double J=-1,E1t=-3*(m*n),E2t=-3*(m*n);

	double l_ow=-2, u_pper=.8, l_ength=.9;
	double small_piece=(u_pper-l_ow-l_ength)/(ntasks-1);

	double E_min=l_ow+small_piece*double(id), E_max=E_min+l_ength;

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

	// specifaclly designed for q=13
	if(id>=int(ntasks/2)+3){
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){
//			s[i][j] = qq-1;
			s[i][j] = rand()%qq;
			}}}
	else{
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){
			s[i][j] = qq-1;
//			s[i][j] = rand()%qq;
			}}}

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

		int i = rand()%m; int j = rand()%n; // the indices should fall within the boundaries
		int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
		s[i][j] = nn;
		E1t=E_nergy(s, cos_theta, m, n, J);
		// need to check what id vs.energy range is stacked at this position
}

		log_gE[int(E1t)+2*m*n]+=log(f);
		H_E[int(E1t)+2*m*n]+=1;

	int c_ount=0;

	while(c_ount<=total_steps){c_ount++;
		for (int ii=0;ii<3*m*n;ii++){ // we assume that we should randomly access pretty much all spins in the array

			int i = rand()%m; int j = rand()%n; // the indices should fall within the boundaries
			int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
			int nn_old=s[i][j];s[i][j]=nn;
									
			E2t=E_nergy(s, cos_theta, m, n, J);

			if((E_min<=E2t/(m*n)) && (E2t/(m*n)<E_max)){
//-----------------------------------------------------------------------------							
				if (exp(log_gE[int(E1t)+2*m*n]-log_gE[int(E2t)+2*m*n])>=rand()/(RAND_MAX + 1.0)){
					E1t=E2t;
					log_gE[int(E2t)+2*m*n]+=log(f);
					H_E[int(E2t)+2*m*n]+=1;
//-----------------------------------------------------------------------------					
				}else{// walker stays. 
					log_gE[int(E1t)+2*m*n]+=log(f);
					H_E[int(E1t)+2*m*n]+=1;
					// no spin configuration update here. Reverse s[][]
					s[i][j]=nn_old;}
//-----------------------------------------------------------------------------											
			}else{// walker stays and do not leave the energy range
					log_gE[int(E1t)+2*m*n]+=log(f);
					H_E[int(E1t)+2*m*n]+=1;//}
					// no spin configuration update here. Reverse s[][]
					s[i][j]=nn_old;}
	 }
			int max_H=1;
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>=max_H){max_H=H_E[i];}}
			
			int min_H=max_H;			
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>1 && H_E[i]<min_H){min_H=H_E[i];}}
			(min_H!=max_H)?(r_atio=double(min_H)/double(max_H)):(r_atio=0);

//*********************************************************************************************
			if(r_atio>=0.90){
				ofstream gfile("log_gE_HE.txt");
				cout<<id<<"\t"<<max_H<<"\t"<<min_H<<"\t"<<E_min<<"\t"<<E_max<<"\t"<<c_ount<<"\t"<<f<<"\n";				
				// normalize data
				// search backward for the first non zero element and subtruct it
				// may cause the appearance of the NEGATIVE values for the pieces with negative slope
				int non_zero=log_gE[n_odes-1];
				for(int i=n_odes-1;i>=0;i--){
					if(log_gE[i]!=0){non_zero=log_gE[i];}}	
				for(int i=0;i<n_odes;i++){
					if(log_gE[i]!=0){log_gE[i]=log_gE[i]-non_zero;}} //need +2 in matlab to normalize	


				MPI_Gather(H_E,n_odes,MPI_DOUBLE,H_E_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Gather(log_gE,n_odes,MPI_DOUBLE,log_gE_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);

				// stitching pieces of different energy range together
//---------------------------------------------------------------------------------------------				
				if(id==0){
					for(int j=0;j<ntasks;j++){
						for(int i=0;i<n_odes;i++){
							gfile<<EpN[i]/(m*n)<<"\t"<<log_gE_total[i+j*n_odes]<<
									     "\t"<<H_E_total[i+j*n_odes]<<"\n";}}}
//---------------------------------------------------------------------------------------------
				f=sqrt(f);
				for(int i=0;i<n_odes;i++){H_E[i]=0;}
				c_ount=0;
				gfile.close();}
//*********************************************************************************************				

			if (f<=1.002){break;}


	}
	
	
	// dynamic arrays clean up
	//************************
	delete[] EpN;
	
	delete[] s_data;
	delete[] s;
	
	delete[] cost_heta;
	delete[] cos_theta;
	
	delete[] log_gE;
	//delete[] loggE_dis;
	delete[] H_E;
	delete[] H_E_total;
	delete[] log_gE_total;

	MPI_Finalize();		
	return 0;
}
