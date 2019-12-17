// December 01, 2016 : Need to figure out. what is going on with double counting in the E_nergy function
// December 05, 2016 : Does not go beyond 0
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include <sys/time.h>

using namespace std;

const double pi = 3.14159;
const int qq=7;
const int m=32, n=32;
int nn;

// timing functions
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){ return 0;}
    return (double)time.tv_sec + (double)time.tv_usec*.000001;}

double get_cpu_time(){return (double)clock() / CLOCKS_PER_SEC;}


// fixed to the conventional treatment of the boundaries
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

    return s_um*J/2;} // need to figure out. what is going on with double counting... I divided by 2
			// seems still working

int main(int argc, char *argv[]){
	
	//  Start Timers
	double wall0 = get_wall_time();
	double cpu0  = get_cpu_time();

	int id, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	const clock_t begin_time = clock();

	cout<<setprecision(2);

	double r_atio=0.1;
	double J=-1,E1t=-3*(m*n),E2t=-3*(m*n);

	double l_ow=-2, u_pper=-0.1, l_ength=.3;
	double small_piece=(u_pper-l_ow-l_ength)/(ntasks-1);

	double E_min=l_ow+small_piece*double(id), E_max=E_min+l_ength;

	long int total_steps=90000000;
	
	double f=2.71828;
	// Introduce arrays to keep log_gE and H_E data. it covers the whole energy range
	int n_odes=m*n*4+1;
	double* log_gE = new double[n_odes];
	double* loggE_dis = new double[n_odes];

	double* log_gE_total = new double[n_odes*ntasks];
	double* H_E    = new double[n_odes];
	double* H_E_total = new double[n_odes*ntasks];
	double* H_E_cumul = new double[n_odes*ntasks];

	if(id==0){for(int i=0;i<n_odes*ntasks;i++){H_E_cumul[i]=0;}}

	double* EpN = new double[n_odes];
	for(int i=0;i<n_odes;i++){
		log_gE[i]=0;H_E[i]=0;EpN[i]=double(i-2*m*n);}


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
			// clock model
			 cos_theta[q][t] = cos(2*pi*double(q-t)/(qq-1));
			// Potts with delta function
			//(q==t)?(cos_theta[q][t]=1):(cos_theta[q][t]=0);
			// Ising
			// (q==t)?(cos_theta[q][t]=1):(cos_theta[q][t]=-1);
	}}
		
	// create the main, basic 2D array of spins
	int* s_data = new int[m*n];
	int** s = new int*[m];
	// initialize it
	for(int i = 0; i < m; i++){
		s[i] = s_data + n*i;}
	
	// see WikiPedia for Potts model description
	// set element values such that cos(2*pi*s[i][j]/qq) = -1 for qq=2 ==> cos(pi), for qq>2 see by yourself..

	// need to replace by the next while loop. it will be the proper way.
	//if(id>ntasks/2){
		// m and n should be EVEN for this scheme to work!!!!
	//	for(int j=0;j<n;j++){
	//		for(int i=0;i<m;i++){
	//			((i+j)%2==0)?(s[i][j]=0):(s[i][j]=(qq-1)/2);}}}
	//else{
	//	for(int j=0;j<n;j++){
	//		for(int i=0;i<m;i++){
	//		s[i][j] = qq-1;}}}

	//if(E_min*E_max<0){
	//	for(int j=0;j<n;j++){
	//		for(int i=0;i<m;i++){
	//		s[i][j] = rand()%qq;}}}

	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++){
			s[i][j] = qq-1;}}


	E1t=E_nergy(s, cos_theta, m, n, J)/(m*n);
	// usefull info for some inital nergy distribution schemes
     //	cout<<"**"<<id<<"\t"<<E_min<<"\t"<<E_max<<E1t<<"\t";

	while((E1t/(m*n)<E_min)||(E1t/(m*n)>E_max)){ // very inefficent. need to push it to the right boundaris
		int i = rand()%m; int j = rand()%n; // the indices should fall within the boundaries
		int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
		s[i][j] = nn;
		E1t=E_nergy(s, cos_theta, m, n, J);}

		cout<<id<<"\t["<<E_min<<" "<<E_max<<"]\t"<<E1t/(m*n)<<"\n";
		MPI_Barrier(MPI_COMM_WORLD);

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
			}else{// walker stays and do not leave the energy range where it was PROPERLY placed earlier
					log_gE[int(E1t)+2*m*n]+=log(f);
					H_E[int(E1t)+2*m*n]+=1;//}
					// no spin configuration update here. Reverse s[][]
					s[i][j]=nn_old;}
	 	}
			int i_max,i_min;
			int max_H=1;
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>=max_H){max_H=H_E[i];i_max=i;}}
			
			int min_H=max_H;			
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>1 && H_E[i]<min_H){min_H=H_E[i];i_min=i;}}
			(min_H!=max_H)?(r_atio=double(min_H)/double(max_H)):(r_atio=0);

//if(id==0){cout<<"\n"<<"\tH["<<i_max<<" "<<i_min<<"]\n";}
//if(id==0){cout<<"\n"<<"\tH["<<max_H<<" "<<min_H<<"]\n";}

//*********************************************************************************************
			if(r_atio>=0.80){

//			int local_counter[1]={c_ount};
//			MPI_Bcast(local_counter, 1 , MPI_INT, 0, MPI_COMM_WORLD);

			cout<<setprecision(2)<<"id="<<id<<"\tH["<<max_H<<" '"<<i_max<<"' "<<min_H<<" '"<<i_min
				<<"']\tE["<<E_min<<" "<<E_max<<"]\t"<<c_ount<<"\t"<<setprecision(5)<<f<<"\n";				
				// normalize data
				// search backward for the first non zero element and subtruct it
				// may cause the appearance of the NEGATIVE values for the pieces with negative slope
				int non_zero=log_gE[n_odes-1];
				for(int i=n_odes-1;i>=0;i--){
					if(log_gE[i]!=0){non_zero=log_gE[i];}}	
				for(int i=0;i<n_odes;i++){
					if(log_gE[i]!=0){log_gE[i]=log_gE[i]-non_zero;}} //still need +2 in matlab to normalize	


				MPI_Gather(H_E,n_odes,MPI_DOUBLE,H_E_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Gather(log_gE,n_odes,MPI_DOUBLE,log_gE_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				//MPI_Barrier(MPI_COMM_WORLD);


				// stitching pieces of different energy range together
//---------------------------------------------------------------------------------------------				
				if(id==0){
					ofstream gfile("log_gE_HE.txt");
					//H_E_total[0]=ntasks; // need to save number of tasks for Matlab reader
					for(int j=0;j<ntasks;j++){
						for(int i=0;i<n_odes;i++){
							H_E_cumul[i+j*n_odes]+=H_E_total[i+j*n_odes];
							gfile<<EpN[i]/(m*n)<<"\t"<<log_gE_total[i+j*n_odes]<<
									     "\t"<<H_E_cumul[i+j*n_odes]<<"\n";}}
				gfile.close();}
//---------------------------------------------------------------------------------------------
				f=sqrt(f);
				for(int i=0;i<n_odes;i++){H_E[i]=0;}

				c_ount=0;}
//*********************************************************************************************				

			if (f<=1.0001){break;}

	}
	
		if(id==0){cout<<"\nqq="<<qq<<", m="<<m<<", n="<<n
		     <<",\nl_ow="<<l_ow<<", u_pper="<<u_pper<<", l_ength="<<l_ength
		     <<",\nntasks="<<ntasks<<", f="<<f<<", r_atio="<<r_atio<<"\n\n";}
	
	
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
	delete[] H_E_cumul;

	MPI_Finalize();	
 
	//  Stop timers
    	double wall1 = get_wall_time();
    	double cpu1  = get_cpu_time();

	cout<<"id="<<id<<"\tt_wall="<<(wall1-wall0)/60<<"\tt_cpu="<<(cpu1-cpu0)/60<<"\n";	

	return 0;
}
