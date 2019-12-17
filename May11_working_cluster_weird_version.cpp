//mpicxx \-o isi_AK Apr08_test.cpp
//mpirun -np 1 isi_AK
// if t>12 add clusster algorithm///
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include <sys/time.h>

using namespace std;

const double pi = 3.14159;
const int qq=2;
const int m=16, n=16;
int nn, u_p, d_own, l_eft, r_ight;

void four_pattern(int i, int j, int m, int n, int *s[], int *visited_added[], int neighbour[2][4], int v_alue){
	
	for(int ii=0;ii<2;ii++){ // make equal to the default value 9999
		for(int jj=0;jj<4;jj++)neighbour[ii][jj]=9999;}
	
	int count_1=0;
	for(int i1=i-1;i1<=i+1;i1++){
		for(int j1=j-1;j1<=j+1;j1++){
			if ((abs(j1+i1-j-i)==1) && (i1>=0) && (i1<m) && (j1>=0) && (j1<n)){
					int f_ound=0;
					for(int s_can=0;s_can<m*n;s_can++){
						if((visited_added[0][s_can]==i1)&&(visited_added[1][s_can]==j1)){f_ound=1;break;}}
						
					if((s[i1][j1]==v_alue)&&(f_ound!=1)){
						neighbour[0][count_1]=i1;
						neighbour[1][count_1]=j1;
						count_1++;}
			}
		}
	}
}

void get_cluster(int i, int j, int m, int n, int *s[], int *visited_added[], int neighbour[2][4]){
int v_alue=s[i][j];

for(int jj=0;jj<m*n;jj++){
	for(int ii=0;ii<2;ii++){visited_added[ii][jj] = 9999;}} 

visited_added[0][0]=i;
visited_added[1][0]=j;

int next_start=0;
int next_end=0;
int h_its=0;

while(next_start<=next_end){
	int s_hift=0;
	for (int l=next_start;l<=next_end;l++){
	 
		int ic=visited_added[0][l];
		int jc=visited_added[1][l];

		four_pattern(ic, jc, m, n, s, visited_added, neighbour, v_alue);
		
		for(int io=0; io<=3;io++){
			if((neighbour[0][io]!=9999)&&((rand()/(RAND_MAX + 1.0))>0.50)){
				h_its=h_its+1;
                s_hift=s_hift+1;

				visited_added[0][h_its]=neighbour[0][io];
				visited_added[1][h_its]=neighbour[1][io];
				//s[visited_added[0][h_its]][visited_added[1][h_its]]=8;
			}
		}
	}
	
	next_start=next_end+1;
    next_end=next_start+s_hift-1;
}
    
}

// timing functions
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){ return 0;}
    return (double)time.tv_sec + (double)time.tv_usec*.000001;}

double get_cpu_time(){return (double)clock() / CLOCKS_PER_SEC;}

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

double dE(int *s[], double* cos_theta[], int m, int n, int i, int j, int nn ,double J){
    double d_e=0;
    int up,down,left,right;

	i==0   ? up=m-1  : up=i-1;
	i==m-1 ? down=0  : down=i+1;
	j==0   ? left=n-1: left=j-1;
	j==n-1 ? right=0 : right=j+1;

	d_e=(cos_theta[nn][s[down][j]]+
	     cos_theta[nn][s[up][j]]+
	     cos_theta[nn][s[i][right]]+
	     cos_theta[nn][s[i][left]])-

	    (cos_theta[s[i][j]][s[down][j]]+
	     cos_theta[s[i][j]][s[up][j]]+
	     cos_theta[s[i][j]][s[i][right]]+
	     cos_theta[s[i][j]][s[i][left]]);

    return d_e*J;} // in dE I've removed division by 2...because we calcuating dE for a single spin

int main(int argc, char *argv[]){

	int id, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	//  Start Timers
	double wall0 = get_wall_time();
	double cpu0  = get_cpu_time();

	cout<<setprecision(3);

	double r_atio=0.1;
	double J=-1;
	double E1t=-3*(m*n),E2t=-3*(m*n); // Just the arbitrary numbers to start

	//double l_ow=-2, u_pper=-0.1, l_ength=1.2*((u_pper-l_ow)/ntasks); // 1.01 just to give a slight overlap
	//double small_piece=(u_pper-l_ow-l_ength)/(ntasks-1);
	//double E_min=l_ow+small_piece*double(id), E_max=E_min+l_ength;


	double E_min=-2, E_max=-.1;


	long int total_steps=90000000;
	
	double f=2.71828;
	// Introduce arrays to keep log_gE and H_E data. it covers the whole energy range
	int n_odes=m*n*4+1;
	double* log_gE = new double[n_odes];

	double* log_gE_total = new double[n_odes*ntasks];
	double* H_E    = new double[n_odes];
	double* H_E_copy    = new double[n_odes];
	double* H_E_total = new double[n_odes*ntasks];

	double* EpN = new double[n_odes];
	for(int i=0;i<n_odes;i++){
		log_gE[i]=0;H_E[i]=0;EpN[i]=double(i-2*m*n);}

//*******************************************************************************************************
	// 2D array of the precomputed cos(theta_j-theta_i) values to makes program a little bit faster
	double* cost_heta = new double[qq*qq];
	double** cos_theta = new double*[qq];
	// initialize it
	for(int i = 0; i < qq; ++i){
		cos_theta[i] = cost_heta + qq*i;}
	// precomputing the array
	for(int t = 0; t < qq; t++){
		for (int q = 0; q < qq; q++){

			// clock model
			// cos_theta[q][t] = cos(2*pi*double(q-t)/(qq-1));

			// standard Pott with delta function
			// does not work properly...in the begining and for qq=2
			// need to check how it's interact with ising...
			//  (q==t)?(cos_theta[q][t]=1):(cos_theta[q][t]=0); 

			// Ising
			 (q==t)?(cos_theta[q][t]=1):(cos_theta[q][t]=-1);
	}}
		
	// create the main, basic 2D array of spins
	int* s_data = new int[m*n];
	int** s = new int*[m];
	// initialize it
	for(int i = 0; i < m; i++){
		s[i] = s_data + n*i;}
	
	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++){
			s[i][j] = qq-1;}} // filling array with MAX value

	// create the array of neighbouring points
	int neighbour[2][4];
	for(int ii=0;ii<2;ii++){
		for(int jj=0;jj<4;jj++)neighbour[ii][jj]=9999;}

	// create the cluster
	int* a_dded = new int[2*m*n];
	int** visited_added = new int*[2];
	// initialize it
	for(int i = 0; i < 2; i++){
		visited_added[i] = a_dded + m*n*i;}
	
	for(int j=0;j<m*n;j++){
		for(int i=0;i<2;i++){
			visited_added[i][j] = 9999;}} 
//*******************************************************************************************************

	E1t=E_nergy(s, cos_theta, m, n, J)/(m*n);
//	cout<<"**"<<id<<"\t"<<E_min<<"\t"<<E_max<<E1t<<"\t";

	while((E1t/(m*n)<E_min)||(E1t/(m*n)>E_max)){ // very inefficent. need to push it to the right boundaris
		int i = rand()%m; int j = rand()%n; // the indices should fall within the boundaries
		int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
		s[i][j] = nn;
		E1t=E_nergy(s, cos_theta, m, n, J);}


		cout<<id<<"\t"<<E_min<<"\t"<<E_max<<"\t"<<E1t/(m*n)<<"\n";

		log_gE[int(E1t)+2*m*n]+=log(f);
		H_E[int(E1t)+2*m*n]+=1;

	int c_ount=0;
	
// TEST DELETE LATER
//	srand (time(NULL));
//	for(int j=0;j<n;j++){
//		for(int i=0;i<m;i++){
//			s[i][j] = int(round(0.2+rand()/(RAND_MAX + 1.0)));}} // filling array with MAX value
// (((((((((((((((((
	

//	for(int i=0;i<m;i++){
//		for(int j=0;j<n;j++){			
//			cout<<s[i][j]<<" ";}cout<<"\n";}
//	cout<<"\n";
	
//		get_cluster(0, 0, m, n, s, visited_added, neighbour);	

//	for(int j=0;j<m*n;j++){
//		if(visited_added[0][j]!=9999){
//			//cout<<visited_added[0][j]<<" "<<visited_added[1][j];cout<<"\n";
//			s[visited_added[0][j]][visited_added[1][j]]=8;}} 
//		cout<<"\n";		
			
//	for(int i=0;i<m;i++){
//		for(int j=0;j<n;j++){			
//			cout<<s[i][j]<<" ";}cout<<"\n";}
//	cout<<"\n";
	
	
	srand (time(NULL)+id);
	while(c_ount<=total_steps){c_ount++;

		
		for (int ii=0;ii<m*n;ii++){ // we assume that we should randomly access pretty much all spins in the array
			int i = rand()%m; int j = rand()%n; // the indices should fall within the boundaries
			int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
			
			// should stay here. before I've change the matrix s
	//		float d_e=dE(s, cos_theta, m, n, i,j,nn,J);
		//	int nn_old=s[i][j];s[i][j]=nn;
		//  E2t=E1t+d_e;


///////////////// Cluster part
			get_cluster(i, j, m, n, s, visited_added, neighbour);	

for(int jl=0;jl<m*n;jl++){
	if(visited_added[0][jl]!=9999){
		int ik=visited_added[0][jl];int jk=visited_added[1][jl];
		float d_e=dE(s, cos_theta, m, n, ik,jk,nn,J);
		int nn_old=s[ik][jk];
		s[ik][jk]=nn; // veeery strange cluster algorithm 00.08 AM
		E2t=E1t+d_e;
				
		
		
		
		//}else{break;}} 

/////////////////



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
					s[ik][jk]=nn_old;}
//-----------------------------------------------------------------------------											
			}else{// walker stays and do not leave the energy range where it was PROPERLY placed earlier
					log_gE[int(E1t)+2*m*n]+=log(f);
					H_E[int(E1t)+2*m*n]+=1;//}
					// no spin configuration update here. Reverse s[][]
					s[ik][jk]=nn_old;}
	 	
	 			}else{break;}} 
	 	
	 	
	 	}


			int max_H=1;
			for(int i=0;i<n_odes;i++){	
				if(H_E[i]>=max_H){max_H=H_E[i];}}

			int min_H=max_H;	
			for(int i=0;i<n_odes;i++){
				if(H_E[i]>1 && H_E[i]<min_H){min_H=H_E[i];}}
			(min_H!=max_H)?(r_atio=double(min_H)/double(max_H)):(r_atio=0);
				
//*********************************************************************************************
//cout<<id<<" "<<r_atio<<"\n";

			if(r_atio>=0.95){
				//  Stop timers
			    	double wall1 = get_wall_time();
			    	double cpu1  = get_cpu_time();

			
				if(wall1-wall0>=60){
				cout<<setprecision(3)<<id<<"\t["<<min_H<<":"<<max_H<<"]\t["<<E_min<<":"<<E_max<<"]\t"<<c_ount<<"\t"
				<<setprecision(5)<<f<<setprecision(3)<<"\t"<<(wall1-wall0)/60<<"m\t"<<(cpu1-cpu0)/60<<"m\n";}
				else{
				cout<<setprecision(3)<<id<<"\t["<<min_H<<":"<<max_H<<"]\t["<<E_min<<":"<<E_max<<"]\t"<<c_ount<<"\t"
				<<setprecision(5)<<f<<setprecision(3)<<"\t"<<(wall1-wall0)<<"s\t"<<(cpu1-cpu0)<<"s\n";}			

				//  Start Timers
//				double wall0 = get_wall_time();
	//			double cpu0  = get_cpu_time();

				f=sqrt(f);
				for(int i=0;i<n_odes;i++){H_E_copy[i]=H_E[i];H_E[i]=0;}

				c_ount=0;
			}
//*********************************************************************************************				

			if (f<=1.001){break;}
	}

				MPI_Gather(H_E_copy,n_odes,MPI_DOUBLE,H_E_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Gather(log_gE,n_odes,MPI_DOUBLE,log_gE_total,n_odes,MPI_DOUBLE,0,MPI_COMM_WORLD);
				// stitching pieces of different energy range together
//---------------------------------------------------------------------------------------------				
				if(id==0){
					ofstream gfile("log_gE_HE.txt");
					H_E_total[0]=ntasks; // need to save number of tasks for Matlab reader
					for(int j=0;j<ntasks;j++){
						for(int i=0;i<n_odes;i++){
							gfile<<EpN[i]/(m*n)<<"\t"<<log_gE_total[i+j*n_odes]<<
									     "\t"<<H_E_total[i+j*n_odes]<<"\n";}}
				gfile.close();}
//---------------------------------------------------------------------------------------------
	
		
	
	// dynamic arrays clean up
	//************************
	delete[] EpN;
	
	delete[] s_data;
	delete[] s;
	
	delete[] cost_heta;
	delete[] cos_theta;
	
	delete[] log_gE;

	delete[] H_E;
	delete[] H_E_total;
	delete[] H_E_copy;
	delete[] log_gE_total;

	delete[] a_dded;
	delete[] visited_added; 

	//MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Finalize();	
 

	return 0;
}
