#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;
int main()
{
  int l;  double J=+1,H=0,T,k=1;
  // read values
  int m=200, n=200;
  //std::cin >> m >> n;

  // create array of spins
  double* s_data = new double[m*n];
  double** s = new double*[m];
  // initialize it
  for (int i = 0; i < m; ++i)
    s[i] = s_data + n*i;
 
  // write element values
  for(int i=0;i<m;i++){
   for(int j=0;j<n;j++){
    s[i][j] = 1.0;}}
  
  ofstream file("ising.txt");
  ofstream mfile("magnetization.txt");
  int total_steps=5100;  int therm_steps=int(0.2*total_steps);
  double dT=0.01;
  for(double T=2.1;T<=2.3;T+=dT){
   double M=0;
   for(int c_ount=1;c_ount<=total_steps;c_ount++){
 
   for (int ii=0;ii<m*n;ii++){
   int i = rand()%m; int j = rand()%n; 

   int l_eft   = i == 0 ? m-1 : i-1;
   int r_ight  = i == m-1 ? 0 : i+1;
   int b_ottom = j == 0 ? n-1 : j-1;
   int u_p =     j == n-1 ? 0 : j+1;

   double dE=2*J*s[i][j]*(s[l_eft][j] + s[r_ight][j] + s[i][b_ottom] + s[i][u_p]);
   
  if (exp(-dE/(k*T))>(rand()/(RAND_MAX + 1.0))){
	s[i][j] = -s[i][j];}
}
  double M_0=0;
  for(int i=0;i<m;i++){
   for(int j=0;j<n;j++){M_0+=s[i][j];}}
 
  if (c_ount>=therm_steps){M+=M_0;}
}
  cout<<T<<"\t"<<M/double(m*n)/double(0.8*total_steps)<<"\n"; //?? cast double to int
  file<<T<<"\t"<<M/double(m*n)/double(0.8*total_steps)<<"\n";
}
  file.close(); 
  // get rid of array
  delete[] s;
  delete[] s_data;
  return 0;
}
