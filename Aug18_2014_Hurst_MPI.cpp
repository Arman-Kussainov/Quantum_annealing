/*  http://www.cplusplus.com/reference/clibrary/cstdio/fgets/   */
/*  http://www.cplusplus.com/reference/clibrary/cstdio/feof/    */
/*  http://www.cplusplus.com/reference/clibrary/cstring/strstr/ */
/*  http://www.cplusplus.com/reference/cstdio/fgetc/ */

#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>            // Step #1 - Needed for file streaming.
#include <cmath>
#include <new>
#include <math.h>
#include <mpi.h>

#include <sys/time.h>

using namespace std;

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}


int main(int argc, char *argv[])
{

    //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

   FILE * pFile;
   FILE * mFile;
   
   int row_count;
   char mystring [200]; // quite arbitrary and big enough number to accomodate one line
   char mean_value [200]; // quite arbitrary and big enough number to accomodate one line   
   int chnls_nmbr=6;
   string m; // data channels string format
   string read_out; // to keep one line   
   double ch;// data from channel in numerical format
   
   int s_um=0.0;
   double m_ax, m_in;
   double Z=0.0;
   int E_count;
   double m_ean,s_td,E;
   int c_ount,point_count=0;
   
   double xy=0.0;
   double x=0.0;
   double y=0.0;
   double x2=0.0;
   double a,b;
 
   int id, ntasks, len;
   char single_channel[]="single_channel_";
   char mean_data[]="mean_data_"; 
   char E_data[]="E_data_";

   char current_id [33];
   char hostname[MPI_MAX_PROCESSOR_NAME];


  MPI_Init(&argc, &argv);                  /* Initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);  /* Get nr of tasks */
  MPI_Comm_rank(MPI_COMM_WORLD, &id);      /* Get id of this process */
  MPI_Get_processor_name(hostname, &len);


  sprintf(current_id,"%d",id);
  strncat(strncat(single_channel, current_id,1),".txt",4);
  strncat(strncat(mean_data, current_id,1),".txt",4);
  strncat(strncat(E_data, current_id,1),".txt",4);

   int clmn=id+1; // selecting a single data channel
  
    pFile = fopen ("izmi!borons!31.7.2014!10.8.2014.txt" , "r");
	ofstream SingleChannel; SingleChannel.open(single_channel);
	if (pFile == NULL) perror ("Error opening data file");
      else{while(!feof(pFile)){
	  fgets (mystring ,200, pFile);
	  if (strlen(mystring)>=4){ read_out=string(mystring);
	  row_count++;
	  read_out.erase(0,20);
	  m.assign(read_out,(clmn-1)*4,3);
	  SingleChannel<<atof(m.c_str())<<"\n";}}c_ount=row_count;}
	  fclose (pFile);SingleChannel.close();
   
   
   point_count=0;
   ofstream Efile; Efile.open(E_data);       
   for(int i=4;i<=900;i++){
   point_count++;
   int n=i*2; // how points are spreading out

   for(int pass=1;pass<=2;pass++){
   /* Calculate and save to a separate file mean vales*/
   if(pass==1){row_count=0;
   	
   pFile = fopen (single_channel, "r");
   ofstream outFile; outFile.open(mean_data);      
   if (pFile == NULL) perror ("Error opening data file");
      else{s_um=0.0;while(!feof(pFile)){
  	  		  	  						   						   
			 fgets (mystring ,200, pFile); // reads with the end of line symbol
			 //test if I'm reading the right, specific, line of text
			 //if (strlen(mystring)==44){ // the exact length of the single correct data row !!! elaborate further
			 	read_out=string(mystring);
				row_count++;


			 	// removing time stamp
			 	m.assign(read_out,0,3);
				s_um=s_um+atof(m.c_str());
				
				if(row_count%n==0 && row_count<=c_ount){
				outFile<<double(s_um)/double(n)<<"\n";s_um=0;}


				
			//	}
	  }}
   fclose (pFile);outFile.close();
   }
  
   /******************************************************************************/  
   else{E=0.0;E_count=0;row_count=0;
   pFile = fopen (single_channel, "r");
   mFile = fopen (mean_data, "r");

   if (mFile == NULL) perror ("Error opening data file");
      else{while(!feof(pFile)){

			 row_count++;

			 if (row_count%n==1 && row_count<=c_ount){fgets (mean_value ,200, mFile);
	   		 					 read_out=string(mean_value);
								 m.assign(read_out,0,8); // 8 because I need more decimal points									 
								 m_ean=atof(m.c_str());}

			 fgets (mystring ,200, pFile);
			 read_out=string(mystring);
			 m.assign(read_out,0,3);
			 ch=double(atof(m.c_str()));		

			 if (row_count%n==1 && row_count<=c_ount){m_ax=ch-m_ean; m_in=ch-m_ean;s_td=0.0;}

			 Z=Z+ch-m_ean;
			 s_td=s_td +(ch-m_ean)*(ch-m_ean);      	
			 if (Z>m_ax){m_ax=Z;};if(Z<m_in){m_in=Z;};

			 if (row_count%n==0 && row_count<=c_ount){E_count++;
    		 	                 E=E+(m_ax-m_in)/(sqrt(s_td/double(n)));
		  	                     s_td=0;Z=0;}


}}
     
   fclose (pFile);fclose(mFile);  	
   }
   

   
   }

///////  cout<<n<<"\t"<<E/double(E_count)<<"\t"<<id<<"\t"<<hostname<<"\n";
  Efile<<n<<"\t"<<E/double(E_count)<<"\n";

  // if E=C*n^H then log(E)=log(C)+H*log(n)
  // y = a + b*x
  // b=[<xy>-<x><y>]/[<x^2>-<x>^2]
  // a=<y>-b*<x>;
  xy+=log(n)*log(E/double(E_count));
  x+=log(n);y+=log(E/double(E_count));
  x2+=log(n)*log(n);
   }
   b=(xy-x*y/double(point_count))/(x2-x*x/double(point_count));
   a=(y-b*x)/double(point_count);

    //  Stop timers
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();

printf("\n %s \t %s \t %s \t %s \t %s \t %s \n","log(C)","H","id", "hostname", "    wall time, sec", "CPU time, sec");
printf("%6.3f \t %s  %6.3f \t %d \t %s \t %6.3f \t %6.3f \n",a,"   ",b,id,hostname,wall1-wall0,cpu1-cpu0);

   Efile<<a<<"\t"<<b<<"\n";
   Efile.close();    
      

    MPI_Finalize();	         /* Terminate MPI */
    if (id==0) printf("Ready\n");
    exit(0);
   }
																						 
																						 
