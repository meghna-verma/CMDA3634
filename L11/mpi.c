#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
 
  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

   
  //need running tallies
  long long int Ntotal;
  long long int Ncircle;
  long long int ntotal;
  long long int ncircle;

  printf("Rank of each MPI process is %d \n", rank);
  printf("Size of the MPI processes is %d \n", size);

  //seed random number generator
  //double seed = 1.0;
  long long int n;
  double seed = rank;
  srand48(seed);
   
  for (n=0; n<1000000000;n++) {
    //gererate two random numbers
    double rand1 = drand48(); //drand48 returns a number between 0 and 1
    double rand2 = drand48();
    
    double x = -1 + 2*rand1; //shift to [-1,1]
    double y = -1 + 2*rand2;

    //check if its in the circle
    if (sqrt(x*x+y*y)<=1) Ncircle++;
    Ntotal++;
}
 MPI_Allreduce(&Ncircle,
	       &ncircle,
	       1,
	       MPI_FLOAT,
               MPI_SUM,
	       MPI_COMM_WORLD);

 MPI_Allreduce(&Ntotal,
              &ntotal,
              1,               // double pi = Ncircle/ (double) Ntotal;
              MPI_FLOAT,
              MPI_SUM,  
	      MPI_COMM_WORLD);  

double pi = 4.0*ncircle/ (double) ntotal;     

int r;

for (r=0; r<size; r++) {
   if (r==rank) {
     printf("Rank %d has value %f after the reduction. \n", rank, pi);
      }
    MPI_Barrier(MPI_COMM_WORLD);
  }

//printf("Our estimate of pi is %f \n", pi);
     MPI_Finalize();
     return 0;
}
