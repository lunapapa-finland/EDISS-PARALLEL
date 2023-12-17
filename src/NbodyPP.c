#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

const double G  = 6.67259e-7;  /* Gravitational constant (should be e-10 but modified to get more action */
const double dt = 1.0;         /* Length of timestep */

/* Writes out positions (x,y) of N particles to the file fn 
   Returns zero if the file couldn't be opened, otherwise 1 */
int write_particles(int N, double *X, double *Y, char *fn) {
  FILE *fp;
  /* Open the file */
  if ((fp=fopen(fn, "w")) == NULL) {
    printf("Couldn't open file %s\n", fn);
    return 0;
  }
  /* Write the positions to the file fn */
  for (int i=0; i<N; i++) {
    fprintf(fp, "%3.2f %3.2f \n", X[i], Y[i] );
  }
  fprintf(fp, "\n" );
  fclose(fp);  /* Close the file */
  return(1);
}

// Distance between points with coordinates (px,py) and (qx,qy)
double dist(double px, double py, double qx, double qy) {
  return sqrt (pow(px-qx,2)+pow(py-qy,2) );
  // Could also be written as sqrt( (px-qx)*(px-qx) + (py-qy)*(py-qy) )
}


/* Parallelly Computes forces between bodies */
void ComputeForce_parallel(int first, int last,int N, double *X, double *Y, double *mass, double *Fx, double *Fy){
    const double mindist  = 0.0001;  /* Minimal distance of two bodies of being in interaction*/
    // GlobalIndex is from [first to last)
    // update local Fx and Fy
    for (int i = 0; i < last - first; i++) {      // Compute the force for all bodies
      Fx[i] = Fy[i] = 0.0;             // Initialize force vector to zero
      for ( int j = 0; j < N; j++) {   // The force on a body i is the sum of forces from all other bodies j
        // Fx[i] and Fy[i] is calculated based on !X[i] and !Y[i]
        // Thus local Fx[i] and Fy[i] => global X[i+first] and Y[i+firstj] need to exclude the Global X[j] and Y[j]
        if (i+first != j) {                  //     but not from it self
          // Distance between points i and j
          double r = dist(X[i+first], Y[i+first], X[j], Y[j]); 
          if (r>mindist) {        // Very near-distance forces are ignored
            double r3 = pow(r,3);     // Could also be written as r3=r*r*r;
            Fx[i] += G*mass[i+first]*mass[j]*(X[j]-X[i+first]) / r3;
            Fy[i] += G*mass[i+first]*mass[j]*(Y[j]-Y[i+first]) / r3;
          }
        }
      }
    }

}

int main(int argc, char * argv[]) {
    // Define Variables   
    int np, id;
    int currentsteps;
    double starttime, endtime;
    const int N=1000;                   // Number of bodies 
    const int timesteps = 1000;        // Number of timesteps
    const double size = 100.0;  
    double *mass;          /* mass of bodies */
    double *X;             /* x-positions of bodies */
    double *Y;             /* y-positions of bodies */
    double *tempX;             /* x-positions of bodies */
    double *tempY;             /* y-positions of bodies */
    double *Vx;            /* velocities on x-axis of bodies */
    double *Vy;            /* velocities on y-axis of bodies */
    double *Fx;            /* forces on x-axis of bodies */
    double *Fy;            /* forces on y-axis of bodies */
    int first, last; 

    // Initialize MPI 
    MPI_Init( & argc, & argv);
    MPI_Comm_size(MPI_COMM_WORLD, & np);
    MPI_Comm_rank(MPI_COMM_WORLD, & id);

    if(N % np != 0){
      if (id == 0){
        printf("Unable to distribute the workload evenly among %d processes\n"\
        "due to the total number of principals being 1000. \n"\
        "Please reconsider the number of processes to ensure that\n"\
        "it can be evenly divided by the number of entities. Say, 10, 20, 25, etc.\n", np); 
      }
      exit(0);
    }

    // Calculate the number of particles per process
    int particles_per_process = N / np;

    
    // Acrossed Variable: Initialize mass and position arrays in all processes
    mass = (double *)malloc(N * sizeof(double));
    X = (double *)malloc(N * sizeof(double));
    Y = (double *)malloc(N * sizeof(double));

    // Use process 0 to generate oroginal data of mass and position arrays
    if(id == 0){
      // Seed the random number generator so that it generates a fixed sequence
      unsigned short int seedval[3] = {7, 7, 7};
      seed48(seedval);
      /* Initialize mass and position of bodies */
      for (int i = 0; i<N; i++){
        mass[i]  = 1000.0*drand48();   // 0 <= mass < 1000
        X[i] = size*drand48();      // 0 <= X < 100
        Y[i] = size*drand48();      // 0 <= Y < 100
        }
      // Write intial particle coordinates to a file
      write_particles(N, X, Y, "initial_pos_PP.txt");
      // Start Timer
      starttime = MPI_Wtime();
    }

    // Synchronize before Broadcasting
    MPI_Barrier(MPI_COMM_WORLD);

    // Broadcast the initial mass and position X and Y of bodies to each process other than 0
    MPI_Bcast( mass, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( X, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( Y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Domestic Variable: calculate the first and last bodies in current process.
    first = id*particles_per_process;
    last = (id+1)*particles_per_process;

    // Domestic Varialble: allocate space for forces, velocities in each process
    Vx = (double *)malloc(particles_per_process * sizeof(double));
    Vy = (double *)malloc(particles_per_process * sizeof(double));
    Fx = (double *)malloc(particles_per_process * sizeof(double));
    Fy = (double *)malloc(particles_per_process * sizeof(double));
    tempX = (double *)malloc(particles_per_process * sizeof(double));
    tempY = (double *)malloc(particles_per_process * sizeof(double));

    // Every process get a copy of initial mass and position
    // Compute the initial forces for local particles_per_process
    ComputeForce_parallel(first, last, N, X, Y, mass, Fx, Fy);

    // local and global Index matching
    // in Each process, set up the velocity vectors caused by initial forces for Leapfrog method
    for(int i = 0; i<particles_per_process; i++){
      Vx[i] = 0.5*dt*Fx[i]/mass[i+first];
      Vy[i] = 0.5*dt*Fy[i]/mass[i+first];
    }
    /* Main loop: 
      - Move the bodies
      - Calculate forces of the bodies with their new position
      - Calculate velocities of the bodies with the new forces
      - Copy the updated positions to the old positions (for use in next timestep)
    */
    int t=0;
    while (t<timesteps) {    // Loop for this many timesteps

      t++;
      // Use Process ID 0 to print steps
      if (id == 0){
        printf("%d ", t); fflush(stdout);  // Print out the timestep
      }

      // Calculate new positions 
      for (int i=0;i<particles_per_process;i++){
        tempX[i] = X[i+first] + Vx[i]*dt;
        tempY[i] = Y[i+first] + Vy[i]*dt;
      }
      
      // MPI_Allgather can be used if comment the MPI_Gather and MPI_Bcast
      MPI_Allgather(tempX, particles_per_process, MPI_DOUBLE, X, particles_per_process, MPI_DOUBLE,
            MPI_COMM_WORLD);
      MPI_Allgather(tempY, particles_per_process, MPI_DOUBLE, Y, particles_per_process, MPI_DOUBLE,
            MPI_COMM_WORLD);

      // need a machenisum to update local X and Y to Others and sycn in different Process.
      // MPI_Gather(tempX, particles_per_process, MPI_DOUBLE, X, particles_per_process, MPI_DOUBLE, 0,
      //       MPI_COMM_WORLD);
      // MPI_Gather(tempY, particles_per_process, MPI_DOUBLE, Y, particles_per_process, MPI_DOUBLE, 0,
      //       MPI_COMM_WORLD);

      // // Bcast X and Y in Process 0
      // MPI_Bcast( X, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      // MPI_Bcast( Y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      // // Synchronize processes after updating X and Y Globally
      // MPI_Barrier(MPI_COMM_WORLD);

      // Compute the initial forces for local particles_per_process
      ComputeForce_parallel(first, last, N, X, Y, mass, Fx, Fy);
      
       /* Update velocities of bodies */ 
      for (int i=0;i<particles_per_process;i++){		
        Vx[i] = Vx[i] + Fx[i]*dt/mass[i+first];
        Vy[i] = Vy[i] + Fy[i]*dt/mass[i+first];
      }	
      
    }  /* end of while-loop */



    // Use Process 0 to print time and write final status to file.
    if(id == 0){
      // End timer
      endtime = MPI_Wtime();

      printf("\n");
      printf("Time: %6.2f seconds\n", (endtime - starttime));
      // Write final particle coordinates to a file
      write_particles(N, X, Y, "final_pos_pp.txt");
    }

    // Clean up allocated memory
    free(X);
    free(Y);
    free(mass);
    free(tempX);
    free(tempY);
    free(Vx);
    free(Vy);
    free(Fx);
    free(Fy);

    MPI_Finalize();
    exit(0);   
}