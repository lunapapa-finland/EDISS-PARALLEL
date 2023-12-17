#include <stdlib.h>

#include <stdio.h>

#include <string.h>

#include <mpi.h>

#include <time.h>

int main(int argc, char * argv[]) {

  // Define Variables   
  int MAXSIZE = 1000000000;
  int np, id;
  double starttime, endtime;
  int inputSize;
  char * buffer;
  // Initialize MPI 
  MPI_Init( & argc, & argv);
  MPI_Comm_size(MPI_COMM_WORLD, & np);
  MPI_Comm_rank(MPI_COMM_WORLD, & id);

//   While loop, break criteria is either input is less than or equal to 0, or larger than 1,000,000,000
  while (1) {
    // get the size of message from process 0 
    if (id == 0) {
      // Repeat until an inputSize <= 0 is entered
      printf("Please give an input size in bytes:\n");
      fflush(stdout);
      scanf("%d", & inputSize);
    }
    // Broadcast the size of the data to the each process other than 0
    MPI_Bcast( & inputSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // if the inputsieze if <=0 or > MAXSIZE, printout hints, break the loop and finalize the MPI.
    if (inputSize > MAXSIZE) {
      if (id == 0) {
        printf("Input size is too large, maximum value is %d\n", MAXSIZE);
      }
      continue;
      // break;
    }else if (inputSize <= 0) break;
    
      /* Allocate a sufficiently large message buffer for each process*/
    buffer = (char * ) malloc(inputSize * sizeof(char));
    if (buffer == NULL) {
      if (id == 0) {
        printf("Mem allocation failed");
      }
      break;
    }

    // if process 0, print message size, initialize necessary data to allocated buffer
    if (id == 0) {
        printf("Message size: %d\n", inputSize);
        fflush(stdout);
        // Initialize the whole message buffer to some values
        for (int i = 0; i < inputSize; i++) {
            buffer[i] = (char) 75;
        }
        // send the Message to (id + 1) % np , i.e. process 1
      MPI_Send(buffer, inputSize, MPI_CHAR, (id + 1) % np, 0,
        MPI_COMM_WORLD);
        // start to record time
      starttime = MPI_Wtime();
        // wait for the message back after passing through the whole ring at the very last
      MPI_Recv(buffer, inputSize, MPI_CHAR, np - 1, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // stop to record time
      endtime = MPI_Wtime();
        //   print out total time
      printf("Total time for sending size of %d bytes message in %f s\n", inputSize, endtime - starttime);
    } else {
        // each process other than 0 will firstly wait for the message sending from its previous node 
      MPI_Recv(buffer, inputSize, MPI_CHAR, id - 1, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // then send the message to the next node, mod is used to calculate process 0 as the next one if the current process is the last node
      MPI_Send(buffer, inputSize, MPI_CHAR, (id + 1) % np, 0,
        MPI_COMM_WORLD);
    }
  }
  free(buffer); // Deallocate the buffer                                                                
  MPI_Finalize(); // Terminate MPI                                                                        
  exit(0); // Quit     
}