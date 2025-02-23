#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

/* For correct usage look at the bottom of this file. */

#define PI 3.14159265358979323846

int N = 16777216;//1024//16384//16777216//*4
bool dgbMode = false;
bool printRes= false;


bool is_pwr_of_two(unsigned int );
unsigned int bitFlip(unsigned int, int, int);
unsigned int FlipBit(unsigned int, unsigned int);
void fft(double complex *, double complex *, int, int);
void communicate(int splitPoint, int dataChunk, int rank, double complex *, double complex *);
void saveDisk(int rank, int currData, int size, double *, FILE *file);
void dispInfo();

int main(int argc, char ** argv)
{


FILE *file = fopen("InputFFT.txt", "w");

  if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }


double start, finish;

  MPI_Init(&argc, &argv);
  
  int rank, size;
  int root=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    /* User input choices. */
    /* No choices = default settings. */
    if (argc == 4) 
    {
      int fftSize,debugMode,printOutput;

      fftSize  = atoi(argv[1]);
      debugMode = atoi(argv[2]);
      printOutput = atoi(argv[3]);

      /* Set the number (globally) of points to calculate the fft. */
      N=fftSize;

      if (debugMode==1)
      dgbMode = true;

      if (printOutput==1)
      printRes= true;
    }
    else
    { 
       if (rank==0)
       printf("\nRunning on the default settings, N=%d\n", N); 
    }

    /* Initializations. */
    double complex * arr;
    double complex * y;
    double complex * output;
    double complex * resultRev;

    double z;
    int currData;

    /* The program expects powers of two, both for processes and input data size. */
    if (is_pwr_of_two(N) && is_pwr_of_two(size) && (N>=size))
    {
     currData = N/size;

     if (rank==0)
      printf("\nCalculating... \n");
    }
    else
    {
      if (rank==0)
      dispInfo();
      remove("InputFFT.txt");
      MPI_Finalize();
      return 0;
    }
    
    /* Random Initialization per process. */
    srand(rank+1);

    arr = (double complex *) malloc(sizeof(double complex)*currData);

/* In case we need to print the input sequence to a file.*/
if(dgbMode)

    {
      double* barr;
      barr = (double*) malloc (currData*sizeof(double));

     for (int i = 0; i < currData; i++)
      {
         z=rand() % 40 - 20;
         arr[i]=z;
         barr[i]=z;
      }
      /* Wait all processes to come to this point. */
      /* Write data to disk. */
      MPI_Barrier(MPI_COMM_WORLD);
      saveDisk(rank, currData, size, barr, file);
      free(barr);
  }

  else
   {
    remove("InputFFT.txt");
      for (int i = 0; i < currData; i++)
      {
         z=rand() % 40 - 20;
         arr[i]=z;
      }
   }
    
    output = (double complex *) malloc(sizeof(double complex)*N);
    y = (double complex *) malloc(sizeof(double complex)*currData);

    resultRev = (double complex *) malloc(sizeof(double complex)*N);
    /* Wait all processes to come to this point.*/
    MPI_Barrier(MPI_COMM_WORLD);
    start=MPI_Wtime(); /*start timer*/
    fft(arr, y, rank, size);

    /* Get the elements from the other processes. */
    MPI_Gather(y, N/size, MPI_C_DOUBLE_COMPLEX, output, N/size, 
               MPI_C_DOUBLE_COMPLEX, root, MPI_COMM_WORLD);

    finish=MPI_Wtime(); /*end timer*/

/* In case we need to print the FFT output sequence 
in a file (with the same order as the input).*/
if (printRes)
{
    int i=0,j=0;
    if(rank == 0)
    {
      for( i = 0; i < N; i++)
      {
        resultRev[FlipBit(i,(int) ceil(log2(N)))]=output[i];
      }

      FILE *fout = fopen("OutputFFT.txt", "w");

      for( j = 0; j < N; j++)
      {
       fprintf(fout,"y[%d] = %lf + %lf*i\n", j, creal(resultRev[j]), cimag(resultRev[j]));
      }
    }//rank==0
}//printRes


    free(arr);
    free(y);
    free(output);
    free(resultRev);

  
  MPI_Finalize();
  printf("Parallel Elapsed time: %f seconds\n", finish-start);

  return 0;
}


void fft(double complex * R,double complex * y, int rank, int size)
{

  int nbits = (int) ceil(log2(N));
  int dataChunk = N/size;
  int start = rank*dataChunk;
  double start1, finish1;
 
  double complex * S = malloc(sizeof(double complex)*dataChunk);
  double complex * Sk = malloc(sizeof(double complex)*dataChunk);



  for(int m = 0; m < nbits; m++)
  {

    for(int i = 0; i < dataChunk; i++)
    {
      /* Update the partial sums. */
      Sk[i] = S[i] = R[i];
    }

    unsigned int bit = 1 << (nbits - m - 1); /* N/2 over 2 over 2 over 2 ... 1 */

    /* Do the butterfly. */

    int splitPoint = size / (1 << (m+1));  /* size/2 , size/4 ... size/(N) */

    if(splitPoint > 0)

    {
      start1=MPI_Wtime(); /*start timer*/
     
      /* Do the communication part. */
      communicate(splitPoint, dataChunk, rank, S, Sk);

      finish1=MPI_Wtime(); /*end timer*/
    }

    else
      /* No communication is needed. */
    {
        for(int i = 0; i < dataChunk; i++)
        Sk[i] = S[i];
    }
    /* Compute the sums given the partial sums. */
    for(int i = start, l = 0; l < dataChunk; i++, l++)
    {

      int j = (i & (~bit)) % dataChunk;
      int k = (i | bit) % dataChunk;

      /* Compute the twiddle factor. */
      int expFactor = bitFlip(i, nbits, m);
      R[l] = S[j] + Sk[k] * cexp( (2*PI*I*expFactor)/N );
    }
  }

  for(int i = 0; i < dataChunk; i++)
  {
    /* Update the results. */
    y[i] = R[i];
  }

  free(S);
  free(Sk);
  printf("Communication Elapsed time: %f seconds\n", finish1-start1);
}


unsigned int bitFlip(unsigned int i, int bits, int m)
{
    unsigned int out = 0;
    /* Computes the twiddle factor. */
    i = i >> (bits-m-1);

    for(int j = 0; j < m+1; j++)
    {
      out |= i & 1;
      i = i >> 1;
      if(j < m)
        out = out << 1;
    }

    out = out << (bits-m-1);

    return out;
}


unsigned int FlipBit(unsigned int n, unsigned int bits)
{   
    /* Bit Reversal */
    unsigned int nrev, L;
    unsigned int count;   
    L = 1<<bits;
    count = bits-1;   /* initialize the count variable */
    nrev = n;
    for(n>>=1; n; n>>=1)
    {
        nrev <<= 1;
        nrev |= n & 1;
        count--;
    }

    nrev <<= count;
    nrev &= L - 1;

    return nrev;
}

void communicate(int splitPoint, int dataChunk, int rank, double complex * S, double complex * Sk)
{
    int destination, source;
    int err;  // Variable to capture error codes from MPI functions

    if ((rank % (splitPoint * 2)) < splitPoint)
    {
        /* My rank is < than the one that i have to communicate with. */
        /* This process sends S and receives into Sk. */
        destination = source = rank + splitPoint;
        err = MPI_Sendrecv(S, dataChunk, MPI_C_DOUBLE_COMPLEX, destination, 0,
                             Sk, dataChunk, MPI_C_DOUBLE_COMPLEX, source, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (err != MPI_SUCCESS)
        {
            char error_string[BUFSIZ];
            int length_of_error_string;
            MPI_Error_string(err, error_string, &length_of_error_string);
            printf("MPI Error in Sendrecv (rank %d): %s\n", rank, error_string);
            MPI_Abort(MPI_COMM_WORLD, err);
        }
    }
    else
    {
        /* My rank is > than the one that i have to communicate with. */
        /* This process sends Sk and receives into S. */
        destination = source = rank - splitPoint;
        err = MPI_Sendrecv(Sk, dataChunk, MPI_C_DOUBLE_COMPLEX, destination, 0,
                             S, dataChunk, MPI_C_DOUBLE_COMPLEX, source, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (err != MPI_SUCCESS)
        {
            char error_string[BUFSIZ];
            int length_of_error_string;
            MPI_Error_string(err, error_string, &length_of_error_string);
            printf("MPI Error in Sendrecv (rank %d): %s\n", rank, error_string);
            MPI_Abort(MPI_COMM_WORLD, err);
        }
    }
}


void saveDisk(int rank, int currData, int size, double *barr, FILE *file)
{

MPI_Status status;

int a; 

if (rank == 0) 
  {

    for( a = 0; a < currData; a++ )
      {
       //printf("I am proc: %d %f \n",rank, barr[a]);
       fprintf(file, "%f\n",barr[a]);
      }
  fclose(file);
  //Send a message to the next process to continue 
  // writting on the disk on append 'a' mode.

  //We care about the message not its content so
  //we use mpi_bottom address.
    if (size > 1) 

    MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, 1, 2373, MPI_COMM_WORLD);

  }//if p==0
else
  {
    //Wait for message from previous process 
    MPI_Recv(MPI_BOTTOM, 0, MPI_DOUBLE, rank - 1, 2373, MPI_COMM_WORLD, &status);
    file=fopen("InputFFT.txt", "a");

    for( a = 0; a < currData; a++ )
      {
       //  printf("I am proc: %d %f \n",rank, barr[a]);
         fprintf(file, "%f\n",barr[a]);
      }


  fclose(file);

  // If not reached end process yet
  // send message to the next process
  if (rank < size-1)
  MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, rank + 1, 2373, MPI_COMM_WORLD);

  } 
}

 /* Check if number is a power of two. */
 bool is_pwr_of_two(unsigned int x) 
 {
   return x && !(x & (x - 1));
 }

 void dispInfo()
 {  printf("\n\n");
    printf("\t\t|~~~~~~~~~~~~~~ Something went wrong! ~~~~~~~~~~~~|\n");
    printf("\t\t|                                                 |\n");
    printf("\t\t|----------------- House Rules: ------------------|\n");
    printf("\t\t|                                                 |\n");
    printf("\t\t| 1) #CPUS and #Inputs are Powers of two.         |\n");
    printf("\t\t| 2) #CPUS < = #Inputs.                           |\n");
    printf("\t\t| 3) To print FFT input  -> dgbMode  = true.      |\n");
    printf("\t\t| 4) To print FFT output -> printRes = true.      |\n");
    printf("\t\t|                                                 |\n");
    printf("\t\t|                Run  Example:                    |\n");
    printf("\t\t|                                                 |\n");
    printf("\t\t| (8 processes, 128 elements, no print, yes print)|\n");
    printf("\t\t|          mpirun -np 8 myfft 128 0 1             |\n");
    printf("\t\t|                                                 |\n");
    printf("\t\t|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|\n\n\n");
 }
