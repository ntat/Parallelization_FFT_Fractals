#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h> 
#include <mpi.h>


// Create A Dynamic Matrix stored on the Heap (especially good when 
// we need to work with big Matrices that don't fit in the stack.) 

int** makeArray(unsigned nrows, unsigned ncols)
{
   unsigned i;
   int** ptr = malloc(nrows * sizeof(int*));  // allocate pointers to rows
   if ( ptr != NULL )
   { 
       int* chunk = malloc(nrows * ncols * sizeof(int)); // allocate chunk of memory
       if ( chunk != NULL )
       {
           for (i = 0; i < nrows; ++i, chunk += ncols )
               ptr[i] = chunk;  // point the row pointers into the chunk
       }
       else
       { 
          free(ptr);  
          ptr = NULL;
       }
   }
   return ptr;
}

//Function to free() the allocated memory. 
void freeArray(int** arr)
{
   free(arr[0]);  // free the chunk of memory
   free(arr);     // free the pointers
}

//Globals
double MaxRe = 2.0;
double MinRe = -2.0;

double MaxIm = 2.0;
double MinIm = -2.0;


// Image dimensions 
int w = 4096; 
int h = 4096;

double newMinRe,newMaxRe,newMinIm,newMaxIm,dx,dy;

int count(double complex d, int b, int N);

// Creates a square around this (x,y) 
// pixel to zoom into the mandelbrot set.
// Factor controls how big is the area.
void cordTransZoom(int x_0, int y_0, int factor);

int main(int argc, char **argv)
{

    FILE *fp;
    double complex d;

    int rank, size, tag, rc, i;
    MPI_Status status;

    //Initialize MPI
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    tag = 132;

    int N = 1024; 
    int b = 2;

    if((fp=fopen("test", "wb"))==NULL) {
    printf("Cannot open file.\n");
    }

    double dreal,dimag;

    //Step
    dx=(float)(MaxRe-MinRe)/(float)(w-1);
    dy=(float)(MaxIm-MinIm)/(float)(h-1);

    newMinRe=MinRe;
    newMinIm=MinIm;

    

    if (argc == 4)//User passes zoom pixel cordinates. 
    {
      //Zoom 1x around pixel(x_0,y_0) from input given
      // (x,y) and scale
     cordTransZoom(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
    }
    
    if (argc > 4) //Just some fancy zoom.
    {
      //Zoom 2x around pixel(x_0,y_0).
       cordTransZoom(1291,1860,3);
       cordTransZoom(945,1866,10);
    }
    //Else we calculate the normal mandelbrot set.


    int wp = w/size;
    int hp = h;
    int xoff = rank*w/size;
    int yoff = 0;

    int r,s,l,k,x,y;

    printf("SubMatrix Size: %d, %d\n", wp,h);

    int **temp = makeArray(wp, hp);

    for ( x = 0; x < wp; x++)

    {
        dreal=(x+xoff)*dx + newMinRe;

        for ( y = 0; y < h; y++)
        {
            dimag=(y+yoff)*dy + newMinIm;

            d=dreal+dimag*I;
            temp[x][y]=count(d, b, N);
        

        }
    
    }


printf("I am process #%d\n", rank);

    if (rank == 0) // If master
    {
        int **img = makeArray(w, h);
          // Firstly, insert results of 0th master
          // process elements into the final image.
        for ( s = 0; s < wp; s++)

                {
                    
                    for ( r = 0; r < h; r++)
                    {
                        img[s][r]=temp[s][r];
                    }
                }



        for (i = 1; i < size; i++)
        {   //Receive the submatrix from each process
            rc = MPI_Recv(&temp[0][0], wp*hp, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

            int cnt=0;
            for ( k = i*wp; k < (i+1)*wp; k++)

                {
                    
                    for ( l = 0; l < h; l++)
                    {
                        img[k][l]=temp[cnt][l];
                       // printf("%d\n", img[k][l] );
                    }
                    cnt++;
                }
        }


        freeArray(temp);
        int i,j;
        //Write output in binary format (readable by matlab).
        for ( i = 0; i < w; i++)
        {
            for (j = 0; j < h; j++)
            {
            fwrite(&img[i][j], sizeof(int), 1, fp);
            }
        }

        fclose(fp); 
        freeArray(img);

    }// END_if_(rank == 0)

    else
    {
        //Send subMatrices to master
        rc = MPI_Send(&temp[0][0], wp*hp, MPI_INT, 0, tag, MPI_COMM_WORLD);
        freeArray(temp);
    }


    //Finish Session
    rc = MPI_Finalize();
   return 0;
}

// Evaluate Mandelbrot function and
// assign corresponding colors.
int count(double complex d, int b, int N) 
{
    int count =1;
    double complex z=0;

   while((cabs(z) < b) && (count < N)) 
   {
      z= cpow(z,2)+d;
      count=count+1;
   }
 
return count; 
}

// Transformation from Pixel coordinate space
// into complex plane.
void cordTransZoom(int x_0, int y_0, int factor) 
{
  //Upper Left Corner of zoom window (x1,y1)
  int upperLeftX=x_0-factor;
  int upperLeftY=y_0-factor;
  //Lower Right Corner of zoom window (x2,y2)
  int lowerRightX=x_0+factor;
  int lowerRightY=y_0+factor;

  //Coord Tranformation
   newMinRe = MinRe + (dx* upperLeftX);
   newMaxRe = MinRe + (dx* lowerRightX);
   newMinIm = MinIm + (dy* upperLeftY);
   newMaxIm = MinIm + (dy* lowerRightY);
    
  dx=(float)(newMaxRe-newMinRe)/(float)(w-1);
  dy=(float)(newMaxIm-newMinIm)/(float)(h-1);

  MinRe=newMinRe;
  MinIm=newMinIm;   
}

