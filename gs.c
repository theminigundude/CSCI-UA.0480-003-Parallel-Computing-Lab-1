#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/*
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;

  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);

    for(j = 0; j < num; j++)
       if( j != i)
	 sum += fabs(a[i][j]);

    if( aii < sum)
    {
      printf("The matrix will not converge.\n");
      exit(1);
    }

    if(aii > sum)
      bigger++;

  }

  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;

  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
  {
	printf("Cannot allocate a!\n");
	exit(1);
  }

 for(i = 0; i < num; i++)
  {
    a[i] = (float *)malloc(num * sizeof(float));
    if( !a[i])
  	{
		printf("Cannot allocate a[%d]!\n",i);
		exit(1);
  	}
  }

 x = (float *) malloc(num * sizeof(float));
 if( !x)
  {
	printf("Cannot allocate x!\n");
	exit(1);
  }


 b = (float *) malloc(num * sizeof(float));
 if( !b)
  {
	printf("Cannot allocate b!\n");
	exit(1);
  }

 /* Now .. Filling the blanks */

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
	fscanf(fp,"%f ", &x[i]);

 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i][j]);

   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
 }

 fclose(fp);

}

void MPI_CALCULATE(int dim, float prec, float *values, float **matrix, float *constants, int begin, int end, float *toSendValues)
{
  // init loop counter variable
  int i = 0;
  // Count how many are within precison
  int numCount = 0;

  for(i = begin; i < (end-begin); i++){
    toSendValues[i-begin] = values[i];
  }

  for(i = begin; i< end; i++){
    float lastNum = constants[i];
    int symbols = 0;
    for(symbols = 0; symbols < dim; symbols++){
	     if (symbols != i){
         lastNum -= (matrix[i][symbols] * values[symbols]);
       }
    }
    toSendValues[i-begin] = lastNum / matrix[i][i];
  }
}

int error (float *x, float *newX)
{
	int counter = 0;
	for(int i = 0; i < num; i++){
		if (fabs( (newX[i] - x[i]) / newX[i] )<= err){
			counter ++;
		}
	}

	return counter;
}

void divideWork (int comm_sz, int *length, int *prev_count)
{
	int dividedWork = num / comm_sz;
	int remainder = num % comm_sz;

	length [0] = dividedWork + remainder;
	prev_count[0] = 0;

	for(int i = 1; i < comm_sz; i++){
		length[i] = dividedWork;
		prev_count[i] = length[i-1] + prev_count[i-1];
	}

}




/************************************************************/


int main(int argc, char *argv[])
{

 int i;
 int nit = 0; /* number of iterations */
 FILE * fp;
 char output[100] ="";

 if( argc != 2)
 {
   printf("Usage: ./gsref filename\n");
   exit(1);
 }

 /* Read the input file and fill the global data structure above */
 get_input(argv[1]);

 /* Check for convergence condition */
 /* This function will exit the program if the coffeicient will never converge to
  * the needed absolute error.
  * This is not expected to happen for this programming assignment.
  */
 check_matrix();

 int comm_sz;
 int my_rank;
 int condition = 0;

 MPI_Init(NULL, NULL);
 MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 int precisionNum = 0;
 int length[comm_sz];
 int prev_count[comm_sz];
 float newX[num];

 //divide work for processors
 divideWork(comm_sz, length, prev_count);

 //calculate
 while (precisionNum != num){
   nit ++;
   precisionNum = 0;
   float tempX[length[my_rank]];

   MPI_CALCULATE(num, err, x, a, b, prev_count[my_rank], prev_count[my_rank] + length[my_rank], tempX);
   MPI_Allgatherv(&tempX, length[my_rank], MPI_FLOAT, &newX, length, prev_count, MPI_FLOAT, MPI_COMM_WORLD);

   precisionNum = error(x, newX);

   for(int i = 0; i < num; i++){
     x[i] = newX[i];
   }

 }

 MPI_Finalize();

 if (my_rank == 0){
   /* Writing results to file */
   sprintf(output,"%d.sol",num);
   fp = fopen(output,"w");
   if(!fp)
   {
     printf("Cannot create the file %s\n", output);
     exit(1);
   }

   for( i = 0; i < num; i++)
     fprintf(fp,"%f\n",x[i]);

   printf("total number of iterations: %d\n", nit);
   fclose(fp);
 }

 exit(0);

}
