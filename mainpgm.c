#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

void matmult (int m, int n, int p, double lhs[][n], double rhs[][p],  double result[][p]);
void matranspose (int m, int n, double mat[][n], double transp[][m]);
double signe (int i, int j);
double determinant (int m, double mat[][m]);
void cofactor (int m, double mat [][m], double cofact [][m]);
void inverse (int m, double mat[][m], double inv[][m]);
int random_training_data_selector (int dataElements, int percentOfData, unsigned char bitmap [dataElements / 8 + 1]);
int SetBit (unsigned char arr[], int position);
int isTrainingDataRecord (int dataElementPos, unsigned char selectorBitmap []);
void freader (FILE *fp, double *valptr);
void solveRegressionEqn (int r, int c, double X[r] [c], double Y [r] [1], double B [r] [1]);
void solve ( int r, int c, double X[r][c], double B [c] [1], double Y [r][1]);


void main (int argc, char *argv[])
{
	int nRows, nCols, bHasNoHeader, nTrainingRows, nTestRows, fret;
	long int fpos;
	unsigned char *bmp;
	char filename [50];
	FILE *fp;
	int TestCounter = 0, TrainingCounter = 0;
	int i, j, k;
	
	strcpy (filename, argv [1]);
	sscanf (argv[2], "%d", &nRows);
	sscanf (argv[3], "%d", &nCols);
	sscanf (argv[4], "%d", &bHasNoHeader);
	
	double inputs[1][nCols];
	if( !bHasNoHeader) nRows--;
	
	bmp = (unsigned char *) calloc (nRows / 8 + 1, sizeof (unsigned char));

	// Choose the percentage of data that is training data (the remaining is test data)
	nTrainingRows =  random_training_data_selector (nRows, 95, bmp);
	nTestRows = nRows - nTrainingRows;

	double TrainingSetX [nTrainingRows] [nCols], TrainingSetY [nTrainingRows] [1], Coefficients [nCols][1];
	double TestSetX [nTestRows][nCols], TestSetY [nTestRows] [1], CalculatedY [nTestRows] [1];
	char header_names [nCols][50];

	fp = fopen ( filename, "r");

	for (i = 0; i < nCols; i++)
		if (bHasNoHeader) 
			if (i == 0)
				sprintf (header_names [i], "y");
			else
				sprintf (header_names [i], "x%d", i);
		else
			fscanf (fp, "%s", header_names [i]);


	for ( i=0; i< nRows; i++)
		{
			if ( isTrainingDataRecord (i, bmp))
			{
				freader (fp, &TrainingSetY [TrainingCounter] [0]);
				for( j=0; j< nCols; j++){
					if( j==0){
						 TrainingSetX [TrainingCounter][j]=1;
						 continue;
						}
					freader (fp, &TrainingSetX [TrainingCounter] [j]);
				}
				TrainingCounter++;
			}
			else
			{
				freader (fp, &TestSetY [TestCounter] [0]);
				for( j=0; j< nCols; j++){
					if(j==0){
						TestSetX [TestCounter][j]=1;
						continue;
					}
					freader (fp, &TestSetX [TestCounter][j]);
				}
				TestCounter++;
			}
		}

	fclose (fp);
	cfree (bmp);

	solveRegressionEqn (nTrainingRows, nCols, TrainingSetX, TrainingSetY, Coefficients);
	solve ( nTestRows, nCols, TestSetX, Coefficients, CalculatedY);


	double diff = 0;
	double SSres = 0;		// residual sum of squares
	double SStot = 0i, avgY = 0;
	double Rsquared;
	double StdErrorRegression = 0;

	for ( i = 0; i < nTestRows; i++) avgY += TestSetY [i][0];
	avgY = avgY / nTestRows;

	for (int i =0; i< nTestRows; i++)
	{
		diff = TestSetY [i][0] - CalculatedY [i][0];
		SSres += diff * diff;

		diff = TestSetY [i][0] - avgY;
		SStot += diff * diff;
	} 

	Rsquared = 1.0 - SSres / SStot;
	StdErrorRegression = sqrt (SSres / nTestRows);

	printf ("Std Error of Regression is: %f\n", StdErrorRegression);

	printf ("\nThe multiple linear regression equation is:\n%s = ", header_names [0]);
	for (int k = 0; k < nCols; k++) {
		printf ("%f", Coefficients [k][0]);
		if (k > 0)
			printf (" * %s ", header_names [k]);
		if (k != nCols - 1)
			printf (" + ");
	}
	printf ("\n\n");

	printf("To calculate %s :\n", header_names [0]);
	inputs[0][0]=1;
	for( i =1; i<nCols; i++)
	{
		printf("Enter %s\n", header_names[i]);
		scanf("%lf", &inputs[0][i]);
	}
	double price[1][1];
	solve (1, nCols, inputs, Coefficients, price);
	printf("The price for the given values is : %lf", price[0][0]); 
}
//---------------------------------------------------------------------------------------
void solve ( int r, int c, double X[r][c], double B [c] [1], double Y [r][1])
{
	matmult(r, c, 1, X, B, Y);
}
//---------------------------------------------------------------------------------------
void solveRegressionEqn (int r, int c, double X[r] [c], double Y [r] [1], double B [c] [1])
{
	double transpose [c][r], mult1 [c][c], inv [c][c], mult2 [c][r];

	matranspose (r, c, X, transpose);
	matmult (c, r, c, transpose, X, mult1);
	inverse ( c, mult1, inv);
	matmult (c, c, r, inv, transpose, mult2);
	matmult ( c, r, 1, mult2, Y, B);	
}
//---------------------------------------------------------------------------------------
void freader (FILE *fp, double *valptr)
{
	long int fpos;
	char tempstr [100], *tstrptr = tempstr;

	fpos = ftell (fp);
	if (!fscanf (fp, "%lf", valptr)) {
		fseek (fp, fpos, SEEK_SET);
		fscanf (fp, "%s", tstrptr);
		if (tstrptr [0] == '"') tstrptr++;
		sscanf (tstrptr, "%lf", valptr);
	}
}
//---------------------------------------------------------------------------------------
void matmult (int m, int n, int p, double lhs[][n], double rhs[][p],  double result[][p])
{
	// lhs must be a m x n matrix and rhs must be a n x p matrix
	// result will be a m x p matrix
	
	int i, j, k;
	for (i = 0; i < m; i++)
	    for (k = 0; k < p; k++)
		for (j = 0; j < n; j++)  {
		    if (j ==0) result [i][k] = 0;
		    result [i][k] += lhs [i][j] * rhs [j][k];
		}
}
//----------------------------------------------------------------------------------------
void matranspose (int m, int n, double mat[][n], double transp[][m])
{
	// Transpose a m x n matrix into a n x m matrix

	int i, j;
	for (i = 0; i < m; i++)
	     for (j = 0; j < n; j++)
		transp[j][i] = mat [i][j];
}
//----------------------------------------------------------------------------------------
void minor_matrix (int which_row, int which_col, int sz, double matrix[][sz], double minor[][sz-1])
{
	// Determine the minor of a sz x sz matrix, eliminating which_row and which_col data

	int i, j;
	int r = -1, c = -1;
	
	for (i = 0; i < sz; i++) {
		if (i == which_row) continue;
		else r++;

		c = -1;

		for (j = 0; j < sz; j++) {
			if (j == which_col) continue;
			else c++;

			minor [r][c] = matrix [i][j];
		}
	}
}
//----------------------------------------------------------------------------------------
double determinant (int m, double mat[][m])
{
	int i;
	double minor [m-1][m-1];
	double determ = 0;

	if (m == 1) return mat [0][0];

	for (i = 0; i < m; i++)   {
		minor_matrix (0, i, m, mat, minor);
		determ += mat [0][i] * signe (0, i) * determinant (m - 1, minor);
		if (m > 10) { putchar ('d');  fflush (stdout); }
	}
	return determ;
}
//----------------------------------------------------------------------------------------
double signe (int i, int j)
{
	// For determinant and cofactor calculations, the sign of a term depends on the row and column
	// If the sum of row no and col no is odd, then the sign is negative else the sign is positive

	return (((i + j) % 2 == 1) ? -1 : 1);
}
//----------------------------------------------------------------------------------------
void cofactor (int m, double mat [][m], double cofact [][m])
{
	int i, j;
	double minnie [m-1][m-1];
	int ii, jj;
	
	for (i = 0; i < m; i++) {	
		for (j = 0; j < m; j++) {
			minor_matrix (i, j, m, mat, minnie);
			cofact [i][j] = signe (i, j) * determinant (m - 1, minnie);
			if (m > 10) { putchar ('C'); fflush (stdout); }
		}
		if (m > 10) printf ("\nCompleted %d\% \n", (i+1)*100 / m);
	}
}
//---------------------------------------------------------------------------------------
void inverse (int m, double mat[][m], double inv[][m])
{
	if (m > 10) printf ("Starting inverse of matrix of size %d ..\n", m);

	int i,j;
	double det;
	double cofact[m][m], adjoint[m][m];

	if (m > 10) printf ("Computing cofactor:\n");
	cofactor( m, mat, cofact);

	matranspose (m, m, cofact, adjoint);

	if (m > 10) printf ("Computing determinant:\n");
	det=determinant( m, mat);

	for(i=0; i<m; i++)
		for(j=0;j<m;j++)
			inv[i][j] = (1.00/det)* adjoint[i][j];
}
//----------------------------------------------------------------------------------------
int random_training_data_selector (int dataElements, int percentOfData, unsigned char bitmap [dataElements / 8 + 1])

// Given dataElements number of elements, and a percentage of data to select (as training data), fill the provided 
// bitmap with 1s and 0s, with 1 representing training data and 0 representing test data. The choice of elements as
// training or test data elements is driven by a random function.

{
	int i, index, targetNo, randomNosCount = 0; 

	// Since dataElements may not be an exact multiple of 8 (for bitwise storage), use 1 additional byte

	for (i = 0; i < (dataElements / 8 + 1); i++)
		bitmap [i] = 0;

	srand ((unsigned) time (NULL));

	targetNo = dataElements * percentOfData / 100;

	// rand () returns values from 0 to RAND_MAX. However, we want values between 0 and dataElements - 1

	while (randomNosCount < targetNo) {
		index = (int) ((double) rand () / RAND_MAX * (dataElements - 1));
		if (!SetBit (bitmap, index))
			randomNosCount++;
	}

	return targetNo;
}
//---------------------------------------------------------------------------------------
int SetBit (unsigned char arr[], int position)

// Set the bit in position in array arr to 1. If it was not already set, return 0 to caller
// else return 1

{
	int arrayIndex, whichBitToSet, flag = 1, currentValue;

	arrayIndex = position / 8;
	whichBitToSet = position % 8;
	currentValue = (arr [arrayIndex] >> whichBitToSet) & 1;

	flag = flag << whichBitToSet;
	arr [arrayIndex] |= flag;
	return currentValue;
	
}
//--------------------------------------------------------------------------------------
int isTrainingDataRecord (int dataElementPos, unsigned char selectorBitmap [])

// A bit mapped representation of all the records read from the file is in selectorBitmap. This 
// function returns whether the row of data represented by dataElementPos is part of training 
// data set or not (else it is part of test data set)

{
	return ((selectorBitmap [dataElementPos / 8] >> (dataElementPos % 8)) & 1);
}
