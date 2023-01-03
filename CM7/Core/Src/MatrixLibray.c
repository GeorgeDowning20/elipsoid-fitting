#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ZERO 0
#define ONES 1
#define RANDOM 2
#define IDENTITY 3




float **createNewMatrix(int mode, int r, int c)
{
    int i, j;
    float **mat = (float **)malloc(r * sizeof(float *));

        for (i = 0; i < r; i++)
        {
            mat[i] = (float *)malloc(c * sizeof(float));

            for (j = 0; j < c; j++)
            {
            	switch (mode){

            	case 0:
            	mat[i][j] = 0;
            	break;

            	case 1:
            	mat[i][j] = 1;
            	break;

            	case 2:
            	mat[i][j] = rand() % 10;
            	break;

            	case 3:
            	mat[i][j] = i == j ? 1 : 0;
            	break;
            	}
            }
        }

    return mat;
}

float **normaliseMatrix (float **Matrix, int r, int c)
{
	float **MATRIX = createNewMatrix(0, r, c);
	float K = -Matrix[r-1][c-1];

    int i, j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
        	MATRIX[i][j] = Matrix[i][j]/K;
        }
    }
    return MATRIX;
}

float **ArraytoMatrix (float *array , int r , int c)
{
	float ** matrix = createNewMatrix(0, r, c);
    int i, j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
        	matrix[i][j] = array[(i*c)+j];
        }
    }
    return matrix;
}

float *MatrixtoArray (float **Matrix , int r , int c, float *ArrayPTR)
{
	//float Array[r][c];
    int i, j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
        	ArrayPTR[(i*c)+j] = Matrix[i][j];
        }
    }
    return ArrayPTR;
}
#
float *MatrixtoArrayB (float **Matrix , int r , int c, float *ArrayPTR)
{
	//float Array[r][c];
    int i, j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
        	ArrayPTR[(i*9)+j] = Matrix[i][j];
        }
    }
    return ArrayPTR;
}



void printArray (float *Array , int r , int c)
{
	printf("\r\n");
	    int i, j;
	    for (i = 0; i < r; i++)
	    {
	        for (j = 0; j < c; j++)
	        {
	            printf("%f ", Array[(i*c)+j]);
	        }
	        printf("\r\n");
	    }
	    printf("\r\n");
}


float **MultiplyMatrix(float **matrix1, int r1, int c1, float **matrix2, int r2, int c2)
{
	if (r2 != c1)
	return NULL;

    int i, j, k;
    float **result = createNewMatrix(ZERO, r1, c2);

    for (i = 0; i < r1; i++)
    {
        for (j = 0; j < c2; j++)
        {
            for (k = 0; k < r2; k++)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

float **TransposeMatrix(float **matrix, int r, int c)
{
    int i, j;
    float **result = createNewMatrix(ONES, c, r);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
                result[j][i] = matrix[i][j];
        }
    }
    return result;
}
float **SquareMatrix(float **matrix, int r, int c)
{
	float **rusult = MultiplyMatrix(TransposeMatrix(matrix,r,c),c,r,matrix,r,c);
    return rusult;
}
float **NegitiveMatrix (float **Matrix,int r,int c)
{
	int i, j;
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			Matrix[i][j] = -Matrix[i][j];
		}
	}

	return Matrix;
}





