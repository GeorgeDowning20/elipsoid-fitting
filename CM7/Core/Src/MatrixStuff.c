#include "main.h"
#include "math.h"


float detrimant123;
float detrminant(float a[9][9], float k);
float **cofactors(float num[9][9], float f);
float a[];
float **inverseMatrix(float ** Matrix,int size);
float ** ArrayElipsoidMat (float **ArrayM, int samples);
float **EigenValues3x3 (float **Matrix);
void printMatrix(float **matrix, int r, int c)
{
	printf("\r\n[");
    int i, j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            printf("%f ", matrix[i][j]);
        }
        if (!(i == (r-1)))
        printf(";\r\n");
    }
    printf("]\r\n");
}

float **get_cubic(float a,float b,float c,float d)
{
float e,f,g,h,i,j,k,l,m,n,p,r,s,t,u,x1,x2,x3;
int w;
printf("\r\n%f\t%f\t%f\t%f\t\r\n",a,b,c,d);


e=2.7182818284590;
f=((3*c/a)-(b*b/(a*a)))/3;
g=((2*b*b*b/(a*a*a))-(9*b*c/(a*a))+(27*d/a))/27;
h=(g*g/4)+(f*f*f/27);
i=sqrtf(((g*g/4)-h));
j=expf(log10(i)/log10(e)/3);
k=acosf((-1)*(g/(2*i)));
l=j*(-1);
m=cosf(k/3);
n=sqrtf(3)*sinf(k/3);
p=(b/3*a)*(-1);
r=(-1)*(g/2)+sqrtf(h);
s=expf(log10(r)/log10(e)/3);
t=(-1)*(g/2)-sqrt(h);
u=expf(log10(t)/log10(e)/3);

if (h>0) w=1;
if (h<=0) w=3;
if ((f==0) && (g==0) && (h==0)) w=2;

float **eigenValues = createNewMatrix(0, 3, 1);

eigenValues[0][0]=2*j*cos(k/3)-(b/3*a);
eigenValues[1][0]=l*(m+n)+p;
eigenValues[2][0]=l*(m-n)+p;

printMatrix(eigenValues, 3, 1);
return eigenValues;


}


float **ArrayElipseMatrix(float *Array, int samples)
{
	float **ElipseM = createNewMatrix(0, samples, 5);
	float **ArrayM = ArraytoMatrix(Array, 2, samples);
	float **test = TransposeMatrix(ArrayM, 2, samples);
	printMatrix(test, samples, 2);
	printf("hi\r\r");

	   for (uint16_t c = 0 ; c < samples; c++)
	                {
		   	   	   	   	    ElipseM[c][0]= ArrayM[1][c]*ArrayM[1][c];
		   	   	   	   	    ElipseM[c][1]= 2*ArrayM[1][c]*ArrayM[0][c];
		   	   	   	   	   	ElipseM[c][2]= ArrayM[0][c]*ArrayM[0][c];
		   	   	   			ElipseM[c][3]= 2*ArrayM[1][c];
		   	   				ElipseM[c][4]= 2*ArrayM[0][c];
	                }


	   return ElipseM;
}

float ** ArrayElipsoidMat (float **ArrayM, int samples)
{

	//x2 y2 z2 2xy 2xz 2yz 2x 2y 2z

	float **ElipseoidM = createNewMatrix(0, samples, 9);

		for (uint16_t c = 0 ; c < samples; c++)
					{
							ElipseoidM[c][0]= ArrayM[0][c]*ArrayM[0][c];
							ElipseoidM[c][1]= ArrayM[1][c]*ArrayM[1][c];
							ElipseoidM[c][2]= ArrayM[2][c]*ArrayM[2][c];

							ElipseoidM[c][3]= 2*ArrayM[0][c]*ArrayM[1][c];
							ElipseoidM[c][4]= 2*ArrayM[0][c]*ArrayM[2][c];
							ElipseoidM[c][5]= 2*ArrayM[1][c]*ArrayM[2][c];

							ElipseoidM[c][6]= 2*ArrayM[0][c];
							ElipseoidM[c][7]= 2*ArrayM[1][c];
							ElipseoidM[c][8]= 2*ArrayM[2][c];
			        }
		printMatrix(ElipseoidM, samples, 9);
			   return ElipseoidM;
}



float ** generateElispoidData(void)
{
#define Samples 30
float** data = createNewMatrix(0, 3, Samples);
//XY
float XStretch = 0.9,YStretch = 1.1,ZStretch = 1.2;

float Xoff= 0.07,Yoff = 0.08,Zoff = -0.02;
int a;
for ( a = 0 ; a < (Samples/3); a++)//Z = 0
     {
    	data[0][a] = (XStretch*(sinf(((2*M_PI)/(Samples/3))*(a)))) + Xoff;//x
    	data[1][a] = (YStretch*(cosf(((2*M_PI)/(Samples/3))*(a+1)))) + Yoff ;//y
    	data[2][a] = Zoff;//z
     }
for ( a = 0 ; a < (Samples/3); a++)//Z = 0
     {
    	data[0][a+(Samples/3)] = (XStretch*(sinf(((2*M_PI)/(Samples/3))*(a)))) + Xoff ;//x
    	data[1][a+(Samples/3)] = Yoff;
    	data[2][a+(Samples/3)] = (ZStretch*(cosf(((2*M_PI)/(Samples/3))*(a)))) + Zoff ;//z
     }
for ( a = 0 ; a < (Samples/3); a++)//Z = 0
     {
    	data[0][a+(2*Samples/3)] = Xoff;//x
    	data[1][a+(2*Samples/3)] = (YStretch*(sinf(((2*M_PI)/(Samples/3))*(a+1)))) + Yoff ;//y
    	data[2][a+(2*Samples/3)] = (ZStretch*(cosf(((2*M_PI)/(Samples/3))*(a)))) + Zoff ;
     }


//fit elipseoid
printMatrix(data, 3, Samples);
float **M = ArrayElipsoidMat(data, 30);
float **MT = TransposeMatrix(M, Samples, 9);
float **MTM = MultiplyMatrix(MT, 9, Samples, M, Samples, 9);




float **INVMTM = inverseMatrix(MTM, 9);
float **MTMINVMT = MultiplyMatrix(INVMTM, 9, 9, MT, 9, Samples);
float **OnesColumVector = createNewMatrix(1, Samples, 1);
float **coefs = MultiplyMatrix(MTMINVMT, 9, Samples, OnesColumVector, Samples, 1);
printMatrix(INVMTM, 9, 9);
printMatrix(coefs, 9, 1);

//center elipseoid
float** A = createNewMatrix(0, 4, 4);
	A[0][0] = coefs[0][0];
	A[0][1] = coefs[3][0];
	A[0][2] = coefs[4][0];
	A[0][3] = coefs[6][0];

	A[1][0] = coefs[3][0];
	A[1][1] = coefs[1][0];
	A[1][2] = coefs[5][0];
	A[1][3] = coefs[7][0];

	A[2][0] = coefs[4][0];
	A[2][1] = coefs[5][0];
	A[2][2] = coefs[2][0];
	A[2][3] = coefs[8][0];

	A[3][0] = coefs[6][0];
	A[3][1] = coefs[7][0];
	A[3][2] = coefs[8][0];
	A[3][3] = -1;

	printMatrix(A, 4, 4);

	float **aINV = inverseMatrix(A, 3);
	float **NEGaINV = NegitiveMatrix(aINV, 3, 3);
	printMatrix(NEGaINV, 3, 3);

	float **CVector = createNewMatrix(0, 3, 1);
	CVector[0][0] = coefs[6][0];
	CVector[1][0] = coefs[7][0];
	CVector[2][0] = coefs[8][0];
	printMatrix(CVector, 3, 1);

	//get offsets
	float ** ofs = MultiplyMatrix(NEGaINV, 3, 3, CVector, 3, 1);
	printMatrix(ofs, 3, 1);

	float **Tmtx = createNewMatrix(3,4,4);
	Tmtx[3][0] = ofs[0][0];
	Tmtx[3][1] = ofs[1][0];
	Tmtx[3][2] = ofs[2][0];
	printMatrix(Tmtx, 4, 4);

	float ** Tmtxa = MultiplyMatrix(Tmtx, 4, 4, A, 4, 4);
	float ** tTmtx = TransposeMatrix(Tmtx, 4, 4);

	// transposed elipseoid:
	float ** AT = MultiplyMatrix(Tmtxa, 4, 4, tTmtx, 4, 4);
	printMatrix(AT, 4, 4);

	// normalise
	float **NAT = normaliseMatrix(AT, 4, 4);
	printMatrix(NAT, 3, 3);

	EigenValues3x3(NAT);



return data;
}

float ** getElipsoicoeff (void)
{
	float ** elipsoid = generateElispoidData();



}
float **elipsecoef (void)
{
#define SAMPLES 15

	float data[2][SAMPLES];


    uint16_t r,c;
     for (uint16_t a = 0 ; a < (SAMPLES); a++)
     {
    	 data[0][a] = 1.3*(cosf((2*M_PI/(SAMPLES))*(a)));//y
    	data[1][a] = 1.2*(sinf((2*M_PI/(SAMPLES))*(a+1)));//x
     }

     float **M = ArrayElipseMatrix(data, SAMPLES);
     float **MT = TransposeMatrix(M, SAMPLES, 5);
     float **MTM = MultiplyMatrix(MT, 5, SAMPLES, M, SAMPLES, 5);
     float **MTMINV = inverseMatrix(MTM, 5);
     float **MTMINVMT = MultiplyMatrix(MTMINV, 5, 5, MT, 5, SAMPLES);
     float **OnesColumVector = createNewMatrix(1, SAMPLES, 1);
     float **coefs = MultiplyMatrix(MTMINVMT, 5, SAMPLES, OnesColumVector, SAMPLES, 1);


if (1)
{
     printArray(data, 2, SAMPLES);
     printMatrix(M, SAMPLES, 5);
     printMatrix(MT, 5, SAMPLES);
     printMatrix(MTM, 5, 5);
      printMatrix(MTMINV, 5, 5);
      printMatrix(MTMINVMT, 5, SAMPLES);
      printMatrix(OnesColumVector, SAMPLES, 1);
     printMatrix(coefs, 5, 1);
}

return coefs;

}

float **EigenValues3x3 (float **Matrix)
{
	static float a,b,c,d;

	a = -1;

	b = Matrix[0][0] + Matrix[1][1] + Matrix[2][2];

	c = - (Matrix[2][2]*Matrix[0][0])
		- (Matrix[2][2]*Matrix[1][1])
		- (Matrix[0][0]*Matrix[1][1])
		+ (Matrix[0][2]*Matrix[2][0])
		+ (Matrix[1][2]*Matrix[2][1])
		+ (Matrix[0][1]*Matrix[1][0]);

	d = +(Matrix[2][2]*Matrix[0][0]*Matrix[1][1])
		+(Matrix[0][2]*Matrix[2][0]*Matrix[1][2])
		+(Matrix[0][2]*Matrix[1][0]*Matrix[2][1])
		-(Matrix[0][2]*Matrix[1][1]*Matrix[2][0])
		-(Matrix[0][0]*Matrix[1][2]*Matrix[2][1])
		-(Matrix[0][1]*Matrix[1][0]*Matrix[2][2]);

	printf("\r\n%f\t%f\t%f\t%f\t\r\n",a,b,c,d);

	float **roots = get_cubic(a,b,c,d);
	return roots;
}


float **quadraticSolveRoots (float a, float b, float c)
{

	float **roots = createNewMatrix(0, 2, 1);
	roots[0][0] = (-b - sqrtf((b*b)-4*a*c))/2*a;
	roots[1][0] = (-b + sqrtf((b*b)-4*a*c))/2*a;

	return roots;
}

float **EigenValues2x2 (float **Matrix)
{
	float a,b,c;
	float **roots = createNewMatrix(0, 2, 1);
	a = 1;
	b = -Matrix[0][0] - Matrix[1][1];
	c = Matrix[0][0]*Matrix[1][1]-Matrix[1][0]*Matrix[0][1];

	roots = quadraticSolveRoots(a,b,c);
	return roots;
}

float **GetEigenVelctors2x2 (float **EValues, float **Matrix)
{
	int i;
	float **Eigenvectors = createNewMatrix(0, 2, 2);
	float Norm;

	Eigenvectors[1][0] = 1;
	Eigenvectors[1][1] = 1;

	for (i=0;i<2;i++)
	{


		if (fabs(Matrix[0][1]) < 0.01)
		{

			printf("\r\n e \r\n");

			Eigenvectors[0][0] = 1;
			Eigenvectors[0][1] = 0;
			Eigenvectors[1][0] = 0;
			Eigenvectors[1][1] = 1;
		}
		else
		{
		Eigenvectors[0][i] = -Matrix[0][1]/(Matrix[0][0]-EValues[i][0]);
		//Eigenvectors[0][i] = (Matrix[0][0]-EValues[i][0]) > 00.1 ? 1 : Eigenvectors[0][i];
		Norm = sqrtf(Eigenvectors[0][i]*Eigenvectors[0][i] + Eigenvectors[1][i]*Eigenvectors[1][i]);
		Eigenvectors[0][i] = Eigenvectors[0][i]/Norm;
		Eigenvectors[1][i] = Eigenvectors[1][i]/Norm;
		}

		//Eigenvectors[0][0] = 1;
		//Eigenvectors[1][0] = 0;

	}



	return Eigenvectors;
}


void centerElipse (float **coefs)
{
	float **a = createNewMatrix(0, 3, 3);
	float **aINV = createNewMatrix(0, 3, 3);
	float **NEGaINV = createNewMatrix(0, 3, 3);
	float **CVector = createNewMatrix(0, 2, 1);
	float **ofs = createNewMatrix(0,2,1);
	float **Tmtx = createNewMatrix(3,3,3);
	float **AT = createNewMatrix(0,3,3);


	a[0][0] = coefs[0][0];
	a[0][1] = coefs[1][0];
	a[0][2] = coefs[3][0];

	a[1][0] = coefs[1][0];
	a[1][1] = coefs[2][0];
	a[1][2] = coefs[4][0];

	a[2][0] = coefs[3][0];
	a[2][1] = coefs[4][0];
	a[2][2] = -1;

	printMatrix(a, 3, 3);
	printMatrix(a, 2, 2);

	aINV = inverseMatrix(a, 2);
	printMatrix(aINV, 2, 2);

	NEGaINV = NegitiveMatrix(aINV, 2, 2);
	printMatrix(NEGaINV, 2, 2);

	CVector[0][0] = coefs[3][0];
	CVector[1][0] = coefs[4][0];

	printMatrix(CVector, 2, 1);

	ofs = MultiplyMatrix(NEGaINV, 2, 2, CVector, 2, 1);

	printMatrix(ofs, 2, 1);
	Tmtx[2][0] = ofs[0][0];
	Tmtx[2][1] = ofs[1][0];
	printMatrix(Tmtx, 3, 3);

	float **Tmtxa = createNewMatrix(0,3,3);
	float **tTmtx = createNewMatrix(0,3,3);


	Tmtxa = MultiplyMatrix(Tmtx, 3, 3, a, 3, 3);

	tTmtx = TransposeMatrix(Tmtx, 3, 3);

	AT = MultiplyMatrix(Tmtxa, 3, 3, tTmtx, 3, 3);
	printMatrix(AT, 3, 3);
	float ** NAT = createNewMatrix(0, 3, 3);
	NAT = normaliseMatrix(AT, 3, 3);
	printMatrix(NAT, 3, 3);

	float **eigenValues = createNewMatrix(0, 2, 1);
	eigenValues = EigenValues2x2(NAT);
	printMatrix(eigenValues, 2, 1);


	float **eigenVec = createNewMatrix(0, 2, 2);
	eigenVec = GetEigenVelctors2x2(eigenValues, NAT);
	printMatrix(eigenVec, 2, 2);

	float **Rt = createNewMatrix(0, 2, 2);
	Rt = TransposeMatrix(eigenVec, 2, 2);

	float **sqrt1 = createNewMatrix(0, 2, 2);
	sqrt1[0][0] = sqrtf(eigenValues[1][0]);
	sqrt1[1][1] = sqrtf(eigenValues[0][0]);
	printMatrix(sqrt1, 2, 2);

	float **rinv = createNewMatrix(0, 2, 2);
	rinv = inverseMatrix(eigenVec, 2);
	//T = inverseMatrix(T, 2);
	printMatrix(rinv, 2, 2);


	float **T = createNewMatrix(0, 2, 2);
	T = MultiplyMatrix(rinv, 2, 2, sqrt1, 2, 2);
		printMatrix(T, 2, 2);

		if (detrimant123>0)
		{
			T[0][0] = -T[0][0];
			T[0][1] = -T[0][1];
				printf ("flipping\r\n");
		}

	float dataa[2][SAMPLES];
	float data[2][SAMPLES];
	float circle[2][SAMPLES];
	    uint16_t r,c;
	     for (uint16_t a = 0 ; a < (SAMPLES); a++)
	     {
	    	 data[0][a] = 1.3*(cosf((2*M_PI/(SAMPLES))*(a+1)));//y
	    	   data[1][a] = 1.2*(sinf((2*M_PI/(SAMPLES))*(a)));//x
	     }
	     printArray(data, 2, SAMPLES);

	     for (uint16_t a = 0 ; a < (SAMPLES); a++)
	     	     {
	    	 circle[0][a] = (cosf((2*M_PI/(SAMPLES))*(a+1)));//y
	    	 circle[1][a] = (sinf((2*M_PI/(SAMPLES))*(a)));//x
	     	     }
	     	     printArray(circle, 2, SAMPLES);

	     	    for (uint16_t a = 0 ; a < (SAMPLES); a++)
	     	    	     {
	     	    	dataa[0][a] = data[0][a]*T[0][0] + data[1][a]*T[1][0];//y
	     	    	dataa[1][a] = data[0][a]*T[0][1] + data[1][a]*T[1][1];//x
	     	    	     }
	     	    	     //printArray(data, 2, SAMPLES);

	     printArray(dataa, 2, SAMPLES);

}




float **inverseMatrix(float ** Matrix,int size) // max size = 5
{
	float a[25][25];

	//a = MatrixtoArray(Matrix, 5, 5,a);
	MatrixtoArrayB(Matrix, size, size,a);

	printArray(a, size, size);

	float d = detrminant(a, size);
	printf("detrimant is %f\r\n",d);
	detrimant123 =d;
	if (d != 0)
	{
		float **INV = cofactors(a, size);
		return INV;
	}
	else
		return NULL;

}






float detrminant(float a[9][9], float k) {

	float s = 1, det = 0;

	float b[9][9];


	int i, j, m, n, c;




	if (k == 1) {
		return (a[0][0]);
	}
	else
		{
		det = 0;
		for (c = 0; c < k; c++)
			{
			m = 0;
			n = 0;
			for (i = 0; i < k; i++)
				{
				for (j = 0; j < k; j++)
					{
					b[i][j] = 0;
					if (i != 0 && j != c)
						{
						b[m][n] = a[i][j];
						if (n < (k - 2))
							{
						    	n++;
							}
						 else
						 	 {
							n = 0;
							m++;
						 	 }
						}
					}
				}

			det = det + s * (a[0][c] * detrminant(b, k - 1));
			s = -1 * s;
		}
	}
	return det;
}

float **cofactors(float num[9][9], float f) {



	//printf("\r\n num: \r\n");
	//printArray(num, 5, 5);
	float b[9][9], fac[9][9];
	float c[9][9], inv[9][9], d;

	int p, q, m, n, i, j;
	for (q = 0; q < f; q++) {
		for (p = 0; p < f; p++) {
			m = 0;
			n = 0;
			for (i = 0; i < f; i++) {
				for (j = 0; j < f; j++) {
					b[i][j] = 0;
					if (i != q && j != p) {
						b[m][n] = num[i][j];
						if (n < (f - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = pow(-1, q + p) * detrminant(b, f - 1);
		}
	}

	//printArray(fac, 5, 5);
	for (i = 0; i < f; i++) {
		for (j = 0; j < f; j++) {
			c[i][j] = fac[j][i];
		}
	}

	d = detrminant(num, f);

	//inv[i][j] = 0;
	for (i = 0; i < f; i++) {
		for (j = 0; j < f; j++) {
			inv[i][j] = c[i][j] / d;
		}
	}
//printf("\r\n invers is: \r\n");
	//printArray(inv, 5, 5);
	float **INVS = ArraytoMatrix(inv, 9, 9);
	return INVS;

}
