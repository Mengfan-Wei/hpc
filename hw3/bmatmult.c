#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define FLOAT_TYPE double
#define n 512

void do_block(int b, FLOAT_TYPE* A, FLOAT_TYPE* B, FLOAT_TYPE* C)
{
	int i, j, k;
	for (i = 0; i < b; i++)
	{ 
		for (j = 0; j < b; j++)
		{
			double Cij = C[i*n+j];
			for (k = 0; k < b; k++)
				Cij += A[i*n+k] * B[k*n+j];
			C[i*n+j] = Cij;	
		}
	}
}

int main(void)
{
	double timeinsec;
	clock_t start, end;
	int i, j, l, k, m, q;
	int b, LOOP_COUNT = 5;
	FLOAT_TYPE* A;
	FLOAT_TYPE* B;
	FLOAT_TYPE* C;

	m = q = n;

	//Partition heap space to store matrices A, B, and C
	A = (FLOAT_TYPE*)malloc(sizeof(double)*m*q);
	B = (FLOAT_TYPE*)malloc(sizeof(double)*q*n);
	C = (FLOAT_TYPE*)malloc(sizeof(double)*m*n);


	//Generate random entries in [0, 50] for matrices A and B
	//srand((size_t)time(NULL)); 
	for (i = 0; i < m*q; i++)
			A[i] = rand() % 51;
	for (i = 0; i < q*n; i++)
			B[i] = rand() % 51;

	//Compute the product of A and B
	for (b=4; b<512; b*=2){
		start = clock();
		for (l = 0; l < LOOP_COUNT; l++)
			for (i = 0; i < n; i+=b)
				for (j = 0; j < n; j+=b)
					for (k = 0; k < n; k+=b)
						do_block(b, A+j*n+k, B+k*n+i, C+j*n+i);	
		end = clock();

		timeinsec = (double)(end-start) / CLOCKS_PER_SEC / LOOP_COUNT;
		double gflop = (2.0*m*q*n + 1.0*m*n)*1E-9;
		printf("When b = %d\n", b);
		printf("Average time: %e sec\n", timeinsec);
		printf("GFlops	    : %.5f\n", gflop);
		printf("GFlops/sec  : %.5f\n", gflop/timeinsec);
		printf("\n");
	}
	return 0;
}
