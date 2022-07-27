#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mkl.h>
#define FLOAT_TYPE double

int main()
{
	double timeinsec;
	clock_t start, end;
	int m, n, k, i, LOOP_COUNT = 5;

	FLOAT_TYPE* A;
	FLOAT_TYPE* B;
	FLOAT_TYPE* C;
	FLOAT_TYPE alpha = 1.0, beta = 0;
	char transa = 'N';
	char transb = 'N';

	for (n = 128; n < 8193; n*=2)
	{
		m = k = n;
		//Partition heap space to store matrices A, B, and C
		A = (FLOAT_TYPE*)malloc(sizeof(double) * m*k);
		B = (FLOAT_TYPE*)malloc(sizeof(double) * k*n);
		C = (FLOAT_TYPE*)malloc(sizeof(double) * m*n);


		//Generate random entries in [0, 50] for matrices A and B
		//srand((size_t)time(NULL)); 
		for (i = 0; i < m*k; i++)
			A[i] = rand() % 51;
		for (i = 0; i < k*n; i++)
			B[i] = rand() % 51;

		dgemm(&transa, &transb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);

		//Compute the product of A and B
		start = clock();
		for (i = 0; i < LOOP_COUNT; i++)
			dgemm(&transa, &transb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
		end = clock();

		timeinsec = (double)(end - start) / CLOCKS_PER_SEC / LOOP_COUNT;
		double gflop = (2.0*m*k*n + 1.0*m*n)*1E-9;
		
		printf("When n = %d\n", n);
		printf("Average time: %e sec\n", timeinsec);
		printf("GFlops	    : %.5f\n", gflop);
		printf("GFlops/sec  : %.5f\n", gflop / timeinsec);
		printf("\n");
	}
	return 0;
}

