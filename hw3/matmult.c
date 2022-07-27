#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FLOAT_TYPE double
int main(void)
{
	double timeinsec;
	clock_t start, end;
	int i, j, l, k, m, q, n, sum;
	int LOOP_COUNT = 5;
	FLOAT_TYPE** a;
	FLOAT_TYPE** b;
	FLOAT_TYPE** c;

	for (n = 128; n < 2049; n*=2)
	{
		m = q = n;
		//Partition heap space to store matrices A, B, and C
		a = (FLOAT_TYPE**)malloc(sizeof(double) * m);
		for (i = 0; i < m; i++)
			a[i] = (FLOAT_TYPE*)malloc(sizeof(double) * q);

		b = (FLOAT_TYPE**)malloc(sizeof(double) * q);
		for (i = 0; i < q; i++)
			b[i] = (FLOAT_TYPE*)malloc(sizeof(double) * n);

		c = (FLOAT_TYPE**)malloc(sizeof(double) * m);
		for (i = 0; i < m; i++)
			c[i] = (FLOAT_TYPE*)malloc(sizeof(double) * n);

		//Generate random entries in [0, 50] for matrices A and B
		//srand((size_t)time(NULL)); 
		for (i = 0; i < m; i++)
			for (j = 0; j < q; j++)
				a[i][j] = rand() % 51;
		for (i = 0; i < q; i++)
			for (j = 0; j < n; j++)
				b[i][j] = rand() % 51;

		//Compute the product of A and B
		start = clock();
		for (l = 0; l < LOOP_COUNT; l++)
			for (i = 0; i < m; i++)
				for (j = 0; j < n; j++)
				{
					sum = 0;
					for (k = 0; k < q; k++)
						sum += a[i][k] * b[k][j];
					c[i][j] = sum;
				}	
		end = clock();

		timeinsec = (double)(end-start) / CLOCKS_PER_SEC / LOOP_COUNT;
		double gflop = (2.0*m*q*n + 1.0*m*n)*1E-9;
		printf("When n = %d\n", n);
		printf("Average time: %e sec\n", timeinsec);
		printf("GFlops      : %.5f\n", gflop);
		printf("GFlops/sec  : %.5f\n", gflop/timeinsec);
		printf("\n");
	}
	return 0;
}
