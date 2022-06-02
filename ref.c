#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
int main()
{
	float dx, dt, t, L = 1.0, a[100000], b[100000], k = 1.0;
	/*dx表示delta(x)，dt表示delta(t)，k为热扩散率。由于数组大小随输入值而变化，所以取得比较大*/
	int i, n1, n2, n = 0, total = 1;
	FILE *F;//文件指针
	F = fopen("data.txt", "w");
	dx = 0.01;
	dt = 0.001;
	t = 10;
	n1 = (int)(L / dx);
	n2 = (int)(t / dt);
	//printf("n1=%d  n2=%d\n",n1,n2);
	for (i = 0; i < n1; i++) {     /*初始赋值*/
		if (0 < i < n1) {
			a[i] = exp(i*dx);
		}
		else a[i] = 20;
		fprintf(F, "%8.4f", a[i]);//写入t=0时刻的温度
	}
	fprintf(F, "\n");    //换行
	while (n < n2) {    /*时间到了，停止循环*/
		for (i = 0; i < n1; i++) {
			/*计算*/
			if (i == 0) {
				b[i] = 20;
				a[i] = b[i];
			}     /*表示左壁面的温度始终保持在20摄氏度*/
			//else if (i == n1 - 1) {
			//	b[i] = 20;
			//	a[i] = b[i];
			//}     /*表示右壁面的温度始终保持在20摄氏度*/
			else {
				b[i] = a[i] + k * dt*(a[i + 1] - 2 * a[i] + a[i - 1]) / (dx*dx) + dt*sin(i*dx);   /*差分方程的表达式*/
				a[i] = b[i];
			}
			n++;
			fprintf(F, "%8.4f ", a[i]);
			//printf("%8.4f",a[i]);
			if (total%n1 == 0) {             /*输出data文件*/
				//printf("\n");
				fprintf(F, "\n");
			}
			total++;
		}
	}
	fclose(F);
	return 0;
}