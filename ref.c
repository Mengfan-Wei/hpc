#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
int main()
{
	float dx, dt, t, L = 1.0, a[100000], b[100000], k = 1.0;
	/*dx��ʾdelta(x)��dt��ʾdelta(t)��kΪ����ɢ�ʡ����������С������ֵ���仯������ȡ�ñȽϴ�*/
	int i, n1, n2, n = 0, total = 1;
	FILE *F;//�ļ�ָ��
	F = fopen("data.txt", "w");
	dx = 0.01;
	dt = 0.001;
	t = 10;
	n1 = (int)(L / dx);
	n2 = (int)(t / dt);
	//printf("n1=%d  n2=%d\n",n1,n2);
	for (i = 0; i < n1; i++) {     /*��ʼ��ֵ*/
		if (0 < i < n1) {
			a[i] = exp(i*dx);
		}
		else a[i] = 20;
		fprintf(F, "%8.4f", a[i]);//д��t=0ʱ�̵��¶�
	}
	fprintf(F, "\n");    //����
	while (n < n2) {    /*ʱ�䵽�ˣ�ֹͣѭ��*/
		for (i = 0; i < n1; i++) {
			/*����*/
			if (i == 0) {
				b[i] = 20;
				a[i] = b[i];
			}     /*��ʾ�������¶�ʼ�ձ�����20���϶�*/
			//else if (i == n1 - 1) {
			//	b[i] = 20;
			//	a[i] = b[i];
			//}     /*��ʾ�ұ�����¶�ʼ�ձ�����20���϶�*/
			else {
				b[i] = a[i] + k * dt*(a[i + 1] - 2 * a[i] + a[i - 1]) / (dx*dx) + dt*sin(i*dx);   /*��ַ��̵ı��ʽ*/
				a[i] = b[i];
			}
			n++;
			fprintf(F, "%8.4f ", a[i]);
			//printf("%8.4f",a[i]);
			if (total%n1 == 0) {             /*���data�ļ�*/
				//printf("\n");
				fprintf(F, "\n");
			}
			total++;
		}
	}
	fclose(F);
	return 0;
}