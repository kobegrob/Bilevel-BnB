

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <cstddef>

#define maxn 10000
#define maxb 1500000


typedef int mat[maxn+1];

int gub, glb,n, b1, b2, nodes;
double avg_a1, avg_width, avg_Prod_Low;// Where n is the number of items, b1 is the size of the upper level knapsack and b2 the size of the lower level knapsack
mat c, a1, Prod_Up, Prod_Low, variable, x, y, z, cvariable; // C is the fixed value of an item, a1 is the weight in the upper level, Prod_Up and Prod_Low
// are the upper and lower bound of an item in the lower level, variable represents the variable value of an item and x represents upper level items
// y represents what lower level items are selected and z how much of the selected items is selected.
// Reading the problem
clock_t global_start;
int INF = 1000000000;
double TIME_LIMIT = 900.0; // seconds
bool time_limit_reached = false;
void input (void) {
	int i; char filename[20]; FILE *fp;
	printf ("\n\n	Give the problem name: "); 
	scanf ("%s", filename); strcat (filename, ".knp");
	fp = fopen (filename, "r"); fscanf (fp, "%d %d %d", &n, &b1, &b2);
	if (n > maxn) {printf ("Too many objects\n\n"); return;}
	if (b1 > maxb || b2 > maxb) { printf("Capacities are too large\n\n"); return; }
	for (i = 0; ++i <= n;) fscanf(fp, "%d %d %d %d %d", a1 + i, c + i, Prod_Up + i, Prod_Low + i, variable + i);
	fclose (fp); return;
}

void output (void) {
	if (time_limit_reached) {
		printf("\nTime limit reached (3600s). Best solution found is not guaranteed optimal.\n");
	}
	int i, j, k;
	printf ("\n	Upper Level Items:\n "); 
	for (i = j = 0; ++i <= n;) if (x[i]) 
	{	printf ("%d, ", i); j += a1[i]; break;}
	while (++i <= n) if (x[i])
	{	printf ("%4d, ", i); j += a1[i];} 
	printf("\n	Lower Level = %d <= %d\n", j, b1);
	printf("\n	Lower Level Items:\n ");
	for (i = j = k = 0; ++i <= n;) if (y[i])
	{
		printf("Item %d: %d\n", i, z[i]); j += z[i]; k += c[i] + z[i]*variable[i]; break;
	}
	while (++i <= n) if (y[i])
	{
		printf("Item %d: %d\n", i, z[i]); j += z[i]; k += c[i] + z[i] * variable[i];
	}
	printf ("\n	Lower Level = %d <= %d\n", j, b2);
	printf ("\n	Value = %d\n", k);
	printf("\n\tBest feasible value = %d\n", gub);
	printf("\n\tBest valid lower bound = %d\n", glb);

	if (gub < INF) {
    	double gap = 100.0 * ((double)gub - (double)glb) / std::max(1.0, fabs((double)gub));
    	if (gap < 0.0) gap = 0.0;
    	printf("\n\tGap = %.4f%%\n", gap);
	}
	return;
}

void procedure (void);

int main () {
	double elapsed_time; clock_t start_time;
	nodes = 0;
	input ();
	start_time= clock();
	procedure ();
	elapsed_time = (float) (clock() - start_time)/ CLOCKS_PER_SEC;
	output ();
	printf ("\n	Computation time: %.2f seconds\n", elapsed_time);
	printf ("\n	Nodes Visited: %d \n", nodes);
	return 0;
}
