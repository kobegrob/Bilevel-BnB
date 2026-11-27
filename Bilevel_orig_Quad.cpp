

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
double totalSolveP1Time = 0.0;
double totalOmegaTime = 0.0;
clock_t global_start;
int INF = 1000000000;
double TIME_LIMIT = 3600.0; // seconds
bool time_limit_reached = false;

// =====================================================================================================
// PERFORMANCE FIX 1 (2026-11-23): throttled time-limit check inside the DP recursions
// -----------------------------------------------------------------------------------------------------
// Problem: solveP1(), omega_dense() and omega_sparse() used to call clock() at the start of EVERY call
//   to test the 3600 s time limit. clock() is a system call (~180 ns on the test machine), and these
//   recursions are called millions of times: the number of calls grows with C^L * (u - l), i.e. roughly
//   quadratically when capacities and item quantities grow. Measured on n_25/new252 with all capacities
//   x200 (C^U = 2000, C^L = 2200): 9.0 s with the per-call clock(), 0.29 s without it -- the clock()
//   calls were ~97% of the run time.
// Fix: the recursions now call dp_time_limit_exceeded(), which reads the clock only once every
//   TIME_CHECK_INTERVAL calls. The DP results, the search tree and the optimal values are unchanged.
// Side effect: the time limit can be overshot by at most TIME_CHECK_INTERVAL DP calls (well under a
//   millisecond at typical speeds). The per-node checks in bnb() still read the clock directly.
// =====================================================================================================
const unsigned long long TIME_CHECK_INTERVAL = 65536; // must be a power of two (used as a bit mask)
unsigned long long dp_time_check_counter = 0;          // DP calls since the start of the run
inline bool dp_time_limit_exceeded() {
	if ((++dp_time_check_counter & (TIME_CHECK_INTERVAL - 1)) != 0) return false; // not time to check yet
	return (double)(clock() - global_start) / CLOCKS_PER_SEC >= TIME_LIMIT;
}

// =====================================================================================================
// PERFORMANCE FIX 2 (2026-11-23): RAM-aware choice between the dense and the sparse omega table
// -----------------------------------------------------------------------------------------------------
// Problem: procedure() allocates the dense omega table, (n+1)*(C^U+1)*(C^L+1)*2 ints, whenever it is
//   below a FIXED 10 GiB (the old comment said 8 GB), whatever the machine's memory. The table grows with
//   C^U * C^L, so on a machine with less free RAM the program crashed with std::bad_alloc (or was killed)
//   before it started searching. Measured: n_25/new252 with capacities x400 needs a 3.4 GiB table and
//   crashed on a 4 GB machine, although fewer than 0.2% of the table's cells are ever used.
// Fix: the dense table is only allowed up to DENSE_OMEGA_RAM_FRACTION of the machine's physical RAM
//   (and never above the old 10 GiB). Above that, the existing sparse (hash map) omega is used, which
//   stores only the states actually visited: the same instance then runs in ~230 MB. The sparse table is
//   about 2-3x slower per lookup than the dense one, so the dense table is still used whenever it fits.
//   The threshold < 0.0025 rule for choosing the sparse table is unchanged.
// =====================================================================================================
#include <unistd.h> // sysconf(): physical RAM size, works on Linux and macOS
const double DENSE_OMEGA_RAM_FRACTION = 0.5;                    // share of physical RAM the dense table may use
const size_t OLD_DENSE_OMEGA_LIMIT = (size_t)10 * 1024 * 1024 * 1024; // previous fixed limit (10 GiB)
// Refinement (2026-11-23): the memory that counts is the smaller of the physical RAM and the process's own
// address-space limit (ulimit -v, as set by job schedulers); without it, a run with such a limit could still pick
// a dense table that the limit forbids and stop at once with std::bad_alloc.
#include <sys/resource.h> // getrlimit(): the process's address-space limit, Linux and macOS
// Memory available to this process: the smaller of the physical RAM and the address-space limit (ulimit -v),
// or 0 if unknown. Used by PERFORMANCE FIX 7 (budget of the omega row cache).
inline double available_memory_bytes() {
	long pages = sysconf(_SC_PHYS_PAGES), page_size = sysconf(_SC_PAGE_SIZE);
	double memory = (pages > 0 && page_size > 0) ? (double)pages * (double)page_size : 0.0;
	struct rlimit rl;
	if (getrlimit(RLIMIT_AS, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY && (memory == 0.0 || (double)rl.rlim_cur < memory))
		memory = (double)rl.rlim_cur;
	return memory;
}
inline size_t dense_omega_byte_limit() {
	long pages = sysconf(_SC_PHYS_PAGES), page_size = sysconf(_SC_PAGE_SIZE);
	double memory = (pages > 0 && page_size > 0) ? (double)pages * (double)page_size : 0.0; // 0 = unknown
	struct rlimit rl;
	if (getrlimit(RLIMIT_AS, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY && (memory == 0.0 || (double)rl.rlim_cur < memory))
		memory = (double)rl.rlim_cur;
	if (memory == 0.0) return OLD_DENSE_OMEGA_LIMIT; // memory size unknown: keep the old behaviour
	size_t ram_limit = (size_t)(memory * DENSE_OMEGA_RAM_FRACTION);
	return std::min(ram_limit, OLD_DENSE_OMEGA_LIMIT);
}
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
		printf("Item %d: %d\n", i, z[i]); j += z[i]; k += c[i] + z[i]*z[i]*variable[i]; break;
	}
	while (++i <= n) if (y[i])
	{
		printf("Item %d: %d\n", i, z[i]); j += z[i]; k += c[i] + z[i]* z[i]* variable[i];
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
