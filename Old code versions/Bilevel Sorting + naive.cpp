#include "Bilevel_orig.cpp"
#include <vector>
using namespace std;
int gub, glb, llb, lub, fractionindex; //These variables represent the global and local bounds as well as where the first fractional value is 
mat nextinline,sorteda,sortedc, initx; //These arrays are used for the initialization of the BnB algorithm 
double ratio[maxn + 1], cratio[maxn + 1]; //This arrays contains the ratio of value/weight
mat xv, yu, zu; //Being able to keep the current best solution saved somewhere
int INF = 1000000000;

vector<vector<vector<int>>> omg;// Note that the maximum size of omega is already predefined

vector<vector<double>> memoBellman;
vector<vector<int>> memoP1; // Store the values of the Bellman algorithm + the value of the p1 algorithm, good for the dynamic memory allocation
int choiceZ[maxn+1][maxb+1];

//Here two quicksort algorithms are presented 1 for the ratio array and 1 for the upper bound arrays 
int partition(double arr[], int arr2[], int low, int high)
{
	double pivot = arr[high];
	int pivotItem = arr2[high];

	int i = (low - 1);
	for (int j = low; j <= high - 1; j++) {
		int itemJ = arr2[j];
		if ((arr[j] < pivot) ||
			(arr[j] == pivot && variable[itemJ] < variable[pivotItem]) ||
			(arr[j] == pivot && variable[itemJ] == variable[pivotItem] && c[itemJ] < c[pivotItem])) {
			i++;
			swap(arr[i], arr[j]);
			swap(arr2[i], arr2[j]);
		}
	}
	swap(arr[i + 1], arr[high]);
	swap(arr2[i + 1], arr2[high]);
	return (i + 1);
}

int partition(int arr[], int arr2[], int low, int high)
{
	int pivot = arr[high];
	int pivotItem = arr2[high];

	int i = (low - 1);
	for (int j = low; j <= high - 1; j++) {
		int itemJ = arr2[j];
		if ((arr[j] < pivot) ||
			(arr[j] == pivot && variable[itemJ] < variable[pivotItem]) ||
			(arr[j] == pivot && variable[itemJ] == variable[pivotItem] && c[itemJ] < c[pivotItem])) {
			i++;
			swap(arr[i], arr[j]);
			swap(arr2[i], arr2[j]);
		}
	}
	swap(arr[i + 1], arr[high]);
	swap(arr2[i + 1], arr2[high]);
	return (i + 1);
}
void quickSort(double arr[], int arr2[], int low, int high)
{
	if (low < high)
	{
		int pi = partition(arr,arr2, low, high);
		quickSort(arr, arr2, low, pi - 1);
		quickSort(arr, arr2, pi + 1, high);
	}
}
void quickSort(int arr[], int arr2[], int low, int high)
{
	if (low < high)
	{
		int pi = partition(arr, arr2, low, high);
		quickSort(arr, arr2, low, pi - 1);
		quickSort(arr, arr2, pi + 1, high);
	}
}
//The Bellman recursion that also saves this optimal solution somewhere (note that we can simply use x[] as it is not used yet

double Bellman(int i, int CupLeft) {
	if (i == 0) {
		return 0;
	}

	// 1. Efficiency: Check if we've already solved this state
	if (memoBellman[i][CupLeft] != -1) {
		return memoBellman[i][CupLeft];
	}

	double result;
	if (a1[nextinline[i]] > CupLeft) {
		result = Bellman(i - 1, CupLeft);
	}
	else {
		double value1 = Bellman(i - 1, CupLeft);
		double value2 = Bellman(i - 1, CupLeft - a1[nextinline[i]]) + ratio[nextinline[i]];

		result = (value1 > value2) ? value1 : value2;
	}

	// Store and return
	return memoBellman[i][CupLeft] = result;
}
// Reconstruct the best solution found so far
void reconstructSolutionB (int n, int capacity) {
	int currentCap = capacity;
	for (int i = n; i > 0; i--) {
		// If the value at this state is different from the value without the item,
		// it means we MUST have included item i to get that result.
		if (Bellman(i, currentCap) != Bellman(i - 1, currentCap)) {
			x[nextinline[i]] = 1;
			currentCap -= a1[nextinline[i]];
		}
		else {
			x[nextinline[i]] = 0;
		}
	}
}

//Here we write a function that we can use in all previous works that is able to solve p1 given a set of already interdicted items

int solveP1(int i, int Cleft, int xu[]) {
	if (i == 0) return 0;
	if (memoP1[i][Cleft] != -1) return memoP1[i][Cleft];

	// Case 1: Item cannot be taken
	if (Prod_Low[nextinline[i]] > Cleft || xu[nextinline[i]] == 1) {
		choiceZ[i][Cleft] = 0; // 0 means not taken/impossible
		return memoP1[i][Cleft] = solveP1(i - 1, Cleft, xu);
	}

	// Case 2: Explore taking the item vs. not taking it
	// Start by assuming "not taking it" is best
	int maxVal = solveP1(i - 1, Cleft, xu);
	int bestW = 0;

	int minval = (Prod_Up[nextinline[i]] < Cleft) ? Prod_Up[nextinline[i]] : Cleft;
	for (int w = Prod_Low[nextinline[i]]; w <= minval; w++) {
		int currentVal = solveP1(i - 1, Cleft - w, xu) + c[nextinline[i]] + (w * variable[nextinline[i]]);
		if (currentVal > maxVal) {
			maxVal = currentVal;
			bestW = w; // Store the specific weight that worked best
		}
	}

	choiceZ[i][Cleft] = bestW;
	return memoP1[i][Cleft] = maxVal;
}


// Reconstruct the solution found so far, should be called whenever we put the solution somewhere
void reconstructSolutionP(int n, int Cleft) {
	int currentCap = Cleft;
	for (int i = n; i > 0; i--) {
		int optimalW = choiceZ[i][currentCap];

		if (optimalW > 0) {
			yu[nextinline[i]] = 1;
			zu[nextinline[i]] = optimalW;
			currentCap -= optimalW;
		}
		else {
			yu[nextinline[i]] = 0;
			zu[nextinline[i]] = 0;
		}
	}
}


// Initial greedy algorithm to compute the first upper bound

int greedyupper(void) {
	// Compute the optimal solution of the upper level knapsack problem using the previously calculated ratio as objective function value --> Might be best to use COMBO
	// or the Bellman recursion for simplicity

	// Using this solution compute the optimal solution of the lower level problem using some sort of DP
	Bellman(n, b1);
	reconstructSolutionB(n, b1);
	// Somehow save this initial solution somewhere as it might already be the optimal solution, this should happen in the procedure function somewhere
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= b2; j++) {
			memoP1[i][j] = -1;
		}
	}
	int tempvalue = solveP1(n, b2, x);
	reconstructSolutionP(n, b2);
	for (int i = 1; i <= n; i++) {
		// Make sure to save the current best solution
		y[i] = yu[i];
		z[i] = zu[i];
	}
	return tempvalue; // the upper bound value
}
// Omega which is used to compute a strong lower bound for the problem. 
// Note that the currently implemented version solves the lower level problem as exactly as possible which potentially takes more time. 

int omega(int i, int Cup, int Clow) {
	// 1. Boundary Checks (Keep these first)
	if (Cup < 0) return INF; // Using a safe high value for "infinity"
	if (Clow < 0) return -INF;
	if (i > n) return 0;

	// 2. THE FIX: Check if we've already solved this specific state
	if (omg[i][Cup][Clow] != -1) {
		return omg[i][Cup][Clow];
	}

	int result; // Local variable to hold the calculation


	if (Prod_Low[nextinline[i]] > Clow) {
		result = omega(i + 1, Cup, Clow);
		return omg[i][Cup][Clow] = result;
	}
	else {
		int minimumVal = omega(i + 1, Cup - a1[nextinline[i]], Clow);
		int maxVal = 0;

		int minval = (Prod_Up[nextinline[i]] < Clow) ? Prod_Up[nextinline[i]] : Clow;
		for (int w = Prod_Low[nextinline[i]]; w <= minval; w++) {
			int currentVal = omega (i + 1, Cup, Clow - w) + c[nextinline[i]] + (w * variable[nextinline[i]]);
			if (currentVal > maxVal) {
				maxVal = currentVal;
			}
		}

		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		return omg[i][Cup][Clow] = result;
	}
}




// The Branch and Bound Function

void bnb(int item, int xu[], int CupLeft, std::vector<std::vector<int>> Tempmemo // all the input required meaning a number to signal what item we are considering + the set of interdicted items x
) {
	// Make sure that the initilization happens each time as the final table depends a lot on the interdicted items, this structure ensures the calculations are done at most once. 
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= b2; j++) {
			memoP1[i][j] = Tempmemo [i][j];
		}
	}
	// Also always check if we have reached the final item as in that case we can compute the actual optimal solution
	if (item > n) {
		int tempmax = solveP1(n, b2, xu);
		reconstructSolutionP(n, b2);
		if (tempmax < gub) { // because we have a minimzation problem
			gub = tempmax;
			for (int i = 1; i <= n; i++) {
				x[i] = xu[i]; // Make sure to save the current best solution
				y[i] = yu[i]; // After the function solveP1(), yu and zu should be filled with the current best solution of the problem
				z[i] = zu[i];
			}
		}
		return;
	}
	// Here we should test the lower bound of the problem and compare it to the upper bound, can phantom the node if solution cannot be improved, this is done in a similar way as
	// as in the original paper. 
	int locallower = 0; 
	for (int w = b2; w >= 0; w--) {
		locallower = solveP1(item - 1, w, xu) + omg[item][CupLeft][b2 - w];
		if (locallower >= gub) {
			return;
		}
	} 
	// should be depth first branch and bound always checking if the upper level problem is still feasible
	if (a1[nextinline[item]] <= CupLeft) {
		xu[nextinline[item]] = 1;
		bnb(item + 1, xu, CupLeft - a1[nextinline[item]], memoP1);
		// Makes sure that the solution above does not influence the solution below
		for (int i = 1; i <= n; i++) {
			for (int j = 0; j <= b2; j++) {
				memoP1[i][j] = Tempmemo[i][j];
			}
		}
	}

	xu[nextinline[item]] = 0;
	bnb(item + 1, xu, CupLeft, memoP1);
	std::vector<vector<int>>().swap(Tempmemo);
	return;
}

void procedure(void) {
	// Need some sort of initialization to make sure we start with a completely empty solution
	memoBellman.assign(n + 1, std::vector<double>(b1 + 1, -1));
	memoP1.assign(n + 1, std::vector<int>(b2 + 1, -1)); // Allocate the correct memory to this global variable and initialize
	omg.assign(n + 1, vector<vector<int>>(b1 + 1, vector<int>(b2 + 1, -1)));
	gub = INF;
	for (int i = 1; i <= n; i++) {
		x[i] = 0;
		y[i] = 0;
		z[i] = 0;
		xv[i] = 0;
		yu[i] = 0;
		zu[i] = 0;
		ratio[i] = ((double)c[i] + (double)variable[i] * (double)Prod_Up[i]) / (double)Prod_Up[i]; //Calculate the ratio 
		cratio[i] = ratio[i];
		nextinline[i] = i; //will store the order based on the ratio! This helps keep track of where each item is located. 
	}
	
	// Some sort of sorting function first that sorts them based on the sorting algoritms above
	// Probably need to use the same tricks to keep track of all the actual values after the sorting has happened --> Maybe we do not need to implement the sorting in the beginning!
	quickSort(cratio, nextinline, 1, n);
	// Insert the tricks here

	// Start of the actual bnb procedure
	gub = greedyupper(); // Note that this also gives a first initial solution to the problem

	// This should actually run for all possible values of b1 and b2 and n
	clock_t start = clock();
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= b1; j++) {
			for (int k = 0; k <= b2; k++) {
				omega(i, j, k);
			}
		}
	}
	clock_t end = clock();
	double loop_time = (double)(end - start) / CLOCKS_PER_SEC;

	printf("Omega loop time: %.6f seconds\n", loop_time);
	// Make sure we start with a clear slate for the memory
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= b2; j++) {
			memoP1[i][j] = -1;
		}
	}
	 // This function computes all omega values before the algorithm is able to start
	 bnb(1, xv, b1, memoP1); // Make sure that the initial x is only 0 values in order for this to work
	 std::vector<vector<double>>().swap(memoBellman);
	 std::vector<vector<int>>().swap(memoP1);
	 std::vector<vector<vector<int>>>().swap(omg);
	return;
}










