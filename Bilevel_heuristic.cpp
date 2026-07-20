#include "Bilevel_orig.cpp"
#include <vector>
#include <unordered_map>


// Keeps track of which Omega we already know
struct OmegaState {
	int i, Cup, Clow, fractional;

	bool operator==(const OmegaState& other) const {
		return i == other.i &&
			Cup == other.Cup &&
			Clow == other.Clow &&
			fractional == other.fractional;
	}
};

// Sparse representation in case of large capacities in both levels
struct OmegaHash {
	size_t operator()(const OmegaState& s) const {
		// Fast integer hash (low collision, cheap)
		size_t h = s.i;
		h = h * 1315423911u + s.Cup;
		h = h * 2654435761u + s.Clow;
		h = h * 97531u + s.fractional;
		return h;
	}
};

// Sparse memoization
std::unordered_map<OmegaState, int, OmegaHash> omg_sparse;
using namespace std;
int  llb, lub; //These variables represent the global and local bounds as well as where the first fractional value is 
mat nextinline,sorteda,sortedc, initx; //These arrays are used for the initialization of the BnB algorithm, it keeps in mind the original item number
double ratio[maxn + 1]; //This arrays contains the ratio of value/weight
mat xv, yu, zu; //Being able to keep the current best solution saved somewhere



// Keeping track of intermediate solutions
vector<vector<double>> memoBellman;
std::vector<std::vector<std::vector<int>>> memoP1;
std::vector<std::vector<std::vector<int>>> choiceZ;
std::vector<int> omg_flat;
struct MemoEntry { int i, j, k, oldVal; };
std::vector<MemoEntry> memoLog;

//Sorting the items based on their variable profit structure 
int partition(int arr[], int arr2[], int low, int high)
{
	int pivotItem = arr2[high];
	int i = low - 1;

	for (int j = low; j <= high - 1; j++) {
		int itemJ = arr2[j];

		if (variable[itemJ] > variable[pivotItem]) {
			i++;
			swap(arr[i], arr[j]);
			swap(arr2[i], arr2[j]);
		}
	}

	swap(arr[i + 1], arr[high]);
	swap(arr2[i + 1], arr2[high]);
	return i + 1;
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
//The Bellman recursion to solve the knapsack problem from the greedy algorithm
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



inline void writeMemo(int i, int j, int k, int val) {
	if (memoP1[i][j][k] == -1) { // only log the first write (was -1 = unset)
		memoLog.push_back({ i, j, k, -1 });
	}
	memoP1[i][j][k] = val;
}

// This function solves the lower-level problem using the DP presented in the paper
int solveP1(int i, int Cleft, int xu[], int fractional) {
	if (time_limit_reached) return 0;

	if ((double)(clock() - global_start) / CLOCKS_PER_SEC >= TIME_LIMIT) {
		time_limit_reached = true;
		return 0;
	}
	int bestW = 0;
	if (Cleft < 0) return -INF;
	if (i == 0) return 0;
	if (memoP1[i][Cleft][fractional] != -1) return memoP1[i][Cleft][fractional];

	if (Prod_Low
		[nextinline[i]] > Cleft || xu[nextinline[i]] == 1)
	{
		choiceZ[i][Cleft][fractional] = 0; // 0 means not taken/impossible
		writeMemo(i, Cleft, fractional, solveP1(i - 1, Cleft, xu, fractional));  return solveP1(i - 1, Cleft, xu, fractional);
	}
	else if (Prod_Low
		[nextinline[i]] == Prod_Up
		[nextinline[i]])
	{
	int maxVal = solveP1(i - 1, Cleft, xu, fractional);
	if (maxVal < solveP1(i - 1, Cleft - Prod_Low[nextinline[i]], xu, fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]])) {
			bestW = Prod_Low[nextinline[i]];
			maxVal = solveP1(i - 1, Cleft - Prod_Low[nextinline[i]], xu, fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]]);
		}
	choiceZ[i][Cleft][fractional] = bestW;
	writeMemo(i, Cleft, fractional, maxVal);  return maxVal;
	}
	else if (fractional == 0) {
		int maxVal = solveP1(i - 1, Cleft, xu, fractional);
		if (maxVal < solveP1(i - 1, Cleft - Prod_Low[nextinline[i]], xu, fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]])) {
			bestW = Prod_Low[nextinline[i]];
			maxVal = solveP1(i - 1, Cleft - Prod_Low[nextinline[i]], xu, fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]]);
		}
		for (int w = Prod_Low[nextinline[i]] + 1; w <= Prod_Up[nextinline[i]]; w++) {
			int currentVal = solveP1(i - 1, Cleft - w, xu, 1) + c[nextinline[i]] + (w * variable[nextinline[i]]);
			if (currentVal > maxVal) {
				maxVal = currentVal;
				bestW = w; // Store the specific weight that worked best

			}
		}
		choiceZ[i][Cleft][fractional] = bestW;
		writeMemo(i, Cleft, fractional, maxVal);  return maxVal;
	}
	else {
		int maxVal = solveP1(i - 1, Cleft, xu, fractional);
		int currentVal = solveP1(i - 1, Cleft - Prod_Up[nextinline[i]], xu, fractional) + c[nextinline[i]] + (Prod_Up[nextinline[i]] * variable[nextinline[i]]);
		if (currentVal > maxVal) {
			maxVal = currentVal;
			bestW = Prod_Up[nextinline[i]]; // Store the specific weight that worked best
		}
		choiceZ[i][Cleft][fractional] = bestW;
		writeMemo(i, Cleft, fractional, maxVal);  return maxVal;
	}
}


// Reconstruct the best lower-level solution found using the above DP procedure
void reconstructSolutionP(int n, int Cleft) {
	int currentCap = Cleft;
	int bestfrac = 0;
	for (int i = n; i > 0; i--) {
		int optimalW = choiceZ[i][currentCap][bestfrac];
		
		if (optimalW > 0) {
			yu[nextinline[i]] = 1;
			zu[nextinline[i]] = optimalW;
			currentCap -= optimalW;
			if (optimalW != Prod_Low[nextinline[i]]) {
				bestfrac = 1;
			}
		}
		else {
			yu[nextinline[i]] = 0;
			zu[nextinline[i]] = 0;
		}
	}
}


// Initial greedy algorithm to compute the first upper bound of the problem

int greedyupper(void) {
	// Compute the optimal solution of the upper level knapsack problem (using the above Bellman algorithm) using the previously calculated ratio as objective function value
	

	
	Bellman(n, b1);
	reconstructSolutionB(n, b1);
	// Somehow save this initial solution somewhere as it might already be the optimal solution, this should happen in the procedure function somewhere
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= b2; j++) {
			for (int k = 0; k <= 1; k++) {
				memoP1[i][j][k] = -1;
			}
		}
	}
	// Using this solution compute the optimal solution of the lower level problem using the DP presented above
	int bestsolution = solveP1(n, b2, x, 0);
	reconstructSolutionP(n, b2);
	for (int i = 1; i <= n; i++) {
		// Make sure to save the current best solution
		y[i] = yu[i];
		z[i] = zu[i];
	}
	return bestsolution; // the upper bound value
}
// Omega which is used to compute a strong lower bound for the problem. 
int (*omega)(int, int, int, int);
// When omega is expected to be sparse
int omega_sparse(int i, int Cup, int Clow, int fractional) {
	if (time_limit_reached) return 0;

	if ((double)(clock() - global_start) / CLOCKS_PER_SEC >= TIME_LIMIT) {
		time_limit_reached = true;
		return 0;
	}
	// 1. Boundary Checks
	if (Cup < 0) return INF; // Using a safe high value for "infinity"
	if (Clow < 0) return -INF;
	if (i > n) return 0;

	// 2. Check if we've already solved this specific state
	OmegaState state{ i, Cup, Clow, fractional };
	auto it = omg_sparse.find(state);
	if (it != omg_sparse.end()) {
		return it->second;
	}

	int result; // Local variable to hold the calculation


	if (Prod_Low[nextinline[i]] > Clow) {
		result = omega_sparse(i + 1, Cup, Clow, fractional);
		omg_sparse[state] = result;
		return result;
	}
	else if (Prod_Low[nextinline[i]] == Prod_Up[nextinline[i]])
	{
		int minimumVal = omega_sparse(i + 1, Cup - a1[nextinline[i]], Clow, fractional);
		int maxVal = 0;
		int val1 = omega_sparse(i + 1, Cup, Clow - Prod_Low[nextinline[i]], fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]]);
		int val2 = omega_sparse(i + 1, Cup, Clow, fractional);
		maxVal = (val1 > val2) ? val1 : val2;
		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		omg_sparse[state] = result;
		return result;
	}
	else if (fractional == 0) {
		int minimumVal = omega_sparse(i + 1, Cup - a1[nextinline[i]], Clow, fractional);
		int maxVal = 0;
		int currentVal = 0;

		int minval = (Prod_Up[nextinline[i]] < Clow) ? Prod_Up[nextinline[i]] : Clow;
		for (int w = Prod_Low[nextinline[i]]; w <= minval; w++) {
			if (w == Prod_Up[nextinline[i]]) {
				currentVal = omega_sparse(i + 1, Cup, Clow - w, 0) + c[nextinline[i]] + (w * variable[nextinline[i]]);
	
			}
			else {
				currentVal = omega_sparse(i + 1, Cup, Clow - w, 1) + c[nextinline[i]] + (w * variable[nextinline[i]]);
			}
			if (currentVal > maxVal) {
				maxVal = currentVal;
			}
		}
		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		omg_sparse[state] = result;
		return result;

	}
	else {
		int minimumVal = omega_sparse(i + 1, Cup - a1[nextinline[i]], Clow, fractional);
		int maxVal = 0;
		int val1 = omega_sparse(i + 1, Cup, Clow - Prod_Low[nextinline[i]], fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]]);
		int val2 = omega_sparse(i + 1, Cup, Clow, fractional);
		maxVal = (val1 > val2) ? val1 : val2;
		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		omg_sparse[state] = result;
		return result;
	}
}
inline int& omg_at(int i, int Cup, int Clow, int f) {
	return omg_flat[
		(((size_t)i * (b1 + 1) + Cup) * (b2 + 1) + Clow) * 2 + f
	];
}
// Omega when no sparse matrix is expected
int omega_dense(int i, int Cup, int Clow, int fractional) {
	if (time_limit_reached) return 0;

	if ((double)(clock() - global_start) / CLOCKS_PER_SEC >= TIME_LIMIT) {
		time_limit_reached = true;
		return 0;
	}
	// 1. Boundary Checks 
	if (Cup < 0) return INF; // Using a safe high value for "infinity"
	if (Clow < 0) return -INF;
	if (i > n) return 0;

	int& memo = omg_at(i, Cup, Clow, fractional);
	// 2. Check if we've already solved this specific state
	if (memo != -1) {
		return memo;
	}

	int result; // Local variable to hold the calculation


	if (Prod_Low[nextinline[i]] > Clow) {
		result = omega_dense(i + 1, Cup, Clow, fractional);
		return memo = result;
	}
	else if (Prod_Low[nextinline[i]] == Prod_Up[nextinline[i]]) 
	{
		int minimumVal = omega_dense(i + 1, Cup - a1[nextinline[i]], Clow, fractional);
		int maxVal = 0;
		int val1 = omega_dense(i + 1, Cup, Clow - Prod_Low[nextinline[i]], fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]]);
		int val2 = omega_dense(i + 1, Cup, Clow, fractional);
		maxVal = (val1 > val2) ? val1 : val2;
		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		return memo = result;
	}
	else if (fractional == 0) {
		int minimumVal = omega_dense(i + 1, Cup - a1[nextinline[i]], Clow, fractional);
		int maxVal = 0;
		int currentVal = 0;

		int minval = (Prod_Up[nextinline[i]] < Clow) ? Prod_Up[nextinline[i]] : Clow;
		for (int w = Prod_Low[nextinline[i]]; w <= minval; w++) {
			if (w == Prod_Up[nextinline[i]]) {
				currentVal = omega_dense(i + 1, Cup, Clow - w, 0) + c[nextinline[i]] + (w * variable[nextinline[i]]);

			}
			else {
				currentVal = omega_dense(i + 1, Cup, Clow - w, 1) + c[nextinline[i]] + (w * variable[nextinline[i]]);
			}
			if (currentVal > maxVal) {
				maxVal = currentVal;
			}
		}
		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		return memo = result;

	}
	else {
		int minimumVal = omega_dense(i + 1, Cup - a1[nextinline[i]], Clow, fractional);
		int maxVal = 0;
		int val1 = omega_dense(i + 1, Cup, Clow - Prod_Low[nextinline[i]], fractional) + c[nextinline[i]] + (Prod_Low[nextinline[i]] * variable[nextinline[i]]);
		int val2 = omega_dense(i + 1, Cup, Clow, fractional);
		maxVal = (val1 > val2) ? val1 : val2;
		result = (minimumVal < maxVal) ? minimumVal : maxVal;
		return memo = result;
	}
}


// The main Branch-and-Bound Function

int bnb(int item, int xu[], int CupLeft, int inheritedLB ) {
	nodes++;
	if (time_limit_reached) return inheritedLB;
	double elapsed = (double)(clock() - global_start) / CLOCKS_PER_SEC;
	if (elapsed >= TIME_LIMIT) {
		time_limit_reached = true;
		return inheritedLB;
	}
	
	// Also always check if we have reached the final item as in that case we can compute the actual optimal solution
	if (item > n) {
		int bestsolution = solveP1(n, b2, xu, 0);
		reconstructSolutionP(n, b2);
		if (bestsolution < gub) { // because we have a minimzation problem
			gub = bestsolution;
			for (int i = 1; i <= n; i++) {
				x[i] = xu[i]; // Make sure to save the current best solution
				y[i] = yu[i]; // After the function solveP1(), yu and zu should be filled with the current best solution of the problem
				z[i] = zu[i];
			}
		}
		return bestsolution;
	}
	// Here we should test the lower bound of the problem and compare it to the upper bound, can phantom the node if solution cannot be improved
	int locallower = 0;
	int nodeLB = inheritedLB;
	for (int w = b2; w >= 0; w--) {
		if (time_limit_reached) return nodeLB;
		elapsed = (double)(clock() - global_start) / CLOCKS_PER_SEC;
		if (elapsed >= TIME_LIMIT) {
			time_limit_reached = true;
			return nodeLB;
		}
		locallower = solveP1(item - 1, w, xu, 0) + omega(item, CupLeft, b2 - w, 0);
		nodeLB = max(nodeLB, locallower);
		if (nodeLB >= gub) return nodeLB;
	}
	// Depth first branch-and-bound always checking if the upper level problem is still feasible
	int bestSubtreeLB = INF;
	if (a1[nextinline[item]] <= CupLeft) {
		xu[nextinline[item]] = 1;

		size_t savepoint = memoLog.size();

		int lbTake = bnb(item + 1, xu, CupLeft - a1[nextinline[item]], nodeLB);
		bestSubtreeLB = min(bestSubtreeLB, lbTake);
		
		while (memoLog.size() > savepoint) {
			auto& e = memoLog.back();
			memoP1[e.i][e.j][e.k] = e.oldVal;
			memoLog.pop_back();	
		}
		if (time_limit_reached) return min(bestSubtreeLB, nodeLB);
		// Makes sure that the solution above does not influence the solution below
	}

	xu[nextinline[item]] = 0;
	int lbSkip = bnb(item + 1, xu, CupLeft, nodeLB);
	bestSubtreeLB = min(bestSubtreeLB, lbSkip);
	return max(nodeLB, bestSubtreeLB);
}

void procedure(void) {
	global_start = clock();
	// Need some sort of initialization to make sure we start with a completely empty solution
	memoBellman.assign(n + 1, std::vector<double>(b1 + 1, -1));
	memoP1.assign(n + 1,
		std::vector<std::vector<int>>(b2 + 1,
			std::vector<int>(2, -1)));// Allocate the correct memory to this global variable and initialize
	choiceZ.assign(n + 1,
		std::vector<std::vector<int>>(b2 + 1,
			std::vector<int>(2, -1)));
	gub = INF;
	avg_a1 = 0.0;
	avg_width = 0.0;
	avg_Prod_Low = 0.0;
	for (int i = 1; i <= n; i++) {
		double convexcomb = ((double)variable[i] * (double)Prod_Up[i]) / ((double)c[i] + (double)variable[i] * (double)Prod_Up[i]);
		convexcomb = (1 - convexcomb)*(double)Prod_Low[i] + convexcomb * (double)Prod_Up[i];
		x[i] = 0;
		y[i] = 0;
		z[i] = 0;
		xv[i] = 0;
		yu[i] = 0;
		zu[i] = 0;
		ratio[i] = ((double)c[i] + (double)variable[i] * (double)convexcomb) / (double)convexcomb;
		cvariable[i] = variable[i];
		nextinline[i] = i; //will store the order based on the variable profit! This helps keep track of where each item is located. 
		avg_a1 = avg_a1 + a1[i];
		avg_width += (double) (Prod_Up[i] - Prod_Low[i] + 1);
		avg_Prod_Low += Prod_Low[i];
	}
	avg_a1 = (double)avg_a1/n;
	avg_a1 = min((double)n, (double)b1 / avg_a1);
	avg_a1 = (double)avg_a1 / n;
	avg_width = (double)avg_width/n;
	avg_Prod_Low = (double)avg_Prod_Low / n;
	avg_Prod_Low = min((double)n, (double)b2 / avg_Prod_Low);
	avg_Prod_Low = (double)avg_Prod_Low / n;

	quickSort(cvariable, nextinline, 1, n);


	// Start of the actual bnb procedure
	gub = greedyupper(); // Note that this also gives a first initial solution to the problem


	// Make sure we start with a clear slate for the memory
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= b2; j++) {
			for (int k = 0; k <= 1; k++) {
				memoP1[i][j][k] = -1;
			}
		}
	}
	size_t states = (size_t)(n + 1) * (b1 + 1) * (b2 + 1) * 2;
	size_t bytes = states * sizeof(int);
	const size_t LIMIT = (size_t)10 * 1024 * 1024 * 1024; // 8 GB
	printf("\n	bytes: %.2f \n", (double)bytes / (1024.0 * 1024.0 * 1024.0));
	printf("\n	Limit: %.2f \n", (double)LIMIT / (1024.0 * 1024.0 * 1024.0));
	double threshold = (double)avg_a1 * (double)avg_Prod_Low;
	
	printf("\n	Treshold: %.10f \n", threshold);
	if (bytes > LIMIT || threshold < 0.0025) {
		omg_sparse.reserve(20'000'000);
		omega = omega_sparse;
	}
	else {
		omg_flat.assign(
			(size_t)(n + 1) * (b1 + 1) * (b2 + 1) * 2,
			-1
		);
		omega = omega_dense;
	}
	
	glb = gub; // Make sure that the initial x is only 0 values in order for this to work
	 std::vector<vector<double>>().swap(memoBellman); // Avoiding memory leaks
	 std::vector<std::vector<std::vector<int>>>().swap(memoP1);
	 std::vector<std::vector<std::vector<int>>>().swap(choiceZ);
	 if (bytes > LIMIT || threshold < 0.0025) {
		 omg_sparse.clear();
		 omg_sparse.rehash(0);
	 }
	 else
	 {
		 std::vector<int>().swap(omg_flat);
	 }
	
	return;
}










