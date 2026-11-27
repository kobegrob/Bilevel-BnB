#include "Bilevel_orig_Log.cpp"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <deque>
#include <list>
#include <memory>
#include <algorithm>
#include <climits>


using namespace std;
// =====================================================================================================
// Bilevel_StrongerLBNLLog.cpp: nonlinear profit c + w*ln(w)*v (logarithmic), with the performance fixes of the linear solver
// (generated on 2026-11-24 by build/analysis/gen_nonlinear.py from the linear source with Fixes 1-8 and
// the original Bilevel_StrongerLBNLLog.cpp). Item order, greedy start and all recurrences are those of the original
// nonlinear code, so the optimal values are the same. See the block "NONLINEAR PORT" below for the one
// fix that had to be done differently (Fix 3's sliding window needs a linear profit).
// =====================================================================================================
double llb, lub; //These variables represent the global and local bounds as well as where the first fractional value is 
mat nextinline,sorteda,sortedc, initx; //These arrays are used for the initialization of the BnB algorithm, it keeps in mind the original item number
double ratio[maxn + 1], ratio_value[maxn + 1]; //This arrays contains the ratio of value/weight
mat xv, yu, zu; //Being able to keep the current best solution saved somewhere
mat UpperW,LowerW, LowerW2, LowerUpper; // Used in the lower bound computations



// Keeping track of intermediate solutions
vector<vector<double>> memoBellman;
// =====================================================================================================
// PERFORMANCE FIX 3 (2026-11-23): sliding-window maximum in the two DPs (solveP1 and omega)
// -----------------------------------------------------------------------------------------------------
// Problem: in every DP state where the follower may pack item i with any quantity w in [l, u], the old
//   code looped over all w ("for (int w = wMin; w <= minval; w++)" in solveP1 and omega). That costs
//   O(u - l) per state, so each DP level cost O(C^L * (u - l)): quadratic when capacities and item
//   quantities grow together (measured: doubling all capacities -> about 4x the run time).
// Idea: packing w units of item i is worth c + w * v, which is linear in w. Writing y = C - w (the
//   capacity left for the other items), the loop computes
//        max over y in [lo(C), hi(C)] of ( P(y) - y * v )   plus   C * v + c,
//   and both window ends lo(C), hi(C) move right by one when C grows by one. A monotone deque gives the
//   maximum of such a sliding window in O(1) amortised time per C (the classic bounded-knapsack trick),
//   so a whole row of states over C costs O(C^L) instead of O(C^L * (u - l)).
// Consequence for the code structure: the window needs the whole row of the next/previous level, so both
//   DPs are now computed row by row (all capacities of one level at once) instead of state by state:
//   * solveP1: rows P(i, C, f), i = 1..n, live in flat arrays. Row i depends only on the interdiction
//     decisions xu[] of the items at levels 1..i, so when bnb() changes xu for the item at level `item`,
//     rows >= item are marked invalid (p1_invalidate_from) and recomputed when they are next needed. This
//     replaces the old memoLog undo log. The flat arrays also use ~7x less memory than the old
//     vector<vector<vector<int>>> memoP1/choiceZ.
//   * omega: a row Omega(i, Cup, ., f) over all lower-level capacities is stored for each (i, Cup) pair
//     that is actually reached (hash map from (i, Cup) to the row). This replaces both the dense table,
//     which reserved (n+1)*(C^U+1)*(C^L+1)*2 ints up front, and the per-state hash map, so the
//     dense/sparse choice of PERFORMANCE FIX 2 is no longer needed in this file.
// The values computed are exactly the same as before (same recurrences and boundary rules; on ties the
//   smallest quantity w is kept, as in the old ascending loop with a strict '>'), so bounds, the search
//   tree and the optimal value are unchanged. Only the linear objective allows this trick: in this
//   nonlinear file the loops use the step-start enumeration instead (see NONLINEAR PORT below).
// =====================================================================================================


// Lower-level DP rows (PERFORMANCE FIX 3; replace the old memoP1, choiceZ and memoLog)
std::vector<double> p1Val;          // p1Val[p1_idx(i, C, f)] = P(i, C, f) for 1 <= i <= n, 0 <= C <= UpperW[i]
std::vector<int> p1Choice;       // quantity of item nextinline[i] in that optimum (0 = not packed), for reconstructSolutionP
int p1RowsValid = 0;             // rows 1..p1RowsValid are valid for the interdiction array p1RowsFor
const int* p1RowsFor = nullptr;  // the xu[] array the rows were computed for
inline size_t p1_idx(int i, int C, int f) { return ((size_t)i * (size_t)(b2 + 1) + (size_t)C) * 2 + f; }
inline void p1_invalidate_from(int level) { if (p1RowsValid > level - 1) p1RowsValid = std::max(0, level - 1); }
inline void p1_invalidate_all() { p1RowsValid = 0; p1RowsFor = nullptr; }
// P(i, C, f) with the boundary rules of the old recursion: -INF for a negative capacity, 0 without items,
// capacities above UpperW[i] behave like UpperW[i]. Only valid for rows that are already computed.
inline double P1(int i, int C, int f) {
	if (C < 0) return -INF;
	if (i == 0) return 0;
	if (C > UpperW[i]) C = UpperW[i];
	return p1Val[p1_idx(i, C, f)];
}

// Profit model (logarithmic), without integer truncation.
inline double profit_at(int j, int w) { return c[j] + (w == 0 ? 0.0 : std::log((double)w) * w * variable[j]); }
inline double add_profit(double base, int j, int w) { return base + profit_at(j, w); }
inline double add_profit2(double base, int j1, int w1, int j2, int w2) { return base + profit_at(j1, w1) + profit_at(j2, w2); }
inline double add_profit_case4_u2(double base, int j) { return add_profit(base, j, Prod_Up[j]); }

// =====================================================================================================
// NONLINEAR PORT (2026-11-24): step-start enumeration instead of the sliding window of PERFORMANCE FIX 3
// -----------------------------------------------------------------------------------------------------
// The loops "for (int w = wMin; w <= minval; w++)" of the original solveP1 and omega pack an item with every
// quantity w strictly between its bounds l < w < u on top of a row read with fractional = 1. With a nonlinear
// profit the sliding-window trick of Fix 3 does not apply, but these loops can still skip almost all w:
//  * A row read with fractional = 1 only allows the other items to take 0, l or u units. As a function of the
//    capacity y = C - w, its value therefore changes only at sums of such quantities (the "steps" of the row):
//    it is constant from one sum to the next. For the lower level the sums are over the items at levels
//    1..i-1 (stepsPrefix), for the bound over the items at levels i+1..n or i+2..n (stepsSuffix).
//  * The profit of the item does not decrease with w (both w*w*v and w*ln(w)*v are non-decreasing for w >= 1),
//    so on a range of y where the row is constant the best w is the largest one, i.e. the smallest y.
//  Hence only y equal to a sum inside the window, plus the lowest y of the window, can be optimal; the other
//  quantities are never better. The maximum is the same as with the full loop (exact for any non-decreasing
//  profit); when several w give the same value a different one may be recorded, which only changes which of the
//  equally good follower solutions is reconstructed. The bounds l and u are evaluated as before.
// The step lists contain every sum of 0, l or u units of the items, capped at C^L: they are computed once per
// instance and are a superset of the steps of every row (interdicted items only remove sums).
// =====================================================================================================
std::vector<std::vector<int>> stepsPrefix, stepsSuffix;
void build_step_lists() {
	stepsPrefix.assign(n + 2, std::vector<int>()); stepsSuffix.assign(n + 3, std::vector<int>());
	std::vector<char> r((size_t)b2 + 1, 0), nr;
	auto collect = [&](std::vector<int>& out) { for (int y = 0; y <= b2; y++) if (r[y]) out.push_back(y); };
	auto add_item = [&](int item) {
		nr = r; int l = Prod_Low[item], u = Prod_Up[item];
		for (int y = 0; y <= b2; y++) if (r[y]) { if (y + l <= b2) nr[y + l] = 1; if (y + u <= b2) nr[y + u] = 1; }
		r.swap(nr);
	};
	r[0] = 1; collect(stepsPrefix[0]);
	for (int i = 1; i <= n; i++) { add_item(nextinline[i]); collect(stepsPrefix[i]); }
	std::fill(r.begin(), r.end(), 0); r[0] = 1; collect(stepsSuffix[n + 1]); stepsSuffix[n + 2] = stepsSuffix[n + 1];
	for (int i = n; i >= 1; i--) { add_item(nextinline[i]); collect(stepsSuffix[i]); }
}
// Calls fn(w) for the quantities w in [wlo, whi] that can be optimal when the row read at C - w is constant
// between consecutive elements of T: w = C - y for every y of T in [C - whi, C - wlo] (largest y first), then
// w = whi (the lowest y of the window).
template <class Fn> inline void for_step_starts(int C, int wlo, int whi, const std::vector<int>& T, Fn fn) {
	if (wlo > whi) return;
	int lo = C - whi, hi = C - wlo;
	auto itHi = std::upper_bound(T.begin(), T.end(), hi), itLo = std::upper_bound(T.begin(), T.end(), lo);
	for (auto it = itHi; it != itLo; ) { --it; fn(C - *it); }
	fn(whi);
}

// =====================================================================================================
// PERFORMANCE FIXES 4-8 (2026-11-24). Fixes 4, 5, 7 and 8 do not change any value the algorithm computes;
// fix 6 only gives the search a better first incumbent, so the optimal value is unchanged but fewer nodes
// may be visited. Profiling of fix 3 on a hard instance (n = 100, capacities ~2,500) showed that 64% of the
// run time was spent in clock() system calls made by the lower-bound loop of bnb() (fix 4), and that the
// omega rows used 2 GB (fixes 5 and 7).
//  FIX 4: bnb() computed its node lower bound with a loop over the follower capacity w that read the clock
//         and called solveP1() and omega() for every w. The omega row Omega(item, CupLeft, ., .) and the
//         lower-level row P(item - 1, ., 0) are the same for all w, so they are fetched once per node and the
//         loop only adds two array entries. The time limit is checked once per node and when rows are built.
//  FIX 5: Omega(i, K, ., .) depends on the leader's remaining capacity K only through which subsets of the
//         weights of items i..n fit in K. Every K is therefore replaced by the largest subset sum of those
//         weights that does not exceed it (computed once by a subset-sum DP), so capacities that allow the
//         same interdictions share one row. This generalises the old cap of K at the total remaining weight.
//  FIX 6: the incumbent is improved by a short local search (add an item to the interdiction set, or swap one
//         interdicted item for another), each candidate evaluated exactly with solveP1(). Bounded by 2n + 20
//         evaluations and 1% of the time limit. Revised the same day: it runs once inside bnb(), after
//         LOCAL_SEARCH_AFTER_NODES nodes, so instances solved with fewer nodes do not pay for it (before the
//         revision it cost 25-30 ms on instances solved in 30-60 ms, with no node saved).
//  FIX 7: omega rows are kept in a least-recently-used cache whose size is limited to a share of the memory
//         available to the process. Rows that are evicted are recomputed if they are needed again (values
//         are deterministic), so a run that would have run out of memory runs slower instead.
//  FIX 8: interdicting an item that the follower can never pack (l_i > C^L) uses leader capacity and cannot
//         lower the follower's value, so bnb() no longer branches on it.
// =====================================================================================================

// Omega rows (PERFORMANCE FIX 3), kept in a cache with a memory budget (PERFORMANCE FIX 7)
typedef std::shared_ptr<std::vector<double>> RowPtr; // a row Omega(i, K, Clow, f) at [Clow * 2 + f]; holding it keeps it alive
struct OmCacheEntry { RowPtr row; std::list<long long>::iterator lru; size_t bytes; };
std::unordered_map<long long, OmCacheEntry> omCache; // key i * (b1 + 1) + K
std::list<long long> omLru;                         // most recently used first
size_t omBytes = 0, omBudget = (size_t)-1, omPeakBytes = 0, omComputed = 0, omEvicted = 0;
const double OMEGA_CACHE_MEMORY_FRACTION = 0.6;     // share of the process's memory that the cache may use
const size_t OMEGA_ROW_OVERHEAD = 160;              // bytes of bookkeeping per cached row (map, list, control block)
RowPtr om_cache_get(long long key) {
	auto it = omCache.find(key);
	if (it == omCache.end()) return nullptr;
	omLru.splice(omLru.begin(), omLru, it->second.lru); // now the most recently used
	return it->second.row;
}
RowPtr om_cache_put(long long key, std::vector<double>&& cells) {
	RowPtr row = std::make_shared<std::vector<double>>(std::move(cells));
	size_t bytes = row->size() * sizeof(double) + OMEGA_ROW_OVERHEAD;
	omLru.push_front(key);
	omCache[key] = OmCacheEntry{ row, omLru.begin(), bytes };
	omBytes += bytes; omComputed++;
	while (omBytes > omBudget && omLru.size() > 1) { // evict the least recently used rows; rows in use stay alive through their RowPtr
		auto victim = omCache.find(omLru.back());
		omBytes -= victim->second.bytes; omCache.erase(victim); omLru.pop_back(); omEvicted++;
	}
	if (omBytes > omPeakBytes) omPeakBytes = omBytes;
	return row;
}

// Canonical leader capacities (PERFORMANCE FIX 5): canonCup[i * (b1 + 1) + K] = largest subset sum of the weights of
// the items at levels i..n that is <= K.
std::vector<int> canonCup;
bool canonOn = false;
const size_t CANON_MAX_BYTES = (size_t)256 * 1024 * 1024; // above this size the table is not built (plain capacities are used)
inline int canonical_cup(int i, int K) { return canonOn ? canonCup[(size_t)i * (size_t)(b1 + 1) + K] : K; }
void build_canonical_capacities() {
	canonOn = false; std::vector<int>().swap(canonCup);
	size_t cells = (size_t)(n + 2) * (size_t)(b1 + 1);
	if (cells * sizeof(int) > CANON_MAX_BYTES) return;
	canonCup.assign(cells, 0);
	std::vector<char> reach((size_t)b1 + 1, 0); reach[0] = 1; // subset sums of the weights of the items at levels > i
	for (int i = n; i >= 1; i--) {
		int w = a1[nextinline[i]];
		for (int K = b1; K >= w; K--) if (reach[K - w]) reach[K] = 1;
		int last = 0; int* row = &canonCup[(size_t)i * (size_t)(b1 + 1)];
		for (int K = 0; K <= b1; K++) { if (reach[K]) last = K; row[K] = last; }
	}
	canonOn = true;
}
// Memory budget of the omega cache: a share of the process's memory, minus the other large tables.
void setup_omega_cache() {
	double mem = available_memory_bytes();
	double others = (double)p1Val.size() * sizeof(double) + (double)(p1Choice.size() + canonCup.size()) * sizeof(int) + (double)(n + 1) * (b1 + 1) * sizeof(double);
	if (mem <= 0.0) omBudget = (size_t)-1;
	else omBudget = (size_t)std::max(64.0 * 1024 * 1024, mem * OMEGA_CACHE_MEMORY_FRACTION - others);
	omBytes = omPeakBytes = omComputed = omEvicted = 0;
}

//Sorting the items based on their variable profit structure 
int partition(double arr[], int arr2[], int low, int high)
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

void quickSort(double arr[], int arr2[], int low, int high)
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
		double value2 = Bellman(i - 1, CupLeft - a1[nextinline[i]]) + ratio[nextinline[i]] * ratio_value[nextinline[i]];

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



// This function solves the lower-level problem using the DP presented in the paper.
// PERFORMANCE FIX 3 (row-based): computes one full row P(i, ., .) from row i - 1, with the recurrence of the original
// nonlinear solveP1; the loop over the quantity w uses the step-start enumeration (NONLINEAR PORT).
void compute_p1_row(int i, const int xu[]) {
	int item = nextinline[i];
	int l = Prod_Low[item], u = Prod_Up[item];
	bool interdicted = (xu[item] == 1);
	const std::vector<int>& T = stepsPrefix[i - 1]; // steps of the rows P(i - 1, ., 1)
	for (int C = 0; C <= UpperW[i]; C++) {
		for (int f = 0; f <= 1; f++) {
			double best; int choice = 0;
			if (l > C || interdicted) { // item cannot be packed
				best = P1(i - 1, C, f);
			}
			else if (l == u) { // fixed quantity: skip or pack l units
				best = P1(i - 1, C, f);
				double t = add_profit(P1(i - 1, C - l, f), item, l);
				if (best < t) { best = t; choice = l; }
			}
			else if (f == 0) { // skip, pack l or u (fractional stays 0), or any w strictly between (then fractional = 1)
				best = P1(i - 1, C, 0);
				double t = add_profit(P1(i - 1, C - l, 0), item, l);
				if (best < t) { best = t; choice = l; }
				t = add_profit(P1(i - 1, C - u, 0), item, u);
				if (best < t) { best = t; choice = u; }
				// old loop: w from max(l + 1, C - UpperW[i - 1]) to min(u - 1, C) reading P(i - 1, C - w, 1)
				for_step_starts(C, std::max(l + 1, C - UpperW[i - 1]), std::min(u - 1, C), T, [&](int w) {
					double v = add_profit(P1(i - 1, C - w, 1), item, w);
					if (v > best) { best = v; choice = w; }
				});
			}
			else { // fractional already used: skip, or pack u or l
				best = P1(i - 1, C, 1);
				double t = add_profit(P1(i - 1, C - u, 1), item, u);
				if (t > best) { best = t; choice = u; }
				t = add_profit(P1(i - 1, C - l, 1), item, l);
				if (t > best) { best = t; choice = l; }
			}
			p1Val[p1_idx(i, C, f)] = best;
			p1Choice[p1_idx(i, C, f)] = choice;
		}
	}
}

// Same interface as before; computes the missing rows (bottom-up) and returns P(i, Cleft, fractional).
double solveP1(int i, int Cleft, int xu[], int fractional) {
	if (time_limit_reached) return 0;

	if (dp_time_limit_exceeded()) { // PERFORMANCE FIX 1: was clock() on every call; now checked every TIME_CHECK_INTERVAL calls (see Bilevel_orig.cpp)
		time_limit_reached = true;
		return 0;
	}
	if (Cleft < 0) return -INF;
	if (i == 0) return 0;
	if (p1RowsFor != xu) { p1RowsFor = xu; p1RowsValid = 0; } // rows were computed for another interdiction array
	while (p1RowsValid < i) {
		if ((double)(clock() - global_start) / CLOCKS_PER_SEC >= TIME_LIMIT) { time_limit_reached = true; return 0; }
		compute_p1_row(p1RowsValid + 1, xu);
		p1RowsValid++;
	}
	return P1(i, Cleft, fractional);
}


// Reconstruct the best lower-level solution found using the above DP procedure
void reconstructSolutionP(int n, int Cleft) {
	int currentCap = Cleft;
	int bestfrac = 0;
	for (int i = n; i > 0; i--) {
		if (currentCap > UpperW[i]) currentCap = UpperW[i];
		int optimalW = p1Choice[p1_idx(i, currentCap, bestfrac)]; // PERFORMANCE FIX 3: was choiceZ[i][currentCap][bestfrac]

		if (optimalW > 0) {
			yu[nextinline[i]] = 1;
			zu[nextinline[i]] = optimalW;
			currentCap -= optimalW;
			if (optimalW != Prod_Low[nextinline[i]] && optimalW != Prod_Up[nextinline[i]]) { // nonlinear code: both bounds keep fractional = 0
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

double greedyupper(void) {
	// Compute the optimal solution of the upper level knapsack problem (using the above Bellman algorithm) using the previously calculated ratio as objective function value
	

	
	Bellman(n, b1);
	reconstructSolutionB(n, b1);
	// Somehow save this initial solution somewhere as it might already be the optimal solution, this should happen in the procedure function somewhere
	p1_invalidate_all(); // PERFORMANCE FIX 3: rows are recomputed for the new interdiction array
	// Using this solution compute the optimal solution of the lower level problem using the DP presented above
	double bestsolution = solveP1(n, b2, x, 0);
	reconstructSolutionP(n, b2);
	for (int i = 1; i <= n; i++) {
		// Make sure to save the current best solution
		y[i] = yu[i];
		z[i] = zu[i];
	}
	return bestsolution; // the upper bound value
}
// PERFORMANCE FIX 6: improve the greedy starting solution x (with value cur) by local search. Moves: add one item
// to the interdiction set, or swap one interdicted item for one that is not; a move is kept only if the follower's
// optimal value (solveP1, exact) strictly decreases. The best solution found is left in x, y, z.
mat xls; // candidate interdiction vector of the local search
#ifdef VALIDATION_NO_LOCAL_SEARCH // only for validation runs: the search tree is then the one of the original code
const bool USE_LOCAL_SEARCH = false;
#else
const bool USE_LOCAL_SEARCH = true;
#endif
// PERFORMANCE FIX 6 (revised 2026-11-24): the local search starts only once the search has visited this many nodes;
// a node count rather than a time, so that runs are reproducible.
const int LOCAL_SEARCH_AFTER_NODES = 1000;
bool localSearchDone = false;
double improve_incumbent(double cur);
mat xGreedy; double greedyValue = 0; // the greedy interdiction and its value: first start of the deferred local search
// PERFORMANCE FIX 6 (second revision): two local searches, from the greedy solution and from the best solution found.
// x, y, z hold the incumbent; improve_incumbent() changes them only when it finds a strictly better value.
double deferred_local_search(double incumbent) {
	std::vector<int> xs(n + 1), ys(n + 1), zs(n + 1);
	auto save = [&]() { for (int j = 1; j <= n; j++) { xs[j] = x[j]; ys[j] = y[j]; zs[j] = z[j]; } };
	auto restore = [&]() { for (int j = 1; j <= n; j++) { x[j] = xs[j]; y[j] = ys[j]; z[j] = zs[j]; } };
	save();
	for (int j = 1; j <= n; j++) x[j] = xGreedy[j];
	double g = improve_incumbent(greedyValue);
	if (g < incumbent) { incumbent = g; save(); } else restore();
	return improve_incumbent(incumbent);
}
double improve_incumbent(double cur) {
	if (!USE_LOCAL_SEARCH || n < 2) return cur;
	double t0 = (double)clock() / CLOCKS_PER_SEC, maxSeconds = 0.01 * TIME_LIMIT;
	long long evals = 0, maxEvals = 2LL * n + 20;
	std::vector<int> levelOf(n + 1), order;
	for (int i = 1; i <= n; i++) levelOf[nextinline[i]] = i;
	int used = 0;
	for (int j = 1; j <= n; j++) { xls[j] = x[j]; if (x[j]) used += a1[j]; if (Prod_Low[j] <= b2) order.push_back(j); }
	// most valuable items for the follower first: they are the most useful to interdict
	std::stable_sort(order.begin(), order.end(), [](int p, int q) {
		return profit_at(p, Prod_Up[p]) > profit_at(q, Prod_Up[q]); });
	p1_invalidate_all();
	auto out_of_budget = [&]() { return evals >= maxEvals || time_limit_reached || (double)clock() / CLOCKS_PER_SEC - t0 > maxSeconds; };
	auto accept = [&]() { // the rows hold the accepted candidate: save it as the incumbent
		reconstructSolutionP(n, b2);
		for (int j = 1; j <= n; j++) { x[j] = xls[j]; y[j] = yu[j]; z[j] = zu[j]; }
	};
	bool improved = true;
	while (improved && !out_of_budget()) {
		improved = false;
		for (int j : order) { // add one item
			if (out_of_budget()) break;
			if (xls[j] || a1[j] > b1 - used) continue;
			xls[j] = 1; p1_invalidate_from(levelOf[j]); evals++;
			double v = solveP1(n, b2, xls, 0);
			if (v < cur) { cur = v; used += a1[j]; improved = true; accept(); }
			else { xls[j] = 0; p1_invalidate_from(levelOf[j]); }
		}
		for (int k = (int)order.size() - 1; k >= 0 && !improved; k--) { // swap: least valuable interdicted item first
			int i = order[k];
			if (!xls[i]) continue;
			for (int j : order) {
				if (out_of_budget()) break;
				if (xls[j] || a1[j] > b1 - used + a1[i]) continue;
				int lv = std::min(levelOf[i], levelOf[j]);
				xls[i] = 0; xls[j] = 1; p1_invalidate_from(lv); evals++;
				double v = solveP1(n, b2, xls, 0);
				if (v < cur) { cur = v; used += a1[j] - a1[i]; improved = true; accept(); break; }
				xls[i] = 1; xls[j] = 0; p1_invalidate_from(lv);
			}
		}
	}
	p1_invalidate_all();
	printf("\n	Local search: %lld evaluations, start value %.15g\n", evals, cur);
	return cur;
}

// One cell Omega(i, Cup, C, f) of the single-item recursion of the original nonlinear omega. keep(y, f) and
// interdict(y, f) read Omega(i + 1, Cup, y, f) and Omega(i + 1, Cup - a_i, y, f) with the old boundary rules; LW is
// LowerW[i + 1] and T the steps of the rows of level i + 1.
template <class KV, class IV> inline double omega_single_cell(int item, int C, int f, int LW, const std::vector<int>& T, KV keep, IV interdict) {
	int l = Prod_Low[item], u = Prod_Up[item];
	if (l > C) return keep(C, f); // item cannot be packed
	double minimumVal = interdict(C, f);
	double maxVal = 0;
	if (l == u || f == 1) { // fixed quantity l, or fractional already used: skip, l or u
		double val1 = add_profit(keep(C - l, f), item, l);
		double val2 = keep(C, f);
		maxVal = (val1 > val2) ? val1 : val2;
		if (l != u) { double val3 = add_profit(keep(C - u, f), item, u); maxVal = (maxVal > val3) ? maxVal : val3; }
	}
	else { // fractional == 0: old loop over w in [max(l, C - LW), min(u, C)], starting from maxVal = 0
		int wlo = std::max(l, C - LW), whi = std::min(u, C);
		if (wlo <= l && l <= whi) maxVal = std::max(maxVal, add_profit(keep(C - l, 0), item, l)); // bounds keep fractional = 0
		if (wlo <= u && u <= whi) maxVal = std::max(maxVal, add_profit(keep(C - u, 0), item, u));
		for_step_starts(C, std::max(wlo, l + 1), std::min(whi, u - 1), T, [&](int w) {
			maxVal = std::max(maxVal, add_profit(keep(C - w, 1), item, w)); });
	}
	return (minimumVal < maxVal) ? minimumVal : maxVal;
}

struct OmegaRowRef { bool inf = false; const double* p = nullptr; RowPtr hold; }; // a row of level i + 2 (inf: its Cup < 0); hold keeps it alive (PERFORMANCE FIX 7)

// Omega which is used to compute a strong lower bound for the problem (two-item look-ahead version).
// PERFORMANCE FIX 3 (row-based, as in the linear code) with the step-start enumeration of the NONLINEAR PORT:
// Omega(i, Cup, ., .) is computed as a whole row over the lower-level capacity, once per (i, Cup) pair reached.
// The cases and boundary rules are those of the original nonlinear omega_sparse / omega_dense.
RowPtr omega_row(int i, int Cup) { // 1 <= i <= n, 0 <= Cup <= LowerUpper[i]; nullptr if the time limit is hit
	Cup = canonical_cup(i, Cup); // PERFORMANCE FIX 5: capacities allowing the same interdictions share one row
	long long key = (long long)i * (long long)(b1 + 1) + Cup;
	if (RowPtr cached = om_cache_get(key)) return cached; // PERFORMANCE FIX 7
	if (time_limit_reached) return nullptr;
	if ((double)(clock() - global_start) / CLOCKS_PER_SEC >= TIME_LIMIT) { time_limit_reached = true; return nullptr; }

	int L = std::min(b2, LowerW[i]);
	std::vector<double> row((size_t)(L + 1) * 2);

	if (i + 1 > n) {
		// ---- last item: single-item recursion; the next level is empty (value 0) ----
		int item = nextinline[i], a = a1[item];
		auto keepVal = [&](int Clow, int) -> double { return (Clow < 0) ? -INF : 0; };
		auto intVal = [&](int Clow, int) -> double { if (Cup - a < 0) return INF; return (Clow < 0) ? -INF : 0; };
		for (int C = 0; C <= L; C++)
			for (int f = 0; f <= 1; f++)
				row[(size_t)C * 2 + f] = omega_single_cell(item, C, f, 0, stepsSuffix[n + 1], keepVal, intVal);
		return om_cache_put(key, std::move(row)); // PERFORMANCE FIX 7
	}

	// ---- two items at once: item1 = level i, item2 = level i + 1; rows of level i + 2 are read ----
	int it1 = nextinline[i], it2 = nextinline[i + 1];
	int l1 = Prod_Low[it1], u1 = Prod_Up[it1], a_1 = a1[it1];
	int l2 = Prod_Low[it2], u2 = Prod_Up[it2], a_2 = a1[it2];
	bool last2 = (i + 2 > n);
	int LW2 = LowerW[i + 2]; // LowerW[n + 1] = 0
	auto getRow = [&](int CupX) -> OmegaRowRef {
		OmegaRowRef r; r.inf = (CupX < 0);
		if (!r.inf && !last2) { r.hold = omega_row(i + 2, std::min(CupX, LowerUpper[i + 2])); r.p = r.hold ? r.hold->data() : nullptr; }
		return r;
	};
	OmegaRowRef rBoth = getRow(Cup - a_1 - a_2); // both items interdicted
	OmegaRowRef rInt1 = getRow(Cup - a_1);       // item1 interdicted
	OmegaRowRef rInt2 = getRow(Cup - a_2);       // item2 interdicted
	OmegaRowRef rKeep = getRow(Cup);             // neither interdicted
	if (time_limit_reached) return nullptr;
	// Omega(i + 2, CupX, Clow, f) with the old boundary rules
	auto rv = [&](const OmegaRowRef& r, int Clow, int f) -> double {
		if (r.inf) return INF;
		if (Clow < 0) return -INF;
		if (last2) return 0;
		if (Clow > LW2) Clow = LW2;
		return r.p[(size_t)Clow * 2 + f];
	};
	const std::vector<int>& T2 = stepsSuffix[i + 2]; // steps of the rows of level i + 2
	const std::vector<int>& T1 = stepsSuffix[i + 1]; // contains T2, T2 + l2 and T2 + u2: steps after item2 is added
	// old loops "for (int w = wMin; w <= minval; w++)" packing one item (bounds: fractional 0, otherwise 1) on top of
	// row r, w in [max(l, C - LW2), min(u, C)]; start is the value without the item
	auto item_loop = [&](const OmegaRowRef& r, int it, int C, double start) -> double {
		int l = Prod_Low[it], u = Prod_Up[it]; double best = start;
		int wlo = std::max(l, C - LW2), whi = std::min(u, C);
		if (wlo <= l && l <= whi) best = std::max(best, add_profit(rv(r, C - l, 0), it, l));
		if (u != l && wlo <= u && u <= whi) best = std::max(best, add_profit(rv(r, C - u, 0), it, u));
		for_step_starts(C, std::max(wlo, l + 1), std::min(whi, u - 1), T2, [&](int w) {
			best = std::max(best, add_profit(rv(r, C - w, 1), it, w)); });
		return best;
	};

	for (int C = 0; C <= L; C++) {
		for (int f = 0; f <= 1; f++) {
			// Case 1: both items interdicted
			double both = rv(rBoth, C, f);
			// Case 2: item1 interdicted, item2 not interdicted
			double m1 = rv(rInt1, C, f);
			if (l2 <= C) {
				if (l2 == u2) m1 = std::max(m1, add_profit(rv(rInt1, C - l2, f), it2, l2));
				else if (f == 1) {
					m1 = std::max(m1, add_profit(rv(rInt1, C - l2, 1), it2, l2));
					m1 = std::max(m1, add_profit(rv(rInt1, C - u2, 1), it2, u2));
				}
				else m1 = item_loop(rInt1, it2, C, m1);
			}
			// Case 3: item1 not interdicted, item2 interdicted
			double m2 = rv(rInt2, C, f);
			if (l1 <= C) {
				if (l1 == u1) m2 = std::max(m2, add_profit(rv(rInt2, C - l1, f), it1, l1));
				else if (f == 1) {
					m2 = std::max(m2, add_profit(rv(rInt2, C - l1, 1), it1, l1));
					m2 = std::max(m2, add_profit(rv(rInt2, C - u1, 1), it1, u1));
				}
				else m2 = item_loop(rInt2, it1, C, m2);
			}
			// Case 4: neither item interdicted
			double m4 = rv(rKeep, C, f);
			if (l2 <= C) { // item2 only (in the original the "l2 == u2" test is not followed by "else")
				if (l2 == u2) m4 = std::max(m4, add_profit(rv(rKeep, C - l2, f), it2, l2));
				if (f == 1) {
					m4 = std::max(m4, add_profit(rv(rKeep, C - l2, 1), it2, l2));
					m4 = std::max(m4, add_profit_case4_u2(rv(rKeep, C - u2, 1), it2));
				}
				else m4 = item_loop(rKeep, it2, C, m4);
			}
			if (l1 <= C) { // item1 only, and item1 + item2
				if (l1 == u1) {
					int rem = C - l1;
					m4 = std::max(m4, add_profit(rv(rKeep, rem, f), it1, l1));
					if (l2 <= rem) m4 = std::max(m4, add_profit2(rv(rKeep, rem - l2, f), it1, l1, it2, l2));
				}
				if (f == 1) { // item1 at l1 or u1, then item2 at l2 or u2
					for (int w1 : { l1, u1 }) {
						int rem = C - w1;
						m4 = std::max(m4, add_profit(rv(rKeep, rem, 1), it1, w1));
						if (l2 <= rem) {
							m4 = std::max(m4, add_profit2(rv(rKeep, rem - l2, 1), it1, w1, it2, l2));
							m4 = std::max(m4, add_profit2(rv(rKeep, rem - u2, 1), it1, w1, it2, u2));
						}
					}
				}
				else { // old loop over w1 in [max(l1, C - LW2), min(u1, C)]
					int wlo = std::max(l1, C - LW2), whi = std::min(u1, C);
					// w1 at a bound (fractional stays 0): item1 alone, then item2 with l2 = u2 or with any w2 (inner loop)
					for (int k = 0; k < (l1 == u1 ? 1 : 2); k++) {
						int w1 = (k == 0) ? l1 : u1;
						if (w1 < wlo || w1 > whi) continue;
						int rem = C - w1;
						m4 = std::max(m4, add_profit(rv(rKeep, rem, 0), it1, w1));
						if (l2 > rem) continue;
						if (l2 == u2) { m4 = std::max(m4, add_profit2(rv(rKeep, rem - l2, 0), it1, w1, it2, l2)); continue; }
						int vlo = std::max(l2, rem - LW2), vhi = std::min(u2, rem);
						if (vlo <= l2 && l2 <= vhi) m4 = std::max(m4, add_profit2(rv(rKeep, rem - l2, 0), it1, w1, it2, l2));
						if (vlo <= u2 && u2 <= vhi) m4 = std::max(m4, add_profit2(rv(rKeep, rem - u2, 0), it1, w1, it2, u2));
						for_step_starts(rem, std::max(vlo, l2 + 1), std::min(vhi, u2 - 1), T2, [&](int w2) {
							m4 = std::max(m4, add_profit2(rv(rKeep, rem - w2, 1), it1, w1, it2, w2)); });
					}
					// w1 strictly between its bounds (fractional becomes 1): item1 alone ...
					int ilo = std::max(wlo, l1 + 1), ihi = std::min(whi, u1 - 1);
					for_step_starts(C, ilo, ihi, T2, [&](int w1) {
						m4 = std::max(m4, add_profit(rv(rKeep, C - w1, 1), it1, w1)); });
					// ... and item1 plus item2 at l2 or u2 (only when l2 <= C - w1); the rows are read at C - w1 - l2 and
					// C - w1 - u2, so their steps as functions of C - w1 are in T1
					for_step_starts(C, ilo, std::min(ihi, C - l2), T1, [&](int w1) {
						int rem = C - w1;
						m4 = std::max(m4, add_profit2(rv(rKeep, rem - l2, 1), it1, w1, it2, l2));
						m4 = std::max(m4, add_profit2(rv(rKeep, rem - u2, 1), it1, w1, it2, u2));
					});
				}
			}
			// the upper-level player minimises over the four interdiction choices
			row[(size_t)C * 2 + f] = std::min(std::min(both, m1), std::min(m2, m4));
		}
	}
	return om_cache_put(key, std::move(row)); // PERFORMANCE FIX 7
}

// Same interface and boundary rules as the old omega_dense/omega_sparse.
double omega(int i, int Cup, int Clow, int fractional) {
	if (time_limit_reached) return 0;
	if (dp_time_limit_exceeded()) { // PERFORMANCE FIX 1
		time_limit_reached = true;
		return 0;
	}
	if (Cup < 0) return INF; // Using a safe high value for "infinity"
	if (Clow < 0) return -INF;
	if (i > n) return 0;
	if (Clow > LowerW[i]) Clow = LowerW[i];
	if (Cup > LowerUpper[i]) Cup = LowerUpper[i];
	RowPtr row = omega_row(i, Cup);
	if (row == nullptr) return 0; // time limit reached
	if (Clow > b2) Clow = b2;      // rows cover Clow <= b2, the largest value bnb() asks for
	return (*row)[(size_t)Clow * 2 + fractional];
}


// The main Branch-and-Bound Function

double bnb(int item, int xu[], int CupLeft, double inheritedLB ) {
	nodes++;
	if (time_limit_reached) return inheritedLB;
	double elapsed = (double)(clock() - global_start) / CLOCKS_PER_SEC;
	if (elapsed >= TIME_LIMIT) {
		time_limit_reached = true;
		return inheritedLB;
	}
	// PERFORMANCE FIX 6 (revised): once the search is not trivial, improve the best solution found so far by local
	// search. It rebuilds the lower-level rows for its own candidates and leaves them invalid, so this node and the
	// nodes after it recompute their rows as usual.
	if (!localSearchDone && nodes >= LOCAL_SEARCH_AFTER_NODES) {
		localSearchDone = true;
		double improved = deferred_local_search(gub);
		if (improved < gub) gub = improved;
	}
	
	// Also always check if we have reached the final item as in that case we can compute the actual optimal solution
	if (item > n) {
		double bestsolution = solveP1(n, b2, xu, 0);
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
	// PERFORMANCE FIX 4: the lower-level row P(item - 1, ., 0) and the omega row Omega(item, CupLeft, ., 0) are the same
	// for every w, so they are fetched once per node (the loop used to read the clock and call solveP1() and omega() for
	// every w). Same values: omega(item, CupLeft, b2 - w, 0) is the row's entry at min(b2 - w, LowerW[item]).
	double locallower = 0;
	double nodeLB = inheritedLB;
	if (UpperW[item - 1] >= LowerW2[item]) {
		solveP1(item - 1, 0, xu, 0); // builds the missing lower-level rows up to item - 1
		RowPtr omRowP = omega_row(item, std::min(CupLeft, LowerUpper[item]));
		if (time_limit_reached || omRowP == nullptr) return nodeLB;
		const double* omRow = omRowP->data();
		int ClowMax = std::min(b2, LowerW[item]);
		for (int w = UpperW[item - 1]; w >= LowerW2[item]; w--) {
			int Clow = b2 - w; if (Clow > ClowMax) Clow = ClowMax;
			locallower = P1(item - 1, w, 0) + omRow[(size_t)Clow * 2];
			nodeLB = max(nodeLB, locallower);
			if (nodeLB >= gub) return nodeLB;
		}
	}
	// Depth first branch-and-bound always checking if the upper level problem is still feasible
	double bestSubtreeLB = INF;
#ifdef VALIDATION_NO_FIX8 // only for validation runs: branch like the original code
	if (a1[nextinline[item]] <= CupLeft) {
#else
	if (a1[nextinline[item]] <= CupLeft && Prod_Low[nextinline[item]] <= b2) { // PERFORMANCE FIX 8: never interdict an item the follower cannot pack
#endif
		xu[nextinline[item]] = 1;
		p1_invalidate_from(item); // PERFORMANCE FIX 3: lower-level rows >= item depend on this decision (replaces the memoLog undo)

		double lbTake = bnb(item + 1, xu, CupLeft - a1[nextinline[item]], nodeLB);
		bestSubtreeLB = min(bestSubtreeLB, lbTake);
		
		if (time_limit_reached) return min(bestSubtreeLB, nodeLB);
		// Makes sure that the solution above does not influence the solution below
	}

	if (xu[nextinline[item]] != 0) {
		xu[nextinline[item]] = 0;
		p1_invalidate_from(item); // PERFORMANCE FIX 3: the take branch changed rows >= item
	}
	double lbSkip = bnb(item + 1, xu, CupLeft, nodeLB);
	bestSubtreeLB = min(bestSubtreeLB, lbSkip);
	return max(nodeLB, bestSubtreeLB);
}

void procedure(void) {
	global_start = clock();
	// Need some sort of initialization to make sure we start with a completely empty solution
	memoBellman.assign(n + 1, std::vector<double>(b1 + 1, -1));
	// PERFORMANCE FIX 3: flat lower-level DP rows (were the nested vectors memoP1 and choiceZ)
	p1Val.assign((size_t)(n + 1) * (b2 + 1) * 2, 0);
	p1Choice.assign((size_t)(n + 1) * (b2 + 1) * 2, 0);
	p1_invalidate_all();
	gub = INF;
	avg_a1 = 0.0;
	avg_width = 0.0;
	avg_Prod_Low = 0.0;
	for (int i = 1; i <= n; i++) {
		ratio[i] = 0; 
		for (int j = std::max(1, Prod_Low[i]); j <= Prod_Up[i]; j++){
			if (ratio[i] < ((double)c[i] + (double)variable[i] * std::log(j) * j) / (double)j) {
				ratio[i] = ((double)c[i] + (double)variable[i] * std::log(j) * j) / (double)j;
				ratio_value[i] = j;
			}
		}
		x[i] = 0;
		y[i] = 0;
		z[i] = 0;
		xv[i] = 0;
		yu[i] = 0;
		zu[i] = 0;
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

	quickSort(ratio, nextinline, 1, n); // as in the original nonlinear code


	// Start of the actual bnb procedure
	
	// Precomputing some bounds for the w used in the lower bound computations
	UpperW[0]=0;
	LowerW[n + 1] = 0;
	LowerUpper[n + 1] = 0;
	for (int i = 1; i <= n; i++) {
		UpperW[i]=UpperW[i-1] + Prod_Up[nextinline[i]];
		if (UpperW[i] > b2) UpperW[i] = b2;
		LowerW[n + 1 - i]= LowerW[n + 2 -i] + Prod_Up[nextinline[n + 1 - i]];
		LowerUpper[n + 1 - i]= LowerUpper[n + 2 -i] + a1[nextinline[n + 1 - i]];
	}
	for (int i = 1; i <= n; i++) {
    LowerW2[i] = max(0, b2 - LowerW[i]);
	}
	build_step_lists();           // NONLINEAR PORT: steps of the rows (sums of 0, l or u units)
	build_canonical_capacities(); // PERFORMANCE FIX 5
	setup_omega_cache();          // PERFORMANCE FIX 7
	gub = greedyupper(); // Note that this also gives a first initial solution to the problem
	localSearchDone = false; // PERFORMANCE FIX 6 (revised): the local search runs inside bnb(), see there
	greedyValue = gub; for (int j = 1; j <= n; j++) xGreedy[j] = x[j];


	// Make sure we start with a clear slate for the memory
	p1_invalidate_all(); // PERFORMANCE FIX 3: rows are recomputed for the new interdiction array
	// PERFORMANCE FIX 3: omega rows are created on demand, so no dense/sparse table has to be chosen here.
	printf("\n	Omega table: rows per (item, upper capacity)\n");
	glb = bnb(1, xv, b1, 0);
	if (glb > gub) glb = gub; // Make sure that the initial x is only 0 values in order for this to work
	 std::vector<vector<double>>().swap(memoBellman); // Avoiding memory leaks
	printf("\n	Omega rows: %zu computed, %zu evicted, cache peak %.1f MB, budget %.1f MB\n", omComputed, omEvicted, omPeakBytes / 1048576.0, omBudget == (size_t)-1 ? -1.0 : omBudget / 1048576.0);
	std::vector<double>().swap(p1Val); std::vector<int>().swap(p1Choice); // PERFORMANCE FIX 3: free the DP rows
	omCache.clear(); omLru.clear(); std::vector<int>().swap(canonCup); // PERFORMANCE FIXES 5 and 7
	
	return;
}










