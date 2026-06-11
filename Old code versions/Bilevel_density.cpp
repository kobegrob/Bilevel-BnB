#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>

#define maxn 10000
#define maxb 1500000

typedef int mat[maxn + 1];

int n, b1, b2;
mat c, a1, Prod_Up, Prod_Low, variable;

const int INF = 1000000000;

std::vector<int> omegaMemo;

inline int& omg_at(int i, int Cup, int Clow, int f) {
    return omegaMemo[
        (((size_t)i * (b1 + 1) + Cup) * (b2 + 1) + Clow) * 2 + f
    ];
}

void input() {
    char filename[256];
    FILE* fp;

    printf("\n\n\tGive the problem name: ");
    scanf("%255s", filename);
    strcat(filename, ".knp");

    fp = fopen(filename, "r");
    if (!fp) {
        printf("Could not open file %s\n", filename);
        exit(1);
    }

    fscanf(fp, "%d %d %d", &n, &b1, &b2);

    if (n > maxn) {
        printf("Too many objects\n");
        exit(1);
    }

    if (b1 > maxb || b2 > maxb) {
        printf("Capacities are too large\n");
        exit(1);
    }

    for (int i = 1; i <= n; i++) {
        fscanf(fp, "%d %d %d %d %d",
               &a1[i],
               &c[i],
               &Prod_Up[i],
               &Prod_Low[i],
               &variable[i]);
    }

    fclose(fp);
}

int omega(int i, int Cup, int Clow, int fractional) {
    if (Cup < 0) return INF;
    if (Clow < 0) return -INF;
    if (i > n) return 0;

    int& memo = omg_at(i, Cup, Clow, fractional);
    if (memo != -1) return memo;

    int result;

    if (Prod_Low[i] > Clow) {
        result = omega(i + 1, Cup, Clow, fractional);
    }
    else if (Prod_Low[i] == Prod_Up[i]) {
        int minimumVal = omega(i + 1, Cup - a1[i], Clow, fractional);

        int val1 = omega(i + 1, Cup, Clow - Prod_Low[i], fractional)
                 + c[i] + Prod_Low[i] * variable[i];

        int val2 = omega(i + 1, Cup, Clow, fractional);

        int maxVal = std::max(val1, val2);
        result = std::min(minimumVal, maxVal);
    }
    else if (fractional == 0) {
        int minimumVal = omega(i + 1, Cup - a1[i], Clow, fractional);
        int maxVal = 0;

        int maxWeight = std::min(Prod_Up[i], Clow);

        for (int w = Prod_Low[i]; w <= maxWeight; w++) {
            int nextFractional = (w == Prod_Up[i]) ? 0 : 1;

            int currentVal = omega(i + 1, Cup, Clow - w, nextFractional)
                           + c[i] + w * variable[i];

            maxVal = std::max(maxVal, currentVal);
        }

        result = std::min(minimumVal, maxVal);
    }
    else {
        int minimumVal = omega(i + 1, Cup - a1[i], Clow, fractional);

        int val1 = omega(i + 1, Cup, Clow - Prod_Low[i], fractional)
                 + c[i] + Prod_Low[i] * variable[i];

        int val2 = omega(i + 1, Cup, Clow, fractional);

        int maxVal = std::max(val1, val2);
        result = std::min(minimumVal, maxVal);
    }

    memo = result;
    return result;
}

void computeFullOmega() {
    for (int i = n; i >= 1; i--) {
        for (int Cup = 0; Cup <= b1; Cup++) {
            for (int Clow = 0; Clow <= b2; Clow++) {
                omega(i, Cup, Clow, 0);
                omega(i, Cup, Clow, 1);
            }
        }
    }
}

long long countMinusOnes() {
    long long count = 0;

    for (int i = 1; i <= n; i++) {
        for (int Cup = 0; Cup <= b1; Cup++) {
            for (int Clow = 0; Clow <= b2; Clow++) {
                for (int f = 0; f <= 1; f++) {
                    if (omg_at(i, Cup, Clow, f) == -1) {
                        count++;
                    }
                }
            }
        }
    }

    return count;
}

void computeDensity() {
    long long equalPairs = 0;
    long long totalPairs = 0;

    for (int i = 1; i <= n; i++) {
        for (int Cup = 0; Cup <= b1; Cup++) {
            for (int Clow = 0; Clow <= b2; Clow++) {
                for (int f = 0; f <= 1; f++) {
                    int v = omg_at(i, Cup, Clow, f);

                    if (Cup + 1 <= b1) {
                        int v2 = omg_at(i, Cup + 1, Clow, f);
                        totalPairs++;
                        if (v == v2) equalPairs++;
                    }

                    if (Clow + 1 <= b2) {
                        int v2 = omg_at(i, Cup, Clow + 1, f);
                        totalPairs++;
                        if (v == v2) equalPairs++;
                    }
                }
            }
        }
    }

    double density = totalPairs == 0
        ? 0.0
        : (double)equalPairs / (double)totalPairs;

    printf("\n\tOmega equal adjacent pairs = %lld / %lld",
           equalPairs, totalPairs);

    printf("\n\tOmega density = %.10f\n", density);
}

int main() {
    input();

    size_t states = (size_t)(n + 1) * (b1 + 1) * (b2 + 1) * 2;
    double gb = states * sizeof(int) / (1024.0 * 1024.0 * 1024.0);

    printf("\n\tOmega states: %zu", states);
    printf("\n\tMemory needed: %.4f GB\n", gb);

    omegaMemo.assign(states, -1);

    clock_t start = clock();

    computeFullOmega();

    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;

    long long missing = countMinusOnes();

    printf("\n\tUncomputed omega entries = %lld", missing);
    printf("\n\tComputation time = %.2f seconds\n", elapsed);

    computeDensity();

    return 0;
}







