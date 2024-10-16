#include <bits/stdc++.h>
using namespace std;

void gaussJacobi(vector<vector<double>>& a, vector<double>& b, int n, int max_iter = 100, double tol = 1e-5) {
    vector<double> x(n, 0), x_new(n, 0);
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += a[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / a[i][i];
        }
        double error = 0;
        for (int i = 0; i < n; i++) {
            error += fabs(x_new[i] - x[i]);
            x[i] = x_new[i];
        }
        if (error < tol) break;
    }
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

int main() {
    int n = 3;
    vector<vector<double>> a = {
        {5, -2, 3},
        {-3, 9, 1},
        {2, -1, -7}
    };
    vector<double> b = {-1, 2, 3};
    gaussJacobi(a, b, n);
    return 0;
}
