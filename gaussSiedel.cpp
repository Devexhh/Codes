#include <bits/stdc++.h>
using namespace std;

void gaussSeidel(vector<vector<float> >& a, vector<float>& b, int n, int max_iter = 100, float tol = 1e-5) {
    vector<float> x(n, 0);
    for (int iter = 0; iter < max_iter; iter++) {
        float error = 0;
        for (int i = 0; i < n; i++) {
            float sum = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += a[i][j] * x[j];
                }
            }
            float x_new = (b[i] - sum) / a[i][i];
            error += fabs(x_new - x[i]);
            x[i] = x_new;
        }
        if (error < tol) break;
    }
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

int main() {
    int n = 3;
    vector<vector<float> > a = {
        {5, -2, 3},
        {-3, 9, 1},
        {2, -1, -7}
    };
    vector<float> b = {-1, 2, 3};

    gaussSeidel(a, b, n);

    return 0;
}

