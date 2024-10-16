#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return x * exp(x) - 1;
}

double df(double x) {
    return exp(x) + x * exp(x);
}

void newtonRaphson(double x0, int maxIterations, double tolerance) {
    double x1;
    for (int i = 0; i < maxIterations; i++) {
        x1 = x0 - f(x0) / df(x0);
        if (fabs(f(x1)) < tolerance) {
            cout << "Root found: " << x1 << endl;
            return;
        }
        x0 = x1;
    }
    cout << "Root after " << maxIterations << " iterations: " << x1 << endl;
}

int main() {
    double x0 = 0.5;
    int maxIterations = 10;
    double tolerance = 1e-6;

    newtonRaphson(x0, maxIterations, tolerance);

    return 0;
}