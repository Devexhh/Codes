#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return cos(x) - x * exp(x);
}

void secant(double x0, double x1, int maxIterations, double tolerance) {
    double x2;
    for (int i = 0; i < maxIterations; i++) {
        x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
        if (fabs(f(x2)) < tolerance) {
            cout << "Root found: " << x2 << endl;
            return;
        }
        x0 = x1;
        x1 = x2;
    }
    cout << "Root after " << maxIterations << " iterations: " << x2 << endl;
}

int main() {
    double x0 = 0.0, x1 = 1.0;
    int maxIterations = 10;
    double tolerance = 1e-6;

    secant(x0, x1, maxIterations, tolerance);

    return 0;
}