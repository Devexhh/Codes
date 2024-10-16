#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return cos(x) - x * exp(x);
}

void regulaFalsi(double a, double b, int maxIterations, double tolerance) {
    double c;
    for (int i = 0; i < maxIterations; i++) {
        c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        if (f(c) == 0 || fabs(f(c)) < tolerance) {
            cout << "Root found: " << c << endl;
            return;
        } else if (f(c) * f(a) < 0) {
            b = c;
        } else {
            a = c;
        }
    }
    cout << "Root after " << maxIterations << " iterations: " << c << endl;
}

int main() {
    double a = 0.0, b = 1.0;
    int maxIterations = 20;
    double tolerance = 1e-6;

    regulaFalsi(a, b, maxIterations, tolerance);

    return 0;
}