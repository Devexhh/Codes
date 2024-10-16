#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return cos(x) - x * exp(x);
}

void bisection(double a, double b, int maxIterations, double tolerance) {
    double mid;
    for (int i = 0; i < maxIterations; i++) {
        mid = (a + b) / 2.0;
        if (f(mid) == 0 || (b - a) / 2 < tolerance) {
            cout << "Root found: " << mid << endl;
            return;
        } else if (f(mid) * f(a) < 0) {
            b = mid;
        } else {
            a = mid;
        }
    }
    cout << "Root after " << maxIterations << " iterations: " << mid << endl;
}

int main() {
    double a = 0.0, b = 1.0;
    int maxIterations = 20;
    double tolerance = 1e-6;

    bisection(a, b, maxIterations, tolerance);

    return 0;
}