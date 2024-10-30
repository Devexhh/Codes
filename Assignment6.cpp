#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return 1.0/(1+x*x);
}

double trapezoidalRule(double a, double b, int n) {
    double h = (b-a)/n, x = (f(a)+f(b))/2.0;
    for (int i = 1; i < n; i++) {
        x += f(a+i*h);
    }
    return x*h;
}

double simpsonsOneThirdRule(double a, double b, int n) {
    if (n % 2 != 0) {
        cout << "Error Even Interval";
        return -1;
    }
    double h = (b-a)/n, x = f(a)+f(b);
    for (int i = 1; i < n; i += 2) {
            x += 4*f(a+i*h);
    }
    for (int i = 2; i < n-1; i += 2) {
        x += 2*f(a+i*h);
    }
    return h/3*x;
}

double simpsonsThreeEighthRule(double a, double b, int n) {
    if (n % 3 != 0) {
        cout << "Error Interval is not a multiple of 3";
        return -1;
    }
    double h = (b-a)/n, x = f(a)+f(b);
    for (int i = 1; i < n; i++) {
        if (i % 3 == 0) {
            x += 2*f(a+i*h);
        } else {
            x += 3*f(a+i*h);
        }
    }
    return 3*h/8*x;
}


int main() {
    double a = 0.0, b = 1.0;
    int intervals[] = {10,12,15,18,20};
    cout << fixed << setprecision(8);
    for (int n: intervals) {
        double approx = trapezoidalRule(a,b,n);
        double exact = atan(1)*4;
        double error = fabs(exact-approx);
        cout << "Trapezoidal Rule with n = " << n << " : " << approx << ", Absolute error = " << error << endl;
    }
    for (int n: intervals) {
        if (n % 2 == 0) {
            double approx = simpsonsOneThirdRule(a,b,n);
            double exact = atan(1)*4;
            double error = fabs(exact-approx);
            cout << "Simpson's 1/3 with n = " << n << " : " << approx << ", Absolute error = " << error << endl;
        }
    }
    for (int n: intervals) {
        if (n % 3 == 0) {
            double approx = simpsonsThreeEighthRule(a,b,n);
            double exact = atan(1)*4;
            double error = fabs(exact-approx);
            cout << "Simpson's 3/8 with n = " << n << " : " << approx << ", Absolute error = " << error << endl;
        }
    }
    return 0;
}
