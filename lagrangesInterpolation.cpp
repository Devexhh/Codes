#include <bits/stdc++.h>
using namespace std;

double lagrangesInterpolation(double x[], double y[], int n, double xp) {
    double yp = 0;
    for (int i = 0; i < n; i++) {
        double p = y[i];
        for (int j = 0; j < n; j++) {
            if (i != j) {
                p = p * (xp - x[j])/(x[i] - x[j]);
            }
        }
        yp += p;
    }
    return yp;
}

int main () {
    int n; cout << "Enter n: "; cin >> n;
    double x[n-1], y[n-1];
    for (int i = 0; i < n-1; i++) {
        cout << "x[" << i << "]: ";
        cin >> x[i];
    }
    for (int i = 0; i < n-1; i++) {
        cout << "y[" << i << "]: ";
        cin >> y[i];
    }
    double xp; cout << "Enter xp: " ; cin >> xp;
    cout << "Interpolated Value at x = " << xp << " is " << lagrangesInterpolation(x,y,n,xp) << endl;
    return 0;
}
