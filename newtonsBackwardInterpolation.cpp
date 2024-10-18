#include <bits/stdc++.h>
using namespace std;

void newtonForwardInterpolation(double x[], double y[], int n, double xp) {
    double h = x[1]-x[0];
    double diff[n][n];
    for (int i = 0; i < n; i++) {
        diff[i][0] = y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = n-1; i >= j; i--) {
            diff[i][j] = diff[i][j-1] - diff[i-1][j-1];
        }
    }
    double s = diff[n-1][0];
    double u = (xp - x[n-1])/h;
    double uterm = u;
    for (int i = 1; i < n; i++) {
        s += (uterm * diff[n-1][i]/i);
        uterm *= (u+i);
    }
    cout << "The interpolated value at x = " << xp << " is " << s << endl;
}

int main() {
    int n = 7;
    double x[] = {1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3};
    double y[] = {5.474, 6.050, 6.686, 7.389, 8.166, 9.025, 9.974};
    double xp = 2.25;
    newtonForwardInterpolation(x,y,n,xp);
    return 0;
}
