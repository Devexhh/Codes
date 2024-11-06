#include <bits/stdc++.h>
using namespace std;

double f(double x, double y) {
    return -2*x*y*y;
}

void eulerMethod(double x0, double y0, double h, double xf) {
    double x = x0;
    double y = y0;
    while (x < xf) {
        y += h*f(x,y);
        x += h;
    }
    cout << "The value of y at x = " << xf << " is " << y << endl;
}

int main () {
    double x0 = 1.0;
    double y0 = 1.0;
    double h = 0.1;
    double xf = 1.3;
    eulerMethod(x0, y0,h,xf);
    return 0;
}
