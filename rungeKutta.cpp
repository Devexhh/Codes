#include <bits/stdc++.h>
using namespace std;

double f(double x, double y) {
    return (x*x)+(y*y);
}

void rungeKutta(double x0, double y0, double h, double xf) {
    double x = x0;
    double y = y0;
    while (x < xf) {
        double k1 = h*f(x,y);
        double k2 = h*f(x+h/2,y+k1/2);
        double k3 = h*f(x+h/2,y+k2/2);
        double k4 = h*f(x+h,y+k3);

        y += (k1+2*k2+2*k3+k4)/6;
        x += h;
    }
    cout << "The value of y at x = " << xf << " is " << y << endl;
}

int main () {
    double x0 = 1.0;
    double y0 = 1.2;
    double h = 0.05;
    double xf = 1.05;
    rungeKutta(x0, y0,h,xf);
    return 0;
}
