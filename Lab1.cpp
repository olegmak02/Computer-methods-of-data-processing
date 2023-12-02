#include <iostream>
#include <cmath> 

using namespace std;

double func(double x[2]) {
    return pow(1 - x[0], 2) + pow(2 - x[1], 2) - x[0] - x[1];
}

double* nextX(double a, double* x, double* grad) {
    x[0] = x[0] - a * grad[0];
    x[1] = x[1] - a * grad[1];
    return x;
}

double stepCompute(double* x, double* grad) {
    return (2 * grad[0] * x[0] + 2 * grad[1] * x[1] - 3 * grad[0] - 5 * grad[1]) / (2 * (pow(grad[0], 2) + pow(grad[1], 2)));
}

double vectorLengthCompute(double* v) {
    return sqrt(pow(v[0], 2) + pow(v[1], 2));
}

double* gradientCompute(double* x) {
    return new double[2]{ 2 * x[0] - 3, 2 * x[1] - 5 };
}

double* cauchyMethod(double* initX) {
    int k = 0;
    double e1 = 0.0001;
    double e2 = 0.0001;
    int steps = 10;
    double a = 0;
    double* x = new double[2]{ initX[0], initX[1] };

    double* grad = new double[2]{ 0, 0 };

    while (k < steps) {
        grad = gradientCompute(x);

        if (vectorLengthCompute(grad) <= e1) {
            return x;
        }

        a = stepCompute(x, grad);

        double* nextStepX = nextX(a, x, grad);

        if ((vectorLengthCompute(nextStepX) - vectorLengthCompute(x) / vectorLengthCompute(x)) <= e2) {
            return nextStepX;
        }

        x = nextStepX;
        k++;
    }

    return x;
}

int main() {
    double* x = new double[2]{ -1, 6 };
    double* xTarg = cauchyMethod(x);
    cout << "x1 = " << xTarg[0] << endl;
    cout << "x2 = " << xTarg[1] << endl;
    cout << "Minimal value of function is " << func(xTarg) << endl;
}
