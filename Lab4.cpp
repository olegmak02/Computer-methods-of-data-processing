#include <iostream>
#include <cmath> 

using namespace std;

double func(double x[2]) {
    return pow(1 - x[0], 2) + pow(2 - x[1], 2) - x[0] - x[1];
}


double* vectorAddition(double* v1, double* v2) {
    double* res = new double[2];
    double e = 0;
    for (int i = 0; i < 2; i++) {
        res[i] = v1[i] + v2[i];
    }
    return res;
}

double* nextX(double* x, double* grad) {
    double step =  -(2 * grad[0] * x[0] + 2 * grad[1] * x[1] - 3 * grad[0] - 5 * grad[1]) / (2 * (pow(grad[0], 2) + pow(grad[1], 2)));
    return vectorAddition(x, new double[] {grad[0] * step, grad[1] * step});
}

double* gradientCompute(double* x) {
    return new double[] { 2 * x[0] - 3, 2 * x[1] - 5 };
}

double vectorLengthCompute(double* v) {
    return sqrt(pow(v[0], 2) + pow(v[1], 2));
}

double* scalarVectorMultipication(double a, double* v) {
    return new double[2]{ a * v[0], a * v[1] };
}

double* nextD(double* nextX, double* x, double* d) {
    return vectorAddition(gradientCompute(nextX), scalarVectorMultipication(pow(vectorLengthCompute(gradientCompute(nextX)) / vectorLengthCompute(gradientCompute(x)), 2), d));
}

double* fletcherRivesMethod(double* initX) {
    int k = 0;
    double e = 0.0001;
    int steps = 2;
    double* x = new double[2]{ initX[0], initX[1] };
    double* d = gradientCompute(x);
    double* grad;

    while (k < steps) {
        grad = gradientCompute(x);

        if (vectorLengthCompute(grad) < e) {
            return x;
        }

        if (d[0]*grad[0]+ d[1] * grad[1] > 0 && k > 0) {
            return x;
        }

        double* nextStepX = nextX(x, d);

        if (abs(func(nextStepX) - func(x)) <= e) {
            return nextStepX;
        }
        
        d = nextD(nextStepX, x, d);
        x = nextStepX;
        k++;
    }
    return x;
}

int main() {
    double* x = new double[] { -1, 6 };
    double* xTarg = fletcherRivesMethod(x);
cout << "x1 = " << xTarg[0] << endl;
cout << "x2 = " << xTarg[1] << endl;
cout << "Minimal value of function is " << func(xTarg) << endl;

}