#include <iostream>
#include <cmath> 

using namespace std;

double func(double x[2]) {
    return pow(1 - x[0], 2) + pow(2 - x[1], 2) - x[0] - x[1];
}

double* matrixVectorMultiply(double matrix[2][2], double* vector) {
    double* res = new double[sizeof(matrix) / sizeof(matrix[0])];
    double e = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            e += matrix[i][j] * vector[j];
        }
        res[i] = e;
        e = 0;
    }
    return res;
}

double* vectorSubtraction(double* v1, double* v2) {
    double* res = new double[2];
    double e = 0;
    for (int i = 0; i < 2; i++) {
        res[i] = v1[i] - v2[i];
    }
    return res;
}

double* nextX(double matrix[2][2], double* x, double* grad) {
    return vectorSubtraction(x, matrixVectorMultiply(matrix, grad));
}

double* gradientCompute(double* x) {
    return new double[] { 2 * x[0] - 3, 2 * x[1] - 5 };
}

double* newtoneMethod(double* initX) {
    int k = 0;
    double e1 = 0.0001;
    double e2 = 0.0001;
    int steps = 10;
    double reverseMatrix[2][2] = {{0.5, 0}, {0, 0.5}};
    double* x = new double[2]{ initX[0], initX[1] };

    double* grad = new double[2]{ 0, 0 };

    while (k < steps) {
        grad = gradientCompute(x);

        double* nextStepX = nextX(reverseMatrix, x, grad);

        if (abs(func(nextStepX) - func(x)) <= e1) {
            return nextStepX;
        }

        x = nextStepX;
        k++;
    }

    return x;
}

int main() {
    double* x = new double[] { -1, 6 };
    double* xTarg = newtoneMethod(x);
    cout << "x1 = " << xTarg[0] << endl;
      cout << "x2 = " << xTarg[1] << endl;
	  cout << "Minimal value of function is " << func(xTarg) << endl;

}