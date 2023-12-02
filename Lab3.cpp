#include <iostream>
#include <cmath> 

using namespace std;

double func(double x[2]) {
    return pow(1 - x[0], 2) + pow(2 - x[1], 2) - x[0] - x[1];
}

double* matrixVectorMultiply(double** matrix, double* vector) {
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

double** reverseMatrix(double** matrix) {
    double** rev = new double* [] {new double[] {0, 0}, new double[] {0, 0}};
    double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            rev[i][j] = pow(-1, i+j) * matrix[1 - i][1 - j] / det;
        }
    }
    return rev;
}

double* nextX(double H[2][2], double* x, double a, double* grad) {
    double** temp = new double*[2]{new double[]{ 0, 0 }, new double[]{0, 0}};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (i == j) {
                temp[i][j] = H[i][j] + a;
            } else {
                temp[i][j] = H[i][j];
            }
        }
    }
    return vectorSubtraction(x, matrixVectorMultiply(reverseMatrix(temp), grad));
}

double* gradientCompute(double* x) {
    return new double[] { 2 * x[0] - 3, 2 * x[1] - 5 };
}

double vectorLengthCompute(double* v) {
    return sqrt(pow(v[0], 2) + pow(v[1], 2));
}

double* marquardtMethod(double* initX) {
    int k = 0;
    double e = 0.0001;
    int steps = 100;
    double a = 10000;
    double H[2][2] = { {2, 0}, {0, 2} };
    double* x = new double[2]{ initX[0], initX[1] };

    double* grad = new double[2]{ 0, 0 };

    while (k < steps) {
        grad = gradientCompute(x);

        if (vectorLengthCompute(grad) < e) {
            return x;
        }

        double* nextStepX = nextX(H, x, a, grad);

        if (abs(func(nextStepX) - func(x)) <= e) {
            return nextStepX;
        }

        if (func(nextStepX) < func(x)) {
            a /= 2;
        } else {
            a *= 2;
        }
        
        x = nextStepX;
        k++;
    }

    return x;
}

int main() {
    double* x = new double[] { -1, 6 };
    double* xTarg = marquardtMethod(x);
	  cout << "x1 = " << xTarg[0] << endl;
	  cout << "x2 = " << xTarg[1] << endl;
	  cout << "Minimal value of function is " << func(xTarg) << endl;

}