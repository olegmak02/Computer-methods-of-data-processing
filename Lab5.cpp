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

double* matrixVectorMultipication(double** m, double* v) {
    double* res = new double[2]{ 0,0 };
    for (int i = 0; i < 2; i++) {
        double a = 0;
        for (int j = 0; j < 2; j++) {
            a += m[i][j] * v[j];
        }
        res[i] = a;
        a = 0;
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

double* gradientCompute(double* x) {
    return new double[] { 2 * x[0] - 3, 2 * x[1] - 5 };
}

double vectorLengthCompute(double* v) {
    return sqrt(pow(v[0], 2) + pow(v[1], 2));
}

double vectorVectorMultipication(double* v1, double* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
}



double** matrixMatrixMultipication(double** m1, double** m2) {
    double** res = new double*[]{ new double[]{0,0}, new double[] {0,0} };
    double temp = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                temp += m1[i][k] * m2[k][j];
            }
            res[i][j] = temp;
            temp = 0;
        }
    }
    return res;
}

double** scalarMatrixMultipication(double a, double** m) {
    double** res = new double* [] { new double[] {0, 0}, new double[] {0, 0} };
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res[i][j] = m[i][j] * a;
        }
    }
    return res;
}

double* scalarVectorMultipication(double a, double* v) {
    return new double[2]{ a * v[0], a * v[1] };
}

double** matrixAddition(double** m1, double** m2) {
    double** res = new double* [] { new double[] {0, 0}, new double[] {0, 0} };
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res[i][j] = m1[i][j]+ m2[i][j];
        }
    }
    return res;
}

double** matrixFromVectorMultiplication(double* v1, double* v2) {
    double** res = new double* [] {new double[] {0, 0}, new double[] {0, 0}};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res[i][j] = v1[i] * v2[j];
        }
    }
    return res;
}

double* vectorMatrixMultiplication(double* v, double** m) {
    double* res = new double[2]{ 0,0 };
    for (int i = 0; i < 2; i++) {
        double temp = 0;
        for (int j = 0; j < 2; j++) {
            temp += v[j] * m[i][j];
        }
        res[i] = temp;
        temp = 0;
    }
    return res;
}

double* nextX(double* x, double** A, double* g) {
    double step = (2 * x[0] - 3) / (2 * (A[0][0] * g[0] + A[0][1] * g[1]));
    return vectorAddition(x, scalarVectorMultipication(-step, matrixVectorMultipication(A, g)));
}

double** nextA(double* nextX, double* x, double** A) {
    double* deltaX = vectorSubtraction(nextX, x);
    double** deltaA = matrixFromVectorMultiplication(deltaX, deltaX);
    double* deltaG = vectorSubtraction(gradientCompute(nextX), gradientCompute(x));
    deltaA = scalarMatrixMultipication(1/vectorVectorMultipication(deltaX, deltaG), deltaA);
    double** temp = matrixMatrixMultipication(matrixMatrixMultipication(A, matrixFromVectorMultiplication(deltaG, deltaG)), A);
    temp = scalarMatrixMultipication(-1/(vectorVectorMultipication(vectorMatrixMultiplication(deltaG, A), deltaG)), temp);
    return matrixAddition(A, matrixAddition(deltaA, temp));
}

double* DFPMethod(double* initX) {
    int k = 0;
    double e = 0.0001;
    double d = 0.0001;
    int steps = 10;
    double* x = new double[2]{ initX[0], initX[1] };
    double** A = new double* [] {new double[] {1,0}, new double[] {0,1}};

    while (k < steps) {
        double* nextStepX = nextX(x, A, gradientCompute(x));

        if ((vectorLengthCompute(nextStepX) - vectorLengthCompute(x)) <= e || vectorLengthCompute(gradientCompute(nextStepX)) <= d) {
            return nextStepX;
        }

        A = nextA(nextStepX, x, A);
        x = nextStepX;
        k++;
    }
    return x;
}

int main() {
    double* x = new double[] { -1, 6 };
    double* xTarg = DFPMethod(x);
cout << "x1 = " << xTarg[0] << endl;
cout << "x2 = " << xTarg[1] << endl;
cout << "Minimal value of function is " << func(xTarg) << endl;

}