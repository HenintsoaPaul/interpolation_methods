#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void getData(int *pDim, float **pXd, float **pYd);

float **allocateMatrix(int dimension);

float *allocateVector(int dimension);

void logError(const char *string);

void gplot(int dim, float *xd, float *yd);

void plotSetup(FILE *GP, int xmin, int xmax, int ymin, int ymax);

float lagrange(float x, const float *xd, const float *yd, int dim);

void lambdaPicoY(int dim, const float *xd, const float *yd, float **ptabLambda, float **ptabPico, float **ptabY);

void createSystem(float ***matrixA, float **vectorB, float **vectorSolution, int dimension, const float *tabLambda,
                  const float *tabPico, const float *tabY);

void plotLinear(FILE *GP, int dim, float *xd, float *yd);

void plotLagrange(FILE *GP, int dim, float *xd, float *yd, float xmin, float xmax, float *pYmin, float *pYmax);

void plotSpline(FILE *GP, int dim, float *xd, float *yd, float xmin, float xmax);

void cholesky(int dim, float **matrixA, float *vectorB, float *vectorS);

void factorisationCholesky(int dim, float **A);

void solveTriangleInf(int dim, float **A, float *y);

void solveTriangleSup(int dim, float **A, float *x, const float *y);

int main() {
    /// Initialisation
    int dim = 7;
    float *xd = NULL, *yd = NULL;

    getData(&dim, &xd, &yd);

    /// Operations
    gplot(dim, xd, yd);

    return 0;
}

void gplot(int dim, float *xd, float *yd) {
    FILE *GP = popen("gnuplot -persist", "w");
    float xmin = xd[0], xmax = xd[dim - 1], ymin = 0, ymax = 0; // setup x range

    if (GP) {
        /// put data in reusable internal data blocks
        plotLinear(GP, dim, xd, yd);
        plotLagrange(GP, dim, xd, yd, xmin, xmax, &ymin, &ymax);
        plotSpline(GP, dim, xd, yd, xmin, xmax);

        /// plot setup
        int Ymin = ymin - 1, Ymax = ceil(ymax), Xmin = xmin, Xmax = xmax;   // setup x and y ranges
        plotSetup(GP, Xmin, Xmax, Ymin, Ymax);

        /// plot now
        /*fprintf(GP, "plot $fileData using 1:2 w linespoints lc'black' pointtype 7 t'Linear' \n");
        fprintf(GP, "plot $fileData using 1:2 w linespoints lc'black' pointtype 7 t'Linear'"
                    ", $lagrangeData using 1:2 w linespoints lc'red' pointtype 7 t'Lagrange'\n");*/
        fprintf(GP, "plot $fileData using 1:2 w linespoints lc'black' pointtype 7 t'Linear'"
                    ", $lagrangeData using 1:2 w linespoints lc'red' pointtype 7 t'Lagrange'"
                    ", $splineData using 1:2 w linespoints lc'blue' pointtype 7 t'Cubic spline'\n");

        /// send command to gnuplot and close pipe
        fflush(GP);
        fclose(GP);
    }
    else {
        printf("gnuplot not found ...\n");
    }
}

void createSystem(float ***matrixA, float **vectorB, float **vectorSolution, int dimension, const float *tabLambda,
                  const float *tabPico, const float *tabY) {
    // allocate space
    *vectorB = allocateVector(dimension);
    *vectorSolution = allocateVector(dimension);
    *matrixA = allocateMatrix(dimension);

    // add elements of the matrix
    (*matrixA)[0][0] = 2;
    for (int i = 1; i < dimension; ++i) {
        (*matrixA)[i][i] = 2;
        (*matrixA)[i - 1][i] = tabLambda[i - 1];
        (*matrixA)[i][i - 1] = tabPico[i - 1];

        for (int j = i + 1; j < dimension; ++j) {
            (*matrixA)[i - 1][j] = 0;   // null en haut a droite
        }
        for (int j = i - 1; j > -1; --j) {
            (*matrixA)[i][j] = 0;   // null en haut a gauche
        }
    }

    // init elements in vectorB and vectorSolution
    for (int i = 0; i < dimension; ++i) {
        vectorSolution[0][i] = 0;
        vectorB[0][i] = tabY[i];
    }
}

void lambdaPicoY(int dim, const float *xd, const float *yd, float **ptabLambda, float **ptabPico, float **ptabY) {
     *ptabLambda = allocateVector(dim - 1),
     *ptabPico = allocateVector(dim - 2),
     *ptabY = allocateVector(dim - 1);



    for (int i = 0; i < dim - 1; ++i) {
        float lambda_j = (*ptabLambda)[i + 1],
                pico_j_1 = 1 - lambda_j;
        (*ptabPico)[i] = pico_j_1;
    }
}

void plotSpline(FILE *GP, int dim, float *xd, float *yd, float xmin, float xmax) {
    // calculate lambda, pico and Y
    float *tabLambda, *tabPico, *tabY;
    lambdaPicoY(dim, xd, yd, &tabLambda, &tabPico, &tabY);
    for (int i = 0; i < dim + 1; ++i) {
        printf("i : %d\n", i);
        printf("Lambda : %f\n", tabLambda[i]);
        printf("Pico : %f\n", tabPico[i]);
        printf("Y : %f\n", tabY[i]);
        printf("---\n");
    }

    // create the system of matrix
    float **matrixA, *vectorB, *vectorS;
    createSystem(&matrixA, &vectorB, &vectorS, dim, tabLambda, tabPico, tabY);

    // resolve by cholesky
    cholesky(dim, matrixA, vectorB, vectorS);

    // add data blocks
    fprintf(GP, "$splineData << EOF\n");

    fprintf(GP, "EOF\n");
}

void cholesky(int dim, float **matrixA, float *vectorB, float *vectorS) {
    // factorisation - calculus of B
    factorisationCholesky(dim, matrixA);

    // resolution system tri inf -- calculus of y
    solveTriangleInf(dim, matrixA, vectorB);

    // resolution system tri sup -- calculus of x
    solveTriangleSup(dim, matrixA, vectorS, vectorB);
}

void solveTriangleSup(int dim, float **A, float *x, const float *y) {
    for (int i = dim - 1; i >= 0; --i) {
        float sum = 0;
        for (int j = i + 1; j < dim; ++j) {
            sum += A[j][i] * x[j];
        }
        x[i] = (y[i] - sum) / A[i][i];
    }
}

void solveTriangleInf(int dim, float **A, float *y) {
    for (int i = 0; i < dim; ++i) {
        float sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += A[i][j] * y[j];
        }
        y[i] = (y[i] - sum) / A[i][i];
    }
}

void factorisationCholesky(int dim, float **A) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < i; ++j) {
            float temp1 = 0;
            for (int k = 0; k < j - 1; ++k) {
                temp1 += A[i][k]  * A[j][k];
            }
            A[i][j] = (A[i][j] - temp1) / A[j][j];
        }
        float temp2 = 0;
        for (int k = 0; k < i - 1; ++k) {
            temp2 += pow(A[i][k], 2);
        }
        A[i][i] = sqrt(A[i][i] / temp2);
    }
}

float lagrange(float x, const float *xd, const float *yd, int dim) {
    float result = 0;
    for (int j = 0; j < dim; ++j) {

        float product = yd[j]; // y_j
        for (int k = 0; k < dim; ++k) {
            if ( j != k ) {

                float x_k = xd[k], x_j = xd[j],
                    P_j = (x - x_k) / (x_j - x_k);

                product *= P_j; // y_j * P_j(x_i)
            }
        }
        result += product;
    }
    return result;
}

void plotLagrange(FILE *GP, int dim, float *xd, float *yd, float xmin, float xmax, float *pYmin, float *pYmax) {
    fprintf(GP, "$lagrangeData << EOF\n");
    int np = 100; // nb of points to show in Lagrange's representation
    float x = xmin, dx = (xmax - xmin) / (np - 1);
    x -= dx;
    for (int i = 0; i < np; ++i) {
        x += dx;
        float y = lagrange(x, xd, yd, dim);
        fprintf(GP, "%f %f\n", x, y);

        if (*pYmin > y) *pYmin = y;
        if (*pYmax < y) *pYmax = y;
    }
    fprintf(GP, "EOF\n");
}

void plotLinear(FILE *GP, int dim, float *xd, float *yd) {
    fprintf(GP, "$fileData << EOF\n");
    for (int i = 0; i < dim; ++i) {
        fprintf(GP, "%f %f\n", xd[i], yd[i]); // data
    }
    fprintf(GP, "EOF\n");
}

void plotSetup(FILE *GP, int xmin, int xmax, int ymin, int ymax) {
    fprintf(GP, "set term wxt size 640,480\n");
    fprintf(GP, "set title 'Interpolation'\n");
    fprintf(GP, "set xlabel 'x'\n");
    fprintf(GP, "set ylabel 'y'\n");
    fprintf(GP, "set xzeroaxis\n");
    fprintf(GP, "set xrange [%d:%d]\n", xmin, xmax);
    fprintf(GP, "set yrange [%d:%d]\n", ymin, ymax);
    fprintf(GP, "set grid\n");
}

void getData(int *pDim, float **pXd, float **pYd) {
    /// open file if exists
    char filePath[] = "../data.txt";                                // must begin with '../'
    FILE *file = NULL;
    file = fopen(filePath, "r");

    if (file == NULL) {
        printf("File %s not found!\n", filePath);
        exit(2);
    }
    else {
        printf("File %s found!\n", filePath);

        /// load data from the file
        int dimension;                                              // load dimension from file
        fscanf(file, "%d", &dimension);

        float *vectorXd = allocateVector(dimension),                // allocate vector of x data
            *vectorYd = allocateVector(dimension);                  // allocate vector of y data
        for (int i = 0; i < dimension; ++i) {                       // load x;y in their vectors
            fscanf(file, "%f", &vectorXd[i]);
            fscanf(file, "%f", &vectorYd[i]);
        }

        /// close file
        fclose(file);

        /// send data from the file to main()
        *pDim = dimension;
        *pXd = vectorXd;
        *pYd = vectorYd;
    }
}

float **allocateMatrix(int dimension) {
    float **A = NULL;
    A = malloc(dimension * sizeof(float*));
    for (int i = 0; i < dimension; ++i) {
        A[i] = allocateVector(dimension);
    }
    if (A == NULL) logError("Cannot allocate matrix!");
    return A;
}

float *allocateVector(int dimension) {
    float *vector = NULL;
    vector = malloc(dimension * sizeof (float));
    if (vector == NULL) logError("Cannot allocate vector!");
    return vector;
}

void logError(const char *string) {
    printf("%s\n", string);
    exit(2);
}
