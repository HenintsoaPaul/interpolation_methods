#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void getData(int *pDim, float **pXd, float **pYd);

float *allocateVector(int dimension);

void logError(const char *string);

void gplot(int dim, float *xd, float *yd);

void plotSetup(FILE *GP, int xmin, int xmax, int ymin, int ymax);

float lagrange(float x, float *xd, float *yd, int dim);

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
        // plot all the data in the file
        fprintf(GP, "$fileData << EOF\n");
        for (int i = 0; i < dim; ++i) {
            fprintf(GP, "%f %f\n", xd[i], yd[i]); // data
        }
        fprintf(GP, "EOF\n");

        // plot Lagrange
        fprintf(GP, "$lagrangeData << EOF\n");
        int np = 100; // nb of points to show in Lagrange's representation
        float x = xmin, dx = (xmax - xmin) / (np - 1);
        x -= dx;
        for (int i = 0; i < np; ++i) {
            x += dx;
            float y = lagrange(x, xd, yd, dim);
            fprintf(GP, "%f %f\n", x, y);

            if (ymin > y) ymin = y;
            if (ymax < y) ymax = y;
        }
        fprintf(GP, "EOF\n");



    /// plot setup
        int Ymin = ymin - 1, Ymax = ceil(ymax), Xmin = xmin, Xmax = xmax;   // setup x and y ranges
        plotSetup(GP, Xmin, Xmax, Ymin, Ymax);



    /// plot now
//        fprintf(GP, "plot $fileData using 1:2 w linespoints lc'black' pointtype 7 t'Linear' \n");
        fprintf(GP, "plot $fileData using 1:2 w linespoints lc'black' pointtype 7 t'Linear'"
                    ", $lagrangeData using 1:2 w linespoints lc'red' pointtype 7 t'Lagrange'\n");



    /// send command to gnuplot and close pipe
        fflush(GP);
        fclose(GP);
    }
    else {
        printf("gnuplot not found ...\n");
    }
}

float lagrange(float x, float *xd, float *yd, int dim) {
    float result = 0;
    for (int i = 0; i < dim; ++i) {

        float product = yd[i]; // y_j
        for (int j = 0; j < dim; ++j) {
            if ( i != j ) {

                float x_k = xd[j], x_j = xd[i],
                    P_j = (x - x_k) / (x_j - x_k);
                product *= P_j; // y_j * P_j(x_i)
            }
        }
        result += product;
    }
    return result;
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

float *allocateVector(int dimension) {
    float *vector = NULL;
    vector = malloc(dimension * sizeof (float ));
    if (vector == NULL) logError("Cannot allocate vector!");
    return vector;
}

void logError(const char *string) {
    printf("%s\n", string);
    exit(2);
}
