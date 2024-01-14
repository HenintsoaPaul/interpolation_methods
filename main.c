#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void getData(int *pDim, float **pXd, float **pYd);

float *allocateVector(int dimension);

void logError(const char *string);

void gplot(int dim, float *xd, float *yd);

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

    /// setup x range
    float xmin = xd[0], xmax = xd[dim - 1],
        ymin = 0, ymax = 0;

    if (GP) {
    /// put data in reusable internal data blocks
        // plot all the data in the file
        fprintf(GP, "$fileData << EOF\n");
        for (int i = 0; i < dim; ++i) {
            fprintf(GP, "%f %f\n", xd[i], yd[i]); // data
        }
        fprintf(GP, "EOF\n");

        // plot Lagrange
        // ...

    /// plot setup
        ymin = -2, ymax = 3;
        fprintf(GP, "set term wxt size 640,480\n");
        fprintf(GP, "set title 'Interpolation'\n");
        fprintf(GP, "set xlabel 'x'\n");
        fprintf(GP, "set ylabel 'y'\n");
        fprintf(GP, "set xzeroaxis\n");
        fprintf(GP, "set xrange [%f:%f]\n", xmin, xmax);
        fprintf(GP, "set yrange [%f:%f]\n", ymin, ymax);
        fprintf(GP, "set grid\n");

    /// plot now
        fprintf(GP, "plot $fileData using 1:2 w linespoints lc'black' pointtype 7 t'Linear' \n");

    /// send commmand to gnuplot and close pipe
        fflush(GP);
        fclose(GP);
    } else {
        printf("gnuplot not found ...\n");
    }
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
