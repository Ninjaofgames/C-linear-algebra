#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mt19937.h"

void identityMatrix(int dimension, double identfierMatrix[dimension][dimension]){
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            identfierMatrix[i][j] = (i == j) ? 1 : 0;
        }
    }
}

void nullMatrix(int x, int y, double nullMatrix[y][x]){
    for (int i = 0; i < y; i++){
        for (int j = 0; j < x; j++){
            nullMatrix[i][j] = 0;
        }
    }
}

void matRand(int x, int y, double matrix[y][x]) {
    static int seeded = 0;
    if (!seeded) {
        init_genrand64((uint64_t)time(NULL));
        seeded = 1;
    }
    for (int i = 0; i < y; i++) {
        for (int j = 0; j < x; j++) {
            matrix[i][j] = genrand64_real1();
        }
    }
}

int equiMatrix(int x_m1, int x_m2, int y_m1, int y_m2, double matrix1[y_m1][x_m1], double matrix2[y_m2][x_m2]){
    if(x_m1 != x_m2 || y_m1 != y_m2){
        return 0;
    }
    for(int i = 0; i < y_m1; i++){
        for(int j = 0; j < x_m1; j++){
            if(matrix1[i][j] != matrix2[i][j]){
                return 0;
            }
        }
    }
    return 1;
}

void addMatrix(int x_m1, int x_m2, int y_m1, int y_m2, double matrix1[y_m1][x_m1], double matrix2[y_m2][x_m2], double matrixResult[][ (x_m1 > x_m2) ? x_m1 : x_m2 ]){
    int maxX = (x_m1 > x_m2) ? x_m1 : x_m2;
    int maxY = (y_m1 > y_m2) ? y_m1 : y_m2;
    for (int i = 0; i < maxY; i++){
        for (int j = 0; j < maxX; j++){
            double val1 = (i < y_m1 && j < x_m1) ? matrix1[i][j] : 0;
            double val2 = (i < y_m2 && j < x_m2) ? matrix2[i][j] : 0;
            matrixResult[i][j] = val1 + val2;
        }
    }
}

void scalarMultiplication(int x, int y, double scalar, double matrix[y][x]){
    for (int i = 0; i < y; i++){
        for (int j = 0; j < x; j++){
            matrix[i][j] *= scalar;
        }
    }
}

int checkTriangularMatrix(int dimension, double matrix[dimension][dimension], char type[]){
    int isLower = 1, isUpper = 1;

    for (int i = 0; i < dimension; i++) {
        for (int j = i + 1; j < dimension; j++) {
            if (matrix[i][j] != 0) {
                isLower = 0;
                break;
            }
        }
    }

    for (int i = 1; i < dimension; i++) {
        for (int j = 0; j < i; j++) {
            if (matrix[i][j] != 0) {
                isUpper = 0;
                break;
            }
        }
    }

    if (isLower && isUpper){
        strcpy(type, "diagonal");
        return 1;
    }else if (isLower){
        strcpy(type, "lower");
        return 1;
    }else if (isUpper){
        strcpy(type, "upper");
        return 1;
    }else{
        strcpy(type, "none");
        return 0;
    }
}

void transpose(int x, int y, double matrix[y][x], double result[x][y]){
    for (int i = 0; i < y; i++){
        for (int j = 0; j < x; j++){
            result[j][i] = matrix[i][j];
        }
    }
}

double matrixTrace(int dimension, double matrix[dimension][dimension]){
    double trace = 0;
    for (int i = 0; i < dimension; i++){
        trace += matrix[i][i];
    }
    return trace;
}

void getMinor(int dimension, double srcMatrix[dimension][dimension], double minMatrix[dimension-1][dimension-1], int excludeX, int excludeY){
    int y = 0;
    for (int i = 0; i < dimension; i++){
        if(i == excludeY) continue;
        int x = 0;
        for (int j = 0; j < dimension; j++){
            if(j == excludeX) continue;
            minMatrix[y][x] = srcMatrix[i][j];
            x++;
        }
        y++;
    }
}

double determinant(int dimension, double matrix[dimension][dimension]){
    if(dimension == 1) return matrix[0][0];
    if (dimension == 2){
        return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
    }
    double det = 0.0;
    double minor[dimension - 1][dimension - 1];
    for (int i = 0; i < dimension; i++){
        getMinor(dimension, matrix, minor, i, 0);
        double sign = (i % 2 == 0) ? 1.0 : -1.0;
        double cofactor = sign * matrix[0][i] * determinant(dimension - 1, minor);
        det += cofactor;
    }
    return det;
}

void mulMatrix(int x_m1, int x_m2, int y_m1, int y_m2, double matrix1[y_m1][x_m1], double matrix2[y_m2][x_m2], double resultMatrix[y_m2][x_m1]){
    if (y_m1 != x_m2){
        return;
    }
    for (int i = 0; i < y_m1; i++){
        for (int j = 0; j < x_m2; j++){
            resultMatrix[i][j] = 0;
            for (int k = 0; k < x_m1; k++){
                resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

void matPrint(int x, int y, double matrix[y][x]){
    for (int i = 0; i < y; i++){
        for (int j = 0; j < x; j++){
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(){
    double matrix[3][3] = {
        {2, 3, 1},
        {1, -1, 2},
        {7, 3, 8}
    };
    printf("%lf", determinant(3, matrix));
}