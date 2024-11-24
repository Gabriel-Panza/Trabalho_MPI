#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "libMatriz/matriz.h"

void gaussElimination(double **U, int n, int rank, int size) {
    for (int k = 0; k < n; k++) {
        if (k % size == rank) {
            for (int j = k + 1; j <= n; j++) {
                U[k][j] /= U[k][k];
            }
            U[k][k] = 1.0;
        }
        MPI_Bcast(U[k], n + 1, MPI_DOUBLE, k % size, MPI_COMM_WORLD);
        for (int i = k + 1; i < n; i++) {
            if (i % size == rank) {
                for (int j = k + 1; j <= n; j++) {
                    U[i][j] -= U[i][k] * U[k][j];
                }
                U[i][k] = 0.0;
            }
        }
    }
}

void backSubstitution(double **U, double *x, int n, int rank, int size) {
    for (int i = n - 1; i >= 0; i--) {
        if (i % size == rank) {
            x[i] = U[i][n];
            for (int j = i + 1; j < n; j++) {
                x[i] -= U[i][j] * x[j];
            }
        }
        MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 2000; // Tamanho da matriz
    double **U = criar_matriz(n);
    double x[n];

    double start_time, end_time;
    start_time = MPI_Wtime();

    gaussElimination(U, n, rank, size);
    backSubstitution(U, x, n, rank, size);


    if (rank == 0) {
        printf("Solução do sistema linear:\n");
        for (int i = 0; i < n; i++) {
            printf("x[%d] = %f\n", i, x[i]);
        }
    }

    end_time = MPI_Wtime();
    if (rank == 0)
        printf("Tempo total de execução: %f segundos\n", end_time - start_time);

    liberarMatriz(U, n);

    MPI_Finalize();
    return 0;
}