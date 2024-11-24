#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "libMatriz/matriz.h"

void gaussElimination(double **U, int n) {
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = U[i][k] / U[k][k];
            for (int j = k; j <= n; j++) { // Incluindo a coluna b
                U[i][j] -= factor * U[k][j];
            }
        }
    }
}

void backSubstitution(double **U, double *x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = U[i][n]; // Pegando o valor de b (última coluna)
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 2000; // Tamanho da matriz
    double **U = criar_matriz(n);
    double *x = (double *)malloc(n * sizeof(double));

    if (x == NULL) {
        fprintf(stderr, "Erro ao alocar memória para o vetor x\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    double start_time, end_time;
    start_time = MPI_Wtime();

    if (rank == 0) {
        gaussElimination(U, n);
        end_time = MPI_Wtime();
        printf("Tempo de execução da eliminação de Gauss: %f segundos\n", end_time - start_time);
    }

    // Convertendo a matriz 2D para uma matriz 1D para MPI_Bcast
    double *U_flat = (double *)malloc(n * (n + 1) * sizeof(double));

    if (U_flat == NULL) {
        fprintf(stderr, "Erro ao alocar memória para a matriz U_flat\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            U_flat[i * (n + 1) + j] = U[i][j];
        }
    }

    MPI_Bcast(U_flat, n * (n + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Convertendo de volta para matriz 2D
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            U[i][j] = U_flat[i * (n + 1) + j];
        }
    }

    backSubstitution(U, x, n);

    if (rank == 0) {
        printf("Solução do sistema linear:\n");
        for (int i = 0; i < n; i++) {
            printf("x[%d] = %f\n", i, x[i]);
        }
    }

    end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Tempo total de execução: %f segundos\n", end_time - start_time);
    }

    liberarMatriz(U, n);
    free(x);
    free(U_flat);

    MPI_Finalize();
    return 0;
}
