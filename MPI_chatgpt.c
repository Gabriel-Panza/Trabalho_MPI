#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 3 // Tamanho da matriz (NxN)

void printMatrix(double matrix[N][N+1], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double matrix[N][N+1] = {
        {0, 2, 5, 1},
        {2, 1, 1, 1},
        {3, 1, 0, 2},
    };

    if (rank == 0) {
        printf("Initial Matrix:\n");
        printMatrix(matrix, N, N + 1);
    }

    for (int k = 0; k < N; k++) {
        // Broadcast da linha pivô
        MPI_Bcast(&matrix[k], N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Redução nas outras linhas
        for (int i = rank; i < N; i += size) {
            if (i > k) {
                double factor = matrix[i][k] / matrix[k][k];
                for (int j = k; j < N + 1; j++) {
                    matrix[i][j] -= factor;
                }
            }
        }
        // Reunir as linhas atualizadas de volta ao processo raiz
        for (int i = rank; i < N; i += size) {
            MPI_Gather(matrix[i], N + 1, MPI_DOUBLE, matrix[i], N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // Processo raiz faz a substituição retroativa
    double result[N];
    if (rank == 0) {
        for (int i = N - 1; i >= 0; i--) {
            result[i] = matrix[i][N];
            for (int j = i + 1; j < N; j++) {
                result[i] -= matrix[i][j] * result[j];
            }
            result[i] /= matrix[i][i];
        }

        printf("\nSolution:\n");
        for (int i = 0; i < N; i++) {
            printf("x[%d] = %8.4f\n", i, result[i]);
        }
    }

    MPI_Finalize();
    return 0;
}