#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "libMatriz/matriz.h"

// Função para trocar linhas
void swapRows(double** matrix, int N, int row1, int row2) {
    for (int j = 0; j < N + 1; j++) {
        double temp = matrix[row1][j];
        matrix[row1][j] = matrix[row2][j];
        matrix[row2][j] = temp;
    }
}

int main(int argc, char **argv) {
    int rank, size, tam, N;
    double **matrix;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        printf("Tamanho para a matriz quadrática = ");
        scanf("%d", &N);
        matrix = criar_matriz(N);
        MPI_Bcast(&matrix[0], N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        printf("\nInitial Matrix:\n");
        imprimirMatriz(matrix, N);
    }

    // Início da medição de tempo
    double start_time = MPI_Wtime();

    for (int k = 0; k < N; k++) {
        // Processo raiz realiza o pivotamento parcial
        if (rank == 0) {
            int maxRow = k;
            for (int i = k + 1; i < N; i++) {
                if (fabs(matrix[i][k]) > fabs(matrix[maxRow][k])) {
                    maxRow = i;
                }
            }

            // Trocar linhas se necessário
            if (maxRow != k) {
                swapRows(matrix, N, k, maxRow);
            }

            // Verifica se o pivô ainda é zero
            if (fabs(matrix[k][k]) < 1e-9) {
                printf("Erro: Matriz singular, pivô na linha %d é zero.\n", k);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }

        // Broadcast do pivô e da linha correspondente
        MPI_Bcast(&matrix[k], N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Cada processo atualiza suas linhas
        for (int i = k + 1 + rank; i < N; i += size) {
            double factor = matrix[i][k] / matrix[k][k];
            for (int j = k; j < N + 1; j++) {
                matrix[i][j] -= factor * matrix[k][j];
            }
        }

        // Sincronizar a matriz entre os processos
        for (int i = k + 1; i < N; i++) {
            MPI_Bcast(&matrix[i], N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // Processo raiz realiza a substituição retroativa
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

    // Fim da medição de tempo
    double end_time = MPI_Wtime();

    // Exibe o tempo total de execução
    if (rank == 0) {
        printf("\nExecution Time: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}