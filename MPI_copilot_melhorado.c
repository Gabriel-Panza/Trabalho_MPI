#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "libMatriz/matriz.h"

void gaussElimination(double **U, int n, int rank, int size) {
    // Seção 1:
    // Quando a linha que o processo é responsável chega, ele normaliza o pivô
    for (int k = 0; k < n; k++) { 
        if (k % size == rank) {
            for (int j = k + 1; j <= n; j++) { 
            // É assumido que os valores à esquerda do pivô são zero 
            // pois, exceto para a linha 0, que não tem elementos à esquerda,
            // as colunas à esquerda do pivô já estarão zeradas pela seção 2.
                U[k][j] /= U[k][k];
            }
            U[k][k] = 1.0;
        }
        // todos os processos recebem a linha anterior normalizada 
        // obs: não houve ganho de tempo por paralelismo relevante até agora, 
        // foi uma etapa de controle para que todos tenham a linha normalizada.
        MPI_Bcast(U[k], n + 1, MPI_DOUBLE, k % size, MPI_COMM_WORLD); 

        // Seção 2:
        // O for e o if abaixo controlam a(s)
        // linha(s) que cada processo é responsável
        // Assim, cada processo terá uma ou mais linhas "locais"
        // e também a linha "global" distribuída para todos
        for (int i = k + 1; i < n; i++) {
            if (i % size == rank) {
                // Subtraimos de cada linha a linha "global" 
                // multiplicada pelo valor da coluna do pivô local.
                // Isso garante que ambas as linhas terão o mesmo valor
                // na coluna do pivô, resultando em zero ao subtrair.
                for (int j = k + 1; j <= n; j++) {
                    U[i][j] -= U[i][k] * U[k][j];
                }
                U[i][k] = 0.0;
            }
        }
    }
}

void backSubstitution(double **U, double *x, int n, int rank, int size) {
    // De baixo para cima, resolvemos o valor da incógnita da diagonal principal
    // e transmitimos o resultado para as linhas acima, o que permite as linhas
    // acima a resolverem suas próprias incógnitas.
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

    gaussElimination(U, n, rank, size); // retorna triangular superior
    backSubstitution(U, x, n, rank, size); // resolve incógnitas de uma triangular superior


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