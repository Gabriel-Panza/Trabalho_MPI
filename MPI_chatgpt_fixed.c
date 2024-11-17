#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include "libMatriz/matriz.h"
#include "gauss_op.c"

int main(int argc, char** argv) {
    int rank, size, n;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double** matriz; // Ponteiro para a matriz

    // Apenas o processo de rank 0 solicita o valor de n
    if (rank == 0) {
        printf("Insira o valor de n: ");
        //fflush(stdout); // Garante que o texto é exibido imediatamente
        scanf("%d", &n);
    }

    // Compartilha o valor de n com todos os processos
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Processo rank 0 cria a matriz e extrai a coluna
    if (rank == 0) {
        matriz = criar_matriz(n);
        printf("Matriz gerada:\n");
        imprimir_matriz(matriz, n);
    }
    pivoteamento_parcial(matriz, n);
    normalizar_matriz(matriz, n);
    // Libera os recursos da matriz no processo rank 0
    if (rank == 0) {
        liberar_matriz(matriz, n); // Libere a matriz
    }
    // printf("Rank %d finalizado.\n", rank);
    MPI_Finalize();
    return 0;
}
