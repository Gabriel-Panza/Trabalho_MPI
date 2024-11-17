#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include "libMatriz/matriz.h"
#include "gauss_op.c"

void main(int argc, char** argv) {
    int rank, size, tam;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n; // Dimensão da matriz
    double** matriz; // Ponteiro para a matriz
    ValorIndice* coluna; // Para armazenar a coluna extraída

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
    for(int i=0; i<n; i++){
        int tamanho_coluna = n-i;
        if(rank == 0){
            coluna = extrair_coluna(&matriz[i], tamanho_coluna, n, i);
            // &matriz[i] retorna uma matriz contendo somente linhas a partir de i (n-i linhas totais)
            printf("\n\n");
            printf("Passando tamanho_coluna = %d, indice_coluna = %d\n", tamanho_coluna, i);
            printf("Coluna extraída:\n");
            imprimir_vetor_valor_indice(coluna, tamanho_coluna);
        }
        int i_pivo_coluna = achar_maior_da_coluna(coluna, tamanho_coluna, 0).indice; 
        // sobre a linha acima:
        // 1. somente a raiz não contém lixo de memória
        // 2. o vetor coluna é espalhado (scatter) dentro da função, não sendo necessário se preocupar se outros ranks tem lixo em coluna.
        // 3. o índice retornado é em relação à coluna, que tem linhas do que a matriz da segunda iteração em diante

        if(rank == 0){
            int i_pivo_matriz = i + i_pivo_coluna; // para obter um índice em relação à matriz, em vez da coluna com linhas a menos.
            printf("Maior valor encontrado: %f\n", coluna[i_pivo_coluna].valor);
            troca_linhas(matriz, n, i, i_pivo_matriz);
            printf("Matriz após a troca:\n");
            imprimir_matriz(matriz, n);
        }
        
    }
    // Libera os recursos da matriz no processo rank 0
    if (rank == 0) {
        liberar_matriz(matriz, n); // Libere a matriz
        free(coluna); // Libere o vetor coluna
    }
    // printf("Rank %d finalizado.\n", rank);
    MPI_Finalize();
}
