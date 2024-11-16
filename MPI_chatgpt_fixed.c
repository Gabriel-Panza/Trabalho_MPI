#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include "libMatriz/matriz.h"

ValorIndice achar_maior_da_coluna(ValorIndice* coluna, int tamanho_coluna, int raiz){
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    int tamanho_local = tamanho_coluna / num_procs;
    ValorIndice* vetor_local = malloc(sizeof(ValorIndice) * tamanho_local);

    // Criar o tipo derivado MPI para ValorIndice
    MPI_Datatype MPI_ValorIndice;
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(ValorIndice, valor);
    offsets[1] = offsetof(ValorIndice, indice);
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_ValorIndice);
    MPI_Type_commit(&MPI_ValorIndice);
    // espalhar os vetores locais para cada processo
    MPI_Scatter(coluna, tamanho_local, MPI_ValorIndice, vetor_local, tamanho_local, MPI_ValorIndice, raiz, MPI_COMM_WORLD);
    // print de teste de como ficou o vetor local
    printf("Vetor local do rank %d:\n", rank);
    imprimir_vetor_valor_indice(vetor_local, tamanho_local);
    // convertendo o vetor local recebido para modulo
    aplicar_modulo(vetor_local, tamanho_local);
    // achando elemento de maior módulo
    ValorIndice maior_local = achar_maior_local(vetor_local, tamanho_local);
    // Declarando operação de comparação e passando resultados para cima na árvore de computação 
    ValorIndice resultado_final;
    MPI_Op op_comparar_maximo;
    MPI_Op_create(comparar_maximo, 1, &op_comparar_maximo);
    MPI_Reduce(&maior_local, &resultado_final, 1, MPI_ValorIndice, op_comparar_maximo, raiz, MPI_COMM_WORLD);
    free(vetor_local);
    MPI_Type_free(&MPI_ValorIndice);
    MPI_Op_free(&op_comparar_maximo);

    return resultado_final;
}

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
        coluna = extrair_coluna(matriz, n, 0);
        printf("Matriz gerada:\n");
        imprimir_matriz(matriz, n);
        printf("Coluna extraída:\n");
        imprimir_vetor_valor_indice(coluna, n);
    }

    int indice_pivo = achar_maior_da_coluna(coluna, n, 0).indice; // somente a raiz não contém lixo de memória
    if(rank == 0){
        printf("Linha %d com maior pivô: \n", indice_pivo);
        imprimir_vetor_double(matriz[indice_pivo], n);
        printf("\n");
        troca_linhas(matriz, n, 0, indice_pivo);
        printf("Matriz com pivô da linha 0 maximizado: \n");
        imprimir_matriz(matriz, n);

    }
    // Libera os recursos da matriz no processo rank 0
    if (rank == 0) {
        liberar_matriz(matriz, n); // Libere a matriz
        free(coluna); // Libere o vetor coluna
    }
    // printf("Rank %d finalizado.\n", rank);
    MPI_Finalize();
}
