#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include "libMatriz/matriz.h"

ValorIndice achar_maior_da_coluna(ValorIndice* coluna, int tamanho_coluna, int raiz) {
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Cálculo do tamanho local para cada processo
    int tamanho_local = tamanho_coluna / num_procs;
    int resto = tamanho_coluna % num_procs;

    // Ajustar o tamanho local do último processo para lidar com o resto
    if (rank == num_procs - 1) {
        tamanho_local += resto;
    }

    // Alocar espaço para o vetor local
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

    // Vetor que contém o tamanho local de cada processo
    int* tamanhos_locais = malloc(num_procs * sizeof(int));
    int* deslocamentos = malloc(num_procs * sizeof(int));

    // Calcular tamanhos e deslocamentos para Scatterv
    for (int i = 0; i < num_procs; i++) {
        tamanhos_locais[i] = tamanho_coluna / num_procs;
        if (i == num_procs - 1) {
            tamanhos_locais[i] += resto;
        }
        deslocamentos[i] = (i > 0) ? deslocamentos[i - 1] + tamanhos_locais[i - 1] : 0;
    }

    // Usar Scatterv para distribuir partes desiguais do vetor
    MPI_Scatterv(
        coluna, tamanhos_locais, deslocamentos, MPI_ValorIndice,
        vetor_local, tamanho_local, MPI_ValorIndice, 
        raiz, MPI_COMM_WORLD
    );

    // Aplicar módulo ao vetor local
    aplicar_modulo(vetor_local, tamanho_local);

    // Achar o maior elemento local
    ValorIndice maior_local = achar_maior_local(vetor_local, tamanho_local);

    // Declarar operação de comparação e reduzir para o maior elemento global
    ValorIndice resultado_final;
    MPI_Op op_comparar_maximo;
    MPI_Op_create(comparar_maximo, 1, &op_comparar_maximo);
    MPI_Reduce(&maior_local, &resultado_final, 1, MPI_ValorIndice, op_comparar_maximo, raiz, MPI_COMM_WORLD);

    // Liberar recursos
    free(vetor_local);
    free(tamanhos_locais);
    free(deslocamentos);
    MPI_Type_free(&MPI_ValorIndice);
    MPI_Op_free(&op_comparar_maximo);

    return resultado_final;
}

void pivoteamento_parcial(double** matriz, int n){
    int rank;
    ValorIndice* coluna;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
}