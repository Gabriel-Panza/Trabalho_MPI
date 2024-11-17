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