#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include "matriz.h"

void troca_linhas(double **matriz, int n, int linha_1, int linha_2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Tamanho do bloco para cada processo
    int bloco = n / size;
    int resto = n % size; // Para lidar com tamanhos não divisíveis por size
    int local_count = bloco + (rank < resto ? 1 : 0);

    // Vetores locais para os sub-blocos das linhas a serem trocadas
    double *sub_linha_1 = malloc(local_count * sizeof(double));
    double *sub_linha_2 = malloc(local_count * sizeof(double));

    // Vetores temporários para distribuir dados no processo 0
    double *linha_1_flat = NULL;
    double *linha_2_flat = NULL;

    if (rank == 0) {
        linha_1_flat = matriz[linha_1];
        linha_2_flat = matriz[linha_2];
    }

    // Scatter para dividir as linhas em pedaços
    int *send_counts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    int offset = 0;

    for (int i = 0; i < size; i++) {
        send_counts[i] = bloco + (i < resto ? 1 : 0);
        displs[i] = offset;
        offset += send_counts[i];
    }

    MPI_Scatterv(linha_1_flat, send_counts, displs, MPI_DOUBLE, sub_linha_1, send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(linha_2_flat, send_counts, displs, MPI_DOUBLE, sub_linha_2, send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Cada processo troca os elementos localmente
    for (int i = 0; i < local_count; i++) {
        double temp = sub_linha_1[i];
        sub_linha_1[i] = sub_linha_2[i];
        sub_linha_2[i] = temp;
    }

    // Processo 0 recolhe os resultados com Reduce
    MPI_Gatherv(sub_linha_1, send_counts[rank], MPI_DOUBLE, linha_1_flat, send_counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sub_linha_2, send_counts[rank], MPI_DOUBLE, linha_2_flat, send_counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // Liberar memória local
    free(sub_linha_1);
    free(sub_linha_2);
    free(send_counts);
    free(displs);
}

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

void somar_linha(double *linha_atual, double *linha_pivo, int tamanho) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int tamanho_local = tamanho / size;
    int resto = tamanho % size;
    // Calcula os tamanhos locais para todos os processos
    int *tamanhos = (int *)malloc(size * sizeof(int));
    int *deslocamentos = (int *)malloc(size * sizeof(int));
    for (int ind = 0; ind < size; ind++) {
        tamanhos[ind] = tamanho_local + (ind < resto ? 1 : 0);
        deslocamentos[ind] = (ind == 0) ? 0 : deslocamentos[ind - 1] + tamanhos[ind - 1];
    }

    double *sublinha = (double *)malloc(tamanhos[rank] * sizeof(double));
    double *sublinha_pivo = (double *)malloc(tamanhos[rank] * sizeof(double));

    // Distribui o vetor para os processos
    MPI_Scatterv(linha_atual, tamanhos, deslocamentos, MPI_DOUBLE, sublinha, tamanhos[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(linha_pivo, tamanhos, deslocamentos, MPI_DOUBLE, sublinha_pivo, tamanhos[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Cada processo zera as linhas atribuídas ao pivô
    for (int j = 0; j < tamanhos[rank]; j++) {
        sublinha[j] -= sublinha_pivo[j];
    }

    // Coletar as linhas atualizadas nos processos
    MPI_Gatherv(sublinha, tamanhos[rank], MPI_DOUBLE, linha_atual, tamanhos, deslocamentos, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Libera memória
    free(tamanhos);
    free(deslocamentos);
    free(sublinha);
    free(sublinha_pivo);
}

void multiplicar_linha(double* linha, int tamanho, double valor) {
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    int tamanho_local = tamanho / num_procs;
    int resto = tamanho % num_procs;
    // Calcula os tamanhos locais para todos os processos
    int *tamanhos = (int *)malloc(num_procs * sizeof(int));
    int *deslocamentos = (int *)malloc(num_procs * sizeof(int));
    for (int i = 0; i < num_procs; i++) {
        tamanhos[i] = tamanho_local + (i < resto ? 1 : 0);
        deslocamentos[i] = (i == 0) ? 0 : deslocamentos[i - 1] + tamanhos[i - 1];
    }

    double *subvetor = (double *)malloc(tamanhos[rank] * sizeof(double));

    // Distribui o vetor para os processos
    MPI_Scatterv(linha, tamanhos, deslocamentos, MPI_DOUBLE,
                 subvetor, tamanhos[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Multiplica localmente
    for (int i = 0; i < tamanhos[rank]; i++) {
        subvetor[i] *= valor;
    }

    // Reúne os resultados
    MPI_Gatherv(subvetor, tamanhos[rank], MPI_DOUBLE,
                linha, tamanhos, deslocamentos, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Libera memória
    free(tamanhos);
    free(deslocamentos);
    free(subvetor);
}

void transformacao_triangular_superior(double **matriz, int n, int rank, int size, int k) {
    // Distribuir o pivô para todos os processos
    double pivo;
    double *linha, *linha_pivo;
    if(rank==0){
        pivo = matriz[k][k];
    } 
    MPI_Bcast(&pivo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Atualizar as linhas abaixo do pivô usando somar_linha
    double fator;
    for (int i = k+1; i < n; i++) {
        if (rank==0){
            fator = matriz[i][k]/pivo;
            linha = matriz[i];
            linha_pivo = matriz[k];
        }
        MPI_Bcast(&fator, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        multiplicar_linha(linha_pivo, n+1, fator);
        somar_linha(linha, linha_pivo, n+1);
        multiplicar_linha(linha_pivo, n+1, 1/fator);
    }
}

void pivoteamento_parcial(double** matriz, int n){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ValorIndice* coluna;
    for(int i=0; i<n-1; i++){
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
        int i_pivo_matriz;
        if (rank == 0) {
            i_pivo_matriz = i + i_pivo_coluna; // para obter um índice em relação à matriz, em vez da coluna com linhas a menos.
            printf("Maior valor encontrado: %f\n", coluna[i_pivo_coluna].valor);
        }   

        MPI_Bcast(&i_pivo_matriz, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // troca a linha atual com a linha do pivô
        troca_linhas(matriz, n+1, i, i_pivo_matriz);
        if(rank == 0){
            printf("Matriz após a troca:\n");
            imprimir_matriz(matriz, n);
        }
        // zero os elementos abaixo do pivô ignorando a ultima interação
        transformacao_triangular_superior(matriz, n, rank, size, i);
        if (rank == 0) {
            printf("Matriz após zerar a coluna %d:\n", i + 1);
            imprimir_matriz(matriz, n);
        }
    }
}

void scatter_linhas(double **matriz, double *linha_local, int n, int i_linha) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Scatter apenas a linha necessária
        MPI_Scatter(matriz[i_linha], n / size, MPI_DOUBLE, linha_local, n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatter(NULL, n / size, MPI_DOUBLE, linha_local, n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void normalizar_linha(double *linha_local, int n_local, double pivo) {
    for (int j = 0; j < n_local; j++) {
        linha_local[j] /= pivo;
    }
}

void gather_linha(double *linha_local, double *linha_global, int n, int i_linha) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Gather(linha_local, n / size, MPI_DOUBLE, linha_global, n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Compartilhar a linha normalizada com todos os processos
    MPI_Bcast(linha_global, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void normalizar_matriz(double **matriz, int n) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *linha_local = malloc((n / size) * sizeof(double)); // Parte da linha que cada processo manipula
    double *linha_global = malloc(n * sizeof(double));         // Linha completa para comunicação

    for (int i = 0; i < n; i++) {
        double pivo;

        // Apenas o rank 0 calcula o pivô
        if (rank == 0) {
            pivo = matriz[i][i];
            if (pivo == 0.0) {
                fprintf(stderr, "Erro: o pivô da linha %d é zero. Não é possível normalizar.\n", i);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Compartilhar o pivô com todos os processos
        MPI_Bcast(&pivo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Scatter da linha
        scatter_linhas(matriz, linha_local, n+1, i);

        // Normalizar localmente
        normalizar_linha(linha_local, n / size, pivo);

        // Gather da linha normalizada
        gather_linha(linha_local, linha_global, n, i);

        // Apenas o rank 0 atualiza a matriz
        if (rank == 0) {
            for (int j = 0; j < n; j++) {
                matriz[i][j] = linha_global[j];
            }
        }
    }

    free(linha_local);
    free(linha_global);
}

// Processo raiz realiza a substituição retroativa
void substituicao_retroativa(double** matriz, int n, int rank, int size) {
    double* result = (double*) malloc(n * sizeof(double));
    double soma_parcial;
    
    // Inicialmente, todos os processos conhecem o valor inicial do vetor result como 0
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
    }

    // Cada processo realiza a substituição retroativa parcial
    for (int i = n - 1; i >= 0; i--) {
        // Processo raiz calcula a soma parcial
        if (rank == 0) {
            soma_parcial = matriz[i][n];
            for (int j = i + 1; j < n; j++) {
                soma_parcial -= matriz[i][j] * result[j];
            }
            result[i] = soma_parcial / matriz[i][i];
        }
        // Processo raiz comunica o valor de result[i] para todos os outros processos
        MPI_Bcast(&result[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Imprimir a solução no processo raiz
    if (rank == 0) {
        printf("\nSolution:\n");
        for (int i = 0; i < n; i++) {
            printf("x[%d] = %8.4f\n", i, result[i]);
        }
    }

    free(result);
}