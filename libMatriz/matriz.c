#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "matriz.h"
#include <mpi.h>
#include <stddef.h>
#include <math.h>


// Função para criar uma matriz n x (n+1) com valores aleatórios (double)
double** criarMatrizAleatoria(int n) {
    // Aloca a matriz com n linhas
    double** matriz = (double**)malloc(n * sizeof(double*));

    for (int i = 0; i < n; i++) {
        // Aloca n + 1 colunas para cada linha (n para A e 1 para B)
        matriz[i] = (double*)malloc((n + 1) * sizeof(double));

        for (int j = 0; j < n + 1; j++) {
            // Preenche os elementos da matriz com valores aleatórios entre 0.0 e 99.99
            matriz[i][j] = (double)(rand() % 100) + ((double)rand() / RAND_MAX);
        }
    }

    return matriz;
}

// Função auxiliar para verificar se um número é primo
bool ehPrimo(int num) {
    if (num < 2) return false;
    for (int i = 2; i * i <= num; i++) {
        if (num % i == 0) return false;
    }
    return true;
}

// Função para gerar um vetor de n números primos
double* gerarPrimos(int n) {
    double* primos = (double*)malloc(n * sizeof(double));
    int contador = 0;
    int candidato = 2;

    while (contador < n) {
        if (ehPrimo(candidato)) {
            primos[contador++] = (double)candidato; // Armazena os primos como double
        }
        candidato++;
    }
    return primos;
}

// Função para multiplicar os elementos de cada linha da matriz pelos primos
void multiplicarPorPrimos(double** matriz, int n, double* primos) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matriz[i][j] *= primos[i]; // Multiplica os elementos da linha i pelo primo correspondente
        }
    }
}

// Função para imprimir a matriz
void imprimir_matriz(double** matriz, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) { // Considera n + 1 colunas para Ax = B
            printf("%10.2f ", matriz[i][j]); // Exibe com duas casas decimais
        }
        printf("\n");
    }
}

// Função para liberar a memória da matriz
void liberar_matriz(double** matriz, int n) {
    for (int i = 0; i < n; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

// Função para imprimir o vetor de primos
void imprimir_vetor_valor_indice(ValorIndice* vetor, int n) {
    printf("[");
    for (int i = 0; i < n; i++) {
        printf("(%.2f, %d)", vetor[i].valor, vetor[i].indice);
        if (i < n - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

// Função para criar a matriz e realizar as operações
double** criar_matriz(int n) {
    // Cria a matriz aleatória
    double** matriz = criarMatrizAleatoria(n);

    // Gera um vetor de n números primos
    double* primos = gerarPrimos(n);

    // Multiplica os elementos da matriz pelos primos
    multiplicarPorPrimos(matriz, n, primos);

    // Libera o vetor de primos
    free(primos);

    // Retorna a matriz resultante
    return matriz;
}

void troca_linhas(double** matriz, int n, int linha1, int linha2) {
    for (int j = 0; j < n + 1; j++) {
        double temp = matriz[linha1][j];
        matriz[linha1][j] = matriz[linha2][j];
        matriz[linha2][j] = temp;
    }
}

void imprimir_vetor_double(double* vetor, int n){
    printf("[");
    for(int i=0; i < n-1; i++){
        printf("%f, ", vetor[i]);
    }
    printf("%f]", vetor[n-1]);
}

ValorIndice achar_maior_local(ValorIndice* vetor, int tamanho) {
    ValorIndice maior = vetor[0];
    for (int i = 1; i < tamanho; i++) {
        if (vetor[i].valor > maior.valor) {
            maior = vetor[i];
        }
    }

    return maior;
}

void aplicar_modulo(ValorIndice* vetor, int tamanho) {
    for (int i = 0; i < tamanho; i++) {
        vetor[i].valor = fabs(vetor[i].valor);
    }
}

void comparar_maximo(void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) {
    ValorIndice* in = (ValorIndice*)invec;
    ValorIndice* inout = (ValorIndice*)inoutvec;

    for (int i = 0; i < *len; i++) {
        if (in[i].valor > inout[i].valor) {
            inout[i] = in[i];
        }
    }
}

ValorIndice* extrair_coluna(double** matriz, int linhas, int colunas, int indiceColuna) {
    // Validar se o índice da coluna é válido
    if (indiceColuna < 0 || indiceColuna >= colunas) {
        printf("Erro: Índice de coluna fora do intervalo: indice = %d, colunas = %d\n", indiceColuna, colunas);
        return NULL;
    }

    // Alocar memória para o array de ValorIndice
    ValorIndice* resultado = (ValorIndice*)malloc(linhas * sizeof(ValorIndice));
    if (resultado == NULL) {
        printf("Erro ao alocar memória.\n");
        return NULL;
    }

    // Preencher o array com os valores e índices da coluna
    for (int i = 0; i < linhas; i++) {
        resultado[i].valor = matriz[i][indiceColuna];
        resultado[i].indice = i;
    }

    // Retornar o array com os valores extraídos
    return resultado;
}
