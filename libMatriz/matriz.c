#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

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
void imprimirMatriz(double** matriz, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) { // Considera n + 1 colunas para Ax = B
            printf("%10.2f ", matriz[i][j]); // Exibe com duas casas decimais
        }
        printf("\n");
    }
}

// Função para liberar a memória da matriz
void liberarMatriz(double** matriz, int n) {
    for (int i = 0; i < n; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

// Função para imprimir o vetor de primos
void imprimirVetor(double* vetor, int n) {
    for (int i = 0; i < n; i++) {
        printf("%.2f ", vetor[i]);
    }
    printf("\n");
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
