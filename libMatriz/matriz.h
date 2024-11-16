#ifndef MATRIZ_PRIMOS_H
#define MATRIZ_PRIMOS_H

#include <stdbool.h>

// Função para criar uma matriz n x n com valores aleatórios
int** criarMatrizAleatoria(int n);

// Função auxiliar para verificar se um número é primo
bool ehPrimo(int num);

// Função para gerar um vetor de n números primos
int* gerarPrimos(int n);

// Função para multiplicar os elementos de cada linha da matriz por cada primo
void multiplicarPorPrimos(int** matriz, int n, int* primos);

// Função para imprimir a matriz
void imprimirMatriz(double** matriz, int n);

// Função para liberar a memória da matriz
void liberarMatriz(int** matriz, int n);

// Função para imprimir o vetor de primos
void imprimirVetor(int* vetor, int n);

// Função que combina todas as etapas
double** criar_matriz(int n);

#endif // MATRIZ_PRIMOS_H
