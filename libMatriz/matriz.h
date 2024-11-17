#ifndef MATRIZ_PRIMOS_H
#define MATRIZ_PRIMOS_H
#include <mpi.h>
#include <stdbool.h>

typedef struct valor_indice{
    double valor;
    int indice;
} ValorIndice;

ValorIndice* extrair_coluna(double** matriz, int linhas, int colunas, int indiceColuna);


// Função para criar uma matriz n x n com valores aleatórios
double** criarMatrizAleatoria(int n);

// Função auxiliar para verificar se um número é primo
bool ehPrimo(int num);

// Função para gerar um vetor de n números primos
double* gerarPrimos(int n);

// Função para multiplicar os elementos de cada linha da matriz por cada primo
void multiplicarPorPrimos(double** matriz, int n, double* primos);

// Função para imprimir a matriz
void imprimir_matriz(double** matriz, int n);

// Função para liberar a memória da matriz
void liberar_matriz(double** matriz, int n);

// Função para imprimir o vetor de primos
void imprimir_vetor_valor_indice(ValorIndice* vetor, int n);
// Função que combina todas as etapas
double** criar_matriz(int n);

void troca_linhas(double** matriz, int n, int linha1, int linha2);

ValorIndice achar_maior_local(ValorIndice* vetor, int tamanho);

void imprimir_vetor_double(double* vetor, int n);

void aplicar_modulo(ValorIndice* vetor, int tamanho);

void comparar_maximo(void* invec, void* inoutvec, int* len, MPI_Datatype* datatype);

#endif // MATRIZ_PRIMOS_H
