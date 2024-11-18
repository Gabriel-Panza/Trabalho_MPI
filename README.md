# Trabalho_MPI
## Codigo ChatGPT-4
Para compilar o código com suporte a MPI:
```bash
mpicc MPI_chatgpt_fixed.c libMatriz/matriz.c -o mpi_chatgpt_fixed -lm
```

Para executar em múltiplos processos (substitua x pelo número de processos):
```bash
mpirun -np x ./mpi_chatgpt_fixed
```

## Código do Grupo
Para compilar o código com suporte a MPI:
```bash
mpicc MPI_gauss.c libMatriz/matriz.c -o MPI_gauss -lm
```

Para executar em múltiplos processos (substitua x pelo número de processos):
```bash
mpirun -np x ./MPI_gauss
```