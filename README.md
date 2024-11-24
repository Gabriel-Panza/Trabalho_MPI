# Trabalho_MPI
## Codigo ChatGPT-4
Para compilar o código com suporte a MPI:
```bash
mpicc MPI_copilot.c libMatriz/matriz.c -o MPI_copilot -lm
```

Para executar em múltiplos processos (substitua x pelo número de processos):
```bash
mpirun -np x ./MPI_copilot
```

## Código do Grupo
Para compilar o código com suporte a MPI:
```bash
mpicc MPI_copilot_melhorado.c libMatriz/matriz.c -o MPI_copilot_melhorado -lm
```

Para executar em múltiplos processos (substitua x pelo número de processos):
```bash
mpirun -np x ./MPI_copilot_melhorado
```