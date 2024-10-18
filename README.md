# Original
Versão original do codigo com flags na MakeFile

# versao1

Versão otimizada do chatgpt (estudar)

# versao2

Troca de Ciclos

# versao3

Merge de versao1 e 2

# versao4

(introduzir o que se fez)

# versao5

(introduzir o que se fez)

# Comandos importantes
correr perf, se for preciso introduzir mais flags
```
perf stat -r 5 -e instructions,cycles ./fluid_sim 
```
criar o ficheiro para criacao de graficos
```
gprof fluid_sim gmon.out > *versaon*
```
criar grafico
```
python3 gprof2dot.py <*versaon* | dot -Tpng -o *versaon*.png
```
criar mais versoes
```
mkdir versaon
cp original/. versaon
```

