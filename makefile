final : mcmc.py rutina3
	python mcmc.py

rutina3 : perfilrho.py rutina2
	python perfilrho.py

rutina2 : rutina1
	time ./exe.out 450 0.1

rutina1 : main.c evolve.c inicial.c
	gcc main.c evolve.c inicial.c -lm -fopenmp -o exe.out
