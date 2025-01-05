CPP = g++ -Wall -fopenmp
SRCS = main.cpp fluid_solver.cpp EventManager.cpp
SEQ_SRCS = main.cpp fluid_solver.cpp EventManager.cpp 

all: phase2 phase2_seq

phase2:
	$(CPP) -Ofast $(SRCS) -o fluid_sim

phase2_seq:
	$(CPP) -Ofast $(SEQ_SRCS) -o fluid_sim_seq

clean:
	@echo Cleaning up...
	@rm fluid_sim fluid_sim_seq
	@echo Done.

runseq: phase2_seq
	./fluid_sim_seq

runpar: phase2
	OMP_NUM_THREADS=20 ./fluid_sim
