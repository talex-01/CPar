CPP = g++ -Wall 
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all:
	$(CPP) -pg -funroll-loops -march=native -Ofast -ftree-vectorize -msse4 $(SRCS) -o fluid_sim

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
