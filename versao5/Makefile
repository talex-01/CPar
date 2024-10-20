CPP = g++ -Wall 
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all:
	$(CPP) -pg -march=native -Ofast -funroll-loops -ftree-vectorize $(SRCS) -o fluid_sim

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
