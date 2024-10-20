CPP = g++ -Wall 
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all:
	$(CPP) -march=native -Ofast -funroll-loops -ffast-math -flto -fstrict-aliasing -ftree-vectorize -floop-block -fprefetch-loop-arrays $(SRCS) -o fluid_sim

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
