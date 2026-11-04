GUROBI_HOME ?= /opt/gurobi1301/linux64
CXXFLAGS = -O2 -std=c++17 -I../knapsack -I$(GUROBI_HOME)/include
LIBS = $(GUROBI_HOME)/lib/libgurobi_c++.a -L$(GUROBI_HOME)/lib -lgurobi130 -lpthread -lm

nested_benders: nested_benders.cpp ../knapsack/GurobiSolve3Stage.cpp
	g++ $(CXXFLAGS) -o $@ nested_benders.cpp ../knapsack/GurobiSolve3Stage.cpp $(LIBS)
