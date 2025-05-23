SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio2211/concert
CPLEXLIBDIR   = $(CPLEXDIR)/cplex/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CLNFLAGS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -lilocplex -lconcert -lcplex -m64 -lm -lpthread

CC_DEBUG = -DDEBUG -g
CC_RELEASE = -DNDEBUG
CC_VALGRIND = -DNDEBUG -g -O0

COPT  = -m64 -O2 -fPIC -fexceptions $(CC_RELEASE) -DIL_STD -DLONG_MAX=0x7FFFFFFFL
GENERALINCDIR   = -I ./codes
CPLEXINCDIR   = -I $(CPLEXDIR)/cplex/include -I $(CPLEXDIR)/concert/include
CFLAGS = $(COPT) $(GENERALINCDIR) -std=c++17
CFLAGS2  = $(COPT) $(CPLEXINCDIR)

CC=ccache g++ -std=c++17

PROG_DIR=codes/src
PROG_BIN=codes/bin

MAIN_SRC=$(PROG_DIR)/main.cpp

ARC_SRC=$(PROG_DIR)/arc.cpp
ARC_H=$(PROG_DIR)/arc.h
ARC_OBJ=$(PROG_BIN)/arc.o

GRAPH_SRC=$(PROG_DIR)/graph.cpp
GRAPH_H=$(PROG_DIR)/graph.h
GRAPH_OBJ=$(PROG_BIN)/graph.o

GRAPH_ALGO_SRC=$(PROG_DIR)/graph_algorithms.cpp
GRAPH_ALGHO_H=$(PROG_DIR)/graph_alhorithms.h
GRAPH_ALGO_OBJ=$(PROG_BIN)/graph_algorithms.o

GENERAL_SRC=$(PROG_DIR)/general.cpp
GENERAL_H=$(PROG_DIR)/general.h
GENERAL_OBJ=$(PROG_BIN)/general.o

USERCUT_SRC=$(PROG_DIR)/user_cut.cpp
USERCUT_H=$(PROG_DIR)/user_cut.h
USERCUT_OBJ=$(PROG_BIN)/user_cut.o

TIMER_SRC=$(PROG_DIR)/timer.cpp
TIMER_H=$(PROG_DIR)/timer.h
TIMER_OBJ=$(PROG_BIN)/timer.o

FEASIBILITY_PUMP_SRC=$(PROG_DIR)/feasibility_pump/feasibility_pump.cpp
FEASIBILITY_PUMP_H=$(PROG_DIR)/feasibility_pump/feasibility_pump.h
FEASIBILITY_PUMP_OBJ=$(PROG_BIN)/feasibility_pump.o

KERNEL_SEARCH_SRC=$(PROG_DIR)/kernel_search/kernel_search.cpp
KERNEL_SEARCH_H=$(PROG_DIR)/kernel_search/kernel_search.h
KERNEL_SEARCH_OBJ=$(PROG_BIN)/kernel_search.o

LOCAL_BRANCHING_SRC=$(PROG_DIR)/local_branching.cpp
LOCAL_BRANCHING_H=$(PROG_DIR)/local_branching.h
LOCAL_BRANCHING_OBJ=$(PROG_BIN)/local_branching.o

ALNS_SRC=$(PROG_DIR)/ALNS/ALNS.cpp
ALNS_H=$(PROG_DIR)/ALNS/ALNS.h
ALNS_OBJ=$(PROG_BIN)/ALNS.o

SIMULATED_ANNEALING_SRC=$(PROG_DIR)/simulated_annealing/simulated_annealing.cpp
SIMULATED_ANNEALING_H=$(PROG_DIR)/simulated_annealing/simulated_annealing.h
SIMULATED_ANNEALING_OBJ=$(PROG_BIN)/simulated_annealing.o

INITIAL_SOL_SRC=$(PROG_DIR)/initial_solution/initial_solution.cpp
INITIAL_SOL_H=$(PROG_DIR)/initial_solution/initial_solution.h
INITIAL_SOL_OBJ=$(PROG_BIN)/initial_solution.o

LOCAL_SEARCHES_SRC=$(PROG_DIR)/local_searches/local_searches.cpp
LOCAL_SEARCHES_H=$(PROG_DIR)/local_searches/local_searches.h
LOCAL_SEARCHES_OBJ=$(PROG_BIN)/local_searches.o

HEURISTIC_SOLUTION_SRC=$(PROG_DIR)/heuristic_solution.cpp
HEURISTIC_SOLUTION_H=$(PROG_DIR)/heuristic_solution.h
HEURISTIC_SOLUTION_OBJ=$(PROG_BIN)/heuristic_solution.o

ROUTE_SRC=$(PROG_DIR)/route.cpp
ROUTE_H=$(PROG_DIR)/route.h
ROUTE_OBJ=$(PROG_BIN)/route.o

SOLUTION_HPP=$(PROG_DIR)/solution.hpp
SOLUTION_OBJ=$(PROG_BIN)/solution.o

MATRIX_HPP=$(PROG_DIR)/matrix.hpp
MATRIX_OBJ=$(PROG_BIN)/matrix.o

INSTANCE_SRC=$(PROG_DIR)/instance.cpp
INSTANCE_H=$(PROG_DIR)/instance.h
INSTANCE_OBJ=$(PROG_BIN)/instance.o

FORM_SRC=$(PROG_DIR)/exact/formulations.cpp
FORM_H=$(PROG_DIR)/exact/formulations.h
FORM_OBJ=$(PROG_BIN)/formulations.o

BENDERS_CALLBACK_SRC=$(PROG_DIR)/benders_generic_callback.cpp
BENDERS_CALLBACK_H=$(PROG_DIR)/benders_generic_callback.h
BENDERS_CALLBACK_OBJ=$(PROG_BIN)/benders_generic_callback.o

arc: $(ARC_SRC) $(ARC_H)
	$(CC) $(CFLAGS) -c $(ARC_SRC) -o $(ARC_OBJ)

graph: $(GRAPH_SRC) $(GRAPH_H)
	$(CC) $(CFLAGS) -c $(GRAPH_SRC) -o $(GRAPH_OBJ)

graph_algo: $(GRAPH_ALGO_SRC) $(GRAPH_ALGO_H)
	$(CC) $(CFLAGS) -c $(GRAPH_ALGO_SRC) -o $(GRAPH_ALGO_OBJ)

general: $(GENERAL_SRC) $(GENERAL_H)
	$(CC) $(CFLAGS) -c $(GENERAL_SRC) -o $(GENERAL_OBJ)

USERCUT_OBJ: $(USERCUT_SRC) $(USERCUT_H)
	$(CC) $(CFLAGS) -c $(USERCUT_SRC) -o $(USERCUT_OBJ)

timer: $(TIMER_SRC) $(TIMER_H)
	$(CC) $(CFLAGS) -c $(TIMER_SRC) -o $(TIMER_OBJ)

feasibility_pump: $(FEASIBILITY_PUMP_SRC) $(FEASIBILITY_PUMP_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(FEASIBILITY_PUMP_SRC) -o $(FEASIBILITY_PUMP_OBJ)

kernel_search: $(KERNEL_SEARCH_SRC) $(KERNEL_SEARCH_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(KERNEL_SEARCH_SRC) -o $(KERNEL_SEARCH_OBJ)

LOCAL_BRANCHING_OBJ: $(LOCAL_BRANCHING_SRC) $(LOCAL_BRANCHING_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(LOCAL_BRANCHING_SRC) -o $(LOCAL_BRANCHING_OBJ)

alns: $(ALNS_SRC) $(ALNS_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(ALNS_SRC) -o $(ALNS_OBJ)

simulated_annealing: $(SIMULATED_ANNEALING_SRC) $(SIMULATED_ANNEALING_H)
	$(CC) $(CFLAGS) -c $(SIMULATED_ANNEALING_SRC) -o $(SIMULATED_ANNEALING_OBJ)

initial_solution: $(INITIAL_SOL_SRC) $(INITIAL_SOL_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(INITIAL_SOL_SRC) -o $(INITIAL_SOL_OBJ)

local_searches: $(LOCAL_SEARCHES_SRC) $(LOCAL_SEARCHES_H)
	$(CC) $(CFLAGS) -c $(LOCAL_SEARCHES_SRC) -o $(LOCAL_SEARCHES_OBJ)

heuristic_solution: $(HEURISTIC_SOLUTION_SRC) $(HEURISTIC_SOLUTION_H)
	$(CC) $(CFLAGS) -c $(HEURISTIC_SOLUTION_SRC) -o $(HEURISTIC_SOLUTION_OBJ)

route: $(ROUTE_SRC) $(ROUTE_H)
	$(CC) $(CFLAGS) -c $(ROUTE_SRC) -o $(ROUTE_OBJ)

instance: $(INSTANCE_SRC) $(INSTANCE_H)
	$(CC) $(CFLAGS) -c $(INSTANCE_SRC) -o $(INSTANCE_OBJ)

solution: $(SOLUTION_HPP)
	$(CC) $(CFLAGS) -c $(SOLUTION_HPP) -o $(SOLUTION_OBJ)

matrix: $(MATRIX_HPP)
	$(CC) $(CFLAGS) -c $(MATRIX_HPP) -o $(MATRIX_OBJ)

formulations: $(FORM_SRC) $(FORM_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(FORM_SRC) -o $(FORM_OBJ)

benders_callback: $(BENDERS_CALLBACK_SRC) $(BENDERS_CALLBACK_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(BENDERS_CALLBACK_SRC) -o $(BENDERS_CALLBACK_OBJ)

stop: arc graph graph_algo general instance route matrix timer heuristic_solution solution benders_callback formulations feasibility_pump local_searches alns simulated_annealing initial_solution kernel_search
	$(CC) $(CFLAGS) $(CFLAGS2) $(ARC_SRC) $(GRAPH_SRC) $(GRAPH_ALGO_SRC) $(GENERAL_SRC) $(USERCUT_SRC) $(ROUTE_SRC) $(MATRIX_HPP) $(INSTANCE_SRC) $(TIMER_SRC) $(HEURISTIC_SOLUTION_SRC) $(SOLUTION_HPP) $(BENDERS_CALLBACK_SRC) $(FORM_SRC) $(FEASIBILITY_PUMP_SRC) $(INITIAL_SOL_SRC) $(LOCAL_SEARCHES_SRC) $(ALNS_SRC) $(SIMULATED_ANNEALING_SRC) $(KERNEL_SEARCH_SRC) $(MAIN_SRC) -o $(PROG_BIN)/stop $(CLNFLAGS)

new: matrix instance timer general formulations solution
	$(CC) $(CFLAGS) $(CFLAGS2) $(MATRIX_HPP) $(INSTANCE_SRC) $(TIMER_SRC) $(GENERAL_SRC) $(FORM_SRC) $(SOLUTION_HPP) $(MAIN_SRC) -o $(PROG_BIN)/stop $(CLNFLAGS)

clean:
	rm $(PROG_BIN)/*
