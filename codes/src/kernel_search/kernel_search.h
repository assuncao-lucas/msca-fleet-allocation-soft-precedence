#pragma once

#include "src/heuristic_solution.h"
#include "src/benders_generic_callback.h"
#include "src/initial_solution/initial_solution.h"
#include <vector>
#include <list>
#include <boost/dynamic_bitset.hpp>

class Instance;

class KernelSearch
{
public:
    explicit KernelSearch(Instance &instance);
    virtual ~KernelSearch();

    KSHeuristicSolution *Run(Formulation formulation, int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool feasibility_emphasis, bool multithreading);

    friend HeuristicSolution *InitalSolutionGenerator::GenerateInitialSolution(Instance &inst, bool multihreading);

private:
    IloEnv *env_ = nullptr;     // Cplex environment
    IloCplex *cplex_ = nullptr; // Cplex solver
    IloModel *model_ = nullptr; // Cplex model

    IloNumVar slack_;
    MasterVariables master_vars_;

    IloNumArray curr_x_values_;
    IloNumArray curr_y_values_;
    IloNumVarArray curr_mip_start_vars_;
    IloNumArray curr_mip_start_vals_;

    boost::dynamic_bitset<> curr_int_y_;
    boost::dynamic_bitset<> curr_int_x_;
    double curr_best_solution_value_ = -1;

    double *R0_ = nullptr;
    double *Rn_ = nullptr;

    const Instance &instance_;

    void BuildKernelAndBuckets(Formulation formulation, KSHeuristicSolution *solution, int ks_max_size_bucket, bool multithreading);
    void InitCplex(bool multithreading);
    void ResetCplex();
    void BuildModel(Formulation formulation, bool linearly_relaxed, bool disable_all_binary_vars, bool export_model);
    void RetrieveSolutionArcVertexValues();
    void BuildHeuristicSolution(KSHeuristicSolution *);
    void PrintKernelAndBuckets();
    void UpdateModelVarBounds(boost::dynamic_bitset<> &vars_entering_kernel, boost::dynamic_bitset<> &vars_leaving_kernel, boost::dynamic_bitset<> &curr_reference_kernel);

    boost::dynamic_bitset<> curr_kernel_bitset_;
    std::vector<boost::dynamic_bitset<>> buckets_bitsets_;

    bool found_int_x_ = false;
};