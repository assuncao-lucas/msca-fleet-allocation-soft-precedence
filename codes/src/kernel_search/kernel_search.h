#pragma once

#include "src/heuristic_solution.h"
#include "src/exact/formulations.h"
#include <vector>
#include <list>
#include <boost/dynamic_bitset.hpp>

class Instance;

class KernelSearch
{
public:
    explicit KernelSearch(Instance &instance);
    virtual ~KernelSearch();

    KSHeuristicSolution *Run(Model &formulation, std::list<UserCut *> *initial_cuts, int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool feasibility_emphasis, bool multithreading);

private:
    IloNumArray curr_x_values_;
    IloNumArray curr_y_values_;
    IloNumArray curr_mip_start_vals_;

    boost::dynamic_bitset<> curr_solution_items_;

    double curr_best_solution_value_ = -1;

    const Instance &instance_;

    void BuildKernelAndBuckets(Model &formulation, std::list<UserCut *> *initial_cuts, KSHeuristicSolution *solution, std::optional<std::unordered_set<int>> initial_kernel_items, int ks_max_size_bucket, bool multithreading);
    void BuildHeuristicSolution(KSHeuristicSolution *);
    void PrintKernelAndBuckets();
    void UpdateModelVarBounds(boost::dynamic_bitset<> &vars_entering_kernel, boost::dynamic_bitset<> &vars_leaving_kernel, boost::dynamic_bitset<> &curr_reference_kernel);

    boost::dynamic_bitset<> curr_kernel_bitset_;
    std::vector<boost::dynamic_bitset<>> buckets_bitsets_;

    bool found_feasible_solution = false;

    static std::pair<std::vector<std::pair<int, std::vector<int>>>, std::unordered_set<int>> findInitialSolution(const Instance &instance);
};