#pragma once

#include <vector>
#include <unordered_map>
#include <ilcplex/ilocplex.h>
#include <boost/functional/hash/hash.hpp>
#include "src/instance.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/matrix.hpp"
#include "src/general.h"
#include "src/heuristic_solution.h"

class VarToIndexMap
{
public:
  VarToIndexMap() = default;
  ~VarToIndexMap() = default;

  int varToIndex(int i, int j)
  {
    return var_to_index_map_[{i, j}];
  }

  std::pair<int, int> indexToVar(int index)
  {
    return index_to_var_map_[index];
  }

  size_t numVars() const
  {
    return index_to_var_map_.size();
  }

  void addEntry(int i, int j)
  {
    size_t counter = numVars();
    var_to_index_map_[{i, j}] = counter;
    index_to_var_map_[counter] = {i, j};
  }

  struct pair_hash
  {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const
    {
      // auto h1 = std::hash<T1>{}(p.first);
      // auto h2 = std::hash<T2>{}(p.second);
      // return h1 ^ h2;
      std::size_t seed = 0;
      boost::hash_combine(seed, p.first);
      boost::hash_combine(seed, p.second);

      return seed;
    }
  };

  const std::unordered_map<std::pair<int, int>, int, pair_hash> &var_to_index_map() const
  {
    return var_to_index_map_;
  }
  const std::unordered_map<int, std::pair<int, int>> &index_to_var_map() const
  {
    return index_to_var_map_;
  }

private:
  std::unordered_map<std::pair<int, int>, int, pair_hash> var_to_index_map_;
  std::unordered_map<int, std::pair<int, int>> index_to_var_map_;
};

struct VehicleSequencingModelVariables
{
  IloNumVarArray x;
  IloNumVarArray z;
  IloNumVarArray w;
  IloNumVarArray U;

  VarToIndexMap x_var_to_index;
  VarToIndexMap z_var_to_index;
  VarToIndexMap w_var_to_index;
  VarToIndexMap U_var_to_index;
};

int a_var_to_index(int vertex, int budget, int num_vertices);
int f_var_to_index(int arc_pos, int budget, int num_arcs);

std::pair<int, int> index_to_a_var(int index, int num_vertices);

template <class T>
void SetSolutionStatus(IloCplex &cplex, Solution<T> &solution, bool solve_relax)
{
  // std::cout << cplex.getCplexStatus() << std::endl;
  if (solve_relax)
  {
    if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
      solution.lp_ = -1;
    else
      solution.lp_ = cplex.getObjValue();
    // if((cplex.getCplexStatus() == IloCplex::Optimal)||(cplex.getCplexStatus() == IloCplex::OptimalTol))
    // {
    //   solution.is_optimal_ =  true;
    // }
  }
  else
  {
    if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
      solution.is_feasible_ = false;
    else
    {
      double cost = cplex.getObjValue();
      if ((cplex.getCplexStatus() == IloCplex::Optimal) || (cplex.getCplexStatus() == IloCplex::OptimalTol))
      {
        solution.lb_ = cost;
        solution.ub_ = cost;
        solution.is_optimal_ = true;
      }
      else
      {
        solution.lb_ = cost;
        solution.ub_ = cplex.getBestObjValue();
        if (!(double_less(solution.lb_, solution.ub_)))
          solution.is_optimal_ = true;
      }
    }
  }
  solution.num_nodes_ = cplex.getNnodes();
}

void allocateVehicleSequencingModelVariables(IloEnv &env, VehicleSequencingModelVariables &vars, const Instance &instance, bool solve_relax, bool disable_all_binary_vars = false);

static void populateByRowVehicleSequencingModel(IloCplex &cplex, IloEnv &env, IloModel &model, VehicleSequencingModelVariables &vars, const Instance &instance, bool export_model);

void vehicleSequencingModel(Instance &inst, double time_limit, bool solve_relax, bool export_model);

void optimize(IloCplex &cplex, IloEnv &env, IloModel &model, std::optional<Formulation> formulation, VehicleSequencingModelVariables &vars, Instance &instance, double total_time_limit, bool solve_relax, bool apply_benders, bool apply_benders_generic_callback, bool combine_feas_op_cuts, bool separate_benders_cuts_relaxation, bool use_valid_inequalities, bool find_root_cuts, double *R0, double *Rn, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &solution);
void optimizeLP(IloCplex &cplex, IloEnv &env, IloModel &model, Instance &instance, double total_time_limit, double *R0, double *Rn, Solution<double> &);