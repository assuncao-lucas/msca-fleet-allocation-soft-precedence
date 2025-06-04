#pragma once

#include <vector>
#include <unordered_map>
#include <string>
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
    auto element = var_to_index_map_.find({i, j});
    if (element == var_to_index_map_.end())
      throw "tried to access var not mapped: " + std::to_string(i) + " " + std::to_string(j);
    else
      return element->second;
  }

  std::pair<int, int> indexToVar(int index)
  {
    auto element = index_to_var_map_.find(index);
    if (element == index_to_var_map_.end())
      throw "tried to access var not mapped: " + std::to_string(index);
    else
      return element->second;
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
    // std::cout << i << " " << j << " -> " << counter << std::endl;
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
    friend std::ostream &operator<<(std::ostream &out, const VarToIndexMap &map);
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

class VehicleSequencingModelVariables
{
public:
  VehicleSequencingModelVariables() = default;
  ~VehicleSequencingModelVariables() = default;

  friend std::ostream &operator<<(std::ostream &out, const VehicleSequencingModelVariables &vars);

  IloNumVarArray x_;
  IloNumVarArray z_;
  IloNumVarArray w_;
  IloNumVarArray U_;
  IloNumVarArray X_;

  IloNumVar &x(int i, int j)
  {
    return x_[x_index(i, j)];
  }

  IloNumVar &z(int i, int h)
  {
    return z_[z_index(i, h)];
  }

  IloNumVar &w(int j, int k)
  {
    return w_[w_index(j, k)];
  }

  IloNumVar &U(int j)
  {
    return U_[U_index(j)];
  }

  IloNumVar &X(int i, int q)
  {
    return X_[X_index(i, q)];
  }

  int x_index(int i, int j)
  {
    return x_var_to_index_.varToIndex(i, j);
  }

  int z_index(int i, int h)
  {
    return z_var_to_index_.varToIndex(i, h);
  }

  int w_index(int j, int k)
  {
    return w_var_to_index_.varToIndex(j, k);
  }

  int U_index(int j)
  {
    return U_var_to_index_.varToIndex(0, j);
  }

  int X_index(int i, int q)
  {
    return X_var_to_index_.varToIndex(i, q);
  }

  VarToIndexMap x_var_to_index_;
  VarToIndexMap z_var_to_index_;
  VarToIndexMap w_var_to_index_;
  VarToIndexMap U_var_to_index_;
  VarToIndexMap X_var_to_index_;
};

class ItemSequencingModelVariables
{
public:
  ItemSequencingModelVariables() = default;
  ~ItemSequencingModelVariables() = default;

  friend std::ostream &operator<<(std::ostream &out, const ItemSequencingModelVariables &vars);

  IloNumVarArray x_;
  IloNumVarArray u_;
  IloNumVarArray U_;
  IloNumVarArray X_;

  IloNumVar &x(int i, int j)
  {
    return x_[x_index(i, j)];
  }

  IloNumVar &u(int j, int k)
  {
    return u_[u_index(j, k)];
  }

  IloNumVar &U(int j)
  {
    return U_[U_index(j)];
  }

  IloNumVar &X(int i, int q)
  {
    return X_[X_index(i, q)];
  }

  int x_index(int i, int j)
  {
    return x_var_to_index_.varToIndex(i, j);
  }

  int u_index(int j, int k)
  {
    return u_var_to_index_.varToIndex(j, k);
  }

  int U_index(int j)
  {
    return U_var_to_index_.varToIndex(0, j);
  }

  int X_index(int i, int q)
  {
    return X_var_to_index_.varToIndex(i, q);
  }

  VarToIndexMap x_var_to_index_;
  VarToIndexMap u_var_to_index_;
  VarToIndexMap U_var_to_index_;
  VarToIndexMap X_var_to_index_;
};

class VehicleSlotsModelVariables
{
public:
  VehicleSlotsModelVariables() = default;
  ~VehicleSlotsModelVariables() = default;

  friend std::ostream &operator<<(std::ostream &out, const VehicleSlotsModelVariables &vars);

  IloNumVarArray y1_;
  IloNumVarArray y2_;
  IloNumVarArray U_;
  IloNumVarArray Y_;

  IloNumVar &y1(int i, int j)
  {
    return y1_[y1_index(i, j)];
  }

  IloNumVar &y2(int i, int h)
  {
    return y2_[y2_index(i, h)];
  }

  IloNumVar &U(int j)
  {
    return U_[U_index(j)];
  }

  IloNumVar &Y(int i, int q)
  {
    return Y_[Y_index(i, q)];
  }

  int y1_index(int i, int j)
  {
    return y1_var_to_index_.varToIndex(i, j);
  }

  int y2_index(int i, int h)
  {
    return y2_var_to_index_.varToIndex(i, h);
  }

  int U_index(int j)
  {
    return U_var_to_index_.varToIndex(0, j);
  }

  int Y_index(int i, int q)
  {
    return Y_var_to_index_.varToIndex(i, q);
  }

  VarToIndexMap y1_var_to_index_;
  VarToIndexMap y2_var_to_index_;
  VarToIndexMap U_var_to_index_;
  VarToIndexMap Y_var_to_index_;
};

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

class Model // abstract class for the models.
{
public:
  Model() = default;
  ~Model();
  virtual void solve(const Instance &instance, double time_limit) = 0;

protected:
  virtual void allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars = false) = 0;
  virtual void populateByRow(const Instance &instance, bool reformulate, bool symmetry_breaking, bool export_model) = 0;
  void optimize(const Instance &instance, double total_time_limit, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCut *> *initial_cuts, std::list<UserCut *> *root_cuts, Solution<double> &solution);
  virtual void addInitialCuts(std::list<UserCut *> *initial_cuts, std::list<UserCut *> *root_cuts, Solution<double> &solution);
  virtual void addCut(UserCut *curr_cut) = 0;
  virtual bool findAndAddValidInqualities(const Instance &instance, Solution<double> &sol, std::list<UserCut *> *root_cuts);

  IloEnv *env_ = nullptr;     // Cplex environment.
  IloCplex *cplex_ = nullptr; // Cplex solver.
  IloModel *model_ = nullptr; // Cplex model.
  bool is_relaxed_ = false;
  bool reformulated_ = false;
};

class VehicleSequencingModel : public Model
{
public:
  explicit VehicleSequencingModel(Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed, bool export_model);
  ~VehicleSequencingModel() = default;
  virtual void solve(const Instance &instance, double time_limit);

private:
  virtual void allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars = false);
  virtual void populateByRow(const Instance &instance, bool reformulate, bool symmetry_breaking, bool export_model);
  virtual void addCut(UserCut *curr_cut);

  VehicleSequencingModelVariables vars_;
};

class ItemSequencingModel : public Model
{
public:
  explicit ItemSequencingModel(Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed, bool export_model);
  ~ItemSequencingModel() = default;
  virtual void solve(const Instance &instance, double time_limit);

private:
  virtual void allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars = false);
  virtual void populateByRow(const Instance &instance, bool reformulate, bool symmetry_breaking, bool export_model);
  virtual void addCut(UserCut *curr_cut);

  ItemSequencingModelVariables vars_;
};

class VehicleSlotsModel : public Model
{
public:
  explicit VehicleSlotsModel(Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed, bool export_model);
  ~VehicleSlotsModel() = default;
  virtual void solve(const Instance &instance, double time_limit);

private:
  virtual void allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars = false);
  virtual void populateByRow(const Instance &instance, bool reformulate, bool symmetry_breaking, bool export_model);
  virtual void addCut(UserCut *curr_cut);

  VehicleSlotsModelVariables vars_;
};