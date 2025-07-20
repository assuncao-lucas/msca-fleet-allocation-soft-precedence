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

  int varToIndex(int i, int j) const
  {
    auto element = var_to_index_map_.find({i, j});
    if (element == var_to_index_map_.end())
      throw "tried to access var not mapped: " + std::to_string(i) + " " + std::to_string(j);
    else
      return element->second;
  }

  std::pair<int, int> indexToVar(int index) const
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

  void fill_all_vars_array(IloNumVarArray &all_vars)
  {
    all_vars.add(x_);
    all_vars.add(z_);
    all_vars.add(w_);
    all_vars.add(U_);
    all_vars.add(X_);
  }

  int num_vars()
  {
    return x_.getSize() + z_.getSize() + w_.getSize() + U_.getSize() + X_.getSize();
  }

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

  void fill_all_vars_array(IloNumVarArray &all_vars)
  {
    all_vars.add(x_);
    all_vars.add(u_);
    all_vars.add(U_);
    all_vars.add(X_);
  }

  int num_vars()
  {
    return x_.getSize() + u_.getSize() + U_.getSize() + X_.getSize();
  }

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

  void fill_all_vars_array(IloNumVarArray &all_vars)
  {
    all_vars.add(y1_);
    all_vars.add(y2_);
    all_vars.add(U_);
    all_vars.add(Y_);
  }

  int num_vars()
  {
    return y1_.getSize() + y2_.getSize() + U_.getSize() + Y_.getSize();
  }

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

static UserCut *GenerateCliqueConflictCuts(Instance &instance, const VarToIndexMap &var_to_index_map, const IloNumVarArray &item_to_vehicle_or_slot_vars, const IloNumArray &item_to_vehicle_or_slot_values, bool root_cuts, Solution<double> &sol, std::list<UserCut *> &cuts);

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
  Model(const Instance &instance);
  virtual ~Model();
  virtual void fillSolution(Solution<double> &solution) = 0;
  bool optimize(double total_time_limit, bool find_root_cuts, std::list<UserCut *> *initial_cuts, std::list<UserCut *> *root_cuts, Solution<double> &solution);
  void set_multithreading(bool multithreading);
  void set_emphasis_feasibility()
  {
    cplex_->setParam(IloCplex::Param::Emphasis::MIP, IloCplex::MIPEmphasisFeasibility);
  }

  void set_time_limit(double time_limit)
  {
    cplex_->setParam(IloCplex::Param::TimeLimit, time_limit);
  }

  void addMIPStart(const IloNumArray &curr_mip_start_vals)
  {
    cplex_->addMIPStart(all_vars_array_, curr_mip_start_vals, IloCplex::MIPStartSolveMIP);
  }

  IloNum getObjValue()
  {
    return cplex_->getObjValue();
  }

  void getSolutionValues(IloNumArray &values)
  {
    values.clear();
    cplex_->getValues(values, all_vars_array_);
  }

  IloEnv &env()
  {
    return *env_;
  }
  IloCplex::CplexStatus status()
  {
    return cplex_->getCplexStatus();
  }

  void exportModel(const char *name) const
  {
    cplex_->exportModel(name);
  }

  virtual Model *getClone(bool relaxed) = 0;

  virtual void getSolutionItems(boost::dynamic_bitset<> &curr_items) = 0;
  virtual int num_vars() = 0;
  virtual IloNumArray get_items_values() = 0;
  virtual IloNumArray get_items_reduced_costs() = 0;
  virtual void UpdateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active) = 0;

protected:
  virtual void allocateVariables(bool reformulate) = 0;
  virtual void populateByRow(bool reformulate, bool symmetry_breaking) = 0;
  virtual void addInitialCuts(std::list<UserCut *> *initial_cuts, Solution<double> &solution);
  virtual void addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts) = 0;
  virtual bool findAndAddValidInqualities(Solution<double> &sol, std::list<UserCut *> *root_cuts);
  virtual UserCut *separateCuts(bool root_cuts, Solution<double> &sol, std::list<UserCut *> &cuts) = 0;

  IloEnv *env_ = nullptr;     // Cplex environment.
  IloCplex *cplex_ = nullptr; // Cplex solver.
  IloModel *model_ = nullptr; // Cplex model.
  IloNumVarArray all_vars_array_;
  bool is_relaxed_ = false;
  bool reformulated_ = false;
  bool symmetry_breaking_ = false;

  const Instance &instance_;
};

class VehicleSequencingModel : public Model
{
public:
  explicit VehicleSequencingModel(const Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed);
  virtual ~VehicleSequencingModel() = default;
  virtual void fillSolution(Solution<double> &solution);
  virtual int num_vars()
  {
    return vars_.num_vars();
  }
  virtual IloNumArray get_items_values();
  virtual IloNumArray get_items_reduced_costs();

  virtual Model *getClone(bool relaxed)
  {
    return new VehicleSequencingModel(instance_, reformulated_, symmetry_breaking_, relaxed);
  }

  void UpdateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active);
  virtual void getSolutionItems(boost::dynamic_bitset<> &curr_items);

private:
  virtual void allocateVariables(bool reformulate);
  virtual void populateByRow(bool reformulate, bool symmetry_breaking);
  virtual void addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts);
  virtual UserCut *separateCuts(bool root_cuts, Solution<double> &sol, std::list<UserCut *> &cuts);

  VehicleSequencingModelVariables vars_;
};

class ItemSequencingModel : public Model
{
public:
  explicit ItemSequencingModel(const Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed);
  virtual ~ItemSequencingModel() = default;
  virtual void fillSolution(Solution<double> &solution);
  virtual int num_vars()
  {
    return vars_.num_vars();
  }

  virtual IloNumArray get_items_values();
  virtual IloNumArray get_items_reduced_costs();

  virtual Model *getClone(bool relaxed)
  {
    return new ItemSequencingModel(instance_, reformulated_, symmetry_breaking_, relaxed);
  }

  void UpdateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active);
  virtual void getSolutionItems(boost::dynamic_bitset<> &curr_items);

private:
  virtual void allocateVariables(bool reformulate);
  virtual void populateByRow(bool reformulate, bool symmetry_breaking);
  virtual void addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts);
  virtual UserCut *separateCuts(bool root_cuts, Solution<double> &sol, std::list<UserCut *> &cuts);

  ItemSequencingModelVariables vars_;
};

class VehicleSlotsModel : public Model
{
public:
  explicit VehicleSlotsModel(const Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed);
  virtual ~VehicleSlotsModel() = default;
  virtual void fillSolution(Solution<double> &solution);
  virtual int num_vars()
  {
    return vars_.num_vars();
  }

  virtual IloNumArray get_items_values();
  virtual IloNumArray get_items_reduced_costs();

  virtual Model *getClone(bool relaxed)
  {
    return new VehicleSlotsModel(instance_, reformulated_, symmetry_breaking_, relaxed);
  }

  void UpdateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active);
  virtual void getSolutionItems(boost::dynamic_bitset<> &curr_items);

private:
  virtual void allocateVariables(bool reformulate);
  virtual void populateByRow(bool reformulate, bool symmetry_breaking);
  virtual void addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts);
  virtual UserCut *separateCuts(bool root_cuts, Solution<double> &sol, std::list<UserCut *> &cuts);

  VehicleSlotsModelVariables vars_;
};