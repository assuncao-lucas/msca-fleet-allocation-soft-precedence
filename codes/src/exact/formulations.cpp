#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/exact/formulations.h"

ILOSTLBEGIN

std::ostream &
operator<<(std::ostream &out, const VarToIndexMap &map)
{
  out << "index to map:" << std::endl;
  for (auto &[key, value] : map.index_to_var_map())
  {
    out << key << " -> (" << value.first << " ," << value.second << ")" << std::endl;

    if (key != map.var_to_index_map().find(value)->second)
      throw "inconsistency map vars";
  }

  out << "map to index:" << std::endl;
  for (auto &[key, value] : map.var_to_index_map())
  {
    out << "(" << key.first << " ," << key.second << ") -> " << value << std::endl;
    auto iter = map.index_to_var_map().find(value);
    if (key.first != iter->second.first || key.second != iter->second.second)
      throw "inconsistency map vars";
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const VehicleSequencingModelVariables &vars)
{
  out << "x vars:" << std::endl
      << vars.x_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "z vars:" << std::endl
      << vars.z_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "w vars:" << std::endl
      << vars.w_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "U vars:" << std::endl
      << vars.U_var_to_index_ << std::endl
      << "***************" << std::endl;

  return out;
}

void allocateVehicleSequencingModelVariables(IloEnv &env, VehicleSequencingModelVariables &vars, const Instance &instance, bool solve_relax, bool disable_all_binary_vars)
{
  const int num_items = instance.num_items();
  const int num_items_for_transport = instance.num_items_for_transport();
  const int num_vehicles = instance.num_vehicles();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();

  double binary_upper_bound = disable_all_binary_vars ? 0.0 : 1.0;
  auto type = solve_relax ? ILOFLOAT : ILOINT;

  // fill x variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance.items_for_transport())
      vars.x_var_to_index_.addEntry(i, j);

  vars.x_ = IloNumVarArray(env, vars.x_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill z variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (int h = 0; h < num_vehicles; ++h)
      if (i != h)
        vars.z_var_to_index_.addEntry(i, h);

  vars.z_ = IloNumVarArray(env, vars.z_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill w variables and index map.
  for (auto i : instance.items_for_transport())
  {
    for (auto &k : successors_for_transport_per_item[i])
      vars.w_var_to_index_.addEntry(i, k);
  }

  vars.w_ = IloNumVarArray(env, vars.w_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill U variables and index map.
  for (int i = 0; i < num_items; ++i)
    if (!(successors_for_transport_per_item[i].empty()))
      vars.U_var_to_index_.addEntry(0, i);

  vars.U_ = IloNumVarArray(env, vars.U_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
}

void vehicleSequencingModel(Instance &inst, double time_limit, bool add_symmetry_breaking, bool solve_relax, bool export_model)
{
  const int num_vertices = inst.num_items();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  VehicleSequencingModelVariables vars{};

  allocateVehicleSequencingModelVariables(env, vars, inst, solve_relax);

  // std::cout << vars << std::endl;

  populateByRowVehicleSequencingModel(cplex, env, model, vars, inst, add_symmetry_breaking, export_model);

  Solution<double> solution;
  optimize(cplex, env, model, time_limit, solve_relax, solution);
  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  if (solve_relax)
  {
    std::cout << "LP: " << solution.lp_ << std::endl;
  }

  if (solution.is_feasible_)
  {

    int num_vehicles = inst.num_vehicles();
    int num_items = inst.num_items();
    const auto &successors_for_transport_per_item = inst.successors_for_transport_per_item();

    IloNumArray x_values(env), z_values(env), w_values(env), U_values(env);
    cplex.getValues(x_values, vars.x_);
    cplex.getValues(z_values, vars.z_);
    cplex.getValues(w_values, vars.w_);
    cplex.getValues(U_values, vars.U_);

    // get values of x variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (auto j : inst.items_for_transport())
      {
        auto value = x_values[vars.x_var_to_index_.varToIndex(i, j)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "x[" << i << "," << j << "] = " << value << std::endl;
          num_items_loaded++;
        }
      }
    }

    // get values of z variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (int h = 0; h < num_vehicles; ++h)
      {
        if (i != h)
        {
          auto value = z_values[vars.z_var_to_index_.varToIndex(i, h)];
          if (!double_equals(value, 0.0))
            std::cout << "z[" << i << "," << h << "] = " << value << std::endl;
        }
      }
    }

    // get values of w variables.
    for (auto i : inst.items_for_transport())
    {
      for (auto &k : successors_for_transport_per_item[i])
      {
        auto value = w_values[vars.z_var_to_index_.varToIndex(i, k)];
        if (!double_equals(value, 0.0))
          std::cout << "w[" << i << "," << k << "] = " << value << std::endl;
      }
    }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = U_values[vars.U_var_to_index_.varToIndex(0, i)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    std::cout << "num_items_loaded: " << num_items_loaded << std::endl
              << "num_unproductive_moves: " << num_unproductive_moves << std::endl;
    x_values.end();
    z_values.end();
    w_values.end();
    U_values.end();
  }
  cplex.end();
  env.end();
}

void populateByRowVehicleSequencingModel(IloCplex &cplex, IloEnv &env, IloModel &model, VehicleSequencingModelVariables &vars, const Instance &instance, bool add_symmetry_breaking, bool export_model)
{

  IloExpr obj(env);
  const int num_vehicles = instance.num_vehicles();
  const auto &items = instance.items();
  const int num_items = instance.num_items();
  const auto &successors_fixed_per_item = instance.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();
  const auto &fleet = instance.fleet();
  int M = 0;

  // compute M.
  for (int j = 0; j < num_items; ++j)
    if (!successors_fixed_per_item[j].empty() || !successors_for_transport_per_item[j].empty())
      ++M;

  // add objective function.
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (item->available_for_transport())
    {
      for (int i = 0; i < num_vehicles; ++i)
      {
        obj += operator*(M, vars.x(i, j));
      }

      if (!successors_fixed_per_item[j].empty() || !successors_for_transport_per_item[j].empty())
        obj -= vars.U(j);
    }
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  // Restrictions (6).
  for (auto &j : instance.items_for_transport())
  {
    IloExpr exp(env);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars.x(i, j);
    model.add(exp <= 1);
    exp.end();
  }

  // Restrictions (7).
  for (auto &j : instance.items_for_transport())
  {
    for (auto &k : instance.items_for_transport())
    {
      if ((j < k) && (items[j]->group() != items[k]->group()))
      {
        for (int i = 0; i < num_vehicles; ++i)
        {
          IloExpr exp(env);
          model.add(vars.x(i, j) + vars.x(i, k) <= 1);
          exp.end();
        }
      }
    }
  }

  // Restrictions (8).
  for (int i = 0; i < num_vehicles; ++i)
  {
    IloExpr exp(env);
    for (auto &j : instance.items_for_transport())
      exp += vars.x(i, j);
    model.add(exp <= fleet[i]->capacity());
    exp.end();
  }

  // Restrictions (9) and (10).
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (!successors_fixed_per_item[j].empty() || !successors_for_transport_per_item[j].empty())
    {
      IloExpr exp(env);
      if (item->available_for_transport())
      {
        for (auto &k : successors_for_transport_per_item[j])
          exp += vars.w(j, k);
      }
      else
      {
        for (auto &k : successors_for_transport_per_item[j])
          for (int i = 0; i < num_vehicles; ++i)
            exp += vars.x(i, k);
      }
      model.add(operator*((int)(successors_for_transport_per_item[j]).size(), vars.U(j)) >= exp);
      exp.end();
    }
  }

  // Restrictions (11).
  for (int i = 0; i < num_vehicles; ++i)
    for (int h = i + 1; h < num_vehicles; ++h)
      model.add(vars.z(i, h) + vars.z(h, i) == 1);

  // Restrictions (12) == (19) and (20).
  for (int i = 0; i < num_vehicles; ++i)
  {
    for (int h = i + 1; h < num_vehicles; ++h)
    {
      for (int s = h + 1; s < num_vehicles; ++s)
      {
        model.add(vars.z(i, h) + vars.z(h, s) + vars.z(s, i) <= 2);
        model.add(vars.z(i, s) + vars.z(s, h) + vars.z(h, i) <= 2);
      }
    }
  }

  // Restrictions (13) and (14).
  for (auto &j : instance.items_for_transport())
  {
    for (auto &k : successors_for_transport_per_item[j])
    {
      IloExpr sum_x_i_j(env), sum_x_i_k(env);
      for (int i = 0; i < num_vehicles; ++i)
      {
        sum_x_i_j += vars.x(i, j);
        sum_x_i_k += vars.x(i, k);
        for (int h = 0; h < num_vehicles; ++h)
          if (i != h)
            model.add(vars.w(j, k) >= vars.x(i, j) + vars.x(h, k) + vars.z(h, i) - 2);
      }

      model.add(vars.w(j, k) >= sum_x_i_k - sum_x_i_j);
      sum_x_i_j.end();
      sum_x_i_k.end();
    }
  }

  if (add_symmetry_breaking)
  {
    // Inequalities (21)-(23).
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (int h = i + 1; h < num_vehicles; ++h)
      {
        IloExpr sum_x_i_j(env), sum_x_h_j(env);
        for (auto &j : instance.items_for_transport())
        {
          sum_x_i_j += vars.x(i, j);
          sum_x_h_j += vars.x(h, j);
        }
        model.add(vars.z(i, h) >= 1 - sum_x_i_j - sum_x_h_j);
        model.add(vars.z(i, h) >= sum_x_i_j - sum_x_h_j);
        model.add(vars.z(i, h) >= sum_x_h_j - sum_x_i_j);

        sum_x_i_j.end();
        sum_x_h_j.end();
      }
    }
  }

  if (export_model)
  {
    // add name to variables.
    for (auto &[key, value] : vars.x_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "x(%d)(%d)", value.first, value.second);
      vars.x_[key].setName(strnum);
    }

    for (auto &[key, value] : vars.z_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "z(%d)(%d)", value.first, value.second);
      vars.z_[key].setName(strnum);
    }

    for (auto &[key, value] : vars.w_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "w(%d)(%d)", value.first, value.second);
      vars.w_[key].setName(strnum);
    }

    for (auto &[key, value] : vars.U_var_to_index_.index_to_var_map())
    {
      char strnum[15];
      sprintf(strnum, "U(%d)", value.second);
      vars.U_[key].setName(strnum);
    }

    cplex.exportModel("vehicle_sequencing_model.lp");
  }
}

void optimize(IloCplex &cplex, IloEnv &env, IloModel &model, double total_time_limit, bool solve_relax, Solution<double> &solution)
{
  Timestamp *ti = NewTimestamp(), *tf = NewTimestamp();
  Timer *timer = GetTimer();
  timer->Clock(ti);
  cplex.setParam(IloCplex::Param::WorkMem, 100000);
  std::cout << "limit of memory 100000MB" << std::endl;
  cplex.setParam(IloCplex::IloCplex::Param::MIP::Strategy::File, 3);
  cplex.setOut(env.getNullStream());

  // if solving relaxed problem, force dual cplex method (multithreading showed to be slower! Dual cplex forces single threading and is also useful when
  // the main LP is solved several times, either on Benders or while separating valid inequalities).
  if (solve_relax)
  {
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);
  }

  if (solve_relax)
    total_time_limit = -1;

  if (!double_equals(total_time_limit, -1))
  {
    // std::cout << total_time_limit - instance.time_spent_in_preprocessing() << std::endl;
    cplex.setParam(IloCplex::Param::ClockType, 2);
    double time_left = total_time_limit;
    cplex.setParam(IloCplex::Param::TimeLimit, time_left);
  }

  // Optimize the problem and obtain solution.
  cplex.solve();

  timer->Clock(tf);
  if (solve_relax)
    solution.root_time_ = timer->ElapsedTime(ti, tf);
  else
    solution.milp_time_ += timer->ElapsedTime(ti, tf);

  // if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
  //   solution.is_feasible_ = false;
  // else
  SetSolutionStatus(cplex, solution, solve_relax);

  delete (ti);
  ti = nullptr;
  delete (tf);
  tf = nullptr;
}