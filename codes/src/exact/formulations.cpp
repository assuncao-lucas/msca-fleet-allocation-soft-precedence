#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/exact/formulations.h"

ILOSTLBEGIN

// Model.
Model::~Model()
{
  if (cplex_)
  {
    cplex_->end();
    delete cplex_;
  }

  if (model_)
  {
    model_->end();
    delete model_;
  }

  if (env_)
  {
    env_->end();
    delete env_;
  }
}

void Model::optimize(const Instance &instance, double total_time_limit, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCut *> *initial_cuts, std::list<UserCut *> *root_cuts, Solution<double> &solution)
{
  Timestamp *ti = NewTimestamp(), *tf = NewTimestamp();
  Timer *timer = GetTimer();
  timer->Clock(ti);
  bool found_cuts = false;

  cplex_->setParam(IloCplex::Param::WorkMem, 100000);
  std::cout << "limit of memory 100000MB" << std::endl;
  cplex_->setParam(IloCplex::IloCplex::Param::MIP::Strategy::File, 3);
  cplex_->setOut(env_->getNullStream());

  addInitialCuts(initial_cuts, root_cuts, solution);

  // if solving relaxed problem, force dual cplex method (multithreading showed to be slower! Dual cplex forces single threading and is also useful when
  // the main LP is solved several times, either on Benders or while separating valid inequalities).
  if (is_relaxed_)
  {
    cplex_->setParam(IloCplex::Param::Threads, 1);
    cplex_->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);
  }

  if (is_relaxed_)
    total_time_limit = -1;

  if (!double_equals(total_time_limit, -1))
  {
    // std::cout << total_time_limit - instance.time_spent_in_preprocessing() << std::endl;
    cplex_->setParam(IloCplex::Param::ClockType, 2);
    double time_left = total_time_limit;
    cplex_->setParam(IloCplex::Param::TimeLimit, time_left);
  }

  double previous_bound = 0.0, curr_bound = std::numeric_limits<double>::infinity();
  do
  {
    previous_bound = curr_bound;
    // Optimize the problem and obtain solution.
    if (!cplex_->solve())
    {
      timer->Clock(tf);

      std::cout << cplex_->getCplexStatus() << std::endl;

      if (is_relaxed_)
        solution.root_time_ = timer->ElapsedTime(ti, tf);
      else
        solution.milp_time_ += timer->ElapsedTime(ti, tf);
      if ((cplex_->getCplexStatus() == IloCplex::Infeasible) || (cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
        solution.is_feasible_ = false;
      // SetSolutionStatus(cplex,solution,solve_relax);

      delete (ti);
      ti = nullptr;
      delete (tf);
      tf = nullptr;
      return;
    }
    // getchar(); getchar();
    // cont++;
    // std::cout << cont << std::endl;
    found_cuts = false;

    if (is_relaxed_)
    {
      curr_bound = cplex_->getObjValue();

      if (find_root_cuts && double_less(curr_bound, previous_bound, K_TAILING_OFF_TOLERANCE))
        found_cuts |= findAndAddValidInqualities(instance, solution, root_cuts);
    }
  } while (found_cuts);

  // if ((cplex_->getCplexStatus() == IloCplex::Infeasible) || (cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
  //   solution.is_feasible_ = false;
  // else
  SetSolutionStatus(*cplex_, solution, is_relaxed_);

  delete (ti);
  ti = nullptr;
  delete (tf);
  tf = nullptr;
}

void Model::addInitialCuts(std::list<UserCut *> *initial_cuts, std::list<UserCut *> *root_cuts, Solution<double> &solution)
{
  if ((initial_cuts != nullptr) && (!(initial_cuts->empty())))
  {
    IloRangeArray root_cuts(*env_);
    for (std::list<UserCut *>::iterator it = initial_cuts->begin(); it != initial_cuts->end(); it++)
    {
      UserCut *curr_user_cut = static_cast<UserCut *>((*it));
      addCut(curr_user_cut);
      solution.set_cut_added(-1, false);

      // delete (*it);
      //*it = NULL;
    }

    is_relaxed_ ? model_->add(root_cuts) : cplex_->addUserCuts(root_cuts);
    root_cuts.endElements();
    root_cuts.end();
  }
}

bool Model::findAndAddValidInqualities(const Instance &instance, Solution<double> &sol, std::list<UserCut *> *root_cuts)
{
  // std::cout << " ***************** looking for cuts" << std::endl;
  bool found_cut = false;
  (sol.num_calls_to_callback_lp_) += 1;

  std::list<UserCut *> cuts;

  UserCut *best_cut = nullptr, *curr_cut = nullptr, *local_best_cut = nullptr;

  best_cut = nullptr;
  // find clique conflict cuts.
  // local_best_cut = GenerateCliqueConflictCuts(instance, visited_nodes, nodes_sum, values_x,
  //                                             visited_nodes_list, subgraph_to_graph_map, graph_to_subgraph_map, g,
  //                                             l_nodes, capacity, g_inv, l_nodes_inv, capacity_inv, sol, cuts, true, instance.active_conflict_cliques());

  if ((local_best_cut != nullptr) && (local_best_cut->isBetterThan(best_cut)))
    best_cut = local_best_cut;

  if (best_cut != nullptr)
  {
    found_cut = true;

    // adds best cut (most violated).
    addCut(best_cut);
    sol.set_cut_added(-1, true);
    // std::cout << "added cut to LP: " << std::endl;
    // for (auto i : best_cut->rhs_nonzero_coefficients_indexes_)
    //   std::cout << i << std::endl;
    // std::cout << std::endl;

    // for (auto [i, j] : best_cut->lhs_nonzero_coefficients_indexes_)
    //   std::cout << i << ", " << j << std::endl;
    // std::cout << std::endl;

    if (root_cuts != nullptr)
    {
      root_cuts->push_back(best_cut);
    }

    // checks the angle between the each Cut and the most violated one
    for (std::list<UserCut *>::iterator it = cuts.begin(); it != cuts.end(); ++it)
    {
      curr_cut = static_cast<UserCut *>(*it);

      if (curr_cut == best_cut)
        continue;

      // std::cout << *curr_cut << std::endl;
      //  adds only the ones sufficiently orthogonal to the best_cut
      // std::cout << (*curr_cut)*(*best_cut) << std::endl;
      // getchar(); getchar();
      bool is_sufficiently_orthogonal = double_less((*curr_cut) * (*best_cut), K_CUTS_ANGLE_COSIN_LIMIT);
      if (is_sufficiently_orthogonal)
      {
        addCut(curr_cut);

        sol.set_cut_added(-1, true);
        // std::cout << "added cut to LP 2: " << std::endl;
        // for (auto i : curr_cut->rhs_nonzero_coefficients_indexes_)
        //   std::cout << i << std::endl;
        // std::cout << std::endl;

        // for (auto [i, j] : curr_cut->lhs_nonzero_coefficients_indexes_)
        //   std::cout << i << ", " << j << std::endl;
        // std::cout << std::endl;

        if (root_cuts != nullptr)
        {
          root_cuts->push_back(curr_cut);
          // sol.set_cut_added(curr_cut->type_,false);
        }
      }

      if ((!is_sufficiently_orthogonal) || (root_cuts == nullptr))
      {
        delete curr_cut;
        *it = nullptr;
      }
    }

    if (root_cuts == nullptr)
    {
      delete best_cut;
      best_cut = nullptr;
    }
    cuts.clear();
  }

  return found_cut;
}

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

// Vehicle Sequencing model.
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
  out << "X vars:" << std::endl
      << vars.X_var_to_index_ << std::endl
      << "***************" << std::endl;
  return out;
}

void VehicleSequencingModel::addCut(UserCut *curr_cut)
{
  IloExpr exp(*env_);
  int vehicle = -1, item = -1;
  for (std::list<std::pair<int, int>>::iterator it = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it != curr_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
  {
    vehicle = (*it).first;
    item = (*it).second;
    exp += vars_.x(vehicle, item);
  }

  model_->add(exp <= 1);
  exp.end();
}

void VehicleSequencingModel::allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars)
{
  const int num_items = instance.num_items();
  const int num_items_for_transport = instance.num_items_for_transport();
  const int num_vehicles = instance.num_vehicles();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();

  double binary_upper_bound = disable_all_binary_vars ? 0.0 : 1.0;
  auto type = is_relaxed_ ? ILOFLOAT : ILOINT;

  // fill x variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance.items_for_transport())
      vars_.x_var_to_index_.addEntry(i, j);

  vars_.x_ = IloNumVarArray(*env_, vars_.x_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill z variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (int h = 0; h < num_vehicles; ++h)
      if (i != h)
        vars_.z_var_to_index_.addEntry(i, h);

  vars_.z_ = IloNumVarArray(*env_, vars_.z_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill w variables and index map.
  for (auto i : instance.items_for_transport())
  {
    for (auto &k : successors_for_transport_per_item[i])
      vars_.w_var_to_index_.addEntry(i, k);
  }

  vars_.w_ = IloNumVarArray(*env_, vars_.w_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill U variables and index map.
  for (int i = 0; i < num_items; ++i)
    if (!(successors_for_transport_per_item[i].empty()))
      vars_.U_var_to_index_.addEntry(0, i);

  vars_.U_ = IloNumVarArray(*env_, vars_.U_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  if (reformulate)
  {
    // fill X variables and index map.
    for (int i = 0; i < num_vehicles; ++i)
      for (auto &[group, _] : instance.items_for_transport_per_group())
        vars_.X_var_to_index_.addEntry(i, group);

    vars_.X_ = IloNumVarArray(*env_, vars_.X_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
  }
}

VehicleSequencingModel::VehicleSequencingModel(Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed, bool export_model)
{
  is_relaxed_ = relaxed;
  reformulated_ = reformulate;

  env_ = new IloEnv();
  model_ = new IloModel(*env_);
  cplex_ = new IloCplex(*env_);

  cplex_->extract(*model_);

  allocateVariables(inst, reformulate);

  std::cout << vars_ << std::endl;

  populateByRow(inst, reformulate, symmetry_breaking, export_model);
}

void VehicleSequencingModel::populateByRow(const Instance &instance, bool reformulate, bool symmetry_breaking, bool export_model)
{
  IloExpr obj(*env_);
  const int num_vehicles = instance.num_vehicles();
  const auto &items = instance.items();
  const int num_items = instance.num_items();
  // const auto &successors_fixed_per_item = instance.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();
  const auto &fleet = instance.fleet();
  int M = 0;

  // compute M.
  for (int j = 0; j < num_items; ++j)
    if (!successors_for_transport_per_item[j].empty())
      ++M;

  // add objective function.
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (item->available_for_transport())
    {
      for (int i = 0; i < num_vehicles; ++i)
      {
        obj += operator*(M, vars_.x(i, j));
      }
    }

    // even if not available for transport.
    if (!successors_for_transport_per_item[j].empty())
      obj -= vars_.U(j);
  }

  model_->add(IloMaximize(*env_, obj));
  obj.end();

  // Restrictions (6).
  for (auto &j : instance.items_for_transport())
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.x(i, j);
    model_->add(exp <= 1);
    exp.end();
  }

  if (reformulate)
  {
    // Constraints (22) and (23).
    for (int i = 0; i < num_vehicles; ++i)
    {
      IloExpr exp(*env_);
      for (auto &[group, items] : instance.items_for_transport_per_group())
      {
        exp += vars_.X(i, group);
        IloExpr sum_items(*env_);

        for (auto j : items)
          sum_items += vars_.x(i, j);

        model_->add(operator*((int)items.size(), vars_.X(i, group)) >= sum_items);
        sum_items.end();
      }

      model_->add(exp <= 1);
      exp.end();
    }
  }
  else
  {
    // Restrictions (7).
    for (auto &j : instance.items_for_transport())
    {
      for (auto &k : instance.items_for_transport())
      {
        if ((j < k) && (items[j]->group() != items[k]->group()))
        {
          for (int i = 0; i < num_vehicles; ++i)
          {
            IloExpr exp(*env_);
            model_->add(vars_.x(i, j) + vars_.x(i, k) <= 1);
            exp.end();
          }
        }
      }
    }
  }

  // Restrictions (8).
  for (int i = 0; i < num_vehicles; ++i)
  {
    IloExpr exp(*env_);
    for (auto &j : instance.items_for_transport())
      exp += vars_.x(i, j);
    model_->add(exp <= fleet[i]->capacity());
    exp.end();
  }

  // Restrictions (9) and (10).
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (!(successors_for_transport_per_item[j].empty()))
    {
      IloExpr exp(*env_);
      if (item->available_for_transport())
      {
        for (auto &k : successors_for_transport_per_item[j])
          exp += vars_.w(j, k);
      }
      else
      {
        for (auto &k : successors_for_transport_per_item[j])
          for (int i = 0; i < num_vehicles; ++i)
            exp += vars_.x(i, k);
      }
      model_->add(operator*((int)(successors_for_transport_per_item[j]).size(), vars_.U(j)) >= exp);
      exp.end();
    }
  }

  // Restrictions (11).
  for (int i = 0; i < num_vehicles; ++i)
    for (int h = i + 1; h < num_vehicles; ++h)
      model_->add(vars_.z(i, h) + vars_.z(h, i) == 1);

  // Restrictions (12) == (19) and (20).
  for (int i = 0; i < num_vehicles; ++i)
  {
    for (int h = i + 1; h < num_vehicles; ++h)
    {
      for (int s = h + 1; s < num_vehicles; ++s)
      {
        model_->add(vars_.z(i, h) + vars_.z(h, s) + vars_.z(s, i) <= 2);
        model_->add(vars_.z(i, s) + vars_.z(s, h) + vars_.z(h, i) <= 2);
      }
    }
  }

  // Restrictions (13) and (14).
  for (auto &j : instance.items_for_transport())
  {
    for (auto &k : successors_for_transport_per_item[j])
    {
      IloExpr sum_x_i_j(*env_), sum_x_i_k(*env_);
      for (int i = 0; i < num_vehicles; ++i)
      {
        sum_x_i_j += vars_.x(i, j);
        sum_x_i_k += vars_.x(i, k);
        for (int h = 0; h < num_vehicles; ++h)
          if (i != h)
            model_->add(vars_.w(j, k) >= vars_.x(i, j) + vars_.x(h, k) + vars_.z(h, i) - 2);
      }

      model_->add(vars_.w(j, k) >= sum_x_i_k - sum_x_i_j);
      sum_x_i_j.end();
      sum_x_i_k.end();
    }
  }

  if (reformulate && symmetry_breaking) // in this case, symmetry breaking constraint rely on X vars_.
  {
    // Inequalities (24)-(26).
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (int h = i + 1; h < num_vehicles; ++h)
      {
        IloExpr sum_X_i_q(*env_), sum_X_h_q(*env_);
        for (auto &[q, _] : instance.items_for_transport_per_group())
        {
          sum_X_i_q += vars_.X(i, q);
          sum_X_h_q += vars_.X(h, q);
        }
        model_->add(vars_.z(i, h) >= 1 - sum_X_i_q - sum_X_h_q);
        model_->add(vars_.z(i, h) >= sum_X_i_q - sum_X_h_q);
        model_->add(vars_.z(h, i) >= sum_X_h_q - sum_X_i_q);

        sum_X_i_q.end();
        sum_X_h_q.end();
      }
    }
  }

  if (export_model)
  {
    // add name to variables.
    for (auto &[key, value] : vars_.x_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "x(%d)(%d)", value.first, value.second);
      vars_.x_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.z_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "z(%d)(%d)", value.first, value.second);
      vars_.z_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.w_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "w(%d)(%d)", value.first, value.second);
      vars_.w_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.U_var_to_index_.index_to_var_map())
    {
      char strnum[15];
      sprintf(strnum, "U(%d)", value.second);
      vars_.U_[key].setName(strnum);
    }

    if (reformulate)
    {
      for (auto &[key, value] : vars_.X_var_to_index_.index_to_var_map())
      {
        char strnum[26];
        sprintf(strnum, "X(%d)(%d)", value.first, value.second);
        vars_.X_[key].setName(strnum);
      }
    }
    cplex_->exportModel("vehicle_sequencing_model.lp");
  }
}

void VehicleSequencingModel::solve(const Instance &inst, double time_limit)
{
  Solution<double> solution;
  optimize(inst, time_limit, false, false, nullptr, nullptr, solution);

  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  if (is_relaxed_)
  {
    std::cout << "LP: " << solution.lp_ << std::endl;
  }

  if (solution.is_feasible_)
  {

    int num_vehicles = inst.num_vehicles();
    int num_items = inst.num_items();
    const auto &successors_for_transport_per_item = inst.successors_for_transport_per_item();

    IloNumArray x_values(*env_), z_values(*env_), w_values(*env_), U_values(*env_);
    cplex_->getValues(x_values, vars_.x_);
    cplex_->getValues(z_values, vars_.z_);
    cplex_->getValues(w_values, vars_.w_);
    cplex_->getValues(U_values, vars_.U_);

    // get values of x variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (auto j : inst.items_for_transport())
      {
        auto value = x_values[vars_.x_index(i, j)];
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
          auto value = z_values[vars_.z_index(i, h)];
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
        auto value = w_values[vars_.w_index(i, k)];
        if (!double_equals(value, 0.0))
          std::cout << "w[" << i << "," << k << "] = " << value << std::endl;
      }
    }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = U_values[vars_.U_index(i)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    if (reformulated_)
    {
      IloNumArray X_values(*env_);
      cplex_->getValues(X_values, vars_.X_);

      // get values of X variables.
      for (int i = 0; i < num_vehicles; ++i)
      {
        for (auto &[group, _] : inst.items_for_transport_per_group())
        {
          auto value = X_values[vars_.X_index(i, group)];
          if (!double_equals(value, 0.0))
            std::cout << "X[" << i << "," << group << "] = " << value << std::endl;
        }
      }
      X_values.end();
    }

    std::cout << "num_items_loaded: " << num_items_loaded << std::endl
              << "num_unproductive_moves: " << num_unproductive_moves << std::endl;
    x_values.end();
    z_values.end();
    w_values.end();
    U_values.end();
  }
}

// Item Sequencing model.
std::ostream &
operator<<(std::ostream &out, const ItemSequencingModelVariables &vars)
{
  out << "x vars:" << std::endl
      << vars.x_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "u vars:" << std::endl
      << vars.u_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "U vars:" << std::endl
      << vars.U_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "X vars:" << std::endl
      << vars.X_var_to_index_ << std::endl
      << "***************" << std::endl;
  return out;
}

void ItemSequencingModel::addCut(UserCut *curr_cut)
{
  IloExpr exp(*env_);
  int vehicle = -1, item = -1;
  for (std::list<std::pair<int, int>>::iterator it = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it != curr_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
  {
    vehicle = (*it).first;
    item = (*it).second;
    exp += vars_.x(vehicle, item);
  }

  model_->add(exp <= 1);
  exp.end();
}

void ItemSequencingModel::allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars)
{
  const int num_items = instance.num_items();
  const int num_items_for_transport = instance.num_items_for_transport();
  const int num_vehicles = instance.num_vehicles();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();

  double binary_upper_bound = disable_all_binary_vars ? 0.0 : 1.0;
  auto type = is_relaxed_ ? ILOFLOAT : ILOINT;

  // fill x variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance.items_for_transport())
      vars_.x_var_to_index_.addEntry(i, j);

  vars_.x_ = IloNumVarArray(*env_, vars_.x_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill u variables and index map.
  for (auto j : instance.items_for_transport())
    for (auto k : instance.items_for_transport())
      if (j != k)
        vars_.u_var_to_index_.addEntry(j, k);

  vars_.u_ = IloNumVarArray(*env_, vars_.u_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill U variables and index map.
  for (int i = 0; i < num_items; ++i)
    if (!(successors_for_transport_per_item[i].empty()))
      vars_.U_var_to_index_.addEntry(0, i);

  vars_.U_ = IloNumVarArray(*env_, vars_.U_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  if (reformulate)
  {
    // fill X variables and index map.
    for (int i = 0; i < num_vehicles; ++i)
      for (auto &[group, _] : instance.items_for_transport_per_group())
        vars_.X_var_to_index_.addEntry(i, group);

    vars_.X_ = IloNumVarArray(*env_, vars_.X_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
  }
}

ItemSequencingModel::ItemSequencingModel(Instance &inst, bool reformulate, bool add_symmetry_breaking, bool solve_relax, bool export_model)
{
  is_relaxed_ = solve_relax;
  reformulated_ = reformulate;
  const int num_vertices = inst.num_items();

  env_ = new IloEnv();
  model_ = new IloModel(*env_);
  cplex_ = new IloCplex(*env_);
  cplex_->extract(*model_);

  allocateVariables(inst, reformulate);

  std::cout << vars_ << std::endl;

  populateByRow(inst, reformulate, add_symmetry_breaking, export_model);
}

void ItemSequencingModel::solve(const Instance &inst, double time_limit)
{
  Solution<double> solution;
  optimize(inst, time_limit, false, false, nullptr, nullptr, solution);

  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  if (is_relaxed_)
  {
    std::cout << "LP: " << solution.lp_ << std::endl;
  }

  if (solution.is_feasible_)
  {

    int num_vehicles = inst.num_vehicles();
    int num_items = inst.num_items();
    const auto &successors_for_transport_per_item = inst.successors_for_transport_per_item();

    IloNumArray x_values(*env_), z_values(*env_), u_values(*env_), U_values(*env_);
    cplex_->getValues(x_values, vars_.x_);
    cplex_->getValues(u_values, vars_.u_);
    cplex_->getValues(U_values, vars_.U_);

    // get values of x variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (auto j : inst.items_for_transport())
      {
        auto value = x_values[vars_.x_index(i, j)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "x[" << i << "," << j << "] = " << value << std::endl;
          num_items_loaded++;
        }
      }
    }

    // get values of u variables.
    for (auto j : inst.items_for_transport())
    {
      for (auto k : inst.items_for_transport())
      {
        if (j != k)
        {
          auto value = u_values[vars_.u_index(j, k)];
          if (!double_equals(value, 0.0))
            std::cout << "u[" << j << "," << k << "] = " << value << std::endl;
        }
      }
    }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = U_values[vars_.U_index(i)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    if (reformulated_)
    {
      IloNumArray X_values(*env_);
      cplex_->getValues(X_values, vars_.X_);

      // get values of X variables.
      for (int i = 0; i < num_vehicles; ++i)
      {
        for (auto &[group, _] : inst.items_for_transport_per_group())
        {
          auto value = X_values[vars_.X_index(i, group)];
          if (!double_equals(value, 0.0))
            std::cout << "X[" << i << "," << group << "] = " << value << std::endl;
        }
      }
      X_values.end();
    }

    std::cout << "num_items_loaded: " << num_items_loaded << std::endl
              << "num_unproductive_moves: " << num_unproductive_moves << std::endl;
    x_values.end();
    z_values.end();
    u_values.end();
    U_values.end();
  }
}

void ItemSequencingModel::populateByRow(const Instance &instance, bool reformulate, bool add_symmetry_breaking, bool export_model)
{
  IloExpr obj(*env_);
  const int num_vehicles = instance.num_vehicles();
  const auto &items = instance.items();
  const int num_items = instance.num_items();
  // const auto &successors_fixed_per_item = instance.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();
  const auto &fleet = instance.fleet();
  const auto &items_for_transport = instance.items_for_transport();
  const int num_items_for_transport = items_for_transport.size();
  int M = 0;

  // compute M.
  for (int j = 0; j < num_items; ++j)
    if (!successors_for_transport_per_item[j].empty())
      ++M;

  // add objective function.
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (item->available_for_transport())
    {
      for (int i = 0; i < num_vehicles; ++i)
      {
        obj += operator*(M, vars_.x(i, j));
      }
    }

    // even if not available for transport.
    if (!successors_for_transport_per_item[j].empty())
      obj -= vars_.U(j);
  }

  model_->add(IloMaximize(*env_, obj));
  obj.end();

  // Restrictions (6).
  for (auto &j : instance.items_for_transport())
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.x(i, j);
    model_->add(exp <= 1);
    exp.end();
  }

  if (reformulate)
  {
    // Constraints (25) and (26).
    for (int i = 0; i < num_vehicles; ++i)
    {
      IloExpr exp(*env_);
      for (auto &[group, items] : instance.items_for_transport_per_group())
      {
        exp += vars_.X(i, group);
        IloExpr sum_items(*env_);

        for (auto j : items)
          sum_items += vars_.x(i, j);

        model_->add(operator*((int)items.size(), vars_.X(i, group)) >= sum_items);
        sum_items.end();
      }

      model_->add(exp <= 1);
      exp.end();
    }
  }
  else
  {
    // Restrictions (7).
    for (auto &j : instance.items_for_transport())
    {
      for (auto &k : instance.items_for_transport())
      {
        if ((j < k) && (items[j]->group() != items[k]->group()))
        {
          for (int i = 0; i < num_vehicles; ++i)
          {
            IloExpr exp(*env_);
            model_->add(vars_.x(i, j) + vars_.x(i, k) <= 1);
            exp.end();
          }
        }
      }
    }
  }

  // Restrictions (8).
  for (int i = 0; i < num_vehicles; ++i)
  {
    IloExpr exp(*env_);
    for (auto &j : instance.items_for_transport())
      exp += vars_.x(i, j);
    model_->add(exp <= fleet[i]->capacity());
    exp.end();
  }

  // Restrictions (9) and (28).
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (!(successors_for_transport_per_item[j].empty()))
    {
      IloExpr exp(*env_);
      if (item->available_for_transport())
      {
        for (auto &k : successors_for_transport_per_item[j])
          exp += vars_.u(k, j);
      }
      else
      {
        for (auto &k : successors_for_transport_per_item[j])
          for (int i = 0; i < num_vehicles; ++i)
            exp += vars_.x(i, k);
      }
      model_->add(operator*((int)(successors_for_transport_per_item[j]).size(), vars_.U(j)) >= exp);
      exp.end();
    }
  }

  // Restrictions (29).
  for (int j = 0; j < num_items_for_transport; ++j)
  {
    int item_pos_j = items_for_transport[j];
    for (int k = j + 1; k < num_items_for_transport; ++k)
    {
      int item_pos_k = items_for_transport[k];
      model_->add(vars_.u(item_pos_j, item_pos_k) + vars_.u(item_pos_k, item_pos_j) <= 1);
    }
  }

  // Restrictions (30) == (35) and (36).
  for (int j = 0; j < num_items_for_transport; ++j)
  {
    int item_pos_j = items_for_transport[j];
    for (int k = j + 1; k < num_items_for_transport; ++k)
    {
      int item_pos_k = items_for_transport[k];
      for (int l = k + 1; l < num_items_for_transport; ++l)
      {
        int item_pos_l = items_for_transport[l];
        model_->add(vars_.u(item_pos_j, item_pos_k) + vars_.u(item_pos_k, item_pos_l) + vars_.u(item_pos_l, item_pos_j) <= 2);
        model_->add(vars_.u(item_pos_j, item_pos_l) + vars_.u(item_pos_l, item_pos_k) + vars_.u(item_pos_k, item_pos_j) <= 2);
      }
    }
  }

  // Restrictions (31)-(35).
  for (int j = 0; j < num_items_for_transport; ++j)
  {
    int item_pos_j = items_for_transport[j];
    for (int k = j + 1; k < num_items_for_transport; ++k)
    {
      int item_pos_k = items_for_transport[k];
      IloExpr sum_x_i_j(*env_), sum_x_i_k(*env_);
      for (int i = 0; i < num_vehicles; ++i)
      {
        sum_x_i_j += vars_.x(i, item_pos_j);
        sum_x_i_k += vars_.x(i, item_pos_k);

        // (31) and (32).
        model_->add(vars_.x(i, item_pos_j) <= vars_.x(i, item_pos_k) + vars_.u(item_pos_j, item_pos_k) + vars_.u(item_pos_k, item_pos_j));
        model_->add(vars_.x(i, item_pos_k) <= vars_.x(i, item_pos_j) + vars_.u(item_pos_j, item_pos_k) + vars_.u(item_pos_k, item_pos_j));

        // (33) and (34).
        model_->add(vars_.u(item_pos_j, item_pos_k) <= vars_.u(item_pos_k, item_pos_j) + 2 - vars_.x(i, item_pos_j) - vars_.x(i, item_pos_k));
        model_->add(vars_.u(item_pos_k, item_pos_j) <= vars_.u(item_pos_j, item_pos_k) + 2 - vars_.x(i, item_pos_j) - vars_.x(i, item_pos_k));
      }

      // (35).
      model_->add(vars_.u(item_pos_k, item_pos_j) >= sum_x_i_k - sum_x_i_j);
      sum_x_i_j.end();
      sum_x_i_k.end();
    }
  }

  for (auto &[group_q, items_q] : instance.items_for_transport_per_group())
  {
    for (auto j : items_q)
    {
      for (auto k : items_q)
      {
        if (j < k)
        {
          for (auto &[group_v, items_v] : instance.items_for_transport_per_group())
          {
            // if (group_v != group_q)
            {
              for (auto l : items_v)
              {
                if ((l != j) && (l != k))
                {
                  model_->add(vars_.u(j, l) <= vars_.u(k, l) + vars_.u(j, k) + vars_.u(k, j));
                  model_->add(vars_.u(k, l) <= vars_.u(j, l) + vars_.u(j, k) + vars_.u(k, j));

                  model_->add(vars_.u(l, j) <= vars_.u(l, k) + vars_.u(j, k) + vars_.u(k, j));
                  model_->add(vars_.u(l, k) <= vars_.u(l, j) + vars_.u(j, k) + vars_.u(k, j));
                }
              }
            }
          }
        }
      }
    }
  }

  if (add_symmetry_breaking)
  {
    // Inequalities (21)-(23).
    for (int j = 0; j < num_items_for_transport; ++j)
    {
      int item_pos_j = items_for_transport[j];
      for (int k = j + 1; k < num_items_for_transport; ++k)
      {
        int item_pos_k = items_for_transport[k];
        IloExpr sum_x_i_j(*env_), sum_x_i_k(*env_);
        for (int i = 0; i < num_vehicles; ++i)
        {
          sum_x_i_j += vars_.x(i, item_pos_j);
          sum_x_i_k += vars_.x(i, item_pos_k);
        }
        model_->add(vars_.u(item_pos_j, item_pos_k) >= 1 - sum_x_i_j - sum_x_i_k);
        model_->add(vars_.u(item_pos_j, item_pos_k) >= sum_x_i_j - sum_x_i_k);
        model_->add(vars_.u(item_pos_k, item_pos_j) >= sum_x_i_k - sum_x_i_j);

        sum_x_i_j.end();
        sum_x_i_k.end();
      }
    }
  }

  if (export_model)
  {
    // add name to variables.
    for (auto &[key, value] : vars_.x_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "x(%d)(%d)", value.first, value.second);
      vars_.x_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.u_var_to_index_.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "u(%d)(%d)", value.first, value.second);
      vars_.u_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.U_var_to_index_.index_to_var_map())
    {
      char strnum[15];
      sprintf(strnum, "U(%d)", value.second);
      vars_.U_[key].setName(strnum);
    }

    if (reformulate)
    {
      for (auto &[key, value] : vars_.X_var_to_index_.index_to_var_map())
      {
        char strnum[26];
        sprintf(strnum, "X(%d)(%d)", value.first, value.second);
        vars_.X_[key].setName(strnum);
      }
    }
    cplex_->exportModel("item_sequencing_model.lp");
  }
}

// Vehicle slots model.
std::ostream &
operator<<(std::ostream &out, const VehicleSlotsModelVariables &vars)
{
  out << "y1 vars:" << std::endl
      << vars.y1_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "y2 vars:" << std::endl
      << vars.y2_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "U vars:" << std::endl
      << vars.U_var_to_index_ << std::endl
      << "***************" << std::endl;
  out << "Y vars:" << std::endl
      << vars.Y_var_to_index_ << std::endl
      << "***************" << std::endl;
  return out;
}

void VehicleSlotsModel::addCut(UserCut *curr_cut)
{
  IloExpr exp(*env_);
  int vehicle = -1, item = -1;
  for (std::list<std::pair<int, int>>::iterator it = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it != curr_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
  {
    vehicle = (*it).first;
    item = (*it).second;
    exp += vars_.y1(vehicle, item); // in this model, the cuts are related to the allocation of item to slots, instead of the actual vehicles.
  }

  model_->add(exp <= 1);
  exp.end();
}

void VehicleSlotsModel::allocateVariables(const Instance &instance, bool reformulate, bool disable_all_binary_vars)
{
  const int num_items = instance.num_items();
  const int num_items_for_transport = instance.num_items_for_transport();
  const int num_vehicles = instance.num_vehicles();
  const int num_slots = num_vehicles;
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();

  double binary_upper_bound = disable_all_binary_vars ? 0.0 : 1.0;
  auto type = is_relaxed_ ? ILOFLOAT : ILOINT;

  // fill y1 variables and index map.
  for (int i = 0; i < num_slots; ++i)
    for (auto j : instance.items_for_transport())
      vars_.y1_var_to_index_.addEntry(i, j);

  vars_.y1_ = IloNumVarArray(*env_, vars_.y1_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill y2 variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (int h = 0; h < num_slots; ++h)
      vars_.y2_var_to_index_.addEntry(i, h);

  vars_.y2_ = IloNumVarArray(*env_, vars_.y2_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill U variables and index map.
  for (int i = 0; i < num_items; ++i)
    if (!(successors_for_transport_per_item[i].empty()))
      vars_.U_var_to_index_.addEntry(0, i);

  vars_.U_ = IloNumVarArray(*env_, vars_.U_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  if (reformulate)
  {
    // fill X variables and index map.
    for (int i = 0; i < num_slots; ++i)
      for (auto &[group, _] : instance.items_for_transport_per_group())
        vars_.Y_var_to_index_.addEntry(i, group);

    vars_.Y_ = IloNumVarArray(*env_, vars_.Y_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
  }
}

VehicleSlotsModel::VehicleSlotsModel(Instance &inst, bool reformulate, bool add_symmetry_breaking, bool solve_relax, bool export_model)
{
  is_relaxed_ = solve_relax;
  reformulated_ = reformulate;

  const int num_vertices = inst.num_items();

  env_ = new IloEnv();
  model_ = new IloModel(*env_);
  cplex_ = new IloCplex(*env_);

  cplex_->extract(*model_);

  allocateVariables(inst, reformulate);

  std::cout << vars_ << std::endl;

  populateByRow(inst, reformulate, add_symmetry_breaking, export_model);
}

void VehicleSlotsModel::solve(const Instance &inst, double time_limit)
{
  Solution<double> solution;
  optimize(inst, time_limit, false, false, nullptr, nullptr, solution);

  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  if (is_relaxed_)
  {
    std::cout << "LP: " << solution.lp_ << std::endl;
  }

  if (solution.is_feasible_)
  {

    int num_vehicles = inst.num_vehicles();
    int num_slots = num_vehicles;
    int num_items = inst.num_items();
    const auto &successors_for_transport_per_item = inst.successors_for_transport_per_item();

    IloNumArray y1_values(*env_), y2_values(*env_), U_values(*env_);
    cplex_->getValues(y1_values, vars_.y1_);
    cplex_->getValues(y2_values, vars_.y2_);
    cplex_->getValues(U_values, vars_.U_);

    // get values of x variables.
    for (int i = 0; i < num_slots; ++i)
    {
      for (auto j : inst.items_for_transport())
      {
        auto value = y1_values[vars_.y1_index(i, j)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "y1[" << i << "," << j << "] = " << value << std::endl;
          num_items_loaded++;
        }
      }
    }

    // get values of y_2 variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (int h = 0; h < num_slots; ++h)
      {
        auto value = y2_values[vars_.y2_index(i, h)];
        if (!double_equals(value, 0.0))
          std::cout << "y2[" << i << "," << h << "] = " << value << std::endl;
      }
    }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = U_values[vars_.U_index(i)];
        if (!double_equals(value, 0.0))
        {
          std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    if (reformulated_)
    {
      IloNumArray Y_values(*env_);
      cplex_->getValues(Y_values, vars_.Y_);

      // get values of X variables.
      for (int i = 0; i < num_slots; ++i)
      {
        for (auto &[group, _] : inst.items_for_transport_per_group())
        {
          auto value = Y_values[vars_.Y_index(i, group)];
          if (!double_equals(value, 0.0))
            std::cout << "X[" << i << "," << group << "] = " << value << std::endl;
        }
      }
      Y_values.end();
    }

    std::cout << "num_items_loaded: " << num_items_loaded << std::endl
              << "num_unproductive_moves: " << num_unproductive_moves << std::endl;
    y1_values.end();
    y2_values.end();
    U_values.end();
  }
}

void VehicleSlotsModel::populateByRow(const Instance &instance, bool reformulate, bool add_symmetry_breaking, bool export_model)
{
  IloExpr obj(*env_);
  const int num_vehicles = instance.num_vehicles();
  const int num_slots = num_vehicles;
  const auto &items = instance.items();
  const int num_items = instance.num_items();
  // const auto &successors_fixed_per_item = instance.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance.successors_for_transport_per_item();
  const auto &fleet = instance.fleet();
  int M = 0;

  // compute M.
  for (int j = 0; j < num_items; ++j)
    if (!successors_for_transport_per_item[j].empty())
      ++M;

  // add objective function.
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    if (item->available_for_transport())
    {
      for (int i = 0; i < num_slots; ++i)
      {
        obj += operator*(M, vars_.y1(i, j));
      }
    }

    // even if not available for trasnport.
    if (!successors_for_transport_per_item[j].empty())
      obj -= vars_.U(j);
  }

  // only in case we wanted to also minimize the number of vehicles used.
  // for (int i = 0; i < num_vehicles; ++i)
  //   for (int r = 0; r < num_slots; ++r)
  //     obj -= vars_.y2(i, r);

  model_->add(IloMaximize(*env_, obj));
  obj.end();

  // Restrictions (47).
  for (auto &j : instance.items_for_transport())
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.y1(i, j);
    model_->add(exp <= 1);
    exp.end();
  }

  if (reformulate)
  {
    // restrictions (58) and (59).
    for (int i = 0; i < num_slots; ++i)
    {
      IloExpr exp(*env_);
      for (auto &[group, items] : instance.items_for_transport_per_group())
      {
        exp += vars_.Y(i, group);
        IloExpr sum_items(*env_);

        for (auto j : items)
          sum_items += vars_.y1(i, j);

        model_->add(operator*((int)items.size(), vars_.Y(i, group)) >= sum_items);
        sum_items.end();
      }

      model_->add(exp <= 1);
      exp.end();
    }
  }
  else
  {
    // Restrictions (48).
    for (auto &j : instance.items_for_transport())
    {
      for (auto &k : instance.items_for_transport())
      {
        if ((j < k) && (items[j]->group() != items[k]->group()))
        {
          for (int i = 0; i < num_slots; ++i)
          {
            IloExpr exp(*env_);
            model_->add(vars_.y1(i, j) + vars_.y1(i, k) <= 1);
            exp.end();
          }
        }
      }
    }
  }

  // Restrictions (49) and (50).
  for (int j = 0; j < num_items; ++j)
  {
    auto &item = items[j];
    for (auto &k : successors_for_transport_per_item[j])
    {
      if (item->available_for_transport())
      {
        for (int r = 0; r < num_slots; ++r)
        {
          IloExpr exp(*env_);
          for (int i = 0; i <= r; ++i)
            exp += vars_.y1(i, j);

          // (50).
          model_->add(vars_.y1(r, k) <= exp + vars_.U(j));
          exp.end();
        }
      }
      else
      {
        IloExpr exp(*env_);
        for (int r = 0; r < num_slots; ++r)
          exp += vars_.y1(r, k);

        // (49).
        model_->add(vars_.U(j) >= exp);
        exp.end();
      }
    }
  }

  // Restrictions (51).
  for (int i = 0; i < num_vehicles; ++i)
  {
    IloExpr exp(*env_);
    for (int r = 0; r < num_slots; ++r)
      exp += vars_.y2(i, r);
    model_->add(exp <= 1);
    exp.end();
  }

  // Restrictions (52).
  for (int r = 0; r < num_slots; ++r)
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.y2(i, r);
    model_->add(exp <= 1);
    exp.end();
  }

  // Restrictions (53).
  for (int r = 0; r < num_slots; ++r)
  {
    IloExpr exp1(*env_);
    for (auto &j : instance.items_for_transport())
      exp1 += vars_.y1(r, j);

    IloExpr exp2(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp2 += operator*(fleet[i]->capacity(), vars_.y2(i, r));

    model_->add(exp1 <= exp2);
    exp1.end();
    exp2.end();
  }

  if (add_symmetry_breaking)
  {
    // Inequalities (56).

    int max_capacity = 0;
    for (int i = 0; i < num_vehicles; ++i)
      max_capacity = std::max(max_capacity, fleet[i]->capacity());

    for (int r = 0; r < num_slots - 1; ++r)
    {
      IloExpr exp1(*env_), exp2(*env_);
      for (auto &j : instance.items_for_transport())
      {
        exp1 += vars_.y1(r, j);
        exp2 += vars_.y1(r + 1, j);
      }
      model_->add(operator*(max_capacity, exp1) >= exp2);

      exp1.end();
      exp2.end();
    }
  }

  if (export_model)
  {
    // add name to variables.
    for (auto &[key, value] : vars_.y1_var_to_index_.index_to_var_map())
    {
      char strnum[29];
      sprintf(strnum, "y1(%d)(%d)", value.first, value.second);
      vars_.y1_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.y2_var_to_index_.index_to_var_map())
    {
      char strnum[29];
      sprintf(strnum, "y2(%d)(%d)", value.first, value.second);
      vars_.y2_[key].setName(strnum);
    }

    for (auto &[key, value] : vars_.U_var_to_index_.index_to_var_map())
    {
      char strnum[15];
      sprintf(strnum, "U(%d)", value.second);
      vars_.U_[key].setName(strnum);
    }

    if (reformulate)
    {
      for (auto &[key, value] : vars_.Y_var_to_index_.index_to_var_map())
      {
        char strnum[26];
        sprintf(strnum, "Y(%d)(%d)", value.first, value.second);
        vars_.Y_[key].setName(strnum);
      }
    }
    cplex_->exportModel("vehicle_slots_model.lp");
  }
}
