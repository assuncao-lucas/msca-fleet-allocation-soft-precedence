#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/exact/formulations.h"

ILOSTLBEGIN

static UserCut *GenerateCliqueConflictCuts(const Instance &instance, const VarToIndexMap &var_to_index_map, const IloNumVarArray &item_to_vehicle_or_slot_vars, const IloNumArray &item_to_vehicle_or_slot_values, bool root_cuts, Solution<double> &sol, std::list<UserCut *> &cuts)
{
  UserCut *best_cut = nullptr, *curr_cut = nullptr;
  const int num_vehicles = instance.num_vehicles();
  const int cut_dimension = item_to_vehicle_or_slot_vars.getSize();
  double curr_cut_violation = 0.0;
  double max_value_in_group_per_vehicle_or_slot = -1.0;
  double curr_value = 0.0;
  int item_max_value_index = -1;
  int var_index = -1;
  double sum_max_values_vehicle = 0.0;

  for (int i = 0; i < num_vehicles; ++i)
  {
    sum_max_values_vehicle = 0.0;
    curr_cut_violation = 0.0;
    curr_cut = new UserCut(cut_dimension, curr_cut_violation, K_TYPE_CLIQUE_CONFLICT_CUT);

    for (auto &[group, items] : instance.items_for_transport_per_group())
    {
      max_value_in_group_per_vehicle_or_slot = -1.0;
      item_max_value_index = -1;
      for (auto &j : items)
      {
        var_index = var_to_index_map.varToIndex(i, j);
        curr_value = item_to_vehicle_or_slot_values[var_index];
        if (double_greater(curr_value, max_value_in_group_per_vehicle_or_slot))
        {
          max_value_in_group_per_vehicle_or_slot = curr_value;
          item_max_value_index = j;
        }
      }

      sum_max_values_vehicle += max_value_in_group_per_vehicle_or_slot;
      curr_cut->AddLhsElement(i, item_max_value_index, var_to_index_map.varToIndex(i, item_max_value_index));
    }

    curr_cut_violation = sum_max_values_vehicle - 1.0;
    bool cut_violated = double_greater(curr_cut_violation, 0.0) ? true : false;

    if (cut_violated)
    {
      // std::cout << "curr_cut_violation: " << curr_cut_violation << std::endl;
      curr_cut->set_curr_abs_violation(curr_cut_violation);

      cuts.push_back(curr_cut);
      sol.set_cut_found(K_TYPE_CLIQUE_CONFLICT_CUT, root_cuts);
      curr_cut->UpdateMeasures();

      if (curr_cut->isBetterThan(best_cut))
        best_cut = curr_cut;
    }
    else
    {
      delete curr_cut;
      curr_cut = nullptr;
    }
  }

  return best_cut;
}

// Model.

Model::Model(const Instance &instance) : instance_(instance)
{
}

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

void Model::set_multithreading(bool multithreading)
{
  if (multithreading)
    cplex_->setParam(IloCplex::Param::Threads, cplex_->getNumCores());
  else
    cplex_->setParam(IloCplex::Param::Threads, 1);
}

bool Model::optimize(double total_time_limit, bool find_root_cuts, std::list<UserCut *> *initial_cuts, std::list<UserCut *> *root_cuts, Solution<double> &solution)
{
  Timestamp *ti = NewTimestamp(), *tf = NewTimestamp();
  Timer *timer = GetTimer();
  timer->Clock(ti);
  bool found_cuts = false;

  cplex_->setParam(IloCplex::Param::ClockType, 2);
  cplex_->setParam(IloCplex::Param::WorkMem, 100000);
  std::cout << "limit of memory 100000MB" << std::endl;
  cplex_->setParam(IloCplex::IloCplex::Param::MIP::Strategy::File, 3);
  // cplex_->setOut(env_->getNullStream());

  addInitialCuts(initial_cuts, solution);

  // if solving relaxed problem, force dual cplex method (multithreading showed to be slower! Dual cplex forces single threading and is also useful when
  // the main LP is solved several times, either on Benders or while separating valid inequalities).
  // if (is_relaxed_)
  // {
  //   // cplex_->setParam(IloCplex::Param::Threads, 1);
  //   // cplex_->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);
  // }

  if (is_relaxed_)
    total_time_limit = -1;

  if (!double_equals(total_time_limit, -1))
  {
    // std::cout << total_time_limit - instance_.time_spent_in_preprocessing() << std::endl;
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
        solution.milp_time_ = timer->ElapsedTime(ti, tf);
      if ((cplex_->getCplexStatus() == IloCplex::Infeasible) || (cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
        solution.is_feasible_ = false;
      // SetSolutionStatus(cplex,solution,solve_relax);

      delete (ti);
      ti = nullptr;
      delete (tf);
      tf = nullptr;
      return false;
    }
    // getchar(); getchar();
    // cont++;
    // std::cout << cont << std::endl;
    found_cuts = false;

    if (is_relaxed_)
    {
      curr_bound = cplex_->getObjValue();

      if (find_root_cuts && double_less(curr_bound, previous_bound, K_TAILING_OFF_TOLERANCE))
        found_cuts |= findAndAddValidInequalities(solution, root_cuts);

      // getchar();
      // getchar();
    }
  } while (found_cuts);

  // if ((cplex_->getCplexStatus() == IloCplex::Infeasible) || (cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
  //   solution.is_feasible_ = false;
  // else

  timer->Clock(tf);
  if (is_relaxed_)
    solution.root_time_ = timer->ElapsedTime(ti, tf);
  else
    solution.milp_time_ = timer->ElapsedTime(ti, tf);

  SetSolutionStatus(*cplex_, solution, is_relaxed_);

  delete (ti);
  ti = nullptr;
  delete (tf);
  tf = nullptr;

  return true;
}

void Model::addInitialCuts(std::list<UserCut *> *initial_cuts, Solution<double> &solution)
{
  if ((initial_cuts != nullptr) && (!(initial_cuts->empty())))
  {
    IloRangeArray root_cuts(*env_);
    for (std::list<UserCut *>::iterator it = initial_cuts->begin(); it != initial_cuts->end(); it++)
    {
      UserCut *curr_user_cut = static_cast<UserCut *>((*it));
      addCut(curr_user_cut, std::reference_wrapper<IloRangeArray>(root_cuts));
      solution.set_cut_added(K_TYPE_CLIQUE_CONFLICT_CUT, false);

      // delete (*it);
      //*it = NULL;
    }

    if (is_relaxed_)
    {
      model_->add(root_cuts); // uses reference to root_cuts.
    }
    else
    {
      cplex_->addUserCuts(root_cuts); // makes copy of root_cuts, then, can delete each element.
      root_cuts.endElements();
    }

    root_cuts.end();
  }
}

bool Model::findAndAddValidInequalities(Solution<double> &sol, std::list<UserCut *> *root_cuts)
{
  // std::cout << " ***************** looking for cuts" << std::endl;
  bool found_cut = false;
  (sol.num_calls_to_callback_lp_) += 1;

  std::list<UserCut *> cuts;

  UserCut *best_cut = nullptr, *curr_cut = nullptr, *local_best_cut = nullptr;

  best_cut = nullptr;
  // find clique conflict cuts.
  local_best_cut = separateCuts(sol, cuts);

  if ((local_best_cut != nullptr) && (local_best_cut->isBetterThan(best_cut)))
    best_cut = local_best_cut;

  if (best_cut != nullptr)
  {
    found_cut = true;

    // adds best cut (most violated).
    addCut(best_cut, nullopt);
    sol.set_cut_added(K_TYPE_CLIQUE_CONFLICT_CUT, is_relaxed_);
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
        addCut(curr_cut, nullopt);

        sol.set_cut_added(K_TYPE_CLIQUE_CONFLICT_CUT, is_relaxed_);
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

void VehicleSequencingModel::addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts)
{
  IloExpr exp(*env_);
  int vehicle = -1, item = -1;
  for (std::list<std::pair<int, int>>::iterator it = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it != curr_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
  {
    vehicle = (*it).first;
    item = (*it).second;
    exp += vars_.x(vehicle, item);
  }

  if (root_cuts.has_value())
    root_cuts->get().add(exp <= 1);
  else
    model_->add(exp <= 1);

  exp.end();
  // cplex_->exportModel("vehicle_sequencing_model.lp");
}

void VehicleSequencingModel::allocateVariables()
{
  const int num_items = instance_.num_items();
  const int num_items_for_transport = instance_.num_items_for_transport();
  const int num_vehicles = instance_.num_vehicles();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

  double binary_upper_bound = 1.0;
  auto type = is_relaxed_ ? ILOFLOAT : ILOINT;

  // fill x variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      vars_.x_var_to_index_.addEntry(i, j);

  vars_.x_ = IloNumVarArray(*env_, vars_.x_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill z variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (int h = 0; h < num_vehicles; ++h)
      if (i != h)
        vars_.z_var_to_index_.addEntry(i, h);

  vars_.z_ = IloNumVarArray(*env_, vars_.z_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill w variables and index map.
  for (auto i : instance_.items_for_transport())
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

  if (reformulated_)
  {
    // fill X variables and index map.
    for (int i = 0; i < num_vehicles; ++i)
      for (auto &[group, _] : instance_.items_for_transport_per_group())
        vars_.X_var_to_index_.addEntry(i, group);

    vars_.X_ = IloNumVarArray(*env_, vars_.X_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
  }
  else
    vars_.X_ = IloNumVarArray(*env_); // just to initialize and avoid seg fault when getting size of X_.

  all_vars_array_ = IloNumVarArray(*env_);
  vars_.fill_all_vars_array(all_vars_array_);
}

void VehicleSequencingModel::updateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active)
{
  const int num_vehicles = instance_.num_vehicles();
  const auto successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  const auto predecessors_for_transport_per_item = instance_.predecessors_for_transport_per_item();
  const auto all_items = instance_.items();

  // add new variables to kernel.
  for (size_t item = items_entering.find_first(); item != boost::dynamic_bitset<>::npos; item = items_entering.find_next(item))
  {
    // activate x variables.
    for (int i = 0; i < num_vehicles; ++i)
      vars_.x(i, item).setUB(1.0);

    // all z variables already active (only indexed by vehicles).

    // activate w variables.
    for (auto &k : predecessors_for_transport_per_item[item]) // activate for all precedessors.
      vars_.w(k, item).setUB(1.0);

    // U variables should be always active (it is always possible to move any item to the buffer).
    // if (!(successors_for_transport_per_item[item].empty()))
    //   vars_.U(item).setUB(1.0);

    if (reformulated_)
    {
      // activate X variables.
      for (int i = 0; i < num_vehicles; ++i)
        vars_.X(i, (all_items[item])->group()).setUB(1.0);
    }
  }

  // remove variables leaving kernel.
  for (size_t item = items_leaving.find_first(); item != boost::dynamic_bitset<>::npos; item = items_leaving.find_next(item))
  {
    // deactivate x variables.
    for (int i = 0; i < num_vehicles; ++i)
      vars_.x(i, item).setUB(0.0);

    // all z variables remain active (only indexed by vehicles).

    // deactivate w variables.
    for (auto &k : predecessors_for_transport_per_item[item]) // deactivate for all precedessors.
      vars_.w(k, item).setUB(0.0);

    // U variables should be always active (it is always possible to move any item to the buffer).
    // if (!(successors_for_transport_per_item[item].empty()))
    //   vars_.U(item).setUB(0.0);

    if (reformulated_)
    {
      // deactivate X variables.
      for (int i = 0; i < num_vehicles; ++i)
      {
        const int group = (all_items[item])->group();
        bool all_items_of_group_inactive = true;

        // only deactivate if all items of group are currently inactive.
        for (auto item_in_group : instance_.items_for_transport_per_group().at(group))
        {
          if (items_active[item_in_group] == 1)
          {
            all_items_of_group_inactive = false;
            break;
          }
        }

        if (all_items_of_group_inactive)
          vars_.X(i, group).setUB(0.0);
      }
    }
  }
}

void VehicleSequencingModel::getSolutionItems(boost::dynamic_bitset<> &curr_items)
{
  curr_items.reset();
  IloNumArray x_values(*env_);
  cplex_->getValues(x_values, vars_.x_);
  const int num_vehicles = instance_.num_vehicles();

  for (auto j : instance_.items_for_transport())
  {
    for (int i = 0; i < num_vehicles; ++i)
    {
      auto index = vars_.x_index(i, j);
      if (double_equals(x_values[index], 1.0))
      {
        curr_items[j] = 1;
        break;
      }
    }
  }

  x_values.end();
}

VehicleSequencingModel::VehicleSequencingModel(const Instance &inst, bool reformulate, bool symmetry_breaking, bool relaxed) : Model(inst)
{
  is_relaxed_ = relaxed;
  reformulated_ = reformulate;
  symmetry_breaking_ = symmetry_breaking;

  env_ = new IloEnv();
  model_ = new IloModel(*env_);
  cplex_ = new IloCplex(*env_);

  cplex_->extract(*model_);

  allocateVariables();

  // std::cout << vars_ << std::endl;

  populateByRow(symmetry_breaking);
}

void VehicleSequencingModel::populateByRow(bool symmetry_breaking)
{
  IloExpr obj(*env_);
  const int num_vehicles = instance_.num_vehicles();
  const auto &items = instance_.items();
  const int num_items = instance_.num_items();
  // const auto &successors_fixed_per_item = instance_.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  const auto &fleet = instance_.fleet();
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
  for (auto &j : instance_.items_for_transport())
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.x(i, j);
    model_->add(exp <= 1);
    exp.end();
  }

  if (reformulated_)
  {
    // Constraints (22) and (23).
    for (int i = 0; i < num_vehicles; ++i)
    {
      IloExpr exp(*env_);
      for (auto &[group, items] : instance_.items_for_transport_per_group())
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
    for (auto &j : instance_.items_for_transport())
    {
      for (auto &k : instance_.items_for_transport())
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
    for (auto &j : instance_.items_for_transport())
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
  for (auto &j : instance_.items_for_transport())
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

  if (reformulated_ && symmetry_breaking) // in this case, symmetry breaking constraint rely on X vars_.
  {
    // Inequalities (24)-(26).
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (int h = i + 1; h < num_vehicles; ++h)
      {
        IloExpr sum_X_i_q(*env_), sum_X_h_q(*env_);
        for (auto &[q, _] : instance_.items_for_transport_per_group())
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

  // if (export_model)
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

    if (reformulated_)
    {
      for (auto &[key, value] : vars_.X_var_to_index_.index_to_var_map())
      {
        char strnum[26];
        sprintf(strnum, "X(%d)(%d)", value.first, value.second);
        vars_.X_[key].setName(strnum);
      }
    }
    // cplex_->exportModel("vehicle_sequencing_model.lp");
  }
}

void VehicleSequencingModel::fillSolution(Solution<double> &solution, std::optional<IloNumArray> var_values_opt)
{
  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  // if (is_relaxed_)
  // {
  //   std::cout << "LP: " << solution.lp_ << std::endl;
  // }

  if (solution.is_feasible_)
  {

    IloNumArray values(*env_);

    if (var_values_opt.has_value())
    {
      values = var_values_opt.value();
    }
    else
    {
      cplex_->getValues(values, all_vars_array_);
    }

    int num_vehicles = instance_.num_vehicles();
    int num_items = instance_.num_items();
    const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

    // IloNumArray x_values(*env_), z_values(*env_), w_values(*env_), U_values(*env_);
    // cplex_->getValues(x_values, vars_.x_);
    // cplex_->getValues(z_values, vars_.z_);
    // cplex_->getValues(w_values, vars_.w_);
    // cplex_->getValues(U_values, vars_.U_);

    int x_size = vars_.x_.getSize();
    int z_size = vars_.z_.getSize();
    int w_size = vars_.w_.getSize();
    int U_size = vars_.U_.getSize();

    // get values of x variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (auto j : instance_.items_for_transport())
      {
        auto value = values[vars_.x_index(i, j)];
        if (!double_equals(value, 0.0))
        {
          // std::cout << "x[" << i << "," << j << "] = " << value << std::endl;
          num_items_loaded++;
        }
      }
    }

    // // get values of z variables.
    // for (int i = 0; i < num_vehicles; ++i)
    // {
    //   for (int h = 0; h < num_vehicles; ++h)
    //   {
    //     if (i != h)
    //     {
    //       auto value = values[x_size + vars_.z_index(i, h)];
    //       // if (!double_equals(value, 0.0))
    //       // std::cout << "z[" << i << "," << h << "] = " << value << std::endl;
    //     }
    //   }
    // }

    // // get values of w variables.
    // for (auto i : instance_.items_for_transport())
    // {
    //   for (auto &k : successors_for_transport_per_item[i])
    //   {
    //     auto value = values[x_size + z_size + vars_.w_index(i, k)];
    //     // if (!double_equals(value, 0.0))
    //     // std::cout << "w[" << i << "," << k << "] = " << value << std::endl;
    //   }
    // }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = values[x_size + z_size + w_size + vars_.U_index(i)];
        if (!double_equals(value, 0.0))
        {
          // std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    // if (reformulated_)
    // {
    //   // IloNumArray X_values(*env_);
    //   // cplex_->getValues(X_values, vars_.X_);

    //   // get values of X variables.
    //   for (int i = 0; i < num_vehicles; ++i)
    //   {
    //     for (auto &[group, _] : instance_.items_for_transport_per_group())
    //     {
    //       auto value = values[x_size + z_size + w_size + U_size + vars_.X_index(i, group)];
    //       // if (!double_equals(value, 0.0))
    //       // std::cout << "X[" << i << "," << group << "] = " << value << std::endl;
    //     }
    //   }

    //   // X_values.end();
    // }

    solution.num_items_loaded_ = num_items_loaded;
    solution.num_unproductive_moves_ = num_unproductive_moves;

    // std::cout << "num_items_loaded: " << num_items_loaded << std::endl
    //           << "num_unproductive_moves: " << num_unproductive_moves << std::endl;
    // x_values.end();
    // z_values.end();
    // w_values.end();
    // U_values.end();

    if (!var_values_opt.has_value())
      values.end();
  }
}

IloNumArray VehicleSequencingModel::get_items_reduced_costs()
{
  IloNumArray x_reduced_costs(*env_);
  IloNumArray items_reduced_costs(*env_, instance_.num_items());
  cplex_->getReducedCosts(x_reduced_costs, vars_.x_);
  int num_vehicles = instance_.num_vehicles();

  // sum x values for each item for transport.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      items_reduced_costs[j] += x_reduced_costs[vars_.x_index(i, j)];

  // // print
  // for (auto j : instance_.items_for_transport())
  // {
  //   for (int i = 0; i < num_vehicles; ++i)
  //   {
  //     std::cout << "x[" << i << "][" << j << "] = " << x_reduced_costs[vars_.x_index(i, j)] << std::endl;
  //   }
  //   std::cout << "sum " << j << " : " << items_reduced_costs[j] << std::endl;
  // }

  return std::move(items_reduced_costs);
}

IloNumArray VehicleSequencingModel::get_items_values()
{
  IloNumArray x_values(*env_);
  IloNumArray items_values(*env_, instance_.num_items());
  cplex_->getValues(x_values, vars_.x_);

  int num_vehicles = instance_.num_vehicles();

  // sum x values for each item for transport.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      items_values[j] += x_values[vars_.x_index(i, j)];

  // // print
  // for (auto j : instance_.items_for_transport())
  // {
  //   for (int i = 0; i < num_vehicles; ++i)
  //   {
  //     std::cout << "x[" << i << "][" << j << "] = " << x_values[vars_.x_index(i, j)] << std::endl;
  //   }
  //   std::cout << "sum " << j << " : " << items_values[j] << std::endl;
  // }

  return std::move(items_values);
}

UserCut *VehicleSequencingModel::separateCuts(Solution<double> &sol, std::list<UserCut *> &cuts)
{
  IloNumArray item_to_vehicle_or_slot_values(*env_);
  cplex_->getValues(vars_.x_, item_to_vehicle_or_slot_values);
  return GenerateCliqueConflictCuts(instance_, vars_.x_var_to_index_, vars_.x_, item_to_vehicle_or_slot_values, is_relaxed_, sol, cuts);
}

std::vector<int> VehicleSequencingModel::fillVarValuesFromSolution(std::vector<std::pair<int, std::vector<int>>> fleet_allocation)
{
  auto items = instance_.items();
  std::vector<int> values(num_vars(), 0);
  // for (int i = 0; i < num_vars(); ++i)
  //   values[i] = 0;

  int x_size = vars_.x_.getSize();
  int z_size = vars_.z_.getSize();
  int w_size = vars_.w_.getSize();
  int U_size = vars_.U_.getSize();

  const auto &precedence = instance_.precedence_matrix();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  std::list<int> already_loaded_vehicles;

  std::unordered_set<int> items_already_unloaded; // keeps items already moved to the buffer area or loaded to a vehicle.

  for (int slot = 0; slot < fleet_allocation.size(); ++slot)
  {
    const auto vehicle = fleet_allocation[slot].first;
    const auto &items_of_vehicle = fleet_allocation[slot].second;
    // std::cout << vehicle << ":";
    for (auto item : items_of_vehicle)
    {
      // std::cout << " " << item;
      // std::cout << "x " << vehicle << " " << item << std::endl;
      values[vars_.x_index(vehicle, item)] = 1;
      auto group = items[item]->group();
      if (reformulated_)
      {
        // std::cout << "X " << vehicle << " " << group << std::endl;
        values[x_size + z_size + w_size + U_size + vars_.X_index(vehicle, group)] = 1;
      }

      // look for blocking items that were not already removed.
      for (int j = 0; j < items.size(); ++j)
      {
        if (precedence[j][item]) // j precedes/blocks item and item is available for transport.
        {
          if (items_already_unloaded.find(j) == items_already_unloaded.end()) // if not yet moved to buffer or unloaded.
          {
            if (items[j]->available_for_transport())
            {
              values[x_size + z_size + vars_.w_index(j, item)] = 1;
              // std::cout << "w " << j << " " << item << std::endl;
            }
            // std::cout << "U " << j << std::endl;
            values[x_size + z_size + w_size + vars_.U_index(j)] = 1; // mark precedence violation.
            items_already_unloaded.insert(j);                        // move item to buffer and mark violation at U_j.
          }
        }
      }
      items_already_unloaded.insert(item);
    }

    // std::cout << std::endl;

    // std::cout << "y2 " << vehicle << " " << slot << std::endl;
    // marks vehicle precedence for all vehicles already loaded.
    for (auto previous_vehicle : already_loaded_vehicles)
    {
      // std::cout << "z " << previous_vehicle << " " << vehicle << std::endl;
      values[x_size + vars_.z_index(previous_vehicle, vehicle)] = 1;
    }

    already_loaded_vehicles.push_back(vehicle);
  }

  // // print
  // all_vars_array_ = IloNumVarArray(*env_);
  // vars_.fill_all_vars_array(all_vars_array_);
  // for (int i = 0; i < num_vars(); ++i)
  // {
  //   if (values[i] > 0)
  //     std::cout << all_vars_array_[i].getName() << ": " << values[i] << std::endl;
  // }
  // all_vars_array_.end();

  // getchar();
  // getchar();

  return values;
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

void ItemSequencingModel::addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts)
{
  IloExpr exp(*env_);
  int vehicle = -1, item = -1;
  for (std::list<std::pair<int, int>>::iterator it = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it != curr_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
  {
    vehicle = (*it).first;
    item = (*it).second;
    exp += vars_.x(vehicle, item);
  }

  if (root_cuts.has_value())
    root_cuts->get().add(exp <= 1);
  else
    model_->add(exp <= 1);

  exp.end();
}

IloNumArray ItemSequencingModel::get_items_reduced_costs()
{
  IloNumArray x_reduced_costs(*env_);
  IloNumArray items_reduced_costs(*env_, instance_.num_items());
  cplex_->getReducedCosts(x_reduced_costs, vars_.x_);
  int num_vehicles = instance_.num_vehicles();

  // sum x values for each item for transport.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      items_reduced_costs[j] += x_reduced_costs[vars_.x_index(i, j)];

  // // print
  // for (auto j : instance_.items_for_transport())
  // {
  //   for (int i = 0; i < num_vehicles; ++i)
  //   {
  //     std::cout << "x[" << i << "][" << j << "] = " << x_reduced_costs[vars_.x_index(i, j)] << std::endl;
  //   }
  //   std::cout << "sum " << j << " : " << items_reduced_costs[j] << std::endl;
  // }

  return std::move(items_reduced_costs);
}

IloNumArray ItemSequencingModel::get_items_values()
{
  IloNumArray x_values(*env_);
  IloNumArray items_values(*env_, instance_.num_items());
  cplex_->getValues(x_values, vars_.x_);
  int num_vehicles = instance_.num_vehicles();

  // sum x values for each item for transport.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      items_values[j] += x_values[vars_.x_index(i, j)];

  // // print
  // for (auto j : instance_.items_for_transport())
  // {
  //   for (int i = 0; i < num_vehicles; ++i)
  //   {
  //     std::cout << "x[" << i << "][" << j << "] = " << x_values[vars_.x_index(i, j)] << std::endl;
  //   }
  //   std::cout << "sum " << j << " : " << items_values[j] << std::endl;
  // }

  return std::move(items_values);
}

void ItemSequencingModel::allocateVariables()
{
  const int num_items = instance_.num_items();
  const int num_items_for_transport = instance_.num_items_for_transport();
  const int num_vehicles = instance_.num_vehicles();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

  double binary_upper_bound = 1.0;
  auto type = is_relaxed_ ? ILOFLOAT : ILOINT;

  // fill x variables and index map.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      vars_.x_var_to_index_.addEntry(i, j);

  vars_.x_ = IloNumVarArray(*env_, vars_.x_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill u variables and index map.
  for (auto j : instance_.items_for_transport())
    for (auto k : instance_.items_for_transport())
      if (j != k)
        vars_.u_var_to_index_.addEntry(j, k);

  vars_.u_ = IloNumVarArray(*env_, vars_.u_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  // fill U variables and index map.
  for (int i = 0; i < num_items; ++i)
    if (!(successors_for_transport_per_item[i].empty()))
      vars_.U_var_to_index_.addEntry(0, i);

  vars_.U_ = IloNumVarArray(*env_, vars_.U_var_to_index_.numVars(), 0.0, binary_upper_bound, type);

  if (reformulated_)
  {
    // fill X variables and index map.
    for (int i = 0; i < num_vehicles; ++i)
      for (auto &[group, _] : instance_.items_for_transport_per_group())
        vars_.X_var_to_index_.addEntry(i, group);

    vars_.X_ = IloNumVarArray(*env_, vars_.X_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
  }
  else
    vars_.X_ = IloNumVarArray(*env_); // just to initialize and avoid seg fault when getting size of X_.

  all_vars_array_ = IloNumVarArray(*env_);
  vars_.fill_all_vars_array(all_vars_array_);
}

void ItemSequencingModel::updateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active)
{
  const int num_vehicles = instance_.num_vehicles();
  const auto successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  const auto all_items = instance_.items();

  // add new variables to kernel.
  for (size_t item = items_entering.find_first(); item != boost::dynamic_bitset<>::npos; item = items_entering.find_next(item))
  {
    // activate x variables.
    for (int i = 0; i < num_vehicles; ++i)
      vars_.x(i, item).setUB(1.0);

    // should ALWAYS allow all orderings as to avoid messing up with the symmetry breaking constraints !!!!!!
    // // activate u variables.
    // for (auto k : instance_.items_for_transport())
    // {
    //   if (item != k)
    //   {
    //     // if (items_active[k] == 1)
    //     // { // if both are active, allow all possibilities of ordering.
    //     // vars_.u(item, k).setUB(1.0);
    //     // vars_.u(k, item).setUB(1.0);
    //     // }
    //     // else // only allow the active one coming before the inactive.
    //     // {
    //     //   vars_.u(item, k).setUB(1.0);
    //     // }
    //   }
    // }

    // U variables should be always active (it is always possible to move any item to the buffer).
    // if (!(successors_for_transport_per_item[item].empty()))
    //   vars_.U(item).setUB(1.0);

    if (reformulated_)
    {
      // activate X variables.
      for (int i = 0; i < num_vehicles; ++i)
        vars_.X(i, (all_items[item])->group()).setUB(1.0);
    }
  }

  // remove variables leaving kernel.
  for (size_t item = items_leaving.find_first(); item != boost::dynamic_bitset<>::npos; item = items_leaving.find_next(item))
  {
    // deactivate x variables.
    for (int i = 0; i < num_vehicles; ++i)
      vars_.x(i, item).setUB(0.0);

    // should ALWAYS allow all orderings as to avoid messing up with the symmetry breaking constraints !!!!!!
    // // deactivate u variables.
    // for (auto k : instance_.items_for_transport())
    // {
    //   if (item != k)
    //   {
    //     if (items_active[k] == 0)
    //     { // if both are inactive, fix ordering with both u == 0.
    //       vars_.u(item, k).setUB(0.0);
    //       vars_.u(k, item).setUB(0.0);
    //     }
    //     else // forbid the inactive one coming before the active.
    //     {
    //       vars_.u(item, k).setUB(0.0);
    //     }
    //   }
    // }

    // U variables should be always active (it is always possible to move any item to the buffer).
    // if (!(successors_for_transport_per_item[item].empty()))
    //   vars_.U(item).setUB(0.0);

    if (reformulated_)
    {
      // deactivate X variables.
      for (int i = 0; i < num_vehicles; ++i)
      {
        const int group = (all_items[item])->group();
        bool all_items_of_group_inactive = true;

        // only deactivate if all items of group are currently inactive.
        for (auto item_in_group : instance_.items_for_transport_per_group().at(group))
        {
          if (items_active[item_in_group] == 1)
          {
            all_items_of_group_inactive = false;
            break;
          }
        }

        if (all_items_of_group_inactive)
          vars_.X(i, group).setUB(0.0);
      }
    }
  }

  // cplex_->exportModel("item_sequencing_model.lp");
}

void ItemSequencingModel::getSolutionItems(boost::dynamic_bitset<> &curr_items)
{
  curr_items.reset();
  IloNumArray x_values(*env_);
  cplex_->getValues(x_values, vars_.x_);
  const int num_vehicles = instance_.num_vehicles();

  for (auto j : instance_.items_for_transport())
  {
    for (int i = 0; i < num_vehicles; ++i)
    {
      auto index = vars_.x_index(i, j);
      if (double_equals(x_values[index], 1.0))
      {
        curr_items[j] = 1;
        break;
      }
    }
  }

  x_values.end();
}

ItemSequencingModel::ItemSequencingModel(const Instance &instance, bool reformulate, bool add_symmetry_breaking, bool solve_relax) : Model(instance)
{
  is_relaxed_ = solve_relax;
  reformulated_ = reformulate;
  symmetry_breaking_ = add_symmetry_breaking;

  const int num_vertices = instance_.num_items();

  env_ = new IloEnv();
  model_ = new IloModel(*env_);
  cplex_ = new IloCplex(*env_);
  cplex_->extract(*model_);

  allocateVariables();

  // std::cout << vars_ << std::endl;

  populateByRow(add_symmetry_breaking);
}

void ItemSequencingModel::fillSolution(Solution<double> &solution, std::optional<IloNumArray> var_values_opt)
{
  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  // if (is_relaxed_)
  // {
  //   std::cout << "LP: " << solution.lp_ << std::endl;
  // }

  if (solution.is_feasible_)
  {

    IloNumArray values(*env_);

    if (var_values_opt.has_value())
    {
      values = var_values_opt.value();
    }
    else
    {
      cplex_->getValues(values, all_vars_array_);
    }

    int num_vehicles = instance_.num_vehicles();
    int num_items = instance_.num_items();
    const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

    // IloNumArray x_values(*env_), z_values(*env_), u_values(*env_), U_values(*env_);
    // cplex_->getValues(x_values, vars_.x_);
    // cplex_->getValues(u_values, vars_.u_);
    // cplex_->getValues(U_values, vars_.U_);

    int x_size = vars_.x_.getSize();
    int u_size = vars_.u_.getSize();
    int U_size = vars_.U_.getSize();

    // get values of x variables.
    for (int i = 0; i < num_vehicles; ++i)
    {
      for (auto j : instance_.items_for_transport())
      {
        auto value = values[vars_.x_index(i, j)];
        if (!double_equals(value, 0.0))
        {
          // std::cout << "x[" << i << "," << j << "] = " << value << std::endl;
          num_items_loaded++;
        }
      }
    }

    // // get values of u variables.
    // for (auto j : instance_.items_for_transport())
    // {
    //   for (auto k : instance_.items_for_transport())
    //   {
    //     if (j != k)
    //     {
    //       auto value = values[x_size + vars_.u_index(j, k)];
    //       // if (!double_equals(value, 0.0))
    //       // std::cout << "u[" << j << "," << k << "] = " << value << std::endl;
    //     }
    //   }
    // }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = values[x_size + u_size + vars_.U_index(i)];
        if (!double_equals(value, 0.0))
        {
          // std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    // if (reformulated_)
    // {
    //   // IloNumArray X_values(*env_);
    //   // cplex_->getValues(X_values, vars_.X_);

    //   // get values of X variables.
    //   for (int i = 0; i < num_vehicles; ++i)
    //   {
    //     for (auto &[group, _] : instance_.items_for_transport_per_group())
    //     {
    //       auto value = values[x_size + u_size + U_size + vars_.X_index(i, group)];
    //       // if (!double_equals(value, 0.0))
    //       // std::cout << "X[" << i << "," << group << "] = " << value << std::endl;
    //     }
    //   }

    //   // X_values.end();
    // }

    // std::cout << "num_items_loaded: " << num_items_loaded << std::endl
    //           << "num_unproductive_moves: " << num_unproductive_moves << std::endl;

    solution.num_items_loaded_ = num_items_loaded;
    solution.num_unproductive_moves_ = num_unproductive_moves;

    // x_values.end();
    // z_values.end();
    // u_values.end();
    // U_values.end();
    if (!var_values_opt.has_value())
      values.end();
  }
}

void ItemSequencingModel::populateByRow(bool add_symmetry_breaking)
{
  IloExpr obj(*env_);
  const int num_vehicles = instance_.num_vehicles();
  const auto &items = instance_.items();
  const int num_items = instance_.num_items();
  // const auto &successors_fixed_per_item = instance_.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  const auto &fleet = instance_.fleet();
  const auto &items_for_transport = instance_.items_for_transport();
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
  for (auto &j : instance_.items_for_transport())
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.x(i, j);
    model_->add(exp <= 1);
    exp.end();
  }

  if (reformulated_)
  {
    // Constraints (25) and (26).
    for (int i = 0; i < num_vehicles; ++i)
    {
      IloExpr exp(*env_);
      for (auto &[group, items] : instance_.items_for_transport_per_group())
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
    for (auto &j : instance_.items_for_transport())
    {
      for (auto &k : instance_.items_for_transport())
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
    for (auto &j : instance_.items_for_transport())
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

  for (auto &[group_q, items_q] : instance_.items_for_transport_per_group())
  {
    for (auto j : items_q)
    {
      for (auto k : items_q)
      {
        if (j < k)
        {
          for (auto &[group_v, items_v] : instance_.items_for_transport_per_group())
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

  // if (export_model)
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

    if (reformulated_)
    {
      for (auto &[key, value] : vars_.X_var_to_index_.index_to_var_map())
      {
        char strnum[26];
        sprintf(strnum, "X(%d)(%d)", value.first, value.second);
        vars_.X_[key].setName(strnum);
      }
    }
  }

  // cplex_->exportModel("item_sequencing_model.lp");
}

UserCut *ItemSequencingModel::separateCuts(Solution<double> &sol, std::list<UserCut *> &cuts)
{
  IloNumArray item_to_vehicle_or_slot_values(*env_);
  cplex_->getValues(vars_.x_, item_to_vehicle_or_slot_values);
  return GenerateCliqueConflictCuts(instance_, vars_.x_var_to_index_, vars_.x_, item_to_vehicle_or_slot_values, is_relaxed_, sol, cuts);
}

std::vector<int> ItemSequencingModel::fillVarValuesFromSolution(std::vector<std::pair<int, std::vector<int>>> fleet_allocation)
{
  auto items = instance_.items();
  std::vector<int> values(num_vars(), 0);
  // for (int i = 0; i < num_vars(); ++i)
  //   values[i] = 0;

  int x_size = vars_.x_.getSize();
  int u_size = vars_.u_.getSize();
  int U_size = vars_.U_.getSize();

  const auto &precedence = instance_.precedence_matrix();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

  std::unordered_set<int> items_already_unloaded; // keeps items already moved to the buffer area or loaded to a vehicle.

  for (int slot = 0; slot < fleet_allocation.size(); ++slot)
  {
    const auto vehicle = fleet_allocation[slot].first;
    const auto &items_of_vehicle = fleet_allocation[slot].second;
    // std::cout << vehicle << ":";
    for (auto item : items_of_vehicle)
    {
      // std::cout << " " << item;
      // std::cout << "x " << vehicle << " " << item << std::endl;
      values[vars_.x_index(vehicle, item)] = 1;
      auto group = items[item]->group();
      if (reformulated_)
      {
        // std::cout << "X " << vehicle << " " << group << std::endl;
        values[x_size + u_size + U_size + vars_.X_index(vehicle, group)] = 1;
      }

      // look for blocking items that were not already removed.
      for (int j = 0; j < items.size(); ++j)
      {
        if (precedence[j][item]) // j precedes/blocks item and item is available for transport.
        {
          if (items_already_unloaded.find(j) == items_already_unloaded.end()) // if not yet moved to buffer or unloaded.
          {
            if (items[j]->available_for_transport())
            {
              values[x_size + vars_.u_index(j, item)] = 1;
              // std::cout << "u " << j << " " << item << std::endl;
            }
            // std::cout << "U " << j << std::endl;
            values[x_size + u_size + vars_.U_index(j)] = 1; // mark precedence violation.
            items_already_unloaded.insert(j);               // move item to buffer and mark violation at U_j.
          }
        }
      }
      items_already_unloaded.insert(item);
    }

    // std::cout << std::endl;
  }

  // // print
  // all_vars_array_ = IloNumVarArray(*env_);
  // vars_.fill_all_vars_array(all_vars_array_);
  // for (int i = 0; i < num_vars(); ++i)
  // {
  //   if (values[i] > 0)
  //     std::cout << all_vars_array_[i].getName() << ": " << values[i] << std::endl;
  // }
  // all_vars_array_.end();

  // getchar();
  // getchar();

  return values;
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

void VehicleSlotsModel::addCut(UserCut *curr_cut, std::optional<std::reference_wrapper<IloRangeArray>> root_cuts)
{
  IloExpr exp(*env_);
  int vehicle = -1, item = -1;
  for (std::list<std::pair<int, int>>::iterator it = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it != curr_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
  {
    vehicle = (*it).first;
    item = (*it).second;
    exp += vars_.y1(vehicle, item); // in this model, the cuts are related to the allocation of item to slots, instead of the actual vehicles.
  }

  if (root_cuts.has_value())
    root_cuts->get().add(exp <= 1);
  else
    model_->add(exp <= 1);

  exp.end();
}

IloNumArray VehicleSlotsModel::get_items_reduced_costs()
{
  IloNumArray y1_reduced_costs(*env_);
  IloNumArray items_reduced_costs(*env_, instance_.num_items());
  cplex_->getReducedCosts(y1_reduced_costs, vars_.y1_);
  int num_vehicles = instance_.num_vehicles();

  // sum x values for each item for transport.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      items_reduced_costs[j] += y1_reduced_costs[vars_.y1_index(i, j)];

  // // print
  // for (auto j : instance_.items_for_transport())
  // {
  //   for (int i = 0; i < num_vehicles; ++i)
  //   {
  //     std::cout << "x[" << i << "][" << j << "] = " << y1_reduced_costs[vars_.y1_index(i, j)] << std::endl;
  //   }
  //   std::cout << "sum " << j << " : " << items_reduced_costs[j] << std::endl;
  // }

  return std::move(items_reduced_costs);
}

IloNumArray VehicleSlotsModel::get_items_values()
{
  IloNumArray y1_values(*env_);
  IloNumArray items_values(*env_, instance_.num_items());
  cplex_->getValues(y1_values, vars_.y1_);
  int num_vehicles = instance_.num_vehicles();

  // sum x values for each item for transport.
  for (int i = 0; i < num_vehicles; ++i)
    for (auto j : instance_.items_for_transport())
      items_values[j] += y1_values[vars_.y1_index(i, j)];

  // // print
  // for (auto j : instance_.items_for_transport())
  // {
  //   for (int i = 0; i < num_vehicles; ++i)
  //   {
  //     std::cout << "x[" << i << "][" << j << "] = " << y1_values[vars_.y1_index(i, j)] << std::endl;
  //   }
  //   std::cout << "sum " << j << " : " << items_values[j] << std::endl;
  // }

  return std::move(items_values);
}

void VehicleSlotsModel::allocateVariables()
{
  const int num_items = instance_.num_items();
  const int num_items_for_transport = instance_.num_items_for_transport();
  const int num_vehicles = instance_.num_vehicles();
  const int num_slots = num_vehicles;
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

  double binary_upper_bound = 1.0;
  auto type = is_relaxed_ ? ILOFLOAT : ILOINT;

  // fill y1 variables and index map.
  for (int i = 0; i < num_slots; ++i)
    for (auto j : instance_.items_for_transport())
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

  if (reformulated_)
  {
    // fill Y variables and index map.
    for (int i = 0; i < num_slots; ++i)
      for (auto &[group, _] : instance_.items_for_transport_per_group())
        vars_.Y_var_to_index_.addEntry(i, group);

    vars_.Y_ = IloNumVarArray(*env_, vars_.Y_var_to_index_.numVars(), 0.0, binary_upper_bound, type);
  }
  else
    vars_.Y_ = IloNumVarArray(*env_); // just to initialize and avoid seg fault when getting size of Y_.

  all_vars_array_ = IloNumVarArray(*env_);
  vars_.fill_all_vars_array(all_vars_array_);
}

void VehicleSlotsModel::updateModelVarBounds(boost::dynamic_bitset<> &items_entering, boost::dynamic_bitset<> &items_leaving, boost::dynamic_bitset<> &items_active)
{
  const int num_vehicles = instance_.num_vehicles();
  const int num_slots = num_vehicles;
  const auto successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  const auto all_items = instance_.items();

  // add new variables to kernel.
  for (size_t item = items_entering.find_first(); item != boost::dynamic_bitset<>::npos; item = items_entering.find_next(item))
  {
    // activate y1 variables.
    for (int i = 0; i < num_slots; ++i)
      vars_.y1(i, item).setUB(1.0);

    // all y2 variables already active (only indexed by vehicles/slots).

    // U variables should be always active (it is always possible to move any item to the buffer).
    // if (!(successors_for_transport_per_item[item].empty()))
    //   vars_.U(item).setUB(1.0);

    if (reformulated_)
    {
      // activate Y variables.
      for (int i = 0; i < num_vehicles; ++i)
        vars_.Y(i, (all_items[item])->group()).setUB(1.0);
    }
  }

  // remove variables leaving kernel.
  for (size_t item = items_leaving.find_first(); item != boost::dynamic_bitset<>::npos; item = items_leaving.find_next(item))
  {

    // deactivate y1 variables.
    for (int i = 0; i < num_slots; ++i)
      vars_.y1(i, item).setUB(0.0);

    // all y2 variables remain active (only indexed by vehicles).

    // U variables should be always active (it is always possible to move any item to the buffer).
    // if (!(successors_for_transport_per_item[item].empty()))
    //   vars_.U(item).setUB(0.0);

    if (reformulated_)
    {
      // deactivate Y variables.
      for (int i = 0; i < num_vehicles; ++i)
      {
        const int group = (all_items[item])->group();
        bool all_items_of_group_inactive = true;

        // only deactivate if all items of group are currently inactive.
        for (auto item_in_group : instance_.items_for_transport_per_group().at(group))
        {
          if (items_active[item_in_group] == 1)
          {
            all_items_of_group_inactive = false;
            break;
          }
        }

        if (all_items_of_group_inactive)
          vars_.Y(i, group).setUB(0.0);
      }
    }
  }
}

void VehicleSlotsModel::getSolutionItems(boost::dynamic_bitset<> &curr_items)
{
  curr_items.reset();
  IloNumArray y1_values(*env_);
  cplex_->getValues(y1_values, vars_.y1_);
  const int num_slots = instance_.num_vehicles();

  for (auto j : instance_.items_for_transport())
  {
    for (int i = 0; i < num_slots; ++i)
    {
      auto index = vars_.y1_index(i, j);
      if (double_equals(y1_values[index], 1.0))
      {
        curr_items[j] = 1;
        break;
      }
    }
  }

  y1_values.end();
}

VehicleSlotsModel::VehicleSlotsModel(const Instance &inst, bool reformulate, bool add_symmetry_breaking, bool solve_relax) : Model(inst)
{
  is_relaxed_ = solve_relax;
  reformulated_ = reformulate;
  symmetry_breaking_ = add_symmetry_breaking;

  const int num_vertices = instance_.num_items();

  env_ = new IloEnv();
  model_ = new IloModel(*env_);
  cplex_ = new IloCplex(*env_);

  cplex_->extract(*model_);

  allocateVariables();

  // std::cout << vars_ << std::endl;

  populateByRow(add_symmetry_breaking);
}

void VehicleSlotsModel::fillSolution(Solution<double> &solution, std::optional<IloNumArray> var_values_opt)
{
  int num_items_loaded = 0;
  int num_unproductive_moves = 0;

  // if (is_relaxed_)
  // {
  //   std::cout << "LP: " << solution.lp_ << std::endl;
  // }

  if (solution.is_feasible_)
  {

    IloNumArray values(*env_);

    if (var_values_opt.has_value())
    {
      values = var_values_opt.value();
    }
    else
    {
      cplex_->getValues(values, all_vars_array_);
    }

    int num_vehicles = instance_.num_vehicles();
    int num_slots = num_vehicles;
    int num_items = instance_.num_items();
    const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

    // IloNumArray y1_values(*env_), y2_values(*env_), U_values(*env_);
    // cplex_->getValues(y1_values, vars_.y1_);
    // cplex_->getValues(y2_values, vars_.y2_);
    // cplex_->getValues(U_values, vars_.U_);

    int y1_size = vars_.y1_.getSize();
    int y2_size = vars_.y2_.getSize();
    int U_size = vars_.U_.getSize();

    // get values of x variables.
    for (int i = 0; i < num_slots; ++i)
    {
      for (auto j : instance_.items_for_transport())
      {
        auto value = values[vars_.y1_index(i, j)];
        if (!double_equals(value, 0.0))
        {
          // std::cout << "y1[" << i << "," << j << "] = " << value << std::endl;
          num_items_loaded++;
        }
      }
    }

    // // get values of y_2 variables.
    // for (int i = 0; i < num_vehicles; ++i)
    // {
    //   for (int h = 0; h < num_slots; ++h)
    //   {
    //     auto value = values[y1_size + vars_.y2_index(i, h)];
    //     // if (!double_equals(value, 0.0))
    //     // std::cout << "y2[" << i << "," << h << "] = " << value << std::endl;
    //   }
    // }

    // get values of U variables.
    for (int i = 0; i < num_items; ++i)
    {
      if (!(successors_for_transport_per_item[i].empty()))
      {
        auto value = values[y1_size + y2_size + vars_.U_index(i)];
        if (!double_equals(value, 0.0))
        {
          // std::cout << "U[" << i << "] = " << value << std::endl;
          num_unproductive_moves++;
        }
      }
    }

    // if (reformulated_)
    // {
    //   // IloNumArray Y_values(*env_);
    //   // cplex_->getValues(Y_values, vars_.Y_);

    //   // get values of X variables.
    //   for (int i = 0; i < num_slots; ++i)
    //   {
    //     for (auto &[group, _] : instance_.items_for_transport_per_group())
    //     {
    //       auto value = values[y1_size + y2_size + U_size + vars_.Y_index(i, group)];
    //       // if (!double_equals(value, 0.0))
    //       // std::cout << "X[" << i << "," << group << "] = " << value << std::endl;
    //     }
    //   }
    //   // Y_values.end();
    // }

    // std::cout << "num_items_loaded: " << num_items_loaded << std::endl
    //           << "num_unproductive_moves: " << num_unproductive_moves << std::endl;

    solution.num_items_loaded_ = num_items_loaded;
    solution.num_unproductive_moves_ = num_unproductive_moves;
    // y1_values.end();
    // y2_values.end();
    // U_values.end();

    if (!var_values_opt.has_value())
      values.end();
  }
}

void VehicleSlotsModel::populateByRow(bool add_symmetry_breaking)
{
  IloExpr obj(*env_);
  const int num_vehicles = instance_.num_vehicles();
  const int num_slots = num_vehicles;
  const auto &items = instance_.items();
  const int num_items = instance_.num_items();
  // const auto &successors_fixed_per_item = instance_.successors_fixed_per_item();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();
  const auto &fleet = instance_.fleet();
  int M = 0;

  // compute M.
  for (int j = 0; j < num_items; ++j)
    if (!successors_for_transport_per_item[j].empty())
      ++M;

  // std::cout << "M: " << M << std::endl;
  // getchar();
  // getchar();

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
  for (auto &j : instance_.items_for_transport())
  {
    IloExpr exp(*env_);
    for (int i = 0; i < num_vehicles; ++i)
      exp += vars_.y1(i, j);
    model_->add(exp <= 1);
    exp.end();
  }

  if (reformulated_)
  {
    // restrictions (58) and (59).
    for (int i = 0; i < num_slots; ++i)
    {
      IloExpr exp(*env_);
      for (auto &[group, items] : instance_.items_for_transport_per_group())
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
    for (auto &j : instance_.items_for_transport())
    {
      for (auto &k : instance_.items_for_transport())
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
    for (auto &j : instance_.items_for_transport())
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
      for (auto &j : instance_.items_for_transport())
      {
        exp1 += vars_.y1(r, j);
        exp2 += vars_.y1(r + 1, j);
      }
      model_->add(operator*(max_capacity, exp1) >= exp2);

      exp1.end();
      exp2.end();
    }
  }

  // if (export_model)
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

    if (reformulated_)
    {
      for (auto &[key, value] : vars_.Y_var_to_index_.index_to_var_map())
      {
        char strnum[26];
        sprintf(strnum, "Y(%d)(%d)", value.first, value.second);
        vars_.Y_[key].setName(strnum);
      }
    }
    // cplex_->exportModel("vehicle_slots_model.lp");
  }
}

UserCut *VehicleSlotsModel::separateCuts(Solution<double> &sol, std::list<UserCut *> &cuts)
{
  IloNumArray item_to_vehicle_or_slot_values(*env_);
  cplex_->getValues(vars_.y1_, item_to_vehicle_or_slot_values);
  return GenerateCliqueConflictCuts(instance_, vars_.y1_var_to_index_, vars_.y1_, item_to_vehicle_or_slot_values, is_relaxed_, sol, cuts);
}

std::vector<int> VehicleSlotsModel::fillVarValuesFromSolution(std::vector<std::pair<int, std::vector<int>>> fleet_allocation)
{
  auto items = instance_.items();
  std::vector<int> values(num_vars(), 0);
  // for (int i = 0; i < num_vars(); ++i)
  //   values[i] = 0;

  int y1_size = vars_.y1_.getSize();
  int y2_size = vars_.y2_.getSize();
  int U_size = vars_.U_.getSize();

  const auto &precedence = instance_.precedence_matrix();
  const auto &successors_for_transport_per_item = instance_.successors_for_transport_per_item();

  std::unordered_set<int> items_already_unloaded; // keeps items already moved to the buffer area or loaded to a vehicle.

  for (int slot = 0; slot < fleet_allocation.size(); ++slot)
  {
    const auto vehicle = fleet_allocation[slot].first;
    const auto &items_of_vehicle = fleet_allocation[slot].second;
    // std::cout << vehicle << ":";
    for (auto item : items_of_vehicle)
    {
      // std::cout << " " << item;
      // std::cout << "y1 " << slot << " " << item << std::endl;
      values[vars_.y1_index(slot, item)] = 1;
      auto group = items[item]->group();
      if (reformulated_)
      {
        // std::cout << "Y " << slot << " " << group << std::endl;
        values[y1_size + y2_size + U_size + vars_.Y_index(slot, group)] = 1;
      }

      // look for blocking items that were not already removed.
      for (int j = 0; j < items.size(); ++j)
      {
        if (precedence[j][item]) // j precedes/blocks item and item is available for transport.
        {
          if (items_already_unloaded.find(j) == items_already_unloaded.end()) // if not yet moved to buffer or unloaded.
          {
            // std::cout << "U " << j << std::endl;
            values[y1_size + y2_size + vars_.U_index(j)] = 1; // mark precedence violation.
            items_already_unloaded.insert(j);                 // move item to buffer and mark violation at U_j.
          }
        }
      }
      items_already_unloaded.insert(item);
    }

    // std::cout << std::endl;

    // std::cout << "y2 " << vehicle << " " << slot << std::endl;
    values[y1_size + vars_.y2_index(vehicle, slot)] = 1;
  }

  // // print
  // all_vars_array_ = IloNumVarArray(*env_);
  // vars_.fill_all_vars_array(all_vars_array_);
  // for (int i = 0; i < num_vars(); ++i)
  // {
  //   if (values[i] > 0)
  //     std::cout << all_vars_array_[i].getName() << ": " << values[i] << std::endl;
  // }
  // all_vars_array_.end();

  return values;
}