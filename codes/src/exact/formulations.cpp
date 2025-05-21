#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/general.h"
#include "src/exact/formulations.h"

ILOSTLBEGIN

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
    for (int j = 0; j < num_items_for_transport; ++i)
      vars.x_var_to_index.addEntry(i, j);

  vars.x = IloNumVarArray(env, vars.x_var_to_index.numVars(), 0.0, binary_upper_bound, type);

  // fill z variables and index map.

  for (int i = 0; i < num_items_for_transport; ++i)
    for (int h = 0; h < num_items_for_transport; ++h)
      if (i != h)
        vars.z_var_to_index.addEntry(i, h);

  vars.z = IloNumVarArray(env, vars.z_var_to_index.numVars(), 0.0, binary_upper_bound, type);

  // fill w variables and index map.

  for (int i = 0; i < num_items_for_transport; ++i)
    for (auto &k : successors_for_transport_per_item[i])
      vars.w_var_to_index.addEntry(i, k);

  vars.w = IloNumVarArray(env, vars.w_var_to_index.numVars(), 0.0, binary_upper_bound, type);

  // fill U0 variables and index map.

  for (int i = 0; i < num_items; ++i)
    if (!(successors_for_transport_per_item[i].empty()))
      vars.U_var_to_index.addEntry(0, i);

  vars.U = IloNumVarArray(env, vars.U_var_to_index.numVars(), 0.0, binary_upper_bound, type);
}

void vehicleSequencingModel(Instance &inst, double time_limit, bool solve_relax, bool export_model)
{
  const int num_vertices = inst.num_items();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  VehicleSequencingModelVariables vars{};

  allocateVehicleSequencingModelVariables(env, vars, inst, solve_relax);

  populateByRowVehicleSequencingModel(cplex, env, model, vars, inst, export_model);

  // optimize(cplex, env, model, vars, inst, time_limit);

  cplex.end();
  env.end();
}

void populateByRowVehicleSequencingModel(IloCplex &cplex, IloEnv &env, IloModel &model, VehicleSequencingModelVariables &vars, const Instance &instance, bool export_model)
{
  // add objective function.
  // IloExpr obj(env);
  // const int budget = instance.uncertainty_budget();

  // for (int j = num_mandatory + 1; j < num_vertices; ++j)
  // {
  //   const auto &vertex_info = vertices_info[j];
  //   obj += operator*(vertex_info.profit_, master_vars.y[j]);
  //   obj -= operator*(vertex_info.decay_ratio_ * route_limit, master_vars.y[j]);

  //   for (const int &i : graph->AdjVerticesIn(j))
  //     obj += operator*(vertex_info.decay_ratio_, f[f_var_to_index(graph->pos(i, j), budget, num_arcs)]);
  // }

  // model.add(IloMaximize(env, obj));
  // obj.end();

  if (export_model)
  {
    // add name to variables.
    for (auto &[key, value] : vars.x_var_to_index.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "x(%d)(%d)", value.first, value.second);
      vars.x[key].setName(strnum);
    }

    for (auto &[key, value] : vars.z_var_to_index.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "z(%d)(%d)", value.first, value.second);
      vars.z[key].setName(strnum);
    }

    for (auto &[key, value] : vars.w_var_to_index.index_to_var_map())
    {
      char strnum[26];
      sprintf(strnum, "w(%d)(%d)", value.first, value.second);
      vars.w[key].setName(strnum);
    }

    for (auto &[key, value] : vars.U_var_to_index.index_to_var_map())
    {
      char strnum[15];
      sprintf(strnum, "U(%d)", value.second);
      vars.U[key].setName(strnum);
    }

    cplex.exportModel("vehicle_sequencing_model.lp");
  }
}