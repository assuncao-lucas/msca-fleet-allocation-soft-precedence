#include <cmath>
#include <algorithm>
#include "src/kernel_search/kernel_search.h"
#include "src/graph_algorithms.h"
#include "src/exact/formulations.h"
#include "src/general.h"

KernelSearch::KernelSearch(Instance &instance) : instance_(instance)
{
    const Graph *graph = instance_.graph();
    int num_vertices = graph->num_vertices();
    int num_arcs = graph->num_arcs();

    R0_ = Dijkstra(instance_.graph(), false, false);
    Rn_ = Dijkstra(instance_.graph(), true, false);
    curr_int_y_ = boost::dynamic_bitset<>(num_vertices, 0);
    curr_int_x_ = boost::dynamic_bitset<>(num_arcs, 0);
    curr_kernel_bitset_ = boost::dynamic_bitset<>(num_vertices, 0);
}

KernelSearch::~KernelSearch()
{
    ResetCplex();
    if (R0_)
    {
        delete[] R0_;
        R0_ = nullptr;
    }
    if (Rn_)
    {
        delete[] Rn_;
        Rn_ = nullptr;
    }
}

void KernelSearch::ResetCplex()
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

void KernelSearch::InitCplex(bool multithreading)
{
    env_ = new IloEnv();
    model_ = new IloModel(*env_);
    cplex_ = new IloCplex(*env_);
    curr_x_values_ = IloNumArray(*env_);
    curr_y_values_ = IloNumArray(*env_);
    cplex_->setOut(env_->getNullStream());
    cplex_->setWarning(env_->getNullStream());
    cplex_->setError(env_->getNullStream());
    cplex_->extract(*model_);
    if (!multithreading)
        cplex_->setParam(IloCplex::Param::Threads, 1);
}

void KernelSearch::BuildModel(Formulation formulation, bool linearly_relaxed, bool disable_all_binary_vars, bool export_model)
{
    const Graph *graph = instance_.graph();
    int num_arcs = graph->num_arcs();
    int num_routes = instance_.num_vehicles();
    int budget = instance_.uncertainty_budget();
    const auto num_vertices = graph->num_vertices();

    IloNumVarArray flow_vars;

    switch (formulation)
    {
    case Formulation::baseline:
    {
        AllocateMasterVariablesBaseline(*env_, master_vars_, instance_, false, linearly_relaxed, disable_all_binary_vars);
        // budget + 1 to consider level 0 of budget! 0,..., budget
        flow_vars = IloNumVarArray(*env_, num_vertices * (budget + 1), 0, IloInfinity, ILOFLOAT);
        break;
    }
    case Formulation::single_commodity:
    {
        AllocateMasterVariablesSingleCommodity(*env_, master_vars_, instance_, false, linearly_relaxed, disable_all_binary_vars);
        // budget + 1 to consider level 0 of budget! 0,..., budget
        flow_vars = IloNumVarArray(*env_, num_arcs * (budget + 1), 0, IloInfinity, ILOFLOAT);
        break;
    }
    default:
        throw "invalid formulation";
    }

    slack_ = IloNumVar(*env_, 0, num_routes, ILOFLOAT);

    curr_mip_start_vars_ = IloNumVarArray(*env_, num_vertices + num_arcs);
    curr_mip_start_vals_ = IloNumArray(*env_, num_vertices + num_arcs);

    for (IloInt i = 0; i < num_vertices; ++i)
        curr_mip_start_vars_[i] = master_vars_.y[i];

    for (IloInt i = 0; i < num_arcs; ++i)
        curr_mip_start_vars_[i + num_vertices] = master_vars_.x[i];

    switch (formulation)
    {
    case Formulation::baseline:
    {
        PopulateByRowCompactBaseline(*cplex_, *env_, *model_, master_vars_, flow_vars, instance_, R0_, Rn_, false, export_model);
        break;
    }
    case Formulation::single_commodity:
    {
        PopulateByRowCompactSingleCommodity(*cplex_, *env_, *model_, master_vars_, flow_vars, instance_, R0_, Rn_, false, export_model);
        break;
    }
    default:
        throw "invalid formulation";
    }
}

void KernelSearch::RetrieveSolutionArcVertexValues()
{
    cplex_->getValues(curr_x_values_, master_vars_.x);
    cplex_->getValues(curr_y_values_, master_vars_.y);
    const Graph *graph = instance_.graph();
    int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
    int v1 = 0;

    // std::cout << curr_y_values_ << std::endl;
    int curr_pos = 0;

    curr_int_x_.reset();
    curr_int_y_.reset();

    // we don't have to check every arc: due to the flow conservation constraints, only perform depth first search.
    // since the solution is necessarily connected and all routes depart from origin, just add origin as starting point.
    std::list<int> q;
    q.push_front(0);
    do
    {
        v1 = q.front();
        q.pop_front();

        // if vertex not yet visited and has value equal to one (Note: since y = sum x, y == 0 implies no arc entering/leaving vertex).
        if ((curr_int_y_[v1] == 0) && (double_equals(curr_y_values_[v1], 1.0)))
        {
            for (int v2 : graph->AdjVerticesOut(v1))
            {
                curr_pos = graph->pos(v1, v2);

                if (double_equals(curr_x_values_[curr_pos], 1.0))
                {
                    assert(v2 == 0 | curr_int_y_[v2] == 0);
                    // just to avoid adding again the origin vertex, which is also the depot.
                    if (curr_int_y_[v2] == 0)
                        q.push_front(v2);

                    curr_int_x_[curr_pos] = 1;
                }
            }
            curr_int_y_[v1] = 1;
        }
    } while (!q.empty());

    // fill current MIP start.
    for (IloInt i = 0; i < num_vertices; ++i)
        curr_mip_start_vals_[i] = curr_y_values_[i];

    for (IloInt i = 0; i < num_arcs; ++i)
        curr_mip_start_vals_[i + num_vertices] = curr_x_values_[i];
}

void KernelSearch::BuildHeuristicSolution(KSHeuristicSolution *solution)
{
    const Graph *graph = instance_.graph();
    int v1 = 0, curr_route_index = 0;
    int last_vertex_added = 0;
    bool has_added_vertex_to_route = false;
    VertexStatus *status = nullptr;
    Route *curr_route = nullptr;
    GArc *curr_arc = nullptr;
    size_t cont = 0;

    std::list<int> q;
    q.push_back(0);

    do
    {
        v1 = q.front();
        q.pop_front();

        // end current route
        if ((v1 == 0) && has_added_vertex_to_route)
        {
            curr_route = &((solution->routes_vec_)[curr_route_index]);
            curr_arc = (*graph)[last_vertex_added][v1];
            if (curr_arc == nullptr)
                throw "Inappropriate addition of vertex to route";
            //   else
            //   {
            //     (curr_route->time_) += curr_arc->dist();
            //   }

            // compute route's profit sum and max duration.
            auto [route_sum_profits, route_max_duration] = instance_.ComputeRouteCosts(curr_route->vertices_);

            curr_route->time_ = route_max_duration;
            solution->profits_sum_ += route_sum_profits;
            curr_route->sum_profits_ = route_sum_profits;

            // auto [route_sum_profits2, route_max_duration2] = instance_.ComputeRouteCostsRec(*curr_route, true);

            // if (!double_equals(route_sum_profits, route_sum_profits2))
            //     std::cout << route_sum_profits << " x " << route_sum_profits2 << std::endl;

            // if (!double_equals(route_max_duration, route_max_duration2))
            //     std::cout << route_max_duration << " x " << route_max_duration2 << std::endl;

            // std::cout << *curr_route << std::endl;

            ++curr_route_index;
            has_added_vertex_to_route = false;
            last_vertex_added = 0;
            continue; // in this case, should avoid the last step of the loop that adds to stack the neighbors of 0.
        }
        else if (v1 != 0)
        {
            // add vertex to current route
            curr_route = &((solution->routes_vec_)[curr_route_index]);
            status = &((solution->vertex_status_vec_)[v1]);
            has_added_vertex_to_route = true;

            // remove from list of unvisited_vertices
            (solution->unvisited_vertices_).erase(status->pos_);

            // adds vertex to route
            status->selected_ = true;
            status->route_ = curr_route_index;
            status->pos_ = (curr_route->vertices_).insert((curr_route->vertices_).end(), v1);

            curr_arc = (*graph)[last_vertex_added][v1];
            assert(curr_arc != nullptr);

            last_vertex_added = v1;
        }

        for (int v2 : graph->AdjVerticesOut(v1))
        {
            if ((curr_int_x_)[graph->pos(v1, v2)])
            {
                ++cont;
                q.push_front(v2);
            }
        }
    } while (!(q.empty()));

    // this is not true because of the decreasing profits...the order in which vertices appear change the profit collected.
    // if (((solution_).unvisited_vertices_).empty())
    // 	(solution_).is_optimal_ = true;

    solution->BuildBitset(instance_);
    // std::cout << solution->bitset_arcs_ << std::endl;
    // std::cout << solution->bitset_vertices_ << std::endl;
    assert((cont == (curr_int_x_).count()));
    assert(double_equals(curr_best_solution_value_, solution->profits_sum_));
    // if (!double_equals(curr_best_solution_value_, solution->profits_sum_))
    //     std::cout << curr_best_solution_value_ << " != " << solution->profits_sum_ << std::endl;
}

void KernelSearch::BuildKernelAndBuckets(Formulation formulation, KSHeuristicSolution *solution, int ks_max_size_bucket, bool multithreading)
{
    InitCplex(multithreading);
    BuildModel(formulation, true, false, false);
    curr_kernel_bitset_.reset();

    cplex_->extract(*model_);
    // cplex_.exportModel("kernel_search_lp.lp");
    // getchar(); getchar();
    if (!(cplex_->solve()))
    {
        if ((cplex_->getCplexStatus() == IloCplex::Infeasible) || (cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
            solution->is_infeasible_ = true;

        return;
    }

    // Sort vertices in non-ascending order of y values. For vertices with y == 0, sort by reduced costs.
    IloNumArray y_values(*env_);
    IloNumArray y_reduced_costs(*env_);

    cplex_->getValues(y_values, master_vars_.y);
    cplex_->getReducedCosts(y_reduced_costs, master_vars_.y);

    struct VertexValueReducedCost
    {
        int vertex;
        double value;
        double reduced_cost;
    };

    const int num_vertices = instance_.graph()->num_vertices();
    const int num_mandatory = instance_.num_mandatory();
    std::vector<VertexValueReducedCost> vertex_value_red_cost;
    for (int i = 0; i < num_vertices; ++i)
        vertex_value_red_cost.push_back(VertexValueReducedCost{.vertex = i, .value = y_values[i], .reduced_cost = y_reduced_costs[i]});

    // std::cout << "Before sorting" << std::endl;
    // for (auto &item : vertex_value_red_cost)
    // {
    //     std::cout << item.vertex << " " << item.value << " " << item.reduced_cost << std::endl;
    // }

    auto compare = [/*num_mandatory*/](const VertexValueReducedCost &a, const VertexValueReducedCost &b)
    {
        // if ((a.vertex <= num_mandatory) && (b.vertex > num_mandatory))
        //     return true;
        // if ((a.vertex > num_mandatory) && (b.vertex <= num_mandatory))
        //     return false;

        return double_equals(a.value, b.value) ? a.reduced_cost > b.reduced_cost : a.value > b.value;
    };

    // make sure that the origin and the mandatory are left at the beginning of the vector, to force their entry in the Kernel.
    auto start = vertex_value_red_cost.begin();
    start += (num_mandatory + 1);
    std::sort(start, vertex_value_red_cost.end(), compare);

    // std::cout << "After sorting" << std::endl;
    // for (auto &item : vertex_value_red_cost)
    // {
    //     std::cout << item.vertex << " " << item.value << " " << item.reduced_cost << std::endl;
    // }

    // start by building kernel with at least all the mandatory vertices (+ origin).
    int size_kernel = std::min(num_vertices, std::max(num_mandatory + 1, ks_max_size_bucket));

    for (int i = 0; i < size_kernel; ++i)
        curr_kernel_bitset_[vertex_value_red_cost[i].vertex] = 1;

    int num_buckets = std::ceil(1.0 * ((num_vertices - size_kernel)) / ks_max_size_bucket);
    buckets_bitsets_ = std::vector<boost::dynamic_bitset<>>(num_buckets, boost::dynamic_bitset<>(num_vertices, 0));

    int vertices_added = size_kernel; // since already added some vertices to the kernel.
    for (int curr_bucket = 0; curr_bucket < num_buckets; ++curr_bucket)
    {
        // std::cout << "bucket " << curr_bucket << std::endl;
        int num_elements_in_bucket = std::min(ks_max_size_bucket, num_vertices - size_kernel - curr_bucket * ks_max_size_bucket);
        // std::cout << "num elements in bucket " << num_elements_in_bucket << std::endl;
        for (int curr_element_in_bucket = 0; curr_element_in_bucket < num_elements_in_bucket; ++curr_element_in_bucket)
        {
            buckets_bitsets_[curr_bucket][vertex_value_red_cost[vertices_added].vertex] = 1;
            ++vertices_added;
        }
    }

    assert(vertices_added == num_vertices);
}

KSHeuristicSolution *KernelSearch::Run(Formulation formulation, int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool feasibility_emphasis, bool multithreading)
{
    Timestamp *ti = NewTimestamp();
    Timer *timer = GetTimer();
    timer->Clock(ti);

    std::cout << std::setprecision(2) << std::fixed;
    const int num_mandatory = instance_.num_mandatory();
    const Graph *graph = instance_.graph();
    const int num_vertices = graph->num_vertices();
    const int num_arcs = graph->num_arcs();
    const int num_routes = instance_.num_vehicles();

    KSHeuristicSolution *solution = new KSHeuristicSolution(num_vertices, num_arcs, num_routes);

    // build Kernel by solving LP of given problem.
    BuildKernelAndBuckets(formulation, solution, ks_max_size_bucket, multithreading);
    solution->time_spent_building_kernel_buckets_ = timer->CurrentElapsedTime(ti);

    if (!solution->is_infeasible_)
    {
        // PrintKernelAndBuckets();

        // delete LP model and create the MILP model.
        ResetCplex();
        InitCplex(multithreading);

        BuildModel(formulation, false, true, false); // initially, set all binary variables to zero.

        cplex_->setParam(IloCplex::Param::ClockType, 2);
        if (feasibility_emphasis)
        {
            // std::cout << "feasibility_emphasis" << std::endl;
            cplex_->setParam(IloCplex::Param::Emphasis::MIP, IloCplex::MIPEmphasisFeasibility);
        }

        double curr_time_limit_iteration = ks_max_time_limit;

        int curr_bucket_index = -1; // starts from kernel.
        int total_num_buckets = buckets_bitsets_.size();

        auto curr_reference_kernel = curr_kernel_bitset_;
        auto curr_vertices_entering_kernel = curr_kernel_bitset_;
        auto curr_vertices_leaving_reference_kernel = boost::dynamic_bitset<>(num_vertices, 0);

        for (int curr_bucket_index = -1; curr_bucket_index < total_num_buckets; ++curr_bucket_index)
        {
            cplex_->setParam(IloCplex::Param::TimeLimit, curr_time_limit_iteration);
            // std::cout << "curr_time_limit_iteration: " << curr_time_limit_iteration << std::endl;
            // update the reference kernel to the current kernel (+ current bucket, if not the first iteration).
            if (curr_bucket_index >= 0)
            {
                curr_reference_kernel = curr_kernel_bitset_ | buckets_bitsets_[curr_bucket_index];
                curr_vertices_entering_kernel |= buckets_bitsets_[curr_bucket_index];
            }

            // std::cout << " bucket index " << curr_bucket_index << std::endl;
            // std::cout << " current bucket ";
            // curr_bucket_index >= 0 ? std::cout << buckets_bitsets_[curr_bucket_index] << std::endl : std::cout << " - " << std::endl;
            // std::cout << " kernel: " << curr_kernel_bitset_ << std::endl;
            // std::cout << " reference kernel: " << curr_reference_kernel << std::endl;
            // std::cout << " IN kernel: " << curr_vertices_entering_kernel << std::endl;
            // std::cout << " OUT ref kernel: " << curr_vertices_leaving_reference_kernel << std::endl;
            // std::cout << " curr sol: " << curr_int_y_ << std::endl;
            // std::cout << " best cost: " << curr_best_solution_value_ << std::endl;
            // cplex_.exportModel("model_before.lp");
            //  Enable in the model the variables that are active in the Kernel.

            UpdateModelVarBounds(curr_vertices_entering_kernel, curr_vertices_leaving_reference_kernel, curr_reference_kernel);

            // cplex_.exportModel("model_updated.lp");
            // getchar();
            // getchar();

            // if already found a feasible solution, use it as warm start of the next iteration.
            if (found_int_x_)
            {
                cplex_->addMIPStart(curr_mip_start_vars_, curr_mip_start_vals_, IloCplex::MIPStartSolveMIP);

                // std::cout << "added warm start" << std::endl;
            }

            bool found_better_solution = false;

            if (cplex_->solve())
            {
                found_int_x_ = true;
                double solution_value = cplex_->getObjValue();
                // std::cout << "found feasible solution with cost " << solution_value << std::endl;

                // only update Kernel if found a solution with strictly better objective function value!
                if (double_greater(solution_value, curr_best_solution_value_))
                {
                    found_better_solution = true;
                    curr_best_solution_value_ = solution_value;
                    RetrieveSolutionArcVertexValues();

                    // update kernel with the possibly new vertices used in the current solution found.
                    curr_vertices_entering_kernel = curr_int_y_ - curr_kernel_bitset_;
                    curr_kernel_bitset_ |= curr_vertices_entering_kernel;

                    curr_vertices_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;

                    // KSHeuristicSolution *new_solution = new KSHeuristicSolution(num_vertices, num_arcs, num_routes);
                    // new_solution->is_feasible_ = new_solution->found_x_integer_ = true;
                    // BuildHeuristicSolution(new_solution);
                    // delete new_solution;
                    // new_solution = nullptr;
                }
            }
            else
            {
                if ((cplex_->getCplexStatus() == IloCplex::Infeasible) || (cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
                {
                    solution->is_infeasible_ = true;
                    assert(found_int_x_ == false);
                    break;
                }
            }

            // if hasn't found a batter solution, nothing changes in the kernel...we only remove from reference kernel the vertices added in this iteration.
            if (!found_better_solution)
            {
                curr_vertices_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;
                curr_vertices_entering_kernel.reset();
            }
            // std::cout << "status " << cplex_->getStatus() << std::endl;
            // std::cout << "current sol: " << curr_int_y_ << std::endl;

            curr_time_limit_iteration = std::max(curr_time_limit_iteration * ks_decay_factor, 1.0 * ks_min_time_limit);
        }
    }

    if (found_int_x_)
    {
        solution->is_feasible_ = solution->found_x_integer_ = true;
        BuildHeuristicSolution(solution);
    }

    // std::cout << "Best solution found: " << curr_best_solution_value_ << " " << curr_int_y_ << std::endl;
    // std::cout << "Elapsed time: " << timer->CurrentElapsedTime(ti) << std::endl;

    solution->total_time_spent_ = timer->CurrentElapsedTime(ti);

    // std::cout << *solution << std::endl;

    delete (ti);
    ti = nullptr;

    return solution;
}

void KernelSearch::UpdateModelVarBounds(boost::dynamic_bitset<> &vars_entering_kernel, boost::dynamic_bitset<> &vars_leaving_reference_kernel, boost::dynamic_bitset<> &curr_reference_kernel)
{
    const Graph *graph = instance_.graph();
    // add new variables to kernel.
    for (size_t vertex = vars_entering_kernel.find_first(); vertex != boost::dynamic_bitset<>::npos; vertex = vars_entering_kernel.find_next(vertex))
    {
        // activate vertex.
        master_vars_.y[vertex].setUB(1.0);

        // activate all arcs entering vertex that come from another active vertex.
        for (auto &vertex_in : graph->AdjVerticesIn(vertex))
            if (curr_reference_kernel[vertex_in] == 1)
                master_vars_.x[graph->pos(vertex_in, vertex)].setUB(1.0);

        // activate all arcs leaving vertex that go to another active vertex.
        for (auto &vertex_out : graph->AdjVerticesOut(vertex))
            if (curr_reference_kernel[vertex_out] == 1)
                master_vars_.x[graph->pos(vertex, vertex_out)].setUB(1.0);
    }

    // remove variables leaving kernel.
    for (size_t vertex = vars_leaving_reference_kernel.find_first(); vertex != boost::dynamic_bitset<>::npos; vertex = vars_leaving_reference_kernel.find_next(vertex))
    {
        // deactivate vertex.
        master_vars_.y[vertex].setUB(0.0);

        // deactivate all arcs entering vertex.
        for (auto &vertex_in : graph->AdjVerticesIn(vertex))
            master_vars_.x[graph->pos(vertex_in, vertex)].setUB(0.0);

        // deactivate all arcs leaving vertex.
        for (auto &vertex_out : graph->AdjVerticesOut(vertex))
            master_vars_.x[graph->pos(vertex, vertex_out)].setUB(0.0);
    }
}

void KernelSearch::PrintKernelAndBuckets()
{
    std::cout << "Kernel: ";
    size_t vertex = curr_kernel_bitset_.find_first();
    while (vertex != boost::dynamic_bitset<>::npos)
    {
        std::cout << vertex << " ";
        vertex = curr_kernel_bitset_.find_next(vertex);
    }

    std::cout << std::endl;

    for (int j = 0; j < buckets_bitsets_.size(); ++j)
    {
        std::cout << "Bucket " << j << ": ";
        vertex = buckets_bitsets_[j].find_first();
        while (vertex != boost::dynamic_bitset<>::npos)
        {
            std::cout << vertex << " ";
            vertex = buckets_bitsets_[j].find_next(vertex);
        }
        std::cout << std::endl;
    }
}