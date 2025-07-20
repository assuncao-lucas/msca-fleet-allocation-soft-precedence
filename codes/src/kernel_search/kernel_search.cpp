#include <cmath>
#include <algorithm>
#include "src/kernel_search/kernel_search.h"
#include "src/graph_algorithms.h"
#include "src/exact/formulations.h"
#include "src/general.h"

KernelSearch::KernelSearch(Instance &instance) : instance_(instance)
{
    int num_items = instance.num_items();

    curr_solution_items_ = boost::dynamic_bitset<>(num_items, 0);
    // curr_int_x_ = boost::dynamic_bitset<>(num_arcs, 0);
    curr_kernel_bitset_ = boost::dynamic_bitset<>(num_items, 0);
}

KernelSearch::~KernelSearch()
{
}

void KernelSearch::BuildHeuristicSolution(KSHeuristicSolution *solution)
{
    // const Graph *graph = instance_.graph();
    // int v1 = 0, curr_route_index = 0;
    // int last_vertex_added = 0;
    // bool has_added_vertex_to_route = false;
    // VertexStatus *status = nullptr;
    // Route *curr_route = nullptr;
    // GArc *curr_arc = nullptr;
    // size_t cont = 0;

    // std::list<int> q;
    // q.push_back(0);

    // do
    // {
    //     v1 = q.front();
    //     q.pop_front();

    //     // end current route
    //     if ((v1 == 0) && has_added_vertex_to_route)
    //     {
    //         curr_route = &((solution->routes_vec_)[curr_route_index]);
    //         curr_arc = (*graph)[last_vertex_added][v1];
    //         if (curr_arc == nullptr)
    //             throw "Inappropriate addition of vertex to route";
    //         //   else
    //         //   {
    //         //     (curr_route->time_) += curr_arc->dist();
    //         //   }

    //         // compute route's profit sum and max duration.
    //         auto [route_sum_profits, route_max_duration] = instance_.ComputeRouteCosts(curr_route->vertices_);

    //         curr_route->time_ = route_max_duration;
    //         solution->profits_sum_ += route_sum_profits;
    //         curr_route->sum_profits_ = route_sum_profits;

    //         // auto [route_sum_profits2, route_max_duration2] = instance_.ComputeRouteCostsRec(*curr_route, true);

    //         // if (!double_equals(route_sum_profits, route_sum_profits2))
    //         //     std::cout << route_sum_profits << " x " << route_sum_profits2 << std::endl;

    //         // if (!double_equals(route_max_duration, route_max_duration2))
    //         //     std::cout << route_max_duration << " x " << route_max_duration2 << std::endl;

    //         // std::cout << *curr_route << std::endl;

    //         ++curr_route_index;
    //         has_added_vertex_to_route = false;
    //         last_vertex_added = 0;
    //         continue; // in this case, should avoid the last step of the loop that adds to stack the neighbors of 0.
    //     }
    //     else if (v1 != 0)
    //     {
    //         // add vertex to current route
    //         curr_route = &((solution->routes_vec_)[curr_route_index]);
    //         status = &((solution->vertex_status_vec_)[v1]);
    //         has_added_vertex_to_route = true;

    //         // remove from list of unvisited_vertices
    //         (solution->unvisited_vertices_).erase(status->pos_);

    //         // adds vertex to route
    //         status->selected_ = true;
    //         status->route_ = curr_route_index;
    //         status->pos_ = (curr_route->vertices_).insert((curr_route->vertices_).end(), v1);

    //         curr_arc = (*graph)[last_vertex_added][v1];
    //         assert(curr_arc != nullptr);

    //         last_vertex_added = v1;
    //     }

    //     for (int v2 : graph->AdjVerticesOut(v1))
    //     {
    //         if ((curr_int_x_)[graph->pos(v1, v2)])
    //         {
    //             ++cont;
    //             q.push_front(v2);
    //         }
    //     }
    // } while (!(q.empty()));

    // // this is not true because of the decreasing profits...the order in which vertices appear change the profit collected.
    // // if (((solution_).unvisited_vertices_).empty())
    // // 	(solution_).is_optimal_ = true;

    // solution->BuildBitset(instance_);
    // // std::cout << solution->bitset_arcs_ << std::endl;
    // // std::cout << solution->bitset_vertices_ << std::endl;
    // assert((cont == (curr_int_x_).count()));
    // assert(double_equals(curr_best_solution_value_, solution->profits_sum_));
    // // if (!double_equals(curr_best_solution_value_, solution->profits_sum_))
    // //     std::cout << curr_best_solution_value_ << " != " << solution->profits_sum_ << std::endl;
}

void KernelSearch::BuildKernelAndBuckets(Model &MIPformulation, std::list<UserCut *> *initial_cuts, KSHeuristicSolution *solution, int ks_max_size_bucket, bool multithreading)
{
    auto formulation = MIPformulation.getClone(true);
    formulation->set_multithreading(multithreading);

    curr_kernel_bitset_.reset();
    Solution<double> sol;

    if (!(formulation->optimize(-1, false, initial_cuts, nullptr, sol)))
    {
        auto status = formulation->status();
        if ((status == IloCplex::Infeasible) || (status == IloCplex::InfOrUnbd))
            solution->is_infeasible_ = true;

        delete formulation;
        formulation = nullptr;
        return;
    }

    formulation->exportModel("ks_model.lp");

    // Sort vertices in non-ascending order of sum of x values. For vertices with sum x == 0, sort by sum of reduced costs.

    auto items_values = formulation->get_items_values();
    auto items_reduced_costs = formulation->get_items_reduced_costs();

    auto env = formulation->env();
    auto num_vars = formulation->num_vars();

    struct ItemValueReducedCost
    {
        int item;
        double value;
        double reduced_cost;
    };

    const int num_items = instance_.num_items();
    const int num_items_for_transport = instance_.num_items_for_transport();
    std::vector<ItemValueReducedCost> item_value_red_cost;
    for (int i : instance_.items_for_transport())
        item_value_red_cost.push_back(ItemValueReducedCost{.item = i, .value = items_values[i], .reduced_cost = items_reduced_costs[i]});

    std::cout << "Before sorting" << std::endl;
    for (auto &item : item_value_red_cost)
    {
        std::cout << item.item << " " << item.value << " " << item.reduced_cost << std::endl;
    }

    auto compare = [/*num_mandatory*/](const ItemValueReducedCost &a, const ItemValueReducedCost &b)
    {
        // if ((a.vertex <= num_mandatory) && (b.vertex > num_mandatory))
        //     return true;
        // if ((a.vertex > num_mandatory) && (b.vertex <= num_mandatory))
        //     return false;

        return double_equals(a.value, b.value) ? a.reduced_cost > b.reduced_cost : a.value > b.value;
    };

    auto start = item_value_red_cost.begin();
    std::sort(start, item_value_red_cost.end(), compare);

    std::cout << "After sorting" << std::endl;
    for (auto &item : item_value_red_cost)
    {
        std::cout << item.item << " " << item.value << " " << item.reduced_cost << std::endl;
    }

    int size_kernel = std::min(num_items_for_transport, ks_max_size_bucket);

    for (int i = 0; i < size_kernel; ++i)
        curr_kernel_bitset_[item_value_red_cost[i].item] = 1;

    int num_buckets = std::ceil(1.0 * ((num_items_for_transport - size_kernel)) / ks_max_size_bucket);
    buckets_bitsets_ = std::vector<boost::dynamic_bitset<>>(num_buckets, boost::dynamic_bitset<>(num_items, 0));

    int items_added = size_kernel; // since already added some vertices to the kernel.
    for (int curr_bucket = 0; curr_bucket < num_buckets; ++curr_bucket)
    {
        // std::cout << "bucket " << curr_bucket << std::endl;
        int num_elements_in_bucket = std::min(ks_max_size_bucket, num_items_for_transport - size_kernel - curr_bucket * ks_max_size_bucket);
        // std::cout << "num elements in bucket " << num_elements_in_bucket << std::endl;
        for (int curr_element_in_bucket = 0; curr_element_in_bucket < num_elements_in_bucket; ++curr_element_in_bucket)
        {
            buckets_bitsets_[curr_bucket][item_value_red_cost[items_added].item] = 1;
            ++items_added;
        }
    }

    assert(items_added == num_items_for_transport);

    items_values.end();
    items_reduced_costs.end();
    delete formulation;
    formulation = nullptr;
}

KSHeuristicSolution *KernelSearch::Run(Model &formulation, std::list<UserCut *> *initial_cuts, int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool feasibility_emphasis, bool multithreading)
{
    Timestamp *ti = NewTimestamp();
    Timer *timer = GetTimer();
    timer->Clock(ti);

    int num_items = instance_.num_items();

    std::cout << std::setprecision(2) << std::fixed;

    KSHeuristicSolution *solution = new KSHeuristicSolution(0, 0, 0);

    // build Kernel by solving LP of given problem.
    BuildKernelAndBuckets(formulation, initial_cuts, solution, ks_max_size_bucket, multithreading);

    solution->time_spent_building_kernel_buckets_ = timer->CurrentElapsedTime(ti);
    if (!solution->is_infeasible_)
    {
        PrintKernelAndBuckets();

        // Create the MILP model.
        auto MIPModel = formulation.getClone(false);
        curr_mip_start_vals_ = IloNumArray(MIPModel->env());

        MIPModel->set_multithreading(multithreading);

        if (feasibility_emphasis)
            MIPModel->set_emphasis_feasibility();

        // initially, set all items (available for transport) to inactive.
        auto null_bitset = boost::dynamic_bitset<>(num_items, 0);
        auto all_items_for_transport_bitset = boost::dynamic_bitset<>(num_items, 0);
        for (auto item : instance_.items_for_transport())
            all_items_for_transport_bitset[item] = 1;

        MIPModel->UpdateModelVarBounds(null_bitset, all_items_for_transport_bitset, null_bitset);

        double curr_time_limit_iteration = ks_max_time_limit;

        int curr_bucket_index = -1; // starts from kernel.
        int total_num_buckets = buckets_bitsets_.size();

        auto curr_reference_kernel = curr_kernel_bitset_;
        auto curr_items_entering_kernel = curr_kernel_bitset_;
        auto curr_items_leaving_reference_kernel = boost::dynamic_bitset<>(num_items, 0);

        for (int curr_bucket_index = -1; curr_bucket_index < total_num_buckets; ++curr_bucket_index)
        {
            MIPModel->set_time_limit(curr_time_limit_iteration);

            // std::cout << "curr_time_limit_iteration: " << curr_time_limit_iteration << std::endl;
            // update the reference kernel to the current kernel (+ current bucket, if not the first iteration).
            if (curr_bucket_index >= 0)
            {
                curr_reference_kernel = curr_kernel_bitset_ | buckets_bitsets_[curr_bucket_index];
                curr_items_entering_kernel |= buckets_bitsets_[curr_bucket_index];
            }

            std::cout << " bucket index " << curr_bucket_index << std::endl;
            std::cout << " current bucket ";
            curr_bucket_index >= 0 ? std::cout << buckets_bitsets_[curr_bucket_index] << std::endl : std::cout << " - " << std::endl;
            std::cout << " kernel: " << curr_kernel_bitset_ << std::endl;
            std::cout << " reference kernel: " << curr_reference_kernel << std::endl;
            std::cout << " IN kernel: " << curr_items_entering_kernel << std::endl;
            std::cout << " OUT ref kernel: " << curr_items_leaving_reference_kernel << std::endl;
            std::cout << " curr sol: " << curr_solution_items_ << std::endl;
            std::cout << " best cost: " << curr_best_solution_value_ << std::endl;
            // cplex_.exportModel("model_before.lp");
            //  Enable in the model the variables that are active in the Kernel.

            MIPModel->UpdateModelVarBounds(curr_items_entering_kernel, curr_items_leaving_reference_kernel, curr_reference_kernel);

            // cplex_.exportModel("model_updated.lp");
            // getchar();
            // getchar();

            // if already found a feasible solution, use it as warm start of the next iteration.
            if (found_feasible_solution)
            {
                MIPModel->addMIPStart(curr_mip_start_vals_);

                // std::cout << "added warm start" << std::endl;
            }

            bool found_better_solution = false;

            Solution<double> sol;
            if (MIPModel->optimize(-1, false, initial_cuts, nullptr, sol))
            {
                found_feasible_solution = true;
                double solution_value = MIPModel->getObjValue();
                // std::cout << "found feasible solution with cost " << solution_value << std::endl;

                // only update Kernel if found a solution with strictly better objective function value!
                if (double_greater(solution_value, curr_best_solution_value_))
                {
                    found_better_solution = true;
                    curr_best_solution_value_ = solution_value;

                    // (re)fill current MIP start and items that compose the current solution.
                    MIPModel->getSolutionValues(curr_mip_start_vals_);
                    MIPModel->getSolutionItems(curr_solution_items_);

                    // update kernel with the possibly new vertices used in the current solution found.
                    curr_items_entering_kernel = curr_solution_items_ - curr_kernel_bitset_;
                    curr_kernel_bitset_ |= curr_items_entering_kernel;

                    curr_items_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;

                    // KSHeuristicSolution *new_solution = new KSHeuristicSolution(num_vertices, num_arcs, num_routes);
                    // new_solution->is_feasible_ = new_solution->found_x_integer_ = true;
                    // BuildHeuristicSolution(new_solution);
                    // delete new_solution;
                    // new_solution = nullptr;
                }
            }
            // else
            // {
            //     auto status = MIPModel->status();
            //     if ((status == IloCplex::Infeasible) || (status == IloCplex::InfOrUnbd))
            //     {
            //         solution->is_infeasible_ = true;
            //         assert(found_feasible_solution == false);
            //         break;
            //     }
            // }

            // if hasn't found a batter solution, nothing changes in the kernel...we only remove from reference kernel the vertices added in this iteration.
            if (!found_better_solution)
            {
                curr_items_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;
                curr_items_entering_kernel.reset();
            }
            // std::cout << "status " << cplex_->getStatus() << std::endl;
            // std::cout << "current sol: " << curr_int_y_ << std::endl;

            curr_time_limit_iteration = std::max(curr_time_limit_iteration * ks_decay_factor, 1.0 * ks_min_time_limit);
        }

        curr_mip_start_vals_.end();
        delete MIPModel;
        MIPModel = nullptr;
    }

    if (found_feasible_solution)
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