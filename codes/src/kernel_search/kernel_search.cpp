#include <cmath>
#include <algorithm>
#include <optional>
#include <queue>
#include <vector>
#include "src/kernel_search/kernel_search.h"
#include "src/graph_algorithms.h"
#include "src/exact/formulations.h"
#include "src/general.h"

std::pair<std::vector<std::pair<int, std::vector<int>>>, std::unordered_set<int>> KernelSearch::findInitialSolution(const Instance &instance)
{
    // Custom comparator for max-heap (non-ascending by value)
    struct MaxCompare
    {
        bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const
        {
            return a.second < b.second; // bigger value = higher priority
        }
    };

    struct MinCompare
    {
        bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const
        {
            return a.second > b.second; // bigger value = higher priority
        }
    };

    class MaxTracker
    {
    private:
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, MaxCompare> pq;

    public:
        void insert(int index, int value)
        {
            pq.push({index, value});
        }

        std::pair<int, int> top() const
        {
            return pq.top();
        }

        // Update the current max (remove + reinsert with new value)
        void updateTop(int newValue)
        {
            if (pq.empty())
                return;
            std::pair<int, int> p = pq.top();
            pq.pop();
            pq.push({p.first, newValue});
        }

        void removeTop()
        {
            if (!pq.empty())
                pq.pop();
        }

        bool empty() const
        {
            return pq.empty();
        }
    };

    class MinTracker
    {
    private:
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, MinCompare> pq;

    public:
        void insert(int index, int value)
        {
            pq.push({index, value});
        }

        std::pair<int, int> top() const
        {
            return pq.top();
        }

        // Update the current max (remove + reinsert with new value)
        void updateTop(int newValue)
        {
            if (pq.empty())
                return;
            std::pair<int, int> p = pq.top();
            pq.pop();
            pq.push({p.first, newValue});
        }

        void removeTop()
        {
            if (!pq.empty())
                pq.pop();
        }

        bool empty() const
        {
            return pq.empty();
        }
    };

    MaxTracker max_cardinality_group_heap;
    MaxTracker max_capacity_vehicle_heap;

    const auto items_for_transport_per_group = instance.items_for_transport_per_group();
    const auto precedence_matrix = instance.precedence_matrix();
    const int num_items = instance.num_items();
    std::unordered_map<int, std::list<int>> items_per_group;
    const auto fleet = instance.fleet();
    const int num_vehicles = instance.num_vehicles();

    std::unordered_map<int, MinTracker> ordered_items_per_group;         // order in terms of number of blocking items.
    std::unordered_map<int, std::pair<int, int>> allocation_per_vehicle; // vehicle -> {group, quantity of items}

    for (auto [group, items] : items_for_transport_per_group)
    {
        max_cardinality_group_heap.insert(group, items.size());
        items_per_group[group] = std::list<int>(items.begin(), items.end());

        MinTracker ordered_items_of_group;
        for (auto item : items)
        {
            int num_precedents = 0;
            for (int i = 0; i < num_items; ++i)
                if (precedence_matrix[i][item])
                    ++num_precedents;

            ordered_items_of_group.insert(item, num_precedents);
        }
        ordered_items_per_group[group] = ordered_items_of_group;
    }

    for (int i = 0; i < num_vehicles; ++i)
    {
        auto vehicle = fleet[i];
        max_capacity_vehicle_heap.insert(i, vehicle->capacity());
    }

    while (!max_capacity_vehicle_heap.empty() && !max_cardinality_group_heap.empty())
    {
        auto max_cap_vehicle_index = max_capacity_vehicle_heap.top().first;
        auto max_cap_vehicle = max_capacity_vehicle_heap.top().second;
        max_capacity_vehicle_heap.removeTop();
        auto max_cardinality_group_index = max_cardinality_group_heap.top().first;
        auto max_cardinality_group = max_cardinality_group_heap.top().second;
        auto num_items_allocated = std::min(max_cap_vehicle, max_cardinality_group);

        allocation_per_vehicle[max_cap_vehicle_index] = {max_cardinality_group_index,
                                                         num_items_allocated};

        (num_items_allocated < max_cardinality_group) ? max_cardinality_group_heap.updateTop(max_cardinality_group - num_items_allocated) : max_cardinality_group_heap.removeTop();
    }

    bool optimize_within_vehicle = true;
    bool optimize_vehicle_ordering = true;

    std::unordered_set<int> all_items_selected;
    std::vector<int> vehicle_where_item_was_allocated(num_items, -1);
    std::vector<std::pair<int, std::vector<int>>> vehicles_allocation;
    std::unordered_map<int, std::vector<int>> vehicles_allocation_naive;

    for (const auto &[vehicle, group_allocation] : allocation_per_vehicle)
    {
        std::vector<int> items_for_vehicle;
        // std::cout << "vehicle " << vehicle << ": " << group_allocation.second << " items of group " << group_allocation.first << std::endl;
        auto group = group_allocation.first;
        auto quantity_to_allocate = group_allocation.second;
        auto &items_group = items_per_group[group];
        auto &ordered_items_group = ordered_items_per_group[group];
        for (int item_iter = 0; item_iter < quantity_to_allocate; ++item_iter)
        {
            int curr_item = -1;
            if (optimize_within_vehicle)
            {
                curr_item = ordered_items_group.top().first;
                ordered_items_group.removeTop();
            }
            else
            {
                curr_item = items_group.back();
                items_group.pop_back();
            }
            all_items_selected.insert(curr_item);
            vehicle_where_item_was_allocated[curr_item] = vehicle;
            items_for_vehicle.push_back(curr_item);
        }
        if (optimize_vehicle_ordering)
            vehicles_allocation_naive[vehicle] = items_for_vehicle;
        else
            vehicles_allocation.push_back({vehicle, items_for_vehicle});

        // std::cout << "vehicle " << vehicle << std::endl;
        // for (auto item : items_for_vehicle)
        //     std::cout << " " << item;
        // std::cout << std::endl;
    }

    MinTracker num_allocated_precedents_per_vehicle;
    for (const auto &[vehicle, vehicle_items] : vehicles_allocation_naive)
    {
        std::unordered_set<int> removed_precedents_from_items_in_vehicle; // keep track of precedents already counted, as to avoid double counting for two or more items of a same vehicle.
        int num_allocated_precedents = 0;
        for (const auto item : vehicle_items)
        {
            for (int i = 0; i < num_items; ++i)
            {
                // only add if precedent item belongs to the solution and was allocated to a different vehicle.
                if ((vehicle_where_item_was_allocated[i] != -1) && (vehicle_where_item_was_allocated[i] != vehicle) && (precedence_matrix[i][item]) && removed_precedents_from_items_in_vehicle.find(i) == removed_precedents_from_items_in_vehicle.end())
                {
                    removed_precedents_from_items_in_vehicle.insert(i);
                    ++num_allocated_precedents;
                }
            }
        }
        num_allocated_precedents_per_vehicle.insert(vehicle, num_allocated_precedents);
        // std::cout << "naive" << std::endl
        //           << "vehicle " << vehicle << "| " << num_allocated_precedents << " allocated precedents" << std::endl;
    }

    while (!num_allocated_precedents_per_vehicle.empty())
    {
        int vehicle = num_allocated_precedents_per_vehicle.top().first;
        num_allocated_precedents_per_vehicle.removeTop();
        vehicles_allocation.push_back({vehicle, vehicles_allocation_naive[vehicle]});

        // std::cout << "optimized" << std::endl
        //           << "vehicle " << vehicle << std::endl;
    }
    // getchar();
    // getchar();

    return {vehicles_allocation, all_items_selected};
}

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

void KernelSearch::BuildKernelAndBuckets(Model &MIPformulation, std::list<UserCut *> *initial_cuts, KSHeuristicSolution *solution, std::optional<std::unordered_set<int>> initial_kernel_items, int ks_max_size_bucket, bool multithreading)
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

    // std::cout << "Before sorting" << std::endl;
    // for (auto &item : item_value_red_cost)
    // {
    //     std::cout << item.item << " " << item.value << " " << item.reduced_cost << std::endl;
    // }

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

    // std::cout << "After sorting" << std::endl;
    // for (auto &item : item_value_red_cost)
    // {
    //     std::cout << item.item << " " << item.value << " " << item.reduced_cost << std::endl;
    // }

    boost::dynamic_bitset<> already_added_items(num_items, 0);
    int size_kernel = 0, iter_item_counter = 0;
    // add to kernal all initial items given.
    if (initial_kernel_items.has_value() && !(initial_kernel_items->empty()))
    {
        // std::cout << "size kernal: " << initial_kernel_items->size() << std::endl;
        // getchar();
        // getchar();
        size_kernel = initial_kernel_items->size();
        for (auto item : *initial_kernel_items)
        {
            curr_kernel_bitset_[item] = 1;
            already_added_items[item] = 1;
            // std::cout << item << std::endl;
        }
    }
    else
    {
        size_kernel = std::min(num_items_for_transport, ks_max_size_bucket);

        for (int i = 0; i < size_kernel; ++i)
        {
            auto item = item_value_red_cost[i].item;
            curr_kernel_bitset_[item] = 1;
            already_added_items[item] = 1;
        }
        iter_item_counter = size_kernel;
    }

    // std::cout << "num_items_for_transport - size_kernel = " << num_items_for_transport << " - " << size_kernel << std::endl;
    int num_buckets = std::ceil(1.0 * ((num_items_for_transport - size_kernel)) / ks_max_size_bucket);
    // std::cout << num_buckets << std::endl;
    buckets_bitsets_ = std::vector<boost::dynamic_bitset<>>(num_buckets, boost::dynamic_bitset<>(num_items, 0));

    // int items_added = size_kernel; // since already added some vertices to the kernel.
    for (int curr_bucket = 0; curr_bucket < num_buckets; ++curr_bucket)
    {
        // std::cout << "bucket " << curr_bucket << std::endl;
        int num_elements_in_bucket = std::min(ks_max_size_bucket, num_items_for_transport - size_kernel - curr_bucket * ks_max_size_bucket);
        // std::cout << "num elements in bucket " << num_elements_in_bucket << std::endl;
        for (int curr_element_in_bucket = 0; curr_element_in_bucket < num_elements_in_bucket;)
        {
            // std::cout << num_items_for_transport << " > " << iter_item_counter << std::endl;
            auto item = item_value_red_cost[iter_item_counter].item;
            // std::cout << "opa" << std::endl;
            if (already_added_items[item] == 0)
            {
                buckets_bitsets_[curr_bucket][item] = 1;
                already_added_items[item] = 1;
                ++curr_element_in_bucket;
            }
            ++iter_item_counter;
        }
    }

    assert(iter_item_counter == num_items_for_transport);

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

    KSHeuristicSolution *solution = new KSHeuristicSolution();

    // build Kernel by solving LP of given problem.
    auto [initial_sol, all_items_solution] = findInitialSolution(instance_);
    auto initial_sol_values = formulation.fillVarValuesFromSolution(initial_sol);
    found_feasible_solution = true;
    // getchar();
    // getchar();
    BuildKernelAndBuckets(formulation, initial_cuts, solution, all_items_solution, ks_max_size_bucket, multithreading);

    solution->time_spent_building_kernel_buckets_ = timer->CurrentElapsedTime(ti);
    if (!solution->is_infeasible_)
    {
        PrintKernelAndBuckets();
        // getchar();
        // getchar();

        // Create the MILP model.
        auto MIPModel = formulation.getClone(false);

        curr_mip_start_vals_ = IloNumArray(MIPModel->env(), formulation.num_vars());
        for (int i = 0; i < curr_mip_start_vals_.getSize(); ++i)
            curr_mip_start_vals_[i] = initial_sol_values[i];

        MIPModel->set_multithreading(multithreading);

        if (feasibility_emphasis)
            MIPModel->set_emphasis_feasibility();

        // initially, set all items (available for transport) to inactive.
        auto null_bitset = boost::dynamic_bitset<>(num_items, 0);
        auto all_items_for_transport_bitset = boost::dynamic_bitset<>(num_items, 0);
        for (auto item : instance_.items_for_transport())
            all_items_for_transport_bitset[item] = 1;

        MIPModel->updateModelVarBounds(null_bitset, all_items_for_transport_bitset, null_bitset);

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

            MIPModel->updateModelVarBounds(curr_items_entering_kernel, curr_items_leaving_reference_kernel, curr_reference_kernel);

            // cplex_.exportModel("model_updated.lp");
            // getchar();
            // getchar();

            // if already found a feasible solution, use it as warm start of the next iteration.
            if (found_feasible_solution)
            {
                MIPModel->addMIPStart(curr_mip_start_vals_);

                std::cout << "added warm start" << std::endl;
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

        if (found_feasible_solution)
        {
            solution->is_feasible_ = solution->found_x_integer_ = true;
            Solution<double> best_solution;

            // compute number of items loaded and unproductive moves of the solution found.
            MIPModel->fillSolution(best_solution, curr_mip_start_vals_);
            solution->num_items_loaded_ = best_solution.num_items_loaded_;
            solution->num_unproductive_moves_ = best_solution.num_unproductive_moves_;
        }

        curr_mip_start_vals_.end();
        delete MIPModel;
        MIPModel = nullptr;
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