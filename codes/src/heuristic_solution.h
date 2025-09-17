#pragma once

#include <list>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "route.h"
#include "graph.h"
#include "instance.h"
#include "general.h"

class KSHeuristicSolution
{
public:
	KSHeuristicSolution() = default;
	~KSHeuristicSolution() = default;
	double time_spent_building_kernel_buckets_ = 0.0;
	double total_time_spent_ = 0.0;
	double lb_ = -1.0;
	bool found_x_integer_ = false;
	bool is_infeasible_ = false;
	bool is_feasible_ = false;
	bool is_optimal_ = false;
	int num_items_loaded_ = 0;
	int num_unproductive_moves_ = 0;

	virtual void Reset();
	static std::string GenerateFileName(bool solve_cutting_plane, bool add_symmetry_breaking, bool vehicle_sequencing_model, bool item_sequencing_model, bool vehicle_slots_model, bool reformulate, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool cluster_buckets_by_item_group);
	void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
};