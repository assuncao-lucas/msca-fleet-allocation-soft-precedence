#include <fstream>
#include <sstream>
#include <iomanip>
#include "heuristic_solution.h"

void KSHeuristicSolution::Reset()
{
	is_infeasible_ = false;
	is_feasible_ = false;
	is_optimal_ = false;
	found_x_integer_ = false;
	time_spent_building_kernel_buckets_ = 0.0;
	total_time_spent_ = 0.0;
	num_items_loaded_ = 0;
	num_unproductive_moves_ = 0;
	lb_ = -1.0;
}

std::string KSHeuristicSolution::GenerateFileName(bool solve_cutting_plane, bool add_symmetry_breaking, bool vehicle_sequencing_model, bool item_sequencing_model, bool vehicle_slots_model, bool reformulate, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool cluster_buckets_by_item_group)
{
	std::string algo = "ks_";
	if (solve_cutting_plane)
	{
		algo += "cb_";
		// if (cccs)
		// 	algo += "CCCs_";
	}

	if (vehicle_sequencing_model)
		algo += "vehc_seq_model";
	if (item_sequencing_model)
		algo += "itm_seq_model";
	if (vehicle_slots_model)
		algo += "vehc_slot_model";

	if (add_symmetry_breaking)
		algo += "_sym_break";

	if (reformulate)
		algo += "_reform";

	std::stringstream ss_decay_factor;
	ss_decay_factor << std::fixed << std::setprecision(2) << ks_decay_factor;
	if (!cluster_buckets_by_item_group)
		algo += /*"_b" + std::to_string(ks_max_size_bucket) +*/ "_[" + std::to_string(ks_min_time_limit) + "," + std::to_string(ks_max_time_limit) + "]_d" + ss_decay_factor.str();
	else
		algo += "_cluster_by_group_[" + std::to_string(ks_min_time_limit) + "," + std::to_string(ks_max_time_limit) + "]_d" + ss_decay_factor.str();

	return algo;
}

void KSHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const
{
	std::fstream file;
	std::string path = "..//solutions//";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out);

	file << std::setprecision(5) << std::fixed;

	if (is_infeasible_)
		file << "STATUS: INFEASIBLE" << std::endl;
	else if (found_x_integer_)
		file << "STATUS: FOUND INTEGER FEASIBLE SOLUTION" << std::endl;
	else
		file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;

	file << "time building kernel and buckets (s): " << time_spent_building_kernel_buckets_ << std::endl;
	file << "total time (s): " << total_time_spent_ << std::endl;
	file << "LB: " << lb_ << std::endl;
	file << "num_items_loaded: " << num_items_loaded_ << std::endl
		 << "num_unproductive_moves: " << num_unproductive_moves_ << std::endl;

	file.close();
}