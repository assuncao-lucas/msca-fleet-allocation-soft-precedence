#include <iostream>
#include <string>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <filesystem>
#include "src/instance.h"
#include "src/graph.h"
#include "src/general.h"
#include "src/exact/formulations.h"
#include "src/timer.h"
#include "src/graph_algorithms.h"
#include "src/heuristic_solution.h"
#include "src/kernel_search/kernel_search.h"

namespace fs = std::filesystem;

void GenerateAlgorithmsLatexTablePerInstance(std::string folder)
{
	std::string curr_file;
	std::vector<std::string> algorithms;
	std::vector<std::string> instances;

	instances.push_back("R25_0.8_v2_d0.50_b5.txt");
	instances.push_back("R25_1_v2_d0.50_b5.txt");
	instances.push_back("R50_0.8_v4_d0.50_b5.txt");
	instances.push_back("R50_3_v2_d0.50_b5.txt");
	instances.push_back("RC50_3_v2_d0.50_b5.txt");

	algorithms.push_back("baseline");
	// algorithms.push_back("baseline_benders_lazy_callback");
	// algorithms.push_back("baseline_benders_combined_lazy_callback");
	// algorithms.push_back("baseline_benders_generic_callback");
	// algorithms.push_back("baseline_benders_combined_generic_callback");
	// algorithms.push_back("baseline_benders_combined_cuts_relaxation_generic_callback");

	algorithms.push_back("csc");
	// algorithms.push_back("csc_benders_lazy_callback");
	// algorithms.push_back("csc_benders_combined_lazy_callback");
	// algorithms.push_back("csc_benders_generic_callback");
	// algorithms.push_back("csc_benders_combined_generic_callback");
	// algorithms.push_back("csc_benders_combined_cuts_relaxation_generic_callback");

	std::fstream output;
	std::string output_name;
	output_name = "..//tables//CSV//table_algorithms_per_instance.csv";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	std::vector<double> total_time(algorithms.size(), 0.0);
	std::vector<double> total_avg_gap(algorithms.size(), 0.0);
	std::vector<int> total_num_optimal(algorithms.size(), 0);

	output << std::setprecision(2) << std::fixed;

	output << "instance;";
	for (size_t algo = 0; algo < algorithms.size(); ++algo)
		for (size_t j = 0; j < 2; ++j)
			output << algorithms[algo] << ";";

	output << std::endl;
	output << " ;";

	for (size_t algo = 0; algo < algorithms.size(); ++algo)
		output << "gap;time;";

	output << std::endl;

	for (auto instance : instances)
	{
		// std::cout << instance << std::endl;
		double original_lp = 0.0;
		output << instance << ";";
		for (size_t algo = 0; algo < algorithms.size(); ++algo)
		{
			curr_file = "..//solutions//" + folder + "//";
			curr_file.append("s_");
			curr_file.append(algorithms[algo]);
			curr_file.append("_");
			curr_file.append(instance);

			// std::cout << curr_file << std::endl;

			std::fstream input;
			input.open(curr_file.c_str(), std::fstream::in);

			if (!input.is_open())
			{
				std::cout << "Could not open file " << curr_file << std::endl;
				throw 4;
			}

			std::stringstream s_lb, s_ub, s_time;
			std::string status;
			double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0;
			std::string line;

			getline(input, line);
			size_t pos = line.find_first_of(":");
			status = line.substr(pos + 2);

			std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
			status.erase(end_pos, status.end());

			getline(input, line);
			getline(input, line);
			pos = line.find_first_of(":");
			s_lb << line.substr(pos + 2);
			if (s_lb.str() == "-inf")
				lb = -1;
			else
				s_lb >> lb;

			getline(input, line);
			pos = line.find_first_of(":");
			s_ub << line.substr(pos + 2);
			if (s_ub.str() == "inf")
				ub = -1;
			else
				s_ub >> ub;

			getline(input, line);
			getline(input, line);
			pos = line.find_first_of(":");
			s_time << line.substr(pos + 2);
			s_time >> time;

			if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
			{
				if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
				{
					if ((!double_greater(ub, lb)))
					{
						++(total_num_optimal[algo]);
					}
					else
						gap = (100.0 * (ub - lb)) / ub;
				}
				else
					gap = 100.0;
			}
			else
			{
				++(total_num_optimal[algo]);
			}

			total_time[algo] += time;

			// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

			if (!double_equals(gap, 0.0))
			{
				total_avg_gap[algo] += gap;
			}
			input.close();

			output << gap << ";" << time << ";";
		}

		output << std::endl;
	}

	output << "Total;";
	for (size_t algo = 0; algo < algorithms.size(); ++algo)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		if (double_equals(total_time[algo], 0.0))
			total_time[algo] = -1;
		total_avg_gap[algo] /= instances.size();

		output << total_avg_gap[algo] << ";" << total_time[algo] << ";";
	}

	output << std::endl;

	output.close();
}

void GenerateAlgorithmsLatexTable(std::string folder)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("vehc_seq_model_sym_break_reform");
	// algorithms.push_back("cb_vehc_seq_model_sym_break_reform");

	algorithms.push_back("vehc_slot_model_sym_break_reform");
	// algorithms.push_back("cb_vehc_slot_model_sym_break_reform");

	// algorithms.push_back("itm_seq_model_sym_break_reform");
	// algorithms.push_back("cb_itm_seq_model_sym_break_reform");

	// algorithms.push_back("relax_item_seq_model_sym_break");
	// algorithms.push_back("relax_item_seq_model_sym_break_reform");
	// algorithms.push_back("relax_cb_item_seq_model_sym_break");
	// algorithms.push_back("relax_cb_item_seq_model_sym_break_reform");

	// std::cout << output_name << std::endl;

	std::vector<std::string> instances = {
		"data10-6-37",
		"data20-20-6",
		"data30-30-32",
		// // "data40-40-3",
		"data5-10-36",
		"data5-4-4",
		"data5-6-10",
		"data5-7-37",
		"data5-9-26",
		"data6-6-19"};
	const std::vector<double> percentage_of_items_serviced_by_fleet_vec = {0.25, 0.5};
	const std::vector<int> number_of_item_groups_vec = {2, 5, 10};

	std::fstream output;
	std::string output_name;
	output_name = "..//tables//latex//table_algorithms.txt";
	output.open(output_name.c_str(), std::fstream::out);

	output << std::setprecision(2) << std::fixed;

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(), 0.0);
	std::vector<double> total_avg_gap(algorithms.size(), 0.0);
	std::vector<int> total_num_optimal(algorithms.size(), 0);

	size_t total_num_instances = instances.size() * percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();

	for (auto instance : instances)
	{
		std::vector<std::vector<double>> time_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<double> avg_time_inst_size(algorithms.size(), 0.0);
		std::vector<double> avg_gap_inst_size(algorithms.size(), 0.0);
		std::vector<double> st_dev_inst_size(algorithms.size(), 0.0);

		std::vector<int> num_optimal_inst_size(algorithms.size(), 0);

		for (auto percentage_of_items_serviced_by_fleet : percentage_of_items_serviced_by_fleet_vec)
		{
			for (auto number_of_item_groups : number_of_item_groups_vec)
			{

				std::cout << std::endl;
				// Create an ostringstream object
				std::ostringstream oss;
				// Set the precision and fixed-point notation
				oss << std::fixed << std::setprecision(2) << percentage_of_items_serviced_by_fleet;
				std::string instance_name = instance + "-c" + oss.str() + "-p" + std::to_string(number_of_item_groups);
				std::cout << instance_name << std::endl;

				// std::cout << instance << std::endl;
				double original_lp = 0.0;
				for (size_t algo = 0; algo < algorithms.size(); ++algo)
				{
					curr_file = "..//solutions//" + folder + "//s_" + algorithms[algo] + "_" + instance_name + ".txt";

					// std::cout << curr_file << std::endl;

					std::fstream input;
					input.open(curr_file.c_str(), std::fstream::in);

					if (!input.is_open())
					{
						std::cout << "Could not open file " << curr_file << std::endl;
						throw 4;
					}

					std::stringstream s_lb, s_ub, s_time;
					std::string status;
					double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0;
					std::string line;

					getline(input, line);
					size_t pos = line.find_first_of(":");
					status = line.substr(pos + 2);

					std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
					status.erase(end_pos, status.end());

					getline(input, line);
					getline(input, line);
					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_lb << line.substr(pos + 2);
					if (s_lb.str() == "-inf")
						lb = -1;
					else
						s_lb >> lb;

					getline(input, line);
					pos = line.find_first_of(":");
					s_ub << line.substr(pos + 2);
					if (s_ub.str() == "inf")
						ub = -1;
					else
						s_ub >> ub;

					getline(input, line);
					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_time << line.substr(pos + 2);
					s_time >> time;

					time = std::min(time, 3600.0);

					if (status == "OPTIMAL")
						std::cout << lb << " ";
					// if ((status == "OPTIMAL") || (status == "INFEASIBLE"))
					// 	std::cout << instance << std::endl;

					if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
					{
						if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
						{
							if ((!double_greater(ub, lb)))
							{
								++(num_optimal_inst_size[algo]);
								++(total_num_optimal[algo]);
								// time_per_algo_inst_size[algo].push_back(time);
								// time_per_algo_quantile[algo].push_back(time);
								// total_time_per_algo[algo].push_back(time);

								// total_avg_time[algo] += time;
								// avg_time_inst_size[algo] += time;
								// avg_time_quantile[algo] += time;
							}
							else
								gap = (100.0 * (ub - lb)) / ub;
						}
						else
							gap = 100.0;
					}
					else
					{
						++(num_optimal_inst_size[algo]);
						++(total_num_optimal[algo]);
						// time_per_algo_inst_size[algo].push_back(time);
						// time_per_algo_quantile[algo].push_back(time);
						// total_time_per_algo[algo].push_back(time);

						// total_avg_time[algo] += time;
						// avg_time_inst_size[algo] += time;
						// avg_time_quantile[algo] += time;
					}

					time_per_algo_inst_size[algo].push_back(time);
					total_time_per_algo[algo].push_back(time);

					total_avg_time[algo] += time;
					avg_time_inst_size[algo] += time;

					// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

					// if (!double_equals(gap, 0.0))
					// {
					gap_per_algo_inst_size[algo].push_back(gap);
					total_gap_per_algo[algo].push_back(gap);

					total_avg_gap[algo] += gap;
					avg_gap_inst_size[algo] += gap;
					// }
					input.close();
				}

				// getchar();getchar();
			}
		}
		output << instance << " & &";

		int num_inst_per_vertex_size = percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();

		for (size_t algo = 0; algo < algorithms.size(); ++algo)
		{
			if ((time_per_algo_inst_size[algo]).size() > 0)
				avg_time_inst_size[algo] /= (1.0 * ((time_per_algo_inst_size[algo]).size()));
			else
				avg_time_inst_size[algo] = -1;
			if ((gap_per_algo_inst_size[algo]).size() > 0)
			{
				avg_gap_inst_size[algo] /= (1.0 * ((gap_per_algo_inst_size[algo]).size()));
				st_dev_inst_size[algo] = StDev(gap_per_algo_inst_size[algo], avg_gap_inst_size[algo]);
			}
			else
				avg_gap_inst_size[algo] = st_dev_inst_size[algo] = -1;
			output << " & & " << num_optimal_inst_size[algo] << "/" << num_inst_per_vertex_size << " & " << avg_time_inst_size[algo] << " & " << avg_gap_inst_size[algo] << " & " << st_dev_inst_size[algo];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total & &";
	for (size_t algo = 0; algo < algorithms.size(); ++algo)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		if ((total_time_per_algo[algo]).size() > 0)
			total_avg_time[algo] /= (1.0 * ((total_time_per_algo[algo]).size()));
		else
			total_avg_time[algo] = -1;
		if ((total_gap_per_algo[algo]).size() > 0)
			total_avg_gap[algo] /= (1.0 * ((total_gap_per_algo[algo]).size()));
		else
			total_avg_gap[algo] = -1;
		output << " & & " << total_num_optimal[algo] << "/" << total_num_instances << " & " << total_avg_time[algo] << " & " << total_avg_gap[algo] << " & " << StDev(total_gap_per_algo[algo], total_avg_gap[algo]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateHeuristicsLatexTable(std::string folder_exact, std::string folder_heuristic, bool add_exact_results)
{
	std::string curr_file;
	std::vector<std::pair<int, std::string>> algorithms;
	std::vector<std::string> exact_algorithms;

	exact_algorithms.push_back("vehc_seq_model_sym_break_reform");
	// algorithms.push_back("cb_vehc_seq_model_sym_break_reform");

	exact_algorithms.push_back("vehc_slot_model_sym_break_reform");
	// algorithms.push_back("cb_vehc_slot_model_sym_break_reform");

	// algorithms.push_back("relax_item_seq_model_sym_break");
	// algorithms.push_back("relax_item_seq_model_sym_break_reform");
	// algorithms.push_back("relax_cb_item_seq_model_sym_break");
	// algorithms.push_back("relax_cb_item_seq_model_sym_break_reform");

	// std::cout << output_name << std::endl;

	std::vector<std::string> instances = {
		"data10-6-37",
		"data20-20-6",
		"data30-30-32",
		// // "data40-40-3",
		"data5-10-36",
		"data5-4-4",
		"data5-6-10",
		"data5-7-37",
		"data5-9-26",
		"data6-6-19"};
	const std::vector<double> percentage_of_items_serviced_by_fleet_vec = {0.25, 0.5};
	const std::vector<int> number_of_item_groups_vec = {2, 5, 10};

	// algorithms.push_back("baseline_ks_b5_[84,19]_d0.96_feas");
	algorithms.push_back(std::pair<int, std::string>(0, "ks_vehc_seq_model_sym_break_reform_[90,45]_d0.90"));
	algorithms.push_back(std::pair<int, std::string>(0, "ks_vehc_slot_model_sym_break_reform_[90,45]_d0.90"));
	algorithms.push_back(std::pair<int, std::string>(0, "ks_vehc_seq_model_sym_break_reform_cluster_by_group_[90,45]_d0.90"));
	algorithms.push_back(std::pair<int, std::string>(0, "ks_vehc_slot_model_sym_break_reform_cluster_by_group_[90,45]_d0.90"));

	std::fstream output;
	std::string output_name;
	output_name = "..//tables//latex//table_heuristics.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo_exact(exact_algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo_exact(exact_algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_time_exact(exact_algorithms.size(), 0.0);
	std::vector<double> total_avg_gap_exact(exact_algorithms.size(), 0.0);
	std::vector<int> total_num_optimal_exact(exact_algorithms.size(), 0);
	size_t total_num_instances_exact = instances.size() * percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(), 0.0);
	std::vector<double> total_avg_gap(algorithms.size(), 0.0);
	std::vector<int> total_num_best_known_bound(algorithms.size(), 0);

	size_t total_num_instances = instances.size() * percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();
	size_t num_inst_per_vertex_size_exact = percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();
	size_t num_inst_per_vertex_size = percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();

	for (auto instance : instances)
	{
		std::vector<std::vector<double>> time_per_algo_inst_size_exact(exact_algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo_inst_size_exact(exact_algorithms.size(), std::vector<double>());
		std::vector<double> avg_time_inst_size_exact(exact_algorithms.size(), 0.0);
		std::vector<double> avg_gap_inst_size_exact(exact_algorithms.size(), 0.0);
		std::vector<double> st_dev_inst_size_exact(exact_algorithms.size(), 0.0);

		std::vector<int> num_optimal_inst_size_exact(exact_algorithms.size(), 0);

		std::vector<std::vector<double>> time_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<double> avg_time_inst_size(algorithms.size(), 0.0);
		std::vector<double> avg_gap_inst_size(algorithms.size(), 0.0);
		std::vector<double> st_dev_inst_size(algorithms.size(), 0.0);

		std::vector<int> num_best_known_bound_inst_size(algorithms.size(), 0);

		for (auto percentage_of_items_serviced_by_fleet : percentage_of_items_serviced_by_fleet_vec)
		{
			for (auto number_of_item_groups : number_of_item_groups_vec)
			{

				// Create an ostringstream object
				std::ostringstream oss;
				// Set the precision and fixed-point notation
				oss << std::fixed << std::setprecision(2) << percentage_of_items_serviced_by_fleet;
				std::string instance_name = instance + "-c" + oss.str() + "-p" + std::to_string(number_of_item_groups);
				std::cout << instance_name << std::endl;

				double best_lb = -1;
				double best_num_unproductive_moves = std::numeric_limits<double>::infinity();
				bool is_infeasible = false;
				bool has_optimal_bound = false;

				// compute the best known primal bound among the exact algorithms.
				for (size_t algo = 0; algo < exact_algorithms.size(); ++algo)
				{
					curr_file = "..//solutions//" + folder_exact + "//s_" + exact_algorithms[algo] + "_" + instance_name + ".txt";

					// std::cout << curr_file << std::endl;

					std::fstream input;
					input.open(curr_file.c_str(), std::fstream::in);

					if (!input.is_open())
					{
						std::cout << "Could not open file " << curr_file << std::endl;
						throw 4;
					}

					std::stringstream s_lb, s_ub, s_time, s_unproductive_moves;
					std::string status;
					double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0, num_unproductive_moves = 0.0;
					std::string line;

					getline(input, line);
					size_t pos = line.find_first_of(":");
					status = line.substr(pos + 2);

					std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
					status.erase(end_pos, status.end());

					if (status == "INFEASIBLE")
					{
						is_infeasible = true;
						break;
					}

					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_unproductive_moves << line.substr(pos + 2);
					s_unproductive_moves >> num_unproductive_moves;
					// std::cout << line << std::endl;

					if (double_less(num_unproductive_moves, best_num_unproductive_moves))
					{
						// std::cout << num_unproductive_moves << " < " << best_num_unproductive_moves << std::endl;
						best_num_unproductive_moves = num_unproductive_moves;
					}

					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_lb << line.substr(pos + 2);
					if (s_lb.str() == "-inf")
						lb = -1;
					else
						s_lb >> lb;

					// std::cout << exact_algorithms[algo] << " " << lb << std::endl;

					getline(input, line);
					pos = line.find_first_of(":");
					s_ub << line.substr(pos + 2);
					if (s_ub.str() == "inf")
						ub = -1;
					else
						s_ub >> ub;

					getline(input, line);
					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_time << line.substr(pos + 2);
					s_time >> time;

					if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
					{
						if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
						{
							if (double_greater(ub, lb))
								gap = (100.0 * (ub - lb)) / ub;
						}
						else
						{
							gap = 100.0;
						}
					}

					// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;
					if (double_equals(gap, 0.0))
					{
						has_optimal_bound = true;
						best_lb = lb;
						break;
					}
					else if (double_greater(lb, best_lb))
						best_lb = lb;

					input.close();
				}

				// std::cout << best_lb << std::endl;

				// compute information of the exact algorithms.
				double original_lp = 0.0;
				for (size_t algo = 0; algo < exact_algorithms.size(); ++algo)
				{
					curr_file = "..//solutions//" + folder_exact + "//s_" + exact_algorithms[algo] + "_" + instance_name + ".txt";

					// std::cout << curr_file << std::endl;

					std::fstream input;
					input.open(curr_file.c_str(), std::fstream::in);

					if (!input.is_open())
					{
						std::cout << "Could not open file " << curr_file << std::endl;
						throw 4;
					}

					std::stringstream s_lb, s_ub, s_time;
					std::string status;
					double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0;
					std::string line;

					getline(input, line);
					size_t pos = line.find_first_of(":");
					status = line.substr(pos + 2);

					std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
					status.erase(end_pos, status.end());

					getline(input, line);
					getline(input, line);
					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_lb << line.substr(pos + 2);
					if (s_lb.str() == "-inf")
						lb = -1;
					else
						s_lb >> lb;

					getline(input, line);
					pos = line.find_first_of(":");
					s_ub << line.substr(pos + 2);
					if (s_ub.str() == "inf")
						ub = -1;
					else
						s_ub >> ub;

					getline(input, line);
					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_time << line.substr(pos + 2);
					s_time >> time;

					// if ((status == "OPTIMAL") || (status == "INFEASIBLE"))
					// 	std::cout << instance << std::endl;

					if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
					{
						if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
						{
							if ((!double_greater(ub, lb)))
							{
								++(num_optimal_inst_size_exact[algo]);
								++(total_num_optimal_exact[algo]);
								// time_per_algo_inst_size_exact[algo].push_back(time);
								// time_per_algo_quantile_exact[algo].push_back(time);
								// total_time_per_algo_exact[algo].push_back(time);

								// total_avg_time_exact[algo] += time;
								// avg_time_inst_size_exact[algo] += time;
								// avg_time_quantile_exact[algo] += time;
							}
							else
								gap = (100.0 * (ub - lb)) / ub;
						}
						else
							gap = 100.0;
					}
					else
					{
						++(num_optimal_inst_size_exact[algo]);
						++(total_num_optimal_exact[algo]);
						// time_per_algo_inst_size_exact[algo].push_back(time);
						// time_per_algo_quantile_exact[algo].push_back(time);
						// total_time_per_algo_exact[algo].push_back(time);

						// total_avg_time_exact[algo] += time;
						// avg_time_inst_size_exact[algo] += time;
						// avg_time_quantile_exact[algo] += time;
					}

					time_per_algo_inst_size_exact[algo].push_back(time);
					total_time_per_algo_exact[algo].push_back(time);

					total_avg_time_exact[algo] += time;
					avg_time_inst_size_exact[algo] += time;

					// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

					// if (!double_equals(gap, 0.0))
					{
						gap_per_algo_inst_size_exact[algo].push_back(gap);
						total_gap_per_algo_exact[algo].push_back(gap);

						total_avg_gap_exact[algo] += gap;
						avg_gap_inst_size_exact[algo] += gap;
					}
					input.close();
				}

				for (size_t algo = 0; algo < algorithms.size(); ++algo)
				{
					int num_seeds = algorithms[algo].first;

					if (num_seeds == 0)
					{
						curr_file = "..//solutions//" + folder_heuristic + "//s_" + algorithms[algo].second + "_" + instance_name + ".txt";

						// std::cout << curr_file << std::endl;

						std::fstream input;
						input.open(curr_file.c_str(), std::fstream::in);

						if (!input.is_open())
						{
							std::cout << "Could not open file " << curr_file << std::endl;
							continue;
							// throw 4;
						}
						else
						{
							std::cout << instance << " " << algorithms[algo].second << std::endl;
						}

						std::stringstream s_lb1, s_t1, s_max_improve_iter, s_num_unproductive_moves;
						std::string line, status;
						double heuristic_lb = -1, heuristic_time = 0, improvement = 0, num_unproductive_moves = 0;
						getline(input, line);
						size_t pos = line.find_first_of(":");
						status = line.substr(pos + 2);

						std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
						status.erase(end_pos, status.end());

						// std::cout << status << std::endl;
						getline(input, line);
						getline(input, line);

						pos = line.find_first_of(":");
						s_t1 << line.substr(pos + 2);
						s_t1 >> heuristic_time;

						if (status != "INFEASIBLE")
						{
							getline(input, line);
							pos = line.find_first_of(":");
							s_lb1 << line.substr(pos + 2);
							s_lb1 >> heuristic_lb;
							heuristic_lb = round_decimals(heuristic_lb, 2); // IMPORTANT because I saved heuristic solutions file without 2 decimal precision.

							getline(input, line);
							getline(input, line);
							pos = line.find_first_of(":");
							s_num_unproductive_moves << line.substr(pos + 2);
							s_num_unproductive_moves >> num_unproductive_moves;
							num_unproductive_moves = round_decimals(num_unproductive_moves, 2); // IMPORTANT because I saved heuristic solutions file without 2 decimal precision.
						}

						// if (double_greater(heuristic_time, 1500))
						// std::cout << heuristic_time << std::endl;

						time_per_algo_inst_size[algo].push_back(heuristic_time);
						total_time_per_algo[algo].push_back(heuristic_time);

						total_avg_time[algo] += heuristic_time;
						avg_time_inst_size[algo] += heuristic_time;

						// if ((status == "INFEASIBLE") && !is_infeasible)
						// 	std::cout << "Heuristica provou inviabilidade nova" << std::endl;

						if (((status == "INFEASIBLE") && is_infeasible) || double_equals(best_lb, heuristic_lb))
						{
							++(num_best_known_bound_inst_size[algo]);
							++(total_num_best_known_bound[algo]);
							improvement = 0;
						}
						else
						{
							if ((status == "INFEASIBLE") && !is_infeasible)
							{
								++(num_best_known_bound_inst_size[algo]);

								++(total_num_best_known_bound[algo]);
								improvement = 100;
							}
							else if (((double_equals(heuristic_lb, -1) || double_equals(heuristic_lb, 0)) && !double_equals(best_lb, -1)) || (is_infeasible && (status != "INFEASIBLE")))
							{
								improvement = -100; // means that the exact found a bound (or proved infeasibility) and the heuristic found nothing or zero bound.
							}
							else if ((double_equals(best_lb, -1) || double_equals(best_lb, 0)) && !double_less(heuristic_lb, 0)) // also make improvement = 100 if best_lb = 0 and heuristic_lb is >= 0.
							{
								improvement = 100.0;
								++(num_best_known_bound_inst_size[algo]);
								++(total_num_best_known_bound[algo]);
							}
							else
							{
								improvement = (100 * (heuristic_lb - best_lb)) / best_lb;
								// improvement = (100 * (num_unproductive_moves - best_num_unproductive_moves)) / best_num_unproductive_moves;
								// std::cout << heuristic_lb << " " << best_lb << " " << status << std::endl;
								// std::cout << num_unproductive_moves << " " << best_num_unproductive_moves << " " << status << std::endl;
								if (!double_less(improvement, 0.0))
								{
									++(num_best_known_bound_inst_size[algo]);
									++(total_num_best_known_bound[algo]);
								}
							}
						}

						// if (!double_equals(improvement, 100))
						// {
						// 	std::cout << improvement << std::endl;
						// 	getchar();
						// 	getchar();
						// }
						gap_per_algo_inst_size[algo].push_back(improvement);
						total_gap_per_algo[algo].push_back(improvement);

						avg_gap_inst_size[algo] += improvement;
						total_avg_gap[algo] += improvement;

						input.close();
					}
					else
					{
						double avg_all_seeds_heuristic_time = 0.0;
						double avg_all_seeds_improvement = 0.0;
						double avg_all_num_best_known_bound_inst_size = 0.0;
						double avg_all_num_best_known_bound_quantile = 0.0;
						double avg_all_total_num_best_known_bound = 0.0;
						double avg_all_improvement = 0.0;

						for (int seed = 1; seed <= num_seeds; ++seed)
						{
							curr_file = "..//solutions//" + folder_heuristic + "//";
							curr_file.append("s_");
							curr_file.append(algorithms[algo].second);
							curr_file.append("_seed_");
							curr_file.append(std::to_string(seed));
							curr_file.append("_");
							curr_file.append(instance);

							// std::cout << curr_file << std::endl;

							std::fstream input;
							input.open(curr_file.c_str(), std::fstream::in);

							if (!input.is_open())
							{
								std::cout << "Could not open file " << curr_file << std::endl;
								continue;
								// throw 4;
							}
							else
							{
								std::cout << instance << " " << algorithms[algo].second << std::endl;
							}

							std::stringstream s_lb1, s_t1, s_max_improve_iter;
							std::string line, status;
							double heuristic_lb = -1, heuristic_time = 0, improvement = 0;

							getline(input, line);
							size_t pos = line.find_first_of(":");
							status = line.substr(pos + 2);

							std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
							status.erase(end_pos, status.end());

							getline(input, line);
							pos = line.find_first_of(":");
							s_lb1 << line.substr(pos + 2);
							s_lb1 >> heuristic_lb;
							heuristic_lb = round_decimals(heuristic_lb, 2); // IMPORTANT because I saved heuristic solutions file without 2 decimal precision.

							// std::cout << status << std::endl;
							getline(input, line);
							getline(input, line);
							getline(input, line);

							pos = line.find_first_of(":");
							s_t1 << line.substr(pos + 2);
							s_t1 >> heuristic_time;

							avg_all_seeds_heuristic_time += heuristic_time;

							if (((status == "INFEASIBLE") && is_infeasible) || double_equals(best_lb, heuristic_lb))
							{
								avg_all_num_best_known_bound_inst_size += 1.0;
								avg_all_num_best_known_bound_quantile += 1.0;
								avg_all_total_num_best_known_bound += 1.0;
								improvement = 0;
							}
							else
							{
								if ((status == "INFEASIBLE") && !is_infeasible)
								{
									avg_all_num_best_known_bound_inst_size += 1.0;
									avg_all_num_best_known_bound_quantile += 1.0;
									avg_all_total_num_best_known_bound += 1.0;
									improvement = 100.0;
								}
								else if (((double_equals(heuristic_lb, -1) || double_equals(heuristic_lb, 0)) && !double_equals(best_lb, -1)) || (is_infeasible && (status != "INFEASIBLE")))
								{
									improvement = -100.0; // means that the exact found a bound (or proved infeasibility) and the heuristic found nothing or zero bound.
								}
								else if ((double_equals(best_lb, -1) || double_equals(best_lb, 0)) && !double_less(heuristic_lb, 0)) // also make improvement = 100 if best_lb = 0 and heuristic_lb is >= 0.
								{
									improvement = 100.0;
									avg_all_num_best_known_bound_inst_size += 1.0;
									avg_all_num_best_known_bound_quantile += 1.0;
									avg_all_total_num_best_known_bound += 1.0;
								}
								else
								{
									improvement = (100 * (heuristic_lb - best_lb)) / best_lb;
									std::cout << best_lb << " " << heuristic_lb << " " << status << std::endl;
									if (!double_less(improvement, 0.0))
									{
										avg_all_num_best_known_bound_inst_size += 1.0;
										avg_all_num_best_known_bound_quantile += 1.0;
										avg_all_total_num_best_known_bound += 1.0;
									}
								}
							}

							avg_all_seeds_improvement += improvement;

							input.close();
						}

						avg_all_seeds_heuristic_time /= num_seeds;
						time_per_algo_inst_size[algo].push_back(avg_all_seeds_heuristic_time);
						total_time_per_algo[algo].push_back(avg_all_seeds_heuristic_time);

						total_avg_time[algo] += avg_all_seeds_heuristic_time;
						avg_time_inst_size[algo] += avg_all_seeds_heuristic_time;

						avg_all_seeds_improvement /= num_seeds;
						gap_per_algo_inst_size[algo].push_back(avg_all_seeds_improvement);
						total_gap_per_algo[algo].push_back(avg_all_seeds_improvement);

						avg_gap_inst_size[algo] += avg_all_seeds_improvement;
						total_avg_gap[algo] += avg_all_seeds_improvement;

						avg_all_num_best_known_bound_inst_size /= num_seeds;
						avg_all_num_best_known_bound_quantile /= num_seeds;
						avg_all_total_num_best_known_bound /= num_seeds;

						num_best_known_bound_inst_size[algo] += 1.0;
						total_num_best_known_bound[algo] += 1.0;
					}

					// std::cout << algorithms[algo] << " " << heuristic_lb << " " << heuristic_time << " " << improvement << std::endl;
				}
			}
		}
		output << instance << " & &";

		if (add_exact_results)
		{
			for (size_t algo = 0; algo < exact_algorithms.size(); ++algo)
			{
				if ((time_per_algo_inst_size_exact[algo]).size() > 0)
					avg_time_inst_size_exact[algo] /= (1.0 * ((time_per_algo_inst_size_exact[algo]).size()));
				else
					avg_time_inst_size_exact[algo] = -1;
				if ((gap_per_algo_inst_size_exact[algo]).size() > 0)
				{
					avg_gap_inst_size_exact[algo] /= (1.0 * ((gap_per_algo_inst_size_exact[algo]).size()));
					st_dev_inst_size_exact[algo] = StDev(gap_per_algo_inst_size_exact[algo], avg_gap_inst_size_exact[algo]);
				}
				else
					avg_gap_inst_size_exact[algo] = st_dev_inst_size_exact[algo] = -1;
				output << " & & " << num_optimal_inst_size_exact[algo] << "/" << num_inst_per_vertex_size_exact << " & " << avg_time_inst_size_exact[algo] << " & " << avg_gap_inst_size_exact[algo] << " & " << st_dev_inst_size_exact[algo];
			}
		}

		for (size_t algo = 0; algo < algorithms.size(); ++algo)
		{
			if ((time_per_algo_inst_size[algo]).size() > 0)
				avg_time_inst_size[algo] /= (1.0 * ((time_per_algo_inst_size[algo]).size()));
			else
				avg_time_inst_size[algo] = -1;
			if ((gap_per_algo_inst_size[algo]).size() > 0)
			{
				avg_gap_inst_size[algo] /= (1.0 * ((gap_per_algo_inst_size[algo]).size()));
				st_dev_inst_size[algo] = StDev(gap_per_algo_inst_size[algo], avg_gap_inst_size[algo]);
			}
			else
				avg_gap_inst_size[algo] = st_dev_inst_size[algo] = -1;
			output << " & & " << num_best_known_bound_inst_size[algo] << "/" << num_inst_per_vertex_size << " & " << avg_time_inst_size[algo] << " & " << avg_gap_inst_size[algo] << " & " << st_dev_inst_size[algo];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total & &";

	if (add_exact_results)
	{
		for (size_t algo = 0; algo < exact_algorithms.size(); ++algo)
		{
			// std::cout << total_improvement_per_algo[j].size() << std::endl;
			if ((total_time_per_algo_exact[algo]).size() > 0)
				total_avg_time_exact[algo] /= (1.0 * ((total_time_per_algo_exact[algo]).size()));
			else
				total_avg_time_exact[algo] = -1;
			if ((total_gap_per_algo_exact[algo]).size() > 0)
				total_avg_gap_exact[algo] /= (1.0 * ((total_gap_per_algo_exact[algo]).size()));
			else
				total_avg_gap_exact[algo] = -1;
			output << "& & " << total_num_optimal_exact[algo] << "/" << total_num_instances_exact << " & " << total_avg_time_exact[algo] << " & " << total_avg_gap_exact[algo] << " & " << StDev(total_gap_per_algo_exact[algo], total_avg_gap_exact[algo]);
		}
	}

	for (size_t algo = 0; algo < algorithms.size(); ++algo)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		if ((total_time_per_algo[algo]).size() > 0)
			total_avg_time[algo] /= (1.0 * ((total_time_per_algo[algo]).size()));
		else
			total_avg_time[algo] = -1;
		if ((total_gap_per_algo[algo]).size() > 0)
			total_avg_gap[algo] /= (1.0 * ((total_gap_per_algo[algo]).size()));
		else
			total_avg_gap[algo] = -1;
		output << "& & " << total_num_best_known_bound[algo] << "/" << total_num_instances << " & " << total_avg_time[algo] << " & " << total_avg_gap[algo] << " & " << StDev(total_gap_per_algo[algo], total_avg_gap[algo]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateLPImprovementsLatexTable(std::string folder)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_vehc_seq_model_sym_break");
	algorithms.push_back("relax_vehc_seq_model_sym_break_reform");
	algorithms.push_back("relax_cb_vehc_seq_model_sym_break");
	algorithms.push_back("relax_cb_vehc_seq_model_sym_break_reform");

	algorithms.push_back("relax_vehc_slot_model_sym_break");
	algorithms.push_back("relax_vehc_slot_model_sym_break_reform");
	algorithms.push_back("relax_cb_vehc_slot_model_sym_break");
	algorithms.push_back("relax_cb_vehc_slot_model_sym_break_reform");

	// algorithms.push_back("relax_itm_seq_model_sym_break");
	// algorithms.push_back("relax_itm_seq_model_sym_break_reform");
	// algorithms.push_back("relax_cb_itm_seq_model_sym_break");
	// algorithms.push_back("relax_cb_itm_seq_model_sym_break_reform");

	std::fstream output;
	std::string output_name = "..//tables//latex//table_LP_improvements.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::string> instances = {
		"data10-6-37",
		"data20-20-6",
		"data30-30-32",
		// // "data40-40-3",
		"data5-10-36",
		"data5-4-4",
		"data5-6-10",
		"data5-7-37",
		"data5-9-26",
		"data6-6-19"};
	const std::vector<double> percentage_of_items_serviced_by_fleet_vec = {0.25, 0.5};
	const std::vector<int> number_of_item_groups_vec = {2, 5, 10};

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(), 0.0);
	int total_num_instances = instances.size() * percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();

	for (auto instance : instances)
	{
		std::vector<std::vector<double>> improvement_per_algo(algorithms.size(), std::vector<double>());
		std::vector<double> avg_improvement(algorithms.size(), 0.0);
		std::vector<double> st_dev(algorithms.size(), 0.0);
		for (auto percentage_of_items_serviced_by_fleet : percentage_of_items_serviced_by_fleet_vec)
		{
			for (auto number_of_item_groups : number_of_item_groups_vec)
			{

				// Create an ostringstream object
				std::ostringstream oss;
				// Set the precision and fixed-point notation
				oss << std::fixed << std::setprecision(2) << percentage_of_items_serviced_by_fleet;
				std::string instance_name = instance + "-c" + oss.str() + "-p" + std::to_string(number_of_item_groups);
				std::cout << instance_name << std::endl;

				double original_lp = 0.0;

				for (size_t algo = 0; algo < algorithms.size(); ++algo)
				{
					double curr_improvement = 0.0;
					curr_file = "..//solutions//" + folder + "//s_" + algorithms[algo] + "_" + instance_name + ".txt";

					// std::cout << curr_file << std::endl;

					std::fstream input;
					input.open(curr_file.c_str(), std::fstream::in);

					if (!input.is_open())
					{
						std::cout << "Could not open file " << curr_file << std::endl;
						continue;
					}

					std::stringstream s_lp;
					std::string status;
					double lp = 0.0;
					std::string line;

					getline(input, line);
					size_t pos = line.find_first_of(":");
					status = line.substr(pos + 2);

					std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
					status.erase(end_pos, status.end());
					status = "FEASIBLE";

					getline(input, line);
					getline(input, line);
					getline(input, line);
					pos = line.find_first_of(":");
					s_lp << line.substr(pos + 2);

					if (s_lp.str() == "inf")
						lp = -1;
					else
						s_lp >> lp;

					if (algo == 0)
					{
						if ((double_equals(lp, -1)) || (status == "INFEASIBLE"))
						{
							original_lp = -1.0;
						}
						else
							original_lp = lp;
						curr_improvement = 0.0;
					}
					else
					{
						if ((double_equals(lp, -1)) || (status == "INFEASIBLE") || (double_equals(original_lp, -1)))
						{
							curr_improvement = -1;
						}
						else
						{
							if (double_equals(original_lp, 0.0))
								curr_improvement = 0.0;
							else
								curr_improvement = (100 * (original_lp - lp)) / original_lp;

							if (double_less(curr_improvement, 0))
								std::cout << original_lp << " - " << lp << std::endl;
							// if (!double_greater(original_lp, lp))
							// 	std::cout << original_lp << " - " << lp << std::endl;
						}
					}

					// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

					if (!double_equals(curr_improvement, -1))
					{
						improvement_per_algo[algo].push_back(curr_improvement);
						avg_improvement[algo] += curr_improvement;
						total_improvement_per_algo[algo].push_back(curr_improvement);
						total_avg_improvement[algo] += curr_improvement;
					}
					input.close();

					// getchar();getchar();
				}
			}
		}
		output << instance << " & & " << percentage_of_items_serviced_by_fleet_vec.size() * number_of_item_groups_vec.size();

		for (size_t algo = 1; algo < algorithms.size(); ++algo)
		{
			if (improvement_per_algo[algo].size() > 0)
				avg_improvement[algo] /= (1.0 * ((improvement_per_algo[algo]).size()));
			else
				avg_improvement[algo] = -1;
			st_dev[algo] = StDev(improvement_per_algo[algo], avg_improvement[algo]);
			output << " & & " << avg_improvement[algo] << " & " << st_dev[algo];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total & & " << total_num_instances;
	for (size_t algo = 1; algo < algorithms.size(); algo++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		if ((total_improvement_per_algo[algo]).size() > 0)
			total_avg_improvement[algo] /= (1.0 * ((total_improvement_per_algo[algo]).size()));
		else
			total_avg_improvement[algo] = -1;
		output << " & & " << total_avg_improvement[algo] << " & " << StDev(total_improvement_per_algo[algo], total_avg_improvement[algo]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

int main()
{
	// std::string folder = "2024-07-26_09:18:48_all_relax_new";
	// std::string folder = "2024-06-23_13:01:07_all_kernel_search_less_time";
	std::string folder_heuristic_sol = "2025-09-22_01:52:09_ks_more_time";
	std::string folder_exact_sol = "2025-09-20_01:13:01_exact";
	std::string folder_relax = "2025-09-19_19:04:13";
	// try
	// {
	// GenerateAlgorithmsLatexTablePerInstance(folder);
	// return 1;
	// GenerateLPImprovementsLatexTable(folder_relax);
	GenerateAlgorithmsLatexTable(folder_exact_sol);
	GenerateHeuristicsLatexTable(folder_exact_sol, folder_heuristic_sol, false);

	// namespace fs = std::filesystem;

	// std::string path = "./"; // Change this to your target folder

	// for (const auto &entry : fs::directory_iterator(path))
	// {
	// 	if (entry.is_regular_file())
	// 	{
	// 		auto filePath = entry.path();
	// 		if (filePath.extension() == ".txt")
	// 		{
	// 			std::string filename = filePath.stem().string(); // stem() removes extension
	// 			std::cout << filename << std::endl;
	// 		}
	// 	}
	// }

	return 0;
}
