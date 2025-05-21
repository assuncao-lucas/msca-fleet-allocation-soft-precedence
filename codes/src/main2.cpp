#include <iostream>
#include <string>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "src/instance.h"
#include "src/graph.h"
#include "src/general.h"
#include "src/exact/formulations.h"
#include "src/timer.h"
#include "src/graph_algorithms.h"
#include "src/heuristic_solution.h"
#include "src/feasibility_pump/feasibility_pump.h"
#include "src/ALNS/ALNS.h"
#include "src/kernel_search/kernel_search.h"
#include "src/initial_solution/initial_solution.h"
#include "src/simulated_annealing/simulated_annealing.h"

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
	std::vector<std::string> instance_sizes;
	std::vector<std::string> instance_types;
	std::vector<std::string> instance_limit_quantiles;

	std::vector<std::string> num_vehicles_vec{"2", "4"};
	std::vector<std::string> service_time_deviation_vec{"0.10", "0.25"};
	std::vector<std::string> uncertainty_budget_vec{"1", "5"};

	std::vector<std::string> instances_terminations;

	for (auto num_vehicles : num_vehicles_vec)
		for (auto service_time_deviation : service_time_deviation_vec)
			for (auto uncertainty_budget : uncertainty_budget_vec)
				instances_terminations.push_back("_v" + num_vehicles + "_d" + service_time_deviation + "_b" + uncertainty_budget + ".txt");

	// algorithms.push_back("baseline");
	algorithms.push_back("cb_AVICs_CCCs_baseline");
	// algorithms.push_back("csc");
	algorithms.push_back("cb_AVICs_CCCs_csc");

	instance_types.push_back("C");
	instance_types.push_back("R");
	instance_types.push_back("RC");

	instance_sizes.push_back("25");
	instance_sizes.push_back("50");
	instance_sizes.push_back("100");

	instance_limit_quantiles.push_back("0.8");
	instance_limit_quantiles.push_back("1");
	instance_limit_quantiles.push_back("2");
	instance_limit_quantiles.push_back("3");
	instance_limit_quantiles.push_back("4");

	std::fstream output;
	std::string output_name;
	output_name = "..//tables//latex//table_algorithms.txt";
	output.open(output_name.c_str(), std::fstream::out);

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
	size_t num_inst_per_vertex_size_type_quantile = num_vehicles_vec.size() * service_time_deviation_vec.size() * uncertainty_budget_vec.size();
	size_t num_inst_per_vertex_size_quantile = instance_types.size() * num_inst_per_vertex_size_type_quantile;
	size_t num_inst_per_vertex_size = num_inst_per_vertex_size_quantile * instance_limit_quantiles.size();
	size_t total_num_instances = num_inst_per_vertex_size * instance_sizes.size();

	for (auto instance_size : instance_sizes)
	{
		std::vector<std::vector<double>> time_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<double> avg_time_inst_size(algorithms.size(), 0.0);
		std::vector<double> avg_gap_inst_size(algorithms.size(), 0.0);
		std::vector<double> st_dev_inst_size(algorithms.size(), 0.0);

		std::vector<int> num_optimal_inst_size(algorithms.size(), 0);

		for (auto instance_type : instance_types)
		{
			for (auto instance_limit_quantile : instance_limit_quantiles)
			{
				std::vector<std::vector<double>> time_per_algo_quantile(algorithms.size(), std::vector<double>());
				std::vector<std::vector<double>> gap_per_algo_quantile(algorithms.size(), std::vector<double>());
				std::vector<double> avg_time_quantile(algorithms.size(), 0.0);
				std::vector<double> avg_gap_quantile(algorithms.size(), 0.0);
				std::vector<double> st_dev_quantile(algorithms.size(), 0.0);

				std::vector<int> num_optimal_quantile(algorithms.size(), 0);

				std::string instance_prefix = instance_type + instance_size + "_" + instance_limit_quantile;

				for (auto instances_termination : instances_terminations)
				{
					std::string instance = instance_prefix + instances_termination;
					// std::cout << instance << std::endl;
					double original_lp = 0.0;
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

						time = std::min(time, 3600.0);

						if ((status == "OPTIMAL") || (status == "INFEASIBLE"))
							std::cout << instance << std::endl;

						if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
						{
							if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
							{
								if ((!double_greater(ub, lb)))
								{
									++(num_optimal_inst_size[algo]);
									++(num_optimal_quantile[algo]);
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
							++(num_optimal_quantile[algo]);
							++(total_num_optimal[algo]);
							// time_per_algo_inst_size[algo].push_back(time);
							// time_per_algo_quantile[algo].push_back(time);
							// total_time_per_algo[algo].push_back(time);

							// total_avg_time[algo] += time;
							// avg_time_inst_size[algo] += time;
							// avg_time_quantile[algo] += time;
						}

						time_per_algo_inst_size[algo].push_back(time);
						time_per_algo_quantile[algo].push_back(time);
						total_time_per_algo[algo].push_back(time);

						total_avg_time[algo] += time;
						avg_time_inst_size[algo] += time;
						avg_time_quantile[algo] += time;

						// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

						// if (!double_equals(gap, 0.0))
						// {
						gap_per_algo_inst_size[algo].push_back(gap);
						gap_per_algo_quantile[algo].push_back(gap);
						total_gap_per_algo[algo].push_back(gap);

						total_avg_gap[algo] += gap;
						avg_gap_quantile[algo] += gap;
						avg_gap_inst_size[algo] += gap;
						// }
						input.close();
					}

					// getchar();getchar();
				}

				output << instance_type + instance_size + "Q" + instance_limit_quantile << " & & " << instance_size;

				for (size_t algo = 0; algo < algorithms.size(); ++algo)
				{
					if ((time_per_algo_quantile[algo]).size() > 0)
						avg_time_quantile[algo] /= (1.0 * ((time_per_algo_quantile[algo]).size()));
					else
						avg_time_quantile[algo] = -1;
					if ((gap_per_algo_quantile[algo]).size() > 0)
					{
						avg_gap_quantile[algo] /= (1.0 * ((gap_per_algo_quantile[algo]).size()));
						st_dev_quantile[algo] = StDev(gap_per_algo_quantile[algo], avg_gap_quantile[algo]);
					}
					else
						avg_gap_quantile[algo] = st_dev_quantile[algo] = -1;
					output << " & & " << num_optimal_quantile[algo] << "/" << num_inst_per_vertex_size_type_quantile << " & " << avg_time_quantile[algo] << " & " << avg_gap_quantile[algo] << " & " << st_dev_quantile[algo];
				}
				output << "\\\\" << std::endl;
			}
		}
		output << "Sub-total & & ";

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

	output << "Total & & ";
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
		output << "& & " << total_num_optimal[algo] << "/" << total_num_instances << " & " << total_avg_time[algo] << " & " << total_avg_gap[algo] << " & " << StDev(total_gap_per_algo[algo], total_avg_gap[algo]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateHeuristicsLatexTable(std::string folder_exact, std::string folder_heuristic, bool add_exact_results)
{
	std::string curr_file;
	std::vector<std::pair<int, std::string>> algorithms;
	std::vector<std::string> exact_algorithms;
	std::vector<std::string> instance_sizes;
	std::vector<std::string> instance_types;
	std::vector<std::string> instance_limit_quantiles;

	std::vector<std::string> num_vehicles_vec{"2", "4"};
	std::vector<std::string> service_time_deviation_vec{"0.10", "0.25"};
	std::vector<std::string> uncertainty_budget_vec{"1", "5"};

	std::vector<std::string> instances_terminations;

	for (auto num_vehicles : num_vehicles_vec)
		for (auto service_time_deviation : service_time_deviation_vec)
			for (auto uncertainty_budget : uncertainty_budget_vec)
				instances_terminations.push_back("_v" + num_vehicles + "_d" + service_time_deviation + "_b" + uncertainty_budget + ".txt");

	exact_algorithms.push_back("baseline");
	// algorithms.push_back("cb_baseline");
	exact_algorithms.push_back("csc");
	// exact_algorithms.push_back("cb_AVICs_CCCs_csc");

	// algorithms.push_back("baseline_ks_b5_[84,19]_d0.96_feas");
	algorithms.push_back(std::pair<int, std::string>(0, "csc_ks_b5_[72,25]_d0.86"));
	algorithms.push_back(std::pair<int, std::string>(10, "sa_MT_df_0.98"));

	instance_types.push_back("C");
	instance_types.push_back("R");
	instance_types.push_back("RC");

	instance_sizes.push_back("25");
	instance_sizes.push_back("50");
	instance_sizes.push_back("100");

	instance_limit_quantiles.push_back("0.8");
	instance_limit_quantiles.push_back("1");
	instance_limit_quantiles.push_back("2");
	instance_limit_quantiles.push_back("3");
	instance_limit_quantiles.push_back("4");

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
	size_t num_inst_per_vertex_size_type_quantile_exact = num_vehicles_vec.size() * service_time_deviation_vec.size() * uncertainty_budget_vec.size();
	size_t num_inst_per_vertex_size_quantile_exact = instance_types.size() * num_inst_per_vertex_size_type_quantile_exact;
	size_t num_inst_per_vertex_size_exact = num_inst_per_vertex_size_quantile_exact * instance_limit_quantiles.size();
	size_t total_num_instances_exact = num_inst_per_vertex_size_exact * instance_sizes.size();

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(), 0.0);
	std::vector<double> total_avg_gap(algorithms.size(), 0.0);
	std::vector<int> total_num_best_known_bound(algorithms.size(), 0);
	size_t num_inst_per_vertex_size_type_quantile = num_vehicles_vec.size() * service_time_deviation_vec.size() * uncertainty_budget_vec.size();
	size_t num_inst_per_vertex_size_quantile = instance_types.size() * num_inst_per_vertex_size_type_quantile;
	size_t num_inst_per_vertex_size = num_inst_per_vertex_size_quantile * instance_limit_quantiles.size();
	size_t total_num_instances = num_inst_per_vertex_size * instance_sizes.size();

	for (auto instance_size : instance_sizes)
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

		for (auto instance_type : instance_types)
		{
			for (auto instance_limit_quantile : instance_limit_quantiles)
			{
				std::vector<std::vector<double>> time_per_algo_quantile_exact(exact_algorithms.size(), std::vector<double>());
				std::vector<std::vector<double>> gap_per_algo_quantile_exact(exact_algorithms.size(), std::vector<double>());
				std::vector<double> avg_time_quantile_exact(exact_algorithms.size(), 0.0);
				std::vector<double> avg_gap_quantile_exact(exact_algorithms.size(), 0.0);
				std::vector<double> st_dev_quantile_exact(exact_algorithms.size(), 0.0);

				std::vector<int> num_optimal_quantile_exact(exact_algorithms.size(), 0);

				std::vector<std::vector<double>> time_per_algo_quantile(algorithms.size(), std::vector<double>());
				std::vector<std::vector<double>> gap_per_algo_quantile(algorithms.size(), std::vector<double>());
				std::vector<double> avg_time_quantile(algorithms.size(), 0.0);
				std::vector<double> avg_gap_quantile(algorithms.size(), 0.0);
				std::vector<double> st_dev_quantile(algorithms.size(), 0.0);

				std::vector<int> num_best_known_bound_quantile(algorithms.size(), 0);

				std::string instance_prefix = instance_type + instance_size + "_" + instance_limit_quantile;

				for (auto instances_termination : instances_terminations)
				{
					std::string instance = instance_prefix + instances_termination;
					// std::cout << instance << std::endl;
					double best_lb = -1;
					bool is_infeasible = false;
					bool has_optimal_bound = false;

					// compute the best known primal bound among the exact algorithms.
					for (size_t algo = 0; algo < exact_algorithms.size(); ++algo)
					{
						curr_file = "..//solutions//" + folder_exact + "//";
						curr_file.append("s_");
						curr_file.append(exact_algorithms[algo]);
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

						if (status == "INFEASIBLE")
						{
							is_infeasible = true;
							break;
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
						curr_file = "..//solutions//" + folder_exact + "//";
						curr_file.append("s_");
						curr_file.append(exact_algorithms[algo]);
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

						// if ((status == "OPTIMAL") || (status == "INFEASIBLE"))
						// 	std::cout << instance << std::endl;

						if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
						{
							if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
							{
								if ((!double_greater(ub, lb)))
								{
									++(num_optimal_inst_size_exact[algo]);
									++(num_optimal_quantile_exact[algo]);
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
							++(num_optimal_quantile_exact[algo]);
							++(total_num_optimal_exact[algo]);
							// time_per_algo_inst_size_exact[algo].push_back(time);
							// time_per_algo_quantile_exact[algo].push_back(time);
							// total_time_per_algo_exact[algo].push_back(time);

							// total_avg_time_exact[algo] += time;
							// avg_time_inst_size_exact[algo] += time;
							// avg_time_quantile_exact[algo] += time;
						}

						time_per_algo_inst_size_exact[algo].push_back(time);
						time_per_algo_quantile_exact[algo].push_back(time);
						total_time_per_algo_exact[algo].push_back(time);

						total_avg_time_exact[algo] += time;
						avg_time_inst_size_exact[algo] += time;
						avg_time_quantile_exact[algo] += time;

						// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

						// if (!double_equals(gap, 0.0))
						{
							gap_per_algo_inst_size_exact[algo].push_back(gap);
							gap_per_algo_quantile_exact[algo].push_back(gap);
							total_gap_per_algo_exact[algo].push_back(gap);

							total_avg_gap_exact[algo] += gap;
							avg_gap_quantile_exact[algo] += gap;
							avg_gap_inst_size_exact[algo] += gap;
						}
						input.close();
					}

					for (size_t algo = 0; algo < algorithms.size(); ++algo)
					{

						int num_seeds = algorithms[algo].first;

						if (num_seeds == 0)
						{
							curr_file = "..//solutions//" + folder_heuristic + "//";
							curr_file.append("s_");
							curr_file.append(algorithms[algo].second);
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

							// std::cout << status << std::endl;
							getline(input, line);
							getline(input, line);

							pos = line.find_first_of(":");
							s_t1 << line.substr(pos + 2);
							s_t1 >> heuristic_time;

							if (status != "INFEASIBLE")
							{
								getline(input, line);
								getline(input, line);
								pos = line.find_first_of(":");
								s_lb1 << line.substr(pos + 2);
								s_lb1 >> heuristic_lb;
								heuristic_lb = round_decimals(heuristic_lb, 2); // IMPORTANT because I saved heuristic solutions file without 2 decimal precision.
							}

							// if (double_greater(heuristic_time, 1500))
							// std::cout << heuristic_time << std::endl;

							time_per_algo_inst_size[algo].push_back(heuristic_time);
							time_per_algo_quantile[algo].push_back(heuristic_time);
							total_time_per_algo[algo].push_back(heuristic_time);

							total_avg_time[algo] += heuristic_time;
							avg_time_inst_size[algo] += heuristic_time;
							avg_time_quantile[algo] += heuristic_time;

							// if ((status == "INFEASIBLE") && !is_infeasible)
							// 	std::cout << "Heuristica provou inviabilidade nova" << std::endl;

							if (((status == "INFEASIBLE") && is_infeasible) || double_equals(best_lb, heuristic_lb))
							{
								++(num_best_known_bound_inst_size[algo]);
								++(num_best_known_bound_quantile[algo]);
								++(total_num_best_known_bound[algo]);
								improvement = 0;
							}
							else
							{
								if ((status == "INFEASIBLE") && !is_infeasible)
								{
									++(num_best_known_bound_inst_size[algo]);
									++(num_best_known_bound_quantile[algo]);
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
									++(num_best_known_bound_quantile[algo]);
									++(total_num_best_known_bound[algo]);
								}
								else
								{
									improvement = (100 * (heuristic_lb - best_lb)) / best_lb;
									std::cout << best_lb << " " << heuristic_lb << " " << status << std::endl;
									if (!double_less(improvement, 0.0))
									{
										++(num_best_known_bound_inst_size[algo]);
										++(num_best_known_bound_quantile[algo]);
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
							gap_per_algo_quantile[algo].push_back(improvement);
							total_gap_per_algo[algo].push_back(improvement);

							avg_gap_quantile[algo] += improvement;
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
							time_per_algo_quantile[algo].push_back(avg_all_seeds_heuristic_time);
							total_time_per_algo[algo].push_back(avg_all_seeds_heuristic_time);

							total_avg_time[algo] += avg_all_seeds_heuristic_time;
							avg_time_inst_size[algo] += avg_all_seeds_heuristic_time;
							avg_time_quantile[algo] += avg_all_seeds_heuristic_time;

							avg_all_seeds_improvement /= num_seeds;
							gap_per_algo_inst_size[algo].push_back(avg_all_seeds_improvement);
							gap_per_algo_quantile[algo].push_back(avg_all_seeds_improvement);
							total_gap_per_algo[algo].push_back(avg_all_seeds_improvement);

							avg_gap_quantile[algo] += avg_all_seeds_improvement;
							avg_gap_inst_size[algo] += avg_all_seeds_improvement;
							total_avg_gap[algo] += avg_all_seeds_improvement;

							avg_all_num_best_known_bound_inst_size /= num_seeds;
							avg_all_num_best_known_bound_quantile /= num_seeds;
							avg_all_total_num_best_known_bound /= num_seeds;

							num_best_known_bound_inst_size[algo] += 1.0;
							num_best_known_bound_quantile[algo] += 1.0;
							total_num_best_known_bound[algo] += 1.0;
						}

						// std::cout << algorithms[algo] << " " << heuristic_lb << " " << heuristic_time << " " << improvement << std::endl;
					}

					// getchar();getchar();
				}

				output << instance_type + instance_size + "Q" + instance_limit_quantile << " & & " << instance_size;

				if (add_exact_results)
				{
					for (size_t algo = 0; algo < exact_algorithms.size(); ++algo)
					{
						if ((time_per_algo_quantile_exact[algo]).size() > 0)
							avg_time_quantile_exact[algo] /= (1.0 * ((time_per_algo_quantile_exact[algo]).size()));
						else
							avg_time_quantile_exact[algo] = -1;
						if ((gap_per_algo_quantile_exact[algo]).size() > 0)
						{
							avg_gap_quantile_exact[algo] /= (1.0 * ((gap_per_algo_quantile_exact[algo]).size()));
							st_dev_quantile_exact[algo] = StDev(gap_per_algo_quantile_exact[algo], avg_gap_quantile_exact[algo]);
						}
						else
							avg_gap_quantile_exact[algo] = st_dev_quantile_exact[algo] = -1;
						output << " & & " << num_optimal_quantile_exact[algo] << "/" << num_inst_per_vertex_size_type_quantile_exact << " & " << avg_time_quantile_exact[algo] << " & " << avg_gap_quantile_exact[algo] << " & " << st_dev_quantile_exact[algo];
					}
				}

				for (size_t algo = 0; algo < algorithms.size(); ++algo)
				{
					if ((time_per_algo_quantile[algo]).size() > 0)
						avg_time_quantile[algo] /= (1.0 * ((time_per_algo_quantile[algo]).size()));
					else
						avg_time_quantile[algo] = -1;
					if ((gap_per_algo_quantile[algo]).size() > 0)
					{
						avg_gap_quantile[algo] /= (1.0 * ((gap_per_algo_quantile[algo]).size()));
						st_dev_quantile[algo] = StDev(gap_per_algo_quantile[algo], avg_gap_quantile[algo]);
					}
					else
						avg_gap_quantile[algo] = st_dev_quantile[algo] = -1;
					output << " & & " << num_best_known_bound_quantile[algo] << "/" << num_inst_per_vertex_size_type_quantile << " & " << avg_time_quantile[algo] << " & " << avg_gap_quantile[algo] << " & " << st_dev_quantile[algo];
				}
				output << "\\\\" << std::endl;
			}
		}
		output << "Sub-total & &";

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
	std::vector<std::string> instance_sizes;
	std::vector<std::string> instance_types;
	std::vector<std::string> instance_limit_quantiles;

	std::vector<std::string> num_vehicles_vec{"2", "4"};
	std::vector<std::string> service_time_deviation_vec{"0.10", "0.25"};
	std::vector<std::string> uncertainty_budget_vec{"1", "5"};

	std::vector<std::string> instances_terminations;

	for (auto num_vehicles : num_vehicles_vec)
		for (auto service_time_deviation : service_time_deviation_vec)
			for (auto uncertainty_budget : uncertainty_budget_vec)
				instances_terminations.push_back("_v" + num_vehicles + "_d" + service_time_deviation + "_b" + uncertainty_budget + ".txt");

	algorithms.push_back("relax_baseline");
	// algorithms.push_back("relax_cb_AVICs_baseline");
	// algorithms.push_back("relax_cb_CCCs_baseline");
	algorithms.push_back("relax_cb_AVICs_CCCs_baseline");
	algorithms.push_back("relax_csc");
	// algorithms.push_back("relax_cb_AVICs_csc");
	// algorithms.push_back("relax_cb_CCCs_csc");
	algorithms.push_back("relax_cb_AVICs_CCCs_csc");

	instance_types.push_back("C");
	instance_types.push_back("R");
	instance_types.push_back("RC");

	instance_sizes.push_back("25");
	instance_sizes.push_back("50");
	instance_sizes.push_back("100");

	instance_limit_quantiles.push_back("0.8");

	instance_limit_quantiles.push_back("1");
	instance_limit_quantiles.push_back("2");
	instance_limit_quantiles.push_back("3");
	instance_limit_quantiles.push_back("4");

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

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(), 0.0);
	int total_num_instances = 0;

	for (auto instance_size : instance_sizes)
	{
		std::vector<std::vector<double>> improvement_per_algo_inst_size(algorithms.size(), std::vector<double>());
		std::vector<double> avg_improvement_per_inst_size(algorithms.size(), 0.0);
		std::vector<double> st_dev_per_inst_size(algorithms.size(), 0.0);
		for (auto instance_type : instance_types)
		{
			for (auto instance_limit_quantile : instance_limit_quantiles)
			{
				std::vector<std::vector<double>> improvement_per_algo(algorithms.size(), std::vector<double>());
				std::vector<double> avg_improvement(algorithms.size(), 0.0);
				std::vector<double> st_dev(algorithms.size(), 0.0);

				std::string instance_prefix = instance_type + instance_size + "_" + instance_limit_quantile;

				for (auto instances_termination : instances_terminations)
				{
					std::string instance = instance_prefix + instances_termination;
					std::cout << instance << std::endl;
					double original_lp = 0.0;
					for (size_t algo = 0; algo < algorithms.size(); ++algo)
					{
						double curr_improvement = 0.0;
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
								if (!double_greater(original_lp, lp))
									std::cout << original_lp << " - " << lp << std::endl;
							}
						}

						// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

						if (!double_equals(curr_improvement, -1))
						{
							improvement_per_algo[algo].push_back(curr_improvement);
							improvement_per_algo_inst_size[algo].push_back(curr_improvement);
							total_improvement_per_algo[algo].push_back(curr_improvement);
							total_avg_improvement[algo] += curr_improvement;
							avg_improvement[algo] += curr_improvement;
							avg_improvement_per_inst_size[algo] += curr_improvement;
						}
						input.close();
					}

					// getchar();getchar();
				}

				output << instance_type + instance_size + "Q" + instance_limit_quantile << " & & " << instance_size << " & " << instances_terminations.size();

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
		}
		int subtotal_num_instances = instance_types.size() * instance_limit_quantiles.size() * instances_terminations.size();
		total_num_instances += subtotal_num_instances;

		output << "Subtotal & & & " << subtotal_num_instances;

		for (size_t algo = 1; algo < algorithms.size(); ++algo)
		{
			if ((improvement_per_algo_inst_size[algo]).size() > 0)
				avg_improvement_per_inst_size[algo] /= (1.0 * ((improvement_per_algo_inst_size[algo]).size()));
			else
				avg_improvement_per_inst_size[algo] = -1;
			st_dev_per_inst_size[algo] = StDev(improvement_per_algo_inst_size[algo], avg_improvement_per_inst_size[algo]);
			output << " & & " << avg_improvement_per_inst_size[algo] << " & " << st_dev_per_inst_size[algo];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total & & & " << total_num_instances;
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
	std::string folder_heuristic_sol = "all_heuristic";
	std::string folder_exact_sol = "all_exact";
	std::string folder_relax = "2024-07-26_09:18:48_all_relax_new";
	// try
	// {
	// GenerateAlgorithmsLatexTablePerInstance(folder);
	// return 1;
	// GenerateLPImprovementsLatexTable(folder_relax);
	GenerateAlgorithmsLatexTable(folder_exact_sol);
	// GenerateHeuristicsLatexTable(folder_exact_sol, folder_heuristic_sol, false);
	return 0;
	int time_limit = -1;
	int num_routes = 5;
	int uncertainty_budget = 10;
	auto dev = 0.5;
	int seed = 10;
	std::string instance_name = "R100_1.txt";
	Instance inst("/home/lucas/Documentos/Research/msca-humanitarian-optimization/instances/R-STOP-DP/", instance_name, num_routes, dev, uncertainty_budget, false);

	// std::cout << inst << std::endl;
	// std::list<int> x{1, 2};
	// for (std::list<int>::iterator it = x.begin(); it != x.end(); ++it)
	// 	std::cout << (*it) << " ";
	// std::cout << std::endl;
	// auto [new_route_sum_profits1, new_route_max_duration1] = inst.ComputeRouteCosts(x);
	// std::cout << "new profit1 " << new_route_sum_profits1 << std::endl;
	// std::cout << "new duration1 " << new_route_max_duration1 << std::endl;

	// auto [new_route_sum_profits2, new_route_max_duration2] = inst.ComputeRouteCostsRec(x, true);
	// std::cout << "new profit2 " << new_route_sum_profits2 << std::endl;
	// std::cout << "new duration2 " << new_route_max_duration2 << std::endl;

	// auto [new_route_sum_profits3, new_route_max_duration3] = inst.ComputeRouteCostsRec(x, false);
	// std::cout << "new profit3 " << new_route_sum_profits3 << std::endl;
	// std::cout << "new duration3 " << new_route_max_duration3 << std::endl;

	// return 0;

	auto graph = inst.graph();

	// MetaHeuristicSolution sol(graph->num_vertices(), graph->num_arcs(), num_routes);

	// std::cout << inst << std::endl;
	// double time_variation = 0.0;
	// double profit_variation = 0.0;
	// std::list<int>::iterator it;
	// std::cout << "Add 1 to route 0" << std::endl;
	// if (sol.PreviewAddVertex(inst, 1, 0, 0, it, profit_variation, time_variation))
	// 	sol.AddVertex(1, 0, it, profit_variation, time_variation);

	// std::cout << "Add 2 to route 0" << std::endl;
	// if (sol.PreviewAddVertex(inst, 2, 0, 0, it, profit_variation, time_variation))
	// 	sol.AddVertex(2, 0, it, profit_variation, time_variation);

	// std::cout << "Add 4 to route 0" << std::endl;
	// if (sol.PreviewAddVertex(inst, 4, 0, 2, it, profit_variation, time_variation))
	// 	sol.AddVertex(4, 0, it, profit_variation, time_variation);

	// std::cout << "Add 6 to route 0" << std::endl;
	// if (sol.PreviewAddVertex(inst, 6, 0, 1, it, profit_variation, time_variation))
	// 	sol.AddVertex(6, 0, it, profit_variation, time_variation);

	// std::cout << "Add 3 to route 1" << std::endl;
	// if (sol.PreviewAddVertex(inst, 3, 1, 0, it, profit_variation, time_variation))
	// 	sol.AddVertex(3, 1, it, profit_variation, time_variation);

	// std::cout << sol << std::endl;

	// // std::cout << "Move 3 to route 0, pos 1" << std::endl;
	// // double time_variation2 = 0;
	// // int profit_variation2 = 0;
	// // if (sol.PreviewInterRouteMoveVertex(inst, 3, 0, 1, it, profit_variation, time_variation, profit_variation2, time_variation2))
	// // 	sol.InterRouteMoveVertex(3, 0, it, profit_variation, time_variation, profit_variation2, time_variation2);

	// // std::cout << sol << std::endl;

	// // std::cout << "Move 3 to route 1, pos 1" << std::endl;
	// // if (sol.PreviewInterRouteMoveVertex(inst, 3, 1, 1, it, profit_variation, time_variation, profit_variation2, time_variation2))
	// // 	sol.InterRouteMoveVertex(3, 1, it, profit_variation, time_variation, profit_variation2, time_variation2);

	// // std::cout << sol << std::endl;

	// // std::cout << "Move 2 to route 1, pos 0" << std::endl;
	// // if (sol.PreviewInterRouteMoveVertex(inst, 2, 1, 0, it, profit_variation, time_variation, profit_variation2, time_variation2))
	// // 	sol.InterRouteMoveVertex(2, 1, it, profit_variation, time_variation, profit_variation2, time_variation2);

	// // std::cout << sol << std::endl;

	// // std::cout << "Move 1 to route 1, pos 0" << std::endl;
	// // if (sol.PreviewInterRouteMoveVertex(inst, 1, 1, 0, it, profit_variation, time_variation, profit_variation2, time_variation2))
	// // 	sol.InterRouteMoveVertex(1, 1, it, profit_variation, time_variation, profit_variation2, time_variation2);

	// // std::cout << time_variation2 << std::endl;

	// std::list<int>::iterator it_i1, it_i2, it_f1, it_f2;

	// double profit_variation1 = 0, profit_variation2 = 0;
	// double time_variation1 = 0.0, time_variation2 = 0.0;
	// int r1 = 0;
	// int r2 = 1;
	// int i1 = 0;
	// int f1 = 3;
	// int i2 = 0;
	// int f2 = 0;

	// if (sol.PreviewInterRouteSwap(inst, r1, i1, f1, r2, i2, f2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2))
	// 	sol.InterRouteSwap(r1, r2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2);

	// std::cout << sol << std::endl;

	// return 0;

	// r1 = 1;
	// i1 = 2;
	// f1 = 2;
	// int unrouted = 10;
	// it = ((sol.vertex_status_vec_)[unrouted]).pos_;
	// if (sol.PreviewInterRouteSwapUnrouted(inst, r1, i1, f1, it_i1, it_f1, profit_variation1, time_variation1, it))
	// 	sol.InterRouteSwapUnrouted(r1, it_i1, it_f1, profit_variation1, time_variation1, it);

	// std::cout << sol << std::endl;

	// auto vertex = 1;
	// auto curr_route = 1;
	// if (sol.PreviewAddVertexToRouteWithinMinimumDistanceIncrease(inst, vertex, curr_route, it, profit_variation, time_variation))
	// {
	// 	std::cout << "Add vertex " << vertex << " to route " << curr_route << " right before vertex " << *it << std::endl;
	// 	sol.AddVertex(vertex, curr_route, it, profit_variation, time_variation);
	// }
	// std::cout << "Remove 2" << std::endl;
	// if (sol.PreviewRemoveVertex(inst, 2, profit_variation, time_variation))
	// 	sol.RemoveVertex(2, profit_variation, time_variation);

	// std::cout << "Remove 1" << std::endl;
	// if (sol.PreviewRemoveVertex(inst, 1, profit_variation, time_variation))
	// 	sol.RemoveVertex(1, profit_variation, time_variation);

	// std::cout << "Remove 3" << std::endl;
	// if (sol.PreviewRemoveVertex(inst, 3, profit_variation, time_variation))
	// 	sol.RemoveVertex(3, profit_variation, time_variation);

	// std::cout << sol << std::endl;

	// vertex = 5;
	// if (sol.PreviewAddVertexWithinMinimumDistanceIncrease(inst, vertex, curr_route, it, profit_variation, time_variation))
	// {
	// 	std::cout << "Add vertex " << vertex << " to route " << curr_route << " right before vertex " << *it << std::endl;
	// 	sol.AddVertex(vertex, curr_route, it, profit_variation, time_variation);
	// }

	// std::cout << sol << std::endl;

	// std::cout << inst << std::endl;

	Timestamp *ti = NewTimestamp();
	Timer *timer = GetTimer();
	timer->Clock(ti);
	FeasibilityPump fp;
	fp.Init(inst);
	std::cout << (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed) << std::endl;
	fp.Run();
	(fp.solution_).WriteToFile(inst, (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), "//", instance_name);

	std::cout << fp.solution_.profits_sum_ << std::endl;
	std::cout << "time: " << timer->CurrentElapsedTime(ti) << std::endl;

	timer->Clock(ti);
	try
	{
		auto pool_size = K_ALNS_SIZE_OF_POOL;
		// find initial solution
		timer->Clock(ti);
		auto initial_solution = InitalSolutionGenerator::GenerateInitialSolution(inst, K_ALNS_MULTI_THREAD);
		std::cout << "found initial solution in " << timer->CurrentElapsedTime(ti) << " s" << std::endl;
		// std::cout << *initial_solution << std::endl;
		ALNS alns(inst, initial_solution, pool_size);
		SimulatedAnnealing sa(inst, initial_solution);

		// seed = time(nullptr);
		srand(time(nullptr));

		delete initial_solution;
		initial_solution = nullptr;

		// FeasibilityPump fp;
		// fp.Init(inst);
		// std::cout << (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed) << std::endl;
		// fp.Run();
		// (fp.solution_).WriteToFile(inst, (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), "", instance_name);
		// std::cout << (fp.solution_) << std::endl;
		// ALNS alns(inst, (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), "", instance_name, pool_size);
		// std::cout << "found initial solution in " << fp.solution_.time_stage1_ + fp.solution_.time_stage2_ << " s" << std::endl;
		// std::cout << "found initial solution in " << (alns.best_solution())->total_time_spent_ << " s" << std::endl;

		timer->Clock(ti);
		alns.Run(K_ALNS_NUM_ITERATIONS, K_ALNS_MULTI_THREAD);
		std::cout << MetaHeuristicSolution::GenerateALNSFileName(K_ALNS_NUM_ITERATIONS, K_ALNS_SIZE_OF_POOL, K_ALNS_MULTI_THREAD) << "_seed_" << seed << std::endl;
		std::cout << "best solution: " << alns.best_solution()->profits_sum_ << std::endl;
		std::cout << "elapsed time: " << timer->CurrentElapsedTime(ti) << " s" << std::endl;
		std::cout << "elapsed time: " << (alns.best_solution())->total_time_spent_ << " s" << std::endl;
		std::cout << "num iter: " << alns.best_solution()->num_iterations_ << std::endl;
		std::cout << "las improve iter: " << alns.best_solution()->last_improve_iteration_ << std::endl;

		// std::cout << *(alns.best_solution()) << std::endl;
		(alns.best_solution())->WriteToFile(inst, MetaHeuristicSolution::GenerateALNSFileName(K_ALNS_NUM_ITERATIONS, K_ALNS_SIZE_OF_POOL, K_ALNS_MULTI_THREAD) + "_seed_" + std::to_string(seed), "", instance_name);

		timer->Clock(ti);
		sa.Run(K_SA_TEMP_DECAY_RATE, K_SA_SAMPLING_SIZE, K_SA_TARGET_ACCEPTANCE_PROBABILITY, K_ALNS_MULTI_THREAD);
		std::cout << MetaHeuristicSolution::GenerateSimulatedAnnealingFileName(K_SA_TEMP_DECAY_RATE, K_ALNS_MULTI_THREAD) << "_seed_" << seed << std::endl;
		std::cout << "best solution: " << sa.best_solution()->profits_sum_ << std::endl;
		std::cout << "elapsed time: " << timer->CurrentElapsedTime(ti) << " s" << std::endl;
		std::cout << "elapsed time: " << (sa.best_solution())->total_time_spent_ << " s" << std::endl;
		std::cout << "num iter: " << sa.best_solution()->num_iterations_ << std::endl;
		std::cout << "las improve iter: " << sa.best_solution()->last_improve_iteration_ << std::endl;

		// std::cout << *(alns.best_solution()) << std::endl;
		(alns.best_solution())->WriteToFile(inst, MetaHeuristicSolution::GenerateSimulatedAnnealingFileName(K_SA_TEMP_DECAY_RATE, K_ALNS_MULTI_THREAD) + "_seed_" + std::to_string(seed), "", instance_name);

		delete ti;
		ti = nullptr;
		DeleteTimer();
	}
	catch (const std::runtime_error &re)
	{
		std::cout << "Runtime error: " << re.what() << std::endl;
	}
	catch (const std::exception &ex)
	{
		std::cout << "Error occurred: " << ex.what() << std::endl;
	}
	catch (const int &error)
	{
		std::cout << "Error occurred: " << error << std::endl;
	}
	catch (IloException &e)
	{
		std::cout << "Concert Exception: " << e << std::endl;
	}
	catch (const char *e)
	{
		std::cout << e << std::endl;
	}
	catch (const std::string &e)
	{
		std::cout << e << std::endl;
	}
	catch (...)
	{
		std::cout << "Unknown failure occurred. Possible memory corruption" << std::endl;
	}

	return 0;
	// FeasibilityPump fp;
	// fp.Init(inst);
	// std::cout << (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed) << std::endl;
	// fp.Run();
	// (fp.solution_).WriteToFile(inst, (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), "//", instance_name);

	// std::cout << fp.solution_.profits_sum_ << std::endl;

	// KernelSearch ks(inst);
	// std::cout << KSHeuristicSolution::GenerateFileName(Formulation::single_commodity, K_KS_MAX_SIZE_BUCKET, K_KS_MIN_TIME_LIMIT, K_KS_MAX_TIME_LIMIT, K_KS_DECAY_FACTOR_TIME_LIMIT, true) << std::endl;
	// auto ks_solution = ks.Run(Formulation::single_commodity, K_KS_MAX_SIZE_BUCKET, K_KS_MIN_TIME_LIMIT, K_KS_MAX_TIME_LIMIT, K_KS_DECAY_FACTOR_TIME_LIMIT, true);

	// ks_solution->WriteToFile(inst, KSHeuristicSolution::GenerateFileName(Formulation::single_commodity, K_KS_MAX_SIZE_BUCKET, K_KS_MIN_TIME_LIMIT, K_KS_MAX_TIME_LIMIT, K_KS_DECAY_FACTOR_TIME_LIMIT, true), "//", instance_name);

	// std::cout << ks_solution->profits_sum_ << std::endl;

	// delete ks_solution;
	// ks_solution = nullptr;

	// return 0;
	// Route route;

	// route.vertices_.push_back(1);
	// route.vertices_.push_back(2);
	// route.vertices_.push_back(10);
	// route.vertices_.push_back(11);

	// // std::cout << inst << std::endl;

	// double route_profits, route_time;
	// tie(route_profits, route_time) = inst.ComputeRouteCosts(route);
	// std::cout << route_profits << " " << route_time << std::endl;

	// tie(route_profits, route_time) = inst.ComputeRouteCostsRec(route, false);
	// std::cout << route_profits << " " << route_time << std::endl;

	// tie(route_profits, route_time) = inst.ComputeRouteCostsRec(route, true);
	// std::cout << route_profits << " " << route_time << std::endl;
	// return 0;

	std::vector<bool> *CALLBACKS_SELECTION = GetCallbackSelection();
	(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = false;
	(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = false;

	// std::cout << inst << std::endl;

	// std::cout << *graph << std::endl;

	double *R0 = Dijkstra(graph, false, false);
	double *Rn = Dijkstra(graph, true, false);

	// for (int i =0; i < graph->num_vertices(); ++i)
	// {
	// 	std::cout << "R0," << i << " = " << R0[i] << std::endl;
	// }

	// inst.WriteToFile("../","teste.txt");
	bool force_use_all_vehicles = false;
	bool export_model = false;
	bool solve_relaxed = false;
	bool use_valid_inequalities = false;
	bool combine_feas_op_cuts = false;
	bool apply_benders_generic_callback = false;
	auto root_cuts = new std::list<UserCutGeneral *>();
	Solution<double> solution(graph->num_vertices());
	Formulation formulation = Formulation::single_commodity;

	if (combine_feas_op_cuts)
		std::cout << "combine cuts" << std::endl;
	else
		std::cout << "classic cuts" << std::endl;
	if (formulation == Formulation::single_commodity)
		std::cout << "single commodity " << std::endl;
	else
		std::cout << "baseline " << std::endl;
	if (apply_benders_generic_callback)
		std::cout << "GENERIC CALLBACK" << std::endl;
	else
		std::cout << "Lazy constraint" << std::endl;

	// // CompactBaseline(inst,R0,Rn,-1,true,false,false,nullptr,nullptr,force_use_all_vehicles,export_model, nullptr,solution);
	// // std::cout <<  " LP: " << solution.lp_ << std::endl;
	// // std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// // solution.reset();
	// // CompactBaseline(inst,R0,Rn,-1,true,use_valid_inequalities,true,nullptr,nullptr,force_use_all_vehicles,export_model, root_cuts,solution);
	// // std::cout <<  " LP: " << solution.lp_ << std::endl;
	// // std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// // solution.reset();
	// // CompactBaseline(inst,R0,Rn,-1,false,use_valid_inequalities,false,root_cuts,nullptr,force_use_all_vehicles,export_model, nullptr,solution);

	// // if(!solution.is_feasible_)
	// // 	std::cout << "Infeasible" << std::endl;
	// // solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// // : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// // std::cout << "num cuts: " << solution.num_cuts_found_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
	// // solution.reset();
	// // DeleteCuts(root_cuts);

	// if(formulation == BendersFormulation::baseline)
	{
		CompactBaseline(inst, R0, Rn, time_limit, solve_relaxed, use_valid_inequalities, false, nullptr, nullptr, force_use_all_vehicles, export_model, nullptr, solution);

		if (solve_relaxed)
		{
			std::cout << " LP: " << solution.lp_ << std::endl;
		}
		else
		{
			solution.is_optimal_ ? std::cout << " optimal: " << solution.lb_ << std::endl
								 : std::cout << " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
		}
		std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

		solution.reset();
	}
	// // // CompactBaseline(inst,R0,Rn,-1,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model, root_cuts,solution);

	// // if (solve_relaxed)
	// // {
	// // 	std::cout <<  " LP: " << solution.lp_ << std::endl;
	// // }
	// // else
	// // {
	// // 	solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// // 	: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// // }
	// // std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// // solution.reset();
	// // auto conflicts = inst.conflicts_list();
	// // std::cout << "conflicts: " << std::endl;
	// // for (auto &conflict: conflicts)
	// // {
	// // 	for( auto& vertex: conflict)
	// // 		std::cout << vertex << " ";

	// // 	std::cout << std::endl;
	// // }

	// if(formulation == BendersFormulation::single_commodity)
	{
		CompactSingleCommodity(inst, R0, Rn, time_limit, solve_relaxed, use_valid_inequalities, false, nullptr, nullptr, force_use_all_vehicles, export_model, nullptr, solution);

		if (solve_relaxed)
		{
			std::cout << " LP: " << solution.lp_ << std::endl;
		}
		else
		{
			solution.is_optimal_ ? std::cout << " optimal: " << solution.lb_ << std::endl
								 : std::cout << " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
		}
		std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

		solution.reset();
	}

	// //getchar();getchar();
	// DeleteCuts(root_cuts);

	// Benders(inst,BendersFormulation::baseline,R0,Rn,time_limit,false,false,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model,nullptr,solution);
	// if (solve_relaxed)
	// {
	// 	std::cout <<  " LP: " << solution.lp_ << std::endl;
	// }
	// else
	// {
	// 	solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// 	: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// }
	// std::cout << "# opt cuts: " << solution.num_benders_opt_cuts_ << std::endl;
	// std::cout << "# feas cuts: " << solution.num_benders_feas_cuts_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();

	// Benders(inst,BendersFormulation::baseline,R0,Rn,time_limit,false,true,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model,nullptr,solution);
	// if (solve_relaxed)
	// {
	// 	std::cout <<  " LP: " << solution.lp_ << std::endl;
	// }
	// else
	// {
	// 	solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// 	: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// }
	// std::cout << "# opt cuts: " << solution.num_benders_opt_cuts_ << std::endl;
	// std::cout << "# feas cuts: " << solution.num_benders_feas_cuts_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
	// solution.reset();

	// Benders(inst, Formulation::baseline, R0, Rn, time_limit, false, false, false, solve_relaxed, use_valid_inequalities, false, nullptr, nullptr, force_use_all_vehicles, export_model, nullptr, solution);
	// if (solve_relaxed)
	// {
	// 	std::cout << " LP: " << solution.lp_ << std::endl;
	// }
	// else
	// {
	// 	solution.is_optimal_ ? std::cout << " optimal: " << solution.lb_ << std::endl
	// 						 : std::cout << " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// }
	// std::cout << "# opt cuts: " << solution.num_benders_opt_cuts_ << std::endl;
	// std::cout << "# feas cuts: " << solution.num_benders_feas_cuts_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();

	// Benders(inst, Formulation::single_commodity, R0, Rn, time_limit, true, true, false, solve_relaxed, use_valid_inequalities, false, nullptr, nullptr, force_use_all_vehicles, export_model, nullptr, solution);
	// if (solve_relaxed)
	// {
	// 	std::cout << " LP: " << solution.lp_ << std::endl;
	// }
	// else
	// {
	// 	solution.is_optimal_ ? std::cout << " optimal: " << solution.lb_ << std::endl
	// 						 : std::cout << " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// }
	// std::cout << "# opt cuts: " << solution.num_benders_opt_cuts_ << std::endl;
	// std::cout << "# feas cuts: " << solution.num_benders_feas_cuts_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();

	// IloEnv env;
	// IloNumArray y_values(env,inst.graph()->num_vertices());
	// // for(int i = 0; i < inst.graph()->num_vertices(); ++i)
	// // {
	// //   y_values[i] = 0;
	// // }
	// y_values[0] = 1;
	// y_values[1] = 0.5;
	// y_values[2] = 0.5;

	// IloNumArray x_values(env,inst.graph()->num_arcs());

	// for(int i  = 0; i < inst.graph()->num_vertices(); ++i)
	// {
	//   for (const auto & j: inst.graph()->AdjVerticesOut(i))
	//   {
	//     x_values[inst.graph()->pos(i,j)] = 0;
	//   }
	// }

	// x_values[inst.graph()->pos(0,1)] = 0.6;
	// x_values[inst.graph()->pos(0,2)] = 0.4;
	// x_values[inst.graph()->pos(2,0)] = 0.4;
	// x_values[inst.graph()->pos(1,0)] = 0.6;

	// // PrimalSubproblemCompactBaseline(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// // std::cout << solution.lp_ << std::endl;
	// // solution.reset();
	// // DualSubproblemCompactBaseline(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// // std::cout << solution.lp_ << std::endl;
	// // solution.reset();

	// PrimalSubproblemCompactSingleCommodity(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// std::cout << solution.lp_ << std::endl;
	// solution.reset();
	// DualSubproblemCompactSingleCommodity(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// std::cout << solution.lp_ << std::endl;

	// env.end();

	delete[] R0;
	R0 = nullptr;
	delete[] Rn;
	Rn = nullptr;

	DeleteTimer();
	DeleteCallbackSelection();
	// }
	// catch (IloException &e)
	// {
	// 	cerr << "Concert Exception: " << e << endl;
	// }
	// catch (...)
	// {
	// 	std::cout << "Other Exception" << endl;
	// }

	return 0;
}
