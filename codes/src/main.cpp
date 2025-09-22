
#include <ctime>
#include <iostream>
#include <vector>
#include <getopt.h>
#include "src/instance.h"
#include "src/exact/formulations.h"
#include "src/kernel_search/kernel_search.h"

void generateInstancesFromBlockRelocationInstancesIter(std::string input_dir_path, std::string output_dir_path, std::string file_name, double percentage_of_items_serviced_by_fleet, int number_of_item_groups)
{
	srand(time(0));
	std::string curr_file = input_dir_path, new_instance_file = output_dir_path + file_name;
	curr_file += (file_name + ".dat");

	// Create an ostringstream object
	std::ostringstream oss;
	// Set the precision and fixed-point notation
	oss << std::fixed << std::setprecision(2) << percentage_of_items_serviced_by_fleet;
	new_instance_file += ("-c" + oss.str() + "-p" + std::to_string(number_of_item_groups) + ".txt");

	std::cout << new_instance_file << std::endl;

	// std::cout << dir_path << " " << file_name << std::endl;

	std::fstream file;
	file.open(curr_file.c_str(), std::fstream::in);
	if (!file.is_open())
	{
		std::cout << "Could not open file" << std::endl;
		throw 1;
		return;
	}

	int num_stacks = 0, num_items = 0, priority = -1, num_items_in_stack = 0;
	int min_priority = std::numeric_limits<int>::max(), max_priority = -1;

	file >> num_stacks >> num_items;
	number_of_item_groups = std::min(num_items, number_of_item_groups);

	std::cout << "num_stacks: " << num_stacks << std::endl;
	std::cout << "num items: " << num_items << std::endl;

	std::vector<int> items_priorities(num_items, -1);
	std::vector<int> items_groups(num_items, -1);
	std::vector<std::vector<int>> stacks(num_stacks);
	// Matrix<int> precedence_matrix(num_items, num_items, 0);

	// std::cout << num_vertices << " " << num_mandatory << std::endl;

	int item_counter = 0;
	double avg_size_stack = 0;
	for (int i = 0; i < num_stacks; ++i)
	{
		file >> num_items_in_stack;
		avg_size_stack += num_items_in_stack;
		for (int j = 0; j < num_items_in_stack; ++j)
		{
			file >> priority;
			items_priorities[item_counter] = priority;
			stacks[i].push_back(item_counter);

			if (min_priority > priority)
				min_priority = priority;

			if (max_priority < priority)
				max_priority = priority;

			// if (j < num_items_in_stack - 1)
			// 	precedence_matrix[item_counter][item_counter + 1] = 1;
			++item_counter;
		}
	}

	file.close();
	avg_size_stack /= num_stacks;

	std::cout << "avg_size_stack: " << avg_size_stack << std::endl;

	// int num_vehicles = std::max(1.0, ceil((num_items * percentage_of_items_serviced_by_fleet) / (0.5 * avg_size_stack)));

	int total_fleet_capacity = ceil(num_items * percentage_of_items_serviced_by_fleet), capacity_serviced = 0;
	// std::cout << "num vehicles: " << num_vehicles << std::endl;

	std::vector<int> vehicles_capacities;
	while (capacity_serviced < total_fleet_capacity)
	{
		auto vehicle_capacity = (0.5 + (rand() % 6) / 10.0) * avg_size_stack;
		vehicles_capacities.push_back(vehicle_capacity);
		capacity_serviced += vehicle_capacity;
	}

	// define the ranges of priotity values that define to which group each item belongs.
	double range_portion = (max_priority - min_priority + 1) / (1.0 * number_of_item_groups);

	std::cout << "min - max: " << min_priority << " - " << max_priority << std::endl;
	std::cout << "range portion: " << range_portion << std::endl;

	for (int k = 0; k < num_items; ++k)
		items_groups[k] = floor((items_priorities[k] - min_priority) / range_portion);

	// int num_fixed_items = floor(num_items * percentage_of_fixed_items);

	// for (int i = 0; i < num_fixed_items; ++i)
	// {
	// 	int pos = rand() % num_items;
	// 	// iterate until find an item that was not yet fixed ( with group != -1).
	// 	while (items_groups[pos] == -1)
	// 		pos = (pos + 1) % num_items;

	// 	items_groups[pos] = -1;
	// }

	// std::cout << "num fixed: " << num_fixed_items << std::endl;

	int num_vehicles = vehicles_capacities.size();

	for (int k = 0; k < num_items; ++k)
		std::cout << "item: " << k << ": " << "priority: " << items_priorities[k] << " | group: " << items_groups[k] << std::endl;

	for (int k = 0; k < num_vehicles; ++k)
		std::cout << "vehicle: " << k << ": " << "capacity: " << vehicles_capacities[k] << std::endl;

	for (int i = 0; i < num_stacks; ++i)
	{
		std::cout << "stack " << i << ":";
		for (auto &j : stacks[i])
		{
			std::cout << " " << j;
		}
		std::cout << std::endl;
	}
	// for (int i = 0; i < num_vertices; ++i)
	// 	for (int j = i + 1; j < num_vertices; ++j)
	// 		graph_->AddEdge(i, j, round_decimals(euclidian_distance(vertices_info[i].coordinates_, vertices_info[j].coordinates_), 2));

	// write new instance file.

	std::fstream file_out;
	file_out.open(new_instance_file.c_str(), std::fstream::out);
	if (!file_out.is_open())
	{
		std::cout << "Could not open file" << std::endl;
		throw 1;
		return;
	}

	file_out << num_stacks << " " << num_vehicles << " " << num_items << " " << number_of_item_groups << std::endl;

	for (int vehicle = 0; vehicle < num_vehicles; ++vehicle)
	{
		file_out << vehicles_capacities[vehicle];
		vehicle == num_vehicles - 1 ? file_out << std::endl : file_out << " ";
	}

	for (int stack = 0; stack < num_stacks; ++stack)
	{
		file_out << stacks[stack].size();
		for (auto item : stacks[stack])
			file_out << " " << items_groups[item];
		file_out << std::endl;
	}

	file_out.close();
}

void generateInstancesFromBlockRelocationInstances()
{
	std::string input_dir_path = "/home/lucas/Downloads/brp-instances-caserta-etal-2012/CRPTestcases_Caserta/selected_instances/";
	std::string output_dir_path = "../instances/";
	std::vector<std::string> instances = {
		"data10-6-37",
		"data20-20-6",
		"data30-30-32",
		// "data40-40-3",
		"data5-10-36",
		"data5-4-4",
		"data5-6-10",
		"data5-7-37",
		"data5-9-26",
		"data6-6-19"};
	const std::vector<double> percentage_of_items_serviced_by_fleet_vec = {0.25, 0.5};
	const std::vector<double> number_of_item_groups_vec = {15};

	for (auto instance_name : instances)
		for (auto percentage : percentage_of_items_serviced_by_fleet_vec)
			for (auto num_groups : number_of_item_groups_vec)
				generateInstancesFromBlockRelocationInstancesIter(input_dir_path, output_dir_path, instance_name, percentage, num_groups);
}

static const struct option longOpts[] = {
	{"solution-dir", required_argument, NULL, 'a'},
	{"solve-exact", no_argument, NULL, 'b'},
	{"ks-cluster-buckets-by-item-group", no_argument, NULL, 'c'},
	{"reformulate", no_argument, NULL, 'd'},
	{"vehicle-sequencing-model", no_argument, NULL, 'e'},
	{"vehicle-slots-model", no_argument, NULL, 'f'},
	{"time-limit", required_argument, NULL, 'g'},
	{"instance", required_argument, NULL, 'h'},
	{"symmetry-breaking", no_argument, NULL, 'i'},
	{"item-sequencing-model", no_argument, NULL, 'j'},
	{"CCCs", no_argument, NULL, 'm'},
	{"solve-relaxed", no_argument, NULL, 'o'},
	{"solve-kernel-search", no_argument, NULL, 'u'},
	{"ks-max-size-bucket", required_argument, NULL, 'v'},
	{"ks-min-time-limit", required_argument, NULL, 'w'},
	{"ks-max-time-limit", required_argument, NULL, 'x'},
	{"ks-decay-factor", required_argument, NULL, 'y'},
	{"seed", required_argument, NULL, 'A'},
	{NULL, no_argument, NULL, 0}};

void ParseArgumentsAndRun(int argc, char *argv[])
{
	std::string instance, folder, file_name, dir_solutions;
	int c;
	int seed = 0;
	bool solve_exact = false, vehicle_sequencing_model = false, vehicle_slots_model = false, item_sequencing_model = false;
	bool solve_relaxed = false, add_valid_inequalities = false;
	double time_limit = -1.0, original_time_limit = -1.0;
	bool force_use_all_vehicles = false;
	bool solve_kernel_search = false;
	int ks_max_size_bucket = K_KS_MAX_SIZE_BUCKET, ks_min_time_limit = K_KS_MIN_TIME_LIMIT, ks_max_time_limit = K_KS_MAX_TIME_LIMIT;
	double ks_decay_factor = K_KS_DECAY_FACTOR_TIME_LIMIT;
	bool reformulate = false;
	bool add_symmetry_breaking = false;
	bool cluster_buckets_by_item_group = false;

	while ((c = getopt_long(argc, argv, "g:h:v:w:x:y:A:", longOpts, NULL)) != -1)
	{
		switch (c)
		{
		case 'a':
			dir_solutions = std::string(optarg);
			break;
		case 'b':
			solve_exact = true;
			break;
		case 'c':
			cluster_buckets_by_item_group = true;
			break;
		case 'd':
			reformulate = true;
			break;
		case 'e':
			vehicle_sequencing_model = true;
			break;
		case 'f':
			vehicle_slots_model = true;
			break;
		case 'g':
			if (optarg)
				original_time_limit = std::atoi(optarg);
			break;
		case 'h':
			if (optarg)
				instance = std::string(optarg);
			break;
		case 'i':
			add_symmetry_breaking = true;
			// std::cout << "add symmetry breaking" << std::endl;
			break;
		case 'j':
			item_sequencing_model = true;
			break;
		case 'm':
			add_valid_inequalities = true;
			// std::cout << "add valid inequalities" << std::endl;
			break;
		case 'o':
			solve_relaxed = true;
			break;
		case 'u':
			solve_kernel_search = true;
			break;
		case 'v':
			ks_max_size_bucket = std::atoi(optarg);
			break;
		case 'w':
			ks_min_time_limit = std::atoi(optarg);
			break;
		case 'x':
			ks_max_time_limit = std::atoi(optarg);
			break;
		case 'y':
			ks_decay_factor = std::atof(optarg);
			break;
		case 'A':
			if (optarg)
				seed = std::atoi(optarg);
			break;
		}
	}

	srand(seed);

	if (!solve_exact && !solve_kernel_search && !solve_relaxed)
		throw 2;
	if (((solve_exact) && !vehicle_sequencing_model && !vehicle_slots_model && !item_sequencing_model) ||
		(solve_exact && solve_kernel_search) ||
		(solve_relaxed && solve_kernel_search))
		throw 3;
	if ((vehicle_sequencing_model && vehicle_slots_model) || (vehicle_sequencing_model && item_sequencing_model) || (item_sequencing_model && vehicle_slots_model))
		throw 4;

	if (!dir_solutions.empty())
		dir_solutions += "//";

	split_file_path(instance, folder, file_name);

	std::cout << "* " << file_name << std::endl;
	Timestamp *ti = NewTimestamp();
	Timer *timer = GetTimer();
	timer->Clock(ti);

	Instance inst(folder, file_name);

	// std::cout << inst << std::endl;

	std::unique_ptr<Model> model;
	// VehicleSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relaxed);
	// ItemSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relaxed);
	// VehicleSlotsModel model(inst, reformulate, add_symmetry_breaking, solve_relaxed);

	if (vehicle_sequencing_model)
		model = std::make_unique<VehicleSequencingModel>(inst, reformulate, add_symmetry_breaking, solve_relaxed);

	if (vehicle_slots_model)
		model = std::make_unique<VehicleSlotsModel>(inst, reformulate, add_symmetry_breaking, solve_relaxed);

	if (item_sequencing_model)
		model = std::make_unique<ItemSequencingModel>(inst, reformulate, add_symmetry_breaking, solve_relaxed);

	Solution<double>
		solution;

	if (solve_relaxed)
	{
		time_limit = -1;
	}
	else
	{
		if (!double_equals(original_time_limit, -1))
			time_limit = std::max(0.0, original_time_limit - timer->CurrentElapsedTime(ti));
	}

	std::list<UserCut *> *root_cuts = nullptr;
	if ((add_valid_inequalities) && (!solve_relaxed))
	{
		// std::cout << "solve relaxed to separate cuts" << std::endl;
		// fill root_cuts by solving relaxed model with separation of cuts.
		auto relaxed_model = model->getClone(true);
		root_cuts = new std::list<UserCut *>();
		relaxed_model->optimize(time_limit, true, nullptr, root_cuts, solution);
		// relaxed_model->exportModel("relaxed.lp");
		delete relaxed_model;
		relaxed_model = nullptr;
	}

	// int max_bucket_size = std::ceil(inst.num_items() / inst.num_stacks());

	solution.pre_processing_time_ = timer->CurrentElapsedTime(ti);
	// std::cout << "ja usou " << solution.pre_processing_time_ << "s" << std::endl;
	if (!double_equals(original_time_limit, -1))
	{
		// std::cout << original_time_limit << " - " << solution.pre_processing_time_ << " = ";
		time_limit = std::max(0.0, original_time_limit - solution.pre_processing_time_);
		// std::cout << time_limit << std::endl;
	}

	// std::cout << "time_limit " << time_limit << "s" << std::endl;
	if (solve_exact || solve_relaxed)
	{
		auto algo = Solution<double>::GenerateAlgorithmName(solve_relaxed, add_valid_inequalities, add_symmetry_breaking, vehicle_sequencing_model, item_sequencing_model, vehicle_slots_model, reformulate);
		// solve the main problem (either relaxed or integer).
		if (model->optimize(time_limit, add_valid_inequalities, root_cuts, nullptr, solution))
			model->fillSolution(solution, std::nullopt);

		if (solve_relaxed)
		{
			std::cout << "Bound: " << solution.lp_ << std::endl;
		}
		else
		{
			std::cout << "Bound: " << solution.lb_ << std::endl
					  << "num items: " << solution.num_items_loaded_ << std::endl
					  << "num unproductive moves: " << solution.num_unproductive_moves_ << std::endl;
		}

		solution.write_to_file(algo, dir_solutions, file_name);
		// std::cout << "cuts added root: " << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
		// std::cout << "opt: " << model->getObjValue() << std::endl;
	}
	else if (solve_kernel_search)
	{
		// solve the main problem (either relaxed or integer).
		auto fleet = inst.fleet();
		int total_capacity = 0;
		for (auto vehicle : fleet)
			total_capacity += vehicle->capacity();

		int max_bucket_size = std::ceil(total_capacity / fleet.size());

		KernelSearch ks(inst);
		auto heuristic_solution = ks.Run(*model, root_cuts, max_bucket_size, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, cluster_buckets_by_item_group, true, time_limit);

		auto algo = KSHeuristicSolution::GenerateFileName(add_valid_inequalities, add_symmetry_breaking, vehicle_sequencing_model, item_sequencing_model, vehicle_slots_model, reformulate, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, cluster_buckets_by_item_group);

		std::cout << "Bound: " << heuristic_solution->lb_ << std::endl
				  << "num items: " << heuristic_solution->num_items_loaded_ << std::endl
				  << "num unproductive moves: " << heuristic_solution->num_unproductive_moves_ << std::endl;

		heuristic_solution->WriteToFile(inst, algo, dir_solutions, file_name);

		delete heuristic_solution;
		heuristic_solution = nullptr;
	}

	// if (export_model)
	// 	model->exportModel("original.lp");

	solution.total_time_ = timer->CurrentElapsedTime(ti);
	DeleteCuts(root_cuts);
	delete root_cuts;
	root_cuts = nullptr;

	delete ti;
	ti = nullptr;
	DeleteTimer();
}

int main(int argc, char *argv[])
{
	std::string input_dir_path = "/home/lucas/Downloads/brp-instances-caserta-etal-2012/CRPTestcases_Caserta/selected_instances/";
	std::string output_dir_path = "../../instances/";

	// generateInstancesFromBlockRelocationInstancesIter(input_dir_path, output_dir_path, "data40-40-4", 0.75, 10);
	generateInstancesFromBlockRelocationInstances();

	return 0;
	try
	{
		ParseArgumentsAndRun(argc, argv);
		// // std::vector<std::shared_ptr<Item>> items;
		// // std::vector<std::shared_ptr<Vehicle>> fleet;

		// // items.push_back(std::make_shared<Item>(0, true));
		// // items.push_back(std::make_shared<Item>(0, true));
		// // items.push_back(std::make_shared<Item>(2, false));
		// // items.push_back(std::make_shared<Item>(0, true));
		// // items.push_back(std::make_shared<Item>(1, true));

		// // items.push_back(std::make_shared<Item>(2, false));
		// // items.push_back(std::make_shared<Item>(1, true));

		// // items.push_back(std::make_shared<Item>(0, true));
		// // items.push_back(std::make_shared<Item>(0, true));
		// // items.push_back(std::make_shared<Item>(1, true));
		// // items.push_back(std::make_shared<Item>(1, true));
		// // items.push_back(std::make_shared<Item>(0, true));
		// // items.push_back(std::make_shared<Item>(2, false));
		// // items.push_back(std::make_shared<Item>(0, true));

		// // fleet.push_back(std::make_shared<Vehicle>(4));
		// // fleet.push_back(std::make_shared<Vehicle>(3));
		// // fleet.push_back(std::make_shared<Vehicle>(3));

		// // Matrix<int> precedence_matrix(items.size(), items.size(), 0);

		// // for (int i = 0; i <= 3; ++i)
		// // 	precedence_matrix[i][i + 1] = 1;

		// // for (int i = 5; i <= 7; ++i)
		// // 	precedence_matrix[i][i + 1] = 1;

		// // for (int i = 9; i <= 10; ++i)
		// // 	precedence_matrix[i][i + 1] = 1;

		// // precedence_matrix[12][13] = 1;

		// // Instance inst(items, precedence_matrix, fleet);

		// Instance inst("../../instances/", "data6-10-24-c0.75-p10.txt");
		// // Instance inst("../../instances/", "instance");

		// std::cout << inst << std::endl;

		// bool reformulate = true;
		// bool add_symmetry_breaking = false;
		// bool solve_relax = false;
		// bool find_root_cuts = false;
		// int time_limit = 1000;
		// // VehicleSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relax);
		// // ItemSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relax);
		// VehicleSlotsModel model(inst, reformulate, add_symmetry_breaking, solve_relax);

		// Solution<double>
		// 	solution;

		// std::list<UserCut *> *root_cuts = nullptr;
		// if ((find_root_cuts) && (!solve_relax))
		// {
		// 	// fill root_cuts by solving relaxed model with separation of cuts.
		// 	auto relaxed_model = model.getClone(true);
		// 	root_cuts = new std::list<UserCut *>();
		// 	relaxed_model->optimize(-1.0, true, nullptr, root_cuts, solution);
		// 	relaxed_model->exportModel("relaxed.lp");
		// 	delete relaxed_model;
		// 	relaxed_model = nullptr;
		// }

		// // int max_bucket_size = std::ceil(inst.num_items() / inst.num_stacks());

		// auto fleet = inst.fleet();
		// int total_capacity = 0;
		// for (auto vehicle : fleet)
		// 	total_capacity += vehicle->capacity();

		// int max_bucket_size = std::ceil(total_capacity / fleet.size());

		// KernelSearch ks(inst);
		// ks.Run(model, root_cuts, max_bucket_size, 30, 60, 0.9, false, true);

		// // solve the main problem (either relaxed or integer).
		// if (model.optimize(time_limit, find_root_cuts, root_cuts, nullptr, solution))
		// 	model.fillSolution(solution);

		// model.exportModel("original.lp");

		// std::cout << "cuts added root: " << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
		// std::cout << "opt: " << model.getObjValue() << std::endl;

		// DeleteCuts(root_cuts);
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
}