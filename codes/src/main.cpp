
#include "src/instance.h"
#include "src/exact/formulations.h"

int main()
{
	try
	{
		std::vector<std::shared_ptr<Item>> items;
		std::vector<std::shared_ptr<Vehicle>> fleet;

		items.push_back(std::make_shared<Item>(0, true));
		items.push_back(std::make_shared<Item>(0, true));
		items.push_back(std::make_shared<Item>(2, false));
		items.push_back(std::make_shared<Item>(0, true));
		items.push_back(std::make_shared<Item>(1, true));

		items.push_back(std::make_shared<Item>(2, false));
		items.push_back(std::make_shared<Item>(1, true));

		items.push_back(std::make_shared<Item>(0, true));
		items.push_back(std::make_shared<Item>(0, true));
		items.push_back(std::make_shared<Item>(1, true));
		items.push_back(std::make_shared<Item>(1, true));
		items.push_back(std::make_shared<Item>(0, true));
		items.push_back(std::make_shared<Item>(2, false));
		items.push_back(std::make_shared<Item>(0, true));

		fleet.push_back(std::make_shared<Vehicle>(4));
		fleet.push_back(std::make_shared<Vehicle>(3));
		fleet.push_back(std::make_shared<Vehicle>(3));

		Matrix<int> precedence_matrix(items.size(), items.size(), 0);

		for (int i = 0; i <= 3; ++i)
			precedence_matrix[i][i + 1] = 1;

		for (int i = 5; i <= 7; ++i)
			precedence_matrix[i][i + 1] = 1;

		for (int i = 9; i <= 10; ++i)
			precedence_matrix[i][i + 1] = 1;

		precedence_matrix[12][13] = 1;

		Instance inst(items, precedence_matrix, fleet);

		std::cout << inst << std::endl;

		bool reformulate = false;
		bool add_symmetry_breaking = true;
		bool solve_relax = true;
		bool find_root_cuts = true;
		bool export_model = true;
		int time_limit = 1000;
		VehicleSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relax, export_model);
		// ItemSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relax, export_model);
		// VehicleSlotsModel model(inst, reformulate, add_symmetry_breaking, solve_relax, export_model);

		Solution<double> solution;

		std::list<UserCut *> *root_cuts = nullptr;
		if ((find_root_cuts) && (!solve_relax))
		{
			// fill root_cuts by solving relaxed model with separation of cuts.
			VehicleSequencingModel relaxed_model(inst, reformulate, add_symmetry_breaking, true, false);
			// ItemSequencingModel relaxed_model(inst, reformulate, add_symmetry_breaking, true, false);
			// VehicleSlotsModel relaxed_model(inst, reformulate, add_symmetry_breaking, true, false);
			root_cuts = new std::list<UserCut *>();
			relaxed_model.optimize(inst, -1.0, true, nullptr, root_cuts, solution);
		}

		// solve the main problem (either relaxed or integer).
		if (model.optimize(inst, time_limit, find_root_cuts, root_cuts, nullptr, solution))
			model.fillSolution(inst, solution);

		std::cout << "cuts added root: " << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

		DeleteCuts(root_cuts);
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