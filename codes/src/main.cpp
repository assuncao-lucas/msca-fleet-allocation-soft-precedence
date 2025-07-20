
#include "src/instance.h"
#include "src/exact/formulations.h"
#include "src/kernel_search/kernel_search.h"

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

		bool reformulate = true;
		bool add_symmetry_breaking = true;
		bool solve_relax = false;
		bool find_root_cuts = true;
		int time_limit = 1000;
		// VehicleSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relax);
		// ItemSequencingModel model(inst, reformulate, add_symmetry_breaking, solve_relax);
		VehicleSlotsModel model(inst, reformulate, add_symmetry_breaking, solve_relax);

		Solution<double>
			solution;

		std::list<UserCut *> *root_cuts = nullptr;
		if ((find_root_cuts) && (!solve_relax))
		{
			// fill root_cuts by solving relaxed model with separation of cuts.
			auto relaxed_model = model.getClone(true);
			root_cuts = new std::list<UserCut *>();
			relaxed_model->optimize(-1.0, true, nullptr, root_cuts, solution);
			relaxed_model->exportModel("relaxed.lp");
			delete relaxed_model;
			relaxed_model = nullptr;
		}

		KernelSearch ks(inst);
		ks.Run(model, root_cuts, 5, 100, 200, 0.1, false, true);

		// solve the main problem (either relaxed or integer).
		if (model.optimize(time_limit, find_root_cuts, root_cuts, nullptr, solution))
			model.fillSolution(solution);

		model.exportModel("original.lp");

		std::cout << "cuts added root: " << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
		std::cout << "opt: " << model.getObjValue() << std::endl;

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