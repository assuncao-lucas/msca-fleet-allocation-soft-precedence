
#include "src/instance.h"

int main()
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

	fleet.push_back(std::make_shared<Vehicle>(3));
	fleet.push_back(std::make_shared<Vehicle>(3));
	fleet.push_back(std::make_shared<Vehicle>(3));

	Matrix<int> precedence_matrix(items.size(), items.size(), 0);

	for (int i = 0; i <= 3; ++i)
		precedence_matrix[i][i + 1] = 1;

	for (int i = 5; i <= 7; ++i)
		precedence_matrix[i][i + 1] = 1;

	for (int i = 9; i <= 10; ++i)
		precedence_matrix[i][i + 1] = 1;

	Instance inst(items, precedence_matrix, fleet);

	std::cout << inst << std::endl;

	return 0;
}