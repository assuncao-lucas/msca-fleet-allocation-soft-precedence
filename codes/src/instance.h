#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include "src/matrix.hpp"
#include "src/graph.h"

class Item
{
private:
	bool available_for_transport_ = false;
	int group_;

public:
	explicit Item(int group, bool available_for_transport) : group_(group), available_for_transport_(available_for_transport)
	{
	}
	virtual ~Item() {}

	int group() const
	{
		return group_;
	}

	bool available_for_transport() const
	{
		return available_for_transport_;
	}
};

class Vehicle
{
private:
	int capacity_;

public:
	explicit Vehicle(int capacity) : capacity_(capacity)
	{
	}

	virtual ~Vehicle() {}

	int capacity() const
	{
		return capacity_;
	}
};

class Instance
{
private:
	std::vector<std::shared_ptr<Item>> items_;										 // complete set of items.
	std::vector<int> items_for_transport_;											 // indexes of items that can be transported.
	std::vector<int> items_fixed_;													 // indexes of items fixed.
	std::unordered_map<int, std::unordered_set<int>> items_for_transport_per_group_; // partition of items_for_transport per group.

	Matrix<int> precedence_matrix_;

	std::vector<std::list<int>> successors_for_transport_per_item_;
	std::vector<std::list<int>> successors_fixed_per_item_;

	std::vector<std::shared_ptr<Vehicle>> fleet_;

	bool PropagatePrecedence(); // guarantees transitivity of precedence. if a -> b and b -> c, then a -> c.
	bool PropagatePrecedenceIter(std::stack<std::pair<int, bool>> &main_stack, std::vector<bool> &visited, std::vector<bool> &in_stack);

public:
	explicit Instance() = default;
	explicit Instance(std::string dir_path, std::string file_name);
	explicit Instance(std::vector<std::shared_ptr<Item>> items, Matrix<int> precedence_matrix, std::vector<std::shared_ptr<Vehicle>> fleet);
	virtual ~Instance() = default;

	void FillInstanceFromFile(std::string dir_path, std::string file_name, double service_time_deviation);
	void WriteToFile(std::string dir_path, std::string curr_file);

	friend std::ostream &operator<<(std::ostream &out, const Instance &instance);

	// getters.
	const std::vector<std::shared_ptr<Item>> &items() const
	{
		return items_;
	}

	const std::vector<int> &items_for_transport() const
	{
		return items_for_transport_;
	}

	const std::vector<int> &items_fixed() const
	{
		return items_fixed_;
	}

	const std::unordered_map<int, std::unordered_set<int>> &items_for_transport_per_group() const
	{
		return items_for_transport_per_group_;
	}

	const Matrix<int> &precedence_matrix() const
	{
		return precedence_matrix_;
	}

	int num_successors_for_transport_per_item(int item) const
	{
		return successors_for_transport_per_item_[item].size();
	}

	int num_successors_fixed_per_item(int item) const
	{
		return successors_fixed_per_item_[item].size();
	}

	const std::vector<std::list<int>> &successors_for_transport_per_item() const
	{
		return successors_for_transport_per_item_;
	}

	const std::vector<std::list<int>> &successors_fixed_per_item() const
	{
		return successors_fixed_per_item_;
	}

	const std::vector<std::shared_ptr<Vehicle>> fleet() const
	{
		return fleet_;
	}

	const int num_items() const
	{
		return items_.size();
	}

	const int num_items_for_transport() const
	{
		return items_for_transport_.size();
	}

	const int num_vehicles() const
	{
		return fleet_.size();
	}

	// setters.
	void set_items(std::vector<std::shared_ptr<Item>> items);
	void set_precedence_matrix(Matrix<int> &precedence_matrix);
	void set_fleet(std::vector<std::shared_ptr<Vehicle>> fleet);
};
