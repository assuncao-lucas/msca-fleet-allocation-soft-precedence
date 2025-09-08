#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include "src/instance.h"
#include "src/general.h"

void Instance::FillInstanceFromFile(std::string dir_path, std::string file_name)
{
    std::string curr_file = dir_path;
    curr_file += file_name;

    // std::cout << dir_path << " " << file_name << std::endl;

    std::fstream file;
    file.open(curr_file.c_str(), std::fstream::in);
    if (!file.is_open())
    {
        std::cout << "Could not open file" << std::endl;
        throw 1;
        return;
    }

    int num_stacks = 0, num_items = 0, group = -1, num_items_in_stack = 0, num_vehicles = 0, capacity = 0, num_groups = 0;
    int min_priority = std::numeric_limits<int>::max(), max_priority = -1;

    file >> num_stacks >> num_vehicles >> num_items >> num_groups;

    std::vector<std::shared_ptr<Item>> items;
    std::vector<std::shared_ptr<Vehicle>> fleet;
    Matrix<int> precedence_matrix(num_items, num_items, 0);

    for (int vehicle = 0; vehicle < num_vehicles; ++vehicle)
    {
        file >> capacity;
        fleet.push_back(std::make_shared<Vehicle>(capacity));
    }

    int item_counter = 0;
    for (int i = 0; i < num_stacks; ++i)
    {
        file >> num_items_in_stack;
        for (int j = 0; j < num_items_in_stack; ++j)
        {
            file >> group;
            // make the last group as fixed (not available for transport).
            items.push_back(std::make_shared<Item>(group, group == num_groups - 1 ? false : true));

            if (j < num_items_in_stack - 1)
                precedence_matrix[item_counter + 1][item_counter] = 1;
            ++item_counter;
        }
    }

    file.close();

    set_items(items);
    set_precedence_matrix(precedence_matrix);
    set_fleet(fleet);
    num_stacks_ = num_stacks;
}

Instance::Instance(std::string dir_path, std::string file_name)
{
    FillInstanceFromFile(dir_path, file_name);
}

Instance::Instance(std::vector<std::shared_ptr<Item>> items, Matrix<int> precedence_matrix, std::vector<std::shared_ptr<Vehicle>> fleet, int num_stacks)
{
    set_items(items);
    set_precedence_matrix(precedence_matrix);
    set_fleet(fleet);
    num_stacks_ = num_stacks;
}

std::ostream &
operator<<(std::ostream &out, const Instance &instance)
{
    const auto &items = instance.items();

    out << "Items:" << std::endl;

    for (int i = 0; i < items.size(); ++i)
    {
        auto &item = items[i];
        out << "id: " << i << " group: " << item->group() << " available for transp.: " << item->available_for_transport() ? "yes" : "no";
        out << std::endl;
    }

    out << "Items for transport per group:" << std::endl;

    for (auto &[group, items_of_group] : instance.items_for_transport_per_group())
    {
        out << "group " << group << ":" << std::endl;
        for (auto item : items_of_group)
            out << " " << item;
        out << std::endl;
    }

    out << "Vehicles:" << std::endl;
    const auto &fleet = instance.fleet();

    for (int i = 0; i < fleet.size(); ++i)
    {
        auto &vehicle = fleet[i];
        out << "id: " << i << " capacity: " << vehicle->capacity() << std::endl;
    }

    out << "Precedence matrix:" << std::endl;
    const auto &precedence = instance.precedence_matrix();

    for (int i = 0; i < items.size(); ++i)
    {
        out << i << ":";
        for (int j = 0; j < items.size(); ++j)
            if (precedence[i][j])
                out << " " << j;

        out << std::endl;
    }

    out << "Successors:" << std::endl;

    for (int i = 0; i < items.size(); ++i)
    {
        out << i << ":";
        for (auto &successor : instance.successors_for_transport_per_item()[i])
            out << " " << successor;
        for (auto &successor : instance.successors_fixed_per_item()[i])
            out << " " << successor;

        out << std::endl;
    }

    out << "Predecessors for transport:" << std::endl;

    for (int i = 0; i < items.size(); ++i)
    {
        out << i << ":";
        for (auto &predecessor : instance.predecessors_for_transport_per_item()[i])
            out << " " << predecessor;

        out << std::endl;
    }

    return out;
}

bool Instance::PropagatePrecedenceIter(std::stack<std::pair<int, bool>> &main_stack, std::vector<bool> &visited, std::vector<bool> &in_stack)
{
    size_t num_items = items_.size();

    while (!main_stack.empty())
    {
        auto [v, return_mark] = main_stack.top();
        main_stack.pop();

        if (return_mark)
        {
            in_stack[v] = false; // Ending visitation of item.
                                 // update precedence relation to all items that precede v. Precisely, if an item is a sucessor of v, then it should also be a sucessor of any precedent of v.
            for (int precedent_v = 0; precedent_v < num_items; ++precedent_v)
            {
                if (precedence_matrix_[precedent_v][v])
                {
                    for (int successor_v = 0; successor_v < num_items; ++successor_v)
                        precedence_matrix_[precedent_v][successor_v] = std::max(precedence_matrix_[precedent_v][successor_v], precedence_matrix_[v][successor_v]);
                }
            }
            continue;
        }

        if (!visited[v])
        {
            visited[v] = true;
            in_stack[v] = true;
            main_stack.push({v, true}); // Set it back to stack for return mark.

            // add to stack its sucessors.
            for (int successor = 0; successor < num_items; ++successor)
            {
                if (precedence_matrix_[v][successor])
                {
                    if (!visited[successor])
                    {
                        main_stack.push({successor, false});
                    }
                    else if (in_stack[successor])
                    {
                        // if successor already in stack, means a cycle was detected.
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

bool Instance::PropagatePrecedence()
{
    size_t num_items = items_.size();
    std::vector<bool> visited(num_items, false);
    std::vector<bool> in_stack(num_items, false);
    std::stack<std::pair<int, bool>> main_stack; // tuples of node number and if it's first visit or marks its return.

    for (int v = 0; v < num_items; ++v)
    {
        if (!visited[v])
        {
            main_stack.push({v, false});
            // if find cycle.
            if (!PropagatePrecedenceIter(main_stack, visited, in_stack))
                return false;
        }
    }

    return true;
}

// setters.
void Instance::set_items(std::vector<std::shared_ptr<Item>> items)
{
    items_ = items;

    // update partition of items per group.
    items_for_transport_per_group_.clear();

    for (int i = 0; i < items.size(); ++i)
    {
        const auto &item = items[i];

        if (item->available_for_transport())
        {
            items_for_transport_.push_back(i);
            items_for_transport_per_group_[item->group()].insert(i);
        }
        else
            items_fixed_.push_back(i);
    }
    successors_for_transport_per_item_.resize(items.size());
    predecessors_for_transport_per_item_.resize(items.size());
    successors_fixed_per_item_.resize(items.size());
}

void Instance::set_precedence_matrix(Matrix<int> &precedence_matrix)
{
    precedence_matrix_ = precedence_matrix;

    if (PropagatePrecedence())
    {
        const auto lines = precedence_matrix.lines();
        const auto columns = precedence_matrix.columns();

        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                if (precedence_matrix_[i][j])
                {
                    if (items_[j]->available_for_transport())
                    {
                        successors_for_transport_per_item_[i].push_back(j);
                    }
                    else
                    {
                        successors_fixed_per_item_[i].push_back(j);
                    }

                    if (items_[i]->available_for_transport())
                    {
                        predecessors_for_transport_per_item_[j].push_back(i);
                    }
                }
            }
        }
    }
    else
        throw "Invalid instance - cycle detected";
}

void Instance::set_fleet(std::vector<std::shared_ptr<Vehicle>> fleet)
{
    fleet_ = fleet;
}