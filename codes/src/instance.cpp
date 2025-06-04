#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include "src/instance.h"
#include "src/general.h"

void Instance::FillInstanceFromFile(std::string dir_path, std::string file_name, double service_time_deviation)
{
    // std::string curr_file = dir_path;
    // curr_file.append(file_name);
    // Graph *graph = nullptr;

    // // std::cout << dir_path << " " << file_name << std::endl;

    // std::fstream file;
    // file.open(curr_file.c_str(), std::fstream::in);
    // if (!file.is_open())
    // {
    //     std::cout << "Could not open file" << std::endl;
    //     throw 1;
    //     return;
    // }

    // int num_vertices = 0, mandatory_iter = 0;

    // file >> num_vertices >> num_mandatory_;

    // // std::cout << num_vertices << " " << num_mandatory << std::endl;

    // for (int i = 0; i < num_mandatory_; ++i)
    // {
    //     file >> mandatory_iter;
    //     mandatory_list_.push_back(mandatory_iter);
    // }

    // std::vector<std::pair<double, double>> coordinates(num_vertices);
    // Graph::VertexInfo *vertices_info = new Graph::VertexInfo[num_vertices];

    // // fill Vertices info.
    // double garbage = 0.0;
    // for (int i = 0; i < num_vertices; ++i)
    // {
    //     file >> vertices_info[i].coordinates_.first >> vertices_info[i].coordinates_.second >> vertices_info[i].profit_ >> vertices_info[i].decay_ratio_ >> vertices_info[i].nominal_service_time_ >> garbage;
    //     vertices_info[i].dev_service_time_ = round_decimals(vertices_info[i].nominal_service_time_ * service_time_deviation, 2);
    // }
    // assert(vertices_info[0].profit_ == 0);
    // assert(double_equals(vertices_info[0].decay_ratio_, 0.0));
    // assert(double_equals(vertices_info[0].nominal_service_time_, 0.0));
    // assert(double_equals(vertices_info[0].dev_service_time_, 0.0));

    // // read route time limit.
    // file >> limit_;

    // file.close();

    // graph_ = new Graph(num_vertices, vertices_info);

    // for (int i = 0; i < num_vertices; ++i)
    //     for (int j = i + 1; j < num_vertices; ++j)
    //         graph_->AddEdge(i, j, round_decimals(euclidian_distance(vertices_info[i].coordinates_, vertices_info[j].coordinates_), 2));
}

Instance::Instance(std::string dir_path, std::string file_name)
{
    // FillInstanceFromFile(dir_path, file_name, service_time_deviation);
}

Instance::Instance(std::vector<std::shared_ptr<Item>> items, Matrix<int> precedence_matrix, std::vector<std::shared_ptr<Vehicle>> fleet)
{
    set_items(items);
    set_precedence_matrix(precedence_matrix);
    set_fleet(fleet);
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

    return out;
}

void Instance::WriteToFile(std::string folder, std::string curr_file)
{
    // std::fstream output;

    // Graph *graph = graph_;
    // int num_vertices = graph->num_vertices();
    // double limit = limit_;
    // const Graph::VertexInfo *vertices_info = graph->vertices_info();

    // std::string file_path = folder + curr_file;
    // output.open(file_path.c_str(), std::fstream::out);

    // // std::cout << folder << " " << curr_file << std::endl;

    // if (!output.is_open())
    // {
    //     std::cout << "Could not open file" << std::endl;
    //     throw 1;
    //     return;
    // }

    // output << num_vertices << "\t" << num_mandatory_ << std::endl;

    // std::list<int>::iterator last_element = mandatory_list_.end();
    // --last_element;
    // for (std::list<int>::iterator it = mandatory_list_.begin(); it != mandatory_list_.end(); ++it)
    //     it != last_element ? output << map_reordered_vertices_to_original_positions_[*it] << "\t" : output << map_reordered_vertices_to_original_positions_[*it] << std::endl;

    // // output << std::setprecision(2) << std::fixed;

    // for (int i = 0; i < num_vertices; ++i)
    //     output
    //         << vertices_info[i].coordinates_.first << "\t"
    //         << vertices_info[i].coordinates_.second << "\t"
    //         << vertices_info[i].profit_ << "\t"
    //         << vertices_info[i].decay_ratio_ << "\t"
    //         << vertices_info[i].nominal_service_time_ << "\t"
    //         << vertices_info[i].nominal_service_time_ << std::endl;

    // output << limit_;
    // output.close();
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
            for (int j = 0; j < columns; ++j)
                if (precedence_matrix_[i][j])
                    items_[j]->available_for_transport() ? successors_for_transport_per_item_[i].push_back(j) : successors_fixed_per_item_[i].push_back(j);
    }
    else
        throw "Invalid instance - cycle detected";
}

void Instance::set_fleet(std::vector<std::shared_ptr<Vehicle>> fleet)
{
    fleet_ = fleet;
}