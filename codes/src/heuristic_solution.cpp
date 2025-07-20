#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include "heuristic_solution.h"

VertexStatus::VertexStatus()
{
	selected_ = false;
	route_ = -1;
}

VertexStatus::~VertexStatus()
{
}

HeuristicSolution::~HeuristicSolution()
{
}

HeuristicSolution::HeuristicSolution(int num_vertices, int num_arcs, int num_routes)
{
	HeuristicSolution::Reset(num_vertices, num_arcs, num_routes);
}

void HeuristicSolution::Reset(int num_vertices, int num_arcs, int num_routes)
{
	// profits_sum_ = 0;
	// num_vertices_ = num_vertices;
	// num_arcs_ = num_arcs;
	// num_routes_ = num_routes;
	// routes_vec_ = std::vector<Route>(num_routes);
	// vertex_status_vec_ = std::vector<VertexStatus>(num_vertices);
	// is_infeasible_ = false;
	// is_feasible_ = false;
	// is_optimal_ = false;
	// (unvisited_vertices_).clear();
	// total_time_spent_ = 0.0;

	// bitset_arcs_ = boost::dynamic_bitset<>(num_arcs, 0);
	// bitset_vertices_ = boost::dynamic_bitset<>(num_vertices, 0);

	// for (int i = 1; i < num_vertices; ++i)
	// {
	// 	(unvisited_vertices_).push_front(i);
	// 	VertexStatus *status = &((vertex_status_vec_)[i]);
	// 	status->selected_ = false;
	// 	status->route_ = -1;
	// 	status->pos_ = (unvisited_vertices_).begin();
	// }
}

bool HeuristicSolution::operator==(HeuristicSolution &other)
{
	return (bitset_arcs_ == other.bitset_arcs_);
}

std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol)
{
	// if (sol.is_infeasible_)
	// 	out << "[INFEASIBLE]";
	// if (sol.is_optimal_)
	// 	out << "[OPTMAL]";
	// out << "[profits_sum: " << sol.profits_sum_ << "]" << std::endl;
	// for (int i = 0; i < sol.num_routes_; ++i)
	// 	out << ((sol.routes_vec_)[i]);
	// out << "unvisited: ";
	// for (auto it = (sol.unvisited_vertices_).begin(); it != (sol.unvisited_vertices_).end(); ++it)
	// 	out << *it << " ";
	// out << std::endl;

	// out << "vertex status: " << std::endl;
	// for (int vertex = 1; vertex < sol.vertex_status_vec_.size(); ++vertex) // skip the zero, cause it's not supposed to be in any list.
	// {
	// 	auto status = &(sol.vertex_status_vec_[vertex]);
	// 	out << vertex << " " << status->selected_ << " " << status->route_ << " " << *status->pos_ << std::endl;
	// }
	return out;
}

void HeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const
{
	// std::fstream file;
	// std::string path = "..//solutions//";
	// path.append(folder);
	// // struct stat sb;
	// // if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	// path.append("s_");
	// path.append(algo);
	// path.append("_");
	// path.append(file_name);
	// // std::cout << path << std::endl;

	// file.open(path.c_str(), std::fstream::out | std::fstream::app);

	// file << std::setprecision(5) << std::fixed;

	// file << K_FILE_DELIMITER << std::endl
	// 	 << "Total profit sum: " << profits_sum_ << std::endl
	// 	 << "Num routes: ";
	// if ((is_infeasible_) || (!(is_feasible_)))
	// 	file << "0" << std::endl;
	// else
	// {
	// 	file << num_routes_ << std::endl;

	// 	for (int i = 0; i < num_routes_; i++)
	// 	{
	// 		const Route &curr_route = (routes_vec_)[i];
	// 		file << curr_route.sum_profits_ << " " << curr_route.time_ << " ";

	// 		for (auto it = (curr_route.vertices_).begin(); it != (curr_route.vertices_).end(); ++it)
	// 		{
	// 			file << instance.getOriginalVertexPosition(*it) << " ";
	// 		}

	// 		// file << num_vertices_ - 1;
	// 		file << std::endl;
	// 	}
	// }

	// file.close();
}

void HeuristicSolution::ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	// const Graph *graph = inst.graph();
	// int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();

	// Reset(num_vertices, num_arcs, inst.num_vehicles());
	// Route *curr_route = nullptr;
	// int num_routes = 0;
	// std::fstream input;
	// std::string path = "..//solutions//";
	// path.append(folder);
	// // struct stat sb;
	// // if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	// path.append("s_");
	// path.append(algo);
	// path.append("_");
	// path.append(file_name);
	// // std::cout << path << std::endl;

	// input.open(path.c_str(), std::fstream::in);
	// std::string line;

	// if (!(input.is_open()))
	// 	throw "Error opening file 1";
	// while (!(input.eof()))
	// {
	// 	getline(input, line);
	// 	if (line == K_FILE_DELIMITER)
	// 	{
	// 		std::stringstream s_profits_sum, s_num_routes;

	// 		getline(input, line);
	// 		size_t pos = line.find_first_of(":");
	// 		s_profits_sum << line.substr(pos + 2);
	// 		s_profits_sum >> profits_sum_;
	// 		// std::cout << profits_sum_ << std::endl;

	// 		getline(input, line);
	// 		pos = line.find_first_of(":");
	// 		s_num_routes << line.substr(pos + 2);
	// 		s_num_routes >> num_routes;

	// 		if ((num_routes == 0) || (num_routes != num_routes_))
	// 			is_feasible_ = false;
	// 		else
	// 		{
	// 			is_feasible_ = true;
	// 			int pre_vertex = 0, curr_vertex = 0, vertex_pos_before_reordering = 0;
	// 			VertexStatus *status = nullptr;
	// 			for (int i = 0; i < num_routes; i++)
	// 			{
	// 				curr_route = &((routes_vec_)[i]);
	// 				getline(input, line);
	// 				std::stringstream s_route(line);
	// 				// std::cout << s_route.str() << std::endl;
	// 				s_route >> curr_route->sum_profits_ >> curr_route->time_;
	// 				// std::cout << curr_route->sum_profits_ << " " << curr_route->time_ << std::endl;
	// 				curr_vertex = pre_vertex = 0;
	// 				while (s_route >> vertex_pos_before_reordering)
	// 				{
	// 					// in the file/solution, is saved with original positions of vertices (before reordering!)
	// 					curr_vertex = inst.getReorderedVertexPosition(vertex_pos_before_reordering);
	// 					status = &((vertex_status_vec_)[curr_vertex]);

	// 					// remove from list of unvisited_vertices
	// 					(unvisited_vertices_).erase(status->pos_);

	// 					// adds vertex to route
	// 					status->selected_ = true;
	// 					status->route_ = i;
	// 					status->pos_ = (curr_route->vertices_).insert((curr_route->vertices_).end(), curr_vertex);

	// 					//(bitset_arcs_)[graph->pos(pre_vertex,curr_vertex)] = 1;
	// 					pre_vertex = curr_vertex;
	// 				}

	// 				curr_vertex = num_vertices - 1;
	// 				//(bitset_arcs_)[graph->pos(pre_vertex,curr_vertex)] = 1;
	// 			}
	// 		}
	// 		break;
	// 	}
	// }

	// if (!(is_infeasible_) && (is_feasible_) && (unvisited_vertices_).empty())
	// 	is_optimal_ = true;

	// input.close();
}

KSHeuristicSolution::KSHeuristicSolution(int num_vertices, int num_arcs, int num_routes) : HeuristicSolution(num_vertices, num_arcs, num_routes)
{
	found_x_integer_ = false;
	time_spent_building_kernel_buckets_ = 0.0;
}

KSHeuristicSolution::KSHeuristicSolution() : HeuristicSolution()
{
	found_x_integer_ = false;
	time_spent_building_kernel_buckets_ = 0.0;
}

KSHeuristicSolution::~KSHeuristicSolution()
{
}

void KSHeuristicSolution::Reset(int num_vertices, int num_arcs, int num_routes)
{
	HeuristicSolution::Reset(num_vertices, num_arcs, num_routes);

	found_x_integer_ = false;
	time_spent_building_kernel_buckets_ = 0.0;
}

std::string KSHeuristicSolution::GenerateFileName(Formulation formulation, int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool feasibility_emphasis)
{
	std::string formulation_name;
	switch (formulation)
	{
	case Formulation::baseline:
		formulation_name = "baseline";
		break;
	case Formulation::single_commodity:
		formulation_name = "csc";
		break;
	default:
		throw "invalid formulation";
	}

	std::stringstream ss_decay_factor;
	ss_decay_factor << std::fixed << std::setprecision(2) << ks_decay_factor;
	std::string file_name;
	file_name += formulation_name + "_ks_b" + std::to_string(ks_max_size_bucket) + "_[" + std::to_string(ks_max_time_limit) + "," + std::to_string(ks_min_time_limit) + "]_d" + ss_decay_factor.str();
	if (feasibility_emphasis)
		file_name += "_feas";
	return file_name;
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
		file << "STATUS: FOUND ARC INTEGER FEASIBLE" << std::endl;
	else
		file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;

	file << "time building kernel and buckets (s): " << time_spent_building_kernel_buckets_ << std::endl;
	file << "total time (s): " << total_time_spent_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
}