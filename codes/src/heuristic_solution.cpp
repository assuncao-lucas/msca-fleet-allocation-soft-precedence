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
	profits_sum_ = 0;
	num_vertices_ = num_vertices;
	num_arcs_ = num_arcs;
	num_routes_ = num_routes;
	routes_vec_ = std::vector<Route>(num_routes);
	vertex_status_vec_ = std::vector<VertexStatus>(num_vertices);
	is_infeasible_ = false;
	is_feasible_ = false;
	is_optimal_ = false;
	(unvisited_vertices_).clear();
	total_time_spent_ = 0.0;

	bitset_arcs_ = boost::dynamic_bitset<>(num_arcs, 0);
	bitset_vertices_ = boost::dynamic_bitset<>(num_vertices, 0);

	for (int i = 1; i < num_vertices; ++i)
	{
		(unvisited_vertices_).push_front(i);
		VertexStatus *status = &((vertex_status_vec_)[i]);
		status->selected_ = false;
		status->route_ = -1;
		status->pos_ = (unvisited_vertices_).begin();
	}
}

bool HeuristicSolution::Do2OptImprovement(const Instance &inst, const int &route)
{
	const auto *graph = inst.graph();
	Route *curr_route = &((routes_vec_)[route]);
	int size_route = (int)((curr_route->vertices_).size());

	// std::cout << "before 2 opt" << std::endl;
	// BuildBitset(inst);
	// CheckCorrectness(inst);

	if (size_route < 2)
		return false;

	std::list<int>::iterator vi_it;
	std::list<int>::iterator vk_it;
	std::list<int>::iterator swap_i, swap_f;
	std::list<int>::iterator pre_vi_it, pos_vk_it, it_init, it_end;
	int vi = 0, vk = 0, pre_vi = 0, pos_vk = 0, v1 = 0, v2 = 0;
	int num_swaps = 0, temp = 0, num_vertices = graph->num_vertices();
	GArc *pre_arc = nullptr, *pos_arc = nullptr, *new_pre_arc = nullptr, *new_pos_arc = nullptr, *curr_arc = nullptr;

	double time_variation = 0.0;

	vi_it = (curr_route->vertices_).begin();
	for (int i = 0; i < size_route - 1; ++i)
	{
		vk_it = vi_it;
		vi = (*vi_it);

		if (i > 0)
		{
			pre_vi_it = vi_it;
			--pre_vi_it;
			pre_vi = (*pre_vi_it);
		}
		else
		{
			pre_vi = 0;
		}

		for (int k = i + 1; k < size_route; ++k)
		{
			++vk_it;
			vk = (*vk_it);

			if (k < size_route - 1)
			{
				pos_vk_it = vk_it;
				++pos_vk_it;
				pos_vk = (*pos_vk_it);
			}
			else
			{
				pos_vk = 0;
			}

			// compute time variation properly!!
			new_pre_arc = (*graph)[pre_vi][vk];
			new_pos_arc = (*graph)[vi][pos_vk];

			if ((new_pre_arc == nullptr) || (new_pos_arc == nullptr))
				continue;

			pre_arc = (*graph)[pre_vi][vi];
			pos_arc = (*graph)[vk][pos_vk];

			// std::cout << "pre_vi: " << pre_vi << " vi: " << vi << " vk:" << vk << " pos_vk: " << pos_vk << std::endl;

			// compute time variation after 2-opt
			time_variation = (new_pre_arc->distance() + new_pos_arc->distance() - pre_arc->distance() - pos_arc->distance());

			it_init = vi_it;
			v2 = *it_init;
			it_end = vk_it;
			bool invalid_swap = false;

			do
			{
				v1 = v2;
				++it_init;
				v2 = *it_init;

				curr_arc = (*graph)[v1][v2];
				time_variation -= curr_arc->distance();
				curr_arc = (*graph)[v2][v1];
				if (curr_arc == nullptr)
				{
					invalid_swap = true;
					break;
				}
				else
					time_variation += curr_arc->distance();
			} while (it_init != it_end);

			// if the resulting route is still feasible, do the actual swap.
			if ((!invalid_swap) && !double_greater(curr_route->time_ + time_variation, inst.limit())) // && double_less(time_variation, 0.0))
			{
				// std::cout << *curr_route << std::endl;
				// std::cout << "após reverse Swap de " << vi << " a " << vk << std::endl;
				//(curr_route->time_) += (time_variation);
				num_swaps = (k - i + 1) / 2;
				swap_i = vi_it;
				swap_f = vk_it;
				for (int cont = 0; cont < num_swaps; ++cont)
				{
					temp = *swap_i;

					*swap_i = *swap_f;
					*swap_f = temp;

					// atualiza os iterators no PoleStatus
					((vertex_status_vec_)[*swap_i]).pos_ = swap_i;
					((vertex_status_vec_)[*swap_f]).pos_ = swap_f;

					++swap_i;
					--swap_f;
				}
				// std::cout << *curr_route << std::endl;
				auto [new_route_sum_profits, new_route_max_duration] = inst.ComputeRouteCosts(curr_route->vertices_);

				// auto [new_route_sum_profits12, new_route_max_duration12] = inst.ComputeRouteCosts(curr_route->vertices_, true);

				// if (!double_equals(new_route_sum_profits, new_route_sum_profits12) || !double_equals(new_route_max_duration, new_route_max_duration12))
				// {
				// 	std::cout << new_route_sum_profits << " x " << new_route_sum_profits12 << std::endl;
				// 	std::cout << new_route_max_duration << " x " << new_route_max_duration12 << std::endl;
				// 	getchar();
				// 	getchar();
				// }
				// std::cout << "new profit " << new_route_sum_profits << std::endl;
				// std::cout << "new duration " << new_route_max_duration << std::endl;

				double profit_variation = new_route_sum_profits - curr_route->sum_profits_;

				// if (!double_equals(new_route_max_duration, curr_route->time_ + time_variation))
				// {
				// 	std::cout << new_route_max_duration << " != " << curr_route->time_ + time_variation << " = " << curr_route->time_ << " + " << time_variation << std::endl;
				// 	getchar();
				// 	getchar();
				// }

				// return true if achieved improvement.
				if (double_greater(profit_variation, 0.0) || (double_equals(profit_variation, 0.0) && double_less(new_route_max_duration, curr_route->time_)))
				{
					curr_route->sum_profits_ = new_route_sum_profits;
					curr_route->time_ = new_route_max_duration;
					profits_sum_ += profit_variation;

					// getchar();
					// getchar();
					return true;
				}

				// otherwise, undo the swap!
				swap_i = vi_it;
				swap_f = vk_it;
				for (int cont = 0; cont < num_swaps; ++cont)
				{
					temp = *swap_i;

					*swap_i = *swap_f;
					*swap_f = temp;

					// atualiza os iterators no PoleStatus
					((vertex_status_vec_)[*swap_i]).pos_ = swap_i;
					((vertex_status_vec_)[*swap_f]).pos_ = swap_f;

					++swap_i;
					--swap_f;
				}

				// std::cout << "after undo swap " << std::endl;

				// std::cout << *curr_route << std::endl;
				// getchar();
				// getchar();
			}
		}
		++vi_it;
	}
	return false;
}

void HeuristicSolution::AddVertex(int vertex, int route, std::list<int>::iterator it, double profit_variation, double time_variation)
{
	VertexStatus *status = &((vertex_status_vec_)[vertex]);

	// remove from list of unvisited_vertices
	unvisited_vertices_.erase(status->pos_);

	// add to route
	Route *curr_route = &((routes_vec_)[route]);

	status->selected_ = true;
	status->route_ = route;
	status->pos_ = (curr_route->vertices_).insert(it, vertex);

	(curr_route->time_) += time_variation;
	(curr_route->sum_profits_) += profit_variation;

	(profits_sum_) += profit_variation;
}

void HeuristicSolution::InterRouteSwap(int r1, int r2, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1,
									   double profit_variation1, double time_variation1, std::list<int>::iterator it_i2, std::list<int>::iterator it_f2, double profit_variation2, double time_variation2)
{
	Route *route1 = &((routes_vec_)[r1]), *route2 = &((routes_vec_)[r2]);

	// change the routes of the status of the vertices moved!!!!!!!!!!
	std::list<int>::iterator it = it_i1;
	do
	{
		((vertex_status_vec_)[*it]).route_ = r2;
		++it;
	} while (it != it_f1);

	it = it_i2;
	do
	{
		((vertex_status_vec_)[*it]).route_ = r1;
		++it;
	} while (it != it_f2);

	(route1->vertices_).splice(it_i1, route2->vertices_, it_i2, it_f2);
	(route2->vertices_).splice(it_f2, route1->vertices_, it_i1, it_f1);

	(route1->time_) += time_variation1;
	(route1->sum_profits_) += profit_variation1;

	(route2->time_) += time_variation2;
	(route2->sum_profits_) += profit_variation2;

	profits_sum_ += (profit_variation1 + profit_variation2);
}

void HeuristicSolution::InterRouteSwapUnrouted(int r1, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1, double profit_variation1, double time_variation1, std::list<int>::iterator it_i2)
{
	Route *route1 = &((routes_vec_)[r1]);

	// change the routes of the status of the vertices moved!!!!!!!!!!
	std::list<int>::iterator it = it_i1;
	do
	{
		((vertex_status_vec_)[*it]).route_ = -1;
		((vertex_status_vec_)[*it]).selected_ = false;
		//(unvisited_vertices_).push_front(*it);
		//((vertex_status_vec_)[*it]).pos_ = (unvisited_vertices_).begin();
		++it;
	} while (it != it_f1);

	((vertex_status_vec_)[*it_i2]).route_ = r1;
	((vertex_status_vec_)[*it_i2]).selected_ = true;

	(route1->vertices_).splice(it_i1, unvisited_vertices_, it_i2);
	(unvisited_vertices_).splice((unvisited_vertices_).begin(), route1->vertices_, it_i1, it_f1);

	(route1->time_) += time_variation1;
	(route1->sum_profits_) += profit_variation1;

	(profits_sum_) += profit_variation1;

	//(route2->time_) += time_variation2;
	//(route2->sum_profits_) += profit_variation2;
}

bool HeuristicSolution::PreviewInterRouteSwap(Instance &inst, int r1, int pos_i1, int pos_f1, int r2, int pos_i2, int pos_f2, std::list<int>::iterator &it_i1, std::list<int>::iterator &it_f1,
											  double &profit_variation1, double &time_variation1, std::list<int>::iterator &it_i2, std::list<int>::iterator &it_f2, double &profit_variation2, double &time_variation2)
{
	const Graph *graph = inst.graph();
	profit_variation1 = profit_variation2 = 0;
	time_variation1 = time_variation2 = 0.0;
	int pre_vertex1 = 0, pos_vertex1 = 0, pre_vertex2 = 0, pos_vertex2 = 0;
	Route *route1 = &((routes_vec_)[r1]), *route2 = &((routes_vec_)[r2]);
	int max_pos1 = (int)((route1->vertices_).size()) - 1, max_pos2 = (int)((route2->vertices_).size()) - 1;
	GArc *new_pre_arc = nullptr, *new_pos_arc = nullptr;
	int begin_segment_vertex1, end_segment_vertex1, begin_segment_vertex2, end_segment_vertex2;

	int num_vertices = graph->num_vertices();

	int intra_route_profit_loss1 = 0, intra_route_profit_loss2 = 0;
	double intra_route_time_variation1 = 0.0, intra_route_time_variation2 = 0.0;

	if ((pos_i1 > pos_f1) || (pos_i1 > max_pos1) || (pos_f1 > max_pos1))
		return false;
	if ((pos_i2 > pos_f2) || (pos_i2 > max_pos2) || (pos_f2 > max_pos2))
		return false;

	// ROUTE 1
	// compute pre_vertex
	it_i1 = (route1->vertices_).begin();
	if (pos_i1 == 0)
		pre_vertex1 = 0;
	else
	{
		if (pos_i1 > 1)
			std::advance(it_i1, pos_i1 - 1);
		pre_vertex1 = *it_i1;
		++it_i1;
	}

	begin_segment_vertex1 = *it_i1;
	// std::cout << "pre_vertex1: " << pre_vertex1 << std::endl;
	// std::cout << "begin_seg_vertex1: " << begin_segment_vertex1 << std::endl;

	it_f1 = it_i1;
	auto dif = pos_f1 - pos_i1;
	if (dif > 0)
		std::advance(it_f1, dif);

	end_segment_vertex1 = *it_f1;
	++it_f1; // always increment one because the slice function works in the interval [it_i1,it_f1).
	// compute pos_vertex
	if (pos_f1 == max_pos1)
		pos_vertex1 = 0;
	else
		pos_vertex1 = *it_f1;

	// std::cout << "pos_vertex1: " << pos_vertex1 << std::endl;
	// std::cout << "end_seg_vertex1: " << end_segment_vertex1 << std::endl;

	// time_variation1 -= ((*graph)[v1][pos_vertex1])->distance();

	// ROUTE 2
	// compute pre_vertex
	it_i2 = (route2->vertices_).begin();
	if (pos_i2 == 0)
		pre_vertex2 = 0;
	else
	{
		if (pos_i2 > 1)
			std::advance(it_i2, pos_i2 - 1);
		pre_vertex2 = *it_i2;
		++it_i2;
	}

	begin_segment_vertex2 = *it_i2;
	// std::cout << "pre_vertex2: " << pre_vertex2 << std::endl;
	// std::cout << "begin_seg_vertex2: " << begin_segment_vertex2 << std::endl;
	it_f2 = it_i2;

	dif = pos_f2 - pos_i2;
	if (dif > 0)
		std::advance(it_f2, dif);

	end_segment_vertex2 = *it_f2;
	++it_f2; // always increment one because the slice function works in the interval [it_i2,it_f2).
	// compute pos_vertex
	if (pos_f2 == max_pos2)
		pos_vertex2 = 0;
	else
		pos_vertex2 = *it_f2;

	// std::cout << "pos_vertex2: " << pos_vertex2 << std::endl;
	// std::cout << "end_seg_vertex2: " << end_segment_vertex2 << std::endl;

	// time_variation2 -= ((*graph)[v1][pos_vertex2])->distance();

	// profit_variation1 = intra_route_profit_loss2 - intra_route_profit_loss1;

	// route 1
	new_pre_arc = ((*graph)[pre_vertex1][begin_segment_vertex2]);
	new_pos_arc = ((*graph)[end_segment_vertex2][pos_vertex1]);

	if ((new_pre_arc == nullptr && pre_vertex1 != 0 && begin_segment_vertex2 != 0) || (new_pos_arc == nullptr && end_segment_vertex2 != 0 && pos_vertex1 != 0))
		return false;
	// time_variation1 += (new_pre_arc->distance() + intra_route_time_variation2 + new_pos_arc->distance());

	// route 2
	// profit_variation2 = intra_route_profit_loss1 - intra_route_profit_loss2;

	new_pre_arc = ((*graph)[pre_vertex2][begin_segment_vertex1]);
	new_pos_arc = ((*graph)[end_segment_vertex1][pos_vertex2]);

	if ((new_pre_arc == nullptr && pre_vertex2 != 0 && begin_segment_vertex1 != 0) || (new_pos_arc == nullptr && end_segment_vertex1 != 0 && pos_vertex2 != 0))
		return false;

	// SIMULATE SWAP and compute new costs of routes.
	(route1->vertices_).splice(it_i1, route2->vertices_, it_i2, it_f2);
	(route2->vertices_).splice(it_f2, route1->vertices_, it_i1, it_f1);

	auto [new_route_sum_profits1, new_route_max_duration1] = inst.ComputeRouteCosts(route1->vertices_);
	auto [new_route_sum_profits2, new_route_max_duration2] = inst.ComputeRouteCosts(route2->vertices_);

	// auto [new_route_sum_profits12, new_route_max_duration12] = inst.ComputeRouteCosts(route1->vertices_, true);
	// auto [new_route_sum_profits22, new_route_max_duration22] = inst.ComputeRouteCosts(route2->vertices_, true);

	// if (!double_equals(new_route_sum_profits1, new_route_sum_profits12) || !double_equals(new_route_sum_profits2, new_route_sum_profits22) || !double_equals(new_route_max_duration1, new_route_max_duration12) || !double_equals(new_route_max_duration2, new_route_max_duration22))
	// {
	// 	std::cout << new_route_sum_profits1 << " x " << new_route_sum_profits12 << std::endl;
	// 	std::cout << new_route_sum_profits2 << " x " << new_route_sum_profits22 << std::endl;
	// 	std::cout << new_route_max_duration1 << " x " << new_route_max_duration12 << std::endl;
	// 	std::cout << new_route_max_duration2 << " x " << new_route_max_duration22 << std::endl;
	// 	getchar();
	// 	getchar();
	// }

	// UNDO SWAP.
	(route1->vertices_).splice(it_f1, route2->vertices_, it_i1, it_f2);
	(route2->vertices_).splice(it_f2, route1->vertices_, it_i2, it_i1);

	if (double_greater(new_route_max_duration1, inst.limit()) || double_greater(new_route_max_duration2, inst.limit()))
		return false;

	profit_variation1 = new_route_sum_profits1 - route1->sum_profits_;
	time_variation1 = new_route_max_duration1 - route1->time_;

	profit_variation2 = new_route_sum_profits2 - route2->sum_profits_;
	time_variation2 = new_route_max_duration2 - route2->time_;

	return true;
}

bool HeuristicSolution::PreviewInterRouteSwapUnrouted(Instance &inst, int r1, int pos_i1, int pos_f1, std::list<int>::iterator &it_i1, std::list<int>::iterator &it_f1,
													  double &profit_variation1, double &time_variation1, std::list<int>::iterator it_i2)
{
	const Graph *graph = inst.graph();
	profit_variation1 = 0;
	time_variation1 = 0.0;
	int pre_vertex1 = 0, pos_vertex1 = 0;
	Route *route1 = &((routes_vec_)[r1]);
	int max_pos1 = (int)((route1->vertices_).size()) - 1;
	GArc *new_pre_arc = nullptr, *new_pos_arc = nullptr;
	int begin_segment_vertex1, end_segment_vertex1;

	int num_vertices = graph->num_vertices();

	int intra_route_profit_loss1 = 0, intra_route_profit_loss2 = 0;
	double intra_route_time_variation1 = 0.0;

	if ((pos_i1 > pos_f1) || (pos_i1 > max_pos1) || (pos_f1 > max_pos1))
		return false;

	// ROUTE 1
	// compute pre_vertex
	it_i1 = (route1->vertices_).begin();
	if (pos_i1 == 0)
		pre_vertex1 = 0;
	else
	{
		if (pos_i1 > 1)
			std::advance(it_i1, pos_i1 - 1);
		pre_vertex1 = *it_i1;
		++it_i1;
	}

	begin_segment_vertex1 = *it_i1;
	// std::cout << "pre_vertex1: " << pre_vertex1 << std::endl;
	// std::cout << "begin_seg_vertex1: " << begin_segment_vertex1 << std::endl;

	it_f1 = it_i1;
	auto dif = pos_f1 - pos_i1;
	if (dif > 0)
		std::advance(it_f1, dif);

	// compute pos_vertex
	++it_f1; // always increment one because the slice function works in the interval [it_i1,it_f1).
	if (pos_f1 == max_pos1)
		pos_vertex1 = 0;
	else
		pos_vertex1 = *it_f1;

	// route 1
	new_pre_arc = ((*graph)[pre_vertex1][*it_i2]);
	new_pos_arc = ((*graph)[*it_i2][pos_vertex1]);

	if ((new_pre_arc == nullptr && pre_vertex1 != 0 && *it_i2 != 0) || (new_pos_arc == nullptr && *it_i2 != 0 && pos_vertex1 != 0))
		return false;

	auto it_f2 = it_i2;
	++it_f2;

	// simulate swap.
	(route1->vertices_).splice(it_i1, unvisited_vertices_, it_i2);
	(unvisited_vertices_).splice(it_f2, route1->vertices_, it_i1, it_f1);

	auto [new_route_sum_profits1, new_route_max_duration1] = inst.ComputeRouteCosts(route1->vertices_);
	// auto [new_route_sum_profits12, new_route_max_duration12] = inst.ComputeRouteCosts(route1->vertices_, true);

	// if (!double_equals(new_route_sum_profits1, new_route_sum_profits12) || !double_equals(new_route_max_duration1, new_route_max_duration12))
	// {
	// 	std::cout << new_route_sum_profits1 << " x " << new_route_sum_profits12 << std::endl;
	// 	std::cout << new_route_max_duration1 << " x " << new_route_max_duration12 << std::endl;
	// 	getchar();
	// 	getchar();
	// }

	// undo swap.
	(route1->vertices_).splice(it_f1, unvisited_vertices_, it_i1, it_f2);
	(unvisited_vertices_).splice(it_f2, route1->vertices_, it_i2, it_i1);

	if (double_greater(new_route_max_duration1, inst.limit()))
		return false;

	profit_variation1 = new_route_sum_profits1 - route1->sum_profits_;
	time_variation1 = new_route_max_duration1 - route1->time_;

	return true;
}

bool HeuristicSolution::PreviewAddVertexWithinMaximumProfitIncrease(Instance &inst, int vertex, int &route, std::list<int>::iterator &it, double &profit_variation, double &time_variation)
{
	profit_variation = 0;
	time_variation = std::numeric_limits<double>::infinity();
	double curr_profit_variation = 0.0;
	double curr_time_variation = 0.0;
	std::list<int>::iterator curr_it;
	bool can_add = false;

	for (int curr_route = 0; curr_route < num_routes_; ++curr_route)
	{
		if (PreviewAddVertexToRouteWithinMaximumProfitIncrease(inst, vertex, curr_route, curr_it, curr_profit_variation, curr_time_variation))
		{
			if (double_greater(curr_profit_variation, profit_variation) || (double_equals(curr_profit_variation, profit_variation) && double_less(curr_time_variation, time_variation)))
			{
				can_add = true;
				route = curr_route;
				it = curr_it;
				profit_variation = curr_profit_variation;
				time_variation = curr_time_variation;
			}
		}
	}

	return can_add;
}

bool HeuristicSolution::PreviewAddVertexToRouteWithinMaximumProfitIncrease(Instance &inst, int vertex, int route, std::list<int>::iterator &it, double &profit_variation, double &time_variation, bool allow_infeasible_routes)
{
	const Graph *graph = inst.graph();
	profit_variation = 0;
	time_variation = std::numeric_limits<double>::infinity();
	double curr_time_variation = 0.0, curr_profit_variation = 0.0;
	int pre_vertex = 0, pos_vertex = 0;
	Route *curr_route = &((routes_vec_)[route]);
	size_t max_pos = (int)((curr_route->vertices_).size());
	VertexStatus *status = &((vertex_status_vec_)[vertex]);
	int num_vertices = graph->num_vertices();
	bool can_add = false;

	if (status->selected_)
		return false;

	std::list<int>::iterator curr_it = (curr_route->vertices_).begin();
	for (int pos = 0; pos <= max_pos; ++pos)
	{
		// compute pos_vertex
		if (pos == max_pos)
			pos_vertex = 0;
		else
			pos_vertex = *curr_it;

		GArc *pre_arc = (*graph)[pre_vertex][vertex], *pos_arc = (*graph)[vertex][pos_vertex], *curr_arc = (*graph)[pre_vertex][pos_vertex];

		if ((pre_arc != nullptr) && (pos_arc != nullptr) && (curr_arc != nullptr || (pre_vertex == 0 && pos_vertex == 0)))
		{
			// simulate addition of vertex.

			auto inserted_vertex_it = curr_route->vertices_.insert(curr_it, vertex);
			auto [new_route_sum_profits, new_route_max_duration] = inst.ComputeRouteCosts(curr_route->vertices_);

			// auto [new_route_sum_profits12, new_route_max_duration12] = inst.ComputeRouteCosts(curr_route->vertices_, true);

			// if (!double_equals(new_route_sum_profits, new_route_sum_profits12) || !double_equals(new_route_max_duration, new_route_max_duration12))
			// {
			// 	std::cout << new_route_sum_profits << " x " << new_route_sum_profits12 << std::endl;
			// 	std::cout << new_route_max_duration << " x " << new_route_max_duration12 << std::endl;
			// 	getchar();
			// 	getchar();
			// }

			// remove back the vertex added.
			curr_route->vertices_.erase(inserted_vertex_it);

			if (!double_greater(new_route_max_duration, inst.limit()) || allow_infeasible_routes)
			{
				curr_time_variation = new_route_max_duration - curr_route->time_;
				curr_profit_variation = new_route_sum_profits - curr_route->sum_profits_;

				if (double_greater(curr_profit_variation, profit_variation) || (double_equals(curr_profit_variation, profit_variation) && double_less(curr_time_variation, time_variation)))
				{
					can_add = true;
					time_variation = curr_time_variation;
					profit_variation = new_route_sum_profits - curr_route->sum_profits_;
					it = curr_it;
				}
			}
		}
		++curr_it;
		pre_vertex = pos_vertex;
	}

	return can_add;
}

bool HeuristicSolution::PreviewAddVertex(Instance &inst, int vertex, int route, int pos, std::list<int>::iterator &it, double &profit_variation, double &time_variation)
{
	const Graph *graph = inst.graph();
	profit_variation = 0;
	time_variation = 0.0;
	int pre_vertex = 0, pos_vertex = 0;
	Route *curr_route = &((routes_vec_)[route]);
	int max_pos = (int)((curr_route->vertices_).size());
	VertexStatus *status = &((vertex_status_vec_)[vertex]);
	int num_vertices = graph->num_vertices();

	if ((pos == -1) || (pos > max_pos))
		pos = max_pos;

	// compute pre_vertex
	it = (curr_route->vertices_).begin();
	if (pos == 0)
		pre_vertex = 0;
	else
	{
		if (pos > 1)
			std::advance(it, pos - 1);
		pre_vertex = *it;
		++it;
	}

	// compute pos_vertex
	if (pos == max_pos)
		pos_vertex = 0;
	else
		pos_vertex = *it;

	GArc *pre_arc = (*graph)[pre_vertex][vertex], *pos_arc = (*graph)[vertex][pos_vertex], *curr_arc = (*graph)[pre_vertex][pos_vertex];

	if ((!(status->selected_)) && (pre_arc != nullptr) && (pos_arc != nullptr) && ((curr_arc != nullptr) || (pre_vertex == 0 && pos_vertex == 0)))
	{
		// simulate addition of vertex to route.
		it = (curr_route->vertices_).begin();
		if (pos > 0)
			std::advance(it, pos);

		auto inserted_vertex_it = curr_route->vertices_.insert(it, vertex);
		auto [new_route_sum_profits, new_route_max_duration] = inst.ComputeRouteCosts(curr_route->vertices_);

		// auto [new_route_sum_profits12, new_route_max_duration12] = inst.ComputeRouteCosts(curr_route->vertices_, true);

		// if (!double_equals(new_route_sum_profits, new_route_sum_profits12) || !double_equals(new_route_max_duration, new_route_max_duration12))
		// {
		// 	std::cout << new_route_sum_profits << " x " << new_route_sum_profits12 << std::endl;
		// 	std::cout << new_route_max_duration << " x " << new_route_max_duration12 << std::endl;
		// 	getchar();
		// 	getchar();
		// }

		profit_variation = new_route_sum_profits - curr_route->sum_profits_;
		time_variation = new_route_max_duration - curr_route->time_;

		// remove back the vertex added.
		curr_route->vertices_.erase(inserted_vertex_it);

		if (!double_greater(new_route_max_duration, inst.limit()))
			return true;
	}
	return false;
}

void HeuristicSolution::RemoveVertex(int vertex, double profit_variation, double time_variation)
{
	VertexStatus *status = &((vertex_status_vec_)[vertex]);
	// std::cout << "removing from route " << status->route_ << std::endl;
	//  remove from route
	Route *curr_route = &((routes_vec_)[status->route_]);
	// std::cout << *curr_route << std::endl;
	(curr_route->vertices_).erase(status->pos_);
	(curr_route->time_) += time_variation;
	(curr_route->sum_profits_) += profit_variation;

	// add to list of unvisited_vertices_
	(unvisited_vertices_).push_front(vertex);
	status->selected_ = false;
	status->route_ = -1;
	status->pos_ = (unvisited_vertices_).begin();

	(profits_sum_) += profit_variation;
}

bool HeuristicSolution::PreviewRemoveVertex(Instance &inst, int vertex, double &profit_variation, double &time_variation, bool allow_infeasible_routes)
{
	const Graph *graph = inst.graph();
	profit_variation = 0;
	time_variation = 0.0;
	VertexStatus *status = &((vertex_status_vec_)[vertex]);
	if (!(status->selected_))
		return false;

	int num_vertices = graph->num_vertices();
	int pre_vertex = 0, pos_vertex = 0;
	Route *curr_route = &((routes_vec_)[status->route_]);
	std::list<int>::iterator pos = status->pos_;
	int max_pos = (int)((curr_route->vertices_).size());

	if (max_pos == 0)
		return false;

	std::list<int>::iterator pre_vertex_it = pos;

	if (pos == (curr_route->vertices_).begin())
		pre_vertex = 0;
	else
	{
		--pre_vertex_it;
		pre_vertex = *pre_vertex_it;
	}
	// std::cout << "pre_vertex: " << pre_vertex << std::endl;

	std::list<int>::iterator pos_vertex_it = pos; // at this point, pos is never the .end(), cause size > 0!!!
	++pos_vertex_it;
	if (pos_vertex_it == (curr_route->vertices_).end())
		pos_vertex = 0;
	else
		pos_vertex = *pos_vertex_it;

	// std::cout << "pos_vertex: " << pos_vertex << std::endl;

	GArc *new_arc = (*graph)[pre_vertex][pos_vertex];

	// only allow if there is a valid connection between pre_vertex and pos_vertex.
	if ((new_arc != nullptr) || ((pre_vertex == 0) && (pos_vertex == 0)))
	{
		// simulate removal from route.
		(curr_route->vertices_).erase(status->pos_);

		auto [new_route_sum_profits, new_route_max_duration] = inst.ComputeRouteCosts(curr_route->vertices_);

		// auto [new_route_sum_profits12, new_route_max_duration12] = inst.ComputeRouteCosts(curr_route->vertices_, true);

		// if (!double_equals(new_route_sum_profits, new_route_sum_profits12) || !double_equals(new_route_max_duration, new_route_max_duration12))
		// {
		// 	std::cout << new_route_sum_profits << " x " << new_route_sum_profits12 << std::endl;
		// 	std::cout << new_route_max_duration << " x " << new_route_max_duration12 << std::endl;
		// 	getchar();
		// 	getchar();
		// }

		profit_variation = new_route_sum_profits - curr_route->sum_profits_;
		time_variation = new_route_max_duration - curr_route->time_;

		// std::cout << "calculou" << std::endl;

		// add vertex back.
		// if (pre_vertex == 0) // if the first vertex in route, update pre_vertex_it to the beginning of the list again, since, when the first element was removed, the previous iterator was also destroyed.
		// 	pre_vertex_it = (curr_route->vertices_).begin();
		status->pos_ = (curr_route->vertices_).insert(pos_vertex_it, vertex);
		// std::cout << "adicionou de volta" << std::endl;
		// std::cout << *(status->pos_) << std::endl;

		if (!double_greater(new_route_max_duration, inst.limit()) || allow_infeasible_routes)
			return true;
	}
	return false;
}

bool HeuristicSolution::PreviewInterRouteMoveVertex(Instance &inst, int vertex, int r2, int pos, std::list<int>::iterator &it, double &profit_variation1, double &time_variation1, double &profit_variation2, double &time_variation2)
{
	const Graph *graph = inst.graph();
	bool can_move = false;
	profit_variation1 = profit_variation2 = 0;
	time_variation1 = time_variation2 = 0.0;
	VertexStatus *status = &((vertex_status_vec_)[vertex]);
	if (status->route_ != r2)
	{
		if (PreviewRemoveVertex(inst, vertex, profit_variation1, time_variation1))
		{
			// to do the preview of addition!!!
			status->selected_ = false;
			if (PreviewAddVertex(inst, vertex, r2, pos, it, profit_variation2, time_variation2))
			{
				can_move = true;
			}
			status->selected_ = true;
		}
	}
	return can_move;
}

void HeuristicSolution::InterRouteMoveVertex(int vertex, int r2, std::list<int>::iterator it, double profit_variation1, double time_variation1, double profit_variation2, double time_variation2)
{
	VertexStatus *status = &((vertex_status_vec_)[vertex]);
	// remove from route r1
	Route *curr_route = &((routes_vec_)[status->route_]);
	(curr_route->vertices_).erase(status->pos_);
	(curr_route->time_) += time_variation1;
	(curr_route->sum_profits_) += profit_variation1;

	// add to route r2
	curr_route = &((routes_vec_)[r2]);
	int max_pos = (int)((curr_route->vertices_).size());
	// std::list<int>::iterator it = (curr_route->vertices_).begin();
	// if((pos == -1)||(pos > max_pos)) pos = max_pos;
	// if(pos > 0) std::advance(it,pos);

	status->route_ = r2;
	status->pos_ = (curr_route->vertices_).insert(it, vertex);

	(curr_route->time_) += time_variation2;
	(curr_route->sum_profits_) += profit_variation2;

	profits_sum_ += (profit_variation1 + profit_variation2);
}

double HeuristicSolution::ComputeSolutionCost(Instance &instance) const
{
	double total_cost = 0.0;
	for (auto route : routes_vec_)
	{
		auto [curr_cost, _] = instance.ComputeRouteCosts(route.vertices_);
		total_cost += curr_cost;
	}
	return total_cost;
}
double HeuristicSolution::ComputeSolutionCostRec(Instance &instance, bool memoization) const
{
	double total_cost = 0.0;
	for (auto route : routes_vec_)
	{
		auto [curr_cost, _] = instance.ComputeRouteCostsRec(route.vertices_, memoization);
		total_cost += curr_cost;
	}
	return total_cost;
}

bool HeuristicSolution::CheckCorrectness(const Instance &instance)
{
	const Graph *graph = instance.graph();
	int num_vertices = graph->num_vertices(), num_mandatory = instance.num_mandatory(), num_vehicles = instance.num_vehicles(), num_arcs = graph->num_arcs();
	Route *curr_route = nullptr;
	VertexStatus *curr_status = nullptr;
	double total_profits = 0.0;
	GArc *curr_arc = nullptr;
	int v1 = 0, v2 = 0;
	int count_mandatory = 0;
	boost::dynamic_bitset<> visited_vertices(num_vertices, 0);
	boost::dynamic_bitset<> visited_arcs(num_arcs, 0);

	visited_vertices[0] = 1;

	if (num_routes_ != num_vehicles)
		throw "Invalid solution 1";

	if (is_infeasible_)
		return true;

	for (int i = 0; i < num_vehicles; i++)
	{
		curr_route = &((routes_vec_)[i]);
		// std::cout << *curr_route << std::endl;
		v1 = 0;
		if (!((curr_route->vertices_).empty()))
		{
			for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			{
				v2 = *it;
				if (visited_vertices[v2])
					throw "Invalid solution 2";
				visited_vertices[v2] = 1;

				if ((v2 > 0) && (v2 <= num_mandatory))
					++count_mandatory;
				curr_status = &((vertex_status_vec_)[v2]);
				if ((!(curr_status->selected_)) || (it != curr_status->pos_) || (curr_status->route_ != i))
				{
					std::cout << curr_status->selected_ << std::endl;
					std::cout << *it << " x " << *(curr_status->pos_) << std::endl;
					std::cout << i << " x " << curr_status->route_ << std::endl;
					throw "Invalid solution 3";
				}
				// if (v2 > num_mandatory)
				// {
				// 	total_profits += (graph->vertices_info())[v2].profit_;
				// 	curr_profits += (graph->vertices_info())[v2].profit_;
				// }

				curr_arc = (*graph)[v1][v2];
				if (curr_arc == nullptr)
					throw "Invalid solution 4";
				// curr_time += (curr_arc->distance());
				visited_arcs[graph->pos(v1, v2)] = 1;
				v1 = v2;
			}

			v2 = 0;
			curr_arc = (*graph)[v1][v2];
			if (curr_arc == nullptr)
				throw "Invalid solution 5";
			// curr_time += (curr_arc->distance());
			visited_arcs[graph->pos(v1, v2)] = 1;
		}

		auto [route_sum_profits, route_max_duration] = instance.ComputeRouteCosts(curr_route->vertices_);

		if (!double_equals(route_max_duration, curr_route->time_))
		{
			std::cout << route_max_duration << " " << curr_route->time_ << std::endl;
			throw "Invalid solution 7";
		}
		if (double_greater(curr_route->time_, instance.limit()))
			throw "Invalid solution 7.2";
		if (!double_equals(route_sum_profits, curr_route->sum_profits_))
		{
			std::cout << route_sum_profits << " x " << curr_route->sum_profits_ << std::endl;
			throw "Invalid solution 8";
		}

		total_profits += route_sum_profits;
	}

	if (!double_equals(total_profits, profits_sum_))
		throw "Invalid solution 9";

	if (visited_vertices.count() != ((size_t)num_vertices - (unvisited_vertices_).size()))
		throw "Invalid solution 10";

	for (auto it = (unvisited_vertices_).begin(); it != (unvisited_vertices_).end(); ++it)
	{
		curr_status = &((vertex_status_vec_)[*it]);
		if ((curr_status->selected_) || (curr_status->pos_ != it) || (curr_status->route_ != -1) || visited_vertices[*it])
			throw "Invalid solution 11";
	}

	if (count_mandatory != num_mandatory)
		throw "invalid solution 12";

	if (visited_arcs != bitset_arcs_)
		throw "Invalid solution 13";
	if (visited_vertices != bitset_vertices_)
		throw "Invalid solution 14";
	return true;
}

void HeuristicSolution::BuildBitset(const Instance &instance)
{
	bitset_arcs_ = boost::dynamic_bitset<>(num_arcs_, 0);
	bitset_vertices_ = boost::dynamic_bitset<>(num_vertices_, 0);
	const Graph *graph = instance.graph();
	int num_vertices = graph->num_vertices(), num_vehicles = instance.num_vehicles();
	Route *curr_route = nullptr;
	int v1 = 0, v2 = 0;

	if (is_infeasible_)
		return;

	(bitset_vertices_)[0] = 1;

	for (int i = 0; i < num_vehicles; i++)
	{
		v1 = 0;
		curr_route = &((routes_vec_)[i]);
		if (!((curr_route->vertices_).empty()))
		{
			for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			{
				v2 = *it;
				(bitset_vertices_)[v2] = 1;

				(bitset_arcs_)[graph->pos(v1, v2)] = 1;
				v1 = v2;
			}

			v2 = 0;
			(bitset_arcs_)[graph->pos(v1, v2)] = 1;
		}
	}
}

bool HeuristicSolution::operator==(HeuristicSolution &other)
{
	return (bitset_arcs_ == other.bitset_arcs_);
}

std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol)
{
	if (sol.is_infeasible_)
		out << "[INFEASIBLE]";
	if (sol.is_optimal_)
		out << "[OPTMAL]";
	out << "[profits_sum: " << sol.profits_sum_ << "]" << std::endl;
	for (int i = 0; i < sol.num_routes_; ++i)
		out << ((sol.routes_vec_)[i]);
	out << "unvisited: ";
	for (auto it = (sol.unvisited_vertices_).begin(); it != (sol.unvisited_vertices_).end(); ++it)
		out << *it << " ";
	out << std::endl;

	out << "vertex status: " << std::endl;
	for (int vertex = 1; vertex < sol.vertex_status_vec_.size(); ++vertex) // skip the zero, cause it's not supposed to be in any list.
	{
		auto status = &(sol.vertex_status_vec_[vertex]);
		out << vertex << " " << status->selected_ << " " << status->route_ << " " << *status->pos_ << std::endl;
	}
	return out;
}

void HeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const
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

	file.open(path.c_str(), std::fstream::out | std::fstream::app);

	file << std::setprecision(5) << std::fixed;

	file << K_FILE_DELIMITER << std::endl
		 << "Total profit sum: " << profits_sum_ << std::endl
		 << "Num routes: ";
	if ((is_infeasible_) || (!(is_feasible_)))
		file << "0" << std::endl;
	else
	{
		file << num_routes_ << std::endl;

		for (int i = 0; i < num_routes_; i++)
		{
			const Route &curr_route = (routes_vec_)[i];
			file << curr_route.sum_profits_ << " " << curr_route.time_ << " ";

			for (auto it = (curr_route.vertices_).begin(); it != (curr_route.vertices_).end(); ++it)
			{
				file << instance.getOriginalVertexPosition(*it) << " ";
			}

			// file << num_vertices_ - 1;
			file << std::endl;
		}
	}

	file.close();
}

void HeuristicSolution::ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	const Graph *graph = inst.graph();
	int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();

	Reset(num_vertices, num_arcs, inst.num_vehicles());
	Route *curr_route = nullptr;
	int num_routes = 0;
	std::fstream input;
	std::string path = "..//solutions//";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	input.open(path.c_str(), std::fstream::in);
	std::string line;

	if (!(input.is_open()))
		throw "Error opening file 1";
	while (!(input.eof()))
	{
		getline(input, line);
		if (line == K_FILE_DELIMITER)
		{
			std::stringstream s_profits_sum, s_num_routes;

			getline(input, line);
			size_t pos = line.find_first_of(":");
			s_profits_sum << line.substr(pos + 2);
			s_profits_sum >> profits_sum_;
			// std::cout << profits_sum_ << std::endl;

			getline(input, line);
			pos = line.find_first_of(":");
			s_num_routes << line.substr(pos + 2);
			s_num_routes >> num_routes;

			if ((num_routes == 0) || (num_routes != num_routes_))
				is_feasible_ = false;
			else
			{
				is_feasible_ = true;
				int pre_vertex = 0, curr_vertex = 0, vertex_pos_before_reordering = 0;
				VertexStatus *status = nullptr;
				for (int i = 0; i < num_routes; i++)
				{
					curr_route = &((routes_vec_)[i]);
					getline(input, line);
					std::stringstream s_route(line);
					// std::cout << s_route.str() << std::endl;
					s_route >> curr_route->sum_profits_ >> curr_route->time_;
					// std::cout << curr_route->sum_profits_ << " " << curr_route->time_ << std::endl;
					curr_vertex = pre_vertex = 0;
					while (s_route >> vertex_pos_before_reordering)
					{
						// in the file/solution, is saved with original positions of vertices (before reordering!)
						curr_vertex = inst.getReorderedVertexPosition(vertex_pos_before_reordering);
						status = &((vertex_status_vec_)[curr_vertex]);

						// remove from list of unvisited_vertices
						(unvisited_vertices_).erase(status->pos_);

						// adds vertex to route
						status->selected_ = true;
						status->route_ = i;
						status->pos_ = (curr_route->vertices_).insert((curr_route->vertices_).end(), curr_vertex);

						//(bitset_arcs_)[graph->pos(pre_vertex,curr_vertex)] = 1;
						pre_vertex = curr_vertex;
					}

					curr_vertex = num_vertices - 1;
					//(bitset_arcs_)[graph->pos(pre_vertex,curr_vertex)] = 1;
				}
			}
			break;
		}
	}

	if (!(is_infeasible_) && (is_feasible_) && (unvisited_vertices_).empty())
		is_optimal_ = true;

	input.close();
}

MetaHeuristicSolution::MetaHeuristicSolution() : HeuristicSolution()
{
}

MetaHeuristicSolution::MetaHeuristicSolution(int num_vertices, int num_arcs, int num_routes) : HeuristicSolution(num_vertices, num_arcs, num_routes)
{
}

MetaHeuristicSolution::MetaHeuristicSolution(const HeuristicSolution *sol) : HeuristicSolution(sol->num_vertices_, sol->num_arcs_, sol->num_routes_)
{
	is_infeasible_ = sol->is_infeasible_;
	is_feasible_ = sol->is_feasible_;
	is_optimal_ = sol->is_optimal_;
	profits_sum_ = sol->profits_sum_;
	total_time_spent_ = sol->total_time_spent_;
	routes_vec_ = sol->routes_vec_;
	vertex_status_vec_ = sol->vertex_status_vec_;
	unvisited_vertices_ = sol->unvisited_vertices_;
	VertexStatus *status = nullptr;
	Route *curr_route = nullptr;

	for (int i = 0; i < num_routes_; ++i)
	{
		curr_route = &((routes_vec_)[i]);
		for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			((vertex_status_vec_)[*it]).pos_ = it;
	}

	for (auto it = (unvisited_vertices_).begin(); it != (unvisited_vertices_).end(); ++it)
		((vertex_status_vec_)[*it]).pos_ = it;
}

MetaHeuristicSolution::MetaHeuristicSolution(const MetaHeuristicSolution *sol)
{
	(*this) = (*sol);

	VertexStatus *status = nullptr;
	Route *curr_route = nullptr;

	for (int i = 0; i < num_routes_; ++i)
	{
		curr_route = &((routes_vec_)[i]);
		for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
			((vertex_status_vec_)[*it]).pos_ = it;
	}

	for (auto it = (unvisited_vertices_).begin(); it != (unvisited_vertices_).end(); ++it)
		((vertex_status_vec_)[*it]).pos_ = it;
}

MetaHeuristicSolution::~MetaHeuristicSolution()
{
}

void MetaHeuristicSolution::Reset(int num_vertices, int num_arcs, int num_routes)
{
	HeuristicSolution::Reset(num_vertices, num_arcs, num_routes);

	initial_solution_profits_sum_ = 0;
	num_iterations_ = 0;
	last_improve_iteration_ = 0;
}

void MetaHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const
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
	else if (is_optimal_)
		file << "STATUS: OPTIMAL" << std::endl;
	else if (is_feasible_)
		file << "STATUS: FEASIBLE" << std::endl;
	else
		file << "STATUS: POSSIBLY FEASIBLE" << std::endl;
	file << "profits sum: " << profits_sum_ << std::endl;
	file << "initial profits sum: " << initial_solution_profits_sum_ << std::endl;
	file << "# iterations: " << num_iterations_ << std::endl;
	file << "time(s): " << total_time_spent_ << std::endl;
	file << "last iteration of improvement: " << last_improve_iteration_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
}

void MetaHeuristicSolution::ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	HeuristicSolution::ReadFromFile(inst, algo, folder, file_name);

	std::fstream input;
	std::string path = "..//solutions//";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	input.open(path.c_str(), std::fstream::in);

	if (!(input.is_open()))
		throw "Error opening file 2";
	std::string line;
	getline(input, line);

	if (line == "STATUS: FOUND ARC INTEGER FEASIBLE")
		is_feasible_ = true;
	if (line == "STATUS: INFEASIBLE")
	{
		is_infeasible_ = true;
		is_feasible_ = false;
	}

	std::stringstream s_time_phase1, s_time_phase2;
	double time_phase1 = 0.0, time_phase2 = 0.0;
	for (int i = 1; i <= 6; ++i)
		getline(input, line);
	size_t pos = line.find_first_of(":");
	s_time_phase1 << line.substr(pos + 2);
	s_time_phase1 >> time_phase1;

	for (int i = 1; i <= 5; ++i)
		getline(input, line);
	pos = line.find_first_of(":");
	s_time_phase2 << line.substr(pos + 2);
	s_time_phase2 >> time_phase2;

	total_time_spent_ = time_phase1 + time_phase2;

	input.close();
}

std::string MetaHeuristicSolution::GenerateALNSFileName(int num_iterations, int pool_size, int multithreading)
{
	std::string file_name = "alns";
	if (multithreading)
		file_name += "_MT";
	file_name += "_i_" + std::to_string(num_iterations);
	file_name += "_ps_" + std::to_string(pool_size);

	return file_name;
}

std::string MetaHeuristicSolution::GenerateSimulatedAnnealingFileName(double temp_decay_rate, int multithreading)
{
	std::stringstream ss_decay_factor;
	ss_decay_factor << std::fixed << std::setprecision(2) << temp_decay_rate;
	std::string file_name = "sa";
	if (multithreading)
		file_name += "_MT";
	file_name += "_df_" + ss_decay_factor.str();

	return file_name;
}

FPHeuristicSolution::FPHeuristicSolution(int num_vertices, int num_arcs, int num_routes) : HeuristicSolution(num_vertices, num_arcs, num_routes)
{
	num_iterations_stage1_ = 0;
	num_iterations_stage2_ = 0;
	num_perturbations_stage1_ = 0;
	num_perturbations_stage2_ = 0;
	num_restarts_stage1_ = 0;
	num_restarts_stage2_ = 0;
	found_y_integer_ = false;
	found_x_integer_ = false;
	time_stage1_ = 0.0;
	time_stage2_ = 0.0;
}

FPHeuristicSolution::FPHeuristicSolution() : HeuristicSolution()
{
	num_iterations_stage1_ = 0;
	num_iterations_stage2_ = 0;
	num_perturbations_stage1_ = 0;
	num_perturbations_stage2_ = 0;
	num_restarts_stage1_ = 0;
	num_restarts_stage2_ = 0;
	found_y_integer_ = false;
	found_x_integer_ = false;
	time_stage1_ = 0.0;
	time_stage2_ = 0.0;
}

FPHeuristicSolution::~FPHeuristicSolution()
{
}

void FPHeuristicSolution::Reset(int num_vertices, int num_arcs, int num_routes)
{
	HeuristicSolution::Reset(num_vertices, num_arcs, num_routes);

	num_iterations_stage1_ = 0;
	num_iterations_stage2_ = 0;
	num_perturbations_stage1_ = 0;
	num_perturbations_stage2_ = 0;
	num_restarts_stage1_ = 0;
	num_restarts_stage2_ = 0;
	found_y_integer_ = false;
	found_x_integer_ = false;
	time_stage1_ = 0.0;
	time_stage2_ = 0.0;
}

std::string FPHeuristicSolution::GenerateFileName()
{
	std::string file_name;
	if (K_STOP)
		file_name += "stop_";
	file_name += "fp";
	if (K_FEASIBILITY_PUMP_ADD_CUTS)
		file_name += "_cuts";
	if (K_SOLVE_STAGE_1)
		file_name += "_alpha1_" + std::to_string(int(100 * initial_alpha_stage1));
	file_name += ("_alpha2_" + std::to_string(int(100 * initial_alpha_stage2)) + "_p_");
	if (FIXED_FLIP_BASIS)
		file_name += ("fixed_" + std::to_string(K_FLIP_BASIS));
	else
		file_name += std::to_string(int(100 * perturbation_flip_percentage));
	return file_name;
}

void FPHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const
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
	else if (found_y_integer_)
		file << "STATUS: FOUND VERTEX INTEGER FEASIBLE" << std::endl;
	else
		file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;
	file << "profit sum: " << profits_sum_ << std::endl;
	file << "STAGE 1: " << std::endl;
	file << "# iterations: " << num_iterations_stage1_ << std::endl;
	file << "# perturbations: " << num_perturbations_stage1_ << std::endl;
	file << "# restarts: " << num_restarts_stage1_ << std::endl;
	file << "time(s): " << time_stage1_ << std::endl;
	file << "STAGE 2: " << std::endl;
	file << "# iterations: " << num_iterations_stage2_ << std::endl;
	file << "# perturbations: " << num_perturbations_stage2_ << std::endl;
	file << "# restarts: " << num_restarts_stage2_ << std::endl;
	file << "time(s): " << time_stage2_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
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

LBHeuristicSolution::LBHeuristicSolution(int num_vertices, int num_arcs, int num_routes) : HeuristicSolution(num_vertices, num_arcs, num_routes)
{
	num_iterations_ = 0;
}

LBHeuristicSolution::LBHeuristicSolution() : HeuristicSolution()
{
	num_iterations_ = 0;
}

LBHeuristicSolution::~LBHeuristicSolution()
{
}

void LBHeuristicSolution::Reset(int num_vertices, int num_arcs, int num_routes)
{
	HeuristicSolution::Reset(num_vertices, num_arcs, num_routes);

	num_iterations_ = 0;
}

std::string LBHeuristicSolution::GenerateFileName()
{
	/*
		std::string file_name;
		if(K_STOP) file_name += "stop_";
		file_name += "fp";
		if(K_FEASIBILITY_PUMP_ADD_CUTS) file_name += "_cuts";
		if(K_SOLVE_STAGE_1) file_name += "_alpha1_" + std::to_string(int(100*initial_alpha_stage1));
		file_name += ("_alpha2_" + std::to_string(int(100*initial_alpha_stage2)) + "_p_");
		if(FIXED_FLIP_BASIS) file_name += ("fixed_" + std::to_string(K_FLIP_BASIS));
		else file_name += std::to_string(int(100*perturbation_flip_percentage));
		return file_name;*/
	return "";
}

void LBHeuristicSolution::WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const
{ /*
	 std::fstream file;
	 std::string path = ".//solutions";
	 path.append(folder);
	 //struct stat sb;
	 //if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	 path.append("s_");
	 path.append(algo);
	 path.append("_");
	 path.append(file_name);
	 //std::cout << path << std::endl;

	 file.open(path.c_str(),std::fstream::out);

	 file << std::setprecision(5) << std::fixed;

		 if(is_infeasible_) file << "STATUS: INFEASIBLE" << std::endl;
	 else if(found_x_integer_) file << "STATUS: FOUND ARC INTEGER FEASIBLE" << std::endl;
		 else if(found_y_integer_) file << "STATUS: FOUND VERTEX INTEGER FEASIBLE" << std::endl;
		 else file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;
		 file << "profit sum: " << profits_sum_ << std::endl;
	 file << "STAGE 1: " << std::endl;
	 file << "# iterations: " << num_iterations_stage1_ << std::endl;
	 file << "# perturbations: " << num_perturbations_stage1_ << std::endl;
		 file << "# restarts: " << num_restarts_stage1_ << std::endl;
	 file << "time(s): " << time_stage1_ << std::endl;
		 file << "STAGE 2: " << std::endl;
	 file << "# iterations: " << num_iterations_stage2_ << std::endl;
	 file << "# perturbations: " << num_perturbations_stage2_ << std::endl;
		 file << "# restarts: " << num_restarts_stage2_ << std::endl;
	 file << "time(s): " << time_stage2_ << std::endl;

		 file.close();
 */
	HeuristicSolution::WriteToFile(instance, algo, folder, file_name);
}
