#include <algorithm>
#include <queue>
#include "src/graph_algorithms.h"
#include "src/general.h"

#define pp std::pair<int, double>

struct pri
{
	int operator()(const std::pair<int, double> &p1, const std::pair<int, double> &p2)
	{
		return double_less(p1.second, p2.second);
	}
} p;

double *Dijkstra(const Graph *graph, bool inverse, bool unit_dist)
{
	std::priority_queue<pp, std::vector<pp>, pri> q;
	int s = 0, u = 0, v = 0;
	double w = 0.0, vertex_nominal_time = 0.0, new_dist = 0.0;
	GArc *curr_arc = NULL;
	int dimension = graph->num_vertices();
	double *d = new double[dimension];
	const auto *vertices_info = graph->vertices_info();

	// if(inverse) s = dimension - 1; // because, here, source == destination.

	// init
	for (int i = 0; i < dimension; ++i)
		d[i] = std::numeric_limits<double>::infinity();

	d[s] = 0;
	q.push(pp(s, d[s]));

	while (!q.empty())
	{
		u = (q.top()).first;
		q.pop();

		if (inverse)
			v = u;

		for (int i = 0; i < dimension; ++i)
		{
			if (inverse)
				u = i;
			else
				v = i;

			curr_arc = (*graph)[u][v];
			if (curr_arc)
			{
				if (unit_dist)
					w = 1;
				else
					w = curr_arc->distance();
				if (inverse)
				{
					vertex_nominal_time = v == 0 ? 0.0 : vertices_info[v].nominal_service_time_;
					new_dist = d[v] + w + vertex_nominal_time;
					if (double_greater(d[u], new_dist))
					{
						d[u] = new_dist;
						q.push(pp(u, d[u]));
					}
				}
				else
				{
					vertex_nominal_time = u == 0 ? 0.0 : vertices_info[u].nominal_service_time_;
					new_dist = d[u] + w + vertex_nominal_time;
					if (double_greater(d[v], new_dist))
					{
						d[v] = new_dist;
						q.push(pp(v, d[v]));
					}
				}
			}
		}
	}

	return d;
}

void BronKerosch(std::vector<std::list<int>> *adj_lists, std::list<int> R, std::unordered_set<int> P, std::unordered_set<int> X, std::vector<std::list<int>> &cliques)
{
	if ((P.empty()) && (X.empty()) && (R.size() >= 1))
		cliques.push_back(R);
	while (!(P.empty()))
	{
		auto it = P.begin();
		std::unordered_set<int> P2, X2;
		std::list<int> R2 = R;
		R2.push_back(*it);
		for (auto it2 = ((*adj_lists)[*it]).begin(); it2 != ((*adj_lists)[*it]).end(); ++it2)
		{
			if (P.find(*it2) != P.end())
				P2.insert(*it2);
			if (X.find(*it2) != X.end())
				X2.insert(*it2);
		}

		BronKerosch(adj_lists, R2, P2, X2, cliques);
		X.insert(*it);
		P.erase(it);
	}
}

void Tomita(std::vector<boost::dynamic_bitset<>> *adj_lists, std::list<int> R, boost::dynamic_bitset<> P, boost::dynamic_bitset<> X, std::vector<std::list<int>> &cliques)
{
	if ((P.none()) && (X.none()) && (R.size() >= 1))
		cliques.push_back(R);

	if (P.none())
		return;
	boost::dynamic_bitset<> PUX = P | X;
	// std::cout << "P: " << P << std::endl;
	// std::cout << "X: " << X << std::endl;
	// std::cout << "PUX: " << PUX << std::endl;
	// getchar();getchar();

	// selects pivot
	int pivot = -1, max_num_neighbours = -1, curr_num_neighbours = 0;
	for (size_t u = 0; u < adj_lists->size(); ++u)
	{
		if (PUX[u] == 1)
		{
			curr_num_neighbours = (P & ((*adj_lists)[u])).count();
			if (curr_num_neighbours > max_num_neighbours)
			{
				pivot = u;
				max_num_neighbours = curr_num_neighbours;
			}
		}
	}

	// std::cout << "u: " << u << std::endl;
	boost::dynamic_bitset<> aux = P - ((*adj_lists)[pivot]);
	// std::cout << "P-N(u): " << aux << std::endl;
	for (size_t v = 0; v < adj_lists->size(); ++v)
	{
		if (aux[v] == 1)
		{
			aux[v] = 0;
			std::list<int> R2 = R;
			R2.push_back(v);

			Tomita(adj_lists, R2, P & ((*adj_lists)[v]), X & ((*adj_lists)[v]), cliques);
			X[v] = 1;
			P[v] = 0;
		}
	}
}

Matrix<double> *FloydWarshall(const Graph *graph)
{
	if (!graph)
		return nullptr;

	int dimension = graph->num_vertices();
	Matrix<double> *d = new Matrix<double>(dimension, dimension, std::numeric_limits<double>::infinity());
	const auto *vertices_info = graph->vertices_info();
	double vertex_nominal_time = 0.0, new_dist = 0.0;

	for (int k = 0; k < dimension; k++)
	{
		(*d)[k][k] = 0.0;
		for (const auto adj_vertex : graph->AdjVerticesOut(k))
			(*d)[k][adj_vertex] = ((*graph)[k][adj_vertex])->distance();
	}

	for (int k = 0; k < dimension; ++k)
	{
		vertex_nominal_time = k == 0 ? 0 : vertices_info[k].nominal_service_time_;
		// std::cout << "* "<<  k << std::endl;
		for (int i = 0; i < dimension; ++i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				new_dist = (*d)[i][k] + vertex_nominal_time + (*d)[k][j];
				if (double_greater((*d)[i][j], new_dist))
					(*d)[i][j] = new_dist;
			}
		}
	}

	return d;
}
