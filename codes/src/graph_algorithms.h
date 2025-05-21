#pragma once

#include <list>
#include <vector>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>
#include "src/graph.h"

double * Dijkstra(const Graph * graph, bool inverse, bool unit_dist);
Matrix<double> * FloydWarshall(const Graph * graph);
void BronKerosch(std::vector<std::list<int>> * adj_lists, std::list<int> R, std::unordered_set<int> P, std::unordered_set<int> X, std::vector<std::list<int>>& cliques); // finds all maximal cliques
void Tomita(std::vector<boost::dynamic_bitset<>> * adj_lists, std::list<int> R, boost::dynamic_bitset<> P, boost::dynamic_bitset<> X, std::vector<std::list<int>>& cliques); // finds all maximal cliques
