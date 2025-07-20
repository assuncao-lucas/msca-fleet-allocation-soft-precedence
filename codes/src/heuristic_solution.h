#pragma once

#include <list>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "route.h"
#include "graph.h"
#include "instance.h"
#include "general.h"

class MetaHeuristicSolution;

class VertexStatus
{
public:
	VertexStatus();
	~VertexStatus();
	bool selected_;
	int route_;
	std::list<int>::iterator pos_;
};

class HeuristicSolution
{
public:
	explicit HeuristicSolution() = default;
	explicit HeuristicSolution(int num_vertices, int num_arcs, int num_routes);
	virtual ~HeuristicSolution();
	virtual void Reset(int num_vertices, int num_arcs, int num_routes);
	bool is_infeasible_ = false;
	bool is_feasible_ = false;
	bool is_optimal_ = false;
	double profits_sum_ = 0.0;
	double total_time_spent_ = 0.0;
	int num_vertices_ = 0;
	int num_arcs_ = 0;
	int num_routes_ = 0;

	bool operator==(HeuristicSolution &);

	friend std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol);

	virtual void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
	virtual void ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name);

	boost::dynamic_bitset<> bitset_arcs_;
	boost::dynamic_bitset<> bitset_vertices_;
};

class KSHeuristicSolution : public HeuristicSolution
{
public:
	KSHeuristicSolution();
	KSHeuristicSolution(int dimension, int dimension2, int num_routes);
	~KSHeuristicSolution();
	double time_spent_building_kernel_buckets_ = 0.0;
	bool found_x_integer_ = false;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	static std::string GenerateFileName(Formulation formulation, int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool feasibility_emphasis);
	void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
};