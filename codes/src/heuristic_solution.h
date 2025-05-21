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

	std::vector<Route> routes_vec_;
	std::vector<VertexStatus> vertex_status_vec_;
	std::list<int> unvisited_vertices_;
	bool PreviewAddVertex(Instance &, int vertex, int route, int pos, std::list<int>::iterator &it, double &profit_variation, double &time_variation);
	void AddVertex(int vertex, int route, std::list<int>::iterator pos_it, double profit_variation, double time_variation);
	bool PreviewRemoveVertex(Instance &, int vertex, double &profit_variation, double &time_variation, bool allow_infeasible_routes = false);
	void RemoveVertex(int vertex, double profit_variation, double time_variation);
	bool PreviewInterRouteMoveVertex(Instance &, int vertex, int r2, int pos, std::list<int>::iterator &it, double &profit_variation1, double &time_variation1, double &profit_variation2, double &time_variation2);
	void InterRouteMoveVertex(int vertex, int r2, std::list<int>::iterator it, double profit_variation1, double time_variation1, double profit_variation2, double time_variation2);
	bool PreviewInterRouteSwap(Instance &, int r1, int pos_i1, int pos_f1, int r2, int pos_i2, int pos_f2, std::list<int>::iterator &it_i1, std::list<int>::iterator &it_f1,
							   double &profit_variation1, double &time_variation1, std::list<int>::iterator &it_i2, std::list<int>::iterator &it_f2, double &profit_variation2, double &time_variation2);
	void InterRouteSwap(int r1, int r2, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1,
						double profit_variation1, double time_variation1, std::list<int>::iterator it_i2, std::list<int>::iterator it_f2, double profit_variation2, double time_variation2);

	bool PreviewInterRouteSwapUnrouted(Instance &inst, int r1, int pos_i1, int pos_f1, std::list<int>::iterator &it_i1, std::list<int>::iterator &it_f1, double &profit_variation1, double &time_variation1, std::list<int>::iterator it_i2);
	void InterRouteSwapUnrouted(int r1, std::list<int>::iterator it_i1, std::list<int>::iterator it_f1, double profit_variation1, double time_variation1, std::list<int>::iterator it_i2);

	bool PreviewAddVertexToRouteWithinMaximumProfitIncrease(Instance &, int vertex, int route, std::list<int>::iterator &it, double &profit_variation, double &time_variation, bool allow_infeasible_routes = false);
	bool PreviewAddVertexWithinMaximumProfitIncrease(Instance &, int vertex, int &route, std::list<int>::iterator &it, double &profit_variation, double &time_variation);

	bool Do2OptImprovement(const Instance &inst, const int &route);
	bool operator==(HeuristicSolution &);

	friend std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol);
	bool CheckCorrectness(const Instance &instance);

	double ComputeSolutionCost(Instance &instance) const;
	double ComputeSolutionCostRec(Instance &instance, bool memoization = false) const;

	virtual void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
	virtual void ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name);

	boost::dynamic_bitset<> bitset_arcs_;
	boost::dynamic_bitset<> bitset_vertices_;

	void BuildBitset(const Instance &);
};

class FPHeuristicSolution : public HeuristicSolution
{
public:
	FPHeuristicSolution();
	FPHeuristicSolution(int dimension, int dimension2, int num_routes);
	~FPHeuristicSolution();
	int num_iterations_stage1_;
	int num_iterations_stage2_;
	int num_perturbations_stage1_;
	int num_perturbations_stage2_;
	int num_restarts_stage1_;
	int num_restarts_stage2_;
	bool found_y_integer_;
	bool found_x_integer_;
	double time_stage1_;
	double time_stage2_;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	static std::string GenerateFileName();
	void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
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

class LBHeuristicSolution : public HeuristicSolution
{
public:
	LBHeuristicSolution();
	LBHeuristicSolution(int dimension, int dimension2, int num_routes);
	~LBHeuristicSolution();
	int num_iterations_;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	static std::string GenerateFileName();
	void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
};

class MetaHeuristicSolution : public HeuristicSolution
{
public:
	MetaHeuristicSolution();
	MetaHeuristicSolution(const MetaHeuristicSolution *);
	MetaHeuristicSolution(const HeuristicSolution *);
	MetaHeuristicSolution(int dimension, int dimension2, int num_routes);
	~MetaHeuristicSolution();
	int num_iterations_ = 0;
	double initial_solution_profits_sum_ = 0.0;
	int last_improve_iteration_ = 0;

	virtual void Reset(int dimension, int dimension2, int num_routes);
	void WriteToFile(Instance &instance, std::string algo, std::string folder, std::string file_name) const;
	void ReadFromFile(Instance &inst, std::string algo, std::string folder, std::string file_name);
	static std::string GenerateALNSFileName(int num_iterations, int pool_size, int multithreading);
	static std::string GenerateSimulatedAnnealingFileName(double temo_decay_rate, int multithreading);
};
