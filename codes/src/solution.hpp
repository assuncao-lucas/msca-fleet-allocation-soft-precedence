#ifndef SOLUTION_PATH_H
#define SOLUTION_PATH_H
#include <iostream>
// #include <sys/stat.h>
// #include <sys/types.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <unordered_map>
#include "src/matrix.hpp"
#include "src/timer.h"
#include "src/general.h"
#include "src/user_cut.h"

template <class T>
class Solution
{
public:
	Solution(void);
	Solution(int);
	~Solution(void);
	Matrix<T> *solution();
	void set_cut_found(int type, bool root);
	void set_cut_added(int type, bool root);

	int dimension();
	double lb_;
	double ub_;
	double lp_;
	int num_calls_to_callback_;
	int num_calls_to_callback_lp_;
	int num_tailing_offs_;
	int num_nodes_;
	int num_conflicts_ = 0;
	int num_maximal_cliques_;
	int num_benders_opt_cuts_;
	int num_benders_feas_cuts_;
	double separation_time_;
	double root_time_;
	double milp_time_;
	double pre_processing_time_;
	bool is_optimal_;
	bool is_feasible_;
	bool out_of_memory_;
	bool is_root_node_;
	bool stop_adding_cuts_;
	Timestamp *ti();
	Timestamp *tf();
	void reset();
	void increment_elaspsed_time(double time);
	void write_to_file(std::string algo, std::string folder, std::string file_name);
	static std::string GenerateAlgorithmName(bool solve_relaxed, bool solve_cutting_plane, bool solve_benders, bool cccs, bool avics, bool baseline, bool capacity_based, bool benders_combine_feas_opt_cuts, bool benders_separate_benders_cuts_relaxation, bool benders_solve_generic_callback);

	std::vector<int> num_cuts_found_;
	std::vector<int> num_cuts_added_;

	std::vector<int> num_cuts_found_lp_;
	std::vector<int> num_cuts_added_lp_;

private:
	Matrix<T> *solution_;
	int dimension_;
	Timestamp *ti_;
	Timestamp *tf_;

	std::unordered_map<int, std::mutex> mutex_cut_added_vec_;
	std::unordered_map<int, std::mutex> mutex_cut_found_vec_;
	std::mutex mutex_time_count_;
};

template <typename T>
Solution<T>::Solution(void)
{
	this->solution_ = NULL;
	this->dimension_ = 0;
	this->is_optimal_ = false;
	this->is_feasible_ = true;
	this->out_of_memory_ = false;
	this->lb_ = -std::numeric_limits<double>::infinity();
	this->ub_ = std::numeric_limits<double>::infinity();
	this->lp_ = std::numeric_limits<double>::infinity();
	this->num_calls_to_callback_ = 0;
	this->num_calls_to_callback_lp_ = 0;
	this->separation_time_ = 0.0;
	this->root_time_ = 0.0;
	this->milp_time_ = 0.0;
	this->pre_processing_time_ = 0.0;
	this->num_nodes_ = 0;
	this->num_conflicts_ = 0;
	this->num_tailing_offs_ = 0;
	this->num_maximal_cliques_ = 0;
	this->num_cuts_added_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);
	this->num_cuts_found_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);

	this->num_cuts_added_lp_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);
	this->num_cuts_found_lp_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);

	this->num_benders_opt_cuts_ = 0;
	this->num_benders_feas_cuts_ = 0;

	this->is_root_node_ = true;
	this->stop_adding_cuts_ = false;
}

template <typename T>
Solution<T>::Solution(int dimension)
{
	this->solution_ = new Matrix<T>(dimension, dimension, 0);
	this->dimension_ = dimension;
	this->is_optimal_ = false;
	this->is_feasible_ = true;
	this->out_of_memory_ = false;
	this->lb_ = -std::numeric_limits<double>::infinity();
	this->ub_ = std::numeric_limits<double>::infinity();
	this->lp_ = std::numeric_limits<double>::infinity();
	this->num_calls_to_callback_ = 0;
	this->num_calls_to_callback_lp_ = 0;
	this->separation_time_ = 0.0;
	this->root_time_ = 0.0;
	this->milp_time_ = 0.0;
	this->pre_processing_time_ = 0.0;
	this->num_nodes_ = 0;
	this->num_conflicts_ = 0;
	this->num_tailing_offs_ = 0;
	this->num_maximal_cliques_ = 0;
	this->num_cuts_added_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);
	this->num_cuts_found_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);

	this->num_cuts_added_lp_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);
	this->num_cuts_found_lp_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);

	this->num_benders_opt_cuts_ = 0;
	this->num_benders_feas_cuts_ = 0;

	this->is_root_node_ = true;
	this->stop_adding_cuts_ = false;
}

template <typename T>
void Solution<T>::reset()
{
	if (this->solution_ != NULL)
	{
		delete (this->solution_);
		this->solution_ = NULL;
	}

	if (this->dimension_ > 0)
		this->solution_ = new Matrix<T>(this->dimension_, this->dimension_, 0);
	else
		this->solution_ = NULL;
	this->is_optimal_ = false;
	this->out_of_memory_ = false;
	this->is_feasible_ = true;
	this->lb_ = -std::numeric_limits<double>::infinity();
	this->ub_ = std::numeric_limits<double>::infinity();
	this->lp_ = std::numeric_limits<double>::infinity();
	this->num_calls_to_callback_ = 0;
	this->num_calls_to_callback_lp_ = 0;
	this->separation_time_ = 0.0;
	this->root_time_ = 0.0;
	this->milp_time_ = 0.0;
	this->pre_processing_time_ = 0.0;
	this->num_cuts_added_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);
	this->num_cuts_found_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);

	this->num_cuts_added_lp_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);
	this->num_cuts_found_lp_ = std::vector<int>(K_NUM_TYPES_CALLBACKS, 0);

	this->num_benders_opt_cuts_ = 0;
	this->num_benders_feas_cuts_ = 0;

	this->num_nodes_ = 0;
	this->num_conflicts_ = 0;
	this->num_tailing_offs_ = 0;
	this->num_maximal_cliques_ = 0;
	this->is_root_node_ = true;
	this->stop_adding_cuts_ = false;
}

template <typename T>
Solution<T>::~Solution(void)
{
	if (this->solution_ != NULL)
	{
		delete (this->solution_);
		this->solution_ = NULL;
	}
}

template <typename T>
Matrix<T> *Solution<T>::solution()
{
	return this->solution_;
}

template <typename T>
int Solution<T>::dimension()
{
	return this->dimension_;
}

template <typename T>
void Solution<T>::set_cut_added(int type, bool root)
{
	((this->mutex_cut_added_vec_)[type]).lock();
	if (!root)
		((this->num_cuts_added_)[type])++;
	else
		((this->num_cuts_added_lp_)[type])++;
	((this->mutex_cut_added_vec_)[type]).unlock();
}

template <typename T>
void Solution<T>::set_cut_found(int type, bool root)
{
	((this->mutex_cut_found_vec_)[type]).lock();
	if (!root)
		((this->num_cuts_found_)[type])++;
	else
		((this->num_cuts_found_lp_)[type])++;
	((this->mutex_cut_found_vec_)[type]).unlock();
}

template <typename T>
void Solution<T>::increment_elaspsed_time(double time)
{
	(this->mutex_time_count_).lock();
	(this->separation_time_) += time;
	(this->mutex_time_count_).unlock();
}

template <typename T>
std::string Solution<T>::GenerateAlgorithmName(bool solve_relaxed, bool solve_cutting_plane, bool solve_benders, bool cccs, bool avics, bool baseline, bool capacity_based, bool benders_combine_feas_opt_cuts, bool benders_separate_benders_cuts_relaxation, bool benders_solve_generic_callback)
{
	std::string algo;
	if (solve_relaxed)
		algo += "relax_";
	if (solve_cutting_plane)
	{
		algo += "cb_";
		if (avics)
			algo += "AVICs_";
		if (cccs)
			algo += "CCCs_";
	}
	if (solve_benders)
	{
		if (benders_combine_feas_opt_cuts)
			algo += "benders_combined_";
		else
			algo += "benders_";

		if (benders_separate_benders_cuts_relaxation)
			algo += "cuts_relaxation_";

		if (benders_solve_generic_callback)
			algo += "generic_callback_";
		else
			algo += "lazy_callback_";
	}
	if (baseline)
		algo += "baseline";
	if (capacity_based)
		algo += "csc";

	return algo;
}

template <typename T>
void Solution<T>::write_to_file(std::string algo, std::string folder, std::string file_name)
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

	file << std::setprecision(2) << std::fixed;

	if (this->out_of_memory_)
		file << "STATUS: OUT_OF_MEMORY" << std::endl;
	else if (!(this->is_feasible_))
		file << "STATUS: INFEASIBLE" << std::endl;
	else if (this->is_optimal_)
		file << "STATUS: OPTIMAL" << std::endl;
	else
		file << "STATUS: FEASIBLE" << std::endl;

	file << "LP: " << this->lp_ << std::endl;
	file << "Lb: " << this->lb_ << std::endl;
	file << "Ub: " << this->ub_ << std::endl;
	file << "Root time(s): " << this->root_time_ << std::endl;
	file << "MILP time(s): " << this->milp_time_ << std::endl;
	file << "Separation time(s): " << this->separation_time_ << std::endl;
	file << "# iterations of cutting plane at root: " << this->num_calls_to_callback_lp_ << std::endl;
	file << "# cuts added/found via callback (LP): " << std::endl;
	for (int i = 0; i < K_NUM_TYPES_CALLBACKS; i++)
		file << "type " << i << ": " << (this->num_cuts_added_lp_)[i] << "/" << (this->num_cuts_found_lp_)[i] << std::endl;
	file << "# calls to callback: " << this->num_calls_to_callback_ << std::endl;
	file << "# cuts added/found via callback (MILP): " << std::endl;
	for (int i = 0; i < K_NUM_TYPES_CALLBACKS; i++)
		file << "type " << i << ": " << (this->num_cuts_added_)[i] << "/" << (this->num_cuts_found_)[i] << std::endl;
	file << "# nodes explored: " << this->num_nodes_ << std::endl;
	file << "# tailing offs treated: " << this->num_tailing_offs_ << std::endl;
	file << "# maximal cliques: " << this->num_maximal_cliques_ << std::endl;
	file << "Pre Processing time (s): " << this->pre_processing_time_ << std::endl;
	file << "# num conflicts: " << this->num_conflicts_ << std::endl;

	file.close();
}

#endif
