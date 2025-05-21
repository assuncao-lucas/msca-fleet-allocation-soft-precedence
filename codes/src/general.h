#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

#include "src/matrix.hpp"

#define K_TYPE_CUT_TREATMENT IloCplex::CutManagement::UseCutForce // "Force" is the default

const bool K_GETOPT = true;
const bool K_STOP = true; // marks if problem considered is STOP (true) or TOP (false)
const bool K_SEPARATE_CUTS_ONLY_AT_ROOT = false;
const double K_PRECISION = 10000;
const double K_PRECISION_COMPARE_DOUBLE = 0.0000001;
const int K_NUM_TYPES_CALLBACKS = 9;

const double K_TAILING_OFF_TOLERANCE = 0.001;
const double K_GCCS_TOLERANCE = 0.05;
const double K_CONFLICT_TOLERANCE = 0.3;
const double K_CLIQUE_CONFLICT_TOLERANCE = 0.01;
const double K_CUTS_ANGLE_COSIN_LIMIT = 0.03;

const bool K_DISABLE_PREPROCESS_HEURISTICS = false;
const bool K_SET_PRIORITY_ORDER = false;
const bool K_MULTI_THREAD = true;
const bool K_ONLY_MAXIMUM_CLIQUES_PER_VERTEX = false;
const bool K_ADD_INITIAL_HEURISTIC_SOLUTION = true;
const bool K_PREPROCESS_REDUCED_COSTS = false;

const int K_TYPE_GCC_CUT = 0;
const int K_TYPE_FLOW_BOUNDS_CUT = 1;
const int K_TYPE_INITIAL_FLOW_BOUNDS_CUT = 2;
const int K_TYPE_CONFLICT_CUT = 3;
const int K_TYPE_COVER_BOUND_CUT = 4;
const int K_TYPE_CLIQUE_CONFLICT_CUT = 5;
const int K_TYPE_CLIQUE_CONFLICT_CUT_EXT = 6;
const int K_TYPE_COVER_BOUND_CUT_EXT = 7;
const int K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT = 8;

const int K_TYPE_PATH_BOUND_CUT = -1;

const double K_PATH_BOUND_CUT_TOLERANCE = 0.00001;
const double K_COVER_BOUND_CUT_TOLERANCE = 0.00001;

const std::string K_FILE_DELIMITER = "***********************************";

// feasibility pump parameters
const double K_PUMP_IMPROVEMENT_TOLERANCE = 0.1;
const int max_iter_stage1 = 1000;
const int max_iter_stage2 = 2000;
const int max_stalls_stage1 = 10;
const int max_stalls_stage2 = 10;
const double alpha_decrease_rate = 0.9;
const double initial_alpha_stage1 = 0.0;
const double initial_alpha_stage2 = 0.0; // if zero, then feasibility pump, if 1.0, then objective feasibility pump.
const bool K_SOLVE_STAGE_1 = false;
const double perturbation_flip_percentage = 0.1;
const int K_FLIP_BASIS = 10;
const bool FIXED_FLIP_BASIS = true;
const double K_APLHA_DECREMENT_PRECISION = 0.005;
const bool K_FEASIBILITY_PUMP_ADD_CUTS = false;

// DEFAULT Kernel Search parameters
const bool K_KERNEL_SEARCH_ADD_CUTS = false;
const int K_KS_MAX_SIZE_BUCKET = 5;
const int K_KS_MAX_TIME_LIMIT = 20;
const int K_KS_MIN_TIME_LIMIT = 5;
const double K_KS_DECAY_FACTOR_TIME_LIMIT = 0.9;

// Local Branching parameters
const int K_LOCAL_BRANCHING_K = 30;
const bool K_LOCAL_BRANCHING_ADD_CUTS = false;
const bool K_LOCAL_BRANCHING_SOLVE_ALNS = true;

// DEFAULT ALNS parameters
const int K_ALNS_SIZE_OF_POOL = 20;
const int K_ALNS_NUM_ITERATIONS = 2000;
const double K_ALNS_PERTURBATION_PERCENTAGE = 0.75;
const bool K_ALNS_MULTI_THREAD = true;

// DEFAULT SIMULATED ANNEALING PARAMETERS
const double K_SA_TEMP_DECAY_RATE = 0.98;
const double K_SA_TARGET_ACCEPTANCE_PROBABILITY = 0.9;
const size_t K_SA_SAMPLING_SIZE = 20;
const size_t K_SA_MIN_SAMPLING = 5;
const double K_SA_DEFAULT_INITIAL_TEMPERATURE = 1000.0;

// Path Relinking parameters
const bool K_PATH_RELINKING = false;
const int K_NON_IMPROVE_ITERATIONS_LIMIT = 100;
const double K_SIMILARITY_THRESHOLD = 0.9;
const bool K_STOP_AT_FIRST_PR_DETECTION = false;
// const int K_NUM_POOL_SOLUTIONS_SELECTED_FOR_PATH_RELINKING = 5;

std::vector<bool> *GetCallbackSelection(void);

double StDev(std::vector<double> &gaps, double b_avg_gap);
void AddInstancesFromDirectory(std::string directory, std::vector<std::string> &instances, bool add_dir = true);
void DeleteCallbackSelection(void);
void split_file_path(std::string directory, std::string &folder, std::string &file_name);
double euclidian_distance(std::pair<double, double> &coordinates_1, std::pair<double, double> &coordinates_2);
bool double_equals(double a, double b, double epsilon = K_PRECISION_COMPARE_DOUBLE);
bool double_greater(double a, double b, double epsilon = K_PRECISION_COMPARE_DOUBLE);
bool double_less(double a, double b, double epsilon = K_PRECISION_COMPARE_DOUBLE);
double round_decimals(double value, int num_decimals);

template <class T>
Matrix<T> **Allocate3Matrix(int k, int lines, int columns)
{
    Matrix<T> **matrix = new Matrix<T> *[k];
    for (int i = 0; i < k; i++)
    {
        matrix[i] = new Matrix<T>(lines, columns, 0);
    }

    return matrix;
}

template <class T>
void Print3Matrix(Matrix<T> **matrix, int k)
{
    if (matrix != NULL)
    {
        for (int i = 0; i < k; i++)
        {
            std::cout << "k = " << i << std::endl;
            if (matrix[i] != NULL)
            {
                matrix[i]->Print();
            }
        }
    }
}

template <class T>
void Delete3Matrix(Matrix<T> **matrix, int k)
{
    if (matrix != NULL)
    {
        for (int i = 0; i < k; i++)
        {
            if (matrix[i] != NULL)
            {
                delete matrix[i];
                matrix[i] = NULL;
            }
        }
        delete[] matrix;
    }
    matrix = NULL;
}

enum class Formulation
{
    baseline,
    single_commodity
};