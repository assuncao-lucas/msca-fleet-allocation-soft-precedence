#include <math.h>
#include <utility>
#include <dirent.h>
#include <algorithm>
#include "general.h"

static std::vector<bool> *CALLBACKS_SELECTION = NULL;

std::vector<bool> *GetCallbackSelection(void)
{
    if (CALLBACKS_SELECTION == NULL)
    {
        CALLBACKS_SELECTION = new std::vector<bool>(K_NUM_TYPES_CALLBACKS, false);
    }

    return CALLBACKS_SELECTION;
}

void AddInstancesFromDirectory(std::string directory, std::vector<std::string> &instances, bool add_dir)
{
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL)
        {
            std::string curr_file(ent->d_name);
            if ((curr_file != ".") && (curr_file != "..") && ((curr_file.size() > 0) && (curr_file[curr_file.size() - 1] != '~')))
            {
                std::string inst_path;
                if (add_dir)
                    inst_path = directory;
                inst_path.append(curr_file);
                instances.push_back(inst_path);
            }
        }
        closedir(dir);
    }
    else
    {
        std::cout << "Invalid directory" << std::endl;
        throw 1;
    }
    // std::sort(instances.begin(),instances.end());
}

double StDev(std::vector<double> &gaps, double b_avg_gap)
{
    double stdev = 0.0;
    for (size_t i = 0; i < gaps.size(); i++)
    {
        stdev += pow(gaps.at(i) - b_avg_gap, 2.0);
    }

    stdev = stdev / (gaps.size() * 1.0 - 1.0);
    return pow(stdev, 0.5);
}

void split_file_path(std::string directory, std::string &folder, std::string &file_name)
{
    std::string reversed = directory;
    std::reverse(reversed.begin(), reversed.end());
    size_t pos = reversed.find_first_of("//");
    // size_t pos2 = reversed.find("PD-POTS-R"); // R-STOP-DP reversed
    folder = reversed.substr(pos);
    file_name = reversed.substr(0, pos);

    std::reverse(folder.begin(), folder.end());
    std::reverse(file_name.begin(), file_name.end());
    // std::cout << reversed << std::endl;
}

void DeleteCallbackSelection(void)
{
    if (CALLBACKS_SELECTION != NULL)
    {
        delete CALLBACKS_SELECTION;
        CALLBACKS_SELECTION = NULL;
    }
}

double euclidian_distance(std::pair<double, double> &coordinates_1, std::pair<double, double> &coordinates_2)
{
    double dist = pow(coordinates_2.second - coordinates_1.second, 2.0) + pow(coordinates_2.first - coordinates_1.first, 2.0);
    dist = pow(dist, 0.5);
    return dist;
}

bool double_equals(double a, double b, double epsilon)
{
    // if (a == std::numeric_limits<double>::infinity() && b == std::numeric_limits<double>::infinity())
    //     return true;
    return fabs(a - b) < epsilon;
}

bool double_greater(double a, double b, double epsilon)
{
    return a - b > epsilon;
}

bool double_less(double a, double b, double epsilon)
{
    return b - a > epsilon;
}

double round_decimals(double value, int num_decimals)
{
    auto factor = pow(10.0, 1.0 * num_decimals);
    return std::round(value * factor) / factor;
}
