#ifndef USERCUT_H_
#define USERCUT_H_

#include <iostream>
#include <vector>
#include <list>
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

class UserCut
{
public:
    UserCut(int lhs_dimension, double abs_violation, int type);
    UserCut();
    ~UserCut();
    void AddLhsElement(int v1, int v2, int pos);
    void AddRhsElement(int pos);

    void UpdateMeasures();
    bool isBetterThan(UserCut *other);
    double norm_;
    double curr_abs_violation_;
    double curr_normalized_violation_;
    double density_;
    int lhs_dimension_;
    // int lhs_num_nonzero_coefs_;
    boost::dynamic_bitset<> lhs_;
    // std::vector<int> lhs_coefficients_;
    std::list<std::pair<int, int>> lhs_nonzero_coefficients_indexes_;
    double operator*(UserCut &other);
    friend std::ostream &operator<<(std::ostream &out, UserCut &cut);
};

void DeleteCuts(std::list<UserCut *> *cuts);

#endif
