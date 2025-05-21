#include "user_cut.h"
#include "general.h"

void DeleteCuts(std::list<UserCutGeneral*> * cuts)
{
  if(cuts != NULL)
  {
    for(auto it = cuts->begin(); it != cuts->end(); ++it)
    {
      delete *it;
      *it = NULL;
    }
    cuts->clear();
  }
}

CoverBoundCut::CoverBoundCut()
{
    this->rhs_ = 0;
    this->type_ = K_TYPE_COVER_BOUND_CUT;
}

CoverBoundCut::CoverBoundCut(int num_vertices, std::list<int> cut_elements, int rhs)
{
    this->cut_exp_coefs_ = std::vector<int>(num_vertices,0);
    this->cut_elements_ = cut_elements;
    this->rhs_ = rhs;
    this->type_ = K_TYPE_COVER_BOUND_CUT;
}

CoverBoundCut::~CoverBoundCut(){}

PathBoundCut::PathBoundCut()
{
    this->rhs_ = 0;
    this->type_ = K_TYPE_PATH_BOUND_CUT;
}

PathBoundCut::PathBoundCut(std::list<std::pair<int,int>> cut_elements, int rhs)
{
    this->cut_elements_ = cut_elements;
    this->rhs_ = rhs;
    this->type_ = K_TYPE_PATH_BOUND_CUT;
}

PathBoundCut::~PathBoundCut(){}

UserCut::UserCut()
{
    this->norm_ = 0.0;
    this->curr_abs_violation_ = 0.0;
    this->curr_normalized_violation_ = 0.0;
	this->density_ = 0;
	this->lhs_dimension_ = 0;
	this->rhs_dimension_ = 0;
	this->type_ = -1;
}

UserCut::UserCut(int lhs_dimension, int rhs_dimension, double abs_violation, int type)
{
    this->norm_ = 0.0;
    this->curr_abs_violation_ = abs_violation;
    this->curr_normalized_violation_ = 0.0;
	this->density_ = 0;
	this->lhs_dimension_ = lhs_dimension;
	this->rhs_dimension_ = rhs_dimension;
	//this->lhs_coefficients_ = std::vector<int>(lhs_dimension, 0);
	this->lhs_ = boost::dynamic_bitset<>(lhs_dimension,0);
	this->rhs_ = boost::dynamic_bitset<>(rhs_dimension,0);
	this->type_ = type;
}

UserCut::~UserCut()
{
}

void UserCut::UpdateMeasures()
{
    // since any coefficient can only be 0, 1 or -1
    this->norm_ = sqrt(1.0*( (int)((this->lhs_nonzero_coefficients_indexes_).size()) + (int)((this->rhs_nonzero_coefficients_indexes_).size()) ));
    this->density_ = (1.0*( (int)((this->lhs_nonzero_coefficients_indexes_).size()) + (int)((this->rhs_nonzero_coefficients_indexes_).size()) ))/(1.0*(this->lhs_dimension_ + this->rhs_dimension_));
    this->curr_normalized_violation_ = (this->curr_abs_violation_)/(this->norm_);
}

double UserCut::operator* (UserCut &other)
{
    double cuts_inner_product = 0.0;
    //int num_repeated_rhs_nonzero_coefficients_indexes = 0;
    if((this->lhs_dimension_ != other.lhs_dimension_) || ((this->rhs_dimension_ != other.rhs_dimension_)) )
    {
        throw 2;
    }
    else
    {
        /*for(int i = 0; i < this->lhs_dimension_; i++)
        {
           cuts_inner_product += (((this->lhs_coefficients_)[i]) * (other.lhs_coefficients_)[i]);
        }

        for(std::list<int>::iterator it = (this->rhs_nonzero_coefficients_indexes_).begin(); it!= (this->rhs_nonzero_coefficients_indexes_).end(); it++)
        {
            for(std::list<int>::iterator it2 = (other.rhs_nonzero_coefficients_indexes_).begin(); it2!= (other.rhs_nonzero_coefficients_indexes_).end(); it2++)
            {
                if(*it == *it2) num_repeated_rhs_nonzero_coefficients_indexes++;
            }
        }*/

	/*for(int i = 0; i < this->lhs_dimension_; i++)
        {
           cuts_inner_product += (((this->lhs_)[i]) * (other.lhs_)[i]);
        }

        for(int i = 0; i < this->rhs_dimension_; i++)
        {
           cuts_inner_product += (((this->rhs_)[i]) * (other.rhs_)[i]);
        }*/

        //cuts_inner_product += 1.0*num_repeated_rhs_nonzero_coefficients_indexes;
	cuts_inner_product = 1.0*(((this->lhs_)&(other.lhs_)).count() + ((this->rhs_)&(other.rhs_)).count());
    }
	//if(double_equals(cuts_inner_product, ((this->lhs_)&(other.lhs_)).count() + ((this->rhs_)&(other.rhs_)).count())) std::cout << "*funfou!!" << std::endl;

    return cuts_inner_product/((this->norm_)*(other.norm_));
}

bool UserCut::isBetterThan(UserCut * other)
{
    if(other == NULL) return true;
    return double_greater(this->curr_normalized_violation_,other->curr_normalized_violation_);
}

void UserCut::AddRhsElement(int pos)
{
	(this->rhs_nonzero_coefficients_indexes_).push_back(pos);
	(this->rhs_)[pos] = 1;
}

void UserCut::AddLhsElement(int v1, int v2, int pos)
{
	(this->lhs_nonzero_coefficients_indexes_).push_back(std::pair<int,int>(v1,v2));
        (this->lhs_)[pos] = 1;
}

std::ostream& operator<< (std::ostream &out, UserCut &cut)
{
    out << std::endl;
    out << "Norm: " << cut.norm_ << std::endl
        << "Abs violation: " << cut.curr_abs_violation_ << std::endl
        << "Normalized violation: " << cut.curr_normalized_violation_ << std::endl
        << "Density: " << cut.density_ << std::endl;

    return out;
}
