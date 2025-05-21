#include <iomanip>
#include "route.h"

Route::Route()
{
	this->sum_profits_ = 0;
	this->time_ = 0.0;
}

Route::Route(const Route &route)
{
	this->sum_profits_ = route.sum_profits_;
	this->time_ = route.time_;
	this->vertices_ = route.vertices_;
}

Route::~Route()
{
}

std::ostream &operator<<(std::ostream &out, Route &route)
{
	out << std::setprecision(2) << std::fixed;
	out << "[size: " << (route.vertices_).size() << "][profits: " << route.sum_profits_ << "][time: " << route.time_ << "]" << std::endl;
	for (std::list<int>::iterator it = (route.vertices_).begin(); it != (route.vertices_).end(); ++it)
		out << (*it) << " ";
	out << std::endl;
	return out;
}
