#pragma once

#include <iostream>
#include <list>

class Route
{
public:
	Route();
	Route(const Route &route);
	~Route();
	double sum_profits_ = 0;
	double time_ = 0.0;
	std::list<int> vertices_;
	friend std::ostream &operator<<(std::ostream &out, Route &route);
};
