#pragma once 

#include <iostream>

class GArc final
{
	private:
		double distance_ = 0.0;
	public:
		explicit GArc() = default;
		explicit GArc(double);
		virtual ~GArc() = default;
		void set_distance(double distance);
		double distance() const {return distance_;}
		friend std::ostream& operator<< (std::ostream &out, const GArc &arc);
};
