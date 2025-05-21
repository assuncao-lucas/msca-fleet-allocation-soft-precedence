#include "src/arc.h"

GArc::GArc(double distance):distance_(distance){}

void GArc::set_distance(double distance)
{
    distance_ = distance;
}

std::ostream& operator<< (std::ostream &out, const GArc &arc)
{
    // Since operator<< is a friend of the GArc class, we can access
    // GArc's members directly.
	out << "[" << arc.distance() << "]";
    return out;
}
