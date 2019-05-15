
#include "myNameClass.h"



myVector::myVector(Vertex initial_point, Vertex terminal_point) {

	this->x = terminal_point.x - initial_point.x;
	this->y = terminal_point.y - initial_point.y;
}