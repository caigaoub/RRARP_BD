#ifndef _MYNAMECLASS_H_
#define _MYNAMECLASS_H_
#include <algorithm>
#include <math.h>
using namespace std;

struct Vertex {
	double x;
	double y;
	Vertex& operator= (Vertex other) {
		swap(x, other.x);
		swap(y, other.y);
		return *this;
	}
};

class myVector {
private:
	double x;
	double y;

public:
	myVector(Vertex, Vertex);
	~myVector() {};
	inline double get_x() { return x; };
	inline double get_y() { return y; };
	inline double get_vecLen() { return sqrt(x*x + y*y); };
};



#endif // !MYNAMECLASS_H_
