#pragma once
#include <utility>
class Edge{
    
public:
	//eid
	unsigned short eid;

    // Endpoints (a.id<b.id)
	unsigned int a;
	unsigned int b;

    // Weight
	float w;

	Edge(){ a = -1; b = -1; eid = -1; } //initial Null Edge

	Edge(unsigned char v1, unsigned char v2, float myw=0){
		a = v1 < v2 ? v1 : v2;
		b = v1 > v2 ? v1 : v2;
		w = myw;
	}

    inline unsigned char getOtherEnd(unsigned char i){
        if (i==a)
            return b;
        else if (i==b)
            return a;
        else
            return -1;
    }

	inline std::pair<unsigned char, unsigned char> id(){
		return std::pair<unsigned char, unsigned char>(a, b);
	}
};