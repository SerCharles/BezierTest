#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Bezier.hpp"
using namespace std;


int main()
{
    Point a(0.0, 0.0);
    Point b(1.0, 1.0);
    Point c(2.0, 0.0);
    vector<Point> point_list;
    point_list.push_back(a);
    point_list.push_back(b);
    point_list.push_back(c);
    Bezier p = Bezier(point_list);
    Polynomial polynomial = p.SwitchToPolynomial();
    Bezier bezier = polynomial.SwitchToBezier();
    int kebab = 0;
    return 0;
}