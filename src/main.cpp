#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Bezier.hpp"
using namespace std;

/*
Generate a random place between [-range, range]
Args:
    range [int]: [the range of the place]

Returns:
    place [double]: [the result place between [-range, range]]
*/
double GenerateRandomPlace(int range)
{
    int raw_place = rand();
    raw_place = raw_place % (2 * range + 1);
    raw_place = raw_place - range;
    double place = double(raw_place);
    return place;
}

/*
Generate a random polynomial curve

Args:
    n [int]: [the n of the curve]
    range [int]: [the range of the points' x and y must be inside [-range, range]]

Returns:
    polynomial [Polynomial]: [the polynomial curve]
*/
Polynomial GenerateRandomPolynomial(int n, int range)
{
    vector<Point> polynomial_points;
    polynomial_points.clear();
    for(int i = 0; i <= n; i ++)
    {
        double x = GenerateRandomPlace(range);
        double y = GenerateRandomPlace(range);
        Point new_point = Point(x, y);\
        polynomial_points.push_back(new_point);
    }
    Polynomial polynomial = Polynomial(polynomial_points);
    return polynomial;
}

/*
Generate a random bezier curve

Args:
    n [int]: [the n of the curve]
    range [int]: [the range of the points' x and y must be inside [-range, range]]

Returns:
    bezier [Bezier]: [the bezier curve]
*/
Bezier GenerateRandomBezier(int n, int range)
{
    vector<Point> bezier_points;
    bezier_points.clear();
    for(int i = 0; i <= n; i ++)
    {
        double x = GenerateRandomPlace(range);
        double y = GenerateRandomPlace(range);
        Point new_point = Point(x, y);\
        bezier_points.push_back(new_point);
    }
    Bezier bezier = Bezier(bezier_points);
    return bezier;
}

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