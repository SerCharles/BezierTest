#include <iostream>
#include <fstream>
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
Generate a random t between [0, 1]
Returns:
    t [double]: [a random number between 0 and 1]
*/
double GenerateRandomParameter()
{
    int raw_place = rand();
    raw_place = raw_place % 1001;
    double t = raw_place / 1000.0;
    return t;
}

/*
Generate a random Polynomial curve

Args:
    n [int]: [the n of the curve]
    range [int]: [the range of the points' x and y must be inside [-range, range]]

Returns:
    curve [Polynomial]: [the Polynomial curve]
*/
Polynomial GenerateRandomPolynomial(int n, int range)
{
    vector<Point> points;
    points.clear();
    for(int i = 0; i <= n; i ++)
    {
        double x = GenerateRandomPlace(range);
        double y = GenerateRandomPlace(range);
        Point new_point = Point(x, y);\
        points.push_back(new_point);
    }
    Polynomial curve = Polynomial(points);
    return curve;
}

/*
Generate a random Bezier curve

Args:
    n [int]: [the n of the curve]
    range [int]: [the range of the points' x and y must be inside [-range, range]]

Returns:
    curve [Bezier]: [the Bezier curve]
*/
Bezier GenerateRandomBezier(int n, int range)
{
    vector<Point> points;
    points.clear();
    for(int i = 0; i <= n; i ++)
    {
        double x = GenerateRandomPlace(range);
        double y = GenerateRandomPlace(range);
        Point new_point = Point(x, y);\
        points.push_back(new_point);
    }
    Bezier curve = Bezier(points);
    return curve;
}

/*
Test whether the curve result is right
Args:
    ground_truth [Point]: [the ground truth point]
    my_result [Point]: [the result point]
    threshold [double]: [the threshold of error]
Returns:
    result [bool]: [1 is right, 0 is wrong]
*/
bool TestValue(Point& ground_truth, Point& my_result, double threshold)
{
    bool result = 1;
    if(abs(ground_truth.x - my_result.x) > threshold)
    {
        result = 0;
    }
    if(abs(ground_truth.y - my_result.y) > threshold)
    {
        result = 0;
    }
    return result;
}

/*
Test whether the curve result is right
Args:
    ground_truth [Bezier or Polynomial]: [the ground truth curve]
    my_result [Bezier or Polynomial]: [the result curve]
    threshold [double]: [the threshold of error]
Returns:
    result [bool]: [1 is right, 0 is wrong]
*/
template <class T>
bool TestCurve(T& ground_truth, T& my_result, double threshold)
{
    bool result = 1;
    int n;
    if(ground_truth.n != my_result.n)    
    {
        return 0;
    }
    else
    {
        n = ground_truth.n;
    }
    for(int i = 0; i <= n; i ++)
    {
        bool the_result = TestValue(ground_truth.ControlPoints[i], my_result.ControlPoints[i], threshold);
        if(the_result == 0)
        {
            result = 0;
            break;
        }
    }
    return result;
}


/*
Test the artificial easy examples of Bezier and Polynomial curves
*/
void TestExamples()
{
    ifstream infile;
	infile.open("example.txt", ios::in);
    int m;
    infile >> m;
    cout << "Starting example test of curve switch, the examples are in example.txt."<< endl;
    cout << "There are " << m << " examples in total." << endl;
    bool result = 1;
    double threshold = 1e-5;
    for(int i = 1; i <= m; i ++)
    {
        int n;
        infile >> n;
        vector<Point> bezier_points;
        vector<Point> polynomial_points;
        bezier_points.clear();
        polynomial_points.clear();
        for(int j = 0; j <= n; j ++)
        {
            double x, y;
            infile >> x >> y;
            Point new_point = Point(x, y);
            bezier_points.push_back(new_point);
        }
        for(int j = 0; j <= n; j ++)
        {
            double x, y;
            infile >> x >> y;
            Point new_point = Point(x, y);
            polynomial_points.push_back(new_point);
        }
        Bezier bezier_gt = Bezier(bezier_points);
        Polynomial polynomial_gt = Polynomial(polynomial_points);
        Polynomial polynomial_mine = bezier_gt.SwitchToPolynomial();
        Bezier bezier_mine = polynomial_mine.SwitchToBezier();
        bool the_result = 1;
        the_result = the_result & TestCurve(polynomial_gt, polynomial_mine, threshold);
        the_result = the_result & TestCurve(bezier_gt, bezier_mine, threshold);
        result = result & the_result;
        if(the_result)
        {
            cout << "Example " << i << " is right." << endl;
        }
        else 
        {
            cout << "Example " << i << " is wrong." << endl;
        }
    }
    if(result)
    {
        cout << "All the examples are right" << endl;
    }
    else 
    {
        cout << "Not all the examples are right" << endl;
    }
    infile.close();
    cout << "Example test of curve switch has ended."<< endl;
}

/*
Randomly test the switch between Bezier and Polynomial
*/
void TestRandomSwitch()
{
    bool result = 1;
    int range = 1000;
    double threshold = 1e-5;
    cout << "Starting random test of curve switch" << endl;
    cout << "There are ranks from 1 to 10, for each rank there are 10 examples, 100 examples in total" << endl;
    for(int n = 1; n <= 10; n ++)
    {
        for(int i = 1; i <= 10; i ++)
        {
            Bezier bezier_gt = GenerateRandomBezier(n, range);
            Polynomial polynomial_gt = bezier_gt.SwitchToPolynomial();
            Bezier bezier_mine = polynomial_gt.SwitchToBezier();
            Polynomial polynomial_mine = bezier_mine.SwitchToPolynomial();
            bool the_result = 1;
            the_result = the_result & TestCurve(polynomial_gt, polynomial_mine, threshold);
            the_result = the_result & TestCurve(bezier_gt, bezier_mine, threshold);
            result = result & the_result;
            if(the_result)
            {
                cout << "Example " << i << " with rank " << n << " is right." << endl;
            }
            else 
            {
                cout << "Example " << i << " with rank " << n << " is wrong." << endl;
            }
        } 
    }
    if(result)
    {
        cout << "All the examples are right" << endl;
    }
    else 
    {
        cout << "Not all the examples are right" << endl;
    }
    cout << "Random test of curve switch has ended."<< endl;
}

int main()
{
    TestExamples();
    TestRandomSwitch();
}