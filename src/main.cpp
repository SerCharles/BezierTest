#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include "Bezier.hpp"
using namespace std;

/*
Generate a random place between [-range, range]
Args:
    range [int]: [the range of the place]

Returns:
    place [double]: [the result place between [-range, range]]
*/
double GenerateRandomPlace(double range)
{
    int raw_place = rand();
    raw_place = raw_place % (20000 + 1);
    raw_place = raw_place - 10000;
    double place = double(raw_place) / 10000.0 * range;
    return place;
}

/*
Generate a random t between [minumum, maximum]
Args:
    minimum [double]: [the minimum of t, between 0 and 1]
    maximum [double]: [the maximum of t, between 0 and 1]

Returns:
    t [double]: [a random number between 0 and 1]
*/
double GenerateRandomParameter(double minimum, double maximum)
{
    int raw_place = rand();
    raw_place = raw_place % 1001;
    double t = raw_place / 1000.0;
    double d = maximum - minimum;
    t = t * d + minimum;
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
Polynomial GenerateRandomPolynomial(int n, double range)
{
    vector<Point> points;
    points.clear();
    for(int i = 0; i <= n; i ++)
    {
        double x = GenerateRandomPlace(range);
        double y = GenerateRandomPlace(range);
        Point new_point = Point(x, y);
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
Bezier GenerateRandomBezier(int n, double range)
{
    vector<Point> points;
    points.clear();
    for(int i = 0; i <= n; i ++)
    {
        double x = GenerateRandomPlace(range);
        double y = GenerateRandomPlace(range);
        Point new_point = Point(x, y);
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
    double dist = ground_truth.dist(my_result);
    if(dist > threshold)
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
    int range = 1000000;
    double threshold = 1e-5;
    cout << "Starting random test of curve switch" << endl;
    cout << "There are ranks from 1 to 20, for each rank there are 10 examples, 200 examples in total" << endl;
    for(int n = 1; n <= 20; n ++)
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

/*
Randomly test the value calculation of Bezier
Args:
    range [int]: [the range of the point coordinates, [-range, range]]
    min_t [double]: [the minimum of t]
    max_t [double]: [the maximum of t]
*/
void TestRandomValueBezier(double range, double min_t, double max_t)
{
    ofstream outfile;
	outfile.open("bezier.txt", ios::out);
    ofstream outmax;
	outmax.open("bezier_max.txt", ios::out);
    int repeat_times = 1000;
    cout << "Starting random test of bezier curve value" << endl;
    cout << "There are ranks from 1 to 20, for each rank there are 100 examples, 2000 examples in total" << endl;
    for(int n = 1; n <= 20; n ++)
    {
        vector<Bezier> bezier_list;
        vector<Point> result_list_brute_force;
        vector<Point> result_list_de_casteljau;
        vector<Point> result_list_polynomial;
        vector<double> t_list;
        bezier_list.clear();
        result_list_brute_force.clear();
        result_list_de_casteljau.clear();
        result_list_polynomial.clear();
        t_list.clear(); 
        for(int i = 0; i < repeat_times; i ++)
        {
            Bezier bezier = GenerateRandomBezier(n, range);
            bezier_list.push_back(bezier);
            double t = GenerateRandomParameter(min_t, max_t);
            t_list.push_back(t);
        }
        int max_de_casteljau_i = 0;
        int max_polynomial_i = 0;
        

        //get results and times
        clock_t start_brute_force = clock();
        for(int i = 0; i < repeat_times; i ++)
        {
            Bezier bezier = bezier_list[i];
            double t = t_list[i];
            Point result = bezier.GetPlace(t);
            result_list_brute_force.push_back(result);
        }
        clock_t end_brute_force = clock();
        double time_brute_force = (double)(end_brute_force - start_brute_force) / CLOCKS_PER_SEC * 1000;

        clock_t start_de_casteljau = clock();
        for(int i = 0; i < repeat_times; i ++)
        {
            Bezier bezier = bezier_list[i];
            double t = t_list[i];
            Point result = bezier.DeCasteljau(t);
            result_list_de_casteljau.push_back(result);
        }
        clock_t end_de_casteljau = clock();
        double time_de_casteljau = (double)(end_de_casteljau - start_de_casteljau) / CLOCKS_PER_SEC * 1000;

        clock_t start_polynomial = clock();
        for(int i = 0; i < repeat_times; i ++)
        {
            Bezier bezier = bezier_list[i];
            double t = t_list[i];
            Polynomial polynomial = bezier.SwitchToPolynomial();
            Point result = polynomial.GetPlace(t);
            result_list_polynomial.push_back(result);
        }
        clock_t end_polynomial = clock();
        double time_polynomial = (double)(end_polynomial - start_polynomial) / CLOCKS_PER_SEC * 1000;

        //test whether right and error
        double dist_de_casteljau = 0.0;
        double max_de_casteljau = 0.0;
        double dist_polynomial = 0.0;
        double max_polynomial = 0.0;
        for(int i = 0; i < repeat_times; i ++)
        {
            Point result_brute_force = result_list_brute_force[i];
            Point result_de_casteljau = result_list_de_casteljau[i];
            Point result_polynomial = result_list_polynomial[i];
            double dist1 = result_brute_force.dist(result_de_casteljau);
            dist_de_casteljau += dist1;
            if(dist1 > max_de_casteljau)
            {
                max_de_casteljau = dist1;
                max_de_casteljau_i = i;
            }
            double dist2 = result_brute_force.dist(result_polynomial);
            dist_polynomial += dist2;
            if(dist2 > max_polynomial)
            {
                max_polynomial = dist2;
                max_polynomial_i = i;
            }
        }
        dist_de_casteljau = dist_de_casteljau / repeat_times;
        dist_polynomial = dist_polynomial / repeat_times;


        outmax << "n = " << n << endl;
        outmax << "de castelijau:" << endl;
        outmax << "t = " << t_list[max_de_casteljau_i] << endl;
        for(int ii = 0; ii < n; ii ++)
        {
            outmax << "[" << bezier_list[max_de_casteljau_i].ControlPoints[ii].x << "," << bezier_list[max_de_casteljau_i].ControlPoints[ii].y << "], " << endl;
        }
        outmax << endl;
        outmax << "polynomial:" << endl;
        outmax << "t = " << t_list[max_polynomial_i] << endl;
        for(int ii = 0; ii < n; ii ++)
        {
            outmax << "[" << bezier_list[max_polynomial_i].ControlPoints[ii].x << "," << bezier_list[max_polynomial_i].ControlPoints[ii].y << "], " << endl;
        }
        outmax << endl;
        outmax << "----------------------------------------------------------------";

        

        //show results
        cout << "Current rank: " << n << endl;
        cout << "The average dist of de castelijau: " << dist_de_casteljau << endl;
        cout << "The average dist of polynomial: " << dist_polynomial << endl;
        cout << "The max dist of de castelijau: " << max_de_casteljau << endl;
        cout << "The max dist of polynomial: " << max_polynomial << endl;   
        cout << "The average time of brute force: " << time_brute_force << " ms" << endl;
        cout << "The average time of de castelijau: " << time_de_casteljau << " ms" << endl;
        cout << "The average time of polynomial: " << time_polynomial << " ms" << endl;
        cout << "------------------------------------------------------------------------" << endl;
        outfile << n << ' ' << dist_de_casteljau << ' ' << dist_polynomial << ' ' << max_de_casteljau << ' ' << max_polynomial << ' ' << time_brute_force << ' ' << time_de_casteljau << ' ' << time_polynomial << endl;

    }
    outfile.close();
    outmax.close();
}

/*
Randomly test the value calculation of Polynomial
Args:
    range [int]: [the range of the point coordinates, [-range, range]]
    min_t [double]: [the minimum of t]
    max_t [double]: [the maximum of t]
*/
void TestRandomValuePolynomial(double range, double min_t, double max_t)
{
    ofstream outfile;
	outfile.open("polynomial.txt", ios::out);
    ofstream outmax;
	outmax.open("polynomial_max.txt", ios::out);
    int repeat_times = 1000;
    cout << "Starting random test of polynomial curve value" << endl;
    cout << "There are ranks from 1 to 20, for each rank there are 100 examples, 2000 examples in total" << endl;
    
    for(int n = 1; n <= 20; n ++)
    {
        vector<Polynomial> polynomial_list;
        vector<Point> result_list_brute_force;
        vector<Point> result_list_de_casteljau;
        vector<Point> result_list_bezier;
        vector<double> t_list;
        polynomial_list.clear();
        result_list_brute_force.clear();
        result_list_de_casteljau.clear();
        result_list_bezier.clear();
        t_list.clear(); 
        for(int i = 0; i < repeat_times; i ++)
        {
            Polynomial polynomial = GenerateRandomPolynomial(n, range);
            polynomial_list.push_back(polynomial);
            double t = GenerateRandomParameter(min_t, max_t);
            t_list.push_back(t);
        }
        int max_de_casteljau_i = 0;
        int max_bezier_i = 0;

        //get results and times
        clock_t start_brute_force = clock();
        for(int i = 0; i < repeat_times; i ++)
        {
            Polynomial polynomial = polynomial_list[i];
            double t = t_list[i];
            Point result = polynomial.GetPlace(t);
            result_list_brute_force.push_back(result);
        }
        clock_t end_brute_force = clock();
        double time_brute_force = (double)(end_brute_force - start_brute_force) / CLOCKS_PER_SEC * 1000;

        clock_t start_de_casteljau = clock();
        for(int i = 0; i < repeat_times; i ++)
        {
            Polynomial polynomial = polynomial_list[i];
            double t = t_list[i];
            Bezier bezier = polynomial.SwitchToBezier();
            Point result = bezier.DeCasteljau(t);
            result_list_de_casteljau.push_back(result);
        }
        clock_t end_de_casteljau = clock();
        double time_de_casteljau = (double)(end_de_casteljau - start_de_casteljau) / CLOCKS_PER_SEC * 1000;

        clock_t start_bezier = clock();
        for(int i = 0; i < repeat_times; i ++)
        {
            Polynomial polynomial = polynomial_list[i];
            double t = t_list[i];
            Bezier bezier = polynomial.SwitchToBezier();
            Point result = bezier.GetPlace(t);
            result_list_bezier.push_back(result);
        }
        clock_t end_bezier = clock();
        double time_bezier = (double)(end_bezier - start_bezier) / CLOCKS_PER_SEC * 1000;

        //test whether right and error
        double dist_de_casteljau = 0.0;
        double max_de_casteljau = 0.0;
        double dist_bezier = 0.0;
        double max_bezier = 0.0;
        for(int i = 0; i < repeat_times; i ++)
        {
            Point result_brute_force = result_list_brute_force[i];
            Point result_de_casteljau = result_list_de_casteljau[i];
            Point result_bezier = result_list_bezier[i];
            double dist1 = result_brute_force.dist(result_de_casteljau);
            dist_de_casteljau += dist1;
            if(dist1 > max_de_casteljau)
            {
                max_de_casteljau = dist1;
                max_de_casteljau_i = i;
            }
            max_de_casteljau = max(max_de_casteljau, dist_de_casteljau);
            double dist2 = result_brute_force.dist(result_bezier);
            dist_bezier += dist2;
            if(dist2 > max_bezier)
            {
                max_bezier = dist2; 
                max_bezier_i = i;
            }
        }
        dist_de_casteljau = dist_de_casteljau / repeat_times;
        dist_bezier = dist_bezier / repeat_times;

        outmax << "n = " << n << endl;
        outmax << "de castelijau:" << endl;
        outmax << "t = " << t_list[max_de_casteljau_i] << endl;
        for(int ii = 0; ii < n; ii ++)
        {
            outmax << "[" << polynomial_list[max_de_casteljau_i].ControlPoints[ii].x << "," << polynomial_list[max_de_casteljau_i].ControlPoints[ii].y << "], " << endl;
        }
        outmax << endl;
        outmax << "bezier:" << endl;
        outmax << "t = " << t_list[max_bezier_i] << endl;
        for(int ii = 0; ii < n; ii ++)
        {
            outmax << "[" << polynomial_list[max_bezier_i].ControlPoints[ii].x << "," << polynomial_list[max_bezier_i].ControlPoints[ii].y << "], " << endl;
        }
        outmax << endl;
        outmax << "----------------------------------------------------------------";


        //show results
        cout << "Current rank: " << n << endl;
        cout << "The average dist of de castelijau: " << dist_de_casteljau << endl;
        cout << "The average dist of bezier: " << dist_bezier << endl;
        cout << "The max dist of de castelijau: " << max_de_casteljau << endl;
        cout << "The max dist of bezier: " << max_bezier << endl;
        cout << "The average time of brute force: " << time_brute_force << " ms" << endl;
        cout << "The average time of de castelijau: " << time_de_casteljau << " ms" << endl;
        cout << "The average time of bezier: " << time_bezier << " ms" << endl;
        cout << "------------------------------------------------------------------------" << endl;
        outfile << n << ' ' << dist_de_casteljau << ' ' << dist_bezier << ' ' << max_de_casteljau << ' ' << max_bezier << ' ' << time_brute_force << ' ' << time_de_casteljau << ' ' << time_bezier << endl;
    }
    outfile.close();
    outmax.close();
}

int main()
{
    TestExamples();
    TestRandomSwitch();
    TestRandomValueBezier(1e6, 0, 1);    
    TestRandomValuePolynomial(1e6, 0, 1);
}