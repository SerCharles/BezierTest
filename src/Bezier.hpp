#include<vector>
using namespace std;

class Point
{
public:
    double x;
    double y;
    Point()
    {
        x = 0.0;
        y = 0.0;
    }
    Point(double a, double b)
    {
        x = a;
        y = b;
    }
    Point operator+(const Point& a) const
	{
        double new_x = x + a.x;
        double new_y = y + a.y;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
    Point operator-(const Point& a) const
	{
        double new_x = x - a.x;
        double new_y = y - a.y;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
	void operator+=(const Point& a)
	{
		x += a.x;
        y += a.y;
	}
	void operator-=(const Point& a)
	{
		x -= a.x;
        y -= a.y;
	}
	Point operator*(const double& a) const
	{
		double new_x = x * a;
        double new_y = y * a;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
	Point operator/(const double& a) const
	{
		double new_x = x / a;
        double new_y = y / a;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
	void operator*=(const double& a)
	{
		x *= a;
        y *= a;	
    }
	void operator/=(const double& a)
	{
        x /= a;
        y /= a;	
    }
    double operator^(const Point& b) const 
    {
        double result = x * b.y - y * b.x;
        return result;
    }
};

class Polynomial;
class Bezier
{
public:
    int n;
    vector<Point> ControlPoints;
    Bezier()
    {
        n = -1;
        ControlPoints.clear();
    }
    Bezier(vector<Point> control_points)
    {
        n = control_points.size() - 1;
        for(int i = 0; i <= n; i ++)
        {
            ControlPoints.push_back(control_points[i]);
        }
    }
    ~Bezier()
    {
        ControlPoints.clear();
    }
    Polynomial SwitchToPolynomial();
    
};

class Polynomial
{
public:
    int n;
    vector<Point> ControlPoints;
    Polynomial()
    {
        n = -1;
        ControlPoints.clear();
    }
    Polynomial(vector<Point> control_points)
    {
        n = control_points.size() - 1;
        for(int i = 0; i <= n; i ++)
        {
            ControlPoints.push_back(control_points[i]);
        }
    }
    ~Polynomial()
    {
        ControlPoints.clear();
    }
    Bezier SwitchToBezier();

};

Polynomial Bezier::SwitchToPolynomial()
{
    //store factorial from 0 to n
    double factorial[n + 1];
    factorial[0] = 1.0;
    for(int i = 1; i <= n; i ++)
    {
        factorial[i] = factorial[i - 1] * double(i);
    }

    //calculate polynomial controlling points
    vector<Point> polynomial_points;
    polynomial_points.clear();
    for(int j = 0; j <= n; j ++)
    {
            //calculate q_j
        Point q_j;
        for(int i = 0; i <= j; i ++)
        {
            //calculate q_ji
            double parameter = 1.0;
            parameter = parameter * factorial[n] / factorial[i] / factorial[n - j] / factorial[j - i];
            if((j - i) % 2 != 0)
            {
                parameter = parameter * -1;
            }
            Point q_ji = ControlPoints[i] * parameter;
            q_j = q_j + q_ji;
        }
        polynomial_points.push_back(q_j);
    }
    Polynomial polynomial = Polynomial(polynomial_points);
    return polynomial;
}

Bezier Polynomial::SwitchToBezier()
{
    vector<Point> bezier_points;
    bezier_points.clear();
    //TODO
    Bezier bezier = Bezier(bezier_points);
    return bezier;
}
