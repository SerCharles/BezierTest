#include <vector>
#include <cmath>
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
    double dist(const Point& b) const 
    {
        double result = sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y));
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
    Point GetPlace(double t);
    Point DeCasteljau(double t);

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
    Point GetPlace(double t);
};

/*
Switch a bezier curve to a polynomial curve

Returns:
    polynomial [Polynomial]: [the polynomial curve]
*/
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

/*
Switch a polynomial curve to a bezier curve

Returns:
    bezier [Bezier]: [the bezier curve]
*/
Bezier Polynomial::SwitchToBezier()
{
    //store factorial from 0 to n
    double factorial[n + 1];
    factorial[0] = 1.0;
    for(int i = 1; i <= n; i ++)
    {
        factorial[i] = factorial[i - 1] * double(i);
    }

    vector<Point> bezier_points;
    bezier_points.clear();
    //initialize the matrix
    double matrix[n + 1][n + 1];
    for(int j = 0; j <= n; j ++)
    {
        for(int i = 0; i <= n; i ++)
        {
            if(i > j)
            {
                matrix[j][i] = 0;
                continue;
            }
            double parameter = 1.0;
            parameter = parameter * factorial[n] / factorial[i] / factorial[n - j] / factorial[j - i];
            if((j - i) % 2 != 0)
            {
                parameter = parameter * -1;
            }
            matrix[j][i] = parameter;
        }
    }

    //solve the linear equation
    for(int i = 0; i <= n; i ++)
    {
        double left_x = 0.0;
        double left_y = 0.0;
        double right_x = ControlPoints[i].x; 
        double right_y = ControlPoints[i].y;
        double parameter = matrix[i][i];
        for(int j = 0; j < i; j ++)
        {
            double new_left_x = matrix[i][j] * bezier_points[j].x;
            double new_left_y = matrix[i][j] * bezier_points[j].y;
            left_x = left_x + new_left_x; 
            left_y = left_y + new_left_y;
        }
        double remain_x = right_x - left_x;
        double remain_y = right_y - left_y;
        double result_x = remain_x / parameter;
        double result_y = remain_y / parameter;
        Point new_point = Point(result_x, result_y);
        bezier_points.push_back(new_point);
    }


    Bezier bezier = Bezier(bezier_points);
    return bezier;
}


/*
Get the place of the point on the polynomial curve given t

Args:
    t [double]: [the parameter t, which is in [0, 1]]

Returns:
    place [Point]: [the place of the result point on the curve given t]
*/
Point Polynomial::GetPlace(double t)
{
    double base = 1.0;
    Point place = Point();
    for(int i = 0; i <= n; i ++)
    {
        Point new_place = ControlPoints[i] * base;
        place = place + new_place;
        base = base * t;
    }
    return place;
}

/*
Get the place of the point on the bezier curve given t

Args:
    t [double]: [the parameter t, which is in [0, 1]]

Returns:
    place [Point]: [the place of the result point on the curve given t]
*/
Point Bezier::GetPlace(double t)
{
    //store factorial from 0 to n
    double factorial[n + 1];
    factorial[0] = 1.0;
    for(int i = 1; i <= n; i ++)
    {
        factorial[i] = factorial[i - 1] * double(i);
    }

    //store the t power from 0 to n
    double t_power[n + 1];
    t_power[0] = 1.0;
    for(int i = 1; i <= n; i ++)
    {
        t_power[i] = t_power[i - 1] * t;
    }

    //store the (1 - t) power from 0 to n
    double mt_power[n + 1];
    mt_power[0] = 1.0;
    for(int i = 1; i <= n; i ++)
    {
        mt_power[i] = mt_power[i - 1] * (1 - t);
    }

    Point place = Point();
    for(int i = 0; i <= n; i ++)
    {
        double bin = factorial[n] / factorial[i] / factorial[n - i] * t_power[i] * mt_power[n - i];
        Point new_place = ControlPoints[i] * bin;
        place = place + new_place;
    }
    return place;
}

/*
Get the place of the point on the bezier curve given t, using de Casteljau Algorithm

Args:
    t [double]: [the parameter t, which is in [0, 1]]

Returns:
    place [Point]: [the place of the result point on the curve given t]
*/
Point Bezier::DeCasteljau(double t)
{
    Point result_matrix[n + 1][n + 1];
    for(int i = 0; i <= n; i ++)
    {
        result_matrix[0][i] = ControlPoints[i];
    }
    for(int i = 1; i <= n; i ++)
    {
        for(int j = 0; j <= n - i; j ++)
        {
            result_matrix[i][j] = result_matrix[i - 1][j] * (1 - t) + result_matrix[i - 1][j + 1] * t;
        }
    }
    return result_matrix[n][0];
}