#include "momentum_assembling.h"

using namespace std;
using namespace INMOST;

vector<double> LocalNVectorAssembling(const vector<vector<double>>& trnodes, int x_or_y)
{
    double x0 = trnodes[0][0];
    double y0 = trnodes[0][1];

    double x1 = trnodes[1][0];
    double y1 = trnodes[1][1];

    double x2 = trnodes[2][0];
    double y2 = trnodes[2][1];

    std::vector<double> coeffs0 = Solve3x3({{x0, x1, x2},
                                            {y0, y1, y2},
                                            {1.0, 1.0, 1.0}}, {1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = Solve3x3({{x0, x1, x2},
                                            {y0, y1, y2},
                                            {1.0, 1.0, 1.0}}, {0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = Solve3x3({{x0, x1, x2},
                                            {y0, y1, y2},
                                            {1.0, 1.0, 1.0}}, {0.0, 0.0, 1.0});
    
    ScalarFunction phi0(FuncType::linear, coeffs0, trnodes);
    ScalarFunction phi1(FuncType::linear, coeffs1, trnodes);
    ScalarFunction phi2(FuncType::linear, coeffs2, trnodes);

    double first, second, third;

    if (x_or_y == 0)
    {
        first = integral_over_triangle(trnodes, d_dx(phi0));
        second = integral_over_triangle(trnodes, d_dx(phi1));
        third = integral_over_triangle(trnodes, d_dx(phi2));
    }
    else
    {
        first = integral_over_triangle(trnodes, d_dy(phi0));
        second = integral_over_triangle(trnodes, d_dy(phi1));
        third = integral_over_triangle(trnodes, d_dy(phi2));
    }
    
    return {first, second, third};
}