#include "metric.h"



metric::metric::~metric()
{

}



metric::metric::metric(int DIM,matrixFunction metFunc,matrixFunction invMetFunc) : Dimension(DIM),_metric(metFunc), _invMetric(invMetFunc)
{
    _invDervMetric = [invMetFunc](int i, int j, int k, std::vector<double> coords)
    {
        return boost::math::tools::finite_difference_derivative([invMetFunc,i,j,k,coords](double xk) {
            auto newCoords = coords;
            newCoords[k] = xk;
            return invMetFunc(i,j,newCoords);
        }, coords[k]);
    };
}


metric::metric::metric(int DIM,matrixFunction metFunc,matrixFunction invMetFunc,derivativeMatrixFunction invDervMet) : Dimension(DIM),_metric(metFunc), _invMetric(invMetFunc),_invDervMetric(invDervMet)
{

}



