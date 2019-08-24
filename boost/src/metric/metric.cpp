#include "metric.h"


metric::metric::metric(int DIM,matrixFunction metFunc) : Dimension(DIM),_metric(metFunc)
{
    if(DIM > 6)
    {
        std::cout << "Unable to invert for dimensions less than six";
        throw;
    }
    //TODO: all the inversions haha
}


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

metric::metric::metric(functionVector diagonalMetric)
{
    isDiagonal = true;
    Dimension = diagonalMetric.size();
    _metric = [diagonalMetric](int i, int j, std::vector<double> coords) 
    { 
        if(i != j) return 0.0;        
        return diagonalMetric[i](coords);
    };
    _invMetric = [diagonalMetric](int i, int j, std::vector<double> coords) 
    {
        if(i != j) return 0.0;
        return 1.0 / (diagonalMetric[i](coords)); //FIXME: divide by zero check
    };
    _invDervMetric = [diagonalMetric](int i, int j, int k, std::vector<double> coords)
    {
        if(i != j) return 0.0;
        return -1.0 * (1.0 / diagonalMetric[i](coords)) * (1.0 / diagonalMetric[i](coords)) * 
        boost::math::tools::finite_difference_derivative([diagonalMetric,i,k,coords](double xk) {
            auto newCoords = coords;
            newCoords[k] = xk;
            return diagonalMetric[i](newCoords);
        }, coords[k]);
    };
    
}