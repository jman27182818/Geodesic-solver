#include "metric.h"


metric::metric::metric(int DIM,matrixFunction metFunc) : Dimension(DIM),_metric(metFunc)
{
    if(DIM > 6)
    {
        std::cout << "Unable to invert for dimensions less than six";
        throw;
    }
    switch (DIM)
    {
    case 2:
        _invMetric = [metFunc](int i, int j, std::vector<double> coords) 
        {
            auto factor = 1.0 / (-1.0* metFunc(0,1,coords) * metFunc(0,1,coords) + metFunc(0,0,coords) * metFunc(1,1,coords));
            if(i == j)
            {
                return metFunc(i+1 % 2, j + 1 % 2 , coords) * factor;
            }
            return -1.0 * metFunc(i,j,coords) * factor;
        };
        break;
    case 3:
        _invMetric = [metFunc](int i, int j, std::vector<double> coords)
        {
            
            return 0.0;
        };
        break;    
    default:
        std::cout << "invalid dimension" << std::endl; // FIXME: better error handling
        exit(-1);
        break;
    }
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


metric::metric::metric(functionVector diagonalMetric , std::function<double(int, int, std::vector<double>)> invDerFunctions)
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
    _invDervMetric = [invDerFunctions](int i, int j, int k, std::vector<double> coords)
    {
        if(i != j) return 0.0;
        return invDerFunctions(i,k,coords);
    };

    //TODO: overload the operators for the metric function 
    //TODO: make a struct for covariant vectors

}
