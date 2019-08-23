#include "metric.h"

metric::metric::metric(int DIM,matrixFunction metFunc) : _dimensions(DIM),_metric(metFunc)
{
    
}

metric::metric::~metric()
{

}

metric::metric::metric(int DIM,matrixFunction metFunc,matrixFunction invMetFunc) : _dimensions(DIM),_metric(metFunc), _invMetric(invMetFunc)
{

}


metric::metric::metric(int DIM,matrixFunction metFunc,matrixFunction invMetFunc,derivativeMatrixFunction invDervMet) : _dimensions(DIM),_metric(metFunc), _invMetric(invMetFunc),_invDervMetric(invDervMet)
{

}

metric::metric::metric(functionVector diagonalMetric)
{
    
    
}