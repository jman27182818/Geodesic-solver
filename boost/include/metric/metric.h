#pragma once
#include <vector>
#include <memory>
#include <functional>


namespace metric
{
    typedef std::function<double(int , int,std::vector<double>)> matrixFunction;
    typedef std::function<double(int, int, int, std::vector<double>)> derivativeMatrixFunction;
    typedef std::vector< std::function<double(std::vector<double>)> > functionVector;

    class metric
    {
    private:
        /// Metrix
        matrixFunction _metric;
        /// Inverse metric 
        matrixFunction _invMetric;
        /// Derivative of the metrix matrix
        derivativeMatrixFunction _invDervMetric;
        /// Total dimension of the system
        int _dimensions;
        /// boolean value designating if the metric is diagonal
        bool isDiagonal = false;

    public:
        /// Metric constructor with the metric function provided
        metric(int DIM, matrixFunction metFunc);

        /// Metric constructor with metric and inverse metric provided
        metric(int DIM, matrixFunction metFunc,matrixFunction invMetFunc);

        /// Metric constructor with metric, inverse metric, and derivative of the metrix provided (last index is the partial derivative index)
        metric(int DIM, matrixFunction metFunc,matrixFunction invMetFunc,derivativeMatrixFunction invDervMet);


        /// Diagonal metric constructor
        metric(functionVector diagonalMetric );

        ~metric();
    };
            
} // namespace metric
