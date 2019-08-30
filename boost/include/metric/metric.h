#pragma once
#include <vector>
#include <memory>
#include <functional>
#include <iostream>
#include <boost/math/tools/numerical_differentiation.hpp>


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
        int Dimension;

    public:
        /// Metric constructor with the metric function provided
        metric(int DIM, matrixFunction metFunc);

        /// Metric constructor with metric and inverse metric provided
        metric(int DIM, matrixFunction metFunc,matrixFunction invMetFunc);

        /// Metric constructor with metric, inverse metric, and derivative of the metrix provided (last index is the partial derivative index)
        metric(int DIM, matrixFunction metFunc,matrixFunction invMetFunc,derivativeMatrixFunction invDervMet);


        /// Diagonal metric constructor
        metric(functionVector diagonalMetric );


        /// Diagonal metric constructor with derivatives passed in as dg^ii(x)/dx^k = invDerFunctions(i,k,x)
        metric(functionVector diagonalMetric , std::function<double(int, int, std::vector<double>)> invDerFunctions);


        ~metric();

        /// boolean value designating if the metric is diagonal
        bool isDiagonal = false;

        int getDimension() const { return Dimension; }

        std::vector<double> raise(const std::vector<double>& coVector) const;
        std::vector<double> lower(const std::vector<double>& contraVector) const;
        


    };
            
} // namespace metric
