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

    // TODO: move these vectors to their own classes and extend their operations
    typedef std::vector< std::function<double(std::vector<double>)> > functionVector;
    typedef std::vector<double> contravariantVector;
    
    /// A class inhereting from a vector double, used to have a different vector type
    class covariantVector : std::vector<double>
    {

    };

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

        /// Metric constructor with metric and inverse metric provided
        metric(int DIM, matrixFunction metFunc,matrixFunction invMetFunc);

        /// Metric constructor with metric, inverse metric, and derivative of the metrix provided (last index is the partial derivative index)
        metric(int DIM, matrixFunction metFunc,matrixFunction invMetFunc,derivativeMatrixFunction invDervMet);


        ~metric();

        /// If input metrix is diagonal, call this function for optimization
        void setDiagonal() { isDiagonal = true;};

        /// boolean value designating if the metric is diagonal
        bool isDiagonal = false;

        int getDimension() const { return Dimension; }
        
        double operator()(int i, int j, contravariantVector coords) const {return _metric(i,j,coords);};//FIXME: add range check here

        double g(int i,int j, contravariantVector coords) const {return _metric(i,j,coords);};//FIXME: add range check here
        double invg(int i, int j, contravariantVector coords) const {return _invMetric(i,j,coords);}; //FIXME: see above


    };
            
} // namespace metric
