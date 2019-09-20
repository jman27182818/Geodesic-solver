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

    // TODO: move these vectors to their own classes
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
        /// raise a covariant vector
        contravariantVector raise(const covariantVector& coVector) const;

        /// lower a contravariantVector
        covariantVector lower(const contravariantVector& contraVector) const;

        friend contravariantVector operator*(const metric& Metric,const covariantVector& vector);
        friend contravariantVector operator*(const covariantVector& vector,const metric& Metric);


        friend covariantVector operator*(const metric& Metric,const contravariantVector& vector);
        friend covariantVector operator*(const contravariantVector& vector,const metric& Metric);
        
        double operator()(int i, int j, contravariantVector coords) const {return _metric(i,j,coords);};//FIXME: add range check here

        double g(int i,int j, contravariantVector coords) const {return _metric(i,j,coords);};//FIXME: add range check here
        double invg(int i, int j, contravariantVector coords) const {return _invMetric(i,j,coords);}; //FIXME: see above


    };
    //TODO default minkowski constructor
            
} // namespace metric
