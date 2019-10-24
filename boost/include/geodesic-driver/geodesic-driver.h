#include "metric.h"


namespace metric
{

    /*
        Geodesic driver which takess in a Metric, contravariant initial position x(0) = xi, contravariant initial velocity dx/dtau = vi, and a final time
        to get x^mu(tf) and p^mu(tf).  The return value is a pair with the first element is the position and the second is the momentum.
    */
    std::pair< contravariantVector, contravariantVector > geodesicDriver(const metric& Metric, const contravariantVector& xi, 
                                                                        const contravariantVector& vi , double tauf)
    {
        /*
            The following has to be used:

            https://www.boost.org/doc/libs/1_64_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/integrate_functions.html
        */
    }

    

}