#include "metric.h"


namespace metric
{
    using namespace metric;

    /*
        Geodesic driver which takess in a Metric, contravariant initial position (xi), contravariant initial momentum (pi), and a final time
        to get x^mu(tf) and p^mu(tf).  The return value is a pair with the first element is the position and the second is the momentum.
    */
    std::pair< std::vector<double>, std::vector<double> > geodesicDriver(const metric& Metric, const std::vector<double>& xi, 
                                                                        const std::vector<double>& pi , double tauf)
    {
        
    }

    

}