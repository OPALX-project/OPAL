#ifndef __SOCIAL_NETWORK_GRAPH__
#define __SOCIAL_NETWORK_GRAPH__

#include <set>
#include <cstdint>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

// Here, we use a simple mapping from processor ID's to 2D grid
// (determining neigborhood) assuming the total number of masters:
//
//    n_m = n_g * n_g
//
// This is fine, because if we have topology info the masters are numbered
// in ascending order (done in comm splitter) and we can map them from 1D to
// 2D by index calculations.

/**
 * \brief Modeling social graph (Cartesian neighbors plus additional random
 *        link) as underlaying master network.
 *
 * @see MasterNode
 *
 * Due to its nice rumor spreading properties this master network is mimicking
 * a social network graph (augmented grid) and solution states are distributed
 * in a rumor fashion.
 */
template < class TopoDiscoveryStrategy_t >
class SocialNetworkGraph : public TopoDiscoveryStrategy_t {

public:

    std::set<size_t> execute(size_t numMasters, size_t dimensions, size_t id,
                             int /*island_id*/) {

        numMasters_ = numMasters;
        dim_        = dimensions;
        myID_       = id;

        chooseRandomNeighbor();
        setNetworkNeighbors();
        return realNetworkNeighborPIDs_;
    }


private:

    boost::random::mt19937 gen_;
    double alpha_;

    size_t numMasters_;
    size_t dim_;
    size_t myID_;

    size_t randomNeighbor_;
    std::set<size_t> realNetworkNeighborPIDs_;

    void setNetworkNeighbors() {
        size_t m = static_cast<size_t>(sqrt(numMasters_));
        size_t north = myID_ + m % numMasters_;
        int64_t south = myID_ - m;
        if(south < 0) south += numMasters_;

        size_t east = myID_ + 1;
        if((myID_ + 1) % m == 0)
            east = myID_ - m + 1;
        size_t west = myID_ - 1;
        if(myID_ % m == 0)
            west = myID_ + m - 1;

        realNetworkNeighborPIDs_.insert(north);
        realNetworkNeighborPIDs_.insert(south);
        realNetworkNeighborPIDs_.insert(east);
        realNetworkNeighborPIDs_.insert(west);
    }


    //XXX: at the moment assumes square number of masters
    //void setNetworkNeighbors() {

        //size_t m = static_cast<size_t>(sqrt(numMasters_));

        //size_t north = myID_ + m;
        //if(north < numMasters_)
            //realNetworkNeighborPIDs_.insert(north);
        //size_t south = myID_ - m;
        //if(south >= 0)
            //realNetworkNeighborPIDs_.insert(south);

        //size_t east = myID_ + 1;
        //if((myID_ + 1) % m != 0)
            //realNetworkNeighborPIDs_.insert(east);
        //size_t west = myID_ - 1;
        //if(myID_ % m == 0)
            //realNetworkNeighborPIDs_.insert(west);
    //}


    double manhattenDistance(size_t from, size_t to) {

        size_t m = static_cast<size_t>(sqrt(numMasters_));
        int x_from = from / m;
        int y_from = from % m;
        int x_to = to / m;
        int y_to = to % m;

        return abs(x_from - x_to) + abs(y_from - y_to);
    }

    /// compute random neighbor using power law distribution with
    /// \f$ \alpha = 2\f$.
    void chooseRandomNeighbor() {

        std::vector<double> probabilities(numMasters_, 0.0);

        double sum = 0.0;
        for(size_t i = 0; i < numMasters_; i++) {
            if(i == myID_) continue;
            sum += std::pow(manhattenDistance(myID_, i), -alpha_);
        }

        for(size_t i = 0; i < numMasters_; i++) {

            if(i == myID_) continue;

            probabilities[i] =
                std::pow(manhattenDistance(myID_, i), -alpha_) / sum;
        }

        boost::random::discrete_distribution<>
            dist(probabilities.begin(), probabilities.end());
        randomNeighbor_ = dist(gen_);
    }

};

#endif
