/**
 * @file MapGenerator.h
 * This class is based on the paper "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons" (2012)
 * of Dr. Christian Baumgarten.
 * It has one template parameter that specifies the type of the variables and containers.
 *
 * @author Matthias Frey
 * @version 1.1
 */
#ifndef MAPGENERATOR_H
#define MAPGENERATOR_H

#include <cmath>
#include <stdexcept>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "matrix_vector_operation.h"

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FMatrix.h"

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

/// @brief This class generates the matrices for the one turn matrix of a cyclotron.
template<typename Value_type, typename Size_type, typename Series_type,
         typename Map_type, typename Hamiltonian_type, typename Space_charge_type>
class MapGenerator
{
    public:
        /// Type of variables
        typedef Value_type value_type;
        /// Type for specifying sizes
        typedef Size_type size_type;
        /// Type of truncated power series
        typedef Series_type series_type;
        /// Type of a map
        typedef Map_type map_type;
        /// Type of the Hamiltonian
        typedef Hamiltonian_type hamiltonian_type;
        /// Type of the Hamiltonian representing the space charge
        typedef Space_charge_type space_charge_type;
        
        /// Type for specifying matrices
        typedef boost::numeric::ublas::matrix<value_type> matrix_type;
        /// Type for specifying vectors
        typedef std::vector<value_type> vector_type;
        
        /// Initialize
        /*!
         * @param nMaps is the number of maps
         */
        MapGenerator(size_type);

        /// Generates a map based on the Hamiltonian for a given angle
        /*!
         * @param H represents the Hamiltonian
         * @param ds is the step size (angle dependent)
         * @param order is the truncation order of the Taylor series of the exponential function
         */
        matrix_type generateMap(const series_type&, value_type, size_type);
        
        /// Combine given maps
        /*!
         * Combines the space charge maps (for each angle one) and the cyclotron maps (for each angle one)
         * to the ont turn map, taking lists of maps
         * @param Mscs is a list of space charge maps (the higher the index, the higher the angle)
         * @param Mcycs is a list of cyclotron maps (the higher the index, the higher the angle)
         */
        void combine(std::vector<matrix_type>&, std::vector<matrix_type>&);
        
        /*!
         * Combine given container of maps.
         * @param maps to be combined.
         */
        matrix_type combine(std::vector<matrix_type>& maps);
        
        /// Returns the one turn map
        matrix_type getMap();
        
        /*!
         * Compute the radial and vertical tune from the map.
         * @param map from where to compute the tunes.
         * @returns the radial and vertical tunes (in this order)
         */
        std::pair<value_type, value_type> computeTunes(const matrix_type& map);
        
    private:
        /// Number of maps
        size_type nMaps_m;
        /// One-turn matrix
        matrix_type Mturn_m;
//         /// Stores the one turn transfer map
//         map_type M_m;
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------


template<typename Value_type,
         typename Size_type,
         typename Series_type,
         typename Map_type,
         typename Hamiltonian_type,
         typename Space_charge_type>
MapGenerator<Value_type,
             Size_type,
             Series_type,
             Map_type,
             Hamiltonian_type,
             Space_charge_type>::MapGenerator(size_type nMaps)
: nMaps_m(nMaps)
{}

template<typename Value_type,
         typename Size_type,
         typename Series_type,
         typename Map_type,
         typename Hamiltonian_type,
         typename Space_charge_type>
typename MapGenerator<Value_type,
                      Size_type,
                      Series_type,
                      Map_type,
                      Hamiltonian_type,
                      Space_charge_type>::matrix_type
MapGenerator<Value_type,
             Size_type,
             Series_type,
             Map_type,
             Hamiltonian_type,
             Space_charge_type>::generateMap/*cyc*/(const series_type& H,
                                                    value_type ds,
                                                    size_type order) {
    
    // expand
    map_type M = ExpMap(-H * ds, order);
    
    // get linear part
    FMatrix<value_type,6,6> matrix = M.linearTerms();
    
    matrix_type map = boost::numeric::ublas::zero_matrix<value_type>(6);
    
    for (size_type i = 0; i < 6; ++i) {
        for (size_type j = 0; j < 6; ++j) {
            map(i,j) = matrix[i][j];
        }
    }
    
    return map;
}

template<typename Value_type,
         typename Size_type,
         typename Series_type,
         typename Map_type,
         typename Hamiltonian_type,
         typename Space_charge_type>
void MapGenerator<Value_type,
                  Size_type,
                  Series_type,
                  Map_type,
                  Hamiltonian_type,
                  Space_charge_type>::combine(std::vector<matrix_type>& Mscs,
                                              std::vector<matrix_type>& Mcycs) {
    
    if (nMaps_m != Mscs.size() || nMaps_m != Mcycs.size())
        throw OpalException("MapGenerator::combine()", "Wrong vector dimensions.");

    Mturn_m = boost::numeric::ublas::identity_matrix<value_type>(6);
    
    for (size_type i = 0; i < nMaps_m; ++i)
        Mturn_m = matt_boost::gemmm<matrix_type>(Mscs[i],Mcycs[i],Mturn_m);
}

template<typename Value_type,
         typename Size_type,
         typename Series_type,
         typename Map_type,
         typename Hamiltonian_type,
         typename Space_charge_type>
typename MapGenerator<Value_type,
                  Size_type,
                  Series_type,
                  Map_type,
                  Hamiltonian_type,
                  Space_charge_type>::matrix_type
MapGenerator<Value_type,
                  Size_type,
                  Series_type,
                  Map_type,
                  Hamiltonian_type,
                  Space_charge_type>::combine(std::vector<matrix_type>& maps)
{
    matrix_type map = boost::numeric::ublas::identity_matrix<value_type>(6);
    
    for (std::size_t i = 0; i < maps.size(); ++i) {
        matrix_type tmp = prod(maps[i], map);
        map = tmp;
    }
    
    return map;
}



template<typename Value_type,
         typename Size_type,
         typename Series_type,
         typename Map_type,
         typename Hamiltonian_type,
         typename Space_charge_type>
typename MapGenerator<Value_type,
                      Size_type,
                      Series_type,
                      Map_type,
                      Hamiltonian_type,
                      Space_charge_type>::matrix_type
MapGenerator<Value_type,
             Size_type,
             Series_type,
             Map_type,
             Hamiltonian_type,
             Space_charge_type>::getMap() {
    return Mturn_m;
}


template<typename Value_type,
         typename Size_type,
         typename Series_type,
         typename Map_type,
         typename Hamiltonian_type,
         typename Space_charge_type>
         
std::pair<Value_type,
          Value_type> MapGenerator<Value_type,
                                   Size_type,
                                   Series_type,
                                   Map_type,
                                   Hamiltonian_type,
                                   Space_charge_type>::computeTunes(const matrix_type& map)
{
    /*
     * M = [ cos(mu) + alpha * sin(mu)  beta * sin(mu)
     *       -gamma * sin(mu)           cos(mu) - alpha * sin(mu)]
     * 
     * i = 0, 2, 4
     * --> cos(mu) = 0.5 * [ M(i, i) + M(i + 1, i + 1) ]
     */
    
    value_type arg = 0.0;
    
    // horizontal phase advance [rad]
    arg = 0.5 * ( map(0, 0) + map(1, 1) );
    
    if ( std::abs( arg ) > 1 )
        throw OpalException("MapGenerator::computeTunes()",
                            "Horizontal phase advance: Acos argument " + std::to_string(arg) + " out of range.");
    
    value_type mux = std::acos( arg );
    
    // vertical phase advance [rad]
    arg = 0.5 * ( map(2, 2) + map(3, 3) );
    
    if ( std::abs( arg ) > 1 )
        throw OpalException("MapGenerator::computeTunes()",
                            "Vertical phase advance: Acos argument " + std::to_string(arg) + " out of range.");
    
    value_type muz = std::acos( arg );
    
    
    
    return std::make_pair<value_type,
                          value_type>( mux * Physics::u_two_pi,
                                       muz * Physics::u_two_pi );
}

#endif
