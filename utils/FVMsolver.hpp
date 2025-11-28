#pragma once

/* 
 * -----------------------------------
 *       FVMsolver struct 
 * -----------------------------------
 */

#include <cstdlib>
#include "GeneralVector.hpp"

template
<
    class RealNumber,
    template< size_t, class, class > class NumericalFlux,  // numerical flux
    class Equation,  // physical flux F(u)
    class BoundaryCondition,
    class DataStorage
>
class FVMsolver {
public:

    
    // constructor
    FVMsolver( RealNumber dx, Equation Flux, BoundaryCondition BC) : spatial_step_(dx), boundary_condition_(BC) {}

    void rhs( const DataStorage &data, DataStorage &rhs_array )
    {

      Vector< RealNumber, data.getSize() > ghost_cell_left = boundary_condition_.left( data );
      Vector< RealNumber, data.getSize() > ghost_cell_right = boundary_condition_.right( data );

      Vector< RealNumber, data.getSize() > value_left, value_center, value_right; // all the values for one cell, so for eulers this is a vector of three values

      for (std::size_t j = 0; j < data.getLength()/2; j++)
      {
        // get the data for the central, left and right cell
        value_left = ( j == 0 ) ? ghost_cell_left : data(2*j);
        value_center = data(2*j+1);
        value_right = ( j == data.getLength()/2 - 1 ) ? ghost_cell_right : data(2*j+2);

        interface_flux_(2*j) = NumericalFlux< data.getSize(), RealNumber, Equation >::numerical_flux( value_left, value_center );
        interface_flux_(2*j+1) = NumericalFlux< data.getSize(), RealNumber, Equation >::numerical_flux( value_center, value_right );
      }

      for ( std::size_t i = 0; i < data.getSize(); i++ )
      {
        for ( std::size_t j = 0; j < data.getLength(); j++ )
        {
          rhs_array(i,j+1) = - ( interface_flux_(i,j) - interface_flux_(i,j+1) ) / spatial_step_;
        }
      }
    }

private:
    RealNumber spatial_step_; // only for regular grids
    BoundaryCondition boundary_condition_;
    DataStorage interface_flux_;
};


template
<
  size_t Size,
  class RealNumber,
  class Equation
>
struct Rusanov {

  static Vector< RealNumber, Size > numerical_flux(const Vector< RealNumber, Size > &value_left, const Vector< RealNumber, Size > &value_right)
  {
    Vector< RealNumber, Size > flux_left = Equation::flux( value_left );
    Vector< RealNumber, Size > flux_right = Equation::flux( value_left );

    // user must define max eigenvalue estimate - depend on the system of equations
    RealNumber a = Equation::max_wave_speed( value_left , value_right );
    // how to pass this function?

    Vector< RealNumber, Size > flux;

    for ( size_t i = 0; i < Size; i++ )
    {
      flux[i] = 0.5*( flux_left[i] + flux_right[i]) - 0.5*a*( value_right[i] - value_left[i] );
    }

    return flux;
  }
};
