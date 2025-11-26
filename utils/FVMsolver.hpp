#pragma once

/* 
 * -----------------------------------
 *       FVMsolver struct 
 * -----------------------------------
 */

#include <cstdlib>
#include "GeneralVector.hpp"
#include "DataStorage.hpp"

template
<
    size_t Size,
    size_t Length,
    class RealNumber,
    template< size_t, class, class > class NumericalFlux,  // numerical flux
    class Equation,  // physical flux F(u)
    class BoundaryCondition
>
class FVMsolver {
public:

  using Data = DataStorage< RealNumber, Size, Length >;
    
    // constructor
    FVMsolver( RealNumber dx, Equation Flux, BoundaryCondition BC) : spatial_step_(dx), boundary_condition_(BC) {}

    void rhs( RealNumber t, const Data &data, Data &rhs_array )
    {

      Vector< RealNumber, Size > ghost_cell_left = boundary_condition_.left( data );
      Vector< RealNumber, Size > ghost_cell_right = boundary_condition_.right( data );

      Vector< RealNumber, Size > value_left, value_center, value_right; // all the values for one cell, so for eulers this is a vector of three values
      RealNumber interface_flux_l, interface_flux_r;

      for (std::size_t j = 0; j < Length/2; j++)
      {
        for (std::size_t k = 0; k < Size; k++)
        {
          value_left[k] = data[k][2*j];
          value_center[k] = data[k][2*j+1];
          value_right[k] = data[k][2*j+2];
        }

        for ( std::size_t k = 0; k < Size; k++)
        {
          interface_flux_[k][2*j] = NumericalFlux< Size, RealNumber, Equation >::numerical_flux( value_left, value_center );
          interface_flux_[k][2*j+1] = NumericalFlux< Size, RealNumber, Equation >::numerical_flux( value_center, value_right );
        }
      }

      for ( std::size_t i = 0; i < Size; i++ )
      {
        for ( std::size_t j = 0; j < Length - 1; j++ )
        {
          rhs_array[i][j+1] = - ( interface_flux_[i][j] - interface_flux_[i][j+1] ) / spatial_step_;
        }
      }
    }

private:
    RealNumber spatial_step_; // only for regular grids
    BoundaryCondition boundary_condition_;
    DataStorage< RealNumber, Size, Length - 1 > interface_flux_;
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
