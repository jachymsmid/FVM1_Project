#pragma once

/* 
 * -----------------------------------
 *       NumericalSolver struct 
 * -----------------------------------
 */

#include <vector>

template
<
  class Array,
  class RealNumber,
  class NumericalScheme
>
class FVMsolver
{
public:

  // Compute RHS for system: vector at each cell
  std::vector< Array > rhs( const std::vector< Array > &data, RealNumber spatial_step )
  {
    const size_t N = data.size();
    std::vector< Array > rhs( N );


    for (size_t i = 1; i < N - 1; ++i) {
      Array flux_left = NumericalScheme::numerical_flux( data[ i ], data[ i + 1 ] );
      Array flux_right = NumericalScheme::numerical_flux( data[ i - 1 ], data[ i ] );
      Array r = flux_plus;

      for (size_t k = 0; k < r.size(); ++k)
      {
        r[k] = -(flux_plus[k] - flux_minus[k]) / spatial_step;
      }
      rhs[i] = r;
    }
  return rhs;
  }
};

template
<
  class RealNumber,
  class Array
>
struct Rusanov
{
  static Array numerical_flux( const Array &value_left, const Array &value_right )
  {
  }

  static RealNumber time_step( const Array &data )
  {
  }
};
