#pragma once

/* 
 * -----------------------------------
 *       NumericalSolver struct 
 * -----------------------------------
 */

#include <cstdlib>
#include <vector>

template
<
  size_t Size,
  class RealNumber
>
using Vector = std::array< RealNumber, Size >;



template
<
    size_t Size,
    class RealNumber,
    template< size_t, class > class NumericalFlux,  // numerical flux
    class FluxFunction,  // physical flux F(u)
    class BoundaryCondition
>
class FVMsolver {
public:

    using VectorS   = Vector< Size, RealNumber >;
    using DataType = std::vector< VectorS >;

    // constructor
    FVMsolver( RealNumber dx, FluxFunction Flux, BoundaryCondition BC) : dx_(dx), Flux_(Flux), BC_(BC) {}

    DataType operator()( RealNumber t, const DataType &data ) const
    {
      std::size_t N = data.size();
      DataType rhs(N);

      VectorS ghost_cell_left = BC_.left( data );
      VectorS ghost_cell_right = BC_.right( data );

      for (std::size_t i = 0; i < N; ++i)
      {
        const VectorS &value_left = data[i];
        const VectorS &value_right = ( i+1 < N ) ? data[i+1] : ghost_cell_right;

        VectorS flux_right = NumericalFlux< Size, RealNumber >::numerical_flux( Flux_, value_left, value_right );

        const VectorS &value_left_2 = (i > 0) ? data[i-1] : ghost_cell_left;

        VectorS flux_left = NumericalFlux< Size, RealNumber >::numerical_flux( Flux_, value_left_2, value_left );

        for (std::size_t k = 0; k < Size; ++k)
        {
          rhs[i][k] = -( flux_right[k] - flux_left[k] ) / dx_;
        }
      }
      return rhs;
    }

private:
    RealNumber dx_;
    FluxFunction Flux_;
    BoundaryCondition BC_;
};


template< size_t Size, class RealNumber >
struct Rusanov {

    template< class FluxFunction >
    static Vector< Size, RealNumber > numerical_flux(const FluxFunction &Flux, const Vector< Size, RealNumber > &value_left, const Vector< Size, RealNumber > &value_right)
    {
        Vector< Size, RealNumber > flux_left = Flux( value_left );
        Vector< Size, RealNumber > flux_right = Flux( value_left );

        // user must define max eigenvalue estimate - depend on the system of equations
        RealNumber a = max_wave_speed( value_left , value_right ); // user-defined function

        Vector< Size, RealNumber > flux;
        for ( size_t i = 0; i < Size; ++i )
            flux[i] = 0.5*( flux_left[i] + flux_right[i]) - 0.5*a*( value_right[i] - value_left[i] );

        return flux;
    }
};
