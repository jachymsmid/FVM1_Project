#include <array>
#include "../utils/GeneralVector.hpp"
#include "../utils/DataStorage.hpp"


template
<
  std::size_t Size,
  std::size_t Length,
  class RealNumber,
  RealNumber gamma_g
>
struct EulersEquations {

    // flux function 
    static void flux(const DataStorage< RealNumber, Size, Length > &data, DataStorage< RealNumber, Size, Length > &flux)
    {
      for ( std::size_t i = 0; i < Length; i++ )
      {
        RealNumber rho = data[0][i];
        RealNumber mom = data[1][i];
        RealNumber E   = data[2][i];

        // I could use the cons to prim here
        RealNumber v = mom / rho;
        RealNumber p = ( gamma_g - 1.0) * (E - 0.5 * rho * v * v);

        flux[0][i] = mom;
        flux[1][i] = mom * v + p;
        flux[2][i] = v * (E + p);
      }
    }

    //// prim to cons
    //RealNumber operator()( const Vector< Size, RealNumber> value_left, const Vector< Size, RealNumber > value_right )
    //{
    //  // the primitive to conserved conversion would come in handy here
    //}
    //// cons to prim
    //RealNumber operator()( const Vector< Size, RealNumber> value_left, const Vector< Size, RealNumber > value_right )
    //{
    //  // the primitive to conserved conversion would come in handy here
    //}

    // also could use the cons to prim here so do it!
    static RealNumber max_speed(){};
};
