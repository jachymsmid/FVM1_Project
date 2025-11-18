#pragma once

#include <array>
template
<
  std::size_t Size,
  class RealNumber
>
using Vector = std::array< RealNumber, Size >;

template
<
  std::size_t Size,
  class RealNumber, 
  class DataType
>
struct Sods_problem
{
    using VectorS = Vector< Size, RealNumber >;

    RealNumber x0;  // position of the diaphragm 
    VectorS left;
    VectorS right;

    // constructor
    Sods_problem( RealNumber x0, VectorS L, VectorS R)
        : x0(x0), left(L), right(R) {}

    void impose( DataType data, RealNumber spatial_step )
    {
      RealNumber x; // coordinate
      for ( std::size_t i = 0; i < Size; i++ )
      {
        for ( std::size_t j = 0; j < Size; j++ )
        {
          x = i * spatial_step;
          if ( x < x0 )
          {
            data[ i ][ j ] = left[ j ];
          }
          else
          {
            data[ i ][ j ] = right[ j ];
          }
          // is ternary operator here any better?
        }
      }

    }
};
