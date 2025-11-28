#pragma once

#include <cstddef>
template
<
  class RealNumber, 
  class DataStorage,
  class Vector
>
struct Sods_problem
{
    static void impose( DataStorage &data,
                        RealNumber spatial_step,
                        RealNumber x0,
                        Vector left,
                        Vector right )
    {
      RealNumber x; // coordinate
      for ( std::size_t i = 0; i < data.getLength(); i++ )
      {
        x = i * spatial_step;

        data( i ) = ( x < x0 ) ? left : right;
        
      }
    }
};
