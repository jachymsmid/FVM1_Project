#pragma once

#include "GeneralVector.hpp"
#include "DataStorage.hpp"

template
<
  std::size_t Size,
  class RealNumber, 
  class DataType
>
struct Sods_problem
{
    static void impose( DataType &data,
                        RealNumber spatial_step,
                        RealNumber x0,
                        Vector< RealNumber, Size > left,
                        Vector< RealNumber, Size > right )
    {
      RealNumber x; // coordinate
      for ( std::size_t i = 0; i < Size; i++ )
      {
        for ( std::size_t j = 0; j < data[0].size(); j++ )
        {
          x = i * spatial_step;
          if ( x < x0 )
          {
            data[ i ][ j ] = left[ i ];
          }
          else
          {
            data[ i ][ j ] = right[ i ];
          }
        }
      }
    }
};
