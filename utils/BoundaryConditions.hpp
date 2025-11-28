#pragma once

#include "GeneralVector.hpp"
#include "DataStorage.hpp"

template
<
  class RealNumber,
  std::size_t Size,
  std::size_t Length
>
struct NeumannBC {

    Vector< RealNumber, Size > left( const DataStorage< RealNumber, Size, Length > &data ) const 
    {
      Vector< RealNumber, Size > vector_left;

      for ( std::size_t i = 0; i < Size; i++ )
      {
        vector_left[i] = data[i][0];
      }

      return vector_left;

    }

    Vector< RealNumber, Size > right( const DataStorage<RealNumber, Size, Length >  &data ) const
    {
      Vector< RealNumber, Size > vector_right;

      for ( std::size_t i = 0; i < Size; i++ )
      {
        vector_right[i] = data[i][Size - 1];
      }

      return vector_right;

    }
};
