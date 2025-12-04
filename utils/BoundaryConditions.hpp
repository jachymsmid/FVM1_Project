#pragma once

#include "GeneralVector.hpp"

template
<
  class RealNumber,
  std::size_t Size,
  class DataStorage
>
struct NeumannBC {

    Vector< RealNumber, Size > left( const DataStorage &data ) const 
    {
      Vector< RealNumber, Size > vector_left;

      for ( std::size_t i = 0; i < data.getSize(); i++ )
      {
        vector_left[i] = data[i][0];
      }

      return vector_left;

    }

    Vector< RealNumber, Size > right( const DataStorage &data ) const
    {
      Vector< RealNumber, Size > vector_right;

      for ( std::size_t i = 0; i < data.getSize(); i++ )
      {
        vector_right[i] = data[i][data.getSize() - 1];
      }

      return vector_right;

    }
};
