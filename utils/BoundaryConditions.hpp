#pragma once

#include <array>
#include <vector>

template
<
  std::size_t Size,
  class RealNumber
>
using Vector = std::array< RealNumber, Size>;

template
<
  std::size_t Size,
  class RealNumber 
>
struct NeumannBC {
    using VectorS = Vector< Size, RealNumber >;
    using DataType = std::vector< VectorS >;

    VectorS left( const DataType &data ) const {
        return data[0];
    }
    VectorS right( const DataType &data ) const {
        return data[data.size() - 1];
    }
};
