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
  class RealNumber 
>
struct Sods_problem
{
    using VectorS = Vector< Size, RealNumber >;

    RealNumber x0;  // interface position
    VectorS left;
    VectorS right;

    // constructor
    Sods_problem( RealNumber x0, VectorS L, VectorS R)
        : x0(x0), left(L), right(R) {}

    VectorS operator()( RealNumber x) const {
        return (x < x0 ? left : right);
    }
};
