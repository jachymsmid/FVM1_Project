#pragma once

#include <array>

template < std::size_t Size, class RealNumber >
using Vector = std::array< RealNumber, Size>;
