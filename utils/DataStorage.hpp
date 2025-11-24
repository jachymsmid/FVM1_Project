#include "GeneralVector.hpp"
template < class RealNumber, std::size_t Size, std::size_t Length >
using DataStorage = std::array< Vector< RealNumber, Length >, Size >;
