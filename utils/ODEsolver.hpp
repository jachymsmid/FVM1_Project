#pragma once

// -----------------------------------------
//          ODEsolver
// -----------------------------------------

#include <cstddef>
template
<
  class RealNumber,
  class DataStorage,
  class Method,
  class RHS_function_type,
  class TimeStep_function_type
>
class ODEsolver
{
public:
  // constructor - needs the rhs() and time_step()
  ODEsolver( RHS_function_type &rhs_function_, TimeStep_function_type &time_step_function_ ) : rhs_function( rhs_function_ ), time_step_function( time_step_function_ ) {}

  // computes the next time_step
  DataStorage next_step(RealNumber time, const DataStorage &data_previous)
  {
    DataStorage data = data_previous; // this shoul be a copy constructor
    RealNumber time_step = time_step_function( data_previous ); // call time-stepping function defined by FVMsolver
    Method::step( rhs_function, time_step, data); // advences data one time_step into future
    return data;
  }

private:
    RHS_function_type rhs_function;
    TimeStep_function_type time_step_function;
};

// ------------------ Euler ---------------------------
template
<
    class RealNumber,
    class DataStorage
>
struct Euler
{
  template<class RHS>
  static void step(RHS &rhs_function, RealNumber time_step, DataStorage &data)
  {
    DataStorage rhs = rhs_function( data );
    for (std::size_t i = 0; i < data.getSize(); i++ )
    {
      for ( std::size_t j = 0; j < data.getLength(); j++ )
      {
        data(i,j) += time_step * rhs(i,j);
      }
    }
  }
};

// ------------------- Heune ---------------------------
// !! not completed !!
template
<
    class DataType,
    class RealNumber,
    class TimeStepType
>
struct Heune 
{
  template<class RHS>
  static void step(RHS &rhs_function, TimeStepType &time_step, DataType &data)
  {
  }
};

