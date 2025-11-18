#pragma once

#include <vector>


// -----------------------------------------
//          ODEsolver
// -----------------------------------------

template
<
  class DataType,
  class RealNumber,
  template< class, class, class > class Method,
  class RHS,
  class TimeStepType
>
class ODEsolver
{
public:
  // constructor - needs the rhs() and time_step()
  ODEsolver( const RHS &rhs_function_, const TimeStepType &time_step_function_ ) : rhs_function( rhs_function_ ), time_step_function( time_step_function_ ) {}

  // computes the next time_step
  DataType next_step(TimeStepType time_step_function, RealNumber time, const DataType &data_previous)
  {
    DataType data = data_previous; // this shoul be a copy constructor
    RealNumber time_step = time_step_function( data_previous ); // call time-stepping function defined by FVMsolver
    Method< DataType, RealNumber, TimeStepType >::step( rhs_function, time, time_step, data); // advences data one time_step into future
    return data;
  }

  // solve the whole equation, return array of data at different time steps
  std::vector< DataType > solve(){}

private:
    RHS rhs_function;
    TimeStepType time_step_function;
};

// ------------------ Euler ---------------------------
template
<
    class DataType,
    class RealNumber
>
struct Euler
{
  template<class RHS>
  static void step(RHS &rhs_function, RealNumber time, RealNumber time_step, DataType &data)
  {
    DataType k1 = rhs_function( time, data );
    for (size_t i = 0; i < data.size(); ++i)
      data[i] += time_step * k1[i];
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
  static void step(RHS &rhs_function, RealNumber time, TimeStepType &time_step, DataType &data)
  {
    DataType k1 = rhs_function( time, data );
    for (size_t i = 0; i < data.size(); ++i)
      data[i] += time_step * k1[i];
  }
};

