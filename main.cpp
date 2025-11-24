#include "utils/ODEsolver.hpp"
#include "utils/FVMsolver.hpp"
#include "equations/Eulers_equation.hpp"
#include "utils/BoundaryConditions.hpp"
#include "utils/InitialConditions.hpp"
#include <iostream>
#include <string>


// define data types
static constexpr std::size_t M = 3; // Euler's equations
using RealNumber = float;
using VectorS  = Vector< M, RealNumber >;
using DataType = std::vector< VectorS >;

// quick print function
void print_data( size_t N, DataType data )
{
  std::array< std::string, 3 > names{ "density", "momentum", "energy" };
  for ( size_t i = 0; i < N; i++)
  {
    std::cout << names[ i ] << std::endl;
    for ( size_t j = 0; j < data[ i ].size(); j++ )
    {
      std::cout << data[ i ][ j ];
    }
    std::cout << std::endl;
  }
}
// ----------------------------------
//            main
// ----------------------------------


int main()
{
  VectorS state_left;
  VectorS state_right;
  // fill the vectors
  // state[0] : density
  // state[1] : momentum
  // state[2] : energy

  RealNumber x0 = 0.5;

  RealNumber dx = 0.01;
  RealNumber time_step = 0.001;

  DataType data;

  // physical flux
  using Equations_ = EulersEquations< M, RealNumber >;
  Equations_ Equations{1.4};

  // choose numerical flux
  auto max_speed = Equations.max_speed();
  using NumericalFlux_ = Rusanov< M, RealNumber, decltype( max_speed ) >;

  // choose boundary conditions
  using BoundaryConditions_ = NeumannBC< M, RealNumber >;
  BoundaryConditions_ BoundaryConditions;

  // choose initial conditioins
  using InitialConditions_ = Sods_problem< M, RealNumber, DataType >;
  InitialConditions_ InitialConditions( x0, state_left, state_right );
  InitialConditions.impose( data, dx );

  // this should later be defined by the numerical flux used
  auto time_stepping_function = [] ()
  {
    return 0.012;
  };


  // FVMsolver intialization and definition
  FVMsolver< M, RealNumber, NumericalFlux_, Equations_, BoundaryConditions_ > fvm_solver( dx, Equations, BoundaryConditions );

  // ODE solver (time integration)
  ODEsolver< DataType , RealNumber, Euler, decltype( fvm_solver ), decltype( time_step_function )> solver(fvm_solver, time_step);



  // solve
  auto sol = solver.next_step(0.0, 0.1, data);
  
  return 0;
}
