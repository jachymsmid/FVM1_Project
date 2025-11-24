#include "utils/ODEsolver.hpp"
#include "utils/GeneralVector.hpp"
#include "utils/DataStorage.hpp"
#include "utils/FVMsolver.hpp"
#include "equations/Eulers_equation.hpp"
#include "utils/BoundaryConditions.hpp"
#include "utils/InitialConditions.hpp"
#include <iostream>
#include <string>

// start using a naming convention!!

// define data types
static constexpr std::size_t N = 3; // Euler's equations
static constexpr std::size_t Num_Cells = 100;
using RealNumber = float;
static constexpr RealNumber gamma_gas = 1.4;

// vector of one variable
using VectorS = Vector< RealNumber, Num_Cells >;
// array of variable vectors
using Data = DataStorage< RealNumber, N, Num_Cells >;
// specific equations to be solved
using Equations_ = EulersEquations< N, Num_Cells, RealNumber, gamma_gas >;

// spcific numerical flux to be used for spatial discretization
auto max_speed = Equations_::max_speed();

using NumericalFlux_ = Rusanov< N, RealNumber, decltype( max_speed ) >;

// boundary conditions
using BoundaryConditions_ = NeumannBC< N, RealNumber >;

// initial conditioins
using InitialConditions_ = Sods_problem< N, RealNumber, Data >;


//// quick print function
//void print_data( size_t N, DataStorage data )
//{
//  std::array< std::string, 3 > names{ "density", "momentum", "energy" };
//  for ( size_t i = 0; i < N; i++)
//  {
//    std::cout << names[ i ] << std::endl;
//    for ( size_t j = 0; j < data[ i ].size(); j++ )
//    {
//      std::cout << data[ i ][ j ];
//    }
//    std::cout << std::endl;
//  }
//}



// ----------------------------------
//            main
// ----------------------------------

int main()
{
  Data data;

  // for sod's problem
  Vector< RealNumber, N > state_left;
  Vector< RealNumber, N > state_right;
  // index 0 = density ( rho )
  // index 1 = momentum ( rho * u )
  // index 2 = energy ( E )

  RealNumber x0 = 0.5; // position of the divide
  RealNumber dx = 0.01; // spatial step

  // impose initial conditions
  InitialConditions_::impose( data, dx, x0, state_left, state_right);


  // FVMsolver intialization and definition
  FVMsolver< N, RealNumber, NumericalFlux_, Equations_, BoundaryConditions_ > fvm_solver( dx );

  // ODE solver (time integration)
  ODEsolver< DataType , RealNumber, Euler, decltype( fvm_solver ), decltype( time_step_function )> solver(fvm_solver, time_step);



  // solve
  auto sol = solver.next_step(0.0, 0.1, data);
  
  return 0;
}
