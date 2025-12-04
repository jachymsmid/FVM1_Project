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
static constexpr float gamma_gas = 1.4;

// array of variable vectors
using Data = DataStorage< RealNumber, N, Num_Cells >;
// specific equations to be solved
using Equations_ = EulersEquations< RealNumber, Data, Vector< RealNumber, N >, gamma_gas >;

// spcific numerical flux to be used for spatial discretization
auto max_speed = Equations_::max_speed;

using NumericalFlux_ = Rusanov< N, RealNumber, Equations_ >;

// boundary conditions
using BoundaryConditions_ = NeumannBC< RealNumber, N, Data >;

// initial conditioins
using InitialConditions_ = Sods_problem< RealNumber, Data, Vector< RealNumber, N > >;

// numerical integration type
using NumericalIntegration_ = Euler< RealNumber, Data >;

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
  using FVMsolver_ = FVMsolver< RealNumber, NumericalFlux_, Equations_, BoundaryConditions_, Data >;

  // ODE solver (time integration)
  ODEsolver< RealNumber, Data, NumericalIntegration_, decltype( FVMsolver_::rhs ), decltype( FVMsolver_::time_step)> solver( FVMsolver_::rhs, FVMsolver_::time_step );

                                                                                                                                                                                                                                                                                                                                                                                                                                                        


  // solve
  auto sol = solver.next_step( 0.1, data);
  
  return 0;
}
