#include "utils/ODEsolver.hpp"
#include "utils/FVMsolver.hpp"
#include "equations/Eulers_equation.hpp"
#include "utils/BoundaryConditions.hpp"
#include "utils/InitialConditions.hpp"
#include <iostream>
#include <string>

/* TODO:
 *    - [ ] use different DataType
 *        - right now im using an array of size 'n' of arrays of  size 3 ( for Euler's equations )
 *          which is pretty horrible
 *        - mby define the Mesh data type?
 *        - or use an array of size 3 of arrays of size 'n' - reuse code from poea1_project
 *        - or start rewriting using TNL
 *    - [x] implement zero gradient boundary condition
 *    - [ ] discretization in time
 *        - [x] forward Euler
 *        - [ ] Heune ( RK2 )
 *        - [ ] look into other solvers
 *    - [x] tie everything into objects
 *    - [ ] write the main body
 *    - [x] CMake
 *    - [ ] VTK
 *    - [ ] code documentation
 *    - [ ] implement conservative to primitive conversion and vice versa - where? Eulers_equations.hpp?
 *    - [ ] max wave speed for rusanov - define in FVMsolver
 *    - [ ] time-stepping function defined in FVMsolver by specific numerical flux
 */

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
  using Equations = EulersEquations< M, RealNumber >;
  Equations physical_flux{1.4};

  // choose numerical flux
  auto max_speed = Equations.max_speed();
  using NumericalFlux = Rusanov< M, RealNumber, max_speed >;

  // choose boundary conditions
  using BoundaryConditions = NeumannBC< M, RealNumber >;
  BoundaryConditions bc;

  // choose initial conditioins
  using InitialConditions = Sods_problem< M, RealNumber, DataType >;
  InitialConditions ic( x0, state_left, state_right );
  ic.impose( data, dx );

  // this should later be defined by the numerical flux used
  RealNumber time_stepping_function = [] ()
  {
    return 0.012;
  }


  // build spatial operator                                                                                                                                                                                                                                                                  
  auto rhs = FVMsolver< M, RealNumber, Rusanov, Equations, BoundaryConditions >(dx, physical_flux, bc);

  // ODE solver (time integration)
  ODEsolver< DataType , RealNumber, Euler, decltype(rhs), decltype( time_step_function )> solver(rhs, time_step);



  // solve
  auto sol = solver.next_step(0.0, 0.1, data);
  
  return 0;
}
