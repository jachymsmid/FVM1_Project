#include "utils/ODEsolver.hpp"
#include "utils/FVMsolver.hpp"
#include "equations/Eulers_equation.hpp"
#include "utils/BoundaryConditions.hpp"
#include "utils/InitialConditions.hpp"

/* TODO:
 *    - [ ] implement zero gradient boundary condition
 *    - [ ] discretization in time
 *        - [ ] forward Euler
 *        - [ ] Heune ( RK2 )
 *        - [ ] look into other solvers
 *    - [ ] tie everything into objects
 *    - [ ] is returning a triple instead of an array a good idea?
 *    - [ ] write the main body
 *    - [ ] CMake
 *    - [ ] VTK
 */

// ----------------------------------
//            main
// ----------------------------------

using RealNumber = float;

int main()
{
  static constexpr std::size_t M = 3; // Euler equations

  using VectorS   = Vector< M, RealNumber >;
  using StateT = std::vector< VectorS >;

  RealNumber dx = 0.01;
  RealNumber time_step = 0.001;

  // physical flux
  using Equations = EulersEquations< M, RealNumber >;
  Equations physical_flux{1.4};

  // choose numerical flux
  using NumericalFlux = Rusanov< M, RealNumber >;

  // choose boundary conditions
  using BoundaryConditions = NeumannBC< M, RealNumber >;
  BoundaryConditions bc;

  // build spatial operator                                                                                                                                                                                                                                                                  
  auto rhs = FVMsolver< M, RealNumber, Rusanov, Equations, BoundaryConditions >(dx, physical_flux, bc);

  // ODE solver (time integration)
  ODEsolver< StateT, RealNumber, Euler, decltype(rhs), RealNumber > solver(rhs, time_step);

  VectorS data;


  // solve
  auto sol = solver.next_step(0.0, 0.1, data);
  
  return 0;
}
