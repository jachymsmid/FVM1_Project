#include "utils.h"
#include <iostream>

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
  // testing
  Mesh< std::vector< RealNumber >, RealNumber > mesh;
  mesh.write_data();
  std::cout << mesh.getSpatialStep() << std::endl;

  return 0;
}
