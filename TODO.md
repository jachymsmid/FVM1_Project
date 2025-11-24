TODO:
   - [ ] use different DataType
       - right now im using an array of size 'n' of arrays of  size 3 ( for Euler's equations )
         which is pretty horrible
       - mby define the Mesh data type?
       - or use an array of size 3 of arrays of size 'n' - reuse code from poea1_project
       - or start rewriting using TNL
   - [x] implement zero gradient boundary condition
   - [ ] discretization in time
       - [x] forward Euler
       - [ ] Heune ( RK2 )
       - [ ] look into other solvers
   - [x] tie everything into objects
   - [ ] write the main body
   - [x] CMake
   - [ ] VTK
   - [ ] code documentation
   - [ ] implement conservative to primitive conversion and vice versa - where? Eulers_equations.hpp?
   - [ ] max wave speed for rusanov - define in FVMsolver
   - [ ] time-stepping function defined in FVMsolver by specific numerical flux

DONE:
