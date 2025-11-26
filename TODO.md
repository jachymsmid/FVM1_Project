TODO:
   - [ ] discretization in time
       - [x] forward Euler
       - [ ] Heune ( RK2 )
       - [ ] some RK4
       - [ ] Svata's scheme
   - [ ] VTK, vector.dump
   - [ ] code documentation
   - [ ] implement conservative to primitive conversion and vice versa - where? Eulers_equations.hpp?
   - [ ] max wave speed for rusanov - define in FVMsolver, or lambda in rusanov?
   - [ ] time-stepping function defined in FVMsolver by specific numerical flux
   - [ ] paralellize the spatial for loops - mby in the parallel fvm project
   - [ ] use static functions instead of operators ( i just dont like them )
   - [ ] write tests ( sod, gau, ... )
        - [ ] write script for exact solutions

DONE:
   - [x] use different DataType
   - [x] implement zero gradient boundary condition
   - [x] CMake
   - [x] tie everything into objects
