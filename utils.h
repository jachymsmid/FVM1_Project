#pragma once

#include <vector>

using RealNumber = float;

extern RealNumber gamma_g;

struct Triple
{
  RealNumber data1;
  RealNumber data2;
  RealNumber data3;

  Triple operator+( Triple &other );
  Triple operator/( Triple &other );
  Triple operator*( Triple &other );
  Triple operator-( Triple &other );

  Triple operator+( RealNumber scalar );
  Triple operator/( RealNumber scalar );
  Triple operator*( RealNumber scalar );
  Triple operator-( RealNumber scalar );

  static Triple primitive_to_conserved( const Triple primitive );
  static Triple conserved_to_primitive( const Triple conserved );
};

Triple flux( Triple conserved );

class Mesh
{
  private:
    std::vector< RealNumber > data1;
    std::vector< RealNumber > data2;
    std::vector< RealNumber > data3;
    RealNumber spatial_step;
  public:
    size_t number_ghost_cells;
    size_t size;
    Triple getValues( size_t i );
    std::vector< RealNumber > getData1();
    std::vector< RealNumber > getData2();
    std::vector< RealNumber > getData3();
};

class NumericalSolver
{
  struct Rusanov
  {
    static Triple numerical_flux( Triple data_left, Triple data_right );
  };

  template< typename T >
  Triple NumericalFlux( Triple data_left, Triple data_right );
  #include "templates/NumericalFlux.tpp"

  static Triple rhs( Mesh &mesh );
};

class ODEsolver
{
  struct Euler
  {
    static void solve( Mesh &mesh );
  };

  struct Heune
  {
    static void solve( Mesh &mesh );
  };

  template< typename T >
  static void Solve_ODE( Mesh &mesh );
  #include "templates/Solve_ODE.tpp"

};

struct Sods_problem
{
  static void impose( Mesh &mesh );
};

template< typename T >
void InitialConditions( Mesh &mesh );
#include "templates/InitialConditions.tpp"

struct Zero_gradient
{
  static void impose( Mesh &mesh );
};

template< typename T >
void BoundaryConditions( Mesh &mesh );
#include "templates/BoundaryConditions.tpp"
