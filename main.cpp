#include <iostream>
#include <cmath>
#include <vector>

using RealNumber = float;

const RealNumber gamma_g = 1.4;

struct SimulationInfo
{
  size_t number_of_cells;
  size_t domain_length;
  RealNumber spatial_step;
  RealNumber time_step;
  RealNumber end_time;
  RealNumber advection_speed;
};

struct Conservative
{
  RealNumber density;
  RealNumber momentum;
  RealNumber energy;
};

struct Primitive
{
  RealNumber density;
  RealNumber speed;
  RealNumber pressure;
};

Conservative prim_to_cons( const Primitive &primitive )
{
  Conservative conservative;
  conservative.density = primitive.density;
  conservative.momentum = primitive.density * primitive.speed;
  conservative.energy = primitive.pressure / ( gamma_g - 1.0 ) + 0.5 * primitive.density * primitive.speed * primitive.speed;

  return conservative;
}

Primitive cons_to_prim( const Conservative &conservative )
{
  Primitive primitive;



  return primitive;
}

class Mesh
{
  std::vector< RealNumber > x_cord;
  std::vector< RealNumber > density;
  std::vector< RealNumber > momentum;
  std::vector< RealNumber > energy;
  size_t size;

public:
  
  // constructor
  Mesh() = default;

  // getters
  const size_t getSize() const
  {
    return size;
  }

  const Conservative getValues( size_t i ) const
  {
    return {density[i], momentum[i], energy[i]};
  }

  const std::vector< RealNumber > getCoordinates() const
  {
    return x_cord;
  }

  const RealNumber getCoordinates( size_t i ) const
  {
    return x_cord[i];
  }

  void constructGrid( SimulationInfo sim_info )
  {
    for ( size_t i = 0; i <  size; i++ )
    {
      x_cord[i] = i * RealNumber(sim_info.domain_length) / ( RealNumber(sim_info.number_of_cells) - 1 );
    }
  }

  void writeData() {}


};

// -----------------------------------
//        numerical solvers
// -----------------------------------

Conservative flux( Conservative data_c )
{
  Primitive data_p = cons_to_prim( data_c );
  Conservative flux;

  flux.density = data_c.momentum;
  flux.momentum = data_c.momentum * data_p.speed + data_p.pressure;
  flux.energy = ( data_c.density * data_c.energy + data_p.pressure ) * data_p.speed;

  return flux;
}

RealNumber max_signal_speed(const Conservative &U_L, const Conservative &U_R)
{
  Primitive U_Lp = cons_to_prim( U_L );
  Primitive U_Rp = cons_to_prim( U_R );

  RealNumber signal_speed_l = fabs( U_Lp.speed ) + sqrt( gamma_g * U_Lp.pressure / U_Lp.density );
  RealNumber signal_speed_r = fabs( U_Rp.speed ) + sqrt( gamma_g * U_Rp.pressure / U_Rp.density );

  return fmax( signal_speed_l, signal_speed_r );
}

RealNumber max_global_speed( const Mesh &mesh )
{
  RealNumber max_speed;
  for ( size_t i = 0; i < mesh.getSize(); i++ )
  {
    Primitive prim = cons_to_prim( mesh.getValues( i ) );
    RealNumber c = sqrt( gamma_g * prim.pressure / prim.density );
    max_speed = fmax( max_speed, c );
  }
  return max_speed;
}

struct Rusanov
{
  Conservative numerical_flux( const Conservative &U_L, const Conservative &U_R )
  {
    Conservative flux_l = flux( U_L );
    Conservative flux_r = flux( U_R );
    double alpha = max_signal_speed( U_L, U_R );
    Conservative numerical_flux;
    numerical_flux.density = 0.5 * ( flux_l.density + flux_r.density ) - 0.5 * alpha * ( U_R.density - U_L.density );
    numerical_flux.momentum = 0.5 * ( flux_l.momentum + flux_r.momentum ) - 0.5 * alpha * ( U_R.momentum - U_L.momentum );
    numerical_flux.energy = 0.5 * ( flux_l.energy + flux_r.energy ) - 0.5 * alpha * ( U_R.energy - U_L.energy );

    return numerical_flux;
  }
};

// ----------------------------------
//      initial conditions
// ----------------------------------

struct MyInitialConditions
{
  static void impose( Mesh &mesh )
  {
    for ( size_t i = 0; i < mesh.getSize(); i++ )
    {
    }

  }
};

template< typename T >
void InitialConditions( Mesh &mesh )
{
  T::impose( mesh );
}

// -----------------------------------
//      boundary conditions
// -----------------------------------

struct Zero_gradient
{
  static void impose( Mesh &mesh ) {};
};

template< typename T >
void BoundaryConditions( Mesh &mesh )
{
  T::impose( mesh );
}

// ----------------------------------
//            main
// ----------------------------------

int main()
{
  const size_t number_of_cells = 100;
  const size_t domain_length = 1;
  RealNumber t = 0.f;
  const RealNumber end_time = 0.25;
  
}
