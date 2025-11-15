#include "utils.h"

RealNumber gamma_g = 1.4;

Triple Triple::operator+( Triple &other )
{
  Triple temp;
  temp.data1 = data1 + other.data1;
  temp.data2 = data2 + other.data2;
  temp.data3 = data3 + other.data3;
  return temp;
}

Triple Triple::operator*( Triple &other )
{
  Triple temp;
  temp.data1 = data1 * other.data1;
  temp.data2 = data2 * other.data2;
  temp.data3 = data3 * other.data3;
  return temp;
}

Triple Triple::operator-( Triple &other )
{
  Triple temp;
  temp.data1 = data1 - other.data1;
  temp.data2 = data2 - other.data2;
  temp.data3 = data3 - other.data3;
  return temp;
}

Triple Triple::operator/( Triple &other )
{
  Triple temp;
  temp.data1 = data1 / other.data1;
  temp.data2 = data2 / other.data2;
  temp.data3 = data3 / other.data3;
  return temp;
}

Triple Triple::operator+( RealNumber scalar )
{
  Triple temp;
  temp.data1 = data1 + scalar;
  temp.data2 = data2 + scalar;
  temp.data3 = data3 + scalar;
  return temp;
}

Triple Triple::operator*( RealNumber scalar )
{
  Triple temp;
  temp.data1 = data1 * scalar;
  temp.data2 = data2 * scalar;
  temp.data3 = data3 * scalar;
  return temp;
}

Triple Triple::operator-( RealNumber scalar )
{
  Triple temp;
  temp.data1 = data1 - scalar;
  temp.data2 = data2 - scalar;
  temp.data3 = data3 - scalar;
  return temp;
}

Triple Triple::operator/( RealNumber scalar )
{
  Triple temp;
  temp.data1 = data1 / scalar;
  temp.data2 = data2 / scalar;
  temp.data3 = data3 / scalar;
  return temp;
}

Triple Triple::primitive_to_conserved( const Triple primitive)
{
  Triple conserved;

  conserved.data1 = primitive.data1;
  conserved.data2 = primitive.data2 * primitive.data3;
  conserved.data3 = primitive.data3 * (gamma_g - 1.0) + 0.5 * primitive.data1 * primitive.data2 * primitive.data2 ;

  return conserved;
}

Triple Triple::conserved_to_primitive( const Triple conserved )
{
  Triple primitive;

  primitive.data1 = conserved.data1;
  primitive.data2 = conserved.data2 / conserved.data1;
  RealNumber internal = conserved.data3 - 0.5 * primitive.data1 * primitive.data2 * primitive.data2;
  primitive.data3 = ( gamma_g - 1 ) * internal;

  return primitive;
}

Triple flux( const Triple &conserved)
{
  Triple temp;
  const Triple primitive = Triple::conserved_to_primitive( conserved );
  temp.data1 = conserved.data2;
  temp.data2 = conserved.data2 * primitive.data2 + primitive.data3;
  temp.data3 = ( conserved.data3 + primitive.data3 ) * primitive.data2;
  return temp;
}

Triple Mesh::getValues( size_t i )
{
  Triple data;

  data.data1 = data1[ i ];
  data.data2 = data2[ i ];
  data.data3 = data3[ i ];

  return data;
}

std::vector< RealNumber > Mesh::getData1()
{
  return data1;
}

std::vector< RealNumber > Mesh::getData2()
{
  return data2;
}

std::vector< RealNumber > Mesh::getData3()
{
  return data3;
}

Triple NumericalSolver::Rusanov::numerical_flux( Triple data_left, Triple data_right )
{
  Triple num_flux =  ( flux( data_left ) + flux( data_right ) ) +  ( data_left - data_right);
  return num_flux;
}

Triple NumericalSolver::rhs( Mesh &mesh )
{
   
}

void ODEsolver::Euler::solve( Mesh &mesh )
{

}

void ODEsolver::Heune::solve( Mesh &mesh )
{
}


void Sods_problem::impose( Mesh &mesh )
{
  size_t ng = mesh.number_ghost_cells;
  for ( size_t i = 0; i < mesh.size - 2*ng; i++ )
  {
    mesh.getData1()[ i + ng] = 1.0;
    mesh.getData1()[ i + ng] = 1.0;
    mesh.getData1()[ i + ng] = 1.0;
  }
}

void Zero_gradient::impose( Mesh &mesh )
{
  size_t ng = mesh.number_ghost_cells;
  for ( size_t i = 0; i <  ng; i++ )
  {
    mesh.getData1()[ (ng - 1) - i ] = mesh.getData1()[ ng - i ];
    mesh.getData1()[ (mesh.size - ng) + i ] = mesh.getData1()[ (mesh.size - ng) + (i - 1) ];
    mesh.getData2()[ (ng - 1) - i ] = mesh.getData2()[ ng - i ];
    mesh.getData2()[ (mesh.size - ng) + i ] = mesh.getData2()[ (mesh.size - ng) + (i - 1) ];
    mesh.getData3()[ (ng - 1) - i ] = mesh.getData3()[ ng - i ];
    mesh.getData3()[ (mesh.size - ng) + i ] = mesh.getData3()[ (mesh.size - ng) + (i - 1) ];
  }
}
