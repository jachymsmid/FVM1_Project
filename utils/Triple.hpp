#pragma once

/* 
 * -----------------------------------
 *        Triple struct
 * -----------------------------------
 */

template< class RealNumber>
struct Triple
{
  RealNumber data1;
  RealNumber data2;
  RealNumber data3;


  Triple operator+( Triple &other )
  {
    Triple temp;
    temp.data1 = data1 + other.data1;
    temp.data2 = data2 + other.data2;
    temp.data3 = data3 + other.data3;
    return temp;
  }

  Triple operator*( Triple &other )
  {
    Triple temp;
    temp.data1 = data1 * other.data1;
    temp.data2 = data2 * other.data2;
    temp.data3 = data3 * other.data3;
    return temp;
  }

  Triple operator-( Triple &other )
  {
    Triple temp;
    temp.data1 = data1 - other.data1;
    temp.data2 = data2 - other.data2;
    temp.data3 = data3 - other.data3;
    return temp;
  }

  Triple operator/( Triple &other )
  {
    Triple temp;
    temp.data1 = data1 / other.data1;
    temp.data2 = data2 / other.data2;
    temp.data3 = data3 / other.data3;
    return temp;
  }

  Triple operator+( RealNumber scalar )
  {
    Triple temp;
    temp.data1 = data1 + scalar;
    temp.data2 = data2 + scalar;
    temp.data3 = data3 + scalar;
    return temp;
  }

  Triple operator*( RealNumber scalar )
  {
    Triple temp;
    temp.data1 = data1 * scalar;
    temp.data2 = data2 * scalar;
    temp.data3 = data3 * scalar;
    return temp;
  }

  Triple operator-( RealNumber scalar )
  {
    Triple temp;
    temp.data1 = data1 - scalar;
    temp.data2 = data2 - scalar;
    temp.data3 = data3 - scalar;
    return temp;
  }

  Triple operator/( RealNumber scalar )
  {
    Triple temp;
    temp.data1 = data1 / scalar;
    temp.data2 = data2 / scalar;
    temp.data3 = data3 / scalar;
    return temp;
  }

  Triple primitive_to_conserved( const Triple primitive, RealNumber gamma_g )
  {
    Triple conserved;

    conserved.data1 = primitive.data1;
    conserved.data2 = primitive.data2 * primitive.data3;
    conserved.data3 = primitive.data3 * (gamma_g - 1.0) + 0.5 * primitive.data1 * primitive.data2 * primitive.data2 ;

    return conserved;
  }

  Triple conserved_to_primitive( const Triple conserved, RealNumber gamma_g )
  {
    Triple primitive;

    primitive.data1 = conserved.data1;
    primitive.data2 = conserved.data2 / conserved.data1;
    RealNumber internal = conserved.data3 - 0.5 * primitive.data1 * primitive.data2 * primitive.data2;
    primitive.data3 = ( gamma_g - 1 ) * internal;

    return primitive;
  }
};



