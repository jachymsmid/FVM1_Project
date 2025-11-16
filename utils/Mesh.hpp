#pragma once

/* 
 * -----------------------------------
 *           Mesh class 
 * -----------------------------------
 */

#include <vector>

template
<
  class DataType,
  class RealNumber
>
class Mesh
{
public:
    Mesh() = default;

    // copy constructor
    Mesh( Mesh &mesh);

    size_t number_ghost_cells;
    size_t size;
    DataType  getData1() { return data1; }
    DataType  getData2() { return data2; }
    DataType  getData3() { return data3; }
    RealNumber getSpatialStep() { return spatial_step; }

    void write_data()
    {
      spatial_step = 1;
    }

private:
    DataType data1;
    DataType data2;
    DataType data3;
    RealNumber spatial_step;
};
