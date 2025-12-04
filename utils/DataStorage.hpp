#pragma once

#include "GeneralVector.hpp"
#include <fstream>
template < class RealNumber, std::size_t Size, std::size_t Length >
class DataStorage
{
public:

  //constructor
  DataStorage() = default;

  DataStorage( DataStorage &other )
  {
    for ( std::size_t i = 0; i < Size; i++ )
    {
      for ( std::size_t j = 0; j < Length; j++ )
      {
        data[i][j] = other[i][j];
      }
    }
  }

  // the data being stored
  std::array< Vector< RealNumber, Length >, Size > data;

  // getSize and getLength methods
  std::size_t getSize() { return Size; }
  std::size_t getLength() { return Length; }

  // acces the element at i,j
  RealNumber operator()(std::size_t i, std::size_t j)
  {
    return data[i][j];
  }

  // acces the cell values at index i
  Vector< RealNumber, Size > operator()(std::size_t i)
  {
    Vector< RealNumber, Size > temp;
    for ( std::size_t j = 0; j < Size; j++)
    {
      temp[j] = data[i][j];
    }
    return temp;
  }

  // print the data into a file
  void print( std::string file_name )
  {
    std::ofstream file;
    file.open( file_name );
    for ( std::size_t i = 0; i < Length; i++ )
    {
      for ( std::size_t j = 0; j < Size; j++ )
      {
        file << data[ j ][ i ] << " ";
      }
      file << std::endl;
    }
    file.close();
  }

};
