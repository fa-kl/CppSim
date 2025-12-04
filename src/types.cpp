/*****************************************************************************************
 * @file: types.cpp
 *
 * @brief: This file provides operators and functions for some basic types.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "types.hpp"

namespace sim
{

bool_t operator==(dimension_t dim1, const dimension_t& dim2)
{
  return dim1.rows == dim2.rows && dim1.cols == dim2.cols;
}

bool_t operator!=(dimension_t dim1, const dimension_t& dim2)
{
  return !(dim1 == dim2);
}

}  // namespace sim