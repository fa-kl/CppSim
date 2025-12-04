/*****************************************************************************************
 * @file: Vector.hpp
 *
 * @brief: This file provides a mathematical vector class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <cmath>
#include <vector>

#include "types.hpp"

namespace sim
{

/**
 * @brief A type for booleans.
 */
using bool_t = bool;

/**
 * @brief A type for integers.
 */
using int_t = long;

/**
 * @brief A type for natural numbers.
 */
using uint_t = unsigned long;

/**
 * @brief A type for sizes.
 */
using size_t = uint_t;

/**
 * @brief A type for indices.
 */
using index_t = int_t;

/**
 * @brief A type for real numbers.
 */
using real_t = double;

/**
 * @brief A type for matrix dimensions.
 */
typedef struct {
  size_t rows;
  size_t cols;
} dimension_t;

/**
 * @brief Checks if two `dimension_t` elements are equal.
 */
bool_t operator==(dimension_t dim1, const dimension_t& dim2);

/**
 * @brief Checks if two `dimension_t` elements are not equal.
 */
bool_t operator!=(dimension_t dim1, const dimension_t& dim2);

}  // namespace sim