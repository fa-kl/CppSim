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

#include <ostream>
#include <vector>

#include "types.hpp"

namespace sim
{

/**
 * @brief A class for n-dimensional real-valued vectors.
 *
 * @details This class provides a mathematical vector implementation with operations
 * for linear algebra, including arithmetic operations, norms, and element-wise functions.
 */
class Vector
{
  protected:
    /**
     * @brief The underlying data storage.
     */
    std::vector<real_t> m_data;

  public:
    /**
     * @brief Creates an empty vector.
     */
    Vector();

    /**
     * @brief Creates a vector with n elements initialized to zero.
     *
     * @param n Number of elements.
     */
    explicit Vector(size_t n);

    /**
     * @brief Creates a vector from an initializer list.
     *
     * @param values Initializer list of values.
     */
    Vector(std::initializer_list<real_t> values);

    /**
     * @brief Copy constructor.
     *
     * @param other Vector to copy.
     */
    Vector(const Vector& other);

    /**
     * @brief Move constructor.
     *
     * @param other Vector to move.
     */
    Vector(Vector&& other) noexcept;

    /**
     * @brief Destructor.
     */
    ~Vector();

    /**
     * @brief Copy assignment operator.
     *
     * @param other Vector to copy.
     * @return Reference to this vector.
     */
    Vector& operator=(const Vector& other);

    /**
     * @brief Move assignment operator.
     *
     * @param other Vector to move.
     * @return Reference to this vector.
     */
    Vector& operator=(Vector&& other) noexcept;

    /**
     * @brief Get the length (number of elements) of the vector.
     *
     * @return Length of the vector.
     */
    size_t length() const;

    /**
     * @brief Get the size dimensions of the vector.
     *
     * @return Size as dimension_t structure.
     */
    dimension_t size() const;

    /**
     * @brief Get the number of rows.
     *
     * @return Number of rows.
     */
    size_t rows() const;

    /**
     * @brief Get the number of columns.
     *
     * @return Number of columns (always 1 for vectors).
     */
    size_t cols() const;

    /**
     * @brief Check if the vector is empty.
     *
     * @return True if empty, false otherwise.
     */
    bool_t isEmpty() const;

    /**
     * @brief Access element at index i (0-based indexing).
     *
     * @param i Index.
     * @return Reference to element at index i.
     */
    real_t& operator[](size_t i);

    /**
     * @brief Access element at index i (0-based indexing).
     *
     * @param i Index.
     * @return Const reference to element at index i.
     */
    const real_t& operator[](size_t i) const;

    /**
     * @brief Access element at index i (1-based indexing, supports negative indices).
     *
     * @param i Index (1-based, negative indices count from end).
     * @return Reference to element at index i.
     */
    real_t& operator()(index_t i);

    /**
     * @brief Access element at index i (1-based indexing, supports negative indices).
     *
     * @param i Index (1-based, negative indices count from end).
     * @return Const reference to element at index i.
     */
    const real_t& operator()(index_t i) const;

    /**
     * @brief Unary plus operator.
     *
     * @return A copy of the vector.
     */
    Vector operator+() const;

    /**
     * @brief Unary minus operator (negate vector).
     *
     * @return A negated copy of the vector.
     */
    Vector operator-() const;

    /**
     * @brief Add two vectors element-wise.
     *
     * @param rhs Right-hand side vector.
     * @return Result of element-wise addition.
     */
    Vector operator+(const Vector& rhs) const;

    /**
     * @brief Add a scalar to all elements of the vector.
     *
     * @param scalar Scalar value.
     * @return Result vector.
     */
    Vector operator+(real_t scalar) const;

    /**
     * @brief Subtract two vectors element-wise.
     *
     * @param rhs Right-hand side vector.
     * @return Result of element-wise subtraction.
     */
    Vector operator-(const Vector& rhs) const;

    /**
     * @brief Subtract a scalar from all elements of the vector.
     *
     * @param scalar Scalar value.
     * @return Result vector.
     */
    Vector operator-(real_t scalar) const;

    /**
     * @brief Compute dot product with another vector.
     *
     * @param rhs Right-hand side vector.
     * @return Dot product (scalar value).
     */
    real_t operator*(const Vector& rhs) const;

    /**
     * @brief Multiply all elements by a scalar.
     *
     * @param scalar Scalar value.
     * @return Result vector.
     */
    Vector operator*(real_t scalar) const;

    /**
     * @brief Divide all elements by a scalar.
     *
     * @param scalar Scalar value (divisor).
     * @return Result vector.
     */
    Vector operator/(real_t scalar) const;

    /**
     * @brief Element-wise equality comparison.
     *
     * @param rhs Right-hand side vector.
     * @return True if all elements are equal, false otherwise.
     */
    bool_t operator==(const Vector& rhs) const;

    /**
     * @brief Element-wise inequality comparison.
     *
     * @param rhs Right-hand side vector.
     * @return True if any elements differ, false otherwise.
     */
    bool_t operator!=(const Vector& rhs) const;

    /**
     * @brief Multiply a scalar by a vector.
     *
     * @param scalar Scalar value.
     * @param vec Vector.
     * @return Result vector.
     */
    friend Vector operator*(real_t scalar, const Vector& vec);

    /**
     * @brief Add a scalar to a vector.
     *
     * @param scalar Scalar value.
     * @param vec Vector.
     * @return Result vector.
     */
    friend Vector operator+(real_t scalar, const Vector& vec);

    /**
     * @brief Subtract a vector from a scalar.
     *
     * @param scalar Scalar value.
     * @param vec Vector.
     * @return Result vector.
     */
    friend Vector operator-(real_t scalar, const Vector& vec);

    /**
     * @brief Divide a scalar by all elements of a vector.
     *
     * @param scalar Scalar value (dividend).
     * @param vec Vector (divisor).
     * @return Result vector.
     */
    friend Vector operator/(real_t scalar, const Vector& vec);

    /**
     * @brief Output stream operator for vectors.
     *
     * @param os Output stream.
     * @param vec Vector to output.
     * @return Output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Vector& vec);
};

/*****************************************************************************************
 * Vector Functions
 ****************************************************************************************/

/**
 * @brief Create a vector of zeros.
 *
 * @param n Length of the vector.
 * @return Vector filled with zeros.
 */
Vector zeros(size_t n);

/**
 * @brief Create a vector of ones.
 *
 * @param n Length of the vector.
 * @return Vector filled with ones.
 */
Vector ones(size_t n);

/**
 * @brief Get the length of a vector.
 *
 * @param vec Input vector.
 * @return Length of the vector.
 */
size_t length(const Vector& vec);

/**
 * @brief Get the size of a vector.
 *
 * @param vec Input vector.
 * @return Dimension of the vector.
 */
dimension_t size(const Vector& vec);

/**
 * @brief Check if a vector is empty.
 *
 * @param vec Input vector.
 * @return True if empty, false otherwise.
 */
bool_t isEmpty(const Vector& vec);

/**
 * @brief Compute the p-norm of a vector.
 *
 * @param vec Input vector.
 * @param p Norm parameter (default: 2 for Euclidean norm).
 * @return p-norm of the vector.
 */
real_t norm(const Vector& vec, real_t p = 2.0);

/**
 * @brief Normalize a vector (divide by its norm).
 *
 * @param vec Input vector.
 * @return Normalized vector.
 */
Vector normalize(const Vector& vec);

/**
 * @brief Compute absolute value of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with absolute values.
 */
Vector abs(const Vector& vec);

/**
 * @brief Compute exponential of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with exponential values.
 */
Vector exp(const Vector& vec);

/**
 * @brief Compute natural logarithm of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with logarithmic values.
 */
Vector log(const Vector& vec);

/**
 * @brief Compute sine of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with sine values.
 */
Vector sin(const Vector& vec);

/**
 * @brief Compute cosine of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with cosine values.
 */
Vector cos(const Vector& vec);

/**
 * @brief Compute tangent of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with tangent values.
 */
Vector tan(const Vector& vec);

/**
 * @brief Compute arcsine of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with arcsine values.
 */
Vector asin(const Vector& vec);

/**
 * @brief Compute arccosine of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with arccosine values.
 */
Vector acos(const Vector& vec);

/**
 * @brief Compute arctangent of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with arctangent values.
 */
Vector atan(const Vector& vec);

/**
 * @brief Raise each element of a vector to a power.
 *
 * @param vec Input vector.
 * @param exponent Power exponent.
 * @return Vector with powered values.
 */
Vector pow(const Vector& vec, real_t exponent);

/**
 * @brief Compute square root of each element in a vector.
 *
 * @param vec Input vector.
 * @return Vector with square root values.
 */
Vector sqrt(const Vector& vec);

/**
 * @brief Compute sum of all elements in a vector.
 *
 * @param vec Input vector.
 * @return Sum of elements.
 */
real_t sum(const Vector& vec);

/**
 * @brief Compute differences between consecutive elements.
 *
 * @param vec Input vector.
 * @return Vector of differences (length n-1).
 */
Vector diff(const Vector& vec);

/**
 * @brief Element-wise multiplication of two vectors.
 *
 * @param lhs Left-hand side vector.
 * @param rhs Right-hand side vector.
 * @return Vector with element-wise products.
 */
Vector elmul(const Vector& lhs, const Vector& rhs);

/**
 * @brief Element-wise division of two vectors.
 *
 * @param lhs Left-hand side vector (dividend).
 * @param rhs Right-hand side vector (divisor).
 * @return Vector with element-wise quotients.
 */
Vector eldiv(const Vector& lhs, const Vector& rhs);

/**
 * @brief Compute dot product of two vectors.
 *
 * @param lhs Left-hand side vector.
 * @param rhs Right-hand side vector.
 * @return Dot product (scalar value).
 */
real_t dot(const Vector& lhs, const Vector& rhs);

/**
 * @brief Compute cross product of two 3D vectors.
 *
 * @param lhs Left-hand side vector (must be 3D).
 * @param rhs Right-hand side vector (must be 3D).
 * @return Cross product vector.
 */
Vector cross(const Vector& lhs, const Vector& rhs);

/**
 * @brief Find maximum element in a vector.
 *
 * @param vec Input vector.
 * @return Maximum value.
 */
real_t max(const Vector& vec);

/**
 * @brief Find minimum element in a vector.
 *
 * @param vec Input vector.
 * @return Minimum value.
 */
real_t min(const Vector& vec);

/**
 * @brief Compute mean (average) of vector elements.
 *
 * @param vec Input vector.
 * @return Mean value.
 */
real_t mean(const Vector& vec);

/**
 * @brief Compute variance of vector elements.
 *
 * @param vec Input vector.
 * @return Variance.
 */
real_t var(const Vector& vec);

/**
 * @brief Compute standard deviation of vector elements.
 *
 * @param vec Input vector.
 * @return Standard deviation.
 */
real_t std(const Vector& vec);

/**
 * @brief Vertically concatenate two vectors.
 *
 * @param lhs First vector.
 * @param rhs Second vector.
 * @return Concatenated vector.
 */
Vector vcat(const Vector& lhs, const Vector& rhs);

}  // namespace sim