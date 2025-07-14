#pragma once

namespace EngineMath {

  /**
    * @brief Pi constant (π).
  */
  float constexpr pi = 3.14159265358979323846f; // Pi constant

  /**
    * @brief Euler's number (e).
  */
  float constexpr e = 2.71828182845904523536f; // Euler's number

  /**
		* @brief Small value for floating-point comparisons. The error margin used in various calculations.
  */
  float constexpr epsilon = 1e-9f; // Small value for floating-point comparisons

  /**
    * @brief Calculates the square root of a number using Newton's method.
    * @param number The number to calculate the square root of.
    * @return The square root of the number, or -1 for negative input.
  */
  float
  sqrt(float number) {
		if (number == 0) return 0; // The square root of 0 is 0
		if (number == 1) return 1; // The square root of 1 is 1
    if (number < 0) return -1; // Return -1 for negative numbers (not defined in real numbers)

    float result = number / 2;
    for (int i = 0; i < 10; i++) { // Newton's method for square root approximation
			result = 0.5f * (result + number / result);
    }

    return result;
  }

  /**
    * @brief Returns the square of a number.
    * @param number The number to square.
    * @return The squared value.
  */
  inline float
  square(float number) {
    return number * number;
  }

  /**
    * @brief Returns the cube of a number.
    * @param number The number to cube.
    * @return The cubed value.
  */
  inline float
  cube(float number) {
    return number * number * number;
  }

  /**
    * @brief Raises a base to an integer exponent.
    * @param base The base value.
    * @param exponent The exponent value.
    * @return The result of base raised to the exponent.
  */
  inline float
  power(float base, int exponent) {
    float result = 1.0f;
		bool expNeg = exponent < 0;
		exponent = absf(exponent);

    for (int i = 0; i < exponent; i++) {
      result *= base;
    }

    return expNeg ? 1.0f / result : result;
  }

	/**
    * @brief Raises a base to a real (floating-point) exponent.
    * @param base The base value.
    * @param exponent The exponent value.
    * @return The result of base raised to the exponent.
  */
	inline double
  powerf(float base, float exponent) {
    if (base <= 0.0f) return 0.0f; // Avoid negative base for real numbers
    if (exponent == 0.0f) return 1.0f; // Any number to the power of 0 is 1
    if (exponent == 1.0f) return base; // Any number to the power of 1 is itself
		
		int intExponent = absf(floor(exponent));
		double floatExponent = absf(exponent) - intExponent;

    if (floatExponent < 0.0f) { intExponent -= 1; floatExponent += 1.0f; }

		double intPart = power(base, intExponent);

		double decimalPart = 1.0f;
    float factor = base;
    double binaryFractions = floatExponent;

    for (int i = 0; i < 32 && binaryFractions > epsilon; i++) {
      factor = sqrt(factor);
			binaryFractions *= 2.0f;

      if (binaryFractions >= 1.0f) {
        decimalPart *= factor;
        binaryFractions -= 1.0f;
      }
    }

		double result = intPart * decimalPart;
		return exponent < 0.0f ? 1 / result : result;
	}

  /**
    * @brief Returns the absolute value of an integer.
    * @param number The integer value.
    * @return The absolute value.
  */
  inline int
  abs(int number) {
    return (number < 0) ? -1 * number : number;
  }

  /**
    * @brief Returns the absolute value of a floating-point number.
    * @param number The floating-point value.
    * @return The absolute value.
  */
  inline float
  absf(float number) {
    return (number < 0.0f) ? number * -1 : number;
  }

  /**
    * @brief Returns the maximum of two floating-point numbers.
    * @param a First value.
    * @param b Second value.
    * @return The maximum value.
  */
  inline float
  eMax(float a, float b) {
    return a > b ? a : b;
  }

  /**
    * @brief Returns the minimum of two floating-point numbers.
    * @param a First value.
    * @param b Second value.
    * @return The minimum value.
  */
  inline float
  eMin(float a, float b) {
    return a < b ? a : b;
  }

  /**
    * @brief Rounds a floating-point number to the nearest integer.
    * @param number The number to round.
    * @return The rounded integer value.
  */
  inline int
  round(float number) {
    int intPart = static_cast<int>(number);
    float fPart = number - intPart;

    return fPart >= 0.5f ? intPart + 1 : intPart;
  }

  /**
    * @brief Returns the largest integer less than or equal to the given number.
    * @param number The number to floor.
    * @return The floored integer value.
  */
  inline int
  floor(float number) {
    int newValue = static_cast<int>(number);
    
		if (number < 0.0f && number != newValue) newValue -= 1; // Adjust for negative numbers
		return newValue;
  }

  /**
    * @brief Returns the smallest integer greater than or equal to the given number.
    * @param number The number to ceil.
    * @return The ceiled integer value.
  */
  inline int
  ceil(float number) {
    return static_cast<int>(number) + 1;
  }

  /**
    * @brief Returns the floating-point remainder of a divided by b.
    * @param a Dividend.
    * @param b Divisor.
    * @return The remainder, or 0 if b is zero.
  */
  inline float
  mod(float a, float b) {
    if (b == 0) {
      return 0; // Avoid division by zero
    }
    return (a / b) - (static_cast<int>(a / b));
  }

  /**
    * @brief Calculates the exponential function e^x.
    * @param exponent The exponent value.
    * @return The result of e raised to the exponent.
  */
  inline float
  exp(float exponent) {
    return powerf(e, exponent); // Using the power function to calculate e^x
  }

  /**
    * @brief Calculates the natural logarithm (ln) of a number.
    * @param number The number to calculate the logarithm of.
    * @return The natural logarithm, or 0 for non-positive input.
  */
  inline double
  ln(float number) {
    if (number <= 0) return 0; // Logarithm is undefined for non-positive numbers
    
		// Normalize number to [0, 2.718] range
    int k = 0;
    while (number >= e) {
      number /= e;
      k++;
		}
    while (number < 1.0f) {
      number *= e;
      k--;
		}

    int n = 0;
    double result = number - 1.0f;
    
    while(absf(exp(result) - number) > epsilon && n < 32) {
			double expResult = exp(result);
			result += (number - expResult) / expResult;
      n++;
    }

		return result + k;
  }

  /**
    * @brief Calculates the base-10 logarithm of a number.
    * @param number The number to calculate the logarithm of.
    * @return The base-10 logarithm, or 0 for non-positive input.
  */
  inline float
  log10(float number) {
    if (number <= 0) return 0; // Logarithm is undefined for non-positive numbers
    
		return ln(number) / ln(10.0f); // Using natural logarithm to calculate log base 10
	}
  
  /**
    * @brief Calculates the sine of an angle (in radians) using Taylor series.
    * @param angle The angle in radians.
    * @return The sine of the angle.
  */
  inline double
  sin(float angle) {
    int n = 0;
    float result = 0.0f;
    float currentTerm = angle;

    while (absf(currentTerm) >= epsilon && n < 32) {
      result += currentTerm;
      n++;
      currentTerm *= (-1.0 * square(angle)) / ((2.0 * n) * (2.0 * n + 1.0));
    }

    return result;
  }

  /**
    * @brief Calculates the cosine of an angle (in radians).
    * @param angle The angle in radians.
    * @return The cosine of the angle.
  */
  inline double
  cos(float angle) {
    return sqrt(1.0f - square(sin(angle)));
  }

  /**
    * @brief Calculates the tangent of an angle (in radians).
    * @param angle The angle in radians.
    * @return The tangent of the angle.
  */
  inline double
  tan(float angle) {
    return sin(angle) / cos(angle);
  }

  /**
    * @brief Calculates the arcsine (inverse sine) of a value (in radians).
    * @param angle The value to calculate arcsine for.
    * @return The arcsine in radians, or 0 if out of range.
  */
  inline double
  aSin(float angle) {
    if (angle < -1.0f || angle > 1.0f) {
      return 0.0f; // Out of range for arcsin
    }

    float result = angle + (power(angle, 3) / 6);
    float i = 0;

    while (absf(sin(result) - angle) > epsilon && i < 32) {
      result -= (sin(result) - angle) / cos(result);
      i++;
    }

    return result;
  }

  /**
    * @brief Calculates the arccosine (inverse cosine) of a value (in radians).
    * @param angle The value to calculate arccosine for.
    * @return The arccosine in radians.
  */
  inline double
  aCos(float angle) {
    return pi / 2.0f - aSin(angle);
  }

  /**
    * @brief Calculates the arctangent (inverse tangent) of a value (in radians).
    * @param angle The value to calculate arctangent for.
    * @return The arctangent in radians.
  */
  inline double
  aTan(float angle) {
    bool invert = false;
    float x = angle;
    if (angle < -1.0f || angle > 1.0f) {
      invert = true;
      x = 1.0f / angle; // Invert if out of range for arctangent
    }

    float currentTerm = x;
    float result = currentTerm;
    int i = 1;

    while (absf(currentTerm) > epsilon && i < 32) {
      currentTerm = power(-1.0f, i) * power(x, 2 * i + 1) / (2 * i + 1);
      result += currentTerm;
      i++;
    }

    if (invert) {
      angle > 0.0f ? result = (pi / 2.0f) - result : result = -(pi / 2.0f) - result; // Adjust for quadrant
    }

    return result;
  }

  /**
    * @brief Calculates the hyperbolic sine of an angle (in radians).
    * @param angle The angle in radians.
    * @return The hyperbolic sine value.
  */
  inline double
  sinH(float angle) {
    return (exp(angle) - exp(-angle)) / 2.0f;
  }
  
  /**
    * @brief Calculates the hyperbolic cosine of an angle (in radians).
    * @param angle The angle in radians.
    * @return The hyperbolic cosine value.
  */
  inline double
  cosH(float angle) {
    return (exp(angle) + exp(-angle)) / 2.0f;
  }
  
  /**
    * @brief Calculates the hyperbolic tangent of an angle (in radians).
    * @param angle The angle in radians.
    * @return The hyperbolic tangent value.
  */
  inline double
  tanh(float angle) {
    return sinH(angle) / cosH(angle);
  }

  /**
    * @brief Converts degrees to radians.
    * @param degrees The angle in degrees.
    * @return The angle in radians.
  */
  inline float
  radians(float degrees) {
    return (degrees * pi) / 180.0f;
  }
  
  /**
    * @brief Converts radians to degrees.
    * @param radian The angle in radians.
    * @return The angle in degrees.
  */
  inline float
  degrees(float radian) {
    return (radian * 180) / pi;
  }

  /**
    * @brief Calculates the area of a circle.
    * @param radius The radius of the circle.
    * @return The area of the circle.
  */
  inline float
  circleArea(float radius) {
    return pi * radius * radius;
  }

  /**
    * @brief Calculates the circumference of a circle.
    * @param radius The radius of the circle.
    * @return The circumference of the circle.
  */
  inline float
  circleCircumference(float radius) {
    return 2.0f * pi * radius;
  }

  /**
    * @brief Calculates the area of a rectangle.
    * @param width The width of the rectangle.
    * @param height The height of the rectangle.
    * @return The area of the rectangle.
  */
  inline float
  rectArea(float width, float height) {
    return width * height;
  }

  /**
    * @brief Calculates the perimeter of a rectangle.
    * @param width The width of the rectangle.
    * @param height The height of the rectangle.
    * @return The perimeter of the rectangle.
  */
  inline float
  rectPerimeter(float width, float height) {
    return 2.0f * (width + height);
  }

  /**
    * @brief Calculates the area of a triangle.
    * @param base The base of the triangle.
    * @param height The height of the triangle.
    * @return The area of the triangle.
  */
  inline float
  triArea(float base, float height) {
    return 0.5f * base * height;
  }

  /**
    * @brief Calculates the perimeter of a triangle given three sides.
    * @param side1 First side.
    * @param side2 Second side.
    * @param side3 Third side.
    * @return The perimeter of the triangle.
  */
  inline float
  triPerimeter(float side1, float side2, float side3) {
    return side1 + side2 + side3;
  }

  /**
    * @brief Calculates the perimeter of an equilateral triangle.
    * @param side The length of one side.
    * @return The perimeter of the triangle.
  */
  inline float
  triPerimeter(float side) {
    return 3 * side;
  }

  /**
    * @brief Calculates the distance between two points in 2D space.
    * @param x1 X coordinate of the first point.
    * @param y1 Y coordinate of the first point.
    * @param x2 X coordinate of the second point.
    * @param y2 Y coordinate of the second point.
    * @return The distance between the two points.
  */
  inline float
  distance(float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;

    return sqrt(dx * dx + dy * dy);
  }

  /**
    * @brief Performs linear interpolation between two values.
    * @param start The start value.
    * @param end The end value.
    * @param t Interpolation factor (0.0 to 1.0).
    * @return The interpolated value.
  */
  inline float
  lerp(float start, float end, float t) {
    return start + (end - start) * t;
  }

  /**
    * @brief Calculates the factorial of an integer.
    * @param number The integer value.
    * @return The factorial of the number.
  */
  inline long
  factorial(int number) {
    int result = 1;
    for (int i = number; i > 0; i--) {
      result *= i;
    }
    return result;
  }

  /**
    * @brief Checks if two floating-point numbers are approximately equal within a given error.
    * @param a First value.
    * @param b Second value.
    * @param error Allowed error margin.
    * @return True if the values are approximately equal, false otherwise.
  */
	inline bool
  approxEqual(float a, float b, float error) {
    return absf(a - b) < error;
  }

  /**
    * @brief Clamps a value between a minimum and maximum.
    * @param value The value to clamp.
    * @param min The minimum value.
    * @param max The maximum value.
    * @return The clamped value.
  */
  inline float
  clamp(float value, float min, float max) {
    return eMax(min, eMin(value, max));
  }

  /**
    * @brief Returns the sign of a floating-point value.
    * @param value The value to check.
    * @return 1 if positive, -1 if negative, 0 if zero.
  */
  inline int
  sign(float value) {
    return (value > 0.0f) - (value < 0.0f);
	}
}