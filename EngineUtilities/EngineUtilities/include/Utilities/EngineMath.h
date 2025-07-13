#pragma once

namespace EngineMath {

  float constexpr pi = 3.14159265358979323846f; // Pi constant
  float constexpr e = 2.71828182845904523536f; // Euler's number

  // Square Root
  float
  sqrt(float number) {
    if (number == 0) return 0;
    if (number == 1) return 1;
    if (number < 0) return -1; // Return -1 for negative numbers (not defined in real numbers)

    float xi = number / 2;
    for (int i = 0; i < 10; i++) { // Newton's method for square root approximation
			xi = 0.5f * (xi + number / xi); // More efficient version
    }

    return xi;
  }

  // Square
  inline float
  square(float number) {
    return number * number;
  }

  // Cube
  inline float
  cube(float number) {
    return number * number * number;
  }

  // Power
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

	// Power for Real Numbers
	inline double
  powerf(float base, float exponent) {
    if (base <= 0.0f) return 0.0f; // Avoid negative base for real numbers
    if (exponent == 0.0f) return 1.0f; // Any number to the power of 0 is 1
    if (exponent == 1.0f) return base; // Any number to the power of 1 is itself
		
		int n = absf(floor(exponent));
		double f = absf(exponent) - n;

    if (f < 0.0f) { n -= 1; f += 1.0f; }

		double intPart = power(base, n);

		double decimalPart = 1.0f;
    float factor = base;
    double binaryFractions = f;

    for (int i = 0; i < 32 && binaryFractions > 1e-9; i++) {
      factor = sqrt(factor);

			binaryFractions *= 2.0f;
      if (binaryFractions >= 1.0f) {
        decimalPart *= factor;
        binaryFractions -= 1.0f;
      }
    }

		double result = intPart * decimalPart;

		return exponent < 0.0f ? 1 / result : result; // Ensure positive result for real numbers
	}

  // Absolute
  inline int
  abs(int number) {
    return (number < 0) ? -1 * number : number;
  }

  // Fabs
  inline float
  absf(float number) {
    return (number < 0.0f) ? number * -1 : number;
  }

  // E Max
  inline float
  eMax(float a, float b) {
    return a > b ? a : b;
  }

  // E Min
  inline float
  eMin(float a, float b) {
    return a < b ? a : b;
  }

  // Round
  inline int
  round(float number) {
    int intPart = static_cast<int>(number);
    float fPart = number - intPart;

    return fPart >= 0.5f ? intPart + 1 : intPart;
  }

  // Floor
  inline int
  floor(float number) {
    int newValue = static_cast<int>(number);
    
		if (number < 0.0f && number != newValue) newValue -= 1; // Adjust for negative numbers
		return newValue;
  }

  // Ceil
  inline int
  ceil(float number) {
    return static_cast<int>(number) + 1;
  }

  // Mod
  inline float
  mod(float a, float b) {
    if (b == 0) {
      return 0; // Avoid division by zero
    }
    return (a / b) - (static_cast<int>(a / b));
  }

  // Exp
  inline float
  exp(float exponent) {
    return powerf(e, exponent); // Using the power function to calculate e^x
  }

  // ln
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
    
    while(absf(exp(result) - number) > 1e-9 && n < 32) {
			double expResult = exp(result);
			result += (number - expResult) / expResult;
      n++;
    }

		return result + k;
  }

  // Log 10
  inline float
  log10(float number) {
    if (number <= 0) return 0; // Logarithm is undefined for non-positive numbers
    
		return ln(number) / ln(10.0f); // Using natural logarithm to calculate log base 10
	}
  
  // Sin
  inline double
  sin(float angle) {
    //angle = mod(angle, 2.0f * pi); // Normalize angle to [0, 2π]
    //if (angle > pi) angle -= 2.0f * pi; // Adjust to [-π, π]
    //if (angle < -pi) angle += 2.0f * pi; // Adjust to [-π, π]
    int n = 0;
    float result = 0.0f;
    float currentTerm = angle;

    while (absf(currentTerm) >= 0.0000001 && n < 25) {
      result += currentTerm;
      n++;
      //currentTerm = (power(-1, n) * power(angle, (2 * n) + 1)) / factorial((2 * n) + 1);
      currentTerm *= (-1.0 * square(angle)) / ((2.0 * n) * (2.0 * n + 1.0));
    }

    return result;
  }

  // Cos
  inline double
  cos(float angle) {
    return sqrt(1.0f - square(sin(angle)));
  }

  // Tan
  inline double
  tan(float angle) {
    return sin(angle) / cos(angle);
  }

  // aSin
  inline double
  aSin(float angle) {
    if (angle < -1.0f || angle > 1.0f) {
      return 0.0f; // Out of range for arcsin
    }

    float result = angle + (power(angle, 3) / 6);
    float i = 0;

    while (absf(sin(result) - angle) > 0.0001 && i < 30) {
      result -= (sin(result) - angle) / cos(result);
      i++;
    }

    return result;
  }

  // aCos
  inline double
  aCos(float angle) {
    return pi / 2.0f - aSin(angle);
  }

  // aTan
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

    while (absf(currentTerm) > 0.0001 && i < 30) {
      currentTerm = power(-1.0f, i) * power(x, 2 * i + 1) / (2 * i + 1);
      result += currentTerm;
      i++;
    }

    if (invert) {
      angle > 0.0f ? result = (pi / 2.0f) - result : result = -(pi / 2.0f) - result; // Adjust for quadrant
    }

    return result;
  }

  // Sinh
  inline double
  sinH(float angle) {
    return (exp(angle) - exp(-angle)) / 2.0f;
  }
  
  // Cosh
  inline double
  cosH(float angle) {
    return (exp(angle) + exp(-angle)) / 2.0f;
  }
  
  // Tanh
  inline double
  tanh(float angle) {
    return sinH(angle) / cosH(angle);
  }

  // Radians
  inline float
  radians(float degrees) {
    return (degrees * pi) / 180.0f;
  }
  
  // Degrees
  inline float
  degrees(float radian) {
    return (radian * 180) / pi;
  }

  // Circle Area
  inline float
  circleArea(float radius) {
    return pi * radius * radius;
  }

  // Circle Circumference
  inline float
  circleCircumference(float radius) {
    return 2.0f * pi * radius;
  }

  // Rectangle Area
  inline float
  rectArea(float width, float height) {
    return width * height;
  }

  // Rectangle Perimeter
  inline float
  rectPerimeter(float width, float height) {
    return 2.0f * (width + height);
  }

  // Triangle Area
  inline float
  triArea(float base, float height) {
    return 0.5f * base * height;
  }

  // Triangle Perimeter
  inline float
  triPerimeter(float side1, float side2, float side3) {
    return side1 + side2 + side3;
  }

  inline float
  triPerimeter(float side) {
    return 3 * side;
  }

  // Distance
  inline float
  distance(float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;

    return sqrt(dx * dx + dy * dy);
  }

  // Lerp
  inline float
  lerp(float start, float end, float t) {
    return start + (end - start) * t;
  }

  // Factorial
  inline long
  factorial(int number) {
    int result = 1;
    for (int i = number; i > 0; i--) {
      result *= i;
    }
    return result;
  }

  // Aprox Equal
	inline bool
  approxEqual(float a, float b, float error) {
    return absf(a - b) < error;
  }

  // Clamp
  inline float
  clamp(float value, float min, float max) {
    return eMax(min, eMin(value, max));
  }

  // Sign
  inline int
  sign(float value) {
    return (value > 0.0f) - (value < 0.0f);
	}
}