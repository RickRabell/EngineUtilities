#pragma once

namespace EngineMath {

  float pi = 3.14159265358979323846f; // Pi constant
  float e = 2.71828182845904523536f; // Euler's number

  // Square Root
  inline float
  sqrt(float number) {
    if (number <= 0) {
      return 0; // Return 0 for negative numbers, as sqrt is not defined
    }

    float xi = number / 2;
    // float epsilon = 0.00001f; // Precision for convergence
    // while (abs(square(xi) - number) > epsilon) {
    
    for (int i = 0; i < 10; i++) { // Newton's method for square root approximation
      xi = xi - (square(xi) - number) / (2 * xi);
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
  power(float base, float exponent) {
    float result = 1.0f;
    for (int i = 0; i < exponent; i++) {
      result *= base;
    }
    return result;
  }

  // Absolute
  inline int
  abs(int number) {
    return (number < 0) ? -1 * number : number;
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
    return static_cast<int>(number);
  }

  // Ceil
  inline int
  ceil(float number) {
    return static_cast<int>(number) + 1;
  }

  // Fabs
  inline float
  fabs(float number) {
    return (number < 0.0f) ? number * -1 : number;
  }

  // Mod
  inline float
  mod(float a, float b) {
    if (b == 0) {
      return 0; // Avoid division by zero
    }
    return (a / b) - (static_cast<int>(a/b));
  }

  // Exp
  inline float
  exp(float exponent) {
    /*float result = 1.0f;
    for (int i = 0; i < 10; i++) { // Using Taylor series expansion for e^x
      result += power(exponent, i) / factorial(i);
    }
    return result;*/

    return power(e, exponent); // Using the power function to calculate e^x
  }

  // Log
  // Log 10
  // Sin
  inline float
  sin(float angle) {
    /*float result = 0.0f;
    float term = angle; // First term is angle itself
    int n = 1; // Counter for terms
    while (fabs(term) > 0.00001f) { // Continue until the term is small enough
      result += term;
      term *= -angle * angle / ((2 * n) * (2 * n + 1)); // Calculate next term in series
      n++;
    }
    return result;*/
    int n = 0;
    float result = 0;
    float currentTerm = angle;

    while (fabs(currentTerm) >= 0.0000001 || n == 0) {
      result += currentTerm;
      n++;
      //currentTerm = (power(-1, n) * power(angle, (2 * n) + 1)) / factorial((2 * n) + 1);
      currentTerm *= (-1.0 * angle * angle) / ((2.0 * n) * (2.0 * n + 1.0));
    }

    return result;
  }
  // Cos
  // Tan
  // aSin
  // aCos
  // aTan
  // Sinh
  // Cosh
  // Tanh
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
}