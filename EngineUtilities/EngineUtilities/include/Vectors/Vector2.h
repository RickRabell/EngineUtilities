namespace EngineUtilities {
	class Vector2 {
	public:
		/**
			* @brief The x component of the vector.
		*/
		float m_x;
		/**
			* @brief The y component of the vector.
		*/
		float m_y;

		/**
			* @brief Default constructor. Initializes both components to zero.
		*/
		Vector2() { 
			m_x = 0; 
			m_y = 0; 
		};

		/**
			* @brief Initializes both components to the same value.
			* @param value Value to assign to both x and y.
		*/
		Vector2(float value) { 
			m_x = value; 
			m_y = value; 
		};

		/**
			* @brief Initializes the vector with specific x and y values.
			* @param x Value for the x component.
			* @param y Value for the y component.
		*/
		Vector2(float x, float y) { 
			m_x = x;
			m_y = y; 
		};

		/**
			* @brief Adds two vectors.
			* @param other The vector to add.
			* @return The result of vector addition.
		*/
		Vector2
		operator+(const Vector2& other) const {
			return Vector2(m_x + other.m_x, 
										 m_y + other.m_y);
		}

		/**
			* @brief Subtracts one vector from another.
			* @param other The vector to subtract.
			* @return The result of vector subtraction.
		*/
		Vector2
		operator-(const Vector2& other) const {
			return Vector2(m_x - other.m_x, 
										 m_y - other.m_y);
		}

		/**
			* @brief Multiplies the vector by a scalar.
			* @param scalar The scalar value.
			* @return The result of scalar multiplication.
		*/
		Vector2
		operator*(float scalar) const {
			return Vector2(m_x * scalar, 
										 m_y * scalar);
		}

		/**
			* @brief Divides the vector by a scalar.
			* @param scalar The scalar value.
			* @return The result of scalar division.
		*/
		Vector2
		operator/(float scalar) const {
			return Vector2(m_x / scalar, 
										 m_y / scalar);
		}

		/**
			* @brief Adds another vector to this vector in place.
			* @param other The vector to add.
			* @return Reference to this vector.
		*/
		Vector2&
		operator+=(const Vector2& other) {
			m_x += other.m_x;
			m_y += other.m_y;
			
			return *this;
		}

		/**
			* @brief Subtracts another vector from this vector in place.
			* @param other The vector to subtract.
			* @return Reference to this vector.
		*/
		Vector2&
		operator-=(const Vector2& other) {
			m_x -= other.m_x;
			m_y -= other.m_y;
			
			return *this;
		}

		/**
			* @brief Multiplies this vector by another vector in place.
			* @param other The vector to multiply.
			* @return Reference to this vector.
		*/
		Vector2&
		operator*=(const Vector2& other) {
			m_x *= other.m_x;
			m_y *= other.m_y;
			
			return *this;
		}

		/**
			* @brief Multiplies this vector by a scalar.
			* @param scalar Scalar value.
			* @return Reference to this vector.
		*/
		Vector2& operator*=(float scalar) {
			m_x *= scalar;
			m_y *= scalar;
			
			return *this;
		}

		/**
			* @brief Checks if two vectors are equal.
			* @param other The vector to compare.
			* @return True if both components are equal, false otherwise.
		*/
		bool
		operator==(const Vector2& other) const {
			return (m_x == other.m_x && 
							m_y == other.m_y);
		}

		/**
			* @brief Checks if two vectors are different.
			* @param other The vector to compare.
			* @return True if any component is different, false otherwise.
		*/
		bool
		operator!=(const Vector2& other) const {
			return (m_x != other.m_x || 
							m_y != other.m_y);
		}

		/**
			* @brief Accesses a component by index.
			* @param index 0 for x, 1 for y.
			* @return Reference to the component.
			* @throws "Index out of range for CVector2" if index is invalid.
		*/
		float&
		operator[](int index) {
			switch (index) {
				case 0: return m_x;
				case 1: return m_y;
				default: throw "Index out of range for Vector2";
			}
		}

		/**
			* @brief Calculates the length (magnitude) of a vector.
			* @param other The vector to calculate the length of.
			* @return The length of the vector.
		*/
		static float
		length(const Vector2& other) {
			return EngineMath::sqrt((other.m_x * other.m_x) + 
															(other.m_y * other.m_y));
		}

		/**
			* @brief Calculates the squared length of this vector.
			* @return The squared length of the vector.
		*/
		float 
		squaredLength() const {
			return (m_x * m_x) + (m_y * m_y);
		}

		/**
			* @brief Calculates the squared length of a vector.
			* @param other The vector to calculate the squared length of.
			* @return The squared length of the vector.
		*/
		static float
		squaredLength(const Vector2& other) {
			return (other.m_x * other.m_x) + 
						 (other.m_y * other.m_y);
		}

		/**
			* @brief Calculates the dot product of two vectors.
			* @param a First vector.
			* @param b Second vector.
			* @return The dot product.
		*/
		static float
		dotProduct(const Vector2& a, const Vector2& b) {
			return (a.m_x * b.m_x) + 
						 (a.m_y * b.m_y);
		}

		/**
			* @brief Calculates the cross product of two vectors.
			* @param a First vector.
			* @param b Second vector.
			* @return The cross product.
		*/
		float
		crossProduct(const Vector2& a, const Vector2& b) const {
			return (a.m_x * b.m_y) - 
						 (a.m_y * b.m_x);
		}

		/**
			* @brief Returns the normalized vector.
			* @param other The vector to normalize.
			* @return The normalized vector.
		*/
		Vector2
		normalize(const Vector2& other) const {
			return Vector2(other.m_x / other.length(other),
										 other.m_y / other.length(other));
		}

		/**
			* @brief Normalizes the vector in place.
			* @param other The vector to normalize.
			* @return Reference to the normalized vector.
		*/
		Vector2&
		normalize(Vector2& other) {
			float len = other.length(other);
			other.m_x /= len;
			other.m_y /= len;
			
			return other;
		}

		/**
			* @brief Calculates the distance between two vectors.
			* @param a First vector.
			* @param b Second vector.
			* @return The distance between the vectors.
		*/
		static float
		distance(const Vector2& a, const Vector2& b) {
			return (b - a).length(b - a);
		}

		/**
			* @brief Performs linear interpolation between two vectors.
			* @param a Start vector.
			* @param b End vector.
			* @param t Interpolation factor.
			* @return The interpolated vector.
		*/
		static Vector2
		lerp(const Vector2& a, const Vector2& b, float t) {
			return a + (b - a) * t;
		}

		/**
			* @brief Returns a zero vector.
			* @return Vector with both components set to zero.
		*/
		static Vector2 
		Zero() {
			return Vector2(0.0f, 0.0f);
		}

		/**
			* @brief Returns a unit vector.
			* @return Vector with both components set to one.
		*/
		static Vector2 
		Unit() {
			return Vector2(1.0f, 1.0f);
		}

		/**
			* @brief Performs spherical linear interpolation (slerp) between two vectors.
			* @param a Start vector.
			* @param b End vector.
			* @param t Interpolation factor.
			* @return The interpolated vector.
		*/
		static Vector2
		slerp(const Vector2& a, const Vector2& b, float t) {
			float dot = a.m_x * b.m_x + a.m_y * b.m_y;
			dot = EngineMath::clamp(dot, -1.0f, 1.0f);
			float theta = EngineMath::aCos(dot) * t;

			Vector2 relative = b - a * dot;
			float relLength = length(relative);

			if (relLength > 0.00001f) {
				relative = relative / relLength;
			} else {
				relative = a;
			}

			return a * EngineMath::cos(theta) + relative * EngineMath::sin(theta);
		}

		/**
			* @brief Reflects a vector off a surface with the given normal.
			* @param direction The incident vector.
			* @param normal The normal vector of the surface.
			* @return The reflected vector.
		*/
		static Vector2
		reflect(const Vector2& direction, const Vector2& normal) {
			float dot = dotProduct(direction, normal);
			return direction - normal * (2.0f * dot);
		}

		/**
			* @brief Projects one vector onto another.
			* @param a The vector to project.
			* @param b The vector to project onto.
			* @return The projected vector.
		*/
		static Vector2
		project(const Vector2& a, const Vector2& b) {
			float dot = dotProduct(a, b);
			float lenSq = squaredLength(b);

			return b * (dot / lenSq);
		}

		/**
			* @brief Clamps a vector between minimum and maximum vectors.
			* @param value The vector to clamp.
			* @param min Minimum vector.
			* @param max Maximum vector.
			* @return The clamped vector.
		*/
		static Vector2
		clamp(const Vector2& value, const Vector2& min, const Vector2& max) {
			float x = EngineMath::clamp(value.m_x, min.m_x, max.m_x);
			float y = EngineMath::clamp(value.m_y, min.m_y, max.m_y);

			return Vector2(x, y);
		}

		/**
			* @brief Calculates the angle between two vectors in radians.
			* @param a First vector.
			* @param b Second vector.
			* @return The angle in radians.
		*/
		static float
		angle(const Vector2& a, const Vector2& b) {
			float dot = dotProduct(a, b);
      float lenA = length(a);
      float lenB = length(b);

			return EngineMath::aCos(dot / (lenA * lenB));
		}
	};
}