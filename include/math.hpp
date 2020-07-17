#ifndef MATH_HPP
#define MATH_HPP

#include <array>
#include <cmath>
#include <iostream>

#include "half.hpp"
using half = half_float::half;

#define PI 3.141592654f
#define INV_PI 0.318309886f
#define PI_OVER_180 0.01745329252f

#define INCH_TO_CM 2.54f
#define MM_TO_CM 0.1f

inline float max(const float a, const float b) {return a > b ? a : b;}
inline float min(const float a, const float b) {return a < b ? a : b;}

inline half max(const half a, const half b) {return a > b ? a : b;}
inline half min(const half a, const half b) {return a < b ? a : b;}

//////// VECTOR ////////

using Vec2i = std::array<int, 2>;
using Vec2f = std::array<float, 2>;
using Vec3i = std::array<int, 3>;
using Vec3f = std::array<float, 3>;
using Vec3h = std::array<half, 3>;

inline const Vec2f operator*	(const float f, const Vec2f v) { return {f*v[0], f*v[1]}; }
inline const Vec2f operator*	(const Vec2f v, const float f) { return {f*v[0], f*v[1]}; }
inline const Vec2f operator+	(const Vec2f a, const Vec2f b) { return {a[0]+b[0], a[1]+b[1]}; }
inline const Vec2f operator+	(const Vec2f v, const float f) { return {v[0]+f, v[1]+f}; }
inline const Vec2f operator-	(const Vec2f a, const Vec2f b) { return {a[0]-b[0], a[1]-b[1]}; }
inline const Vec2f operator-	(const Vec2f v, const float f) { return {v[0]-f, v[1]-f}; }
inline const Vec2f operator*	(const Vec2f a, const Vec2f b) { return {a[0]*b[0], a[1]*b[1]}; }
inline const Vec2f operator/	(const Vec2f a, const Vec2f b) { return {a[0]/b[0], a[1]/b[1]}; }
inline const Vec2f operator/	(const Vec2f v, const float f) { return {v[0]/f, v[1]/f}; }

inline std::ostream& operator<<(std::ostream& os, const Vec2f& v) {
	return os << "[" << v[0] << ", " << v[1] << "]";
}

inline std::ostream& operator<<(std::ostream& os, Vec2f& v) {
	return os << "[" << v[0] << ", " << v[1] << "]";
}

inline       float dot			(const Vec3f& a, const Vec3f& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline const Vec3f cross		(const Vec3f& a, const Vec3f& b) { return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}; }
inline const Vec3f operator-	(const Vec3f& a, const Vec3f& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
inline const Vec3f operator+	(const Vec3f& a, const Vec3f& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
inline const Vec3f operator*	(const Vec3f& a, const Vec3f& b) { return {a[0]*b[0], a[1]*b[1], a[2]*b[2]}; }
inline const Vec3f operator*	(const float f,  const Vec3f& v) { return {f*v[0], f*v[1], f*v[2]}; }
inline const Vec3f operator/	(const Vec3f& v, const float f) { return {v[0]/f, v[1]/f, v[2]/f}; }
inline       float length		(const Vec3f& v) { return std::sqrt(dot(v, v)); }
inline const Vec3f normalize	(const Vec3f& v) { return (1 / length(v))*v; }

inline std::ostream& operator<<(std::ostream& os, const Vec3f& v) {
	return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}

inline std::ostream& operator<<(std::ostream& os, Vec3f& v) {
	return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}

inline const half  dot			(const Vec3h& a, const Vec3h& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline const Vec3h cross		(const Vec3h& a, const Vec3h& b) { return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}; }
inline const Vec3h operator-	(const Vec3h& a, const Vec3h& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
inline const Vec3h operator+	(const Vec3h& a, const Vec3h& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
inline const Vec3h operator*	(const Vec3h& a, const Vec3h& b) { return {a[0]*b[0], a[1]*b[1], a[2]*b[2]}; }
inline const Vec3h operator*	(const half f,   const Vec3h& v) { return {f*v[0], f*v[1], f*v[2]}; }
inline const Vec3h operator/	(const Vec3h& v, const half   f) { return {v[0]/f, v[1]/f, v[2]/f}; }
inline const half  length		(const Vec3h& v) { return half_float::sqrt(dot(v, v)); }
inline const Vec3h normalize	(const Vec3h& v) { return (half)(1.0 / length(v))*v; }

inline const Vec3h fromVec3f(const Vec3f& v) { 
	return
	{
		(half_float::half)v[0], 
		(half_float::half)v[1], 
		(half_float::half)v[2]
	}; 
}

inline std::ostream& operator<<(std::ostream& os, const Vec3h& v) {
	return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}

inline std::ostream& operator<<(std::ostream& os, Vec3h& v) {
	return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}

//////// MATRIX ////////

// Stored column-wise -- as OpenGL does.
class Mat4f {
	std::array<float, 16> data;
public:
	Mat4f() : data {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} {}

	Mat4f(	float a00, float a10, float a20, float a30,
			float a01, float a11, float a21, float a31,
			float a02, float a12, float a22, float a32,
			float a03, float a13, float a23, float a33)
			: data{	a00, a10, a20, a30,
					a01, a11, a21, a31,
					a02, a12, a22, a32,
					a03, a13, a23, a33} {}

	// This assumes vectors are rows
	Mat4f(	const Vec3f& x, const Vec3f& y, const Vec3f& z, const Vec3f& w)
			: data{	x[0], y[0], z[0], w[0],
					x[1], y[1], z[1], w[1],
					x[2], y[2], z[2], w[2],
					   0,    0,    0,    1} {}

	const 	float& operator[](int idx) 		const	{ return data[idx]; }
			float& operator()(int i, int j)			{ return data[i+4*j]; }
	const 	float& operator()(int i, int j) const	{ return data[i+4*j]; }
};

inline const Vec3f transformPoint(const Mat4f& m, const Vec3f& p) {
	return {
		m(0,0)*p[0] + m(1,0)*p[1] + m(2,0)*p[2] + m(3,0), 
		m(0,1)*p[0] + m(1,1)*p[1] + m(2,1)*p[2] + m(3,1), 
		m(0,2)*p[0] + m(1,2)*p[1] + m(2,2)*p[2] + m(3,2)
	};
}

inline const Vec3f transformVector(const Mat4f& m, const Vec3f& v) { 
	return {
		m(0,0)*v[0] + m(1,0)*v[1] + m(2,0)*v[2], 
		m(0,1)*v[0] + m(1,1)*v[1] + m(2,1)*v[2], 
		m(0,2)*v[0] + m(1,2)*v[1] + m(2,2)*v[2]
	};
}

inline const Mat4f refFromVec (const Vec3f& v) {
	Vec3f v2{};
	if (std::abs(v[0]) > std::abs(v[1]))
		v2 = Vec3f{-v[2], 0, v[0]} / std::sqrt(v[0] * v[0] + v[2] * v[2]);
	else
		v2 = Vec3f{0, v[2], -v[1]} / std::sqrt(v[1] * v[1] + v[2] * v[2]);
	
	const auto v3 = cross(v, v2);
	return {v, v2, v3, {}};
}

inline const Mat4f transpose (const Mat4f& m) {
	return {	m[0], m[4],  m[8], m[12],
				m[1], m[5],  m[9], m[13],
				m[2], m[6], m[10], m[14],
				m[3], m[7], m[11], m[15]};
}

inline const Mat4f operator*(const Mat4f& a, const Mat4f& b)
{
	return {
		a[0]*b[0]  + a[4]*b[1]  + a[8]*b[2]  + a[12]*b[3],
		a[1]*b[0]  + a[5]*b[1]  + a[9]*b[2]  + a[13]*b[3],
		a[2]*b[0]  + a[6]*b[1]  + a[10]*b[2]  + a[14]*b[3],
		a[3]*b[0]  + a[7]*b[1]  + a[11]*b[2]  + a[15]*b[3],

		a[0]*b[4]  + a[4]*b[5]  + a[8]*b[6]  + a[12]*b[7],
		a[1]*b[4]  + a[5]*b[5]  + a[9]*b[6]  + a[13]*b[7],
		a[2]*b[4]  + a[6]*b[5]  + a[10]*b[6]  + a[14]*b[7],
		a[3]*b[4]  + a[7]*b[5]  + a[11]*b[6]  + a[15]*b[7],

		a[0]*b[8]  + a[4]*b[9]  + a[8]*b[10] + a[12]*b[11],
		a[1]*b[8]  + a[5]*b[9]  + a[9]*b[10] + a[13]*b[11],
		a[2]*b[8]  + a[6]*b[9]  + a[10]*b[10] + a[14]*b[11],
		a[3]*b[8]  + a[7]*b[9]  + a[11]*b[10] + a[15]*b[11],

		a[0]*b[12] + a[4]*b[13] + a[8]*b[14] + a[12]*b[15],
		a[1]*b[12] + a[5]*b[13] + a[9]*b[14] + a[13]*b[15],
		a[2]*b[12] + a[6]*b[13] + a[10]*b[14] + a[14]*b[15],
		a[3]*b[12] + a[7]*b[13] + a[11]*b[14] + a[15]*b[15]
	};
}

class Mat2 {
	std::array<float, 4> data;
public:
	Mat2() : data {1,0, 0,1} {}

	Mat2(	float a, float b,
			float c, float d)
			: data{ a, b, c, d} {}

	const 	float& operator[](int idx) const { return data[idx]; }

	const	Mat2 operator*	(const float f) {
		return {f*data[0], f*data[1], f*data[2], f*data[3]}; 
	}

};

inline float det(const Mat2& m) {
	return m[0]*m[3] - m[1]*m[2];
}

inline const Mat2 inv(const Mat2& m) {
	return Mat2(m[3],-m[1], -m[2],m[0]) * (1 / det(m));
}

inline const Vec2f operator*(const Mat2& m, const Vec2f& v) {
	return {m[0]*v[0] + m[1]*v[1],
			m[2]*v[0] + m[3]*v[1]};
} 

//////// TRIANGLE ////////

inline float triarea(const Vec2f p0, const Vec2f p1, const Vec2f p2) {
	const auto v01 = p1 - p0;
	const auto v02 = p2 - p0;
	return 0.5 * (v01[0]*v02[1] - v01[1]*v02[0]);
}


#endif