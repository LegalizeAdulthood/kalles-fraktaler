/*
Kalles Fraktaler 2
Copyright (C) 2013-2017 Karl Runmo
Copyright (C) 2017-2018 Claude Heiland-Allen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef KF_FLOATEXP_V16_H
#define KF_FLOATEXP_V16_H

#include "floatexp.h"

typedef int32_t     int32_t_v16     __attribute__ ((vector_size ( 64)));
typedef int64_t     int64_t_v16     __attribute__ ((vector_size (128)));
typedef double      double_v16      __attribute__ ((vector_size (128)));
typedef long double long_double_v16 __attribute__ ((vector_size (256)));

// https://stackoverflow.com/questions/40730815/gnu-c-native-vectors-how-to-broadcast-a-scalar-like-x86s-mm-set1-epi16
template <typename V, typename S> V broadcast(S x) {
  return x - (V){};
}

// https://stackoverflow.com/questions/40730815/gnu-c-native-vectors-how-to-broadcast-a-scalar-like-x86s-mm-set1-epi16
inline int64_t all(const int64_t_v16 &x) {
  return x[0] & x[1] & x[2] & x[3] & x[4] & x[5] & x[6] & x[7] & x[8] & x[9] & x[10] & x[11] & x[12] & x[13] & x[14] & x[15];
}

inline int64_t any(const int64_t_v16 &x) {
  return x[0] | x[1] | x[2] | x[3] | x[4] | x[5] | x[6] | x[7] | x[8] | x[9] | x[10] | x[11] | x[12] | x[13] | x[14] | x[15];
}

class floatexp_v16
{
public:
	double_v16 val;
	int64_t_v16 exp;

  floatexp operator[](int i) const
  {
		assert(0 <= i);
		assert(i < 16);
		return floatexp(val[i], exp[i], 0);
	}

	inline void align()
	{
		union { double_v16 d; int64_t_v16 i; } u;
		u.d = val;
		exp += ((u.i & 0x7FF0000000000000LL) >> 52) - 1023;
		u.i = (u.i & 0x800FFFFFFFFFFFFFLL) | 0x3FF0000000000000LL;
		exp = (val != 0.0) ? exp : EXP_MIN;
		val = (val != 0.0) ? u.d : 0.0;
	}

	inline explicit operator double_v16() const
	{
		union { double_v16 d; int64_t_v16 i; } u;
		u.d = val;
		u.i = (u.i & 0x800FFFFFFFFFFFFFLL) | ((exp + 1023) << 52);
		return ((val != 0.0) && (-MAX_PREC <= exp) && (exp <= MAX_PREC)) ? u.d : 0.0;
	}

	inline floatexp_v16 &abs()
	{
		val = (val < 0.0) ? -val : val;
		return *this;
	}

	inline floatexp_v16()
	{
		val = broadcast<double_v16>(0.0);
		exp = broadcast<int64_t_v16>(EXP_MIN);
	}

	inline floatexp_v16(const floatexp &f)
	{
		val = broadcast<double_v16>(f.val);
		exp = broadcast<int64_t_v16>(f.exp);
	}

	inline floatexp_v16(const floatexp_v16 &f)
	{
		val = f.val;
		exp = f.exp;
	}

	inline floatexp_v16(const double_v16 &a)
	{
		val = a;
		exp = broadcast<int64_t_v16>(0);
		align();
	}

	inline floatexp_v16(const double_v16 &a, const int64_t_v16 &e)
	{
		val = a;
		exp = e;
		align();
	}

	inline floatexp_v16(const double_v16 &a, const int64_t_v16 &e, const int &dummy)
	{
		(void) dummy;
		val = a;
		exp = e;
	}

	inline floatexp_v16(const long_double_v16 &a)
	{
		for (int i = 0; i < 16; ++i)
		{
			floatexp f(a[i]);
			val[i] = f.val;
			exp[i] = f.exp;
		}
	}

	inline floatexp_v16(int a)
	{
		floatexp f(a);
		val = broadcast<double_v16>(f.val);
		exp = broadcast<int64_t_v16>(f.exp);
	}

	inline floatexp_v16(double a)
	{
		floatexp f(a);
		val = broadcast<double_v16>(f.val);
		exp = broadcast<int64_t_v16>(f.exp);
	}

	inline floatexp_v16(long double a)
	{
		floatexp f(a);
		val = broadcast<double_v16>(f.val);
		exp = broadcast<int64_t_v16>(f.exp);
	}

	inline floatexp_v16 &operator=(const floatexp_v16 &a)
	{
		val = a.val;
		exp = a.exp;
		return *this;
	}

	inline floatexp_v16 mul2() const
	{
		return floatexp_v16(val, exp + 1, 0);
	}

	inline floatexp_v16 mul4() const
	{
		return floatexp_v16(val, exp + 2, 0);
	}
	
	inline floatexp_v16 operator+(const floatexp_v16 &a) const
	{
		floatexp_v16 r;
		auto e = exp > a.exp;
		int64_t_v16 diff = e ? exp - a.exp : a.exp - exp;
		r.exp = e ? exp : a.exp;
		r.val = e
		      ? (diff > MAX_PREC ? val : val + double_v16(floatexp_v16(a.val, -diff, 0)))
		      : (diff > MAX_PREC ? a.val : a.val + double_v16(floatexp_v16(val, -diff, 0)));             
		r.align();
		return r;
	}

	inline floatexp_v16 operator -() const
	{
		return floatexp_v16(-val, exp, 0);
	}

	inline floatexp_v16 operator-(const floatexp_v16 &a) const
	{
		return *this + (- a);
	}

	inline floatexp_v16 &operator +=(const floatexp_v16 &a)
	{
		*this = *this + a;
		return *this;
	}

	inline floatexp_v16 &operator-=(const floatexp_v16 &a)
	{
		*this = *this - a;
		return *this;
	}

	inline floatexp_v16 &operator *=(int a)
	{
		val *= a;
		align();
		return *this;
	}

	inline floatexp_v16 &operator *=(double a)
	{
		val *= a;
		align();
		return *this;
	}

	inline floatexp_v16 &operator /=(int a)
	{
		val /= a;
		align();
		return *this;
	}


	inline floatexp_v16 &operator /=(double a)
	{
		val /= a;
		align();
		return *this;
	}

	inline int64_t_v16 operator ==(const floatexp_v16 &a) const
	{
		return (exp != a.exp) ? 0 : (val == a.val);
	}

	inline int64_t_v16 operator>(const floatexp_v16 &a) const
	{
		return
		  val > 0.0 ? (a.val < 0.0 ? -1 : exp > a.exp ? -1 : exp < a.exp ?  0 : val > a.val)
		            : (a.val > 0.0 ?  0 : exp > a.exp ?  0 : exp < a.exp ? -1 : val > a.val);
	}

	inline int64_t_v16 operator<(const floatexp_v16 &a) const
	{
		return a > *this;
	}

	inline int64_t_v16 operator <=(const floatexp_v16 &a) const
	{
		return *this < a || *this == a;
	}

};

inline floatexp_v16 operator *(const floatexp_v16 &a, const floatexp_v16 &b)
{
	return floatexp_v16(a.val * b.val, a.exp + b.exp);
}

inline floatexp_v16 operator /(const floatexp_v16 &a, const floatexp_v16 &b)
{
	return floatexp_v16(a.val / b.val, a.exp - b.exp);
}

inline floatexp_v16 operator*(double a, const floatexp_v16 &b)
{
	return floatexp_v16(a) * b;
}

inline floatexp_v16 operator*(const floatexp_v16 &b, double a)
{
	return floatexp_v16(a) * b;
}

inline floatexp_v16 operator+(double a, const floatexp_v16 &b)
{
	return floatexp_v16(a) + b;
}

inline floatexp_v16 operator+(const floatexp_v16 &b, double a)
{
	return floatexp_v16(a) + b;
}

inline floatexp_v16 operator*(int a, const floatexp_v16 &b)
{
	return double(a) * b;
}

inline floatexp_v16 abs(floatexp_v16 a)
{
	return a.abs();
}

#endif
