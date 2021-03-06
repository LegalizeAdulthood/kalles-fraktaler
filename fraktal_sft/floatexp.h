/*
Kalles Fraktaler 2
Copyright (C) 2013-2017 Karl Runmo
Copyright (C) 2017-2020 Claude Heiland-Allen

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

#ifndef KF_FLOATEXP_H
#define KF_FLOATEXP_H

#include <math.h>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "CFixedFloat.h"

#define MAX_PREC 1020
// this has two fewer 0 than you might expect, this is to give headroom for
// avoiding overflow in + and other functions. it is the exponent for 0.0
#define EXP_MIN (-0x80000000000000LL)

struct floatexp
{
public:
	floatexp() = default; // POD
	floatexp(const floatexp &a) = default; // POD
	floatexp(floatexp &&a) = default; // POD
	floatexp &operator=(const floatexp &a) = default; // POD
	floatexp &operator=(floatexp &&a) = default; // POD
	~floatexp() = default; // POD;

	double val;
	int64_t exp;

	inline void align() noexcept
	{
		if (val != 0)
		{
			union { double d; int64_t i; } u;
			u.d = val;
			exp += ((u.i & 0x7FF0000000000000LL) >> 52) - 1023;
			u.i = (u.i & 0x800FFFFFFFFFFFFFLL) | 0x3FF0000000000000LL;
			val = u.d;
		}
		else
		{
			val = 0;
			exp = EXP_MIN;
		}
	}
	inline floatexp &abs() noexcept
	{
		if(val<0)
			val=-val;
		return *this;
	}
	inline void initFromDouble(double a) noexcept
	{
		val=a;
		exp=0;
		align();
	}
	inline void initFromLongDouble(long double a) noexcept
	{
		using std::frexp;
		int e = 0;
		a = frexp(a, &e);
		val=double(a);
		exp=e;
		align();
	}
	inline double setExp(double newval,int64_t newexp) const noexcept
	{
//		int64_t tmpval = (*((int64_t*)&newval) & 0x800FFFFFFFFFFFFF) | ((newexp+1023)<<52);
//		memcpy(&newval,&tmpval,sizeof(double));
//		return newval;
		union { double d; int64_t i; } u;
		u.d = newval;
		u.i = (u.i & 0x800FFFFFFFFFFFFFLL) | ((newexp + 1023) << 52);
		newval = u.d;
		return newval;
	}

	inline floatexp(int a) noexcept
	{
		initFromDouble(a);
	}
	inline floatexp(int64_t a) noexcept
	{
		initFromDouble(a);
	}
	inline floatexp(double a) noexcept
	{
		initFromDouble(a);
	}
	inline floatexp(double a, int64_t e) noexcept
	{
		val = a;
		exp = e;
		align();
	}
	inline floatexp(double a, int64_t e, int dummy) noexcept
	{
		(void) dummy;
		val = a;
		exp = e;
	}
	inline floatexp(long double a) noexcept
	{
		initFromLongDouble(a);
	}

	inline floatexp &operator =(int a) noexcept
	{
		initFromDouble((double)a);
		return *this;
	}
	inline floatexp &operator =(double a) noexcept
	{
		initFromDouble(a);
		return *this;
	}
	inline floatexp &operator =(long double a) noexcept
	{
		initFromLongDouble(a);
		return *this;
	}
	inline floatexp operator *(const floatexp &a) const noexcept
	{
		floatexp r;
		r.val = a.val*val;
		r.exp = a.exp+exp;
		r.align();
		return r;
	}
	inline floatexp operator /(const floatexp &a) const noexcept
	{
		floatexp r;
		r.val = val/a.val;
		r.exp = exp - a.exp;
		r.align();
		return r;
	}
	__attribute__ ((warn_unused_result))
	inline floatexp mul2() const noexcept
	{
		floatexp r;
		r.val = val;
		r.exp = exp + 1;
		return r;
	}
	__attribute__ ((warn_unused_result))
	inline floatexp mul4() const noexcept
	{
		floatexp r;
		r.val = val;
		r.exp = exp + 2;
		return r;
	}
	inline floatexp operator +(const floatexp &a) const noexcept
	{
		floatexp r;
		int64_t diff;
		if(exp>a.exp){
			diff = exp-a.exp;
			r.exp = exp;
			if(diff>MAX_PREC)
				r.val=val;
			else{
				double aval = setExp(a.val,-diff);
				r.val = val+aval;
			}
		}
		else{
			diff = a.exp-exp;
			r.exp = a.exp;
			if(diff>MAX_PREC)
				r.val=a.val;
			else{
				double aval = setExp(val,-diff);
				r.val = a.val+aval;
			}
		}
		r.align();
		return r;
	}
	inline floatexp operator -() const noexcept
	{
		floatexp r=*this;
		r.val=-r.val;
		return r;
	}
	inline floatexp &operator +=(const floatexp &a) noexcept
	{
		*this = *this+a;
		return *this;
	}
	inline floatexp operator -(const floatexp &a) const noexcept
	{
		floatexp r;
		int64_t diff;
		if(exp>a.exp){
			diff = exp-a.exp;
			r.exp = exp;
			if(diff>MAX_PREC)
				r.val = val;
			else{
				double aval = setExp(a.val,-diff);
				r.val = val-aval;
			}
		}
		else{
			diff = a.exp-exp;
			r.exp = a.exp;
			if(diff>MAX_PREC)
				r.val=-a.val;
			else{
				double aval = setExp(val,-diff);
				r.val = aval-a.val;
			}
		}
		r.align();
		return r;
	}
	inline floatexp &operator -=(const floatexp &a) noexcept
	{
		*this = *this-a;
		return *this;
	}
	inline bool operator >(const floatexp &a) const noexcept
	{
		if(val>0){
			if(a.val<0)
				return true;
			if(exp>a.exp)
				return true;
			else if(exp<a.exp)
				return false;
			return val>a.val;
		}
		else{
			if(a.val>0)
				return false;
			if(exp>a.exp)
				return false;
			else if(exp<a.exp)
				return true;
			return val>a.val;
		}
	}
	inline bool operator <(const floatexp &a) const noexcept
	{
		if(val>0){
			if(a.val<0)
				return false;
			if(exp>a.exp)
				return false;
			else if(exp<a.exp)
				return true;
			return val<a.val;
		}
		else{
			if(a.val>0)
				return true;
			if(exp>a.exp)
				return true;
			else if(exp<a.exp)
				return false;
			return val<a.val;
		}
	}
	inline bool operator <=(const floatexp &a) const noexcept
	{
		return (*this<a || *this==a);
	}
	inline bool operator >=(const floatexp &a) const noexcept
	{
		return (*this>a || *this==a);
	}
	inline bool operator <=(const int a) const noexcept
	{
		return (*this<a || *this==a);
	}
	inline bool operator ==(const floatexp &a) const noexcept
	{
		if(exp!=a.exp)
			return false;
		return val==a.val;
	}
	inline bool iszero() const noexcept
	{
		return (val==0 && exp==0);
	}
	inline double todouble() const noexcept
	{
		if(exp<-MAX_PREC || exp>MAX_PREC)
			return 0;
		return setExp(val,exp);
	}
	inline explicit operator double () const noexcept
	{
		return todouble();
	}
	inline double todouble(int nScaling) const noexcept
	{
		if(!nScaling)
			return todouble();
		floatexp ret = *this;
		while(nScaling>9){
			ret.val*=1e10;
			ret.align();
			nScaling-=10;
		}
		while(nScaling>2){
			ret.val*=1e3;
			ret.align();
			nScaling-=3;
		}
		while(nScaling>0){
			ret.val*=1e1;
			ret.align();
			nScaling--;
		}
		while(nScaling<-9){
			ret.val/=1e10;
			ret.align();
			nScaling+=10;
		}
		while(nScaling<-2){
			ret.val/=1e3;
			ret.align();
			nScaling+=3;
		}
		while(nScaling<0){
			ret.val/=1e1;
			ret.align();
			nScaling++;
		}
		if(ret.exp<-MAX_PREC || ret.exp>MAX_PREC)
			return 0;
		return setExp(ret.val,ret.exp);
	}
	inline floatexp &operator /=(double a) noexcept
	{
		val/=a;
		align();
		return *this;
	}
	inline floatexp &operator *=(double a) noexcept
	{
		val*=a;
		align();
		return *this;
	}
	inline floatexp &operator *=(floatexp a) noexcept
	{
		return *this = *this * a;
	}
	inline floatexp &operator *=(long double a) noexcept
	{
		return *this *= floatexp(a);
	}

	inline floatexp &operator =(const CFixedFloat &a) noexcept
	{
		signed long int e = 0;
		val = mpfr_get_d_2exp(&e, a.m_f.backend().data(), MPFR_RNDN);
		exp = e;
		align();
		return *this;
	}
	inline floatexp(const CFixedFloat &a) noexcept
	{
		*this = a;
	}
	inline floatexp(const CDecNumber &a) noexcept
	{
		signed long int e = 0;
		val = mpfr_get_d_2exp(&e, a.m_dec.backend().data(), MPFR_RNDN);
		exp = e;
		align();
	}
	inline void ToFixedFloat(CFixedFloat &a) const noexcept
	{
		a = val;
		if (exp > int64_t(UINT_MAX))
		{
			a = 1.0 / 0.0;
		}
		else if (exp < -int64_t(UINT_MAX))
		{
			a = 0.0;
		}
		else
		{
			a = val;
			if (exp >= 0)
				mpfr_mul_2ui(a.m_f.backend().data(), a.m_f.backend().data(), exp, MPFR_RNDN);
			else
				mpfr_div_2ui(a.m_f.backend().data(), a.m_f.backend().data(), -exp, MPFR_RNDN);
		}
	}
	inline explicit operator CFixedFloat() const noexcept
	{
		CFixedFloat a;
		ToFixedFloat(a);
		return a;
	}

	inline floatexp setLongDouble(long double a) noexcept
	{
		int e = 0;
		val = double(std::frexp(a, &e));
		exp = e;
		align();
		return *this;
	}
	inline long double toLongDouble() const noexcept
	{
		if (val == 0.0L)
			return 0.0L;
		if (exp >= INT_MAX)
			return (val / 0.0L); // infinity
		if (exp <= INT_MIN)
			return (val * 0.0L); // zero
		return std::ldexp((long double) val, exp);
	}
	inline long double toLongDouble(int nScaling) const noexcept
	{
		if(!nScaling)
			return toLongDouble();
		floatexp ret = *this;
		// FIXME risky to go higher than this? 1e300 might be ok?
		while(nScaling>99){
			ret.val*=1e100;
			ret.align();
			nScaling-=100;
		}
		while(nScaling>29){
			ret.val*=1e30;
			ret.align();
			nScaling-=30;
		}
		while(nScaling>9){
			ret.val*=1e10;
			ret.align();
			nScaling-=10;
		}
		while(nScaling>2){
			ret.val*=1e3;
			ret.align();
			nScaling-=3;
		}
		while(nScaling){
			ret.val*=1e1;
			ret.align();
			nScaling--;
		}
		return ret.toLongDouble();
	}
	inline explicit operator long double () const noexcept
	{
		return toLongDouble();
	}

  inline std::string toString(int digits = 0) const noexcept
  {
		/*
		  f = val 2^exp
		  log10 f = log10 val + exp log10 2
		  e10 = floor(log10 f)
		  d10 = 10^(log10 f - e10)
		  d10 \in [1, 10)
		*/
		double lf = std::log10(std::abs(val)) + exp * std::log10(2.0);
		int64_t e10 = int64_t(std::floor(lf));
		double d10 = std::pow(10, lf - e10) * ((val > 0) - (val < 0));
		if (val == 0) { d10 = 0; e10 = 0; }
		std::ostringstream os; os
		  << std::setprecision(digits ? digits : (std::numeric_limits<double>::digits10 + 1))
		  << std::fixed
		  << d10 << 'E' << e10;
		return os.str();
	}
} __attribute__((packed));

inline std::ostream& operator<<(std::ostream& a, const floatexp& b) noexcept
{
	return a << b.toString();
}

inline floatexp operator*(double a, floatexp b) noexcept
{
	return floatexp(a) * b;
}

inline floatexp operator*(floatexp b, double a) noexcept
{
	return floatexp(a) * b;
}

inline floatexp operator*(long double a, floatexp b) noexcept
{
	return floatexp(a) * b;
}

inline floatexp operator*(floatexp b, long double a) noexcept
{
	return floatexp(a) * b;
}

inline floatexp operator+(double a, floatexp b) noexcept
{
	return floatexp(a) + b;
}

inline floatexp operator+(floatexp b, double a) noexcept
{
	return floatexp(a) + b;
}

inline floatexp operator*(int a, floatexp b) noexcept
{
	return double(a) * b;
}

inline floatexp operator*(floatexp b, int a) noexcept
{
	return double(a) * b;
}

inline floatexp operator-(int a, floatexp b) noexcept
{
	return floatexp(a) - b;
}

inline floatexp abs(floatexp a) noexcept
{
	return a.abs();
}

inline floatexp sqrt(floatexp a) noexcept
{
  return floatexp
    ( std::sqrt((a.exp & 1) ? 2.0 * a.val : a.val)
    , (a.exp & 1) ? (a.exp - 1) / 2 : a.exp / 2
    );
}

inline floatexp log(floatexp a) noexcept
{
	return floatexp(std::log(a.val) + std::log(2.0) * a.exp);
}

inline floatexp log2(floatexp a) noexcept
{
	return floatexp(std::log2(a.val) + a.exp);
}

inline double sqr(double a) noexcept
{
	return a * a;
}

inline long double sqr(long double a) noexcept
{
	return a * a;
}

inline floatexp sqr(floatexp a) noexcept
{
	return a * a;
}

template <typename T> T pow(T x, uint64_t n) noexcept
{
	switch (n)
	{
		case 0: return T(1);
		case 1: return x;
		case 2: return sqr(x);
		case 3: return x * sqr(x);
		case 4: return sqr(sqr(x));
		case 5: return x * sqr(sqr(x));
		case 6: return sqr(x * sqr(x));
		case 7: return x * sqr(x * sqr(x));
		case 8: return sqr(sqr(sqr(x)));
		default:
		{
		  T y(1);
		  while (n > 1)
		  {
		    if (n & 1)
		      y *= x;
		    x = sqr(x);
		    n >>= 1;
			}
		  return x * y;
		}
	}
}

inline bool isnan(const floatexp &a) noexcept
{
	return isnan(a.val);
}

inline bool isinf(const floatexp &a) noexcept
{
	return isinf(a.val);
}

inline floatexp infnan_to_zero(const floatexp &a) noexcept
{
	return isinf(a.val) ? floatexp(copysign(1e30, a.val)) : isnan(a.val) ? floatexp(0) : a;
}

inline floatexp exp(floatexp a) noexcept
{
	using std::exp;
	using std::ldexp;
	using std::pow;
  if (-53 <= a.exp && a.exp <= 8) return floatexp(exp(ldexp(a.val, a.exp)));
  if (61 <= a.exp) return floatexp(a.val > 0.0 ? a.val / 0.0 : 0.0);
  if (a.exp < -53) return floatexp(1.0);
  return pow(floatexp(exp(a.val)), 1ULL << a.exp);
}

inline floatexp expm1(floatexp a) noexcept
{
	using std::expm1;
	using std::ldexp;
  if (a.exp <= -1020) return a;
  if (8 <= a.exp) return exp(a) - 1;
  return floatexp(expm1(ldexp(a.val, a.exp)));
}

inline floatexp sin(floatexp a) noexcept
{
	using std::sin;
	if (a.exp <= -1020) return a;
	return floatexp(sin(ldexp(a.val, a.exp)));
}

inline floatexp cos(floatexp a) noexcept
{
	using std::cos;
	if (a.exp <= -1020) return floatexp(1.0);
	return floatexp(cos(ldexp(a.val, a.exp)));
}

inline floatexp diffabs(const floatexp &c, const floatexp &d) noexcept
{
  const floatexp cd = c + d;
  const floatexp c2d = c.mul2() + d;
  return c.val >= 0.0 ? cd.val >= 0.0 ? d : -c2d : cd.val > 0.0 ? c2d : -d;
}

inline floatexp mpfr_get_fe(const mpfr_t value) noexcept
{
	signed long int e = 0;
	double l = mpfr_get_d_2exp(&e, value, MPFR_RNDN);
	return floatexp(l, e);
}

inline void mpfr_set_fe(mpfr_t value, floatexp fe) noexcept
{
	mpfr_set_d(value, fe.val, MPFR_RNDN);
	if (fe.exp >= 0)
	{
		mpfr_mul_2ui(value, value, fe.exp, MPFR_RNDN);
	}
	else
	{
		mpfr_div_2ui(value, value, -fe.exp, MPFR_RNDN);
	}
}

inline long double mpfr_get_ld(const mpfr_t value) noexcept
{
	using std::ldexp;
	signed long int e = 0;
	long double l = mpfr_get_ld_2exp(&e, value, MPFR_RNDN);
	if (l == 0.0L)
		return 0.0L;
	if (e >= INT_MAX)
		return l / 0.0L;
	if (e <= INT_MIN)
		return l * 0.0L;
	l = ldexp(l, e);
	return l;
}

#endif
