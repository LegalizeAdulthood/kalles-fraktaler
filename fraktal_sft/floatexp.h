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

#ifndef KF_FLOATEXP_H
#define KF_FLOATEXP_H

#include <math.h>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#define MAX_PREC 1020
// this has two fewer 0 than you might expect, this is to give headroom for
// avoiding overflow in + and other functions. it is the exponent for 0.0
#define EXP_MIN (-0x80000000000000LL)

#define isnanl isnan
#define isinfl isinf

#include "../formula/cfloatexp.h"

struct floatexp;
floatexp operator+(const floatexp &a, const floatexp &b); 
floatexp operator-(const floatexp &a, const floatexp &b); 
floatexp operator*(const floatexp &a, const floatexp &b);
floatexp operator/(const floatexp &a, const floatexp &b);


struct floatexp : public cfloatexp
{
	// construct
  inline floatexp(double a, int64_t e, int dummy) { val = a; exp = e; (void) dummy; };
  inline floatexp() : floatexp(0, 0, 0) { };
  inline floatexp(const cfloatexp &x) : floatexp(x.val, x.exp, 0) { };
  inline floatexp(const floatexp &x) : floatexp(x.val, x.exp, 0) { };
  inline floatexp(int x) : floatexp(et_e_set_j(x)) { };
  inline floatexp(double x) : floatexp(et_e_set_d(x)) { };
  inline floatexp(long double x) : floatexp(et_e_set_l(x)) { };
  inline floatexp(double a, int64_t e) : floatexp(et_e_set_dj(a, e)) { };
  // convert
	inline operator cfloatexp () const { return *this; };
	inline explicit operator int () const { return et_i_set_e(*this); };
	inline explicit operator double () const { return et_d_set_e(*this); };
	inline explicit operator long double () const { return et_l_set_e(*this); };
	inline explicit operator std::string () const 
  {
		/*
		  f = val 2^exp
		  log10 f = log10 val + exp log10 2
		  e10 = floor(log10 f)
		  d10 = 10^(log10 f - e10)
		  d10 \in [1, 10)
		*/
		if (isnan(val)) return "nan";
		if (isinf(val) && val > 0) return "inf";
		if (isinf(val) && val < 0) return "-inf";
		double lf = std::log10(std::abs(val)) + exp * std::log10(2.0);
		int64_t e10 = std::floor(lf);
		double d10 = std::pow(10, lf - e10) * ((val > 0) - (val < 0));
		if (val == 0) { d10 = 0; e10 = 0; }
		std::ostringstream os; os
		  << std::setprecision(std::numeric_limits<double>::digits10 + 1)
		  << std::fixed
		  << d10 << 'E' << e10;
		return os.str();
	}
	// assign
	inline floatexp &operator=(const floatexp &x) { val = x.val; exp = x.exp; return *this; };
	inline floatexp &operator+=(const floatexp &x) { return *this = *this + x; };
	inline floatexp &operator-=(const floatexp &x) { return *this = *this - x; };
	inline floatexp &operator*=(const floatexp &x) { return *this = *this * x; };
	inline floatexp &operator/=(const floatexp &x) { return *this = *this / x; };
  // convert with scaling
	inline double todouble(int nScaling) const
	{
		if(!nScaling)
			return double(*this);
		floatexp ret = *this;
		while(nScaling>9){
			ret*=1e10;
			nScaling-=10;
		}
		while(nScaling>2){
			ret*=1e3;
			nScaling-=3;
		}
		while(nScaling>0){
			ret*=1e1;
			nScaling--;
		}
		while(nScaling<-9){
			ret/=1e10;
			nScaling+=10;
		}
		while(nScaling<-2){
			ret/=1e3;
			nScaling+=3;
		}
		while(nScaling<0){
			ret/=1e1;
			nScaling++;
		}
		return double(ret);
	}
	inline long double toLongDouble(int nScaling) const
	{
		if(!nScaling)
			return (long double)(*this);
		floatexp ret = *this;
		// FIXME risky to go higher than this? 1e300 might be ok?
		while(nScaling>99){
			ret*=1e100;
			nScaling-=100;
		}
		while(nScaling>29){
			ret*=1e30;
			nScaling-=30;
		}
		while(nScaling>9){
			ret*=1e10;
			nScaling-=10;
		}
		while(nScaling>2){
			ret*=1e3;
			nScaling-=3;
		}
		while(nScaling){
			ret*=1e1;
			nScaling--;
		}
		return (long double)(ret);
	}

};

inline floatexp operator-(const floatexp &a) { return et_e_neg_e(a); } 

inline floatexp operator+(const floatexp &a, const floatexp &b) { return et_e_add_ee(a, b); } 
inline floatexp operator-(const floatexp &a, const floatexp &b) { return et_e_sub_ee(a, b); } 
inline floatexp operator*(const floatexp &a, const floatexp &b) { return et_e_mul_ee(a, b); } 
inline floatexp operator/(const floatexp &a, const floatexp &b) { return et_e_div_ee(a, b); } 

inline floatexp operator+(const floatexp &a, int b) { return et_e_add_ei(a, b); } 
inline floatexp operator-(const floatexp &a, int b) { return et_e_sub_ei(a, b); } 
inline floatexp operator*(const floatexp &a, int b) { return et_e_mul_ei(a, b); } 
inline floatexp operator/(const floatexp &a, int b) { return et_e_div_ei(a, b); } 

inline floatexp operator+(int a, const floatexp &b) { return et_e_add_ie(a, b); } 
inline floatexp operator-(int a, const floatexp &b) { return et_e_sub_ie(a, b); } 
inline floatexp operator*(int a, const floatexp &b) { return et_e_mul_ie(a, b); } 
inline floatexp operator/(int a, const floatexp &b) { return et_e_div_ie(a, b); } 

inline floatexp operator+(const floatexp &a, double b) { return et_e_add_ed(a, b); } 
inline floatexp operator-(const floatexp &a, double b) { return et_e_sub_ed(a, b); } 
inline floatexp operator*(const floatexp &a, double b) { return et_e_mul_ed(a, b); } 
inline floatexp operator/(const floatexp &a, double b) { return et_e_div_ed(a, b); } 

inline floatexp operator+(double a, const floatexp &b) { return et_e_add_de(a, b); } 
inline floatexp operator-(double a, const floatexp &b) { return et_e_sub_de(a, b); } 
inline floatexp operator*(double a, const floatexp &b) { return et_e_mul_de(a, b); } 
inline floatexp operator/(double a, const floatexp &b) { return et_e_div_de(a, b); } 

inline floatexp operator+(const floatexp &a, long double b) { return et_e_add_el(a, b); } 
inline floatexp operator-(const floatexp &a, long double b) { return et_e_sub_el(a, b); } 
inline floatexp operator*(const floatexp &a, long double b) { return et_e_mul_el(a, b); } 
inline floatexp operator/(const floatexp &a, long double b) { return et_e_div_el(a, b); } 

inline floatexp operator+(long double a, const floatexp &b) { return et_e_add_le(a, b); } 
inline floatexp operator-(long double a, const floatexp &b) { return et_e_sub_le(a, b); } 
inline floatexp operator*(long double a, const floatexp &b) { return et_e_mul_le(a, b); } 
inline floatexp operator/(long double a, const floatexp &b) { return et_e_div_le(a, b); } 

inline int cmp(const floatexp &a, const floatexp &b) { return et_i_cmp_ee(a, b); }
inline bool operator< (const floatexp &a, const floatexp &b) { return cmp(a, b) <  0; }
inline bool operator<=(const floatexp &a, const floatexp &b) { return cmp(a, b) <= 0; }
inline bool operator> (const floatexp &a, const floatexp &b) { return cmp(a, b) >  0; }
inline bool operator>=(const floatexp &a, const floatexp &b) { return cmp(a, b) >= 0; }
inline bool operator==(const floatexp &a, const floatexp &b) { return cmp(a, b) == 0; }
inline bool operator!=(const floatexp &a, const floatexp &b) { return cmp(a, b) != 0; }

inline int cmp(const floatexp &a, int b) { return et_i_cmp_ei(a, b); }
inline bool operator< (const floatexp &a, int b) { return cmp(a, b) <  0; }
inline bool operator<=(const floatexp &a, int b) { return cmp(a, b) <= 0; }
inline bool operator> (const floatexp &a, int b) { return cmp(a, b) >  0; }
inline bool operator>=(const floatexp &a, int b) { return cmp(a, b) >= 0; }
inline bool operator==(const floatexp &a, int b) { return cmp(a, b) == 0; }
inline bool operator!=(const floatexp &a, int b) { return cmp(a, b) != 0; }

inline int cmp(int a, const floatexp &b) { return et_i_cmp_ie(a, b); }
inline bool operator< (int a, const floatexp &b) { return cmp(a, b) <  0; }
inline bool operator<=(int a, const floatexp &b) { return cmp(a, b) <= 0; }
inline bool operator> (int a, const floatexp &b) { return cmp(a, b) >  0; }
inline bool operator>=(int a, const floatexp &b) { return cmp(a, b) >= 0; }
inline bool operator==(int a, const floatexp &b) { return cmp(a, b) == 0; }
inline bool operator!=(int a, const floatexp &b) { return cmp(a, b) != 0; }

inline int cmp(const floatexp &a, double b) { return et_i_cmp_ed(a, b); }
inline bool operator< (const floatexp &a, double b) { return cmp(a, b) <  0; }
inline bool operator<=(const floatexp &a, double b) { return cmp(a, b) <= 0; }
inline bool operator> (const floatexp &a, double b) { return cmp(a, b) >  0; }
inline bool operator>=(const floatexp &a, double b) { return cmp(a, b) >= 0; }
inline bool operator==(const floatexp &a, double b) { return cmp(a, b) == 0; }
inline bool operator!=(const floatexp &a, double b) { return cmp(a, b) != 0; }

inline int cmp(double a, const floatexp &b) { return et_i_cmp_de(a, b); }
inline bool operator< (double a, const floatexp &b) { return cmp(a, b) <  0; }
inline bool operator<=(double a, const floatexp &b) { return cmp(a, b) <= 0; }
inline bool operator> (double a, const floatexp &b) { return cmp(a, b) >  0; }
inline bool operator>=(double a, const floatexp &b) { return cmp(a, b) >= 0; }
inline bool operator==(double a, const floatexp &b) { return cmp(a, b) == 0; }
inline bool operator!=(double a, const floatexp &b) { return cmp(a, b) != 0; }

inline int cmp(const floatexp &a, long double b) { return et_i_cmp_el(a, b); }
inline bool operator< (const floatexp &a, long double b) { return cmp(a, b) <  0; }
inline bool operator<=(const floatexp &a, long double b) { return cmp(a, b) <= 0; }
inline bool operator> (const floatexp &a, long double b) { return cmp(a, b) >  0; }
inline bool operator>=(const floatexp &a, long double b) { return cmp(a, b) >= 0; }
inline bool operator==(const floatexp &a, long double b) { return cmp(a, b) == 0; }
inline bool operator!=(const floatexp &a, long double b) { return cmp(a, b) != 0; }

inline int cmp(long double a, const floatexp &b) { return et_i_cmp_le(a, b); }
inline bool operator< (long double a, const floatexp &b) { return cmp(a, b) <  0; }
inline bool operator<=(long double a, const floatexp &b) { return cmp(a, b) <= 0; }
inline bool operator> (long double a, const floatexp &b) { return cmp(a, b) >  0; }
inline bool operator>=(long double a, const floatexp &b) { return cmp(a, b) >= 0; }
inline bool operator==(long double a, const floatexp &b) { return cmp(a, b) == 0; }
inline bool operator!=(long double a, const floatexp &b) { return cmp(a, b) != 0; }

// misc
inline floatexp abs(const floatexp &a) { return et_e_abs_e(a); }
inline floatexp sqrt(const floatexp &a) { return et_e_sqrt_e(a); }
inline floatexp log(const floatexp &a) { return et_e_log_e(a); }
inline floatexp mul2(const floatexp &a) { return et_e_mul2_e(a); }

inline int sgn(float a) { return (a >= 0) - (a < 0); }
inline int sgn(double a) { return (a >= 0) - (a < 0); }
inline int sgn(long double a) { return (a >= 0) - (a < 0); }
inline int sgn(const floatexp &a) { return sgn(a.val); }

inline float diffabs(float a, float b) { return et_f_diffabs_ff(a, b); }
inline double diffabs(double a, double b) { return et_d_diffabs_dd(a, b); }
inline long double diffabs(long double a, long double b) { return et_l_diffabs_ll(a, b); }
inline floatexp diffabs(const floatexp &a, const floatexp &b) { return et_e_diffabs_ee(a, b); }

#endif
