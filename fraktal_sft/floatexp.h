#ifndef __FLOATEXP_H__
#define __FLOATEXP_H__

#include <math.h>
#include <stdint.h>
#include "CFixedFloat.h"

#define MAX_PREC 1020
// this has two fewer 0 than you might expect, this is to give headroom for
// avoiding overflow in + and other functions. it is the exponent for 0.0
#define EXP_MIN (-0x80000000000000LL)

#define _ALIGN_(val,exp) exp += ((*((int64_t*)&val) & 0x7FF0000000000000LL)>>52) - 1023; *((int64_t*)&val) = (*((int64_t*)&val) & 0x800FFFFFFFFFFFFFLL) | 0x3FF0000000000000LL;
class floatexp
{
public:
	double val;
	int64_t exp;
	inline void align()
	{
		if (val != 0)
		{
			_ALIGN_(val,exp)
		}
		else
		{
			val = 0;
			exp = EXP_MIN;
		}
	}
	inline floatexp &abs()
	{
		if(val<0)
			val=-val;
		return *this;
	}
	inline void initFromDouble(double a)
	{
		val=a;
		exp=0;
		align();
	}
	inline double setExp(double newval,int64_t newexp) const
	{
//		int64_t tmpval = (*((int64_t*)&newval) & 0x800FFFFFFFFFFFFF) | ((newexp+1023)<<52);
//		memcpy(&newval,&tmpval,sizeof(double));
//		return newval;
		*((int64_t*)&newval) = (*((int64_t*)&newval) & 0x800FFFFFFFFFFFFFLL) | ((newexp+1023)<<52);
		return newval;
	}
	inline floatexp()
	{
		val = 0;
		exp = EXP_MIN;
	}
	inline floatexp(double a)
	{
		initFromDouble(a);
	}
	inline floatexp(double a, int64_t e)
	{
		val = a;
		exp = e;
		align();
	}
	inline floatexp(double a, int64_t e, int dummy)
	{
		val = a;
		exp = e;
	}
	inline floatexp &operator =(const floatexp &a)
	{
		val=a.val;
		exp=a.exp;
		return *this;
	}
	inline floatexp &operator =(int a)
	{	
		initFromDouble((double)a);
		return *this;
	}
	inline floatexp &operator =(double a)
	{
		initFromDouble(a);
		return *this;
	}
	inline floatexp operator *(const floatexp &a) const
	{
		floatexp r;
		r.val = a.val*val;
		r.exp = a.exp+exp;
		r.align();
		return r;
	}
	inline floatexp operator /(const floatexp &a) const
	{
		floatexp r;
		r.val = val/a.val;
		r.exp = exp - a.exp;
		r.align();
		return r;
	}
	inline floatexp &mul2()
	{
		exp++;
		return *this;
	}
	inline floatexp &mul4()
	{
		exp+=2;
		return *this;
	}
	inline floatexp operator +(const floatexp &a) const
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
	inline floatexp operator -() const
	{
		floatexp r=*this;
		r.val=-r.val;
		return r;
	}
	inline floatexp &operator +=(const floatexp &a)
	{
		*this = *this+a;
		return *this;
	}
	inline floatexp operator -(const floatexp &a) const
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
	inline floatexp &operator -=(const floatexp &a)
	{
		*this = *this-a;
		return *this;
	}
	inline bool operator >(const floatexp &a) const
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
	inline bool operator <(const floatexp &a) const
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
	inline bool operator <=(const int a) const
	{
		floatexp aa(a);
		return (*this<a || *this==a);
	}
	inline bool operator ==(const floatexp &a) const
	{
		if(exp!=a.exp)
			return false;
		return val==a.val;
	}
	inline bool iszero() const
	{
		return (val==0 && exp==0);
	}
	inline double todouble() const
	{
		if(exp<-MAX_PREC || exp>MAX_PREC)
			return 0;
		return setExp(val,exp);
	}
	inline double todouble(int nScaling) const
	{
		if(!nScaling)
			return todouble();
		floatexp ret = *this;
		while(nScaling>9){
			ret.val*=1e10;
			_ALIGN_(ret.val,ret.exp)
			nScaling-=10;
		}
		while(nScaling>2){
			ret.val*=1e3;
			_ALIGN_(ret.val,ret.exp)
			nScaling-=3;
		}
		while(nScaling){
			ret.val*=1e1;
			_ALIGN_(ret.val,ret.exp)
			nScaling--;
		}
		if(ret.exp<-MAX_PREC || ret.exp>MAX_PREC)
			return 0;
		return setExp(ret.val,ret.exp);
	}
	inline floatexp &operator /=(double a)
	{
		val/=a;
		align();
		return *this;
	}
	inline floatexp &operator *=(double a)
	{
		val*=a;
		align();
		return *this;
	}

	inline floatexp &operator =(const CFixedFloat &a)
	{
		signed long int e = 0;
		val = mpf_get_d_2exp(&e, a.m_f.backend().data());
		if ((mpf_sgn(a.m_f.backend().data()) >= 0) != (val >= 0)) val = -val; // workaround GMP bug
		exp = e;
		align();
		return *this;
	}
	inline void ToFixedFloat(CFixedFloat &a) const
	{
		a = val;
		if (exp >= 0)
			mpf_mul_2exp(a.m_f.backend().data(), a.m_f.backend().data(), exp);
		else
			mpf_div_2exp(a.m_f.backend().data(), a.m_f.backend().data(), -exp);
	}

	inline floatexp setLongDouble(long double a)
	{
		int e = 0;
		val = std::frexp(a, &e);
		exp = e;
		align();
		return *this;
	}
	inline long double toLongDouble() const
	{
		return std::ldexp((long double) val, exp);
	}

};

inline floatexp operator*(double a, floatexp b)
{
	return floatexp(a) * b;
}

inline floatexp operator*(int a, floatexp b)
{
	return double(a) * b;
}

inline floatexp abs(floatexp a)
{
	return a.abs();
}

inline floatexp sqrt(floatexp a)
{
  return floatexp
    ( std::sqrt((a.exp & 1) ? 2.0 * a.val : a.val)
    , (a.exp & 1) ? (a.exp - 1) / 2 : a.exp / 2
    );
}

inline floatexp mpf_get_fe(const mpf_t value)
{
	signed long int e = 0;
	double l = mpf_get_d_2exp(&e, value);
	if ((mpf_sgn(value) >= 0) != (l >= 0)) l = -l; // workaround GMP bug
	return floatexp(l, e);
}

inline void mpf_set_fe(mpf_t value, floatexp fe)
{
	mpf_set_d(value, fe.val);
	if (fe.exp >= 0)
	{
		mpf_mul_2exp(value, value, fe.exp);
	}
	else
	{
		mpf_div_2exp(value, value, -fe.exp);
	}
}

#endif //__FLOATEXP_H__
