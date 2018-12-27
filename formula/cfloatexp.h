#ifndef CFLOATEXP_H
#define CFLOATEXP_H 1

#include <math.h>
#include <stdbool.h>
#include <stdint.h>

static inline int64_t et_j_max_jj(int64_t a, int64_t b)
{
  return a > b ? a : b;
}

static inline int et_i_cmp_ff(float x, float y)
{
  return (x > y) - (x < y);
}

static inline int et_i_cmp_dd(double x, double y)
{
  return (x > y) - (x < y);
}

static inline int et_i_cmp_ll(long double x, long double y)
{
  return (x > y) - (x < y);
}

static inline int et_i_sgn_f(float x)
{
  return (x >= 0) - (x < 0);
}

static inline int et_i_sgn_d(double x)
{
  return (x >= 0) - (x < 0);
}

static inline int et_i_sgn_l(long double x)
{
  return (x >= 0) - (x < 0);
}

static inline float et_f_diffabs_ff(float c, float d) {
  if (c >= 0.0f) {
    if (c + d >= 0.0f) { return d; }
    else { return -(2.0f * c + d); }
  } else {
    if (c + d > 0.0f) { return 2.0f * c + d; }
    else { return -d; }
  }
}

static inline double et_d_diffabs_dd(double c, double d) {
  if (c >= 0.0) {
    if (c + d >= 0.0) { return d; }
    else { return -(2.0 * c + d); }
  } else {
    if (c + d > 0.0) { return 2.0 * c + d; }
    else { return -d; }
  }
}

static inline long double et_l_diffabs_ll(long double c, long double d) {
  if (c >= 0.0L) {
    if (c + d >= 0.0L) { return d; }
    else { return -(2.0L * c + d); }
  } else {
    if (c + d > 0.0L) { return 2.0L * c + d; }
    else { return -d; }
  }
}

typedef struct
{
  double val;
  int64_t exp;
} cfloatexp;

static inline cfloatexp et_e_zero(void)
{
  cfloatexp ret = { 0, 0 };
  return ret;
}

static inline cfloatexp et_e_inf(void)
{
  cfloatexp ret = { 1.0/0.0, 0 };
  return ret;
}

static inline cfloatexp et_e_ninf(void)
{
  cfloatexp ret = { -1.0/0.0, 0 };
  return ret;
}

static inline cfloatexp et_e_nan(void)
{
  cfloatexp ret = { 0.0/0.0, 0 };
  return ret;
}

static inline cfloatexp et_e_set_dj(double val, int64_t exp)
{
  if (val == 0 || isnan(val) || isinf(val))
  {
    exp = 0;
  }
  else
  {
    int e = 0;
    val = frexp(val, &e);
    exp += e;
  }
  cfloatexp ret = { val, exp };
  return ret;
}

static inline cfloatexp et_e_set_lj(long double val, int64_t exp)
{
  if (val == 0 || isnan(val) || isinf(val))
  {
    exp = 0;
  }
  else
  {
    int e = 0;
    val = frexpl(val, &e);
    exp += e;
  }
  cfloatexp ret = { (double) val, exp };
  return ret;
}

static inline cfloatexp et_e_set_d(double val)
{
  return et_e_set_dj(val, 0);
}

static inline cfloatexp et_e_set_l(long double val)
{
  return et_e_set_lj(val, 0);
}

static inline cfloatexp et_e_set_i(int val)
{
  return et_e_set_dj(val, 0);
}

static inline cfloatexp et_e_set_j(int64_t val)
{
  return et_e_set_dj(val, 0);
}

static inline double et_d_set_e(cfloatexp l)
{
  if (l.exp < INT_MIN)
    return l.val * 0.0; // underflow
  if (l.exp > INT_MAX)
    return l.val / 0.0; // overflow
  return ldexp(l.val, l.exp);
}

static inline int et_i_set_e(cfloatexp l)
{
  return et_d_set_e(l);
}

static inline long double et_l_set_e(cfloatexp l)
{
  if (l.exp < INT_MIN)
    return l.val * 0.0; // underflow
  if (l.exp > INT_MAX)
    return l.val / 0.0; // overflow
  return ldexpl(l.val, l.exp);
}

static inline double et_i_sgn_e(cfloatexp l)
{
  return et_i_sgn_d(l.val);
}

static inline cfloatexp et_e_abs_e(cfloatexp l)
{
  cfloatexp ret = { l.val < 0 ? -l.val : l.val, l.exp };
  return ret;
}

static inline cfloatexp et_e_neg_e(cfloatexp l)
{
  cfloatexp ret = { -l.val, l.exp };
  return ret;
}

static inline cfloatexp et_e_mul2_e(cfloatexp l)
{
  cfloatexp ret = { l.val, l.exp + 1 };
  return ret;
}

static inline cfloatexp et_e_mul_ee(cfloatexp l, cfloatexp r)
{
  double val = l.val * r.val;
  int64_t exp = l.exp + r.exp;
  return et_e_set_dj(val, exp);
}

static inline cfloatexp et_e_mul_ie(int l, cfloatexp r)
{
  return et_e_mul_ee(et_e_set_i(l), r);
}

static inline cfloatexp et_e_mul_ei(cfloatexp l, int r)
{
  return et_e_mul_ee(l, et_e_set_i(r));
}

static inline cfloatexp et_e_mul_de(double l, cfloatexp r)
{
  return et_e_mul_ee(et_e_set_d(l), r);
}

static inline cfloatexp et_e_mul_ed(cfloatexp l, double r)
{
  return et_e_mul_ee(l, et_e_set_d(r));
}

static inline cfloatexp et_e_mul_le(long double l, cfloatexp r)
{
  return et_e_mul_ee(et_e_set_l(l), r);
}

static inline cfloatexp et_e_mul_el(cfloatexp l, long double r)
{
  return et_e_mul_ee(l, et_e_set_l(r));
}

static inline cfloatexp et_e_div_ee(cfloatexp l, cfloatexp r)
{
  double val = l.val / r.val;
  int64_t exp = l.exp - r.exp;
  return et_e_set_dj(val, exp);
}

static inline cfloatexp et_e_div_ie(int l, cfloatexp r)
{
  return et_e_div_ee(et_e_set_i(l), r);
}

static inline cfloatexp et_e_div_ei(cfloatexp l, int r)
{
  return et_e_div_ee(l, et_e_set_i(r));
}

static inline cfloatexp et_e_div_de(double l, cfloatexp r)
{
  return et_e_div_ee(et_e_set_d(l), r);
}

static inline cfloatexp et_e_div_ed(cfloatexp l, double r)
{
  return et_e_div_ee(l, et_e_set_d(r));
}

static inline cfloatexp et_e_div_le(long double l, cfloatexp r)
{
  return et_e_div_ee(et_e_set_l(l), r);
}

static inline cfloatexp et_e_div_el(cfloatexp l, long double r)
{
  return et_e_div_ee(l, et_e_set_l(r));
}

static inline cfloatexp et_e_inv_e(cfloatexp l)
{
  double val = 1 / l.val;
  int64_t exp = -l.exp;
  return et_e_set_dj(val, exp);
}

static inline cfloatexp et_e_add_ee(cfloatexp l, cfloatexp r)
{
  if (l.val == 0) return r;
  if (r.val == 0) return l;
  int64_t e = et_j_max_jj(l.exp, r.exp);
  cfloatexp a = { l.val, l.exp - e };
  cfloatexp b = { r.val, r.exp - e };
  return et_e_set_dj(et_d_set_e(a) + et_d_set_e(b), e);
}

static inline cfloatexp et_e_add_ie(int l, cfloatexp r)
{
  return et_e_add_ee(et_e_set_i(l), r);
}

static inline cfloatexp et_e_add_ei(cfloatexp l, int r)
{
  return et_e_add_ee(l, et_e_set_i(r));
}

static inline cfloatexp et_e_add_de(double l, cfloatexp r)
{
  return et_e_add_ee(et_e_set_d(l), r);
}

static inline cfloatexp et_e_add_ed(cfloatexp l, double r)
{
  return et_e_add_ee(l, et_e_set_d(r));
}

static inline cfloatexp et_e_add_le(long double l, cfloatexp r)
{
  return et_e_add_ee(et_e_set_l(l), r);
}

static inline cfloatexp et_e_add_el(cfloatexp l, long double r)
{
  return et_e_add_ee(l, et_e_set_l(r));
}

static inline cfloatexp et_e_sub_ee(cfloatexp l, cfloatexp r)
{
  if (l.val == 0) return et_e_neg_e(r);
  if (r.val == 0) return l;
  int64_t e = et_j_max_jj(l.exp, r.exp);
  cfloatexp a = { l.val, l.exp - e };
  cfloatexp b = { r.val, r.exp - e };
  return et_e_set_dj(et_d_set_e(a) - et_d_set_e(b), e);
}

static inline cfloatexp et_e_sub_ie(int l, cfloatexp r)
{
  return et_e_sub_ee(et_e_set_i(l), r);
}

static inline cfloatexp et_e_sub_ei(cfloatexp l, int r)
{
  return et_e_sub_ee(l, et_e_set_i(r));
}

static inline cfloatexp et_e_sub_de(double l, cfloatexp r)
{
  return et_e_sub_ee(et_e_set_d(l), r);
}

static inline cfloatexp et_e_sub_ed(cfloatexp l, double r)
{
  return et_e_sub_ee(l, et_e_set_d(r));
}

static inline cfloatexp et_e_sub_le(long double l, cfloatexp r)
{
  return et_e_sub_ee(et_e_set_l(l), r);
}

static inline cfloatexp et_e_sub_el(cfloatexp l, long double r)
{
  return et_e_sub_ee(l, et_e_set_l(r));
}

static inline int et_i_cmp_ee(cfloatexp l, cfloatexp r)
{
  if (l.val == 0) return et_i_cmp_dd(0, r.val);
  if (r.val == 0) return et_i_cmp_dd(l.val, 0);
  int64_t e = et_j_max_jj(l.exp, r.exp);
  cfloatexp a = { l.val, l.exp - e };
  cfloatexp b = { r.val, r.exp - e };
  return et_i_cmp_dd(et_d_set_e(a), et_d_set_e(b));
}

static inline int et_i_cmp_ie(int l, cfloatexp r)
{
  return et_i_cmp_ee(et_e_set_i(l), r);
}

static inline int et_i_cmp_ei(cfloatexp l, int r)
{
  return et_i_cmp_ee(l, et_e_set_i(r));
}

static inline int et_i_cmp_de(double l, cfloatexp r)
{
  return et_i_cmp_ee(et_e_set_d(l), r);
}

static inline int et_i_cmp_ed(cfloatexp l, double r)
{
  return et_i_cmp_ee(l, et_e_set_d(r));
}

static inline int et_i_cmp_le(long double l, cfloatexp r)
{
  return et_i_cmp_ee(et_e_set_l(l), r);
}

static inline int et_i_cmp_el(cfloatexp l, long double r)
{
  return et_i_cmp_ee(l, et_e_set_l(r));
}

static inline cfloatexp et_e_sqrt_e(cfloatexp l)
{
  bool even = 0 == (l.exp & 1);
  return et_e_set_dj(sqrt(even ? l.val : 2 * l.val), l.exp / 2);
}

static inline cfloatexp et_e_exp_e(cfloatexp l)
{
  if (l.exp == 0)
  {
    return et_e_set_dj(exp(l.val), 0);
  }
  else if (-53 <= l.exp && l.exp <= 9)
  {
    return et_e_set_dj(exp(et_d_set_e(l)), 0);
  }
  else if (l.exp >= 61)
  {
    return l.val > 0 ? et_e_inf() : et_e_zero();
  }
  else if (l.exp < -53)
  {
    return et_e_set_dj(1, 0);
  }
  else
  {
    // 9 < e < 61
    // exp m ^ (2 ^ e) = (exp m ^ 2)^2)^2 .. e times
    cfloatexp ret = et_e_set_dj(exp(l.val), 0);
    for (int64_t i = 0; i < l.exp; ++i)
    {
      ret = et_e_mul_ee(ret, ret);
    }
    return ret;
  }
}

static inline cfloatexp et_e_log_e(cfloatexp l)
{
  static const double log_2 = 0.6931471805599453;
  return et_e_set_dj(log(l.val) + l.exp * log_2, 0);
}

static inline cfloatexp et_e_expm1_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp > 9)
  {
    return et_e_sub_ee(et_e_exp_e(l), et_e_set_d(1));
  }
  else
  {
    return et_e_set_d(expm1(et_d_set_e(l)));
  }
}

static inline cfloatexp et_e_log1p_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp >= 1024)
  {
    return et_e_log_e(et_e_add_ee(et_e_set_d(1), l));
  }
  else
  {
    return et_e_set_d(log1p(et_d_set_e(l)));
  }
}

static inline cfloatexp et_e_sin_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp < 1024)
  {
    return et_e_set_d(sin(et_d_set_e(l)));
  }
  else
  {
    // sin(2 x) = 2 sin(x) cos(x)
    // cos(2 x) = cos(x)^2 - sin(x)^2
    double s = sin(l.val);
    double c = cos(l.val);
    {
      for (int64_t e = 0; e < l.exp; ++e)
      {
        double u = 2 * s * c;
        double v = c * c - s * s;
        s = u;
        c = v;
      }
    }
    return et_e_set_d(s);
  }
}

static inline cfloatexp et_e_cos_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return et_e_set_d(1);
  }
  else if (l.exp < 1024)
  {
    return et_e_set_d(cos(et_d_set_e(l)));
  }
  else
  {
    // sin(2 x) = 2 sin(x) cos(x)
    // cos(2 x) = cos(x)^2 - sin(x)^2
    double s = sin(l.val);
    double c = cos(l.val);
    {
      for (int64_t e = 0; e < l.exp; ++e)
      {
        double u = 2 * s * c;
        double v = c * c - s * s;
        s = u;
        c = v;
      }
    }
    return et_e_set_d(c);
  }
}

static inline cfloatexp et_e_tan_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp < 1024)
  {
    return et_e_set_d(tan(et_d_set_e(l)));
  }
  else
  {
    // sin(2 x) = 2 sin(x) cos(x)
    // cos(2 x) = cos(x)^2 - sin(x)^2
    double s = sin(l.val);
    double c = cos(l.val);
    {
      for (int64_t e = 0; e < l.exp; ++e)
      {
        double u = 2 * s * c;
        double v = c * c - s * s;
        s = u;
        c = v;
      }
    }
    return et_e_div_ee(et_e_set_d(s), et_e_set_d(c));
  }
}

static inline cfloatexp et_e_sinh_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp <= 0)
  {
    return et_e_set_d(sinh(et_d_set_e(l)));
  }
  else
  {
    cfloatexp ret = et_e_sub_ee(et_e_exp_e(l), et_e_exp_e(et_e_neg_e(l)));
    ret.exp -= 1;
    return ret;
  }
}

static inline cfloatexp et_e_cosh_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return et_e_set_d(1);
  }
  else if (l.exp <= 0)
  {
    return et_e_set_d(cosh(et_d_set_e(l)));
  }
  else
  {
    cfloatexp ret = et_e_add_ee(et_e_exp_e(l), et_e_exp_e(et_e_neg_e(l)));
    ret.exp -= 1;
    return ret;
  }
}

static inline cfloatexp et_e_tanh_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp < 1024)
  {
    return et_e_set_d(tanh(et_d_set_e(l)));
  }
  else
  {
    return et_e_set_i(et_i_sgn_d(l.val));
  }
}

static inline cfloatexp et_e_asin_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else
  {
    return et_e_set_d(asin(et_d_set_e(l)));
  }
}

static inline cfloatexp et_e_acos_e(cfloatexp l)
{
  return et_e_set_d(acos(et_d_set_e(l)));
}

static inline cfloatexp et_e_atan_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.val < 0)
  {
    return et_e_neg_e(et_e_atan_e(et_e_neg_e(l)));
  }
  else if (l.exp >= 1024)
  {
    // pi/2
    return et_e_sub_ee(et_e_set_d(1.5707963267948966), et_e_atan_e(et_e_inv_e(l)));
  }
  else
  {
    return et_e_set_d(atan(et_d_set_e(l)));
  }
}

static inline cfloatexp et_e_asinh_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else if (l.exp < 1024)
  {
    return et_e_set_d(asinh(et_d_set_e(l)));
  }
  else if (l.val < 0)
  {
    return et_e_neg_e(et_e_asinh_e(et_e_neg_e(l)));
  }
  else
  {
    cfloatexp r = l;
    r.exp += 1;
    return et_e_log_e(r);
  }
}

static inline cfloatexp et_e_acosh_e(cfloatexp l)
{
  if (l.exp < 1024)
  {
    return et_e_set_d(acosh(et_d_set_e(l)));
  }
  else if (l.val < 0)
  {
    return et_e_acosh_e(et_e_neg_e(l));
  }
  else
  {
    cfloatexp r = l;
    r.exp += 1;
    return et_e_log_e(r);
  }
}

static inline cfloatexp et_e_atanh_e(cfloatexp l)
{
  if (l.exp <= -1021)
  {
    return l;
  }
  else
  {
    return et_e_set_d(atanh(et_d_set_e(l)));
  }
}

static inline cfloatexp mpfr_get_fe(const mpfr_t value, mpfr_rnd_t rnd)
{
	signed long int e = 0;
	double l = mpfr_get_d_2exp(&e, value, rnd);
	return et_e_set_dj(l, e);
}

static inline void mpfr_set_fe(mpfr_t value, cfloatexp fe, mpfr_rnd_t rnd)
{
	mpfr_set_d(value, fe.val, rnd);
	if (fe.exp >= 0)
	{
		mpfr_mul_2ui(value, value, fe.exp, rnd);
	}
	else
	{
		mpfr_div_2ui(value, value, -fe.exp, rnd);
	}
}

static inline int et_i_cmp_ri(const mpfr_t l, int r)
{
  return mpfr_cmp_si(l, r);
}

static inline cfloatexp et_e_diffabs_ee(cfloatexp c, cfloatexp d)
{
  if (c.val >= 0) if (et_e_add_ee(c, d).val >= 0) return d; else return et_e_neg_e(et_e_add_ee(et_e_mul2_e(c), d));
  else if (et_e_add_ee(c, d).val > 0) return et_e_add_ee(et_e_mul2_e(c), d); else return et_e_neg_e(d);
}

static inline void mpfr_mul_fe(mpfr_t o, const mpfr_t l, cfloatexp r, mpfr_rnd_t rnd)
{
  mpfr_mul_d(o, l, r.val, rnd);
	if (r.exp >= 0)
	{
		mpfr_mul_2ui(o, o, r.exp, rnd);
	}
	else
	{
		mpfr_div_2ui(o, o, -r.exp, rnd);
	}
}

static inline void mpfr_mul_ld(mpfr_t o, const mpfr_t l, long double r, mpfr_rnd_t rnd)
{
  return mpfr_mul_fe(o, l, et_e_set_l(r), rnd);
}

#define b_and_bb(o,l,r) ((o)=(l)&&(r))
#define b_eq_id(o,l,r) ((o)=(l)==(r))
#define b_eq_ie(o,l,r) ((o)=et_i_cmp_ie((l),(r))==0)
#define b_eq_if(o,l,r) ((o)=(l)==(r))
#define b_eq_ii(o,l,r) ((o)=(l)==(r))
#define b_eq_il(o,l,r) ((o)=(l)==(r))
#define b_le_dd(o,l,r) ((o)=(l)<=(r))
#define b_le_ee(o,l,r) ((o)=et_i_cmp_ee((l),(r))<=0)
#define b_le_ff(o,l,r) ((o)=(l)<=(r))
#define b_le_id(o,l,r) ((o)=(l)<=(r))
#define b_le_ie(o,l,r) ((o)=et_i_cmp_ie((l),(r))<=0)
#define b_le_if(o,l,r) ((o)=(l)<=(r))
#define b_le_il(o,l,r) ((o)=(l)<=(r))
#define b_le_ll(o,l,r) ((o)=(l)<=(r))
#define b_lt_dd(o,l,r) ((o)=(l)<(r))
#define b_lt_ee(o,l,r) ((o)=et_i_cmp_ee((l),(r))<0)
#define b_lt_ff(o,l,r) ((o)=(l)<(r))
#define b_lt_id(o,l,r) ((o)=(l)<(r))
#define b_lt_ie(o,l,r) ((o)=et_i_cmp_ie((l),(r))<0)
#define b_lt_if(o,l,r) ((o)=(l)<(r))
#define b_lt_ii(o,l,r) ((o)=(l)<(r))
#define b_lt_il(o,l,r) ((o)=(l)<(r))
#define b_lt_ll(o,l,r) ((o)=(l)<(r))
#define b_lt_ri(o,l,r) ((o)=et_i_cmp_ri((l),(r))<0)
#define b_or_bb(o,l,r) ((o)=(l)||(r))
#define b_read_Vi(o,l,r) ((o)=(l)[(r)])
#define b_set_b(o,l) ((o)=(l))
#define b_set_i(o,l) ((o)=(l))
#define d_add_dd(o,l,r) ((o)=(l)+(r))
#define d_inv_d(o,l) ((o)=1/(l))
#define d_mul_dd(o,l,r) ((o)=(l)*(r))
#define d_mul_id(o,l,r) ((o)=(l)*(r))
#define d_neg_d(o,l) ((o)=-(l))
#define d_read_Di(o,l,r) ((o)=(l)[(r)])
#define d_set_d(o,l) ((o)=(l))
#define d_set_r(o,l) ((o)=mpfr_get_d((l),MPFR_RNDN))
#define d_sqrt_d(o,l) ((o)=sqrt((l)))
#define v_write_Did(p,l,r) ((p)[(l)]=(r))
#define e_add_ee(o,l,r) ((o)=et_e_add_ee((l),(r)))
#define e_inv_e(o,l) ((o)=et_e_inv_e((l)))
#define e_mul_ee(o,l,r) ((o)=et_e_mul_ee((l),(r)))
#define e_mul_ie(o,l,r) ((o)=et_e_mul_ie((l),(r)))
#define e_neg_e(o,l) ((o)=et_e_neg_e((l)))
#define e_read_Ei(o,l,r) ((o)=(l)[(r)])
#define e_set_d(o,l) ((o)=et_e_set_d((l)))
#define e_set_e(o,l) ((o)=(l))
#define e_set_i(o,l) ((o)=et_e_set_j((l)))
#define e_set_r(o,l) ((o)=mpfr_get_fe((l),MPFR_RNDN))
#define e_sqrt_e(o,l) ((o)=et_e_sqrt_e((l)))
#define v_write_Eie(p,l,r) ((p)[(l)]=(r))
#define f_add_ff(o,l,r) ((o)=(l)+(r))
#define f_inv_f(o,l) ((o)=1/(l))
#define f_mul_ff(o,l,r) ((o)=(l)*(r))
#define f_mul_if(o,l,r) ((o)=(l)*(r))
#define f_neg_f(o,l) ((o)=-(l))
#define f_read_Fi(o,l,r) ((o)=(l)[(r)])
#define f_set_f(o,l) ((o)=(l))
#define f_set_r(o,l) ((o)=mpfr_get_flt((l),MPFR_RNDN))
#define f_sqrt_f(o,l) ((o)=sqrtf((l)))
#define v_write_Fif(p,l,r) ((p)[(l)]=(r))
#define i_add_ii(o,l,r) ((o)=(l)+(r))
#define i_eq_ii(o,l,r) ((o)=(l)==(r))
#define i_ifte_bii(o,l,r,x) ((o)=(l)?(r):(x))
#define i_neg_i(o,l) ((o)=-(l))
#define i_read_Ii(o,l,r) ((o)=(l)[(r)])
#define i_read_Vi(o,l,r) ((o)=(l)[(r)])
#define i_set_i(o,l) ((o)=(l))
#define I_set_V(o,l) ((o)=(int*)(l))
#define i_sgn_i(o,l) ((o)=(l)>=0?1:-1)
#define i_sgn_r(o,l) (i_sgn_i((o),mpfr_sgn(l)))
#define v_write_Iii(p,l,r) ((p)[(l)]=(r))
#define l_add_ll(o,l,r) ((o)=(l)+(r))
#define l_inv_l(o,l) ((o)=1/(l))
#define l_mul_il(o,l,r) ((o)=(l)*(r))
#define l_mul_ll(o,l,r) ((o)=(l)*(r))
#define l_neg_l(o,l) ((o)=-(l))
#define l_read_Li(o,l,r) ((o)=(l)[(r)])
#define l_set_l(o,l) ((o)=(l))
#define l_set_r(o,l) ((o)=mpfr_get_ld((l),MPFR_RNDN))
#define l_sqrt_l(o,l) ((o)=sqrtl((l)))
#define v_write_Lil(p,l,r) ((p)[(l)]=(r))
#define r_abs_r(o,l) (mpfr_abs((o),(l),MPFR_RNDN))
#define r_add_ir(o,l,r) (mpfr_add_si((o),(r),(l),MPFR_RNDN))
#define r_add_ri(o,l,r) (mpfr_add_si((o),(l),(r),MPFR_RNDN))
#define r_add_rr(o,l,r) (mpfr_add((o),(l),(r),MPFR_RNDN))
#define r_sub_rr(o,l,r) (mpfr_sub((o),(l),(r),MPFR_RNDN))
#define r_expm1_r(o,l) (mpfr_expm1((o),(l),MPFR_RNDN))
#define r_inv_r(o,l) (mpfr_si_div((o),1,(l),MPFR_RNDN))
#define r_log1p_r(o,l) (mpfr_log1p((o),(l),MPFR_RNDN))
#define r_mul_dr(o,l,r) (mpfr_mul_d((o),(r),(l),MPFR_RNDN))
#define r_mul_ir(o,l,r) (mpfr_mul_si((o),(r),(l),MPFR_RNDN))
#define f_mul2_f(o,l) ((o)=(l)+(l))
#define d_mul2_d(o,l) ((o)=(l)+(l))
#define l_mul2_l(o,l) ((o)=(l)+(l))
#define e_mul2_e(o,l) ((o)=et_e_mul2_e((l)))
#define r_mul2_r(o,l) (mpfr_mul_2ui((o),(l),1,MPFR_RNDN))
#define r_mul_rd(o,l,r) (mpfr_mul_d((o),(l),(r),MPFR_RNDN))
#define r_mul_rr(o,l,r) ((l)==(r)?mpfr_sqr((o),(l),MPFR_RNDN): mpfr_mul((o),(l),(r),MPFR_RNDN))
#define f_sqr_f(o,l) ((o)=(l)*(l))
#define d_sqr_d(o,l) ((o)=(l)*(l))
#define l_sqr_l(o,l) ((o)=(l)*(l))
#define e_sqr_e(o,l) ((o)=et_e_mul_ee((l),(l)))
#define r_sqr_r(o,l) (mpfr_sqr((o),(l),MPFR_RNDN))
#define r_neg_r(o,l) (mpfr_neg((o),(l),MPFR_RNDN))
#define r_set_d(o,l) (mpfr_set_d((o),(l),MPFR_RNDN))
#define r_set_e(o,l) (mpfr_set_fe((o),(l),MPFR_RNDN))
#define r_set_f(o,l) (mpfr_set_flt((o),(l),MPFR_RNDN))
#define r_set_l(o,l) (mpfr_set_ld((o),(l),MPFR_RNDN))
#define r_set_r(o,l) (mpfr_set((o),(l),MPFR_RNDN))
#define r_sqrt_r(o,l) (mpfr_sqrt((o),(l),MPFR_RNDN))
#define f_abs_f(o,l) ((o)=((l)<0?-(l):(l)))
#define i_sgn_f(o,l) ((o)=((l)>=0?1:-1))
#define f_mul_fi(o,l,r) ((o)=(l)*(r))
#define d_abs_d(o,l) ((o)=((l)<0?-(l):(l)))
#define i_sgn_d(o,l) ((o)=((l)>=0?1:-1))
#define d_mul_di(o,l,r) ((o)=(l)*(r))
#define l_abs_l(o,l) ((o)=((l)<0?-(l):(l)))
#define i_sgn_l(o,l) ((o)=((l)>=0?1:-1))
#define l_mul_li(o,l,r) ((o)=(l)*(r))
#define e_abs_e(o,l) ((o)=et_e_abs_e((l)))
#define i_sgn_e(o,l) ((o)=et_i_sgn_e((l)))
#define e_mul_ei(o,l,r) ((o)=et_e_mul_ei((l),(r)))
#define f_diffabs_ff(o,l,r) ((o)=et_f_diffabs_ff((l),(r)))
#define d_diffabs_dd(o,l,r) ((o)=et_d_diffabs_dd((l),(r)))
#define l_diffabs_ll(o,l,r) ((o)=et_l_diffabs_ll((l),(r)))
#define e_diffabs_ee(o,l,r) ((o)=et_e_diffabs_ee((l),(r)))
#define r_mul_ri(o,l,r) (mpfr_mul_si((o),(l),(r),MPFR_RNDN))
#define r_mul_fr(o,l,r) (mpfr_mul_d((o),(r),(l),MPFR_RNDN))
#define r_mul_lr(o,l,r) (mpfr_mul_ld((o),(r),(l),MPFR_RNDN))
#define r_mul_er(o,l,r) (mpfr_mul_fe((o),(r),(l),MPFR_RNDN))

#endif
