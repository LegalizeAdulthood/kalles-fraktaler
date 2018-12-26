/*
Kalles Fraktaler 2
Copyright (C) 2013-2017 Karl Runmo
Copyright (C) 2017-2018 Claude Heiland-Allen

incorporating components derived from:

et -- escape time fractals
Copyright (C) 2018 Claude Heiland-Allen

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

#ifndef KF_FORMULA_H
#define KF_FORMULA_H 1

#include <cassert>
#include <cmath>

#include <mpfr.h>

#include "../../fraktal_sft/floatexp.h"

using std::abs;
using std::floor;
using std::sqrt;
using std::exp;
using std::log;
using std::sin;
using std::cos;
using std::tan;
using std::sinh;
using std::cosh;
using std::tanh;
using std::asin;
using std::acos;
using std::atan;
using std::asinh;
using std::acosh;
using std::atanh;

typedef int f_plainf(int,float,float*,float,float,float,float,float*,volatile int*);
typedef int f_plain(int,double,double*,double,double,double,double,double*,volatile int*);
typedef int f_plainl(int,long double,long double*,long double,long double,long double,long double,long double*,volatile int*);
typedef int f_plainfe(int,floatexp,floatexp*,floatexp,floatexp,floatexp,floatexp,floatexp*,volatile int*);

typedef int f_referencef(int,int,float*,float*,float*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,float,float,float,float,float,float,float,int,int*,float*,float*);
typedef int f_reference(int,int,double*,double*,double*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,double,double,double,double,double,double,double,int,int*,double*,double*);
typedef int f_referencel(int,int,long double*,long double*,long double*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,long double,long double,long double,long double,long double,long double,long double,int,int*,long double*,long double*);
typedef int f_referencefe(int,int,floatexp*,floatexp*,floatexp*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,floatexp,floatexp,floatexp,floatexp,floatexp,floatexp,floatexp,int,int*,floatexp*,floatexp*);

typedef int f_referenceDf(int,int,float*,float*,float*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,float,float,float,float,float,float,float,int,int*,float*,float*,float*,float*);
typedef int f_referenceD(int,int,double*,double*,double*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,double,double,double,double,double,double,double,int,int*,double*,double*,double*,double*);
typedef int f_referenceDl(int,int,long double*,long double*,long double*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,long double,long double,long double,long double,long double,long double,long double,int,int*,long double*,long double*,long double*,long double*);
typedef int f_referenceDfe(int,int,floatexp*,floatexp*,floatexp*,volatile int*,int*,int*,int*,mpfr_t,mpfr_t,floatexp,floatexp,floatexp,floatexp,floatexp,floatexp,floatexp,int,int*,floatexp*,floatexp*,floatexp*,floatexp*);

typedef int f_perturbationf(int,int,float*,float*,float*,int*,float*,float*,int*,float,int,int,float,float,float,float,float*,float*,float,float);
typedef int f_perturbation(int,int,double*,double*,double*,int*,double*,double*,int*,double,int,int,double,double,double,double,double*,double*,double,double);
typedef int f_perturbationl(int,int,long double*,long double*,long double*,int*,long double*,long double*,int*,long double,int,int,long double,long double,long double,long double,long double*,long double*,long double,long double);
typedef int f_perturbationfe(int,int,floatexp*,floatexp*,floatexp*,int*,floatexp*,floatexp*,int*,floatexp,int,int,floatexp,floatexp,floatexp,floatexp,floatexp*,floatexp*,floatexp,floatexp);

typedef int f_perturbationDf(int,int,float*,float*,float*,int*,float*,float*,int*,float,int,int,float,float,float,float,float*,float*,float,float,float*,float*);
typedef int f_perturbationD(int,int,double*,double*,double*,int*,double*,double*,int*,double,int,int,double,double,double,double,double*,double*,double,double,double*,double*);
typedef int f_perturbationDl(int,int,long double*,long double*,long double*,int*,long double*,long double*,int*,long double,int,int,long double,long double,long double,long double,long double*,long double*,long double,long double,long double*,long double*);
typedef int f_perturbationDfe(int,int,floatexp*,floatexp*,floatexp*,int*,floatexp*,floatexp*,int*,floatexp,int,int,floatexp,floatexp,floatexp,floatexp,floatexp*,floatexp*,floatexp,floatexp,floatexp*,floatexp*);

typedef int f_period_tri(int,double,double,double,mpfr_t,mpfr_t,mpfr_t,volatile int*);
typedef int f_period_jsk(int,double,double,double,mpfr_t,mpfr_t,mpfr_t,double*,volatile int*);
typedef int f_newton(int,int,double,double,mpfr_t,mpfr_t,volatile int*);
typedef int f_size(int,double,double,mpfr_t,mpfr_t,mpfr_t,double*,volatile int*);
typedef int f_skew(int,double,double,mpfr_t,mpfr_t,int,double*,volatile int*);
typedef int f_domain_size(int,double,double,mpfr_t,mpfr_t,mpfr_t,volatile int*);

#ifndef KF_MAIN
static f_plainf plainf;
static f_plain plain;
static f_plainl plainl;
static f_plainfe plainfe;
static f_referencef referencef;
static f_reference reference;
static f_referencel referencel;
static f_referencefe referencefe;
static f_referenceDf referenceDf;
static f_referenceD referenceD;
static f_referenceDl referenceDl;
static f_referenceDfe referenceDfe;
static f_perturbationf perturbationf;
static f_perturbation perturbation;
static f_perturbationl perturbationl;
static f_perturbationfe perturbationfe;
static f_perturbationDf perturbationDf;
static f_perturbationD perturbationD;
static f_perturbationDl perturbationDl;
static f_perturbationDfe perturbationDfe;
static f_period_tri period_tri;
static f_period_jsk period_jsk;
static f_newton newton;
static f_size size;
static f_skew skew;
static f_domain_size domain_size;
//static const char name[];
//static const char source[];
#endif

#define MAGIC ((int)(0xC01dCaf3))
#define SIZE ((int)(sizeof(struct formula)))
#define VERSION 7

struct formula
{
  int magic;
  int ssize;
  int version;
  f_reference *reference;
  f_referencel *referencel;
  f_referencefe *referencefe;
  f_referenceD *referenceD;
  f_referenceDl *referenceDl;
  f_referenceDfe *referenceDfe;
  f_perturbation *perturbation;
  f_perturbationl *perturbationl;
  f_perturbationfe *perturbationfe;
  f_perturbationD *perturbationD;
  f_perturbationDl *perturbationDl;
  f_perturbationDfe *perturbationDfe;
  f_period_tri *period_tri;
  f_period_jsk *period_jsk;
  f_newton *newton;
  f_size *size;
  f_skew *skew;
};

#define FORMULA(name,source,degree) \
extern "C" { __declspec(dllexport) struct formula et = \
{ MAGIC, SIZE, VERSION \
, &reference, &referencel, &referencefe \
, &referenceD, &referenceDl, &referenceDfe \
, &perturbation, &perturbationl, &perturbationfe \
, &perturbationD, &perturbationDl, &perturbationDfe \
, &period_tri, &period_jsk, &newton, &size, &skew \
}; }

#endif
