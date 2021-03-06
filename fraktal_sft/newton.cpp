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

#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <float.h>
#include "resource.h"
#include "complex.h"
#include "fraktal_sft.h"
#include "CDecNumber.h"
#include "main.h"
#include "../common/barrier.h"
#include "../common/StringVector.h"
#include "newton.h"
#include "hybrid.h"
#include <iostream>
#include <string>

#define KF_MAIN 1
#include "../formula/generated/formula.h"
#undef abs
#undef sgn

static const struct formula *get_formula(int type, int power)
{
  std::string name = "formula_" + std::to_string(type) + "_" + std::to_string(power);
  const struct formula *f = (const struct formula *) GetProcAddress(GetModuleHandle(nullptr), name.c_str());
  if (f) return f;
  name = "formula_" + std::to_string(type);
  f = (const struct formula *) GetProcAddress(GetModuleHandle(nullptr), name.c_str());
  return f;
}

extern CFraktalSFT g_SFT;
extern HICON g_hIcon;

static volatile int running = 0;
bool g_bNewtonRunning = false;
bool g_bNewtonStop = false;
static bool g_bNewtonExit = false;
bool g_bJustDidNewton = false;
bool g_bUseBallPeriod = true;

static floatexp g_fNewtonDelta2[2];
static int g_nNewtonETA = 0;
static int64_t g_iterations = 0;

// progress reporting threads

struct progress_t
{
	int counters[4];
	volatile bool running;
	HWND hWnd;
	HANDLE hDone;
};

static DWORD WINAPI ThPeriodProgress(progress_t *progress)
{
	while (progress->running)
	{
		int limit = progress->counters[0];
		int iter = progress->counters[1];
		char status[100];
		snprintf(status, 100, "Period %d (%d%%)", iter, (int) (iter * 100.0 / limit));
		SetDlgItemText(progress->hWnd, IDC_EDIT1, status);
		Sleep(250);
	}
	SetEvent(progress->hDone);
	return 0;
}

static DWORD WINAPI ThNewtonProgress(progress_t *progress)
{
	while (progress->running)
	{
		int eta = progress->counters[0];
		int step = progress->counters[1];
		int period = progress->counters[2];
		int iter = progress->counters[3];
		char status[100];
		snprintf(status, 100, "NR %d/%d (%d%%)", step, eta > 0 ? step + eta : 0, (int) (iter * 100.0 / period));
		SetDlgItemText(progress->hWnd, IDC_EDIT1, status);
		Sleep(250);
	}
	SetEvent(progress->hDone);
	return 0;
}

static DWORD WINAPI ThSizeProgress(progress_t *progress)
{
	while (progress->running)
	{
		int period = progress->counters[0];
		int iter = progress->counters[1];
		char status[100];
		snprintf(status, 100, "Size %d%%", (int) (iter * 100.0 / period));
		SetDlgItemText(progress->hWnd, IDC_EDIT1, status);
		Sleep(250);
	}
	SetEvent(progress->hDone);
	return 0;
}

static DWORD WINAPI ThSkewProgress(progress_t *progress)
{
	while (progress->running)
	{
		int period = progress->counters[0];
		int iter = progress->counters[1];
		char status[100];
		snprintf(status, 100, "Skew %d%%", (int) (iter * 100.0 / period));
		SetDlgItemText(progress->hWnd, IDC_EDIT1, status);
		Sleep(250);
	}
	SetEvent(progress->hDone);
	return 0;
}





static int newtonETA(const floatexp &delta0, const floatexp &delta1, const floatexp &epsilon)
{
  floatexp l0 = log(delta0);
  floatexp l1 = log(delta1);
  floatexp e  = log(epsilon);
  // (e - l0) = 1 (l1 - l0) + 2 (l1 - l0) + 4 (l1 - l0) + ... 2^N (l1 - l0)
  // (e - l0) / (l1 - l0) = sum_0^N 2^i
  floatexp N = log2((e - l0) / (l1 - l0));
  return double(N);
}

static std::string g_szRe;
static std::string g_szIm;
static std::string g_szZoom;
static char g_szProgress[128];
static int g_nMinibrotPos=0;
static double g_nMinibrotFactor=1.0;
static std::string g_sMinibrotSourceZoom;

#define flyttyp CDecNumber

const complex<flyttyp> _2(2,0);
const complex<flyttyp> _1(1,0);

static inline int sgn(const flyttyp &z) {
  if (z > 0) { return  1; }
  if (z < 0) { return -1; }
  return 0;
}

static inline bool odd(int a) {
  return a & 1;
}

static inline floatexp cabs2(const complex<floatexp> &z) {
  return z.m_r * z.m_r + z.m_i * z.m_i;
}
static inline flyttyp cabs2(const complex<flyttyp> &z) {
  return z.m_r * z.m_r + z.m_i * z.m_i;
}
static inline bool isfinite(const flyttyp &a){
	(void) a;
	return true;
/*	if (a <= DBL_MAX && a >= -DBL_MAX)
		return true;
	return false;*/
}
static inline bool cisfinite(const complex<flyttyp> &z) {
  return isfinite(z.m_r) && isfinite(z.m_i);
}

struct BallPeriodCommon
{
  HWND hWnd;
  barrier *barrier;
  volatile bool *stop;
  bool *haveperiod;
  int64_t *period;
  int64_t maxperiod;
  mpfr_t cr, ci, zr, zi, zr2, zi2, zri, t;
  floatexp zradius;
};

struct BallPeriod
{
  int threadid;
  HANDLE hDone;
  BallPeriodCommon *c;
};

static DWORD WINAPI ThBallPeriod(BallPeriod *b)
{
  SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
  int t = b->threadid;
  barrier *barrier = b->c->barrier;
  volatile bool *stop = b->c->stop;
  bool *haveperiod = b->c->haveperiod;
  int64_t *period = b->c->period;
  char szStatus[300];
  floatexp r = b->c->zradius;
  complex<floatexp> z(0.0, 0.0);
  complex<floatexp> dz(0.0, 0.0);
  floatexp rz = 0.0;
  floatexp rdz = 0.0;
  floatexp Ei = 0.0;
  floatexp rr = 0.0;
  uint32_t last = 0;
  if (t == 0)
  {
    wsprintf(szStatus,"Finding period, 0...");
    SetDlgItemText(b->c->hWnd,IDC_EDIT1,szStatus);
    last = GetTickCount();
  }
  for (int64_t i = 1; i < b->c->maxperiod; ++i)
  {
    if (t == 0)
    {
      if(i%100==0){
	uint32_t now = GetTickCount();
	if (now - last > 250){
	  snprintf(szStatus,200,"Finding period, %lld...",i);
	  SetDlgItemText(b->c->hWnd,IDC_EDIT1,szStatus);
	  last = now;
	}
      }
      Ei = rdz * rdz + (2 * rz + r * (2 * rdz + r * Ei)) * Ei;
      dz = 2 * z * dz + complex<floatexp>(1.0, 0.0);
      z = complex<floatexp>(mpfr_get_fe(b->c->zr), mpfr_get_fe(b->c->zi));
      rdz = abs(dz);
      rz = abs(z);
      rr = r * (rdz + r * Ei);
      if (rz - rr > 2) // outside everything
      {
	*haveperiod = true;
	*period = 0;
      }
      if (rz - rr <= 0) // surrounds 0
      {
	*haveperiod = true;
	*period = i;
      }
    }
    switch (t)
    {
      case 0:
      {
        mpfr_sqr(b->c->zr2, b->c->zr, MPFR_RNDN);
        mpfr_sqr(b->c->zi2, b->c->zi, MPFR_RNDN);
	break;
      }
      case 1:
      {
        mpfr_mul(b->c->zri, b->c->zr, b->c->zi, MPFR_RNDN);
	break;
      }
    }
    if (barrier->wait(stop)) break;
    if (*haveperiod) break;
    switch (t)
    {
      case 0:
        mpfr_sub(b->c->t, b->c->zr2, b->c->zi2, MPFR_RNDN);
        mpfr_add(b->c->zr, b->c->t, b->c->cr, MPFR_RNDN);
	break;
      case 1:
        mpfr_mul_2ui(b->c->zri, b->c->zri, 1, MPFR_RNDN);
        mpfr_add(b->c->zi, b->c->zri, b->c->ci, MPFR_RNDN);
	break;
    }
    if (barrier->wait(stop)) break;
  }
  if (t == 0)
  {
    if (*period == 0)
    {
      *haveperiod = false;
    }
  }
  SetEvent(b->hDone);
  mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
  return 0;
}

static int64_t ball_period_do(const complex<flyttyp> &center, flyttyp radius, int64_t maxperiod,int &steps,HWND hWnd) {
  radius = flyttyp(4)/radius;
  mp_bitcnt_t bits = mpfr_get_prec(center.m_r.m_dec.backend().data());
  barrier bar(2);
  bool haveperiod = false;
  int64_t period = 0;
  HANDLE hDone[2];
  // prepare threads
  BallPeriod ball[2];
  BallPeriodCommon c;
  c.maxperiod = maxperiod;
  c.hWnd = hWnd;
  c.barrier = &bar;
  c.stop = &g_bNewtonStop;
  c.haveperiod = &haveperiod;
  c.period = &period;
  c.maxperiod = maxperiod;
  mpfr_init2(c.cr, bits); mpfr_set(c.cr, center.m_r.m_dec.backend().data(), MPFR_RNDN);
  mpfr_init2(c.ci, bits); mpfr_set(c.ci, center.m_i.m_dec.backend().data(), MPFR_RNDN);
  mpfr_init2(c.zr, bits); mpfr_set(c.zr, center.m_r.m_dec.backend().data(), MPFR_RNDN);
  mpfr_init2(c.zi, bits); mpfr_set(c.zi, center.m_i.m_dec.backend().data(), MPFR_RNDN);
  mpfr_init2(c.zr2, bits);
  mpfr_init2(c.zi2, bits);
  mpfr_init2(c.zri, bits);
  mpfr_init2(c.t, bits);
  c.zradius = mpfr_get_fe(radius.m_dec.backend().data());
  for (int t = 0; t < 2; ++t)
  {
    ball[t].threadid = t;
    ball[t].hDone = hDone[t] = CreateEvent(NULL, 0, 0, NULL);
    ball[t].c = &c;
  }
  // spawn threads
  for (int i = 0; i < 2; i++)
  {
    DWORD dw;
    HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)ThBallPeriod, (LPVOID)&ball[i], 0, &dw);
    CloseHandle(hThread);
  }
  // wait for threads to complete
  WaitForMultipleObjects(2, hDone, TRUE, INFINITE);
  for (int i = 0; i < 2; i++)
  {
    CloseHandle(hDone[i]);
  }
  mpfr_clear(c.cr);
  mpfr_clear(c.ci);
  mpfr_clear(c.zr);
  mpfr_clear(c.zi);
  mpfr_clear(c.zr2);
  mpfr_clear(c.zi2);
  mpfr_clear(c.zri);
  mpfr_clear(c.t);
  steps = period;
  return haveperiod ? period : 0;
}

#if 0
struct BoxPeriod
{
  int threadid;
  barrier *barrier;
  volatile bool *stop;
  int maxperiod;
  mpfr_t zr, zi, zr2, zi2, zri, cr, ci, t, *z2r, *z2i;
  int *crossing[4];
  int *haveperiod;
  int *period;
  HANDLE hDone;
  HWND hWnd;
};
static DWORD WINAPI ThBoxPeriod(BoxPeriod *b)
{
  SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
  int t = b->threadid;
  barrier *barrier = b->barrier;
  volatile bool *stop = b->stop;
  int *haveperiod = b->haveperiod;
  int *period = b->period;
  char szStatus[300];
  uint32_t last = 0;
  if (t == 0)
  {
    wsprintf(szStatus,"Finding period, 0...");
    SetDlgItemText(b->hWnd,IDC_EDIT1,szStatus);
    last = GetTickCount();
  }
  for (int i = 1; i < b->maxperiod; ++i)
  {
    // *b->crossing[t] = crosses_positive_real_axis(&b->zr, &b->zi, b->z2r, b->z2i);
    bool crossing = false;
    if (mpfr_sgn(b->zi) != mpfr_sgn(*b->z2i))
    {
      // d = b - a;
      mpfr_sub(b->zr2, *b->z2r, b->zr, MPFR_RNDN);
      mpfr_sub(b->zi2, *b->z2i, b->zi, MPFR_RNDN);
      int fs = mpfr_sgn(b->zi2);
      // t = cross(d, a);
      mpfr_mul(b->zri, b->zi2, b->zr, MPFR_RNDN);
      mpfr_mul(b->zi2, b->zr2, b->zi, MPFR_RNDN);
      mpfr_sub(b->t, b->zri, b->zi2, MPFR_RNDN);
      int ft = mpfr_sgn(b->t);
      crossing = fs == ft;
    }
    *b->crossing[t] = crossing;
    if (barrier->wait(stop)) break;
    if (t == 0)
    {
      if(i%100==0){
	uint32_t now = GetTickCount();
	if (now - last > 250){
	  wsprintf(szStatus,"Finding period, %d...",i);
	  SetDlgItemText(b->hWnd,IDC_EDIT1,szStatus);
	  last = now;
	}
      }
      int surround = 0;
      for (int s = 0; s < 4; ++s)
      {
	surround += *b->crossing[s];
      }
      if (surround & 1)
      {
	*haveperiod = true;
	*period = i;
      }
    }
    // z = z * z + c
    mpfr_sqr(b->zr2, b->zr, MPFR_RNDN);
    mpfr_sqr(b->zi2, b->zi, MPFR_RNDN);
    mpfr_mul(b->zri, b->zr, b->zi, MPFR_RNDN);
    mpfr_sub(b->t, b->zr2, b->zi2, MPFR_RNDN);
    mpfr_add(b->zr, b->t, b->cr, MPFR_RNDN);
    mpfr_mul_2ui(b->zri, b->zri, 1, MPFR_RNDN);
    mpfr_add(b->zi, b->zri, b->ci, MPFR_RNDN);
    if (barrier->wait(stop)) break;
    if (*haveperiod) break;
  }
  SetEvent(b->hDone);
  mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
  return 0;
}

static int m_d_box_period_do(const complex<flyttyp> &center, flyttyp radius, int maxperiod,int &steps,HWND hWnd) {
	 radius = flyttyp(4)/radius;

  mp_bitcnt_t bits = mpfr_get_prec(center.m_r.m_dec.backend().data());
  barrier bar(4);
  complex<flyttyp> c[4];
  c[0] = center + complex<flyttyp>(-radius, -radius);
  c[1] = center + complex<flyttyp>(radius, -radius);
  c[2] = center + complex<flyttyp>(radius, radius);
  c[3] = center + complex<flyttyp>(-radius, radius);
  int crossing[4] = { 0, 0, 0, 0 };
  int haveperiod = false;
  int period = 0;
  HANDLE hDone[4];
  // prepare threads
  BoxPeriod box[4];
  for (int t = 0; t < 4; ++t)
  {
    box[t].threadid = t;
    box[t].barrier = &bar;
    box[t].maxperiod = maxperiod;
    mpfr_init2(box[t].cr, bits); mpfr_set(box[t].cr, c[t].m_r.m_dec.backend().data(), MPFR_RNDN);
    mpfr_init2(box[t].ci, bits); mpfr_set(box[t].ci, c[t].m_i.m_dec.backend().data(), MPFR_RNDN);
    mpfr_init2(box[t].zr, bits); mpfr_set(box[t].zr, c[t].m_r.m_dec.backend().data(), MPFR_RNDN);
    mpfr_init2(box[t].zi, bits); mpfr_set(box[t].zi, c[t].m_i.m_dec.backend().data(), MPFR_RNDN);
    mpfr_init2(box[t].zr2, bits);
    mpfr_init2(box[t].zi2, bits);
    mpfr_init2(box[t].zri, bits);
    mpfr_init2(box[t].t, bits);
    box[t].z2r = &box[(t + 1) % 4].zr;
    box[t].z2i = &box[(t + 1) % 4].zi;
    box[t].stop = &g_bNewtonStop;
    box[t].crossing[0] = &crossing[0];
    box[t].crossing[1] = &crossing[1];
    box[t].crossing[2] = &crossing[2];
    box[t].crossing[3] = &crossing[3];
    box[t].haveperiod = &haveperiod;
    box[t].period = &period;
    box[t].hDone = hDone[t] = CreateEvent(NULL, 0, 0, NULL);
    box[t].hWnd = hWnd;
  }
  // spawn threads
  for (int i = 0; i < 4; i++)
  {
    DWORD dw;
    HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)ThBoxPeriod, (LPVOID)&box[i], 0, &dw);
    CloseHandle(hThread);
  }
  // wait for threads to complete
  WaitForMultipleObjects(4, hDone, TRUE, INFINITE);
  for (int i = 0; i < 4; i++)
  {
    CloseHandle(hDone[i]);
    mpfr_clear(box[i].cr);
    mpfr_clear(box[i].ci);
    mpfr_clear(box[i].zr);
    mpfr_clear(box[i].zi);
    mpfr_clear(box[i].zr2);
    mpfr_clear(box[i].zi2);
    mpfr_clear(box[i].zri);
    mpfr_clear(box[i].t);
  }
  steps = period;
  return haveperiod ? period : 0;

}
#endif

struct STEP_STRUCT_COMMON
{
	HWND hWnd;
	barrier *barrier;
	volatile bool *stop;
	int newtonStep;
	int64_t period;
	mpfr_t zr, zi, zrn, zin, zr2, zi2, cr, ci, dcr, dci, dcrn, dcin, dcrzr, dcrzi, dcizr, dcizi;
};
struct STEP_STRUCT
{
	int nType;
	HANDLE hDone;
	STEP_STRUCT_COMMON *common;
};
static DWORD WINAPI ThStep(STEP_STRUCT *t0)
{
  SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
  struct STEP_STRUCT_COMMON *t = t0->common;
  char szStatus[256];
  uint32_t last = GetTickCount();
  switch (t0->nType)
  {
    case 0: // zr
      for (int64_t i = 0; i < t->period; ++i)
      {
	mpfr_sqr(t->zr2, t->zr, MPFR_RNDN);
	mpfr_sqr(t->zi2, t->zi, MPFR_RNDN);
	mpfr_sub(t->zrn, t->zr2, t->zi2, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
	mpfr_add(t->zr, t->zrn, t->cr, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
      }
      break;
    case 1: // zi
      for (int64_t i = 0; i < t->period; ++i)
      {
	if(i%100==0){
		uint32_t now = GetTickCount();
		if (now - last > 250)
		{
			wsprintf(szStatus,"NR %d/%d (%d%%) %s",t->newtonStep + 1,g_nNewtonETA ? t->newtonStep + g_nNewtonETA : 0,i*100/t->period,g_szProgress);
			SetDlgItemText(t->hWnd,IDC_EDIT1,szStatus);
			last = now;
		}
	}
	mpfr_mul(t->zin, t->zr, t->zi, MPFR_RNDN);
	mpfr_mul_2ui(t->zin, t->zin, 1, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
	mpfr_add(t->zi, t->zin, t->ci, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
      }
      break;
    case 2: // dcr
      for (int64_t i = 0; i < t->period; ++i)
      {
	mpfr_mul(t->dcrzr, t->dcr, t->zr, MPFR_RNDN);
	mpfr_mul(t->dcizi, t->dci, t->zi, MPFR_RNDN);
	mpfr_sub(t->dcrn, t->dcrzr, t->dcizi, MPFR_RNDN);
	mpfr_mul_2ui(t->dcrn, t->dcrn, 1, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
	mpfr_add_ui(t->dcr, t->dcrn, 1, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
      }
      break;
    case 3: // dci
      for (int64_t i = 0; i < t->period; ++i)
      {
	mpfr_mul(t->dcrzi, t->dcr, t->zi, MPFR_RNDN);
	mpfr_mul(t->dcizr, t->dci, t->zr, MPFR_RNDN);
	mpfr_add(t->dcin, t->dcrzi, t->dcizr, MPFR_RNDN);
	mpfr_mul_2ui(t->dcin, t->dcin, 1, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
	mpfr_set(t->dci, t->dcin, MPFR_RNDN);
	if (t->barrier->wait(t->stop)) break;
      }
      break;
  }
  SetEvent(t0->hDone);
  mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
  return 0;
}

static int m_d_nucleus_step(complex<flyttyp> *c_out, const complex<flyttyp> &c_guess, const int64_t period,const flyttyp &epsilon2,HWND hWnd,int newtonStep, const flyttyp &radius2) {
  complex<flyttyp> z(0,0);
  complex<flyttyp> zr(0,0);
  complex<flyttyp> dc(0,0);
  int i;

	int threads = 4;
	mp_bitcnt_t bits = mpfr_get_prec(c_guess.m_r.m_dec.backend().data());
	barrier bar(4);
	STEP_STRUCT_COMMON m;
	m.hWnd = hWnd;
	m.barrier = &bar;
	m.stop = &g_bNewtonStop;
	m.newtonStep = newtonStep;
	m.period = period;
	mpfr_init2(m.zr, bits); mpfr_set_ui(m.zr, 0, MPFR_RNDN);
	mpfr_init2(m.zi, bits); mpfr_set_ui(m.zi, 0, MPFR_RNDN);
	mpfr_init2(m.zrn, bits);
	mpfr_init2(m.zin, bits);
	mpfr_init2(m.zr2, bits);
	mpfr_init2(m.zi2, bits);
	mpfr_init2(m.cr, bits); mpfr_set(m.cr, c_guess.m_r.m_dec.backend().data(), MPFR_RNDN);
	mpfr_init2(m.ci, bits); mpfr_set(m.ci, c_guess.m_i.m_dec.backend().data(), MPFR_RNDN);
	mpfr_init2(m.dcr, bits); mpfr_set_ui(m.dcr, 0, MPFR_RNDN);
	mpfr_init2(m.dci, bits); mpfr_set_ui(m.dci, 0, MPFR_RNDN);
	mpfr_init2(m.dcrn, bits);
	mpfr_init2(m.dcin, bits);
	mpfr_init2(m.dcrzr, bits);
	mpfr_init2(m.dcrzi, bits);
	mpfr_init2(m.dcizr, bits);
	mpfr_init2(m.dcizi, bits);
	STEP_STRUCT mc[4];
	HANDLE hDone[4];
	for (i = 0; i<4; i++){
		mc[i].nType =i;
		hDone[i] = mc[i].hDone = CreateEvent(NULL, 0, 0, NULL);
		mc[i].common = &m;
	}

	HANDLE hThread;
	DWORD dw;
	for (i = 0; i<threads; i++){
		hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)ThStep, (LPVOID)&mc[i], 0, &dw);
		CloseHandle(hThread);
	}

	WaitForMultipleObjects(threads, hDone, TRUE, INFINITE);
	for (i = 0; i<threads; i++){
		CloseHandle(hDone[i]);
	}

	mpfr_set(z.m_r.m_dec.backend().data(), m.zr, MPFR_RNDN);
	mpfr_set(z.m_i.m_dec.backend().data(), m.zi, MPFR_RNDN);
	mpfr_set(dc.m_r.m_dec.backend().data(), m.dcr, MPFR_RNDN);
	mpfr_set(dc.m_i.m_dec.backend().data(), m.dci, MPFR_RNDN);
	mpfr_clear(m.zr);
	mpfr_clear(m.zi);
	mpfr_clear(m.zrn);
	mpfr_clear(m.zin);
	mpfr_clear(m.zr2);
	mpfr_clear(m.zi2);
	mpfr_clear(m.cr);
	mpfr_clear(m.ci);
	mpfr_clear(m.dcr);
	mpfr_clear(m.dci);
	mpfr_clear(m.dcrn);
	mpfr_clear(m.dcin);
	mpfr_clear(m.dcrzr);
	mpfr_clear(m.dcrzi);
	mpfr_clear(m.dcizr);
	mpfr_clear(m.dcizi);

  SetDlgItemText(hWnd,IDC_EDIT4,"");
  flyttyp ad;
#if 0
  = 1/cabs2(dc);
  if (ad < epsilon2) {
    *c_out = c_guess;
    return 0;
  }
#endif
  if(dc.m_r==0 && dc.m_i==0)
	  return -1;
  complex<flyttyp> c_new = c_guess - z / dc;
  complex<flyttyp> d = c_new - c_guess;
  ad = cabs2(d);

  // progress reporting
  floatexp delta = floatexp(CFixedFloat(ad.m_dec));
  floatexp epsilon = floatexp(CFixedFloat(epsilon2.m_dec));
  std::string elast = sqrt(floatexp(4.0)/delta).toString(1);
  std::string etarget = sqrt(floatexp(4.0)/epsilon).toString(1);
  snprintf(g_szProgress, sizeof(g_szProgress) - 1, "%s %s %s", elast.c_str(), ad < epsilon2 ? ">" : "<", etarget.c_str());
  g_fNewtonDelta2[0] = g_fNewtonDelta2[1];
  g_fNewtonDelta2[1] = delta;
  if (g_fNewtonDelta2[0] > 0)
    g_nNewtonETA = newtonETA(g_fNewtonDelta2[0], g_fNewtonDelta2[1], epsilon);

  if (ad < epsilon2) {
    *c_out = c_new;
    return 0;
  }
  if (cisfinite(d) && cabs2(c_new) < 16 && ad < radius2) {
    *c_out = c_new;
    return 1;
  } else {
    *c_out = c_guess;
    return -1;
  }
}

bool SaveNewtonBackup(const std::string &szFile, const std::string &re, const std::string &im, const std::string &zoom, int64_t period)
{
	if (! g_SFT.GetSaveNewtonProgress())
		return true;
	bool overwrite = true; // FIXME
	CStringTable stSave;
	stSave.AddRow();
	stSave.AddString(stSave.GetCount() - 1, "Re");
	stSave.AddString(stSave.GetCount() - 1, re);
	stSave.AddRow();
	stSave.AddString(stSave.GetCount() - 1, "Im");
	stSave.AddString(stSave.GetCount() - 1, im);
	stSave.AddRow();
	stSave.AddString(stSave.GetCount() - 1, "Zoom");
	stSave.AddString(stSave.GetCount() - 1, zoom);
	stSave.AddRow();
	stSave.AddString(stSave.GetCount() - 1, "Period");
	stSave.AddInt   (stSave.GetCount() - 1, period);
	char *szData = stSave.ToText(": ", "\r\n");
	HANDLE hFile = CreateFile(szFile.c_str(), GENERIC_WRITE, 0, NULL, overwrite ? CREATE_ALWAYS : CREATE_NEW, 0, NULL);
	if (hFile == INVALID_HANDLE_VALUE)
	{
		stSave.DeleteToText(szData);
		return false;
	}
	DWORD dw;
	WriteFile(hFile, szData, strlen(szData), &dw, NULL);
	CloseHandle(hFile);
	stSave.DeleteToText(szData);
	return true;
}

bool SaveNewtonBackup(const complex<flyttyp> &c_new, const complex<flyttyp> &c_old, int64_t period, int step)
{
  complex<flyttyp> delta = c_new - c_old;
  complex<floatexp> delta_lo = complex<floatexp>(delta);
  std::string re = c_new.m_r.ToText();
  std::string im = c_new.m_i.ToText();
  std::string zoom = sqrt(floatexp(4.0) / cabs2(delta_lo)).toString();
  char extension[100];
  snprintf(extension, sizeof(extension) - 1, "newton-%04d.kfr", step);
  return SaveNewtonBackup(g_szFile == "" ? extension : replace_path_extension(g_szFile, extension), re, im, zoom, period);
}

static int m_d_nucleus(complex<flyttyp> *c_out, const complex<flyttyp> &c_guess, int64_t period, int maxsteps,int &steps,const flyttyp &radius,HWND hWnd) {
  int result = -1, i;
  complex<flyttyp> c = c_guess;
  complex<flyttyp> c_new;

  flyttyp epsilon2 = flyttyp(1)/(radius*radius*radius);
  flyttyp radius2 = radius * radius;
  for (i = 0; i < maxsteps && !g_bNewtonStop && !g_bNewtonExit; ++i) {
    result = m_d_nucleus_step(&c_new, c, period,epsilon2,hWnd,i,radius2);
    SaveNewtonBackup(c_new, c, period, i + 1);
    c = c_new;
    if (result != 1)
      break;
  }
  steps = i;
  *c_out = c;
  if(g_bNewtonExit)
	  result=0;
  return result;
}

static inline complex<floatexp> fec(const complex<flyttyp> &z)
{
  return complex<floatexp>(mpfr_get_fe(z.m_r.m_dec.backend().data()), mpfr_get_fe(z.m_i.m_dec.backend().data()));
}

static complex<floatexp> m_d_size(const complex<flyttyp> &nucleus, int64_t period,HWND hWnd)
{
  complex<floatexp> fec1(1,0);
  complex<floatexp> fec2(2,0);
  complex<floatexp> l(1,0);
  complex<floatexp> b(1,0);
  complex<flyttyp> z(0,0);
  char szStatus[256];
  uint32_t last = GetTickCount();
  for (int64_t i = 1; i < period && !g_bNewtonStop; ++i) {
	  if(i%100==0){
		uint32_t now = GetTickCount();
		if (now - last > 250)
		{
		  wsprintf(szStatus,"Determine size %d%%...",100*i/period);
		  SetDlgItemText(hWnd,IDC_EDIT1,szStatus);
		  last = now;
		}
	  }
    z = z * z + nucleus;
    complex<floatexp> zlo = fec(z);
    l = fec2 * zlo * l;
    b = b + fec1 / l;
  }
  return fec1 / (b * l * l);
}

static int g_useDZ = 1;
static double g_skew[4];
int64_t g_period = 0;

static int WINAPI ThSkew(HWND hWnd)
{
  SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
  const int type = g_SFT.GetFractalType();
  const int power = g_SFT.GetPower();
  const struct formula *f = get_formula(type, power);

  const int64_t iters = g_SFT.GetIterations();
  g_szRe = g_SFT.GetRe();
  g_szIm = g_SFT.GetIm();
  g_szZoom = g_SFT.GetZoom();
  flyttyp radius = g_szZoom;
  radius*=g_SFT.GetZoomSize();
  const char *e = strstr(g_szZoom.c_str(),"E");
  if(!e) e = strstr(g_szZoom.c_str(),"e");
  int exp = (e?atoi(e+1):0);
  unsigned uprec = exp + 6;
  Precision prec(uprec);
  complex<flyttyp> center(g_szRe,g_szIm);

  g_skew[0] = 1;
  g_skew[1] = 0;
  g_skew[2] = 0;
  g_skew[3] = 1;
  // fork progress updater
  progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
  HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThSkewProgress, (LPVOID) &progress, 0, NULL);
  CloseHandle(hThread);
  int ret = 0;
  if (g_SFT.GetUseHybridFormula())
  {
    if (hybrid_skew(g_SFT.GetHybridFormula(), iters, center.m_r.m_dec, center.m_i.m_dec, g_useDZ, &g_skew[0], &running, &progress.counters[0]))
    {
      ret = 2;
    }
    else
    {
      ret = -2;
    }
  }
  else
  if (f && f->skew && f->skew(iters, g_FactorAR, g_FactorAI, center.m_r.m_dec.backend().data(), center.m_i.m_dec.backend().data(), g_useDZ, &g_skew[0], &running, &progress.counters[0]))
  {
    ret = 2;
  }
  else
  {
    ret = -2;
  }
  // join progress updater
  progress.running = false;
  WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
  CloseHandle(progress.hDone);
  PostMessage(hWnd,WM_USER+2,0,ret);
  mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
  return 0;
}


static int WINAPI ThNewton(HWND hWnd)
{
  SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
  const int type = g_SFT.GetFractalType();
  const int power = g_SFT.GetPower();
  const struct formula *f = get_formula(type, power);
  g_skew[0] = 1;
  g_skew[1] = 0;
  g_skew[2] = 0;
  g_skew[3] = 1;

	char szStatus[300];

	flyttyp radius = g_szZoom;
	radius*=g_SFT.GetZoomSize();
	const char *e = strstr(g_szZoom.c_str(),"E");
	if(!e)
		e = strstr(g_szZoom.c_str(),"e");
	int exp = (e?atoi(e+1):0);
	unsigned uprec = exp + 6;
	Precision prec(uprec);
	complex<flyttyp> center(g_szRe,g_szIm);

	char szVal[25];
	int i;
	for(i=0;i<24 && g_szZoom[i] && g_szZoom[i]!='e' && g_szZoom[i]!='E' && g_szZoom[i]!='+' && g_szZoom[i]!='-';i++)
		szVal[i]=g_szZoom[i];
	szVal[i]=0;
	e = strstr(g_szZoom.c_str(),"E");
	if(!e)
		e = strstr(g_szZoom.c_str(),"e");
	int startZooms = (e?atof(e+1)/0.30103:0) + log10(atof(szVal));

	int steps = 0;
	if(SendDlgItemMessage(hWnd,IDC_CHECK1,BM_GETCHECK,0,0)){
		g_period = GetDlgItemInt(hWnd,IDC_EDIT3,NULL,0);
		uprec *= 3;
	}
	else{
	  if (g_SFT.GetUseHybridFormula())
	  {
		  flyttyp r = flyttyp(4) / radius;
		  // fork progress updater
		  progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
		  HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThPeriodProgress, (LPVOID) &progress, 0, NULL);
		  CloseHandle(hThread);
		  g_period = hybrid_period(g_SFT.GetHybridFormula(), INT_MAX, center.m_r, center.m_i, r, &g_skew[0], &running, &progress.counters[0]);
		  // join progress updater
		  progress.running = false;
		  WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
		  CloseHandle(progress.hDone);
		  if (g_period < 0) g_period = 0;
	  }
	  else if (type == 0 && power == 2)
	  {
		int64_t maxperiod = INT_MAX; // FIXME
		g_period= ball_period_do(center,radius,maxperiod,steps,hWnd);
	  }
	  else 
	  {
		if (f)
		{
		  flyttyp r = flyttyp(4) / radius;
		  // fork progress updater
		  progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
		  HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThPeriodProgress, (LPVOID) &progress, 0, NULL);
		  CloseHandle(hThread);
		  if (g_bUseBallPeriod)
		  {
		    g_period = f->period_jsk(INT_MAX, 1e50, g_FactorAR, g_FactorAI, center.m_r.m_dec.backend().data(), center.m_i.m_dec.backend().data(), r.m_dec.backend().data(), &g_skew[0], &running, &progress.counters[0]);
		  }
		  else
		  {
		    g_period = f->period_tri(INT_MAX, 1e50, g_FactorAR, g_FactorAI, center.m_r.m_dec.backend().data(), center.m_i.m_dec.backend().data(), r.m_dec.backend().data(), &running, &progress.counters[0]);
		  }
		  // join progress updater
		  progress.running = false;
		  WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
		  CloseHandle(progress.hDone);
		  if (g_period < 0) g_period = 0;
		}
		else
		{
		  g_period = 0;
		}
	  }
	  uprec *= 2;
	}
	// std::cerr << "period " << g_period << std::endl;
	Precision prec2(uprec);

	SetDlgItemInt(hWnd,IDC_EDIT3,g_period,0);
	BOOL bOK=-1;
	if(g_period){
		sprintf(szStatus,"period=%lld (steps:%d)\n",g_period,steps);
		SetDlgItemText(hWnd,IDC_EDIT1,szStatus);
		complex<flyttyp> c;
		int test = 1;
		if (g_SFT.GetUseHybridFormula())
		{
		    // fork progress updater
		    progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
		    HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThNewtonProgress, (LPVOID) &progress, 0, NULL);
		    CloseHandle(hThread);
		    // c = center; // seems to copy precision too, so do it low-level...
		    mpfr_set(c.m_r.m_dec.backend().data(), center.m_r.m_dec.backend().data(), MPFR_RNDN);
		    mpfr_set(c.m_i.m_dec.backend().data(), center.m_i.m_dec.backend().data(), MPFR_RNDN);
		    flyttyp epsilon2 = flyttyp(1)/(radius*radius*radius);
		    test = hybrid_newton(g_SFT.GetHybridFormula(), 100, g_period, c.m_r, c.m_i, epsilon2, &running, &progress.counters[0]) ? 0 : 1;
		    flyttyp r = flyttyp(4) / radius;
		    if (! (cabs2(c - center) < r * r))
		      test = 1;
		    steps = 1;
		    // join progress updater
		    progress.running = false;
		    WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
		    CloseHandle(progress.hDone);
		}
		else if (type == 0 && power == 2)
		{
		  int maxsteps = INT_MAX; // FIXME
		  test = m_d_nucleus(&c,center,g_period,maxsteps,steps,radius,hWnd);
		}
		else
		{
		  if (f)
		  {
		    // fork progress updater
		    progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
		    HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThNewtonProgress, (LPVOID) &progress, 0, NULL);
		    CloseHandle(hThread);
		    // c = center; // seems to copy precision too, so do it low-level...
		    mpfr_set(c.m_r.m_dec.backend().data(), center.m_r.m_dec.backend().data(), MPFR_RNDN);
		    mpfr_set(c.m_i.m_dec.backend().data(), center.m_i.m_dec.backend().data(), MPFR_RNDN);
		    flyttyp epsilon2 = flyttyp(1)/(radius*radius*radius);
		    test = f->newton(100, g_period, g_FactorAR, g_FactorAI, c.m_r.m_dec.backend().data(), c.m_i.m_dec.backend().data(), epsilon2.m_dec.backend().data(), &running, &progress.counters[0]) ? 0 : 1;
		    flyttyp r = flyttyp(4) / radius;
		    if (! (cabs2(c - center) < r * r))
		      test = 1;
		    steps = 1;
		    // join progress updater
		    progress.running = false;
		    WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
		    CloseHandle(progress.hDone);
		  }
		}
		// std::cerr << "newton " << test << " " << steps << std::endl;

		if(test==0 && steps){
			g_szRe = c.m_r.ToText();
			g_szIm = c.m_i.ToText();

			Precision prec3(exp + 6);
			flyttyp msize = 0;
			if (g_SFT.GetUseHybridFormula())
			{
			    // fork progress updater
			    progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
			    HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThSizeProgress, (LPVOID) &progress, 0, NULL);
			    CloseHandle(hThread);
			    hybrid_size(g_SFT.GetHybridFormula(), g_period, c.m_r, c.m_i, msize, &g_skew[0], &running, &progress.counters[0]);
			    // join progress updater
			    progress.running = false;
			    WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
			    CloseHandle(progress.hDone);
			    msize = flyttyp(.25) / msize;
			}
			else if (type == 0 && power == 2)
			{
			  complex<floatexp> size = m_d_size(c,g_period,hWnd);
			  floatexp msizefe = floatexp(.25)/sqrt(cabs2(size));
			  mpfr_set_fe(msize.m_dec.backend().data(), msizefe);
			}
			else
			{
			  if (f)
			  {
			    // fork progress updater
			    progress_t progress = { { 0, 0, 0, 0 }, true, hWnd, CreateEvent(NULL, 0, 0, NULL) };
			    HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) ThSizeProgress, (LPVOID) &progress, 0, NULL);
			    CloseHandle(hThread);
			    f->size(g_period, g_FactorAR, g_FactorAI, c.m_r.m_dec.backend().data(), c.m_i.m_dec.backend().data(), msize.m_dec.backend().data(), &g_skew[0], &running, &progress.counters[0]);
			    // join progress updater
			    progress.running = false;
			    WaitForMultipleObjects(1, &progress.hDone, TRUE, INFINITE);
			    CloseHandle(progress.hDone);
			    msize = flyttyp(.25) / msize;
			  }
			}
			// std::cerr << "size " << msize.m_dec << std::endl;

			double zooms1;
			{
			  std::ostringstream oss;
			  oss << std::scientific << msize.m_dec;
			  std::string sszSize = oss.str();
			  char *szSize0 = strdup(sszSize.c_str());
			  char *szSize = szSize0;
			  char szTmpSize[200];
			  if(!strstr(szSize,"e") && !strstr(szSize,"E")){
				  char *szP = strstr(szSize, ".");
				  if (szP)
					  *szP = 0;
				  int exp = strlen(szSize) - 1;
				  if (exp>2){
					  int end = 12;
					  if (end>exp)
						  end = exp;
					  szSize[end + 1] = 0;
					  while (end>1){
						  szSize[end] = szSize[end - 1];
						  end--;
					  }
					  szSize[end] = '.';
					  strcpy(szTmpSize,szSize);
					  char szNum[20];
					  sprintf(szNum, "E%d", exp);
					  strcat(szTmpSize, szNum);
					  szSize = szTmpSize;
				  }
			  }
			  if(strstr(szSize,"e") || strstr(szSize,"E")){
				  int i;
				  for(i=0;i<24 && szSize[i] && szSize[i]!='e' && szSize[i]!='E' && szSize[i]!='+' && szSize[i]!='-';i++)
					  szVal[i]=szSize[i];
				  szVal[i]=0;
				  e = strstr(szSize,"E");
				  if(!e)
					  e = strstr(szSize,"e");
#define LOG_10_2 0.30102999566398114
				  zooms1 = (e?atof(e+1)/LOG_10_2:0) + log10(atof(szVal));
			  }
			  else
				  zooms1 = log10(atof(szSize))/LOG_10_2;
#undef LOG_10_2
			  g_szZoom = szSize;
			  free(szSize0);
			}

			double zooms0;
			{
			  char *szSize0 = strdup(g_sMinibrotSourceZoom.c_str());
			  char *szSize = szSize0;
			  char szTmpSize[200];
			  if(!strstr(szSize,"e") && !strstr(szSize,"E")){
				  char *szP = strstr(szSize, ".");
				  if (szP)
					  *szP = 0;
				  int exp = strlen(szSize) - 1;
				  if (exp>2){
					  int end = 12;
					  if (end>exp)
						  end = exp;
					  szSize[end + 1] = 0;
					  while (end>1){
						  szSize[end] = szSize[end - 1];
						  end--;
					  }
					  szSize[end] = '.';
					  strcpy(szTmpSize,szSize);
					  char szNum[20];
					  sprintf(szNum, "E%d", exp);
					  strcat(szTmpSize, szNum);
					  szSize = szTmpSize;
				  }
			  }
			  if(strstr(szSize,"e") || strstr(szSize,"E")){
				  int i;
				  for(i=0;i<24 && szSize[i] && szSize[i]!='e' && szSize[i]!='E' && szSize[i]!='+' && szSize[i]!='-';i++)
					  szVal[i]=szSize[i];
				  szVal[i]=0;
				  e = strstr(szSize,"E");
				  if(!e)
					  e = strstr(szSize,"e");
#define LOG_10_2 0.30102999566398114
				  zooms0 = (e?atof(e+1)/LOG_10_2:0) + log10(atof(szVal));
			  }
			  else
				  zooms0 = log10(atof(szSize))/LOG_10_2;
#undef LOG_10_2
			  free(szSize0);
			}

			double zooms = zooms1;
			double start_iterations = g_SFT.GetIterations();
			double minibrot_iterations = 50.0 * g_period;

			if(g_nMinibrotPos==0)
			{
				g_iterations = minibrot_iterations;
				zooms = zooms1;
			}
			else if(g_nMinibrotPos==1)
			{
				g_iterations = 2 * start_iterations;
				zooms = zooms0 * 0.5 + 0.5 * zooms1;
			}
			else if(g_nMinibrotPos==2)
			{
				g_iterations = 3 * start_iterations;
				zooms = zooms0 * 0.25 + 0.75 * zooms1;
			}
			else if(g_nMinibrotPos==3)
			{
				g_iterations = 4 * start_iterations;
				zooms = zooms0 * 0.125 + 0.875 * zooms1;
			}
			else if(g_nMinibrotPos==4)
			{
				double f = 1 - (1 - g_nMinibrotFactor) * 2;
				g_iterations = fmin(fmax(start_iterations * pow(2, -log(1 - f)), start_iterations), minibrot_iterations);
				zooms = zooms0 * (1 - f) + f * zooms1;
			}
			if(4 * g_period > zooms && zooms>startZooms)
			{
				radius = flyttyp(2)^zooms;
				g_szZoom = radius.ToText();
				bOK=1;
			}
		}
	}
	g_bNewtonRunning=FALSE;
	PostMessage(hWnd,WM_USER+2,0,bOK);
	mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
	return 0;
}
static __int64 t1, t2;
static SYSTEMTIME st1, st2;
extern int WINAPI NewtonProc(HWND hWnd,UINT uMsg,WPARAM wParam,LPARAM lParam)
{
	static int g_AutoSkew = 0;
	if(uMsg==WM_INITDIALOG){
		SendMessage(hWnd, WM_SETICON, ICON_SMALL, LPARAM(g_hIcon));
		SendMessage(hWnd, WM_SETICON, ICON_BIG, LPARAM(g_hIcon));
//		InitToolTip(hWnd,GetModuleHandle(NULL),GetToolText,1);
		if(g_nMinibrotPos==1)
			SendDlgItemMessage(hWnd,IDC_RADIO2,BM_SETCHECK,1,0);
		else if(g_nMinibrotPos==2)
			SendDlgItemMessage(hWnd,IDC_RADIO3,BM_SETCHECK,1,0);
		else if(g_nMinibrotPos==3)
			SendDlgItemMessage(hWnd,IDC_RADIO4,BM_SETCHECK,1,0);
		else if(g_nMinibrotPos==4)
			SendDlgItemMessage(hWnd,IDC_RADIO5,BM_SETCHECK,1,0);
		else
			SendDlgItemMessage(hWnd,IDC_RADIO1,BM_SETCHECK,1,0);
		if (g_AutoSkew)
			SendDlgItemMessage(hWnd,IDC_AUTOSKEW,BM_SETCHECK,1,0);
		if (g_bUseBallPeriod)
			SendDlgItemMessage(hWnd,IDC_BALL_PERIOD,BM_SETCHECK,1,0);
		SendDlgItemMessage(hWnd, IDC_NEWTON_PROGRESS, BM_SETCHECK, g_SFT.GetSaveNewtonProgress() ? 1 : 0,0);
		std::string z = g_SFT.GetZoom();
		SetDlgItemText(hWnd, IDC_EDIT4, z.c_str());
		std::ostringstream s;
		s << g_nMinibrotFactor;
		SetDlgItemText(hWnd, IDC_EDIT2, s.str().c_str());
		return 1;
	}
	if(uMsg==WM_COMMAND && wParam==IDCANCEL){
		if(g_bNewtonRunning){
			running = 0;
			g_bNewtonStop=TRUE;
			while(g_bNewtonRunning){
				MSG msg;
				while(PeekMessage(&msg,NULL,0,0,PM_REMOVE)){
					TranslateMessage(&msg);
					DispatchMessage(&msg);
				}
				Sleep(1);
			}
			SetDlgItemText(hWnd,IDCANCEL,"Close");
		}
		else{
			HWND hMain = GetParent(hWnd);
			PostMessage(hMain,WM_COMMAND,ID_SPECIAL_NEWTON,0);
			DestroyWindow(hWnd);
		}
	}
	if(uMsg==WM_COMMAND && wParam==IDCANCEL2){
		if(g_bNewtonRunning){
			g_bNewtonExit=TRUE;
			EnableWindow(GetDlgItem(hWnd,IDCANCEL2),FALSE);
		}
		else
			MessageBeep((UINT)-1);
	}
	if(uMsg==WM_COMMAND && wParam==IDC_BUTTON2){
		std::string z = g_SFT.GetZoom();
		SetDlgItemText(hWnd, IDC_EDIT4, z.c_str());
	}
	if(uMsg==WM_COMMAND && wParam==IDC_AUTOSKEW_BUTTON){
		if(!g_bNewtonRunning){
			DWORD dw;
			g_bNewtonStop=FALSE;
			g_bNewtonExit=FALSE;
			*g_szProgress=0;
			g_nNewtonETA = 0;
			g_fNewtonDelta2[0] = 0;
			g_fNewtonDelta2[1] = 0;
			running = 1;
			g_useDZ = SendDlgItemMessage(hWnd,IDC_AUTOSKEW_USEDZ,BM_GETCHECK,0,0);
			HANDLE hThread = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)ThSkew,hWnd,0,&dw);
			CloseHandle(hThread);
			g_bNewtonRunning=TRUE;
			SetDlgItemText(hWnd,IDCANCEL,"Stop");
		}
	}
	if(uMsg==WM_USER+1){
		if(!g_bNewtonRunning){
			mat2 m = g_SFT.GetTransformMatrix();
			g_skew[0] = m[0][0];
			g_skew[1] = m[0][1];
			g_skew[2] = m[1][0];
			g_skew[3] = m[1][1];
			RECT r = *(RECT*)lParam;
			g_szRe = g_SFT.GetRe(r.left,r.top,r.right,r.bottom);
			g_szIm = g_SFT.GetIm(r.left,r.top,r.right,r.bottom);
			g_szZoom = g_SFT.GetZoom();
			g_nMinibrotPos=0;
			if(SendDlgItemMessage(hWnd,IDC_RADIO2,BM_GETCHECK,0,0))
				g_nMinibrotPos=1;
			else if(SendDlgItemMessage(hWnd,IDC_RADIO3,BM_GETCHECK,0,0))
				g_nMinibrotPos=2;
			else if(SendDlgItemMessage(hWnd,IDC_RADIO4,BM_GETCHECK,0,0))
				g_nMinibrotPos=3;
			else if(SendDlgItemMessage(hWnd,IDC_RADIO5,BM_GETCHECK,0,0))
				g_nMinibrotPos=4;
			g_bUseBallPeriod = SendDlgItemMessage(hWnd,IDC_BALL_PERIOD,BM_GETCHECK,0,0);
			g_SFT.SetSaveNewtonProgress(SendDlgItemMessage(hWnd, IDC_NEWTON_PROGRESS, BM_GETCHECK, 0, 0));
			DWORD dw;
			g_bNewtonStop=FALSE;
			g_bNewtonExit=FALSE;
			char szText[256];
			GetDlgItemText(hWnd,IDC_EDIT2,szText,sizeof(szText));
			g_nMinibrotFactor = atof(szText);
			GetDlgItemText(hWnd,IDC_EDIT4,szText,sizeof(szText));
			g_sMinibrotSourceZoom = std::string(szText);
			*g_szProgress=0;
			g_nNewtonETA = 0;
			g_fNewtonDelta2[0] = 0;
			g_fNewtonDelta2[1] = 0;
			GetLocalTime(&st1);
			running = 1;
			HANDLE hThread = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)ThNewton,hWnd,0,&dw);
			CloseHandle(hThread);
			g_bNewtonRunning=TRUE;
			SetDlgItemText(hWnd,IDCANCEL,"Stop");
		}
	}
	if(uMsg==WM_USER+2){
		SetDlgItemText(hWnd,IDCANCEL,"Close");
		g_bNewtonRunning=FALSE;
		if(lParam == 1){
			// newton success
			g_SFT.UndoStore();
			if((g_AutoSkew = SendDlgItemMessage(hWnd,IDC_AUTOSKEW,BM_GETCHECK,0,0)))
				g_SFT.SetTransformMatrix(mat2(g_skew[0], g_skew[1], g_skew[2], g_skew[3]));
			g_SFT.SetPosition(g_szRe,g_szIm,g_szZoom);
			g_SFT.SetIterations(g_iterations);
			g_bJustDidNewton = true;
			PostMessage(GetParent(hWnd),WM_KEYDOWN,VK_F5,0);
			if(SendDlgItemMessage(hWnd,IDC_CHECK9,BM_GETCHECK,0,0)){
				GetLocalTime(&st2);
				SystemTimeToFileTime(&st1,(LPFILETIME)&t1);
				SystemTimeToFileTime(&st2,(LPFILETIME)&t2);
				__int64 nT2 = t2-t1;
				FileTimeToSystemTime((LPFILETIME)&nT2,&st1);
				char szResult[256];
				wsprintf(szResult,"Time: %02d:%02d:%02d.%03d.%05d\nPeriod: %d",st1.wHour,st1.wMinute,st1.wSecond,st1.wMilliseconds,(int)(nT2%10000),g_period);
				MessageBox(hWnd,szResult,"Result",MB_OK);
			}
		}
		if(lParam == 2){
			// auto skew success
			g_SFT.UndoStore();
			g_SFT.SetTransformMatrix(mat2(g_skew[0], g_skew[1], g_skew[2], g_skew[3]));
			PostMessage(GetParent(hWnd),WM_KEYDOWN,VK_F5,0);
		}
		PostMessage(GetParent(hWnd),WM_COMMAND,ID_SPECIAL_NEWTON,0);
		if(lParam == -1 && !g_bNewtonStop)
			MessageBox(GetParent(hWnd),"Could not apply Newton-Raphson\nYou may zoom in a little and try again","Error",MB_OK|MB_ICONSTOP);
		if(lParam == -2 && !g_bNewtonStop)
			MessageBox(GetParent(hWnd),"Could not auto find skew\n","Error",MB_OK|MB_ICONSTOP);
		if((lParam == 0 || lParam < -2 || lParam > 2))
			MessageBox(GetParent(hWnd),"Unexpected stop message parameter (internal error)\n","Error",MB_OK|MB_ICONSTOP);
	}
	return 0;
}
