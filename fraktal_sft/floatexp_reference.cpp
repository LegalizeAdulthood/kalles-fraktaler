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

#include "fraktal_sft.h"
#include <float.h>
#include "complex.h"

#include "../formula/formula.h"

extern double g_real;
extern double g_imag;
DWORD WINAPI ThMC2(MC2 *pMC);
DWORD WINAPI ThMC(MC *pMC);
BOOL ISFLOATOK(double a);
extern double g_SeedR;
extern double g_SeedI;
extern double g_FactorAR;
extern double g_FactorAI;

#include "../common/barrier.h"

struct mcthread_common
{
	barrier *barrier;
	mpfr_t xr, xi, xrn, xin, xrn1, xin1, xrxid, xrxid1, sr, si, cr, ci;
	floatexp *m_dxr, *m_dxi;
	double *m_db_z, *terminate, *glitch_threshold;
	int64_t *m_nMaxIter, *m_nGlitchIter, *nMaxIter, *m_nRDone;
	volatile bool *stop;
};

struct mcthread
{
	int nType;
	HANDLE hDone;
	mcthread_common *common;
};

static DWORD WINAPI mcthreadfunc(mcthread *p0)
{
	SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
	mcthread_common *p = p0->common;
	double glitch_threshold = *p->glitch_threshold;
	int64_t i = 0;
	switch (p0->nType)
	{
		case 0:
		{
			for (i = 0; i < *p->nMaxIter; i++)
			{
				mpfr_sub(p->xrn1, p->sr, p->si, MPFR_RNDN);
				mpfr_add(p->xrn, p->cr, p->xrn1, MPFR_RNDN);
				if (p->barrier->wait(p->stop)) break;
				mpfr_set(p->xr, p->xrn, MPFR_RNDN);
				mpfr_sqr(p->sr, p->xrn, MPFR_RNDN);
				p->m_dxr[i] = mpfr_get_fe(p->xrn);
				if (p->barrier->wait(p->stop)) break;
			}
		}
		break;
		case 1:
		{
			for (i = 0; i < *p->nMaxIter; i++)
			{
				mpfr_add(p->xin, p->sr, p->si, MPFR_RNDN);
				mpfr_sub(p->xin1, p->ci, p->xin, MPFR_RNDN);
				mpfr_add(p->xin, p->xin1, p->xrxid, MPFR_RNDN);
				if (p->barrier->wait(p->stop)) break;
				mpfr_set(p->xi, p->xin, MPFR_RNDN);
				mpfr_sqr(p->si, p->xin, MPFR_RNDN);
				p->m_dxi[i] = mpfr_get_fe(p->xin);
				if (p->barrier->wait(p->stop)) break;
			}
		}
		break;
		case 2:
		{
			for (i = 0; i < *p->nMaxIter; i++)
			{
				if (i > 0)
				{
					const floatexp lr = p->m_dxr[i-1];
					const floatexp li = p->m_dxi[i-1];
					const double abs_val = (lr * lr + li * li).todouble();
					p->m_db_z[i-1] = abs_val * glitch_threshold;
					if (abs_val >= *p->terminate){
						if (*p->nMaxIter == *p->m_nMaxIter)
						{
							*p->nMaxIter = i-1 + 3;
							if (*p->nMaxIter > *p->m_nMaxIter)
								*p->nMaxIter = *p->m_nMaxIter;
							*p->m_nGlitchIter = *p->nMaxIter;
						}
					}
					(*p->m_nRDone)++;
				}
				if (p->barrier->wait(p->stop)) break;
				mpfr_add(p->xrxid1, p->xrn, p->xin, MPFR_RNDN);
				mpfr_sqr(p->xrxid, p->xrxid1, MPFR_RNDN);
				if (p->barrier->wait(p->stop)) break;
			}
		}
		break;
	}
	if (p->barrier->wait(p->stop))
	{
		SetEvent(p0->hDone);
		return 0;
	}
	if (p0->nType == 2)
	{
		if (i > 0)
		{
			const floatexp lr = p->m_dxr[i-1];
			const floatexp li = p->m_dxi[i-1];
			const double abs_val = (lr * lr + li * li).todouble();
			p->m_db_z[i-1] = abs_val * glitch_threshold;
			if (abs_val >= *p->terminate){
				if (*p->nMaxIter == *p->m_nMaxIter)
				{
					*p->nMaxIter = i-1 + 3;
					if (*p->nMaxIter > *p->m_nMaxIter)
						*p->nMaxIter = *p->m_nMaxIter;
					*p->m_nGlitchIter = *p->nMaxIter;
				}
			}
			(*p->m_nRDone)++;
		}
		floatexp xr; xr = mpfr_get_fe(p->xr);
		floatexp xi; xi = mpfr_get_fe(p->xi);
		for (; i < *p->nMaxIter && !*p->stop; i++)
		{
			p->m_dxr[i] = xr;
			p->m_dxi[i] = xi;
		}
	}
	SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_END);
	SetEvent(p0->hDone);
	mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
	return 0;
}


void CFraktalSFT::CalculateReferenceEXP()
{
	Precision prec(m_rref.m_f.precision());

	int64_t i;
	if (m_dxr)
		delete[] m_dxr;
	m_dxr = new floatexp[m_nMaxIter];
	if (m_dxi)
		delete[] m_dxi;
	m_dxi = new floatexp[m_nMaxIter];
	if (m_db_z)
		delete[] m_db_z;
	m_db_z = new double [m_nMaxIter];

	double terminate = SMOOTH_BAILOUT*SMOOTH_BAILOUT;
	m_nGlitchIter = m_nMaxIter + 1;
	int64_t nMaxIter = m_nMaxIter;
	if (GetUseHybridFormula())
	{
		bool ok = reference(GetHybridFormula(), m_dxr, m_dxi, m_db_z, m_bStop, m_nRDone, m_nGlitchIter, m_nMaxIter, m_rref, m_iref, g_SeedR, g_SeedI, terminate, GetGlitchLowTolerance());
		assert(ok && "reference_long_double_hybrid");
	}
	else if (m_nFractalType == 0 && m_nPower == 2 && GetThreadedReference()) // FIXME matrix derivatives, option to disable derivatives
	{
		double glitch_threshold = 0.0000001;
		if (GetGlitchLowTolerance()) {
			glitch_threshold = sqrt(glitch_threshold);
		}

		mcthread mc[3];
		barrier barrier(3);
		HANDLE hDone[3];
		mcthread_common co;
	  co.barrier = &barrier;
		mp_bitcnt_t bits = mpfr_get_prec(m_rref.m_f.backend().data());
		mpfr_init2(co.xr, bits);
		mpfr_init2(co.xi, bits);
		mpfr_init2(co.xrn, bits);
		mpfr_init2(co.xin, bits);
		mpfr_init2(co.xrn1, bits);
		mpfr_init2(co.xin1, bits);
		mpfr_init2(co.xrxid, bits);
		mpfr_init2(co.xrxid1, bits);
		mpfr_init2(co.sr, bits);
		mpfr_init2(co.si, bits);
		mpfr_init2(co.cr, bits);
		mpfr_init2(co.ci, bits);
		mpfr_set(co.cr, m_rref.m_f.backend().data(), MPFR_RNDN);
		mpfr_set(co.ci, m_iref.m_f.backend().data(), MPFR_RNDN);
		mpfr_set_d(co.xr, g_SeedR, MPFR_RNDN);
		mpfr_set_d(co.xi, g_SeedI, MPFR_RNDN);
		mpfr_sqr(co.sr, co.xr, MPFR_RNDN);
		mpfr_sqr(co.si, co.xi, MPFR_RNDN);
		mpfr_set_d(co.xrxid, 0, MPFR_RNDN);
		co.m_dxr = m_dxr;
		co.m_dxi = m_dxi;
		co.m_db_z = m_db_z;
		co.terminate = &terminate;
		co.glitch_threshold = &glitch_threshold;
		co.m_nMaxIter = &m_nMaxIter;
		co.m_nGlitchIter = &m_nGlitchIter;
		co.nMaxIter = &nMaxIter;
		co.m_nRDone = &m_nRDone;
		co.stop = &m_bStop;
		// spawn threads
		for (i = 0; i < 3; i++)
		{
			mc[i].nType = i;
			hDone[i] = mc[i].hDone = CreateEvent(NULL, 0, 0, NULL);
			mc[i].common = &co;
			HANDLE hThread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) mcthreadfunc, (LPVOID)&mc[i], 0, NULL);
			CloseHandle(hThread);
		}
		// wait for completion
		WaitForMultipleObjects(3, hDone, TRUE, INFINITE);
		for (i = 0; i < 3; i++){
			CloseHandle(hDone[i]);
		}
		mpfr_clear(co.xr);
		mpfr_clear(co.xi);
		mpfr_clear(co.xrn);
		mpfr_clear(co.xin);
		mpfr_clear(co.xrn1);
		mpfr_clear(co.xin1);
		mpfr_clear(co.xrxid);
		mpfr_clear(co.xrxid1);
		mpfr_clear(co.sr);
		mpfr_clear(co.si);
		mpfr_clear(co.cr);
		mpfr_clear(co.ci);
	}
	else if (m_nFractalType == 0 && m_nPower > 10) // FIXME matrix derivatives, option to disable derivatives
	{
		CFixedFloat xr = g_SeedR, xi = g_SeedI;
		double threashold = 0.0001;
		for (i = 7; i <= m_nPower; i += 2)
			threashold *= 10;
		if (GetGlitchLowTolerance()) {
			threashold = sqrt(threashold);
		}
		if (threashold>.5)
			threashold = .5;
		for (i = 0; i<nMaxIter && !m_bStop; i++){
			complex<CFixedFloat> X(xr, xi), r(m_rref, m_iref);
			complex<CFixedFloat> Xn = (X^m_nPower) + r;
			xr = Xn.m_r;
			xi = Xn.m_i;
			m_dxr[i] = xr;
			m_dxi[i] = xi;
			const double abs_val = (m_dxr[i] * m_dxr[i] + m_dxi[i] * m_dxi[i]).todouble();
			m_db_z[i] = abs_val*threashold;
			if (abs_val >= terminate){
				if (nMaxIter == m_nMaxIter){
					nMaxIter = i + 3;
					if (nMaxIter>m_nMaxIter)
						nMaxIter = m_nMaxIter;
					m_nGlitchIter = nMaxIter;
				}
			}
			m_nRDone++;
		}

	}
	else
	{
    bool ok = reference(m_nFractalType, m_nPower, m_dxr, m_dxi, m_db_z, m_bStop, m_nRDone, m_nGlitchIter, m_nMaxIter, m_rref, m_iref, g_SeedR, g_SeedI, g_FactorAR, g_FactorAI, terminate, GetGlitchLowTolerance());
    assert(ok && "reference_floatexp");
	}
}
