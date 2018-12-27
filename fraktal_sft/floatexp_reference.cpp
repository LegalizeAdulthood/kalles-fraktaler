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

#include "fraktal_sft.h"
#include <float.h>
#include "complex.h"

#include "../formula/generated/formula.h"

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
	floatexp *m_dz, *terminate, *glitch_threshold;
	int *m_nMaxIter, *m_nGlitchIter, *nMaxIter, *m_nRDone;
	int *antal;
	floatexp *test1;
	floatexp *test2;
	volatile BOOL *stop;
	floatexp dr, di;
};

struct mcthread
{
	int nType;
	HANDLE hDone;
	mcthread_common *common;
};

static DWORD WINAPI mcthreadfunc(mcthread *p0)
{
	bool stored = false;
	floatexp old_absval = 0;
	floatexp abs_val = 0;
	int antal = 0;
	floatexp test1 = 0;
	floatexp test2 = 0;
	const floatexp real(g_real);
	const floatexp imag(g_imag);
	mcthread_common *p = p0->common;
	floatexp dr = p->dr;
	floatexp di = p->di;
	floatexp glitch_threshold = *p->glitch_threshold;
	int i = 0;
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
				p->m_dxr[i] = mpfr_get_fe(p->xrn, MPFR_RNDN);
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
				p->m_dxi[i] = mpfr_get_fe(p->xin, MPFR_RNDN);
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
					floatexp drn = 2 * (lr * dr - li * di) + 1;
					floatexp din = 2 * (lr * di + li * dr);
					dr = drn;
					di = din;
					old_absval = abs_val;
					abs_val = real * lr * lr + imag * li * li;
					p->m_dz[i-1] = abs_val * glitch_threshold;
					if (abs_val >= 4)
					{
						if (*p->terminate == 4 && !stored)
						{
							stored = true;
							antal = i;
							test1 = abs_val;
							test2 = old_absval;
						}
					}
					if (abs_val >= *p->terminate){
						if (*p->terminate > 4 && !stored)
						{
							stored = true;
							antal = i;
							test1 = abs_val;
							test2 = old_absval;
						}
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
			floatexp drn = 2 * (lr * dr - li * di) + 1;
			floatexp din = 2 * (lr * di + li * dr);
			dr = drn;
			di = din;
			old_absval = abs_val;
			abs_val = real * lr * lr + imag * li * li;
			p->m_dz[i-1] = abs_val * glitch_threshold;
			if (abs_val >= 4)
			{
				if (*p->terminate == 4 && !stored)
				{
					stored = true;
					antal = i;
					test1 = abs_val;
					test2 = old_absval;
				}
			}
			if (abs_val >= *p->terminate){
				if (*p->terminate > 4 && !stored)
				{
					stored = true;
					antal = i;
					test1 = abs_val;
					test2 = old_absval;
				}
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
		floatexp xr; xr = mpfr_get_fe(p->xr, MPFR_RNDN);
		floatexp xi; xi = mpfr_get_fe(p->xi, MPFR_RNDN);
		for (; i < *p->nMaxIter && !*p->stop; i++)
		{
			p->m_dxr[i] = xr;
			p->m_dxi[i] = xi;
		}
		*p->antal = antal;
		*p->test1 = test1;
		*p->test2 = test2;
		p->dr = dr;
		p->di = di;
	}
	SetEvent(p0->hDone);
	mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
	return 0;
}


void CFraktalSFT::CalculateReferenceEXP()
{
	Precision prec(m_rref.m_f.precision());

	int i;

	DeleteReferenceOrbit();
	m_dxr = new floatexp [m_nMaxIter];
	m_dxi = new floatexp [m_nMaxIter];
	m_dz  = new floatexp [m_nMaxIter];

	floatexp real(g_real);
	floatexp imag(g_imag);

	int antal = 0;
	floatexp test1 = 0;
	floatexp test2 = 0;

	floatexp dr = 1, di = 0;

	floatexp terminate = SMOOTH_BAILOUT*SMOOTH_BAILOUT;
	m_nGlitchIter = m_nMaxIter + 1;
	int nMaxIter = m_nMaxIter;
	if (m_nFractalType == 0 && m_nPower == 2) // FIXME matrix derivatives, option to disable derivatives
	{
		floatexp glitch_threshold = 0.0000001;
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
		co.m_dz = m_dz;
		co.terminate = &terminate;
		co.glitch_threshold = &glitch_threshold;
		co.m_nMaxIter = &m_nMaxIter;
		co.m_nGlitchIter = &m_nGlitchIter;
		co.nMaxIter = &nMaxIter;
		co.m_nRDone = &m_nRDone;
		co.antal = &antal;
		co.test1 = &test1;
		co.test2 = &test2;
		co.stop = &m_bStop;
		co.dr = dr;
		co.di = di;
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
    dr = co.dr;
    di = co.di;
		floatexp pixel_spacing = m_fPixelSpacing;
		dr = dr * pixel_spacing;
		di = di * pixel_spacing;
	}
	else if (m_nFractalType == 0 && m_nPower > 10) // FIXME matrix derivatives, option to disable derivatives
	{
		bool stored = false;
		floatexp old_absval = 0;
		floatexp abs_val = 0;
		CFixedFloat xr = g_SeedR, xi = g_SeedI;
		floatexp threashold = 0.0001;
		for (i = 7; i <= m_nPower; i += 2)
			threashold *= 10.0;
		if (GetGlitchLowTolerance()) {
			threashold = sqrt(threashold);
		}
		if (threashold>.5)
			threashold = .5;
		complex<floatexp> d(1.0, 0.0);
		for (i = 0; i<nMaxIter && !m_bStop; i++){
			complex<CFixedFloat> X(xr, xi), r(m_rref, m_iref);
			complex<CFixedFloat> Xn = (X^m_nPower) + r;
			floatexp xrf(xr);
			floatexp xif(xi);
			complex<floatexp> x(xrf, xif);
			d = m_nPower * d * (x ^ (m_nPower - 1)) + 1;
			xr = Xn.m_r;
			xi = Xn.m_i;
			m_dxr[i] = floatexp(xr);
			m_dxi[i] = floatexp(xi);
			old_absval = abs_val;
			abs_val = real * m_dxr[i] * m_dxr[i] + imag * m_dxi[i] * m_dxi[i];
			m_dz[i] = abs_val*threashold;
			if (abs_val >= 4)
			{
				if (terminate == 4 && !stored)
				{
					stored = true;
					antal = i;
					test1 = abs_val;
					test2 = old_absval;
				}
			}
			if (abs_val >= terminate){
				if (terminate > 4 && !stored)
				{
					stored = true;
					antal = i;
					test1 = abs_val;
					test2 = old_absval;
				}
				if (nMaxIter == m_nMaxIter){
					nMaxIter = i + 3;
					if (nMaxIter>m_nMaxIter)
						nMaxIter = m_nMaxIter;
					m_nGlitchIter = nMaxIter;
				}
			}
			m_nRDone++;
		}
		dr = d.m_r;
		di = d.m_i;
		floatexp pixel_spacing = m_fPixelSpacing;
		dr = dr * pixel_spacing;
		di = di * pixel_spacing;

	}
	else
	{

		floatexp _x, _y, daa, dab, dba, dbb;
		GetPixelCoordinates(g_nAddRefX, g_nAddRefY, _x, _y, daa, dab, dba, dbb);
		floatexp test1f, test2f;
		bool ok;
		if (GetDerivatives())
		{
			floatexp dzc[2] = { 0, 0 };
			floatexp dci[4] = { daa, dab, dba, dbb };
			ok = current_formula->referenceDfe(m_nFractalType, m_nPower, m_dxr, m_dxi, m_dz, &m_bStop, &m_nRDone, &m_nGlitchIter, &m_nMaxIter, m_rref.m_f.backend().data(), m_iref.m_f.backend().data(), floatexp(g_SeedR), floatexp(g_SeedI), floatexp(g_FactorAR), floatexp(g_FactorAI), floatexp(terminate), floatexp(g_real), floatexp(g_imag), GetGlitchLowTolerance(), &antal, &test1f, &test2f, &dzc[0], &dci[0]);
			dr = dzc[0] * m_fPixelSpacing;
			di = dzc[1] * m_fPixelSpacing;
		}
		else
		{
			ok = current_formula->referencefe(m_nFractalType, m_nPower, m_dxr, m_dxi, m_dz, &m_bStop, &m_nRDone, &m_nGlitchIter, &m_nMaxIter, m_rref.m_f.backend().data(), m_iref.m_f.backend().data(), floatexp(g_SeedR), floatexp(g_SeedI), floatexp(g_FactorAR), floatexp(g_FactorAI), floatexp(terminate), floatexp(g_real), floatexp(g_imag), GetGlitchLowTolerance(), &antal, &test1f, &test2f);
	  }
	  test1 = test1f;
	  test2 = test2f;
    assert(ok && "reference_floatexp");

	}

	double de = GetDerivatives()
	  ? double(sqrt(test1) * log(test1) / sqrt(dr * dr + di * di))
	  : 0
	  ;

	if (0 <= g_nAddRefX && g_nAddRefX < m_nX && 0 <= g_nAddRefY && g_nAddRefY < m_nY)
		OutputIterationData(g_nAddRefX, g_nAddRefY, 1, 1, false, antal ? antal + 1 : m_nMaxIter, double(test1), double(test2), de);

}
