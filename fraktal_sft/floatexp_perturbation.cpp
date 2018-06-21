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

#include "../formula/formula.h"

BOOL ISFLOATOK(double a);
extern double g_real;
extern double g_imag;
extern double g_FactorAR;
extern double g_FactorAI;
#define _abs(a) ((_abs_val=(a))>0?_abs_val:-_abs_val)
floatexp lb_abs_exp(const floatexp &c, const floatexp &d)
{
	floatexp abs_val, _abs_val, _2=2;
	if (c>0){
		if (c + d>0)
			abs_val = d;
		else if (d == -c)
			abs_val = d;
		else if (d<-c)
			abs_val = -d - _2 * c;
	}
	else if (c == 0)
		abs_val = _abs(d);
	else if (c < 0){
		if (c + d>0)
			abs_val = d + _2 * c;
		else if (d == -c)
			abs_val = -d;
		else if (d < -c)
			abs_val = -d;
	}
	return abs_val;
}

void CFraktalSFT::MandelCalcEXP()
{
	m_bIterChanged = TRUE;
	floatexp Dnr, Dni, yr, yi;
	int x, y, w, h;
	floatexp real(g_real), imag(g_imag), _abs_val;

  int k = 0;
  bool partial = false;
  int64_t_v16 xs, ys, ws, hs;
	while (!m_bStop && !partial)
  {
	  if (m_P.GetPixel(x, y, w, h, m_bMirrored))
	  {
			int nIndex = x * 3 + (m_bmi->biHeight - 1 - y)*m_row;
			if (m_nPixels[x][y] != -1){
				SetColor(nIndex, m_nPixels[x][y], m_nTrans[x][y], x, y);
				if (m_bMirrored)
					Mirror(x, y);
				continue;
			}
	    if (GuessPixel(x, y, w, h))
	      continue;
      xs[k] = x;
      ys[k] = y;
      ws[k] = w;
      hs[k] = h;
      ++k;
		}
		else
		{
			partial = true;
		}
		if (k == 16 || partial)
		{

		// Series approximation
		floatexp_v16 D0r(0);
		floatexp_v16 D0i(0);
		floatexp_v16 daa(1);
		floatexp_v16 dab(0);
		floatexp_v16 dba(0);
		floatexp_v16 dbb(1);
		for (int j = 0; j < k; ++j)
		{
			floatexp D0rj(0);
			floatexp D0ij(0);
			floatexp daaj(1);
			floatexp dabj(0);
			floatexp dbaj(0);
			floatexp dbbj(1);
			GetPixelCoordinates(xs[j], ys[j], D0rj, D0ij, daaj, dabj, dbaj, dbbj);
      D0r.val[j] = D0rj.val; D0r.exp[j] = D0rj.exp;
      D0i.val[j] = D0ij.val; D0i.exp[j] = D0ij.exp;
      daa.val[j] = daaj.val; daa.exp[j] = daaj.exp;
      dab.val[j] = dabj.val; dab.exp[j] = dabj.exp;
      dba.val[j] = dbaj.val; dba.exp[j] = dbaj.exp;
      dbb.val[j] = dbbj.val; dbb.exp[j] = dbbj.exp;
		}
		int64_t_v16 antal = {};
		floatexp_v16 Dr;
		floatexp_v16 Di;
		floatexp_v16 dr;
		floatexp_v16 di;
		for (int j = 0; j < k; ++j)
		{
			int antalj;
			floatexp Drj;
			floatexp Dij;
			floatexp drj;
			floatexp dij;
			DoApproximation(antalj, D0r[j], D0i[j], Drj, Dij, drj, dij);
			antal[j] = antalj;
      Dr.val[j] = Drj.val; Dr.exp[j] = Drj.exp;
      Di.val[j] = Dij.val; Di.exp[j] = Dij.exp;
      dr.val[j] = drj.val; dr.exp[j] = drj.exp;
      di.val[j] = dij.val; di.exp[j] = dij.exp;
		}
		for (int j = k; j < 16; ++j)
		{
			antal[j] = 0;
			Dr.val[j] = 0.0; Dr.exp[j] = EXP_MIN;
			Di.val[j] = 0.0; Di.exp[j] = EXP_MIN;
			dr.val[j] = 0.0; dr.exp[j] = EXP_MIN;
			di.val[j] = 0.0; di.exp[j] = EXP_MIN;
		}

		// initialized to 0
		double_v16 test1 = {}, test2 = {};
		int64_t_v16 bGlitch = {};
		int nMaxIter = (m_nGlitchIter<m_nMaxIter ? m_nGlitchIter : m_nMaxIter);

		if (m_nFractalType == 0 && m_nPower == 2)
		{
			bool ok = perturbation_0_2_v16(m_nFractalType, m_nPower, m_dxr, m_dxi, m_db_z, antal, test1, test2, bGlitch, m_nBailout2, nMaxIter, Dr, Di, D0r, D0i, k);
			assert(ok && "perturbation_floatexp_v16()");
		}
		else
		{
			for (int j = 0; j < k; ++j)
			{

    if (m_nFractalType == 0 && m_nPower > 10) // FIXME matrix derivatives
		{
			if (GetDerivatives())
			{
			complex<floatexp> d(dr[j], di[j]);
			if (antal[j]<nMaxIter && test1[j] <= m_nBailout2){
				for (; antal[j]<nMaxIter; antal[j]++){
					yr = m_dxr[antal[j]] + Dr[j];
					yi = m_dxi[antal[j]] + Di[j];
					test2[j] = test1[j];
					test1[j] = (real*yr*yr + imag*yi*yi).todouble();
					if (test1[j]<m_db_z[antal[j]]){
						bGlitch[j] = TRUE;
						if (! m_bNoGlitchDetection)
							break;
					}
					if (test1[j] > m_nBailout2)
					{
						break;
					}
					complex<floatexp> y(yr, yi);
					d = m_nPower * d * (y ^ (m_nPower - 1)) + 1;
					complex<floatexp> X(m_dxr[antal[j]], m_dxi[antal[j]]);
					complex<floatexp> D(Dr[j], Di[j]);
					complex<floatexp> D0(D0r[j], D0i[j]);
					complex<floatexp> c(m_pnExpConsts[0], 0);
					int nXExp = m_nPower - 2, nDExp = 2, ci = 1;
					complex<floatexp> Dn = c*(X^(m_nPower - 1))*D;
					while (nXExp){
						c.m_r = m_pnExpConsts[ci++];
						Dn += c*(X^nXExp)*(D^nDExp);
						nXExp--;
						nDExp++;
					}
					Dn += (D^m_nPower) + D0;
					Dr.val[j] = Dn.m_r.val; Dr.exp[j] = Dn.m_r.exp;
					Di.val[j] = Dn.m_i.val; Di.exp[j] = Dn.m_i.exp;
				}
			}
			} else {
			if (antal[j]<nMaxIter && test1[j] <= m_nBailout2){
				for (; antal[j]<nMaxIter; antal[j]++){
					yr = m_dxr[antal[j]] + Dr[j];
					yi = m_dxi[antal[j]] + Di[j];
					test2[j] = test1[j];
					test1[j] = (real*yr*yr + imag*yi*yi).todouble();
					if (test1[j]<m_db_z[antal[j]]){
						bGlitch[j] = TRUE;
						if (! m_bNoGlitchDetection)
							break;
					}
					if (test1[j] > m_nBailout2)
					{
						break;
					}
					complex<floatexp> y(yr, yi);
					complex<floatexp> X(m_dxr[antal[j]], m_dxi[antal[j]]);
					complex<floatexp> D(Dr[j], Di[j]);
					complex<floatexp> D0(D0r[j], D0i[j]);
					complex<floatexp> c(m_pnExpConsts[0], 0);
					int nXExp = m_nPower - 2, nDExp = 2, ci = 1;
					complex<floatexp> Dn = c*(X^(m_nPower - 1))*D;
					while (nXExp){
						c.m_r = m_pnExpConsts[ci++];
						Dn += c*(X^nXExp)*(D^nDExp);
						nXExp--;
						nDExp++;
					}
					Dn += (D^m_nPower) + D0;
					Dr.val[j] = Dn.m_r.val; Dr.exp[j] = Dn.m_r.exp;
					Di.val[j] = Dn.m_i.val; Di.exp[j] = Dn.m_i.exp;
				}
			}
			}
		}
    else
    {
			int antalj = antal[j];
			double test1j = test1[j];
			double test2j = test2[j];
			int bGlitchj = bGlitch[j];
			floatexp drj = dr[j];
			floatexp dij = di[j];
			floatexp Drj = Dr[j];
			floatexp Dij = Di[j];
			floatexp D0rj = D0r[j];
			floatexp D0ij = D0i[j];
			floatexp daaj = daa[j];
			floatexp dabj = dab[j];
			floatexp dbaj = dba[j];
			floatexp dbbj = dbb[j];
			drj *= m_fPixelSpacing;
			dij *= m_fPixelSpacing;
			bool ok = GetDerivatives()
			  ? perturbation_floatexp(m_nFractalType, m_nPower, m_dxr, m_dxi, m_db_z, antalj, test1j, test2j, bGlitchj, m_nBailout2, nMaxIter, m_bNoGlitchDetection, g_real, g_imag, g_FactorAR, g_FactorAI, Drj, Dij, D0rj, D0ij, drj, dij, m_epsilon, m_fPixelSpacing, daaj, dabj, dbaj, dbbj)
			  : perturbation_floatexp(m_nFractalType, m_nPower, m_dxr, m_dxi, m_db_z, antalj, test1j, test2j, bGlitchj, m_nBailout2, nMaxIter, m_bNoGlitchDetection, g_real, g_imag, g_FactorAR, g_FactorAI, Drj, Dij, D0rj, D0ij)
			  ;
			assert(ok && "perturbation_floatexp()");
			antal[j] = antalj;
			test1[j] = test1j;
			test2[j] = test2j;
			bGlitch[j] = bGlitchj;
			dr.val[j] = drj.val; dr.exp[j] = drj.exp;
			di.val[j] = dij.val; di.exp[j] = dij.exp;
			Dr.val[j] = Drj.val; Dr.exp[j] = Drj.exp;
			Di.val[j] = Dij.val; Di.exp[j] = Dij.exp;
			D0r.val[j] = D0rj.val; D0r.exp[j] = D0rj.exp;
			D0i.val[j] = D0ij.val; D0i.exp[j] = D0ij.exp;
			daa.val[j] = daaj.val; daa.exp[j] = daaj.exp;
			dab.val[j] = dabj.val; dab.exp[j] = dabj.exp;
			dba.val[j] = dbaj.val; dba.exp[j] = dbaj.exp;
			dbb.val[j] = dbbj.val; dbb.exp[j] = dbbj.exp;
		}
			}
		}

		for (int j = 0; j < k; ++j)
		{

		double de = GetDerivatives()
		  ? double(sqrt(test1[j]) * log(test1[j]) / sqrt(dr[j] * dr[j] + di[j] * di[j]).todouble())
		  : 0
		  ;

		OutputIterationData(xs[j], ys[j], bGlitch[j], antal[j], test1[j], test2[j], de);

		InterlockedIncrement((LPLONG)&m_nDone);
    OutputPixelData(xs[j], ys[j], ws[j], hs[j], bGlitch[j]);
    
		}

			k = 0;
		}

	}
}
