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
#include "complex.h"
#include "../formula/formula.h"

void CFraktalSFT::CalculateReference()
{
	Precision prec(m_rref.m_f.precision());

	int i;
	if (m_db_dxr)
		delete[] m_db_dxr;
	m_db_dxr = new double [m_nMaxIter];
	if (m_db_dxi)
		delete[] m_db_dxi;
	m_db_dxi = new double [m_nMaxIter];
	if (m_db_z)
		delete[] m_db_z;
	m_db_z = new double [m_nMaxIter];

	CFixedFloat xr = g_SeedR, xi = g_SeedI, xin, xrn, sr = xr.Square(), si = xi.Square(), xrxid = 0;
	double terminate = SMOOTH_BAILOUT*SMOOTH_BAILOUT;

	if (m_Inflections.size() > 0)
	{
		int inf;
		complex<CFixedFloat> c(m_rref,m_iref);
		for(inf=m_Inflections.size()-1;inf>=0;inf--){
			complex<CFixedFloat> d = c-m_Inflections[inf];
			c=m_Inflections[inf]+d*d;
		}
		m_rref=c.m_r;
		m_iref=c.m_i;
	}

	m_nGlitchIter = m_nMaxIter + 1;
	int nMaxIter = m_nMaxIter;

	int antal = 0;
	double test1 = 0;
	double test2 = 0;

  if (m_nFractalType == 0 && m_nPower > 10)
	{

		double threashold = 0.0001;
		for (i = 7; i <= m_nPower; i += 2)
			threashold *= 10;
		if (GetGlitchLowTolerance()) {
			threashold = sqrt(threashold);
		}
		if (threashold>.5)
			threashold = .5;
		complex<CFixedFloat> r(m_rref, m_iref);
		complex<CFixedFloat> X(g_SeedR, g_SeedI);
		bool stored = false;
		double old_absval = 0;
		double abs_val = 0;
		for (i = 0; i<nMaxIter && !m_bStop; i++){
			X = (X^m_nPower) + r;
			m_db_dxr[i] = X.m_r.ToDouble();
			m_db_dxi[i] = X.m_i.ToDouble();
			old_absval = abs_val;
			abs_val = (g_real * m_db_dxr[i] * m_db_dxr[i] + g_imag * m_db_dxi[i] * m_db_dxi[i]);
			m_db_z[i] = abs_val*threashold;
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

	}
	else
	{

		bool ok = reference_double(m_nFractalType, m_nPower, m_db_dxr, m_db_dxi, m_db_z, m_bStop, m_nRDone, m_nGlitchIter, m_nMaxIter, m_rref, m_iref, g_SeedR, g_SeedI, g_FactorAR, g_FactorAI, terminate, g_real, g_imag, GetGlitchLowTolerance(), antal, test1, test2);
		assert(ok && "reference_double");

	}

	if (0 <= g_nAddRefX && g_nAddRefX < m_nX && 0 <= g_nAddRefY && g_nAddRefY < m_nY)
		OutputIterationData(g_nAddRefX, g_nAddRefY, false, antal ? antal + 1 : m_nMaxIter, test1, test2);
}

