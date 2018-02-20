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
	int antal, x, y, w, h;
	floatexp real(g_real), imag(g_imag), _abs_val;


	while (!m_bStop && m_P.GetPixel(x, y, w, h, m_bMirrored)){
		int nIndex = x * 3 + (m_bmi->biHeight - 1 - y)*m_row;
		if (m_nPixels[x][y] != -1){
			SetColor(nIndex, m_nPixels[x][y], m_nTrans[x][y], x, y);
			if (m_bMirrored)
				Mirror(x, y);
			continue;
		}
    if (GuessPixel(x, y, w, h))
      continue;
		// Series approximation
		floatexp dbD0r = 0;
		floatexp dbD0i = 0;
		GetPixelCoordinates(x, y, dbD0r, dbD0i);

		floatexp D0r = dbD0r;//(cr-rref);
		floatexp D0i = dbD0i;
		floatexp Dr = D0r;
		floatexp Di = D0i;
		if (m_nMaxApproximation){
			antal = m_nMaxApproximation - 1;
			Dnr = m_APr[0] * D0r - m_APi[0] * D0i;
			Dni = m_APr[0] * D0i + m_APi[0] * D0r;
			floatexp D_r = D0r*D0r - D0i*D0i;
			floatexp D_i = (D0r*D0i).mul2();
			Dnr += m_APr[1] * D_r - m_APi[1] * D_i;
			Dni += m_APr[1] * D_i + m_APi[1] * D_r;
			int k;
			int m_nTerms = GetApproxTerms();
			for (k = 2; k<m_nTerms; k++){
				floatexp  t = D_r*D0r - D_i*D0i;
				D_i = D_r*D0i + D_i*D0r;
				D_r = t;
				Dnr += m_APr[k] * D_r - m_APi[k] * D_i;
				Dni += m_APr[k] * D_i + m_APi[k] * D_r;
			}
			Dr = Dnr;
			Di = Dni;
		}
		else{
			antal = 0;
			Dr = D0r;
			Di = D0i;
		}


		double test1 = 0, test2 = 0;
		BOOL bGlitch = FALSE;
		int nMaxIter = (m_nGlitchIter<m_nMaxIter ? m_nGlitchIter : m_nMaxIter);

    if (m_nFractalType == 0 && m_nPower > 10)
		{
			if (antal<nMaxIter && test1 <= m_nBailout2){
				for (; antal<nMaxIter && test1 <= m_nBailout2; antal++){
					yr = m_dxr[antal] + Dr;
					yi = m_dxi[antal] + Di;
					test2 = test1;
					test1 = (real*yr*yr + imag*yi*yi).todouble();
					if (test1<m_db_z[antal]){
						bGlitch = TRUE;
						if (!m_bNoGlitchDetection)
							break;
					}
					complex<floatexp> X(m_dxr[antal], m_dxi[antal]);
					complex<floatexp> D(Dr, Di);
					complex<floatexp> D0(D0r, D0i);
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
					Di = Dn.m_i;
					Dr = Dn.m_r;
				}
			}

		}
    else
    {

			bool ok = perturbation_floatexp(m_nFractalType, m_nPower, m_dxr, m_dxi, m_db_z, antal, test1, test2, bGlitch, m_nBailout2, nMaxIter, m_bNoGlitchDetection, g_real, g_imag, g_FactorAR, g_FactorAI, Dr, Di, D0r, D0i);
			assert(ok && "perturbation_floatexp()");

		}

		OutputIterationData(x, y, bGlitch, antal, test1, test2);

		InterlockedIncrement((LPLONG)&m_nDone);
    OutputPixelData(x, y, w, h, bGlitch);
	}
}
