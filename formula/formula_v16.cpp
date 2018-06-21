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

#include "formula.h"

#define FORMULA2(function,type,power) function ## _ ## type ## _ ## power ## _v16
#define FORMULA(a,b,c) FORMULA2(a,b,c)

bool FORMULA(perturbation,0,2)
  ( int m_nFractalType, int m_nPower
  , const floatexp *m_db_dxr, const floatexp *m_db_dxi, const double *m_db_z
  , int64_t_v16 &antal0, double_v16 &test10, double_v16 &test20, int64_t_v16 &bGlitch
  , double m_nBailout2, const int nMaxIter
  , floatexp_v16 &xr0, floatexp_v16 &xi0
  , const floatexp_v16 &cr, const floatexp_v16 &ci
  , int k
  )
{
  if (m_nFractalType == 0 && m_nPower == 2)
  {
    int64_t_v16 antal = antal0;
    double_v16 test1 = test10;
    double_v16 test2 = test20;
    floatexp_v16 xr = xr0;
    floatexp_v16 xi = xi0;
    for (; all(antal < nMaxIter); antal = antal + 1)
    {
      // FIXME assumes antal[i] == antal[0] for all i
      const floatexp_v16 Xr(m_db_dxr[antal[0]]);
      const floatexp_v16 Xi(m_db_dxi[antal[0]]);
      const double       Xz(m_db_z  [antal[0]]);
      const floatexp_v16 Xxr = Xr + xr;
      const floatexp_v16 Xxi = Xi + xi;
      test2 = test1;
      test1 = double_v16(Xxr * Xxr + Xxi * Xxi);
      bGlitch |= test1 < Xz;
      if (any(bGlitch)) break;
      if (any(test1 > m_nBailout2)) break;
      floatexp_v16 xrn, xin;
      {
        xrn=(((((Xr.mul2())+xr)*xr)-(((Xi.mul2())+xi)*xi))+cr);
        xin=((((Xxr*xi)+(Xi*xr)).mul2())+ci);
      }
      xr = xrn;
      xi = xin;
    }
    for (int i = 0; i < k; ++i)
    {
      for (; antal[i] < nMaxIter; antal[i] = antal[i] + 1)
      {
        const floatexp Xr = m_db_dxr[antal[i]];
        const floatexp Xi = m_db_dxi[antal[i]];
        const double   Xz = m_db_z[antal[i]];
        const floatexp Xxr = Xr + xr[i];
        const floatexp Xxi = Xi + xi[i];
        test2[i] = test1[i];
        test1[i] = double(Xxr * Xxr + Xxi * Xxi);
        bGlitch[i] |= test1[i] < Xz;
        if (bGlitch[i]) break;
        if (test1[i] > m_nBailout2) break;
        floatexp xrn, xin;
        {
          xrn=(((((Xr.mul2())+xr[i])*xr[i])-(((Xi.mul2())+xi[i])*xi[i]))+cr[i]);
          xin=((((Xxr*xi[i])+(Xi*xr[i])).mul2())+ci[i]);
        }
        xr.val[i] = xrn.val;
        xr.exp[i] = xrn.exp;
        xi.val[i] = xin.val;
        xi.exp[i] = xin.exp;
      }
    }
    antal0 = antal;
    test10 = test1;
    test20 = test2;
    xr0 = xr;
    xi0 = xi;
    return true;
  }
  return false;
}

// EOF
