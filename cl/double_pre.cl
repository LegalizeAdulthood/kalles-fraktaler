void perturbation_double_loop
( __global const p_config   *g
, __global const double     *m_db_dxr
, __global const double     *m_db_dxi
, __global const double     *m_db_z
,                p_status_d *l
)
{
  const double Ar = g->g_FactorAR;
  const double Ai = g->g_FactorAI;
  bool no_g = g->g_real == 1.0 && g->g_imag == 1.0;
  const double cr = l->cr;
  const double ci = l->ci;
  double xr = l->xr;
  double xi = l->xi;
  // type C derivatives
  double dr = l->dr;
  double di = l->di;
  // type M derivatives
  double daa = l->daa;
  double dab = l->dab;
  double dba = l->dba;
  double dbb = l->dbb;
  double dxa = dr;
  double dxb = -di;
  double dya = di;
  double dyb = dr;
  // conditions
  double test1 = l->test1;
  double test2 = l->test2;
  long antal = l->antal;
  double Xxr = 0;
  double Xxi = 0;
  for (; antal < g->nMaxIter; antal++)
  {
    const double Xr = m_db_dxr[antal - g->nMinIter];
    const double Xi = m_db_dxi[antal - g->nMinIter];
    const double Xz = m_db_z  [antal - g->nMinIter];
    Xxr = Xr + xr;
    Xxi = Xi + xi;
    test2 = test1;
    if (no_g)
    {
      test1 = Xxr * Xxr + Xxi * Xxi;
    }
    else
    {
      test1 = g->g_real * Xxr * Xxr + g->g_imag * Xxi * Xxi;
    }
    if (test1 < Xz)
    {
      l->bGlitch = true;
      if (! g->m_bNoGlitchDetection)
        break;
    }
    if (test1 > g->m_nBailout2)
    {
      break;
    }
    const dcomplex X = {Xr, Xi}, x = {xr, xi}, Xx = {Xxr, Xxi};
    double xrn = 0;
    double xin = 0;
    dcomplex xn = {0, 0};
    const dcomplex d = {dr, di}, d0 = {daa, dba}; // FIXME matrix derivatives
    double drn = 0, din = 0;
    dcomplex dn = {0, 0};
    double dxan = 0, dxbn = 0, dyan = 0, dybn = 0;
// continued...
