/*
Kalles Fraktaler 2
Copyright (C) 2013-2017 Karl Runmo
Copyright (C) 2017-2021 Claude Heiland-Allen

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

#ifndef KF_SETTINGS_H
#define KF_SETTINGS_H 1

#include <string>

#include "exr.h"

// this sets the maximum number of references per image
#define OLD_GLITCH 10000

// this sets the range of approximation terms
// storage is O(terms^2) for R2 fractals, O(terms) for C fractals
#define MIN_APPROX_TERMS 3
#define MAX_APPROX_TERMS 63

class Settings
{

private:

  double m_ZoomSize;
  int64_t m_MaxReferences;
  bool m_GlitchLowTolerance;
  bool m_ApproxLowTolerance;
  bool m_AutoApproxTerms;
  int64_t m_ApproxTerms;
  int64_t m_WindowWidth;
  int64_t m_WindowHeight;
  int64_t m_WindowTop;
  int64_t m_WindowLeft;
  int64_t m_WindowBottom;
  int64_t m_WindowRight;
  int64_t m_ImageWidth;
  int64_t m_ImageHeight;
  double m_ThreadsPerCore;
  int64_t m_ThreadsReserveCore;
  bool m_AnimateZoom;
  bool m_ArbitrarySize;
  bool m_ReuseReference;
  bool m_AutoSolveGlitches;
  bool m_Guessing;
  bool m_SolveGlitchNear;
  bool m_NoApprox;
  bool m_Mirror;
  bool m_LongDoubleAlways;
  bool m_FloatExpAlways;
  bool m_AutoIterations;
  bool m_ShowGlitches;
  bool m_NoReuseCenter;
  int64_t m_IsolatedGlitchNeighbourhood;
  int64_t m_JitterSeed;
  int64_t m_JitterShape;
  double m_JitterScale;
  bool m_Derivatives;
  bool m_ShowCrossHair;
  bool m_UseNanoMB1;
  bool m_UseNanoMB2;
  int64_t m_OrderM;
  int64_t m_OrderN;
  bool m_InteriorChecking;
  double m_RadiusScale;
  int64_t m_Shrink;
  bool m_HalfColour;
  bool m_SaveOverwrites;
  bool m_ThreadedReference;
  int64_t m_SIMDVectorSize;
  int64_t m_SIMDChunkSize;
  int64_t m_GlitchCenterMethod;
  bool m_UseOpenCL;
  int64_t m_OpenCLPlatform;
  int64_t m_EXRChannels;
  bool m_EXRParallel;
  bool m_SaveNewtonProgress;
  bool m_ExponentialMap;
  bool m_UseOpenGL;

public:

  Settings()
  : m_ZoomSize(2.0)
  , m_MaxReferences(10000)
  , m_GlitchLowTolerance(false)
  , m_ApproxLowTolerance(false)
  , m_AutoApproxTerms(true)
  , m_ApproxTerms(10)
  , m_WindowWidth(640)
  , m_WindowHeight(360)
  , m_WindowTop(-1)
  , m_WindowLeft(-1)
  , m_WindowBottom(-1)
  , m_WindowRight(-1)
  , m_ImageWidth(640)
  , m_ImageHeight(360)
  , m_ThreadsPerCore(1)
  , m_ThreadsReserveCore(0)
  , m_AnimateZoom(true)
  , m_ArbitrarySize(true)
  , m_ReuseReference(false)
  , m_AutoSolveGlitches(true)
  , m_Guessing(true)
  , m_SolveGlitchNear(false)
  , m_NoApprox(false)
  , m_Mirror(false)
  , m_LongDoubleAlways(false)
  , m_FloatExpAlways(false)
  , m_AutoIterations(true)
  , m_ShowGlitches(true)
  , m_NoReuseCenter(true)
  , m_IsolatedGlitchNeighbourhood(4)
  , m_JitterSeed(0)
  , m_JitterShape(0)
  , m_JitterScale(1)
  , m_Derivatives(false)
  , m_ShowCrossHair(false)
  , m_UseNanoMB1(false)
  , m_UseNanoMB2(false)
  , m_OrderM(16)
  , m_OrderN(16)
  , m_InteriorChecking(false)
  , m_RadiusScale(0.1)
  , m_Shrink(1)
  , m_HalfColour(false)
  , m_SaveOverwrites(false)
  , m_ThreadedReference(true)
  , m_SIMDVectorSize(2)
  , m_SIMDChunkSize(64)
  , m_GlitchCenterMethod(0)
  , m_UseOpenCL(false)
  , m_OpenCLPlatform(0)
  , m_EXRChannels(~0)
  , m_EXRParallel(true)
  , m_SaveNewtonProgress(false)
  , m_ExponentialMap(false)
  , m_UseOpenGL(false)
  { };

  bool FromText(const std::string &text);
  std::string ToText() const;

  bool OpenFile(const std::string &filename);
  bool SaveFile(const std::string &filename, bool overwrite) const;

  inline double GetZoomSize() const { return m_ZoomSize; };
  inline void   SetZoomSize(double z) { m_ZoomSize = z; };

  inline int64_t    GetMaxReferences() const { return m_MaxReferences; };
  inline void   SetMaxReferences(int64_t r) {
    if (r < 0) r = 0;
    if (r > OLD_GLITCH) r = OLD_GLITCH;
    m_MaxReferences = r;
  };

  inline bool   GetGlitchLowTolerance() const { return m_GlitchLowTolerance; };
  inline void   SetGlitchLowTolerance(bool b) { m_GlitchLowTolerance = b; };

  inline bool   GetApproxLowTolerance() const { return m_ApproxLowTolerance; };
  inline void   SetApproxLowTolerance(bool b) { m_ApproxLowTolerance = b; }

  inline bool   GetAutoApproxTerms() const { return m_AutoApproxTerms; };
  inline void   SetAutoApproxTerms(bool b) { m_AutoApproxTerms = b; };

  inline int64_t    GetApproxTerms() const { return m_ApproxTerms; };
  inline void   SetApproxTerms(int64_t t) {
    if (t < MIN_APPROX_TERMS) t = MIN_APPROX_TERMS;
    if (t > MAX_APPROX_TERMS) t = MAX_APPROX_TERMS;
    m_ApproxTerms = t;
  };

  inline int64_t    GetWindowWidth() const { return m_WindowWidth; };
  inline void   SetWindowWidth(int64_t w) { m_WindowWidth = w; };

  inline int64_t    GetWindowHeight() const { return m_WindowHeight; };
  inline void   SetWindowHeight(int64_t h) { m_WindowHeight = h; };

  inline int64_t    GetWindowTop() const { return m_WindowTop; };
  inline void   SetWindowTop(int64_t w) { m_WindowTop = w; };

  inline int64_t    GetWindowLeft() const { return m_WindowLeft; };
  inline void   SetWindowLeft(int64_t w) { m_WindowLeft = w; };

  inline int64_t    GetWindowBottom() const { return m_WindowBottom; };
  inline void   SetWindowBottom(int64_t w) { m_WindowBottom = w; };

  inline int64_t    GetWindowRight() const { return m_WindowRight; };
  inline void   SetWindowRight(int64_t w) { m_WindowRight = w; };

  inline int64_t    GetImageWidth() const { return m_ImageWidth; };
  inline void   SetImageWidth(int64_t w) { m_ImageWidth = w; };

  inline int64_t    GetImageHeight() const { return m_ImageHeight; };
  inline void   SetImageHeight(int64_t h) { m_ImageHeight = h; };

  inline double GetThreadsPerCore() const { return m_ThreadsPerCore; };
  inline void   SetThreadsPerCore(double t) { m_ThreadsPerCore = t; };

  inline int64_t    GetThreadsReserveCore() const { return m_ThreadsReserveCore; };
  inline void   SetThreadsReserveCore(int64_t t) { m_ThreadsReserveCore = t; };

  inline bool   GetAnimateZoom() const { return m_AnimateZoom; };
  inline void   SetAnimateZoom(bool b) { m_AnimateZoom = b; };

  inline bool   GetArbitrarySize() const { return m_ArbitrarySize; };
  inline void   SetArbitrarySize(bool b) { m_ArbitrarySize = b; };

  inline bool   GetReuseReference() const { return m_ReuseReference; };
  inline void   SetReuseReference(bool b) { m_ReuseReference = b; };

  inline bool   GetAutoSolveGlitches() const { return m_AutoSolveGlitches; };
  inline void   SetAutoSolveGlitches(bool b) { m_AutoSolveGlitches = b; };

  inline bool   GetGuessing() const { return m_Guessing; };
  inline void   SetGuessing(bool b) { m_Guessing = b; };

  inline bool   GetSolveGlitchNear() const { return m_SolveGlitchNear; };
  inline void   SetSolveGlitchNear(bool b) { m_SolveGlitchNear = b; };

  inline bool   GetNoApprox() const { return m_NoApprox; };
  inline void   SetNoApprox(bool b) { m_NoApprox = b; };

  inline bool   GetMirror() const { return m_Mirror; };
  inline void   SetMirror(bool b) { m_Mirror = b; };

  inline bool   GetLongDoubleAlways() const { return m_LongDoubleAlways; };
  inline void   SetLongDoubleAlways(bool b) { m_LongDoubleAlways = b; };

  inline bool   GetFloatExpAlways() const { return m_FloatExpAlways; };
  inline void   SetFloatExpAlways(bool b) { m_FloatExpAlways = b; };

  inline bool   GetAutoIterations() const { return m_AutoIterations; };
  inline void   SetAutoIterations(bool b) { m_AutoIterations = b; };

  inline bool   GetShowGlitches() const { return m_ShowGlitches; };
  inline void   SetShowGlitches(bool b) { m_ShowGlitches = b; };

  inline bool   GetNoReuseCenter() const { return m_NoReuseCenter; };
  inline void   SetNoReuseCenter(bool b) { m_NoReuseCenter = b; };

  inline int64_t    GetIsolatedGlitchNeighbourhood() const { return m_IsolatedGlitchNeighbourhood; };
  inline void   SetIsolatedGlitchNeighbourhood(int64_t n) { m_IsolatedGlitchNeighbourhood = n; };

  inline int64_t    GetJitterSeed() const { return m_JitterSeed; };
  inline void   SetJitterSeed(int64_t n) { m_JitterSeed = n; };

  inline int64_t    GetJitterShape() const { return m_JitterShape; };
  inline void   SetJitterShape(int64_t n) { m_JitterShape = 0 <= n && n <= 1 ? n : 0; };

  inline double GetJitterScale() const { return m_JitterScale; };
  inline void   SetJitterScale(double n) { m_JitterScale = n; };

  inline bool   GetDerivatives() const { return m_Derivatives; };
  inline void   SetDerivatives(bool b) { m_Derivatives = b; };

  inline bool   GetShowCrossHair() const { return m_ShowCrossHair; };
  inline void   SetShowCrossHair(bool b) { m_ShowCrossHair = b; };

  inline bool   GetUseNanoMB1() const { return m_UseNanoMB1; };
  inline void   SetUseNanoMB1(bool b) { m_UseNanoMB1 = b; };

  inline bool   GetUseNanoMB2() const { return m_UseNanoMB2; };
  inline void   SetUseNanoMB2(bool b) { m_UseNanoMB2 = b; };

  inline int64_t    GetOrderM() const { return m_OrderM; };
  inline void   SetOrderM(int64_t t)
  {
    if (t < 1) t = 1;
    if (t > MAX_APPROX_TERMS) t = MAX_APPROX_TERMS;
    m_OrderM = t;
  };

  inline int64_t    GetOrderN() const { return m_OrderN; };
  inline void   SetOrderN(int64_t t)
  {
    if (t < 1) t = 1;
    if (t > MAX_APPROX_TERMS) t = MAX_APPROX_TERMS;
    m_OrderN = t;
  };

  inline bool   GetInteriorChecking() const { return m_InteriorChecking; };
  inline void   SetInteriorChecking(bool b) { m_InteriorChecking = b; };

  inline double GetRadiusScale() const { return m_RadiusScale; };
  inline void   SetRadiusScale(double b) { m_RadiusScale = b; };

  inline int64_t    GetShrink() const { return m_Shrink; };
  inline void   SetShrink(int64_t n) { m_Shrink = 0 <= n && n <= 2 ? n : 1; };

  inline bool   GetHalfColour() const { return m_HalfColour; };
  inline void   SetHalfColour(bool b) { m_HalfColour = b; };

  inline bool   GetSaveOverwrites() const { return m_SaveOverwrites; };
  inline void   SetSaveOverwrites(bool b) { m_SaveOverwrites = b; };

  inline bool   GetThreadedReference() const { return m_ThreadedReference; };
  inline void   SetThreadedReference(bool b) { m_ThreadedReference = b; };

  inline int64_t    GetSIMDChunkSize() const { return m_SIMDChunkSize; };
  inline void   SetSIMDChunkSize(int64_t n) { m_SIMDChunkSize = 0 < n ? n : 1; };

  inline int64_t    GetSIMDVectorSize() const { return m_SIMDVectorSize; };
  inline void   SetSIMDVectorSize(int64_t n) {
    if (n > 1 << KF_SIMD) n = 1 << KF_SIMD;
    switch (n)
    {
      case 2: case 4: case 8: case 16: m_SIMDVectorSize = n; break;
      default: m_SIMDVectorSize = 1; break;
    }
  }

  inline int64_t GetGlitchCenterMethod() const { return m_GlitchCenterMethod; };
  inline void    SetGlitchCenterMethod(int64_t b) { m_GlitchCenterMethod = b; };
  inline bool    GetUseArgMinAbsZAsGlitchCenter() const { return m_GlitchCenterMethod == 1; }

  inline bool   GetUseOpenCL() const { return m_UseOpenCL; };
  inline void   SetUseOpenCL(bool b) { m_UseOpenCL = b; };

  inline int64_t GetOpenCLPlatform() const { return m_OpenCLPlatform; };
  inline void    SetOpenCLPlatform(int64_t n) { m_OpenCLPlatform = n; };

  inline EXRChannels GetEXRChannels() const { return unpack_exr_channels(m_EXRChannels); };
  inline void SetEXRChannels(const EXRChannels n) { m_EXRChannels = pack_exr_channels(n); };

  inline bool   GetEXRParallel() const { return m_EXRParallel; };
  inline void   SetEXRParallel(bool b) { m_EXRParallel = b; };

  inline bool   GetSaveNewtonProgress() const { return m_SaveNewtonProgress; };
  inline void   SetSaveNewtonProgress(bool b) { m_SaveNewtonProgress = b; };

  inline bool   GetExponentialMap() const { return m_ExponentialMap; };
  inline void   SetExponentialMap(bool b) { m_ExponentialMap = b; };

  inline bool   GetUseOpenGL() const { return m_UseOpenGL; };
  inline void   SetUseOpenGL(bool b) { m_UseOpenGL = b; };

};

#endif
