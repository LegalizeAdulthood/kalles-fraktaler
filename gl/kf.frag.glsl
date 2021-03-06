// implement this
vec3 colour(void);

///=====================================================================
/// public API

uniform ivec2 KFP_ImageSize;

uniform sampler1D KFP_Palette;

uniform vec3 KFP_InteriorColor;

uniform bool KFP_ShowGlitches;

uniform uvec2 KFP_Iterations;
uniform uvec2 KFP_IterationsMin;
uniform uvec2 KFP_IterationsMax;

uniform uint KFP_JitterSeed;
uniform int KFP_JitterShape;
uniform float KFP_JitterScale;

uniform float KFP_IterDiv;
uniform float KFP_ColorOffset;

uniform int KFP_ColorMethod;
#define ColorMethod_Standard 0
#define ColorMethod_SquareRoot 1
#define ColorMethod_CubicRoot 2
#define ColorMethod_Logarithm 3
#define ColorMethod_Stretched 4
#define ColorMethod_DistanceLinear 5
#define ColorMethod_DEPlusStandard 6
#define ColorMethod_DistanceLog 7
#define ColorMethod_DistanceSqrt 8
#define ColorMethod_LogLog 9
#define ColorMethod_ATan 10
#define ColorMethod_FourthRoot 11

uniform int KFP_Differences;
#define Differences_Traditional 0
#define Differences_Forward3x3 1
#define Differences_Central3x3 2
#define Differences_Diagonal2x2 3
#define Differences_LeastSquares2x2 4
#define Differences_LeastSquares3x3 5
#define Differences_Laplacian3x3 6
#define Differences_Analytic 7

uniform float KFP_PhaseColorStrength;

uniform bool KFP_Smooth;
uniform bool KFP_Flat;
uniform bool KFP_InverseTransition;

uniform bool KFP_MultiWavesEnabled;
uniform bool KFP_MultiWavesBlend;
uniform int KFP_MultiWavesCount;
#define KFP_MultiWavesCountMax 32
uniform ivec3 KFP_MultiWaves[KFP_MultiWavesCountMax];

uniform bool KFP_Slopes;
uniform float KFP_SlopePower;
uniform float KFP_SlopeRatio;
uniform vec2 KFP_SlopeDir;

uniform sampler2D KFP_Texture;
uniform bool KFP_TextureEnabled;
uniform float KFP_TextureMerge;
uniform float KFP_TexturePower;
uniform float KFP_TextureRatio;

uniform bool KFP_sRGB;

uniform float KFP_ZoomLog2;

/// end of public API
///=====================================================================

#if __VERSION__ >= 330
layout(location = 0, index = 0) out vec4 Internal_Colour;
#else
#define Internal_Colour gl_FragColor
#endif

uniform usampler2D Internal_N1;
uniform usampler2D Internal_N0;
uniform sampler2D Internal_NF;
uniform sampler2D Internal_T;
uniform sampler2D Internal_DEX;
uniform sampler2D Internal_DEY;

uniform ivec2 Internal_TilePadding;
uniform ivec2 Internal_TileOrigin;
uniform ivec2 Internal_TileSize;

// hack to force explicit evaluation order
uniform float Internal_Zero;
float Internal_One = 1.0;
float EXACT(float a) { return Internal_One * a; }

#if __VERSION__ >= 400
float _builtin_ldexp(float a, int b) { return ldexp(a, b); }
float ldexp(float a, int b) { return _builtin_ldexp(a, b); }
#else
float ldexp(float a, int b) { return a * exp2(b); }
#endif

///=====================================================================
/// overload abs()

int _builtin_abs(int a) { return abs(a); }
ivec2 _builtin_abs(ivec2 a) { return abs(a); }
ivec3 _builtin_abs(ivec3 a) { return abs(a); }
ivec4 _builtin_abs(ivec4 a) { return abs(a); }
float _builtin_abs(float a) { return abs(a); }
vec2 _builtin_abs(vec2 a) { return abs(a); }
vec3 _builtin_abs(vec3 a) { return abs(a); }
vec4 _builtin_abs(vec4 a) { return abs(a); }
#if __VERSION__ >= 400
double _builtin_abs(double a) { return abs(a); }
dvec2 _builtin_abs(dvec2 a) { return abs(a); }
dvec3 _builtin_abs(dvec3 a) { return abs(a); }
dvec4 _builtin_abs(dvec4 a) { return abs(a); }
#endif

int abs(int a) { return _builtin_abs(a); }
ivec2 abs(ivec2 a) { return _builtin_abs(a); }
ivec3 abs(ivec3 a) { return _builtin_abs(a); }
ivec4 abs(ivec4 a) { return _builtin_abs(a); }
float abs(float a) { return _builtin_abs(a); }
vec2 abs(vec2 a) { return _builtin_abs(a); }
vec3 abs(vec3 a) { return _builtin_abs(a); }
vec4 abs(vec4 a) { return _builtin_abs(a); }
#if __VERSION__ >= 400
double abs(double a) { return _builtin_abs(a); }
dvec2 abs(dvec2 a) { return _builtin_abs(a); }
dvec3 abs(dvec3 a) { return _builtin_abs(a); }
dvec4 abs(dvec4 a) { return _builtin_abs(a); }
#endif

///=====================================================================
/// overload sqrt()

float _builtin_sqrt(float a) { return sqrt(a); }
vec2 _builtin_sqrt(vec2 a) { return sqrt(a); }
vec3 _builtin_sqrt(vec3 a) { return sqrt(a); }
vec4 _builtin_sqrt(vec4 a) { return sqrt(a); }
#if __VERSION__ >= 400
double _builtin_sqrt(double a) { return sqrt(a); }
dvec2 _builtin_sqrt(dvec2 a) { return sqrt(a); }
dvec3 _builtin_sqrt(dvec3 a) { return sqrt(a); }
dvec4 _builtin_sqrt(dvec4 a) { return sqrt(a); }
#endif

float sqrt(float a) { return _builtin_sqrt(a); }
vec2 sqrt(vec2 a) { return _builtin_sqrt(a); }
vec3 sqrt(vec3 a) { return _builtin_sqrt(a); }
vec4 sqrt(vec4 a) { return _builtin_sqrt(a); }
#if __VERSION__ >= 400
double sqrt(double a) { return _builtin_sqrt(a); }
dvec2 sqrt(dvec2 a) { return _builtin_sqrt(a); }
dvec3 sqrt(dvec3 a) { return _builtin_sqrt(a); }
dvec4 sqrt(dvec4 a) { return _builtin_sqrt(a); }
#endif

///=====================================================================
/// overload floor()

float _builtin_floor(float a) { return floor(a); }
vec2 _builtin_floor(vec2 a) { return floor(a); }
vec3 _builtin_floor(vec3 a) { return floor(a); }
vec4 _builtin_floor(vec4 a) { return floor(a); }
#if __VERSION__ >= 400
double _builtin_floor(double a) { return floor(a); }
dvec2 _builtin_floor(dvec2 a) { return floor(a); }
dvec3 _builtin_floor(dvec3 a) { return floor(a); }
dvec4 _builtin_floor(dvec4 a) { return floor(a); }
#endif

float floor(float a) { return _builtin_floor(a); }
vec2 floor(vec2 a) { return _builtin_floor(a); }
vec3 floor(vec3 a) { return _builtin_floor(a); }
vec4 floor(vec4 a) { return _builtin_floor(a); }
#if __VERSION__ >= 400
double floor(double a) { return _builtin_floor(a); }
dvec2 floor(dvec2 a) { return _builtin_floor(a); }
dvec3 floor(dvec3 a) { return _builtin_floor(a); }
dvec4 floor(dvec4 a) { return _builtin_floor(a); }
#endif

///=====================================================================
/// overload ceil()

float _builtin_ceil(float a) { return ceil(a); }
vec2 _builtin_ceil(vec2 a) { return ceil(a); }
vec3 _builtin_ceil(vec3 a) { return ceil(a); }
vec4 _builtin_ceil(vec4 a) { return ceil(a); }
#if __VERSION__ >= 400
double _builtin_ceil(double a) { return ceil(a); }
dvec2 _builtin_ceil(dvec2 a) { return ceil(a); }
dvec3 _builtin_ceil(dvec3 a) { return ceil(a); }
dvec4 _builtin_ceil(dvec4 a) { return ceil(a); }
#endif

float ceil(float a) { return _builtin_ceil(a); }
vec2 ceil(vec2 a) { return _builtin_ceil(a); }
vec3 ceil(vec3 a) { return _builtin_ceil(a); }
vec4 ceil(vec4 a) { return _builtin_ceil(a); }
#if __VERSION__ >= 400
double ceil(double a) { return _builtin_ceil(a); }
dvec2 ceil(dvec2 a) { return _builtin_ceil(a); }
dvec3 ceil(dvec3 a) { return _builtin_ceil(a); }
dvec4 ceil(dvec4 a) { return _builtin_ceil(a); }
#endif

///=====================================================================
/// overload exp()

float _builtin_exp(float a) { return exp(a); }
vec2 _builtin_exp(vec2 a) { return exp(a); }
vec3 _builtin_exp(vec3 a) { return exp(a); }
vec4 _builtin_exp(vec4 a) { return exp(a); }

float exp(float a) { return _builtin_exp(a); }
vec2 exp(vec2 a) { return _builtin_exp(a); }
vec3 exp(vec3 a) { return _builtin_exp(a); }
vec4 exp(vec4 a) { return _builtin_exp(a); }

///=====================================================================
/// overload log()

float _builtin_log(float a) { return log(a); }
vec2 _builtin_log(vec2 a) { return log(a); }
vec3 _builtin_log(vec3 a) { return log(a); }
vec4 _builtin_log(vec4 a) { return log(a); }

float log(float a) { return _builtin_log(a); }
vec2 log(vec2 a) { return _builtin_log(a); }
vec3 log(vec3 a) { return _builtin_log(a); }
vec4 log(vec4 a) { return _builtin_log(a); }

///=====================================================================
/// overload sin()

float _builtin_sin(float a) { return sin(a); }
vec2 _builtin_sin(vec2 a) { return sin(a); }
vec3 _builtin_sin(vec3 a) { return sin(a); }
vec4 _builtin_sin(vec4 a) { return sin(a); }

float sin(float a) { return _builtin_sin(a); }
vec2 sin(vec2 a) { return _builtin_sin(a); }
vec3 sin(vec3 a) { return _builtin_sin(a); }
vec4 sin(vec4 a) { return _builtin_sin(a); }

///=====================================================================
/// overload cos()

float _builtin_cos(float a) { return cos(a); }
vec2 _builtin_cos(vec2 a) { return cos(a); }
vec3 _builtin_cos(vec3 a) { return cos(a); }
vec4 _builtin_cos(vec4 a) { return cos(a); }

float cos(float a) { return _builtin_cos(a); }
vec2 cos(vec2 a) { return _builtin_cos(a); }
vec3 cos(vec3 a) { return _builtin_cos(a); }
vec4 cos(vec4 a) { return _builtin_cos(a); }

///=====================================================================
/// overload tan()

float _builtin_tan(float a) { return tan(a); }
vec2 _builtin_tan(vec2 a) { return tan(a); }
vec3 _builtin_tan(vec3 a) { return tan(a); }
vec4 _builtin_tan(vec4 a) { return tan(a); }

float tan(float a) { return _builtin_tan(a); }
vec2 tan(vec2 a) { return _builtin_tan(a); }
vec3 tan(vec3 a) { return _builtin_tan(a); }
vec4 tan(vec4 a) { return _builtin_tan(a); }

///=====================================================================
/// overload sinh()

float _builtin_sinh(float a) { return sinh(a); }
vec2 _builtin_sinh(vec2 a) { return sinh(a); }
vec3 _builtin_sinh(vec3 a) { return sinh(a); }
vec4 _builtin_sinh(vec4 a) { return sinh(a); }

float sinh(float a) { return _builtin_sinh(a); }
vec2 sinh(vec2 a) { return _builtin_sinh(a); }
vec3 sinh(vec3 a) { return _builtin_sinh(a); }
vec4 sinh(vec4 a) { return _builtin_sinh(a); }

///=====================================================================
/// overload cosh()

float _builtin_cosh(float a) { return cosh(a); }
vec2 _builtin_cosh(vec2 a) { return cosh(a); }
vec3 _builtin_cosh(vec3 a) { return cosh(a); }
vec4 _builtin_cosh(vec4 a) { return cosh(a); }

float cosh(float a) { return _builtin_cosh(a); }
vec2 cosh(vec2 a) { return _builtin_cosh(a); }
vec3 cosh(vec3 a) { return _builtin_cosh(a); }
vec4 cosh(vec4 a) { return _builtin_cosh(a); }

///=====================================================================
/// overload tanh()

float _builtin_tanh(float a) { return tanh(a); }
vec2 _builtin_tanh(vec2 a) { return tanh(a); }
vec3 _builtin_tanh(vec3 a) { return tanh(a); }
vec4 _builtin_tanh(vec4 a) { return tanh(a); }

float tanh(float a) { return _builtin_tanh(a); }
vec2 tanh(vec2 a) { return _builtin_tanh(a); }
vec3 tanh(vec3 a) { return _builtin_tanh(a); }
vec4 tanh(vec4 a) { return _builtin_tanh(a); }

///=====================================================================
/// overload asin()

float _builtin_asin(float a) { return asin(a); }
vec2 _builtin_asin(vec2 a) { return asin(a); }
vec3 _builtin_asin(vec3 a) { return asin(a); }
vec4 _builtin_asin(vec4 a) { return asin(a); }

float asin(float a) { return _builtin_asin(a); }
vec2 asin(vec2 a) { return _builtin_asin(a); }
vec3 asin(vec3 a) { return _builtin_asin(a); }
vec4 asin(vec4 a) { return _builtin_asin(a); }

///=====================================================================
/// overload acos()

float _builtin_acos(float a) { return acos(a); }
vec2 _builtin_acos(vec2 a) { return acos(a); }
vec3 _builtin_acos(vec3 a) { return acos(a); }
vec4 _builtin_acos(vec4 a) { return acos(a); }

float acos(float a) { return _builtin_acos(a); }
vec2 acos(vec2 a) { return _builtin_acos(a); }
vec3 acos(vec3 a) { return _builtin_acos(a); }
vec4 acos(vec4 a) { return _builtin_acos(a); }

///=====================================================================
/// overload atan()

float _builtin_atan(float a) { return atan(a); }
vec2 _builtin_atan(vec2 a) { return atan(a); }
vec3 _builtin_atan(vec3 a) { return atan(a); }
vec4 _builtin_atan(vec4 a) { return atan(a); }

float _builtin_atan(float a, float b) { return atan(a, b); }
vec2 _builtin_atan(vec2 a, vec2 b) { return atan(a, b); }
vec3 _builtin_atan(vec3 a, vec3 b) { return atan(a, b); }
vec4 _builtin_atan(vec4 a, vec4 b) { return atan(a, b); }

float atan(float a) { return _builtin_atan(a); }
vec2 atan(vec2 a) { return _builtin_atan(a); }
vec3 atan(vec3 a) { return _builtin_atan(a); }
vec4 atan(vec4 a) { return _builtin_atan(a); }

float atan(float a, float b) { return _builtin_atan(a, b); }
vec2 atan(vec2 a, vec2 b) { return _builtin_atan(a, b); }
vec3 atan(vec3 a, vec3 b) { return _builtin_atan(a, b); }
vec4 atan(vec4 a, vec4 b) { return _builtin_atan(a, b); }

///=====================================================================
/// overload asinh()

float _builtin_asinh(float a) { return asinh(a); }
vec2 _builtin_asinh(vec2 a) { return asinh(a); }
vec3 _builtin_asinh(vec3 a) { return asinh(a); }
vec4 _builtin_asinh(vec4 a) { return asinh(a); }

float asinh(float a) { return _builtin_asinh(a); }
vec2 asinh(vec2 a) { return _builtin_asinh(a); }
vec3 asinh(vec3 a) { return _builtin_asinh(a); }
vec4 asinh(vec4 a) { return _builtin_asinh(a); }

///=====================================================================
/// overload acosh()

float _builtin_acosh(float a) { return acosh(a); }
vec2 _builtin_acosh(vec2 a) { return acosh(a); }
vec3 _builtin_acosh(vec3 a) { return acosh(a); }
vec4 _builtin_acosh(vec4 a) { return acosh(a); }

float acosh(float a) { return _builtin_acosh(a); }
vec2 acosh(vec2 a) { return _builtin_acosh(a); }
vec3 acosh(vec3 a) { return _builtin_acosh(a); }
vec4 acosh(vec4 a) { return _builtin_acosh(a); }

///=====================================================================
/// overload atanh()

float _builtin_atanh(float a) { return atanh(a); }
vec2 _builtin_atanh(vec2 a) { return atanh(a); }
vec3 _builtin_atanh(vec3 a) { return atanh(a); }
vec4 _builtin_atanh(vec4 a) { return atanh(a); }

float atanh(float a) { return _builtin_atanh(a); }
vec2 atanh(vec2 a) { return _builtin_atanh(a); }
vec3 atanh(vec3 a) { return _builtin_atanh(a); }
vec4 atanh(vec4 a) { return _builtin_atanh(a); }

///=====================================================================
/// overload pow()

float _builtin_pow(float a, float b) { return pow(a, b); }
vec2 _builtin_pow(vec2 a, vec2 b) { return pow(a, b); }
vec3 _builtin_pow(vec3 a, vec3 b) { return pow(a, b); }
vec4 _builtin_pow(vec4 a, vec4 b) { return pow(a, b); }

float pow(float a, float b) { return _builtin_pow(a, b); }
vec2 pow(vec2 a, vec2 b) { return _builtin_pow(a, b); }
vec3 pow(vec3 a, vec3 b) { return _builtin_pow(a, b); }
vec4 pow(vec4 a, vec4 b) { return _builtin_pow(a, b); }

///=====================================================================
/// overload max()

float _builtin_max(float a, float b) { return max(a, b); }
vec2 _builtin_max(vec2 a, vec2 b) { return max(a, b); }
vec3 _builtin_max(vec3 a, vec3 b) { return max(a, b); }
vec4 _builtin_max(vec4 a, vec4 b) { return max(a, b); }
vec2 _builtin_max(vec2 a, float b) { return max(a, b); }
vec3 _builtin_max(vec3 a, float b) { return max(a, b); }
vec4 _builtin_max(vec4 a, float b) { return max(a, b); }
int _builtin_max(int a, int b) { return max(a, b); }
ivec2 _builtin_max(ivec2 a, ivec2 b) { return max(a, b); }
ivec3 _builtin_max(ivec3 a, ivec3 b) { return max(a, b); }
ivec4 _builtin_max(ivec4 a, ivec4 b) { return max(a, b); }
ivec2 _builtin_max(ivec2 a, int b) { return max(a, b); }
ivec3 _builtin_max(ivec3 a, int b) { return max(a, b); }
ivec4 _builtin_max(ivec4 a, int b) { return max(a, b); }
uint _builtin_max(uint a, uint b) { return max(a, b); }
uvec2 _builtin_max(uvec2 a, uvec2 b) { return max(a, b); }
uvec3 _builtin_max(uvec3 a, uvec3 b) { return max(a, b); }
uvec4 _builtin_max(uvec4 a, uvec4 b) { return max(a, b); }
uvec2 _builtin_max(uvec2 a, uint b) { return max(a, b); }
uvec3 _builtin_max(uvec3 a, uint b) { return max(a, b); }
uvec4 _builtin_max(uvec4 a, uint b) { return max(a, b); }
#if __VERSION__ >= 400
double _builtin_max(double a, double b) { return max(a, b); }
dvec2 _builtin_max(dvec2 a, dvec2 b) { return max(a, b); }
dvec3 _builtin_max(dvec3 a, dvec3 b) { return max(a, b); }
dvec4 _builtin_max(dvec4 a, dvec4 b) { return max(a, b); }
dvec2 _builtin_max(dvec2 a, double b) { return max(a, b); }
dvec3 _builtin_max(dvec3 a, double b) { return max(a, b); }
dvec4 _builtin_max(dvec4 a, double b) { return max(a, b); }
#endif

float max(float a, float b) { return _builtin_max(a, b); }
vec2 max(vec2 a, vec2 b) { return _builtin_max(a, b); }
vec3 max(vec3 a, vec3 b) { return _builtin_max(a, b); }
vec4 max(vec4 a, vec4 b) { return _builtin_max(a, b); }
vec2 max(vec2 a, float b) { return _builtin_max(a, b); }
vec3 max(vec3 a, float b) { return _builtin_max(a, b); }
vec4 max(vec4 a, float b) { return _builtin_max(a, b); }
int max(int a, int b) { return _builtin_max(a, b); }
ivec2 max(ivec2 a, ivec2 b) { return _builtin_max(a, b); }
ivec3 max(ivec3 a, ivec3 b) { return _builtin_max(a, b); }
ivec4 max(ivec4 a, ivec4 b) { return _builtin_max(a, b); }
ivec2 max(ivec2 a, int b) { return _builtin_max(a, b); }
ivec3 max(ivec3 a, int b) { return _builtin_max(a, b); }
ivec4 max(ivec4 a, int b) { return _builtin_max(a, b); }
uint max(uint a, uint b) { return _builtin_max(a, b); }
uvec2 max(uvec2 a, uvec2 b) { return _builtin_max(a, b); }
uvec3 max(uvec3 a, uvec3 b) { return _builtin_max(a, b); }
uvec4 max(uvec4 a, uvec4 b) { return _builtin_max(a, b); }
uvec2 max(uvec2 a, uint b) { return _builtin_max(a, b); }
uvec3 max(uvec3 a, uint b) { return _builtin_max(a, b); }
uvec4 max(uvec4 a, uint b) { return _builtin_max(a, b); }
#if __VERSION__ >= 400
double max(double a, double b) { return _builtin_max(a, b); }
dvec2 max(dvec2 a, dvec2 b) { return _builtin_max(a, b); }
dvec3 max(dvec3 a, dvec3 b) { return _builtin_max(a, b); }
dvec4 max(dvec4 a, dvec4 b) { return _builtin_max(a, b); }
dvec2 max(dvec2 a, double b) { return _builtin_max(a, b); }
dvec3 max(dvec3 a, double b) { return _builtin_max(a, b); }
dvec4 max(dvec4 a, double b) { return _builtin_max(a, b); }
#endif

///=====================================================================
/// BEGIN qd-2.3.22+dfsg.1
/// via apt-get source qd on Debian Bullseye
///=====================================================================

#if __VERSION__ >= 400
#define QD_FMS(a,b,s) fma(a,b,-s)
#endif

#define QD_IEEE_ADD

struct float49 { float x[2]; };
float49 float49_() { return float49(float[2](0.0, 0.0)); }
float49 float49_(float hi) { return float49(float[2](hi, 0.0)); }
float49 float49_(float hi, float lo) { return float49(float[2](hi, lo)); }

///=====================================================================
/// qd-2.3.22+dfsg.1/COPYING
///=====================================================================

/*
This work was supported by the Director, Office of Science, Division
of Mathematical, Information, and Computational Sciences of the
U.S. Department of Energy under contract numbers DE-AC03-76SF00098 and
DE-AC02-05CH11231.

Copyright (c) 2003-2009, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from U.S. Dept. of Energy) All rights reserved.

By downloading or using this software you are agreeing to the modified
BSD license that is in file "BSD-LBNL-License.doc" in the main ARPREC
directory. If you wish to use the software for commercial purposes
please contact the Technology Transfer Department at TTD@lbl.gov or
call 510-286-6457."
*/

///=====================================================================
/// qd-2.3.22+dfsg.1/BSD-LBNL-License.doc
///=====================================================================

/*
1. Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

(1) Redistributions of source code must retain the copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the copyright notice,
this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

3. You are under no obligation whatsoever to provide any bug fixes,
patches, or upgrades to the features, functionality or performance of
the source code ("Enhancements") to anyone; however, if you choose to
make your Enhancements available either publicly, or directly to
Lawrence Berkeley National Laboratory, without imposing a separate
written license agreement for such Enhancements, then you hereby grant
the following license: a non-exclusive, royalty-free perpetual license
to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements
or derivative works thereof, in binary and source code form.
*/

///=====================================================================
/// qd-2.3.22+dfsg.1/include/qd/inline.h
///=====================================================================

/*
 * include/inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file contains the basic functions used both by double-double
 * and quad-double package.  These are declared as inline functions as
 * they are the smallest building blocks of the double-double and
 * quad-double arithmetic.
 */

#define _QD_SPLITTER     4097.0       // = 2^12 + 1
#define _QD_SPLIT_THRESH 4.1538375e34 // 2^115

const float _d_nan = 0.0 / 0.0;
const float _d_inf = 1.0 / 0.0;

/*
ff :: Double -> (Float, Float)
ff d =
  let hi :: Float
      hi = realToFrac d
  in (hi, realToFrac (d - realToFrac hi))
*/
const float49 f49_nan = float49(float[2](_d_nan, _d_nan));
const float49 f49_inf = float49(float[2](_d_inf, _d_inf));
const float49 f49_0 = float49(float[2](0.0, 0.0));
const float49 f49_1 = float49(float[2](1.0, 0.0));
const float49 f49_e = float49(float[2](2.7182817,8.2548404e-8));
const float49 f49_log2 = float49(float[2](0.6931472,-1.9046542e-9));
const float49 f49_log10 = float49(float[2](2.3025851,-3.1975436e-8));
const float49 f49_2pi = float49(float[2](6.2831855,-1.7484555e-7));
const float49 f49_pi = float49(float[2](3.1415927,-8.742278e-8));
const float49 f49_3pi4 = float49(float[2](2.3561945,-5.9624403e-9));
const float49 f49_pi2 = float49(float[2](1.5707964,-4.371139e-8));
const float49 f49_pi4 = float49(float[2](0.7853982,-2.1855694e-8));
const float49 f49_pi16 = float49(float[2](0.19634955,-5.4639235e-9));

const float f49_eps = 1.4210855e-14; // 2^-46
const float f49_min_normalized = 1.9721523e-31; // 2^(-126 + 24)
const float49 f49_max = float49(float[2](3.4028235e38, 2.0282408e31));
const float49 f49_safe_max = float49(float[2](3.401993e38, 2.0282408e31));
const int f49_ndigits = 13;

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
float quick_two_sum(float a, float b, out float err) {
  float s = EXACT(a + b);
  err = EXACT(b - EXACT(s - a));
  return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
float quick_two_diff(float a, float b, out float err) {
  float s = EXACT(a - b);
  err = EXACT(EXACT(a - s) - b);
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
float two_sum(float a, float b, out float err) {
  float s = EXACT(a + b);
  float bb = EXACT(s - a);
  err = EXACT(EXACT(a - EXACT(s - bb)) + EXACT(b - bb));
  return s;
}

/* Computes fl(a-b) and err(a-b).  */
float two_diff(float a, float b, out float err) {
  float s = EXACT(a - b);
  float bb = EXACT(s - a);
  err = EXACT(EXACT(a - EXACT(s - bb)) - EXACT(b + bb));
  return s;
}

#ifndef QD_FMS
/* Computes high word and lo word of a */
void split(float a, out float hi, out float lo) {
  float temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = EXACT(temp - EXACT(temp - a));
    lo = EXACT(a - hi);
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a;
    hi = EXACT(temp - EXACT(temp - a));
    lo = EXACT(a - hi);
  }
}
#endif

/* Computes fl(a*b) and err(a*b). */
float two_prod(float a, float b, out float err) {
#ifdef QD_FMS
  float p = EXACT(a * b);
  err = QD_FMS(a, b, p);
  return p;
#else
  float a_hi, a_lo, b_hi, b_lo;
  float p = EXACT(a * b);
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = EXACT(EXACT(EXACT(a_hi * b_hi) - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
float two_sqr(float a, out float err) {
#ifdef QD_FMS
  float p = EXACT(a * a);
  err = QD_FMS(a, a, p);
  return p;
#else
  float hi, lo;
  float q = EXACT(a * a);
  split(a, hi, lo);
  err = EXACT(EXACT(EXACT(hi * hi) - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

/* Computes the nearest integer to d. */
float nint(float d) {
  if (d == floor(d))
    return d;
  return floor(d + 0.5);
}

/* Computes the truncated integer. */
float aint(float d) {
  return (d >= 0.0) ? floor(d) : ceil(d);
}

/* These are provided to give consistent
   interface for float with float-float and quad-float. */
void sincosh(float t, out float sinh_t, out float cosh_t) {
  sinh_t = sinh(t);
  cosh_t = cosh(t);
}

float sqr(float t) {
  return t * t;
}

float to_float(float a) { return a; }
int    to_int(float a) { return int(a); }

///=====================================================================
/// qd-2.3.22+dfsg.1/src/dd_const.cpp
///=====================================================================

/*
 * src/dd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 */

///=====================================================================
/// qd-2.3.22+dfsg.1/include/qd/dd_inline.h
///=====================================================================

/*
 * include/dd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the float-float
 * arithmetic package.
 */

/*********** Additions ************/
/* float-float = float + float */
float49 f49_add(float a, float b) {
  float s, e;
  s = two_sum(a, b, e);
  return float49_(s, e);
}

/* float-float + float */
float49 add(float49 a, float b) {
  float s1, s2;
  s1 = two_sum(a.x[0], b, s2);
  s2 = EXACT(s2 + a.x[1]);
  s1 = quick_two_sum(s1, s2, s2);
  return float49_(s1, s2);
}

/* float-float + float-float */
float49 ieee_add(float49 a, float49 b) {
  /* This one satisfies IEEE style error bound,
     due to K. Briggs and W. Kahan.                   */
  float s1, s2, t1, t2;

  s1 = two_sum(a.x[0], b.x[0], s2);
  t1 = two_sum(a.x[1], b.x[1], t2);
  s2 = EXACT(s2 + t1);
  s1 = quick_two_sum(s1, s2, s2);
  s2 = EXACT(s2 + t2);
  s1 = quick_two_sum(s1, s2, s2);
  return float49_(s1, s2);
}

float49 sloppy_add(float49 a, float49 b) {
  /* This is the less accurate version ... obeys Cray-style
     error bound. */
  float s, e;

  s = two_sum(a.x[0], b.x[0], e);
  e = EXACT(e + EXACT(a.x[1] + b.x[1]));
  s = quick_two_sum(s, e, e);
  return float49_(s, e);
}

float49 add(float49 a, float49 b) {
#ifndef QD_IEEE_ADD
  return sloppy_add(a, b);
#else
  return ieee_add(a, b);
#endif
}

/* float + float-float */
float49 add(float a, float49 b) {
  return add(b, a);
}


/*********** Self-Additions ************/
/* float-float += float */
void add_set(inout float49 self, float a) {
  float s1, s2;
  s1 = two_sum(self.x[0], a, s2);
  s2 = EXACT(s2 + self.x[1]);
  self.x[0] = quick_two_sum(s1, s2, self.x[1]);
}

/* float-float += float-float */
void add_set(inout float49 self, float49 a) {
#ifndef QD_IEEE_ADD
  float s, e;
  s = two_sum(self.x[0], a.x[0], e);
  e = EXACT(e + self.x[1]);
  e = EXACT(e + a.x[1]);
  self.x[0] = quick_two_sum(s, e, self.x[1]);
#else
  float s1, s2, t1, t2;
  s1 = two_sum(self.x[0], a.x[0], s2);
  t1 = two_sum(self.x[1], a.x[1], t2);
  s2 = EXACT(s2 + t1);
  s1 = quick_two_sum(s1, s2, s2);
  s2 = EXACT(s2 + t2);
  self.x[0] = quick_two_sum(s1, s2, self.x[1]);
#endif
}

/*********** Subtractions ************/
/* float-float = float - float */
float49 f49_sub(float a, float b) {
  float s, e;
  s = two_diff(a, b, e);
  return float49_(s, e);
}

/* float-float - float */
float49 sub(float49 a, float b) {
  float s1, s2;
  s1 = two_diff(a.x[0], b, s2);
  s2 = EXACT(s2 + a.x[1]);
  s1 = quick_two_sum(s1, s2, s2);
  return float49_(s1, s2);
}

/* float-float - float-float */
float49 sub(float49 a, float49 b) {
#ifndef QD_IEEE_ADD
  float s, e;
  s = two_diff(a.x[0], b.x[0], e);
  e = EXACT(e + a.x[1]);
  e = EXACT(e - b.x[1]);
  s = quick_two_sum(s, e, e);
  return float49_(s, e);
#else
  float s1, s2, t1, t2;
  s1 = two_diff(a.x[0], b.x[0], s2);
  t1 = two_diff(a.x[1], b.x[1], t2);
  s2 = EXACT(s2 + t1);
  s1 = quick_two_sum(s1, s2, s2);
  s2 = EXACT(s2 + t2);
  s1 = quick_two_sum(s1, s2, s2);
  return float49_(s1, s2);
#endif
}

/* float - float-float */
float49 sub(float a, float49 b) {
  float s1, s2;
  s1 = two_diff(a, b.x[0], s2);
  s2 = EXACT(s2 - b.x[1]);
  s1 = quick_two_sum(s1, s2, s2);
  return float49_(s1, s2);
}

/*********** Self-Subtractions ************/
/* float-float -= float */
void sub_set(inout float49 self, float a) {
  float s1, s2;
  s1 = two_diff(self.x[0], a, s2);
  s2 = EXACT(s2 + self.x[1]);
  self.x[0] = quick_two_sum(s1, s2, self.x[1]);
}

/* float-float -= float-float */
void sub_set(inout float49 self, float49 a) {
#ifndef QD_IEEE_ADD
  float s, e;
  s = two_diff(self.x[0], a.x[0], e);
  e = EXACT(e + self.x[1]);
  e = EXACT(e - a.x[1];
  self.x[0] = quick_two_sum(s, e, self.x[1]);
#else
  float s1, s2, t1, t2;
  s1 = two_diff(self.x[0], a.x[0], s2);
  t1 = two_diff(self.x[1], a.x[1], t2);
  s2 = EXACT(s2 + t1);
  s1 = quick_two_sum(s1, s2, s2);
  s2 = EXACT(s2 + t2);
  self.x[0] = quick_two_sum(s1, s2, self.x[1]);
#endif
}

/*********** Unary Minus ***********/
float49 neg(float49 a) {
  return float49_(-a.x[0], -a.x[1]);
}

/*********** Multiplications ************/
/* float-float = float * float */
float49 f49_mul(float a, float b) {
  float p, e;
  p = two_prod(a, b, e);
  return float49_(p, e);
}

/* float-float * (2.0 ^ exp) */
float49 ldexp(float49 a, int exp) {
  return float49_(ldexp(a.x[0], exp), ldexp(a.x[1], exp));
}

/* float-float * float,  where float is a power of 2. */
float49 mul_pwr2(float49 a, float b) {
  return float49_(a.x[0] * b, a.x[1] * b);
}

/* float-float * float */
float49 mul(float49 a, float b) {
  float p1, p2;

  p1 = two_prod(a.x[0], b, p2);
  p2 = EXACT(p2 + EXACT(a.x[1] * b));
  p1 = quick_two_sum(p1, p2, p2);
  return float49_(p1, p2);
}

/* float-float * float-float */
float49 mul(float49 a, float49 b) {
  float p1, p2;

  p1 = two_prod(a.x[0], b.x[0], p2);
  p2 = EXACT(p2 + EXACT(a.x[0] * b.x[1] + a.x[1] * b.x[0]));
  p1 = quick_two_sum(p1, p2, p2);
  return float49_(p1, p2);
}

/* float * float-float */
float49 mul(float a, float49 b) {
  return mul(b, a);
}

/*********** Self-Multiplications ************/
/* float-float *= float */
void mul_set(inout float49 self, float a) {
  float p1, p2;
  p1 = two_prod(self.x[0], a, p2);
  p2 = EXACT(p2 + EXACT(self.x[1] * a));
  self.x[0] = quick_two_sum(p1, p2, self.x[1]);
}

/* float-float *= float-float */
void mul_set(inout float49 self, float49 a) {
  float p1, p2;
  p1 = two_prod(self.x[0], a.x[0], p2);
  p2 = EXACT(p2 + EXACT(a.x[1] * self.x[0]));
  p2 = EXACT(p2 + EXACT(a.x[0] * self.x[1]));
  self.x[0] = quick_two_sum(p1, p2, self.x[1]);
}

/*********** Divisions ************/
float49 f49_div(float a, float b) {
  float q1, q2;
  float p1, p2;
  float s, e;

  q1 = EXACT(a / b);

  /* Compute  a - q1 * b */
  p1 = two_prod(q1, b, p2);
  s = two_diff(a, p1, e);
  e = EXACT(e - p2);

  /* get next approximation */
  q2 = EXACT(EXACT(s + e) / b);

  s = quick_two_sum(q1, q2, e);

  return float49_(s, e);
}

/* float-float / float */
float49 div(float49 a, float b) {

  float q1, q2;
  float p1, p2;
  float s, e;
  float49 r;

  q1 = EXACT(a.x[0] / b);   /* approximate quotient. */

  /* Compute  this - q1 * d */
  p1 = two_prod(q1, b, p2);
  s = two_diff(a.x[0], p1, e);
  e = EXACT(e + a.x[1]);
  e = EXACT(e - p2);

  /* get next approximation. */
  q2 = EXACT(EXACT(s + e) / b);

  /* renormalize */
  r.x[0] = quick_two_sum(q1, q2, r.x[1]);

  return r;
}

float49 sloppy_div(float49 a, float49 b) {
  float s1, s2;
  float q1, q2;
  float49 r;

  q1 = EXACT(a.x[0] / b.x[0]);  /* approximate quotient */

  /* compute  this - q1 * dd */
  r = mul(b, q1);
  s1 = two_diff(a.x[0], r.x[0], s2);
  s2 = EXACT(s2 - r.x[1]);
  s2 = EXACT(s2 + a.x[1]);

  /* get next approximation */
  q2 = EXACT(EXACT(s1 + s2) / b.x[0]);

  /* renormalize */
  r.x[0] = quick_two_sum(q1, q2, r.x[1]);
  return r;
}

float49 accurate_div(float49 a, float49 b) {
  float q1, q2, q3;
  float49 r;

  q1 = EXACT(a.x[0] / b.x[0]);  /* approximate quotient */

  r = sub(a, mul(q1, b));

  q2 = EXACT(r.x[0] / b.x[0]);
  sub_set(r, mul(q2, b));

  q3 = EXACT(r.x[0] / b.x[0]);

  q1 = quick_two_sum(q1, q2, q2);
  r = add(float49_(q1, q2), q3);
  return r;
}

/* float-float / float-float */
float49 div(float49 a, float49 b) {
#ifdef QD_SLOPPY_DIV
  return sloppy_div(a, b);
#else
  return accurate_div(a, b);
#endif
}

/* float / float-float */
float49 div(float a, float49 b) {
  return div(float49_(a), b);
}

float49 inv(float49 a) {
  return div(1.0, a);
}

/*********** Self-Divisions ************/
/* float-float /= float */
void div_set(inout float49 self, float a) {
  self = div(self, a);
}

/* float-float /= float-float */
void div_set(inout float49 self, float49 a) {
  self = div(self, a);
}

/*********** Squaring **********/
float49 sqr(float49 a) {
  float p1, p2;
  float s1, s2;
  p1 = two_sqr(a.x[0], p2);
  p2 = EXACT(p2 + EXACT(2.0 * a.x[0] * a.x[1]));
  p2 = EXACT(p2 + EXACT(a.x[1] * a.x[1]));
  s1 = quick_two_sum(p1, p2, s2);
  return float49_(s1, s2);
}

float49 f49_sqr(float a) {
  float p1, p2;
  p1 = two_sqr(a, p2);
  return float49_(p1, p2);
}

/*********** Micellaneous ************/
/*  this == 0 */
bool is_zero(float49 self) {
  return (self.x[0] == 0.0);
}

/*  this == 1 */
bool is_one(float49 self) {
  return (self.x[0] == 1.0 && self.x[1] == 0.0);
}

/*  this > 0 */
bool is_positive(float49 self) {
  return (self.x[0] > 0.0);
}

/* this < 0 */
bool is_negative(float49 self) {
  return (self.x[0] < 0.0);
}

/* Absolute value */
float49 abs(float49 a) {
  return (a.x[0] < 0.0) ? neg(a) : a;
}

float49 fabs(float49 a) {
  return abs(a);
}


/* Computes the n-th power of a float-float number.
   NOTE:  0^0 causes an error.                         */
float49 npwr(float49 a, int n) {

  if (n == 0) {
    if (is_zero(a)) {
      return f49_nan;
    }
    return f49_1;
  }

  float49 r = a;
  float49 s = f49_1;
  int N = abs(n);

  if (N > 1) {
    /* Use binary exponentiation */
    while (N > 0) {
      if (N % 2 == 1) {
        mul_set(s, r);
      }
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }
  } else {
    s = r;
  }

  /* Compute the reciprocal if n is negative. */
  if (n < 0)
    return inv(s);

  return s;
}

/********** Exponentiation **********/
float49 pow(float49 a, int n) {
  return npwr(a, n);
}


/*********** Assignments ************/
/* float-float = float */
void set(out float49 self, float a) {
  self.x[0] = a;
  self.x[1] = 0.0;
}

/*********** Equality Comparisons ************/
/* float-float == float */
bool eq(float49 a, float b) {
  return (a.x[0] == b && a.x[1] == 0.0);
}

/* float-float == float-float */
bool eq(float49 a, float49 b) {
  return (a.x[0] == b.x[0] && a.x[1] == b.x[1]);
}

/* float == float-float */
bool eq(float a, float49 b) {
  return (a == b.x[0] && b.x[1] == 0.0);
}

/*********** Greater-Than Comparisons ************/
/* float-float > float */
bool gt(float49 a, float b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] > 0.0));
}

/* float-float > float-float */
bool gt(float49 a, float49 b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] > b.x[1]));
}

/* float > float-float */
bool gt(float a, float49 b) {
  return (a > b.x[0] || (a == b.x[0] && b.x[1] < 0.0));
}

/*********** Less-Than Comparisons ************/
/* float-float < float */
bool lt(float49 a, float b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] < 0.0));
}

/* float-float < float-float */
bool lt(float49 a, float49 b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] < b.x[1]));
}

/* float < float-float */
bool lt(float a, float49 b) {
  return (a < b.x[0] || (a == b.x[0] && b.x[1] > 0.0));
}

/*********** Greater-Than-Or-Equal-To Comparisons ************/
/* float-float >= float */
bool ge(float49 a, float b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] >= 0.0));
}

/* float-float >= float-float */
bool ge(float49 a, float49 b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] >= b.x[1]));
}

/*********** Less-Than-Or-Equal-To Comparisons ************/
/* float-float <= float */
bool le(float49 a, float b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] <= 0.0));
}

/* float-float <= float-float */
bool le(float49 a, float49 b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] <= b.x[1]));
}

/* float <= float-float */
bool le(float a, float49 b) {
  return ge(b, a);
}

/* float >= float-float */
bool ge(float a, float49 b) {
  return le(b, a);
}

/*********** Not-Equal-To Comparisons ************/
/* float-float != float */
bool ne(float49 a, float b) {
  return (a.x[0] != b || a.x[1] != 0.0);
}

/* float-float != float-float */
bool ne(float49 a, float49 b) {
  return (a.x[0] != b.x[0] || a.x[1] != b.x[1]);
}

/* float != float-float */
bool ne(float a, float49 b) {
  return (a != b.x[0] || b.x[1] != 0.0);
}


float49 max(float a, float49 b)
{
  if (gt(a, b))
  {
    return float49_(a);
  }
  else
  {
    return b;
  }
}

float49 max(float49 b, float a)
{
  if (gt(a, b))
  {
    return float49_(a);
  }
  else
  {
    return b;
  }
}

float49 max(float49 a, float49 b)
{
  if (gt(a, b))
  {
    return a;
  }
  else
  {
    return b;
  }
}

/* Round to Nearest integer */
float49 nint(float49 a) {
  float hi = nint(a.x[0]);
  float lo;

  if (hi == a.x[0]) {
    /* High word is an integer already.  Round the low word.*/
    lo = nint(a.x[1]);

    /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
    hi = quick_two_sum(hi, lo, lo);
  } else {
    /* High word is not an integer. */
    lo = 0.0;
    if (abs(hi-a.x[0]) == 0.5 && a.x[1] < 0.0) {
      /* There is a tie in the high word, consult the low word
         to break the tie. */
      hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
    }
  }

  return float49_(hi, lo);
}

float49 floor(float49 a) {
  float hi = floor(a.x[0]);
  float lo = 0.0;

  if (hi == a.x[0]) {
    /* High word is integer already.  Round the low word. */
    lo = floor(a.x[1]);
    hi = quick_two_sum(hi, lo, lo);
  }

  return float49_(hi, lo);
}

float49 ceil(float49 a) {
  float hi = ceil(a.x[0]);
  float lo = 0.0;

  if (hi == a.x[0]) {
    /* High word is integer already.  Round the low word. */
    lo = ceil(a.x[1]);
    hi = quick_two_sum(hi, lo, lo);
  }

  return float49_(hi, lo);
}

float49 aint(float49 a) {
  return (a.x[0] >= 0.0) ? floor(a) : ceil(a);
}

/********** Remainder **********/
float49 drem(float49 a, float49 b) {
  float49 n = nint(div(a, b));
  return sub(a, mul(n, b));
}

float49 divrem(float49 a, float49 b, out float49 r) {
  float49 n = nint(div(a, b));
  r = sub(a, mul(n, b));
  return n;
}

/* Cast to float. */
float to_float(float49 a) {
  return a.x[0];
}

/* Cast to int. */
int to_int(float49 a) {
  return int(a.x[0]);
}

///=====================================================================
/// qd-2.3.22+dfsg.1/src/dd_real.cpp
///=====================================================================

/*
 * src/dd_real.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Contains implementation of non-inlined functions of float-float
 * package.  Inlined functions are found in dd_inline.h (in include directory).
 */

/* Computes the square root of the float-float number dd.
   NOTE: dd must be a non-negative number.                   */
float49 sqrt(float49 a) {
  /* Strategy:  Use Karp's trick:  if x is an approximation
     to sqrt(a), then

        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision.
  */

  if (is_zero(a))
    return f49_0;

  if (is_negative(a)) {
    return f49_nan;
  }

  float x = EXACT(1.0 / sqrt(a.x[0]));
  float ax = EXACT(a.x[0] * x);
  return f49_add(ax, sub(a, f49_sqr(ax)).x[0] * (x * 0.5));
}

/* Computes the square root of a float in float-float precision.
   NOTE: d must not be negative.                                   */
float49 f49_sqrt(float d) {
  return sqrt(float49_(d));
}

/* Computes the n-th root of the float-float number a.
   NOTE: n must be a positive integer.
   NOTE: If n is even, then a must not be negative.       */
float49 nroot(float49 a, int n) {
  /* Strategy:  Use Newton iteration for the function

          f(x) = x^(-n) - a

     to find its root a^{-1/n}.  The iteration is thus

          x' = x + x * (1 - a * x^n) / n

     which converges quadratically.  We can then find
    a^{1/n} by taking the reciprocal.
  */

  if (n <= 0) {
    return f49_nan;
  }

  if (n%2 == 0 && is_negative(a)) {
    return f49_nan;
  }

  if (n == 1) {
    return a;
  }
  if (n == 2) {
    return sqrt(a);
  }

  if (is_zero(a))
    return f49_0;

  /* Note  a^{-1/n} = exp(-log(a)/n) */
  float49 r = abs(a);
  float49 x; set(x, exp(-log(r.x[0]) / n));

  /* Perform Newton's iteration. */
  add_set(x, div(mul(x, sub(1.0, mul(r, npwr(x, n)))), float(n)));
  if (a.x[0] < 0.0)
    x = neg(x);
  return inv(x);
}

// mapM_ putStrLn [ "  , float49(float[2](" ++ show (ff (1 / factorial (2 + m))) ++ "))" | m <- [1..15] ]
const int n_inv_fact = 15;
const float49 inv_fact[n_inv_fact] = float49[n_inv_fact]
  ( float49(float[2](0.16666667,-4.967054e-9))
  , float49(float[2](4.1666668e-2,-1.2417635e-9))
  , float49(float[2](8.333334e-3,-4.346172e-10))
  , float49(float[2](1.3888889e-3,-3.3631094e-11))
  , float49(float[2](1.984127e-4,-2.7255969e-12))
  , float49(float[2](2.4801588e-5,-3.406996e-13))
  , float49(float[2](2.7557319e-6,3.7935712e-14))
  , float49(float[2](2.755732e-7,-7.575112e-15))
  , float49(float[2](2.5052108e-8,4.417623e-16))
  , float49(float[2](2.0876756e-9,1.108284e-16))
  , float49(float[2](1.6059044e-10,-5.3525265e-18))
  , float49(float[2](1.1470745e-11,2.3722077e-19))
  , float49(float[2](7.6471636e-13,1.22007105e-20))
  , float49(float[2](4.7794773e-14,7.625444e-22))
  , float49(float[2](2.8114574e-15,-1.0462085e-22))
  );

/* Exponential.  Computes exp(x) in float-float precision. */
float49 exp(float49 a) {
  /* Strategy:  We first reduce the size of x by noting that

          exp(kr + m * log(2)) = 2^m * exp(r)^k

     where m and k are integers.  By choosing m appropriately
     we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
     evaluated using the familiar Taylor series.  Reducing the
     argument substantially speeds up the convergence.       */

  const float k = 512.0;
  const float inv_k = 1.0 / k;

  if (a.x[0] <= -709.0)
    return f49_0;

  if (a.x[0] >=  709.0)
    return f49_inf;

  if (is_zero(a))
    return f49_1;

  if (is_one(a))
    return f49_e;

  float m = floor(a.x[0] / f49_log2.x[0] + 0.5);
  float49 r = mul_pwr2(sub(a, mul(f49_log2, m)), inv_k);
  float49 s, t, p;

  p = sqr(r);
  s = add(r, mul_pwr2(p, 0.5));
  mul_set(p, r);
  t = mul(p, inv_fact[0]);
  int i = 0;
  do {
    add_set(s, t);
    mul_set(p, r);
    ++i;
    t = mul(p, inv_fact[i]);
  } while (abs(to_float(t)) > inv_k * f49_eps && i < 5);

  add_set(s, t);

  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  s = add(mul_pwr2(s, 2.0), sqr(s));
  add_set(s, 1.0);

  return ldexp(s, int(m));
}

/* Logarithm.  Computes log(x) in float-float precision.
   This is a natural logarithm (i.e., base e).            */
float49 log(float49 a) {
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x)
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.

     Only one iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

  if (is_one(a)) {
    return f49_0;
  }

  if (a.x[0] <= 0.0) {
    return f49_nan;
  }

  float49 x;
  set(x, log(a.x[0]));   /* Initial approximation */

  x = sub(add(x, mul(a, exp(neg(x)))), 1.0);
  return x;
}

float49 log10(float49 a) {
  return div(log(a), f49_log10);
}

float49 pow(float49 a, float49 b) {
  return exp(mul(b, log(a)));
}

/* Table of sin(k * pi/16) and cos(k * pi/16). */
const float49 sin_table[4] = float49[4]
  ( float49(float[2](0.19509032,-1.6704715e-9))
  , float49(float[2](0.38268343,6.2233507e-9))
  , float49(float[2](0.55557024,-1.1769521e-8))
  , float49(float[2](0.70710677,1.21016175e-8))
  );
const float49 cos_table[4] = float49[4]
  ( float49(float[2](0.98078525,2.9739473e-8))
  , float49(float[2](0.9238795,2.830749e-8))
  , float49(float[2](0.8314696,1.6870263e-8))
  , float49(float[2](0.70710677,1.21016175e-8))
  );

/* Computes sin(a) using Taylor series.
   Assumes |a| <= pi/32.                           */
float49 sin_taylor(float49 a) {
  float thresh = 0.5 * abs(to_float(a)) * f49_eps;
  float49 r, s, t, x;

  if (is_zero(a)) {
    return f49_0;
  }

  int i = 0;
  x = neg(sqr(a));
  s = a;
  r = a;
  do {
    mul_set(r, x);
    t = mul(r, inv_fact[i]);
    add_set(s, t);
    i += 2;
  } while (i < n_inv_fact && abs(to_float(t)) > thresh);

  return s;
}

float49 cos_taylor(float49 a) {
  const float thresh = 0.5 * f49_eps;
  float49 r, s, t, x;

  if (is_zero(a)) {
    return f49_1;
  }

  x = neg(sqr(a));
  r = x;
  s = add(1.0, mul_pwr2(r, 0.5));
  int i = 1;
  do {
    mul_set(r, x);
    t = mul(r, inv_fact[i]);
    add_set(s, t);
    i += 2;
  } while (i < n_inv_fact && abs(to_float(t)) > thresh);

  return s;
}

void sincos_taylor(float49 a, out float49 sin_a, out float49 cos_a) {
  if (is_zero(a)) {
    set(sin_a, 0.0);
    set(cos_a, 1.0);
    return;
  }

  sin_a = sin_taylor(a);
  cos_a = sqrt(sub(1.0, sqr(sin_a)));
}


float49 sin(float49 a) {

  /* Strategy.  To compute sin(x), we choose integers a, b so that

       x = s + a * (pi/2) + b * (pi/16)

     and |s| <= pi/32.  Using the fact that

       sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))

     we can compute sin(x) from sin(s), cos(s).  This greatly
     increases the convergence of the sine Taylor series. */

  if (is_zero(a)) {
    return f49_0;
  }

  // approximately reduce modulo 2*pi
  float49 z = nint(div(a, f49_2pi));
  float49 r = sub(a, mul(f49_2pi, z));

  // approximately reduce modulo pi/2 and then modulo pi/16.
  float49 t;
  float q = floor(r.x[0] / f49_pi2.x[0] + 0.5);
  t = sub(r, mul(f49_pi2, q));
  int j = int(q);
  q = floor(t.x[0] / f49_pi16.x[0] + 0.5);
  sub_set(t, mul(f49_pi16, q));
  int k = int(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    return f49_nan;
  }

  if (abs_k > 4) {
    return f49_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return sin_taylor(t);
      case 1:
        return cos_taylor(t);
      case -1:
        return neg(cos_taylor(t));
      default:
        return neg(sin_taylor(t));
    }
  }

  float49 u = cos_table[abs_k-1];
  float49 v = sin_table[abs_k-1];
  float49 sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  if (j == 0) {
    if (k > 0) {
      r = add(mul(u, sin_t), mul(v, cos_t));
    } else {
      r = sub(mul(u, sin_t), mul(v, cos_t));
    }
  } else if (j == 1) {
    if (k > 0) {
      r = sub(mul(u, cos_t), mul(v, sin_t));
    } else {
      r = add(mul(u, cos_t), mul(v, sin_t));
    }
  } else if (j == -1) {
    if (k > 0) {
      r = sub(mul(v, sin_t), mul(u, cos_t));
    } else if (k < 0) {
      r = neg(add(mul(u, cos_t), mul(v, sin_t)));
    }
  } else {
    if (k > 0) {
      r = neg(add(mul(u, sin_t), mul(v, cos_t)));
    } else {
      r = sub(mul(v, cos_t), mul(u, sin_t));
    }
  }

  return r;
}

float49 cos(float49 a) {

  if (is_zero(a)) {
    return f49_1;
  }

  // approximately reduce modulo 2*pi
  float49 z = nint(div(a, f49_2pi));
  float49 r = sub(a, mul(z, f49_2pi));

  // approximately reduce modulo pi/2 and then modulo pi/16
  float49 t;
  float q = floor(r.x[0] / f49_pi2.x[0] + 0.5);
  t = sub(r, mul(f49_pi2, q));
  int j = int(q);
  q = floor(t.x[0] / f49_pi16.x[0] + 0.5);
  sub_set(t, mul(f49_pi16, q));
  int k = int(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    return f49_nan;
  }

  if (abs_k > 4) {
    return f49_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return cos_taylor(t);
      case 1:
        return neg(sin_taylor(t));
      case -1:
        return sin_taylor(t);
      default:
        return neg(cos_taylor(t));
    }
  }

  float49 sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  float49 u = cos_table[abs_k-1];
  float49 v = sin_table[abs_k-1];

  if (j == 0) {
    if (k > 0) {
      r = sub(mul(u, cos_t), mul(v, sin_t));
    } else {
      r = add(mul(u, cos_t), mul(v, sin_t));
    }
  } else if (j == 1) {
    if (k > 0) {
      r = neg(add(mul(u, sin_t), mul(v, cos_t)));
    } else {
      r = sub(mul(v, cos_t), mul(u, sin_t));
    }
  } else if (j == -1) {
    if (k > 0) {
      r = add(mul(u, sin_t), mul(v, cos_t));
    } else {
      r = sub(mul(u, sin_t), mul(v, cos_t));
    }
  } else {
    if (k > 0) {
      r = sub(mul(v, sin_t), mul(u, cos_t));
    } else {
      r = neg(add(mul(u, cos_t), mul(v, sin_t)));
    }
  }

  return r;
}

void sincos(float49 a, out float49 sin_a, out float49 cos_a) {

  if (is_zero(a)) {
    set(sin_a, 0.0);
    set(cos_a, 1.0);
    return;
  }

  // approximately reduce modulo 2*pi
  float49 z = nint(div(a, f49_2pi));
  float49 r = sub(a, mul(f49_2pi, z));

  // approximately reduce module pi/2 and pi/16
  float49 t;
  float q = floor(r.x[0] / f49_pi2.x[0] + 0.5);
  t = sub(r, mul(f49_pi2, q));
  int j = int(q);
  int abs_j = abs(j);
  q = floor(t.x[0] / f49_pi16.x[0] + 0.5);
  sub_set(t, mul(f49_pi16, q));
  int k = int(q);
  int abs_k = abs(k);

  if (abs_j > 2) {
    cos_a = sin_a = f49_nan;
    return;
  }

  if (abs_k > 4) {
    cos_a = sin_a = f49_nan;
    return;
  }

  float49 sin_t, cos_t;
  float49 s, c;

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    float49 u = cos_table[abs_k-1];
    float49 v = sin_table[abs_k-1];

    if (k > 0) {
      s = add(mul(u, sin_t), mul(v, cos_t));
      c = sub(mul(u, cos_t), mul(v, sin_t));
    } else {
      s = sub(mul(u, sin_t), mul(v, cos_t));
      c = add(mul(u, cos_t), mul(v, sin_t));
    }
  }

  if (abs_j == 0) {
    sin_a = s;
    cos_a = c;
  } else if (j == 1) {
    sin_a = c;
    cos_a = neg(s);
  } else if (j == -1) {
    sin_a = neg(c);
    cos_a = s;
  } else {
    sin_a = neg(s);
    cos_a = neg(c);
  }

}

float49 atan(float49 y, float49 x) {
  /* Strategy: Instead of using Taylor series to compute
     arctan, we instead use Newton's iteration to solve
     the equation

        sin(z) = y/r    or    cos(z) = x/r

     where r = sqrt(x^2 + y^2).
     The iteration is given by

        z' = z + (y - sin(z)) / cos(z)          (for equation 1)
        z' = z - (x - cos(z)) / sin(z)          (for equation 2)

     Here, x and y are normalized so that x^2 + y^2 = 1.
     If |x| > |y|, then first iteration is used since the
     denominator is larger.  Otherwise, the second is used.
  */

  if (is_zero(x)) {

    if (is_zero(y)) {
      /* Both x and y is zero. */
      return f49_nan;
    }

    return (is_positive(y)) ? f49_pi2 : neg(f49_pi2);
  } else if (is_zero(y)) {
    return (is_positive(x)) ? f49_0 : f49_pi;
  }

  if (eq(x, y)) {
    return (is_positive(y)) ? f49_pi4 : neg(f49_3pi4);
  }

  if (eq(x, neg(y))) {
    return (is_positive(y)) ? f49_3pi4 : neg(f49_pi4);
  }

  float49 r = sqrt(add(sqr(x), sqr(y)));
  float49 xx = div(x, r);
  float49 yy = div(y, r);

  /* Compute float precision approximation to atan. */
  float49 z; set(z, atan(to_float(y), to_float(x)));
  float49 sin_z, cos_z;

  if (abs(xx.x[0]) > abs(yy.x[0])) {
    /* Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)  */
    sincos(z, sin_z, cos_z);
    add_set(z, div(sub(yy, sin_z), cos_z));
  } else {
    /* Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)  */
    sincos(z, sin_z, cos_z);
    sub_set(z, div(sub(xx, cos_z), sin_z));
  }

  return z;
}

float49 atan(float49 a) {
  return atan(a, f49_1);
}

float49 tan(float49 a) {
  float49 s, c;
  sincos(a, s, c);
  return div(s, c);
}

float49 asin(float49 a) {
  float49 abs_a = abs(a);

  if (gt(abs_a, 1.0)) {
    return f49_nan;
  }

  if (is_one(abs_a)) {
    return (is_positive(a)) ? f49_pi2 : neg(f49_pi2);
  }

  return atan(a, sqrt(sub(1.0, sqr(a))));
}

float49 acos(float49 a) {
  float49 abs_a = abs(a);

  if (gt(abs_a, 1.0)) {
    return f49_nan;
  }

  if (is_one(abs_a)) {
    return (is_positive(a)) ? f49_0 : f49_pi;
  }

  return atan(sqrt(sub(1.0, sqr(a))), a);
}

float49 sinh(float49 a) {
  if (is_zero(a)) {
    return f49_0;
  }

  if (gt(abs(a), 0.05)) {
    float49 ea = exp(a);
    return mul_pwr2(sub(ea, inv(ea)), 0.5);
  }

  /* since a is small, using the above formula gives
     a lot of cancellation.  So use Taylor series.   */
  float49 s = a;
  float49 t = a;
  float49 r = sqr(t);
  float m = 1.0;
  float thresh = abs((to_float(a)) * f49_eps);

  do {
    m += 2.0;
    mul_set(t, r);
    div_set(t, (m-1) * m);

    add_set(s, t);
  } while (gt(abs(t), thresh));

  return s;

}

float49 cosh(float49 a) {
  if (is_zero(a)) {
    return f49_1;
  }

  float49 ea = exp(a);
  return mul_pwr2(add(ea, inv(ea)), 0.5);
}

float49 tanh(float49 a) {
  if (is_zero(a)) {
    return f49_0;
  }

  if (abs(to_float(a)) > 0.05) {
    float49 ea = exp(a);
    float49 inv_ea = inv(ea);
    return div(sub(ea, inv_ea), add(ea, inv_ea));
  } else {
    float49 s, c;
    s = sinh(a);
    c = sqrt(add(1.0, sqr(s)));
    return div(s, c);
  }
}

void sincosh(float49 a, out float49 s, out float49 c) {
  if (abs(to_float(a)) <= 0.05) {
    s = sinh(a);
    c = sqrt(add(1.0, sqr(s)));
  } else {
    float49 ea = exp(a);
    float49 inv_ea = inv(ea);
    s = mul_pwr2(sub(ea, inv_ea), 0.5);
    c = mul_pwr2(add(ea, inv_ea), 0.5);
  }
}

float49 asinh(float49 a) {
  return log(add(a, sqrt(add(sqr(a), 1.0))));
}

float49 acosh(float49 a) {
  if (lt(a, 1.0)) {
    return f49_nan;
  }

  return log(add(a, sqrt(sub(sqr(a), 1.0))));
}

float49 atanh(float49 a) {
  if (ge(abs(a), 1.0)) {
    return f49_nan;
  }

  return mul_pwr2(log(div(add(1.0, a), sub(1.0, a))), 0.5);
}

float49 fmod(float49 a, float49 b) {
  float49 n = aint(div(a, b));
  return sub(a, mul(b, n));
}

#if 0

///=====================================================================
/// qd-2.3.22+dfsg.1/src/qd_const.h
///=====================================================================

/*
 * src/qd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines constants used in quad-float package.
 */

#if 0
/* Some useful constants. */
const qd_real qd_2pi = qd_real(6.283185307179586232e+00,
                                      2.449293598294706414e-16,
                                      -5.989539619436679332e-33,
                                      2.224908441726730563e-49);
const qd_real qd_pi = qd_real(3.141592653589793116e+00,
                                     1.224646799147353207e-16,
                                     -2.994769809718339666e-33,
                                     1.112454220863365282e-49);
const qd_real qd_pi2 = qd_real(1.570796326794896558e+00,
                                      6.123233995736766036e-17,
                                      -1.497384904859169833e-33,
                                      5.562271104316826408e-50);
const qd_real qd_pi4 = qd_real(7.853981633974482790e-01,
                                      3.061616997868383018e-17,
                                      -7.486924524295849165e-34,
                                      2.781135552158413204e-50);
const qd_real qd_3pi4 = qd_real(2.356194490192344837e+00,
                                       9.1848509936051484375e-17,
                                       3.9168984647504003225e-33,
                                      -2.5867981632704860386e-49);
const qd_real qd_e = qd_real(2.718281828459045091e+00,
                                    1.445646891729250158e-16,
                                    -2.127717108038176765e-33,
                                    1.515630159841218954e-49);
const qd_real qd_log2 = qd_real(6.931471805599452862e-01,
                                       2.319046813846299558e-17,
                                       5.707708438416212066e-34,
                                       -3.582432210601811423e-50);
const qd_real qd_log10 = qd_real(2.302585092994045901e+00,
                                        -2.170756223382249351e-16,
                                        -9.984262454465776570e-33,
                                        -4.023357454450206379e-49);
const qd_real qd_nan = qd_real(_d_nan, _d_nan, _d_nan, _d_nan);
const qd_real qd_inf = qd_real(_d_inf, _d_inf, _d_inf, _d_inf);

const float qd_eps = 1.21543267145725e-63; // = 2^-209
const float qd_min_normalized = 1.6259745436952323e-260; // = 2^(-1022 + 3*53)
const qd_real qd_max = qd_real(
    1.79769313486231570815e+308, 9.97920154767359795037e+291,
    5.53956966280111259858e+275, 3.07507889307840487279e+259);
const qd_real qd_safe_max = qd_real(
    1.7976931080746007281e+308,  9.97920154767359795037e+291,
    5.53956966280111259858e+275, 3.07507889307840487279e+259);
const int qd_ndigits = 62;
#endif

///=====================================================================
/// qd-2.3.22+dfsg.1/include/qd/qd_inline.h
///=====================================================================

/*
 * include/qd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the quad-float
 * arithmetic package.
 */

/********** Constructors **********/
qd_real qd_real(float x0, float x1, float x2, float x3) {
  qd_real r;
  r.x[0] = x0;
  r.x[1] = x1;
  r.x[2] = x2;
  r.x[3] = x3;
  return r;
}

qd_real qd_real(const float xx[4]) {
  qd_real r;
  r.x[0] = xx[0];
  r.x[1] = xx[1];
  r.x[2] = xx[2];
  r.x[3] = xx[3];
  return r;
}

qd_real qd_real(float x0) {
  qd_real r;
  r.x[0] = x0;
  r.x[1] = r.x[2] = r.x[3] = 0.0;
  return r;
}

qd_real qd_real() {
  qd_real r;
  r.x[0] = 0.0;
  r.x[1] = 0.0;
  r.x[2] = 0.0;
  r.x[3] = 0.0;
  return r;
}

qd_real qd_real(float49 a) {
  qd_real r;
  r.x[0] = a.x[0];
  r.x[1] = a.x[1];
  r.x[2] = r.x[3] = 0.0;
  return r;
}

qd_real qd_real(int i) {
  qd_real r;
  r.x[0] = float(i);
  r.x[1] = r.x[2] = r.x[3] = 0.0;
  return r;
}

/********** Accessors **********/
bool isnan(qd_real x) const {
  return QD_ISNAN(x.x[0]) || QD_ISNAN(x.x[1]) || QD_ISNAN(x.x[2]) || QD_ISNAN(x.x[3]);
}

/********** Renormalization **********/
void quick_renorm(inout float c0, inout float c1,
                         inout float c2, inout float c3, inout float c4) {
  float t0, t1, t2, t3;
  float s;
  s  = quick_two_sum(c3, c4, t3);
  s  = quick_two_sum(c2, s , t2);
  s  = quick_two_sum(c1, s , t1);
  c0 = quick_two_sum(c0, s , t0);

  s  = quick_two_sum(t2, t3, t2);
  s  = quick_two_sum(t1, s , t1);
  c1 = quick_two_sum(t0, s , t0);

  s  = quick_two_sum(t1, t2, t1);
  c2 = quick_two_sum(t0, s , t0);

  c3 = t0 + t1;
}

void renorm(inout float c0, inout float c1,
                   inout float c2, inout float c3) {
  float s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = quick_two_sum(c2, c3, c3);
  s0 = quick_two_sum(c1, s0, c2);
  c0 = quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != 0.0) {
    s1 = quick_two_sum(s1, c2, s2);
    if (s2 != 0.0)
      s2 = quick_two_sum(s2, c3, s3);
    else
      s1 = quick_two_sum(s1, c3, s2);
  } else {
    s0 = quick_two_sum(s0, c2, s1);
    if (s1 != 0.0)
      s1 = quick_two_sum(s1, c3, s2);
    else
      s0 = quick_two_sum(s0, c3, s1);
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}

void renorm(inout float c0, inout float c1,
                   inout float c2, inout float c3, inout float c4) {
  float s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = quick_two_sum(c3, c4, c4);
  s0 = quick_two_sum(c2, s0, c3);
  s0 = quick_two_sum(c1, s0, c2);
  c0 = quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;

  if (s1 != 0.0) {
    s1 = quick_two_sum(s1, c2, s2);
    if (s2 != 0.0) {
      s2 = quick_two_sum(s2, c3, s3);
      if (s3 != 0.0)
        s3 += c4;
      else
        s2 = quick_two_sum(s2, c4, s3);
    } else {
      s1 = quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = quick_two_sum(s2, c4, s3);
      else
        s1 = quick_two_sum(s1, c4, s2);
    }
  } else {
    s0 = quick_two_sum(s0, c2, s1);
    if (s1 != 0.0) {
      s1 = quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = quick_two_sum(s2, c4, s3);
      else
        s1 = quick_two_sum(s1, c4, s2);
    } else {
      s0 = quick_two_sum(s0, c3, s1);
      if (s1 != 0.0)
        s1 = quick_two_sum(s1, c4, s2);
      else
        s0 = quick_two_sum(s0, c4, s1);
    }
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}

void renorm(inout qd_real self) {
  renorm(self.x[0], self.x[1], self.x[2], self.x[3]);
}

void renorm(inout qd_real self, inout float e) {
  renorm(self.x[0], self.x[1], self.x[2], self.x[3], e);
}


/********** Additions ************/
void three_sum(inout float a, inout float b, inout float c) {
  float t1, t2, t3;
  t1 = two_sum(a, b, t2);
  a  = two_sum(c, t1, t3);
  b  = two_sum(t2, t3, c);
}

void three_sum2(inout float a, inout float b, inout float c) {
  float t1, t2, t3;
  t1 = two_sum(a, b, t2);
  a  = two_sum(c, t1, t3);
  b = t2 + t3;
}

/* quad-float + float */
qd_real add(qd_real a, float b) {
  float c0, c1, c2, c3;
  float e;

  c0 = two_sum(a[0], b, e);
  c1 = two_sum(a[1], e, e);
  c2 = two_sum(a[2], e, e);
  c3 = two_sum(a[3], e, e);

  renorm(c0, c1, c2, c3, e);

  return qd_real(c0, c1, c2, c3);
}

/* quad-float + float-float */
qd_real add(qd_real a, float49 b) {

  float s0, s1, s2, s3;
  float t0, t1;

  s0 = two_sum(a[0], b.x[0], t0);
  s1 = two_sum(a[1], b.x[1], t1);

  s1 = two_sum(s1, t0, t0);

  s2 = a[2];
  three_sum(s2, t0, t1);

  s3 = two_sum(t0, a[3], t0);
  t0 += t1;

  renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3);
}


/* float + quad-float */
qd_real add(float a, qd_real b) {
  return (b + a);
}

/* float-float + quad-float */
qd_real add(float49 a, qd_real b) {
  return (b + a);
}

/* s = quick_three_accum(a, b, c) adds c to the dd-pair (a, b).
 * If the result does not fit in two floats, then the sum is
 * output into s and (a,b) contains the remainder.  Otherwise
 * s is zero and (a,b) contains the sum. */
float quick_three_accum(inout float a, inout float b, float c) {
  float s;
  bool za, zb;

  s = two_sum(b, c, b);
  s = two_sum(a, s, a);

  za = (a != 0.0);
  zb = (b != 0.0);

  if (za && zb)
    return s;

  if (!zb) {
    b = a;
    a = s;
  } else {
    a = s;
  }

  return 0.0;
}

qd_real ieee_add(qd_real a, qd_real b) {
  int i, j, k;
  float s, t;
  float u, v;   /* float-length accumulator */
  float x[4] = {0.0, 0.0, 0.0, 0.0};

  i = j = k = 0;
  if (abs(a[i]) > abs(b[j]))
    u = a[i++];
  else
    u = b[j++];
  if (abs(a[i]) > abs(b[j]))
    v = a[i++];
  else
    v = b[j++];

  u = quick_two_sum(u, v, v);

  while (k < 4) {
    if (i >= 4 && j >= 4) {
      x[k] = u;
      if (k < 3)
        x[++k] = v;
      break;
    }

    if (i >= 4)
      t = b[j++];
    else if (j >= 4)
      t = a[i++];
    else if (abs(a[i]) > abs(b[j])) {
      t = a[i++];
    } else
      t = b[j++];

    s = quick_three_accum(u, v, t);

    if (s != 0.0) {
      x[k++] = s;
    }
  }

  /* add the rest. */
  for (k = i; k < 4; k++)
    x[3] += a[k];
  for (k = j; k < 4; k++)
    x[3] += b[k];

  renorm(x[0], x[1], x[2], x[3]);
  return qd_real(x[0], x[1], x[2], x[3]);
}

qd_real sloppy_add(qd_real a, qd_real b) {
  /*
  float s0, s1, s2, s3;
  float t0, t1, t2, t3;

  s0 = two_sum(a[0], b[0], t0);
  s1 = two_sum(a[1], b[1], t1);
  s2 = two_sum(a[2], b[2], t2);
  s3 = two_sum(a[3], b[3], t3);

  s1 = two_sum(s1, t0, t0);
  three_sum(s2, t0, t1);
  three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;

  renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3, t0);
  */

  /* Same as above, but addition re-organized to minimize
     data dependency ... unfortunately some compilers are
     not very smart to do this automatically */
  float s0, s1, s2, s3;
  float t0, t1, t2, t3;

  float v0, v1, v2, v3;
  float u0, u1, u2, u3;
  float w0, w1, w2, w3;

  s0 = a[0] + b[0];
  s1 = a[1] + b[1];
  s2 = a[2] + b[2];
  s3 = a[3] + b[3];

  v0 = s0 - a[0];
  v1 = s1 - a[1];
  v2 = s2 - a[2];
  v3 = s3 - a[3];

  u0 = s0 - v0;
  u1 = s1 - v1;
  u2 = s2 - v2;
  u3 = s3 - v3;

  w0 = a[0] - u0;
  w1 = a[1] - u1;
  w2 = a[2] - u2;
  w3 = a[3] - u3;

  u0 = b[0] - v0;
  u1 = b[1] - v1;
  u2 = b[2] - v2;
  u3 = b[3] - v3;

  t0 = w0 + u0;
  t1 = w1 + u1;
  t2 = w2 + u2;
  t3 = w3 + u3;

  s1 = two_sum(s1, t0, t0);
  three_sum(s2, t0, t1);
  three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;

  /* renormalize */
  renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3);
}

/* quad-float + quad-float */
qd_real add(qd_real a, qd_real b) {
#ifndef QD_IEEE_ADD
  return sloppy_add(a, b);
#else
  return ieee_add(a, b);
#endif
}



/********** Self-Additions ************/
/* quad-float += float */
void add_set(inout qd_real self, float a) {
  self = self + a;
}

/* quad-float += float-float */
void add_set(inout qd_real self, float49 a) {
  self = self + a;
}

/* quad-float += quad-float */
void add_set(inout qd_real self, qd_real a) {
  self = self + a;
}

/********** Unary Minus **********/
qd_real neg(qd_real a) {
  return qd_real(-a.x[0], -a.x[1], -a.x[2], -a.x[3]);
}

/********** Subtractions **********/
qd_real sub(qd_real a, float b) {
  return (a + (-b));
}

qd_real sub(float a, qd_real b) {
  return (a + (-b));
}

qd_real sub(qd_real a, float49 b) {
  return (a + (-b));
}

qd_real sub(float49 a, qd_real b) {
  return (a + (-b));
}

qd_real sub(qd_real a, qd_real b) {
  return (a + (-b));
}

/********** Self-Subtractions **********/
void sub_set(inout qd_real self, float a) {
  return ((self) += (-a));
}

void sub_set(inout qd_real self, float49 a) {
  return ((self) += (-a));
}

void sub_set(inout qd_real self, qd_real a) {
  return ((self) += (-a));
}


qd_real mul(float a, qd_real b) {
  return (b * a);
}

qd_real mul(float49 a, qd_real b) {
  return (b * a);
}

qd_real mul_pwr2(qd_real a, float b) {
  return qd_real(a[0] * b, a[1] * b, a[2] * b, a[3] * b);
}

/********** Multiplications **********/
qd_real mul(qd_real a, float b) {
  float p0, p1, p2, p3;
  float q0, q1, q2;
  float s0, s1, s2, s3, s4;

  p0 = two_prod(a[0], b, q0);
  p1 = two_prod(a[1], b, q1);
  p2 = two_prod(a[2], b, q2);
  p3 = a[3] * b;

  s0 = p0;

  s1 = two_sum(q0, p1, s2);

  three_sum(s2, q1, p2);

  three_sum2(q1, q2, p3);
  s3 = q1;

  s4 = q2 + p2;

  renorm(s0, s1, s2, s3, s4);
  return qd_real(s0, s1, s2, s3);

}

/* quad-float * float-float */
/* a0 * b0                        0
        a0 * b1                   1
        a1 * b0                   2
             a1 * b1              3
             a2 * b0              4
                  a2 * b1         5
                  a3 * b0         6
                       a3 * b1    7 */
qd_real mul(qd_real a, float49 b) {
  float p0, p1, p2, p3, p4;
  float q0, q1, q2, q3, q4;
  float s0, s1, s2;
  float t0, t1;

  p0 = two_prod(a[0], b.x[0], q0);
  p1 = two_prod(a[0], b.x[1], q1);
  p2 = two_prod(a[1], b.x[0], q2);
  p3 = two_prod(a[1], b.x[1], q3);
  p4 = two_prod(a[2], b.x[0], q4);

  three_sum(p1, p2, q0);

  /* Five-Three-Sum */
  three_sum(p2, p3, p4);
  q1 = two_sum(q1, q2, q2);
  s0 = two_sum(p2, q1, t0);
  s1 = two_sum(p3, q2, t1);
  s1 = two_sum(s1, t0, t0);
  s2 = t0 + t1 + p4;
  p2 = s0;

  p3 = a[2] * b.x[0] + a[3] * b.x[1] + q3 + q4;
  three_sum2(p3, q0, s1);
  p4 = q0 + s2;

  renorm(p0, p1, p2, p3, p4);
  return qd_real(p0, p1, p2, p3);
}

/* quad-float * quad-float */
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a0 * b2          3
             a1 * b1          4
             a2 * b0          5
                  a0 * b3     6
                  a1 * b2     7
                  a2 * b1     8
                  a3 * b0     9  */
qd_real sloppy_mul(qd_real a, qd_real b) {
  float p0, p1, p2, p3, p4, p5;
  float q0, q1, q2, q3, q4, q5;
  float t0, t1;
  float s0, s1, s2;

  p0 = two_prod(a[0], b[0], q0);

  p1 = two_prod(a[0], b[1], q1);
  p2 = two_prod(a[1], b[0], q2);

  p3 = two_prod(a[0], b[2], q3);
  p4 = two_prod(a[1], b[1], q4);
  p5 = two_prod(a[2], b[0], q5);

  /* Start Accumulation */
  three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  three_sum(p2, q1, q2);
  three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  s0 = two_sum(p2, p3, t0);
  s1 = two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  s1 += a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0] + q0 + q3 + q4 + q5;
  renorm(p0, p1, s0, s1, s2);
  return qd_real(p0, p1, s0, s1);
}

qd_real accurate_mul(qd_real a, qd_real b) {
  float p0, p1, p2, p3, p4, p5;
  float q0, q1, q2, q3, q4, q5;
  float p6, p7, p8, p9;
  float q6, q7, q8, q9;
  float r0, r1;
  float t0, t1;
  float s0, s1, s2;

  p0 = two_prod(a[0], b[0], q0);

  p1 = two_prod(a[0], b[1], q1);
  p2 = two_prod(a[1], b[0], q2);

  p3 = two_prod(a[0], b[2], q3);
  p4 = two_prod(a[1], b[1], q4);
  p5 = two_prod(a[2], b[0], q5);

  /* Start Accumulation */
  three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  three_sum(p2, q1, q2);
  three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  s0 = two_sum(p2, p3, t0);
  s1 = two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  p6 = two_prod(a[0], b[3], q6);
  p7 = two_prod(a[1], b[2], q7);
  p8 = two_prod(a[2], b[1], q8);
  p9 = two_prod(a[3], b[0], q9);

  /* Nine-Two-Sum of q0, s1, q3, q4, q5, p6, p7, p8, p9. */
  q0 = two_sum(q0, q3, q3);
  q4 = two_sum(q4, q5, q5);
  p6 = two_sum(p6, p7, p7);
  p8 = two_sum(p8, p9, p9);
  /* Compute (t0, t1) = (q0, q3) + (q4, q5). */
  t0 = two_sum(q0, q4, t1);
  t1 += (q3 + q5);
  /* Compute (r0, r1) = (p6, p7) + (p8, p9). */
  r0 = two_sum(p6, p8, r1);
  r1 += (p7 + p9);
  /* Compute (q3, q4) = (t0, t1) + (r0, r1). */
  q3 = two_sum(t0, r0, q4);
  q4 += (t1 + r1);
  /* Compute (t0, t1) = (q3, q4) + s1. */
  t0 = two_sum(q3, s1, t1);
  t1 += q4;

  /* O(eps^4) terms -- Nine-One-Sum */
  t1 += a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + q6 + q7 + q8 + q9 + s2;

  renorm(p0, p1, s0, t0, t1);
  return qd_real(p0, p1, s0, t0);
}

qd_real mul(qd_real a, qd_real b) {
#ifdef QD_SLOPPY_MUL
  return sloppy_mul(a, b);
#else
  return accurate_mul(a, b);
#endif
}

/* quad-float ^ 2  = (x0 + x1 + x2 + x3) ^ 2
                    = x0 ^ 2 + 2 x0 * x1 + (2 x0 * x2 + x1 ^ 2)
                               + (2 x0 * x3 + 2 x1 * x2)           */
qd_real sqr(qd_real a) {
  float p0, p1, p2, p3, p4, p5;
  float q0, q1, q2, q3;
  float s0, s1;
  float t0, t1;

  p0 = two_sqr(a[0], q0);
  p1 = two_prod(2.0 * a[0], a[1], q1);
  p2 = two_prod(2.0 * a[0], a[2], q2);
  p3 = two_sqr(a[1], q3);

  p1 = two_sum(q0, p1, q0);

  q0 = two_sum(q0, q1, q1);
  p2 = two_sum(p2, p3, p3);

  s0 = two_sum(q0, p2, t0);
  s1 = two_sum(q1, p3, t1);

  s1 = two_sum(s1, t0, t0);
  t0 += t1;

  s1 = quick_two_sum(s1, t0, t0);
  p2 = quick_two_sum(s0, s1, t1);
  p3 = quick_two_sum(t1, t0, q0);

  p4 = 2.0 * a[0] * a[3];
  p5 = 2.0 * a[1] * a[2];

  p4 = two_sum(p4, p5, p5);
  q2 = two_sum(q2, q3, q3);

  t0 = two_sum(p4, q2, t1);
  t1 = t1 + p5 + q3;

  p3 = two_sum(p3, t0, p4);
  p4 = p4 + q0 + t1;

  renorm(p0, p1, p2, p3, p4);
  return qd_real(p0, p1, p2, p3);

}

/********** Self-Multiplication **********/
/* quad-float *= float */
void mul_set(inout qd_real self, float a) {
  self = (self * a);
}

/* quad-float *= float-float */
void mul_set(inout qd_real self, float49 a) {
  self = (self * a);
}

/* quad-float *= quad-float */
void mul_set(inout qd_real self, qd_real a) {
  self = self * a;
}

qd_real div(qd_real a, float49 b) {
#ifdef QD_SLOPPY_DIV
  return sloppy_div(a, b);
#else
  return accurate_div(a, b);
#endif
}

qd_real div(qd_real a, qd_real b) {
#ifdef QD_SLOPPY_DIV
  return sloppy_div(a, b);
#else
  return accurate_div(a, b);
#endif
}

/* float / quad-float */
qd_real div(float a, qd_real b) {
  return qd_real(a) / b;
}

/* float-float / quad-float */
qd_real div(float49 a, qd_real b) {
  return qd_real(a) / b;
}

/********** Self-Divisions **********/
/* quad-float /= float */
void div_set(inout qd_real self, float a) {
  self = (self / a);
}

/* quad-float /= float-float */
void div_set(inout qd_real self, float49 a) {
  self = (self / a);
}

/* quad-float /= quad-float */
void div_set(inout qd_real self, qd_real a) {
  self = (self / a);
}


/********** Miscellaneous **********/
qd_real abs(qd_real a) {
  return (a[0] < 0.0) ? -a : a;
}

qd_real fabs(qd_real a) {
  return abs(a);
}

/* Quick version.  May be off by one when qd is very close
   to the middle of two integers.                         */
qd_real quick_nint(qd_real a) {
  qd_real r = qd_real(nint(a[0]), nint(a[1]),
      nint(a[2]), nint(a[3]));
  r.renorm();
  return r;
}

/*********** Assignments ************/
/* quad-float = float */
void set(inout qd_real self, float a) {
  self.x[0] = a;
  self.x[1] = self.x[2] = self.x[3] = 0.0;
}

/* quad-float = float-float */
void set(inout qd_real self, float49 a) {
  self.x[0] = a.x[0];
  self.x[1] = a.x[1];
  self.x[2] = self.x[3] = 0.0;
}

/********** Equality Comparison **********/
bool eq(qd_real a, float b) {
  return (a[0] == b && a[1] == 0.0 && a[2] == 0.0 && a[3] == 0.0);
}

bool eq(float a, qd_real b) {
  return (b == a);
}

bool eq(qd_real a, float49 b) {
  return (a[0] == b.x[0] && a[1] == b.x[1] &&
          a[2] == 0.0 && a[3] == 0.0);
}

bool eq(float49 a, qd_real b) {
  return (b == a);
}

bool eq(qd_real a, qd_real b) {
  return (a[0] == b[0] && a[1] == b[1] &&
          a[2] == b[2] && a[3] == b[3]);
}


/********** Less-Than Comparison ***********/
bool lt(qd_real a, float b) {
  return (a[0] < b || (a[0] == b && a[1] < 0.0));
}

bool lt(float a, qd_real b) {
  return (b > a);
}

bool lt(qd_real a, float49 b) {
  return (a[0] < b.x[0] ||
          (a[0] == b.x[0] && (a[1] < b.x[1] ||
                            (a[1] == b.x[1] && a[2] < 0.0))));
}

bool lt(float49 a, qd_real b) {
  return (b > a);
}

bool lt(qd_real a, qd_real b) {
  return (a[0] < b[0] ||
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && (a[2] < b[2] ||
                                              (a[2] == b[2] && a[3] < b[3]))))));
}

/********** Greater-Than Comparison ***********/
bool gt(qd_real a, float b) {
  return (a[0] > b || (a[0] == b && a[1] > 0.0));
}

bool gt(float a, qd_real b) {
  return (b < a);
}

bool gt(qd_real a, float49 b) {
  return (a[0] > b.x[0] ||
          (a[0] == b.x[0] && (a[1] > b.x[1] ||
                            (a[1] == b.x[1] && a[2] > 0.0))));
}

bool gt(float49 a, qd_real b) {
  return (b < a);
}

bool gt(qd_real a, qd_real b) {
  return (a[0] > b[0] ||
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && (a[2] > b[2] ||
                                              (a[2] == b[2] && a[3] > b[3]))))));
}


/********** Less-Than-Or-Equal-To Comparison **********/
bool le(qd_real a, float b) {
  return (a[0] < b || (a[0] == b && a[1] <= 0.0));
}

bool le(float a, qd_real b) {
  return (b >= a);
}

bool le(qd_real a, float49 b) {
  return (a[0] < b.x[0] ||
          (a[0] == b.x[0] && (a[1] < b.x[1] ||
                            (a[1] == b.x[1] && a[2] <= 0.0))));
}

bool le(float49 a, qd_real b) {
  return (b >= a);
}

bool le(qd_real a, qd_real b) {
  return (a[0] < b[0] ||
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && (a[2] < b[2] ||
                                              (a[2] == b[2] && a[3] <= b[3]))))));
}

/********** Greater-Than-Or-Equal-To Comparison **********/
bool ge(qd_real a, float b) {
  return (a[0] > b || (a[0] == b && a[1] >= 0.0));
}

bool ge(float a, qd_real b) {
  return (b <= a);
}

bool ge(qd_real a, float49 b) {
  return (a[0] > b.x[0] ||
          (a[0] == b.x[0] && (a[1] > b.x[1] ||
                            (a[1] == b.x[1] && a[2] >= 0.0))));
}

bool ge(float49 a, qd_real b) {
  return (b <= a);
}

bool ge(qd_real a, qd_real b) {
  return (a[0] > b[0] ||
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && (a[2] > b[2] ||
                                              (a[2] == b[2] && a[3] >= b[3]))))));
}



/********** Not-Equal-To Comparison **********/
bool ne(qd_real a, float b) {
  return !(a == b);
}

bool ne(float a, qd_real b) {
  return !(a == b);
}

bool ne(qd_real a, float49 b) {
  return !(a == b);
}

bool ne(float49 a, qd_real b) {
  return !(a == b);
}

bool ne(qd_real a, qd_real b) {
  return !(a == b);
}



qd_real aint(qd_real a) {
  return (a[0] >= 0) ? floor(a) : ceil(a);
}

bool is_zero(qd_real self) {
  return (self.x[0] == 0.0);
}

bool is_one(qd_real self) {
  return (self.x[0] == 1.0 && self.x[1] == 0.0 && self.x[2] == 0.0 && self.x[3] == 0.0);
}

bool is_positive(qd_real self) {
  return (self.x[0] > 0.0);
}

bool is_negative(qd_real self) {
  return (self.x[0] < 0.0);
}

float49 to_float49_(qd_real a) {
  return float49_(a[0], a[1]);
}

float to_float(qd_real a) {
  return a[0];
}

int to_int(qd_real a) {
  return int(a[0]);
}

qd_real inv(qd_real qd) {
  return 1.0 / qd;
}

qd_real max(qd_real a, qd_real b) {
  return (a > b) ? a : b;
}

qd_real max(qd_real a, qd_real b,
                   qd_real c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

qd_real min(qd_real a, qd_real b) {
  return (a < b) ? a : b;
}

qd_real min(qd_real a, qd_real b,
                   qd_real c) {
  return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

qd_real ldexp(qd_real a, int n) {
  return qd_real(ldexp(a[0], n), ldexp(a[1], n),
                 ldexp(a[2], n), ldexp(a[3], n));
}

///=====================================================================
/// qd-2.3.22+dfsg.1/src/qd_real.h
///=====================================================================
/*
 * src/qd_real.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Contains implementation of non-inlined functions of quad-float
 * package.  Inlined functions are found in qd_inline.h (in include directory).
 */

/********** Multiplications **********/

qd_real nint(qd_real a) {
  float x0, x1, x2, x3;

  x0 = nint(a[0]);
  x1 = x2 = x3 = 0.0;

  if (x0 == a[0]) {
    /* First float is already an integer. */
    x1 = nint(a[1]);

    if (x1 == a[1]) {
      /* Second float is already an integer. */
      x2 = nint(a[2]);

      if (x2 == a[2]) {
        /* Third float is already an integer. */
        x3 = nint(a[3]);
      } else {
        if (abs(x2 - a[2]) == 0.5 && a[3] < 0.0) {
          x2 -= 1.0;
        }
      }

    } else {
      if (abs(x1 - a[1]) == 0.5 && a[2] < 0.0) {
          x1 -= 1.0;
      }
    }

  } else {
    /* First float is not an integer. */
      if (abs(x0 - a[0]) == 0.5 && a[1] < 0.0) {
          x0 -= 1.0;
      }
  }

  renorm(x0, x1, x2, x3);
  return qd_real(x0, x1, x2, x3);
}

qd_real floor(qd_real a) {
  float x0, x1, x2, x3;
  x1 = x2 = x3 = 0.0;
  x0 = floor(a[0]);

  if (x0 == a[0]) {
    x1 = floor(a[1]);

    if (x1 == a[1]) {
      x2 = floor(a[2]);

      if (x2 == a[2]) {
        x3 = floor(a[3]);
      }
    }

    renorm(x0, x1, x2, x3);
    return qd_real(x0, x1, x2, x3);
  }

  return qd_real(x0, x1, x2, x3);
}

qd_real ceil(qd_real a) {
  float x0, x1, x2, x3;
  x1 = x2 = x3 = 0.0;
  x0 = ceil(a[0]);

  if (x0 == a[0]) {
    x1 = ceil(a[1]);

    if (x1 == a[1]) {
      x2 = ceil(a[2]);

      if (x2 == a[2]) {
        x3 = ceil(a[3]);
      }
    }

    renorm(x0, x1, x2, x3);
    return qd_real(x0, x1, x2, x3);
  }

  return qd_real(x0, x1, x2, x3);
}



/********** Divisions **********/
/* quad-float / float */
qd_real div(qd_real a, float b) {
  /* Strategy:  compute approximate quotient using high order
     floats, and then correct it 3 times using the remainder.
     (Analogous to long division.)                             */
  float t0, t1;
  float q0, q1, q2, q3;
  qd_real r;

  q0 = a[0] / b;  /* approximate quotient */

  /* Compute the remainder  a - q0 * b */
  t0 = two_prod(q0, b, t1);
  r = a - float49_(t0, t1);

  /* Compute the first correction */
  q1 = r[0] / b;
  t0 = two_prod(q1, b, t1);
  r -= float49_(t0, t1);

  /* Second correction to the quotient. */
  q2 = r[0] / b;
  t0 = two_prod(q2, b, t1);
  r -= float49_(t0, t1);

  /* Final correction to the quotient. */
  q3 = r[0] / b;

  renorm(q0, q1, q2, q3);
  return qd_real(q0, q1, q2, q3);
}

/* Computes  qd^n, where n is an integer. */
qd_real pow(qd_real a, int n) {
  if (n == 0)
    return 1.0;

  qd_real r = a;   /* odd-case multiplier */
  qd_real s = 1.0;  /* current answer */
  int N = abs(n);

  if (N > 1) {

    /* Use binary exponentiation. */
    while (N > 0) {
      if (N % 2 == 1) {
        /* If odd, multiply by r. Note eventually N = 1, so this
         eventually executes. */
        s *= r;
      }
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }

  } else {
    s = r;
  }

  if (n < 0)
    return (1.0 / s);

  return s;
}

qd_real pow(qd_real a, qd_real b) {
  return exp(b * log(a));
}

qd_real npwr(qd_real a, int n) {
  return pow(a, n);
}

/* Divisions */
/* quad-float / float-float */
qd_real sloppy_div(qd_real a, float49 b) {
  float q0, q1, q2, q3;
  qd_real r;
  qd_real qd_b(b);

  q0 = a[0] / b.x[0];
  r = a - q0 * qd_b;

  q1 = r[0] / b.x[0];
  r -= (q1 * qd_b);

  q2 = r[0] / b.x[0];
  r -= (q2 * qd_b);

  q3 = r[0] / b.x[0];

  renorm(q0, q1, q2, q3);
  return qd_real(q0, q1, q2, q3);
}

qd_real accurate_div(qd_real a, float49 b) {
  float q0, q1, q2, q3, q4;
  qd_real r;
  qd_real qd_b(b);

  q0 = a[0] / b.x[0];
  r = a - q0 * qd_b;

  q1 = r[0] / b.x[0];
  r -= (q1 * qd_b);

  q2 = r[0] / b.x[0];
  r -= (q2 * qd_b);

  q3 = r[0] / b.x[0];
  r -= (q3 * qd_b);

  q4 = r[0] / b.x[0];

  renorm(q0, q1, q2, q3, q4);
  return qd_real(q0, q1, q2, q3);
}

/* quad-float / quad-float */
qd_real sloppy_div(qd_real a, qd_real b) {
  float q0, q1, q2, q3;

  qd_real r;

  q0 = a[0] / b[0];
  r = a - (b * q0);

  q1 = r[0] / b[0];
  r -= (b * q1);

  q2 = r[0] / b[0];
  r -= (b * q2);

  q3 = r[0] / b[0];

  renorm(q0, q1, q2, q3);

  return qd_real(q0, q1, q2, q3);
}

qd_real accurate_div(qd_real a, qd_real b) {
  float q0, q1, q2, q3;

  qd_real r;

  q0 = a[0] / b[0];
  r = a - (b * q0);

  q1 = r[0] / b[0];
  r -= (b * q1);

  q2 = r[0] / b[0];
  r -= (b * q2);

  q3 = r[0] / b[0];

  r -= (b * q3);
  float q4 = r[0] / b[0];

  renorm(q0, q1, q2, q3, q4);

  return qd_real(q0, q1, q2, q3);
}

qd_real sqrt(qd_real a) {
  /* Strategy:

     Perform the following Newton iteration:

       x' = x + (1 - a * x^2) * x / 2;

     which converges to 1/sqrt(a), starting with the
     float precision approximation to 1/sqrt(a).
     Since Newton's iteration more or less doubles the
     number of correct digits, we only need to perform it
     twice.
  */

  if (is_zero(a))
    return 0.0;

  if (is_negative(a)) {
    return qd_nan;
  }

  qd_real r = (1.0 / sqrt(a[0]));
  qd_real h = mul_pwr2(a, 0.5);

  r += ((0.5 - h * sqr(r)) * r);
  r += ((0.5 - h * sqr(r)) * r);
  r += ((0.5 - h * sqr(r)) * r);

  r *= a;
  return r;
}


/* Computes the n-th root of a */
qd_real nroot(qd_real a, int n) {
  /* Strategy:  Use Newton's iteration to solve

        1/(x^n) - a = 0

     Newton iteration becomes

        x' = x + x * (1 - a * x^n) / n

     Since Newton's iteration converges quadratically,
     we only need to perform it twice.

   */
  if (n <= 0) {
    return qd_nan;
  }

  if (n % 2 == 0 && is_negative(a)) {
    return qd_nan;
  }

  if (n == 1) {
    return a;
  }
  if (n == 2) {
    return sqrt(a);
  }
  if (is_zero(a)) {
    return qd_real(0.0);
  }


  /* Note  a^{-1/n} = exp(-log(a)/n) */
  qd_real r = abs(a);
  qd_real x = exp(-log(r.x[0]) / n);

  /* Perform Newton's iteration. */
  float dbl_n = float(n);
  x += x * (1.0 - r * npwr(x, n)) / dbl_n;
  x += x * (1.0 - r * npwr(x, n)) / dbl_n;
  x += x * (1.0 - r * npwr(x, n)) / dbl_n;
  if (a[0] < 0.0){
    x = -x;
  }
  return 1.0 / x;
}

#if 0
const int n_inv_fact = 15;
const qd_real inv_fact[n_inv_fact] = {
  qd_real( 1.66666666666666657e-01,  9.25185853854297066e-18,
           5.13581318503262866e-34,  2.85094902409834186e-50),
  qd_real( 4.16666666666666644e-02,  2.31296463463574266e-18,
           1.28395329625815716e-34,  7.12737256024585466e-51),
  qd_real( 8.33333333333333322e-03,  1.15648231731787138e-19,
           1.60494162032269652e-36,  2.22730392507682967e-53),
  qd_real( 1.38888888888888894e-03, -5.30054395437357706e-20,
          -1.73868675534958776e-36, -1.63335621172300840e-52),
  qd_real( 1.98412698412698413e-04,  1.72095582934207053e-22,
           1.49269123913941271e-40,  1.29470326746002471e-58),
  qd_real( 2.48015873015873016e-05,  2.15119478667758816e-23,
           1.86586404892426588e-41,  1.61837908432503088e-59),
  qd_real( 2.75573192239858925e-06, -1.85839327404647208e-22,
           8.49175460488199287e-39, -5.72661640789429621e-55),
  qd_real( 2.75573192239858883e-07,  2.37677146222502973e-23,
          -3.26318890334088294e-40,  1.61435111860404415e-56),
  qd_real( 2.50521083854417202e-08, -1.44881407093591197e-24,
           2.04267351467144546e-41, -8.49632672007163175e-58),
  qd_real( 2.08767569878681002e-09, -1.20734505911325997e-25,
           1.70222792889287100e-42,  1.41609532150396700e-58),
  qd_real( 1.60590438368216133e-10,  1.25852945887520981e-26,
          -5.31334602762985031e-43,  3.54021472597605528e-59),
  qd_real( 1.14707455977297245e-11,  2.06555127528307454e-28,
           6.88907923246664603e-45,  5.72920002655109095e-61),
  qd_real( 7.64716373181981641e-13,  7.03872877733453001e-30,
          -7.82753927716258345e-48,  1.92138649443790242e-64),
  qd_real( 4.77947733238738525e-14,  4.39920548583408126e-31,
          -4.89221204822661465e-49,  1.20086655902368901e-65),
  qd_real( 2.81145725434552060e-15,  1.65088427308614326e-31,
          -2.87777179307447918e-50,  4.27110689256293549e-67)
};
#endif

qd_real exp(qd_real a) {
  /* Strategy:  We first reduce the size of x by noting that

          exp(kr + m * log(2)) = 2^m * exp(r)^k

     where m and k are integers.  By choosing m appropriately
     we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
     evaluated using the familiar Taylor series.  Reducing the
     argument substantially speeds up the convergence.       */

  const float k = ldexp(1.0, 16);
  const float inv_k = 1.0 / k;

  if (a[0] <= -709.0)
    return qd_0;

  if (a[0] >=  709.0)
    return qd_inf;

  if (is_zero(a))
    return qd_1;

  if (is_one(a))
    return qd_e;

  float m = floor(a.x[0] / qd_log2.x[0] + 0.5);
  qd_real r = mul_pwr2(a - qd_log2 * m, inv_k);
  qd_real s, p, t;
  float thresh = inv_k * qd_eps;

  p = sqr(r);
  s = r + mul_pwr2(p, 0.5);
  int i = 0;
  do {
    p *= r;
    t = p * inv_fact[i++];
    s += t;
  } while (abs(to_float(t)) > thresh && i < 9);

  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s += 1.0;
  return ldexp(s, int(m));
}

/* Logarithm.  Computes log(x) in quad-float precision.
   This is a natural logarithm (i.e., base e).            */
qd_real log(qd_real a) {
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x)
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.

     Two iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

  if (is_one(a)) {
    return qd_0;
  }

  if (a[0] <= 0.0) {
    return qd_nan;
  }

  if (a[0] == 0.0) {
    return -qd_inf;
  }

  qd_real x = log(a[0]);   /* Initial approximation */

  x = x + a * exp(-x) - 1.0;
  x = x + a * exp(-x) - 1.0;
  x = x + a * exp(-x) - 1.0;

  return x;
}

qd_real log10(qd_real a) {
  return log(a) / qd_log10;
}

#if 0
const qd_real _pi1024 = qd_real(
    3.067961575771282340e-03, 1.195944139792337116e-19,
   -2.924579892303066080e-36, 1.086381075061880158e-52);

/* Table of sin(k * pi/1024) and cos(k * pi/1024). */
const qd_real sin_table [] = {
  qd_real( 3.0679567629659761e-03, 1.2690279085455925e-19,
       5.2879464245328389e-36, -1.7820334081955298e-52),
  qd_real( 6.1358846491544753e-03, 9.0545257482474933e-20,
       1.6260113133745320e-37, -9.7492001208767410e-55),
  qd_real( 9.2037547820598194e-03, -1.2136591693535934e-19,
       5.5696903949425567e-36, 1.2505635791936951e-52),
  qd_real( 1.2271538285719925e-02, 6.9197907640283170e-19,
       -4.0203726713435555e-36, -2.0688703606952816e-52),
  qd_real( 1.5339206284988102e-02, -8.4462578865401696e-19,
       4.6535897505058629e-35, -1.3923682978570467e-51),
  qd_real( 1.8406729905804820e-02, 7.4195533812833160e-19,
       3.9068476486787607e-35, 3.6393321292898614e-52),
  qd_real( 2.1474080275469508e-02, -4.5407960207688566e-19,
       -2.2031770119723005e-35, 1.2709814654833741e-51),
  qd_real( 2.4541228522912288e-02, -9.1868490125778782e-20,
       4.8706148704467061e-36, -2.8153947855469224e-52),
  qd_real( 2.7608145778965743e-02, -1.5932358831389269e-18,
       -7.0475416242776030e-35, -2.7518494176602744e-51),
  qd_real( 3.0674803176636626e-02, -1.6936054844107918e-20,
       -2.0039543064442544e-36, -1.6267505108658196e-52),
  qd_real( 3.3741171851377587e-02, -2.0096074292368340e-18,
       -1.3548237016537134e-34, 6.5554881875899973e-51),
  qd_real( 3.6807222941358832e-02, 6.1060088803529842e-19,
       -4.0448721259852727e-35, -2.1111056765671495e-51),
  qd_real( 3.9872927587739811e-02, 4.6657453481183289e-19,
       3.4119333562288684e-35, 2.4007534726187511e-51),
  qd_real( 4.2938256934940820e-02, 2.8351940588660907e-18,
       1.6991309601186475e-34, 6.8026536098672629e-51),
  qd_real( 4.6003182130914630e-02, -1.1182813940157788e-18,
       7.5235020270378946e-35, 4.1187304955493722e-52),
  qd_real( 4.9067674327418015e-02, -6.7961037205182801e-19,
       -4.4318868124718325e-35, -9.9376628132525316e-52),
  qd_real( 5.2131704680283324e-02, -2.4243695291953779e-18,
       -1.3675405320092298e-34, -8.3938137621145070e-51),
  qd_real( 5.5195244349689941e-02, -1.3340299860891103e-18,
       -3.4359574125665608e-35, 1.1911462755409369e-51),
  qd_real( 5.8258264500435759e-02, 2.3299905496077492e-19,
       1.9376108990628660e-36, -5.1273775710095301e-53),
  qd_real( 6.1320736302208578e-02, -5.1181134064638108e-19,
       -4.2726335866706313e-35, 2.6368495557440691e-51),
  qd_real( 6.4382630929857465e-02, -4.2325997000052705e-18,
       3.3260117711855937e-35, 1.4736267706718352e-51),
  qd_real( 6.7443919563664065e-02, -6.9221796556983636e-18,
       1.5909286358911040e-34, -7.8828946891835218e-51),
  qd_real( 7.0504573389613870e-02, -6.8552791107342883e-18,
       -1.9961177630841580e-34, 2.0127129580485300e-50),
  qd_real( 7.3564563599667426e-02, -2.7784941506273593e-18,
       -9.1240375489852821e-35, -1.9589752023546795e-51),
  qd_real( 7.6623861392031492e-02, 2.3253700287958801e-19,
       -1.3186083921213440e-36, -4.9927872608099673e-53),
  qd_real( 7.9682437971430126e-02, -4.4867664311373041e-18,
       2.8540789143650264e-34, 2.8491348583262741e-51),
  qd_real( 8.2740264549375692e-02, 1.4735983530877760e-18,
       3.7284093452233713e-35, 2.9024430036724088e-52),
  qd_real( 8.5797312344439894e-02, -3.3881893830684029e-18,
       -1.6135529531508258e-34, 7.7294651620588049e-51),
  qd_real( 8.8853552582524600e-02, -3.7501775830290691e-18,
       3.7543606373911573e-34, 2.2233701854451859e-50),
  qd_real( 9.1908956497132724e-02, 4.7631594854274564e-18,
       1.5722874642939344e-34, -4.8464145447831456e-51),
  qd_real( 9.4963495329639006e-02, -6.5885886400417564e-18,
       -2.1371116991641965e-34, 1.3819370559249300e-50),
  qd_real( 9.8017140329560604e-02, -1.6345823622442560e-18,
       -1.3209238810006454e-35, -3.5691060049117942e-52),
  qd_real( 1.0106986275482782e-01, 3.3164325719308656e-18,
       -1.2004224885132282e-34, 7.2028828495418631e-51),
  qd_real( 1.0412163387205457e-01, 6.5760254085385100e-18,
       1.7066246171219214e-34, -4.9499340996893514e-51),
  qd_real( 1.0717242495680884e-01, 6.4424044279026198e-18,
       -8.3956976499698139e-35, -4.0667730213318321e-51),
  qd_real( 1.1022220729388306e-01, -5.6789503537823233e-19,
       1.0380274792383233e-35, 1.5213997918456695e-52),
  qd_real( 1.1327095217756435e-01, 2.7100481012132900e-18,
       1.5323292999491619e-35, 4.9564432810360879e-52),
  qd_real( 1.1631863091190477e-01, 1.0294914877509705e-18,
       -9.3975734948993038e-35, 1.3534827323719708e-52),
  qd_real( 1.1936521481099137e-01, -3.9500089391898506e-18,
       3.5317349978227311e-34, 1.8856046807012275e-51),
  qd_real( 1.2241067519921620e-01, 2.8354501489965335e-18,
       1.8151655751493305e-34, -2.8716592177915192e-51),
  qd_real( 1.2545498341154623e-01, 4.8686751763148235e-18,
       5.9878105258097936e-35, -3.3534629098722107e-51),
  qd_real( 1.2849811079379317e-01, 3.8198603954988802e-18,
       -1.8627501455947798e-34, -2.4308161133527791e-51),
  qd_real( 1.3154002870288312e-01, -5.0039708262213813e-18,
       -1.2983004159245552e-34, -4.6872034915794122e-51),
  qd_real( 1.3458070850712620e-01, -9.1670359171480699e-18,
       1.5916493007073973e-34, 4.0237002484366833e-51),
  qd_real( 1.3762012158648604e-01, 6.6253255866774482e-18,
       -2.3746583031401459e-34, -9.3703876173093250e-52),
  qd_real( 1.4065823933284924e-01, -7.9193932965524741e-18,
       6.0972464202108397e-34, 2.4566623241035797e-50),
  qd_real( 1.4369503315029444e-01, 1.1472723016618666e-17,
       -5.1884954557576435e-35, -4.2220684832186607e-51),
  qd_real( 1.4673047445536175e-01, 3.7269471470465677e-18,
       3.7352398151250827e-34, -4.0881822289508634e-51),
  qd_real( 1.4976453467732151e-01, 8.0812114131285151e-18,
       1.2979142554917325e-34, 9.9380667487736254e-51),
  qd_real( 1.5279718525844344e-01, -7.6313573938416838e-18,
       5.7714690450284125e-34, -3.7731132582986687e-50),
  qd_real( 1.5582839765426523e-01, 3.0351307187678221e-18,
       -1.0976942315176184e-34, 7.8734647685257867e-51),
  qd_real( 1.5885814333386145e-01, -4.0163200573859079e-18,
       -9.2840580257628812e-35, -2.8567420029274875e-51),
  qd_real( 1.6188639378011183e-01, 1.1850519643573528e-17,
       -5.0440990519162957e-34, 3.0510028707928009e-50),
  qd_real( 1.6491312048996992e-01, -7.0405288319166738e-19,
       3.3211107491245527e-35, 8.6663299254686031e-52),
  qd_real( 1.6793829497473117e-01, 5.4284533721558139e-18,
       -3.3263339336181369e-34, -1.8536367335123848e-50),
  qd_real( 1.7096188876030122e-01, 9.1919980181759094e-18,
       -6.7688743940982606e-34, -1.0377711384318389e-50),
  qd_real( 1.7398387338746382e-01, 5.8151994618107928e-18,
       -1.6751014298301606e-34, -6.6982259797164963e-51),
  qd_real( 1.7700422041214875e-01, 6.7329300635408167e-18,
       2.8042736644246623e-34, 3.6786888232793599e-51),
  qd_real( 1.8002290140569951e-01, 7.9701826047392143e-18,
       -7.0765920110524977e-34, 1.9622512608461784e-50),
  qd_real( 1.8303988795514095e-01, 7.7349918688637383e-18,
       -4.4803769968145083e-34, 1.1201148793328890e-50),
  qd_real( 1.8605515166344666e-01, -1.2564893007679552e-17,
       7.5953844248530810e-34, -3.8471695132415039e-51),
  qd_real( 1.8906866414980622e-01, -7.6208955803527778e-18,
       -4.4792298656662981e-34, -4.4136824096645007e-50),
  qd_real( 1.9208039704989244e-01, 4.3348343941174903e-18,
       -2.3404121848139937e-34, 1.5789970962611856e-50),
  qd_real( 1.9509032201612828e-01, -7.9910790684617313e-18,
       6.1846270024220713e-34, -3.5840270918032937e-50),
  qd_real( 1.9809841071795359e-01, -1.8434411800689445e-18,
       1.4139031318237285e-34, 1.0542811125343809e-50),
  qd_real( 2.0110463484209190e-01, 1.1010032669300739e-17,
       -3.9123576757413791e-34, 2.4084852500063531e-51),
  qd_real( 2.0410896609281687e-01, 6.0941297773957752e-18,
       -2.8275409970449641e-34, 4.6101008563532989e-51),
  qd_real( 2.0711137619221856e-01, -1.0613362528971356e-17,
       2.2456805112690884e-34, 1.3483736125280904e-50),
  qd_real( 2.1011183688046961e-01, 1.1561548476512844e-17,
       6.0355905610401254e-34, 3.3329909618405675e-50),
  qd_real( 2.1311031991609136e-01, 1.2031873821063860e-17,
       -3.4142699719695635e-34, -1.2436262780241778e-50),
  qd_real( 2.1610679707621952e-01, -1.0111196082609117e-17,
       7.2789545335189643e-34, -2.9347540365258610e-50),
  qd_real( 2.1910124015686980e-01, -3.6513812299150776e-19,
       -2.3359499418606442e-35, 3.1785298198458653e-52),
  qd_real( 2.2209362097320354e-01, -3.0337210995812162e-18,
       6.6654668033632998e-35, 2.0110862322656942e-51),
  qd_real( 2.2508391135979283e-01, 3.9507040822556510e-18,
       2.4287993958305375e-35, 5.6662797513020322e-52),
  qd_real( 2.2807208317088573e-01, 8.2361837339258012e-18,
       6.9786781316397937e-34, -6.4122962482639504e-51),
  qd_real( 2.3105810828067111e-01, 1.0129787149761869e-17,
       -6.9359234615816044e-34, -2.8877355604883782e-50),
  qd_real( 2.3404195858354343e-01, -6.9922402696101173e-18,
       -5.7323031922750280e-34, 5.3092579966872727e-51),
  qd_real( 2.3702360599436720e-01, 8.8544852285039918e-18,
       1.3588480826354134e-34, 1.0381022520213867e-50),
  qd_real( 2.4000302244874150e-01, -1.2137758975632164e-17,
       -2.6448807731703891e-34, -1.9929733800670473e-51),
  qd_real( 2.4298017990326390e-01, -8.7514315297196632e-18,
       -6.5723260373079431e-34, -1.0333158083172177e-50),
  qd_real( 2.4595505033579462e-01, -1.1129044052741832e-17,
       4.3805998202883397e-34, 1.2219399554686291e-50),
  qd_real( 2.4892760574572018e-01, -8.1783436100020990e-18,
       5.5666875261111840e-34, 3.8080473058748167e-50),
  qd_real( 2.5189781815421697e-01, -1.7591436032517039e-17,
       -1.0959681232525285e-33, 5.6209426020232456e-50),
  qd_real( 2.5486565960451457e-01, -1.3602299806901461e-19,
       -6.0073844642762535e-36, -3.0072751311893878e-52),
  qd_real( 2.5783110216215899e-01, 1.8480038630879957e-17,
       3.3201664714047599e-34, -5.5547819290576764e-51),
  qd_real( 2.6079411791527551e-01, 4.2721420983550075e-18,
       5.6782126934777920e-35, 3.1428338084365397e-51),
  qd_real( 2.6375467897483140e-01, -1.8837947680038700e-17,
       1.3720129045754794e-33, -8.2763406665966033e-50),
  qd_real( 2.6671275747489837e-01, 2.0941222578826688e-17,
       -1.1303466524727989e-33, 1.9954224050508963e-50),
  qd_real( 2.6966832557291509e-01, 1.5765657618133259e-17,
       -6.9696142173370086e-34, -4.0455346879146776e-50),
  qd_real( 2.7262135544994898e-01, 7.8697166076387850e-18,
       6.6179388602933372e-35, -2.7642903696386267e-51),
  qd_real( 2.7557181931095814e-01, 1.9320328962556582e-17,
       1.3932094180100280e-33, 1.3617253920018116e-50),
  qd_real( 2.7851968938505312e-01, -1.0030273719543544e-17,
       7.2592115325689254e-34, -1.0068516296655851e-50),
  qd_real( 2.8146493792575800e-01, -1.2322299641274009e-17,
       -1.0564788706386435e-34, 7.5137424251265885e-51),
  qd_real( 2.8440753721127182e-01, 2.2209268510661475e-17,
       -9.1823095629523708e-34, -5.2192875308892218e-50),
  qd_real( 2.8734745954472951e-01, 1.5461117367645717e-17,
       -6.3263973663444076e-34, -2.2982538416476214e-50),
  qd_real( 2.9028467725446239e-01, -1.8927978707774251e-17,
       1.1522953157142315e-33, 7.4738655654716596e-50),
  qd_real( 2.9321916269425863e-01, 2.2385430811901833e-17,
       1.3662484646539680e-33, -4.2451325253996938e-50),
  qd_real( 2.9615088824362384e-01, -2.0220736360876938e-17,
       -7.9252212533920413e-35, -2.8990577729572470e-51),
  qd_real( 2.9907982630804048e-01, 1.6701181609219447e-18,
       8.6091151117316292e-35, 3.9931286230012102e-52),
  qd_real( 3.0200594931922808e-01, -1.7167666235262474e-17,
       2.3336182149008069e-34, 8.3025334555220004e-51),
  qd_real( 3.0492922973540243e-01, -2.2989033898191262e-17,
       -1.4598901099661133e-34, 3.7760487693121827e-51),
  qd_real( 3.0784964004153487e-01, 2.7074088527245185e-17,
       1.2568858206899284e-33, 7.2931815105901645e-50),
  qd_real( 3.1076715274961147e-01, 2.0887076364048513e-17,
       -3.0130590791065942e-34, 1.3876739009935179e-51),
  qd_real( 3.1368174039889146e-01, 1.4560447299968912e-17,
       3.6564186898011595e-34, 1.1654264734999375e-50),
  qd_real( 3.1659337555616585e-01, 2.1435292512726283e-17,
       1.2338169231377316e-33, 3.3963542100989293e-50),
  qd_real( 3.1950203081601569e-01, -1.3981562491096626e-17,
       8.1730000697411350e-34, -7.7671096270210952e-50),
  qd_real( 3.2240767880106985e-01, -4.0519039937959398e-18,
       3.7438302780296796e-34, 8.7936731046639195e-51),
  qd_real( 3.2531029216226293e-01, 7.9171249463765892e-18,
       -6.7576622068146391e-35, 2.3021655066929538e-51),
  qd_real( 3.2820984357909255e-01, -2.6693140719641896e-17,
       7.8928851447534788e-34, 2.5525163821987809e-51),
  qd_real( 3.3110630575987643e-01, -2.7469465474778694e-17,
       -1.3401245916610206e-33, 6.5531762489976163e-50),
  qd_real( 3.3399965144200938e-01, 2.2598986806288142e-17,
       7.8063057192586115e-34, 2.0427600895486683e-50),
  qd_real( 3.3688985339222005e-01, -4.2000940033475092e-19,
       -2.9178652969985438e-36, -1.1597376437036749e-52),
  qd_real( 3.3977688440682685e-01, 6.6028679499418282e-18,
       1.2575009988669683e-34, 2.5569067699008304e-51),
  qd_real( 3.4266071731199438e-01, 1.9261518449306319e-17,
       -9.2754189135990867e-34, 8.5439996687390166e-50),
  qd_real( 3.4554132496398904e-01, 2.7251143672916123e-17,
       7.0138163601941737e-34, -1.4176292197454015e-50),
  qd_real( 3.4841868024943456e-01, 3.6974420514204918e-18,
       3.5532146878499996e-34, 1.9565462544501322e-50),
  qd_real( 3.5129275608556715e-01, -2.2670712098795844e-17,
       -1.6994216673139631e-34, -1.2271556077284517e-50),
  qd_real( 3.5416352542049040e-01, -1.6951763305764860e-17,
       1.2772331777814617e-33, -3.3703785435843310e-50),
  qd_real( 3.5703096123343003e-01, -4.8218191137919166e-19,
       -4.1672436994492361e-35, -7.1531167149364352e-52),
  qd_real( 3.5989503653498817e-01, -1.7601687123839282e-17,
       1.3375125473046791e-33, 7.9467815593584340e-50),
  qd_real( 3.6275572436739723e-01, -9.1668352663749849e-18,
       -7.4317843956936735e-34, -2.0199582511804564e-50),
  qd_real( 3.6561299780477385e-01, 1.6217898770457546e-17,
       1.1286970151961055e-33, -7.1825287318139010e-50),
  qd_real( 3.6846682995337232e-01, 1.0463640796159268e-17,
       2.0554984738517304e-35, 1.0441861305618769e-51),
  qd_real( 3.7131719395183754e-01, 3.4749239648238266e-19,
       -7.5151053042866671e-37, -2.8153468438650851e-53),
  qd_real( 3.7416406297145799e-01, 8.0114103761962118e-18,
       5.3429599813406052e-34, 1.0351378796539210e-50),
  qd_real( 3.7700741021641826e-01, -2.7255302041956930e-18,
       6.3646586445018137e-35, 8.3048657176503559e-52),
  qd_real( 3.7984720892405116e-01, 9.9151305855172370e-18,
       4.8761409697224886e-34, 1.4025084000776705e-50),
  qd_real( 3.8268343236508978e-01, -1.0050772696461588e-17,
       -2.0605316302806695e-34, -1.2717724698085205e-50),
  qd_real( 3.8551605384391885e-01, 1.5177665396472313e-17,
       1.4198230518016535e-33, 5.8955167159904235e-50),
  qd_real( 3.8834504669882630e-01, -1.0053770598398717e-17,
       7.5942999255057131e-34, -3.1967974046654219e-50),
  qd_real( 3.9117038430225387e-01, 1.7997787858243995e-17,
       -1.0613482402609856e-33, -5.4582148817791032e-50),
  qd_real( 3.9399204006104810e-01, 9.7649241641239336e-18,
       -2.1233599441284617e-34, -5.5529836795340819e-51),
  qd_real( 3.9680998741671031e-01, 2.0545063670840126e-17,
       6.1347058801922842e-34, 1.0733788150636430e-50),
  qd_real( 3.9962419984564684e-01, -1.5065497476189372e-17,
       -9.9653258881867298e-34, -5.7524323712725355e-50),
  qd_real( 4.0243465085941843e-01, 1.0902619339328270e-17,
       7.3998528125989765e-34, 2.2745784806823499e-50),
  qd_real( 4.0524131400498986e-01, 9.9111401942899884e-18,
       -2.5169070895434648e-34, 9.2772984818436573e-53),
  qd_real( 4.0804416286497869e-01, -7.0006015137351311e-18,
       -1.4108207334268228e-34, 1.5175546997577136e-52),
  qd_real( 4.1084317105790397e-01, -2.4219835190355499e-17,
       -1.1418902925313314e-33, -2.0996843165093468e-50),
  qd_real( 4.1363831223843456e-01, -1.0393984940597871e-17,
       -1.1481681174503880e-34, -2.0281052851028680e-51),
  qd_real( 4.1642956009763721e-01, -2.5475580413131732e-17,
       -3.4482678506112824e-34, 7.1788619351865480e-51),
  qd_real( 4.1921688836322396e-01, -4.2232463750110590e-18,
       -3.6053023045255790e-34, -2.2209673210025631e-50),
  qd_real( 4.2200027079979968e-01, 4.3543266994128527e-18,
       3.1734310272251190e-34, -1.3573247980738668e-50),
  qd_real( 4.2477968120910881e-01, 2.7462312204277281e-17,
       -4.6552847802111948e-34, 6.5961781099193122e-51),
  qd_real( 4.2755509343028208e-01, 9.4111898162954726e-18,
       -1.7446682426598801e-34, -2.2054492626480169e-51),
  qd_real( 4.3032648134008261e-01, 2.2259686974092690e-17,
       8.5972591314085075e-34, -2.9420897889003020e-50),
  qd_real( 4.3309381885315196e-01, 1.1224283329847517e-17,
       5.3223748041075651e-35, 5.3926192627014212e-51),
  qd_real( 4.3585707992225547e-01, 1.6230515450644527e-17,
       -6.4371449063579431e-35, -6.9102436481386757e-51),
  qd_real( 4.3861623853852766e-01, -2.0883315831075090e-17,
       -1.4259583540891877e-34, 6.3864763590657077e-52),
  qd_real( 4.4137126873171667e-01, 2.2360783886964969e-17,
       1.1864769603515770e-34, -3.8087003266189232e-51),
  qd_real( 4.4412214457042926e-01, -2.4218874422178315e-17,
       2.2205230838703907e-34, 9.2133035911356258e-51),
  qd_real( 4.4686884016237421e-01, -1.9222136142309382e-17,
       -4.4425678589732049e-35, -1.3673609292149535e-51),
  qd_real( 4.4961132965460660e-01, 4.8831924232035243e-18,
       2.7151084498191381e-34, -1.5653993171613154e-50),
  qd_real( 4.5234958723377089e-01, -1.4827977472196122e-17,
       -7.6947501088972324e-34, 1.7656856882031319e-50),
  qd_real( 4.5508358712634384e-01, -1.2379906758116472e-17,
       5.5289688955542643e-34, -8.5382312840209386e-51),
  qd_real( 4.5781330359887723e-01, -8.4554254922295949e-18,
       -6.3770394246764263e-34, 3.1778253575564249e-50),
  qd_real( 4.6053871095824001e-01, 1.8488777492177872e-17,
       -1.0527732154209725e-33, 3.3235593490947102e-50),
  qd_real( 4.6325978355186020e-01, -7.3514924533231707e-18,
       6.7175396881707035e-34, 3.9594127612123379e-50),
  qd_real( 4.6597649576796618e-01, -3.3023547778235135e-18,
       3.4904677050476886e-35, 3.4483855263874246e-51),
  qd_real( 4.6868882203582796e-01, -2.2949251681845054e-17,
       -1.1364757641823658e-33, 6.8840522501918612e-50),
  qd_real( 4.7139673682599764e-01, 6.5166781360690130e-18,
       2.9457546966235984e-34, -6.2159717738836630e-51),
  qd_real( 4.7410021465055002e-01, -8.1451601548978075e-18,
       -3.4789448555614422e-34, -1.1681943974658508e-50),
  qd_real( 4.7679923006332214e-01, -1.0293515338305794e-17,
       -3.6582045008369952e-34, 1.7424131479176475e-50),
  qd_real( 4.7949375766015301e-01, 1.8419999662684771e-17,
       -1.3040838621273312e-33, 1.0977131822246471e-50),
  qd_real( 4.8218377207912277e-01, -2.5861500925520442e-17,
       -6.2913197606500007e-36, 4.0802359808684726e-52),
  qd_real( 4.8486924800079112e-01, -1.8034004203262245e-17,
       -3.5244276906958044e-34, -1.7138318654749246e-50),
  qd_real( 4.8755016014843594e-01, 1.4231090931273653e-17,
       -1.8277733073262697e-34, -1.5208291790429557e-51),
  qd_real( 4.9022648328829116e-01, -5.1496145643440404e-18,
       -3.6903027405284104e-34, 1.5172940095151304e-50),
  qd_real( 4.9289819222978404e-01, -1.0257831676562186e-18,
       6.9520817760885069e-35, -2.4260961214090389e-51),
  qd_real( 4.9556526182577254e-01, -9.4323241942365362e-18,
       3.1212918657699143e-35, 4.2009072375242736e-52),
  qd_real( 4.9822766697278187e-01, -1.6126383830540798e-17,
       -1.5092897319298871e-33, 1.1049298890895917e-50),
  qd_real( 5.0088538261124083e-01, -3.9604015147074639e-17,
       -2.2208395201898007e-33, 1.3648202735839417e-49),
  qd_real( 5.0353838372571758e-01, -1.6731308204967497e-17,
       -1.0140233644074786e-33, 4.0953071937671477e-50),
  qd_real( 5.0618664534515534e-01, -4.8321592986493711e-17,
       9.2858107226642252e-34, 4.2699802401037005e-50),
  qd_real( 5.0883014254310699e-01, 4.7836968268014130e-17,
       -1.0727022928806035e-33, 2.7309374513672757e-50),
  qd_real( 5.1146885043797041e-01, -1.3088001221007579e-17,
       4.0929033363366899e-34, -3.7952190153477926e-50),
  qd_real( 5.1410274419322177e-01, -4.5712707523615624e-17,
       1.5488279442238283e-33, -2.5853959305521130e-50),
  qd_real( 5.1673179901764987e-01, 8.3018617233836515e-18,
       5.8251027467695202e-34, -2.2812397190535076e-50),
  qd_real( 5.1935599016558964e-01, -5.5331248144171145e-17,
       -3.1628375609769026e-35, -2.4091972051188571e-51),
  qd_real( 5.2197529293715439e-01, -4.6555795692088883e-17,
       4.6378980936850430e-34, -3.3470542934689532e-51),
  qd_real( 5.2458968267846895e-01, -4.3068869040082345e-17,
       -4.2013155291932055e-34, -1.5096069926700274e-50),
  qd_real( 5.2719913478190139e-01, -4.2202983480560619e-17,
       8.5585916184867295e-34, 7.9974339336732307e-50),
  qd_real( 5.2980362468629472e-01, -4.8067841706482342e-17,
       5.8309721046630296e-34, -8.9740761521756660e-51),
  qd_real( 5.3240312787719801e-01, -4.1020306135800895e-17,
       -1.9239996374230821e-33, -1.5326987913812184e-49),
  qd_real( 5.3499761988709726e-01, -5.3683132708358134e-17,
       -1.3900569918838112e-33, 2.7154084726474092e-50),
  qd_real( 5.3758707629564551e-01, -2.2617365388403054e-17,
       -5.9787279033447075e-34, 3.1204419729043625e-51),
  qd_real( 5.4017147272989285e-01, 2.7072447965935839e-17,
       1.1698799709213829e-33, -5.9094668515881500e-50),
  qd_real( 5.4275078486451589e-01, 1.7148261004757101e-17,
       -1.3525905925200870e-33, 4.9724411290727323e-50),
  qd_real( 5.4532498842204646e-01, -4.1517817538384258e-17,
       -1.5318930219385941e-33, 6.3629921101413974e-50),
  qd_real( 5.4789405917310019e-01, -2.4065878297113363e-17,
       -3.5639213669362606e-36, -2.6013270854271645e-52),
  qd_real( 5.5045797293660481e-01, -8.3319903015807663e-18,
       -2.3058454035767633e-34, -2.1611290432369010e-50),
  qd_real( 5.5301670558002758e-01, -4.7061536623798204e-17,
       -1.0617111545918056e-33, -1.6196316144407379e-50),
  qd_real( 5.5557023301960218e-01, 4.7094109405616768e-17,
       -2.0640520383682921e-33, 1.2290163188567138e-49),
  qd_real( 5.5811853122055610e-01, 1.3481176324765226e-17,
       -5.5016743873011438e-34, -2.3484822739335416e-50),
  qd_real( 5.6066157619733603e-01, -7.3956418153476152e-18,
       3.9680620611731193e-34, 3.1995952200836223e-50),
  qd_real( 5.6319934401383409e-01, 2.3835775146854829e-17,
       1.3511793173769814e-34, 9.3201311581248143e-51),
  qd_real( 5.6573181078361323e-01, -3.4096079596590466e-17,
       -1.7073289744303546e-33, 8.9147089975404507e-50),
  qd_real( 5.6825895267013160e-01, -5.0935673642769248e-17,
       -1.6274356351028249e-33, 9.8183151561702966e-51),
  qd_real( 5.7078074588696726e-01, 2.4568151455566208e-17,
       -1.2844481247560350e-33, -1.8037634376936261e-50),
  qd_real( 5.7329716669804220e-01, 8.5176611669306400e-18,
       -6.4443208788026766e-34, 2.2546105543273003e-50),
  qd_real( 5.7580819141784534e-01, -3.7909495458942734e-17,
       -2.7433738046854309e-33, 1.1130841524216795e-49),
  qd_real( 5.7831379641165559e-01, -2.6237691512372831e-17,
       1.3679051680738167e-33, -3.1409808935335900e-50),
  qd_real( 5.8081395809576453e-01, 1.8585338586613408e-17,
       2.7673843114549181e-34, 1.9605349619836937e-50),
  qd_real( 5.8330865293769829e-01, 3.4516601079044858e-18,
       1.8065977478946306e-34, -6.3953958038544646e-51),
  qd_real( 5.8579785745643886e-01, -3.7485501964311294e-18,
       2.7965403775536614e-34, -7.1816936024157202e-51),
  qd_real( 5.8828154822264533e-01, -2.9292166725006846e-17,
       -2.3744954603693934e-33, -1.1571631191512480e-50),
  qd_real( 5.9075970185887428e-01, -4.7013584170659542e-17,
       2.4808417611768356e-33, 1.2598907673643198e-50),
  qd_real( 5.9323229503979980e-01, 1.2892320944189053e-17,
       5.3058364776359583e-34, 4.1141674699390052e-50),
  qd_real( 5.9569930449243336e-01, -1.3438641936579467e-17,
       -6.7877687907721049e-35, -5.6046937531684890e-51),
  qd_real( 5.9816070699634227e-01, 3.8801885783000657e-17,
       -1.2084165858094663e-33, -4.0456610843430061e-50),
  qd_real( 6.0061647938386897e-01, -4.6398198229461932e-17,
       -1.6673493003710801e-33, 5.1982824378491445e-50),
  qd_real( 6.0306659854034816e-01, 3.7323357680559650e-17,
       2.7771920866974305e-33, -1.6194229649742458e-49),
  qd_real( 6.0551104140432555e-01, -3.1202672493305677e-17,
       1.2761267338680916e-33, -4.0859368598379647e-50),
  qd_real( 6.0794978496777363e-01, 3.5160832362096660e-17,
       -2.5546242776778394e-34, -1.4085313551220694e-50),
  qd_real( 6.1038280627630948e-01, -2.2563265648229169e-17,
       1.3185575011226730e-33, 8.2316691420063460e-50),
  qd_real( 6.1281008242940971e-01, -4.2693476568409685e-18,
       2.5839965886650320e-34, 1.6884412005622537e-50),
  qd_real( 6.1523159058062682e-01, 2.6231417767266950e-17,
       -1.4095366621106716e-33, 7.2058690491304558e-50),
  qd_real( 6.1764730793780398e-01, -4.7478594510902452e-17,
       -7.2986558263123996e-34, -3.0152327517439154e-50),
  qd_real( 6.2005721176328921e-01, -2.7983410837681118e-17,
       1.1649951056138923e-33, -5.4539089117135207e-50),
  qd_real( 6.2246127937414997e-01, 5.2940728606573002e-18,
       -4.8486411215945827e-35, 1.2696527641980109e-52),
  qd_real( 6.2485948814238634e-01, 3.3671846037243900e-17,
       -2.7846053391012096e-33, 5.6102718120012104e-50),
  qd_real( 6.2725181549514408e-01, 3.0763585181253225e-17,
       2.7068930273498138e-34, -1.1172240309286484e-50),
  qd_real( 6.2963823891492698e-01, 4.1115334049626806e-17,
       -1.9167473580230747e-33, 1.1118424028161730e-49),
  qd_real( 6.3201873593980906e-01, -4.0164942296463612e-17,
       -7.2208643641736723e-34, 3.7828920470544344e-50),
  qd_real( 6.3439328416364549e-01, 1.0420901929280035e-17,
       4.1174558929280492e-34, -1.4464152986630705e-51),
  qd_real( 6.3676186123628420e-01, 3.1419048711901611e-17,
       -2.2693738415126449e-33, -1.6023584204297388e-49),
  qd_real( 6.3912444486377573e-01, 1.2416796312271043e-17,
       -6.2095419626356605e-34, 2.7762065999506603e-50),
  qd_real( 6.4148101280858316e-01, -9.9883430115943310e-18,
       4.1969230376730128e-34, 5.6980543799257597e-51),
  qd_real( 6.4383154288979150e-01, -3.2084798795046886e-17,
       -1.2595311907053305e-33, -4.0205885230841536e-50),
  qd_real( 6.4617601298331639e-01, -2.9756137382280815e-17,
       -1.0275370077518259e-33, 8.0852478665893014e-51),
  qd_real( 6.4851440102211244e-01, 3.9870270313386831e-18,
       1.9408388509540788e-34, -5.1798420636193190e-51),
  qd_real( 6.5084668499638088e-01, 3.9714670710500257e-17,
       2.9178546787002963e-34, 3.8140635508293278e-51),
  qd_real( 6.5317284295377676e-01, 8.5695642060026238e-18,
       -6.9165322305070633e-34, 2.3873751224185395e-50),
  qd_real( 6.5549285299961535e-01, 3.5638734426385005e-17,
       1.2695365790889811e-33, 4.3984952865412050e-50),
  qd_real( 6.5780669329707864e-01, 1.9580943058468545e-17,
       -1.1944272256627192e-33, 2.8556402616436858e-50),
  qd_real( 6.6011434206742048e-01, -1.3960054386823638e-19,
       6.1515777931494047e-36, 5.3510498875622660e-52),
  qd_real( 6.6241577759017178e-01, -2.2615508885764591e-17,
       5.0177050318126862e-34, 2.9162532399530762e-50),
  qd_real( 6.6471097820334490e-01, -3.6227793598034367e-17,
       -9.0607934765540427e-34, 3.0917036342380213e-50),
  qd_real( 6.6699992230363747e-01, 3.5284364997428166e-17,
       -1.0382057232458238e-33, 7.3812756550167626e-50),
  qd_real( 6.6928258834663612e-01, -5.4592652417447913e-17,
       -2.5181014709695152e-33, -1.6867875999437174e-49),
  qd_real( 6.7155895484701844e-01, -4.0489037749296692e-17,
       3.1995835625355681e-34, -1.4044414655670960e-50),
  qd_real( 6.7382900037875604e-01, 2.3091901236161086e-17,
       5.7428037192881319e-34, 1.1240668354625977e-50),
  qd_real( 6.7609270357531592e-01, 3.7256902248049466e-17,
       1.7059417895764375e-33, 9.7326347795300652e-50),
  qd_real( 6.7835004312986147e-01, 1.8302093041863122e-17,
       9.5241675746813072e-34, 5.0328101116133503e-50),
  qd_real( 6.8060099779545302e-01, 2.8473293354522047e-17,
       4.1331805977270903e-34, 4.2579030510748576e-50),
  qd_real( 6.8284554638524808e-01, -1.2958058061524531e-17,
       1.8292386959330698e-34, 3.4536209116044487e-51),
  qd_real( 6.8508366777270036e-01, 2.5948135194645137e-17,
       -8.5030743129500702e-34, -6.9572086141009930e-50),
  qd_real( 6.8731534089175916e-01, -5.5156158714917168e-17,
       1.1896489854266829e-33, -7.8505896218220662e-51),
  qd_real( 6.8954054473706694e-01, -1.5889323294806790e-17,
       9.1242356240205712e-34, 3.8315454152267638e-50),
  qd_real( 6.9175925836415775e-01, 2.7406078472410668e-17,
       1.3286508943202092e-33, 1.0651869129580079e-51),
  qd_real( 6.9397146088965400e-01, 7.4345076956280137e-18,
       7.5061528388197460e-34, -1.5928000240686583e-50),
  qd_real( 6.9617713149146299e-01, -4.1224081213582889e-17,
       -3.1838716762083291e-35, -3.9625587412119131e-51),
  qd_real( 6.9837624940897280e-01, 4.8988282435667768e-17,
       1.9134010413244152e-33, 2.6161153243793989e-50),
  qd_real( 7.0056879394324834e-01, 3.1027960192992922e-17,
       9.5638250509179997e-34, 4.5896916138107048e-51),
  qd_real( 7.0275474445722530e-01, 2.5278294383629822e-18,
       -8.6985561210674942e-35, -5.6899862307812990e-51),
  qd_real( 7.0493408037590488e-01, 2.7608725585748502e-17,
       2.9816599471629137e-34, 1.1533044185111206e-50),
  qd_real( 7.0710678118654757e-01, -4.8336466567264567e-17,
       2.0693376543497068e-33, 2.4677734957341755e-50)
};

const qd_real cos_table [] = {
  qd_real( 9.9999529380957619e-01, -1.9668064285322189e-17,
       -6.3053955095883481e-34, 5.3266110855726731e-52),
  qd_real( 9.9998117528260111e-01, 3.3568103522895585e-17,
       -1.4740132559368063e-35, 9.8603097594755596e-52),
  qd_real( 9.9995764455196390e-01, -3.1527836866647287e-17,
       2.6363251186638437e-33, 1.0007504815488399e-49),
  qd_real( 9.9992470183914450e-01, 3.7931082512668012e-17,
       -8.5099918660501484e-35, -4.9956973223295153e-51),
  qd_real( 9.9988234745421256e-01, -3.5477814872408538e-17,
       1.7102001035303974e-33, -1.0725388519026542e-49),
  qd_real( 9.9983058179582340e-01, 1.8825140517551119e-17,
       -5.1383513457616937e-34, -3.8378827995403787e-50),
  qd_real( 9.9976940535121528e-01, 4.2681177032289012e-17,
       1.9062302359737099e-33, -6.0221153262881160e-50),
  qd_real( 9.9969881869620425e-01, -2.9851486403799753e-17,
       -1.9084787370733737e-33, 5.5980260344029202e-51),
  qd_real( 9.9961882249517864e-01, -4.1181965521424734e-17,
       2.0915365593699916e-33, 8.1403390920903734e-50),
  qd_real( 9.9952941750109314e-01, 2.0517917823755591e-17,
       -4.7673802585706520e-34, -2.9443604198656772e-50),
  qd_real( 9.9943060455546173e-01, 3.9644497752257798e-17,
       -2.3757223716722428e-34, -1.2856759011361726e-51),
  qd_real( 9.9932238458834954e-01, -4.2858538440845682e-17,
       3.3235101605146565e-34, -8.3554272377057543e-51),
  qd_real( 9.9920475861836389e-01, 9.1796317110385693e-18,
       5.5416208185868570e-34, 8.0267046717615311e-52),
  qd_real( 9.9907772775264536e-01, 2.1419007653587032e-17,
       -7.9048203318529618e-34, -5.3166296181112712e-50),
  qd_real( 9.9894129318685687e-01, -2.0610641910058638e-17,
       -1.2546525485913485e-33, -7.5175888806157064e-50),
  qd_real( 9.9879545620517241e-01, -1.2291693337075465e-17,
       2.4468446786491271e-34, 1.0723891085210268e-50),
  qd_real( 9.9864021818026527e-01, -4.8690254312923302e-17,
       -2.9470881967909147e-34, -1.3000650761346907e-50),
  qd_real( 9.9847558057329477e-01, -2.2002931182778795e-17,
       -1.2371509454944992e-33, -2.4911225131232065e-50),
  qd_real( 9.9830154493389289e-01, -5.1869402702792278e-17,
       1.0480195493633452e-33, -2.8995649143155511e-50),
  qd_real( 9.9811811290014918e-01, 2.7935487558113833e-17,
       2.4423341255830345e-33, -6.7646699175334417e-50),
  qd_real( 9.9792528619859600e-01, 1.7143659778886362e-17,
       5.7885840902887460e-34, -9.2601432603894597e-51),
  qd_real( 9.9772306664419164e-01, -2.6394475274898721e-17,
       -1.6176223087661783e-34, -9.9924942889362281e-51),
  qd_real( 9.9751145614030345e-01, 5.6007205919806937e-18,
       -5.9477673514685690e-35, -1.4166807162743627e-54),
  qd_real( 9.9729045667869021e-01, 9.1647695371101735e-18,
       6.7824134309739296e-34, -8.6191392795543357e-52),
  qd_real( 9.9706007033948296e-01, 1.6734093546241963e-17,
       -1.3169951440780028e-33, 1.0311048767952477e-50),
  qd_real( 9.9682029929116567e-01, 4.7062820708615655e-17,
       2.8412041076474937e-33, -8.0006155670263622e-50),
  qd_real( 9.9657114579055484e-01, 1.1707179088390986e-17,
       -7.5934413263024663e-34, 2.8474848436926008e-50),
  qd_real( 9.9631261218277800e-01, 1.1336497891624735e-17,
       3.4002458674414360e-34, 7.7419075921544901e-52),
  qd_real( 9.9604470090125197e-01, 2.2870031707670695e-17,
       -3.9184839405013148e-34, -3.7081260416246375e-50),
  qd_real( 9.9576741446765982e-01, -2.3151908323094359e-17,
       -1.6306512931944591e-34, -1.5925420783863192e-51),
  qd_real( 9.9548075549192694e-01, 3.2084621412226554e-18,
       -4.9501292146013023e-36, -2.7811428850878516e-52),
  qd_real( 9.9518472667219693e-01, -4.2486913678304410e-17,
       1.3315510772504614e-33, 6.7927987417051888e-50),
  qd_real( 9.9487933079480562e-01, 4.2130813284943662e-18,
       -4.2062597488288452e-35, 2.5157064556087620e-51),
  qd_real( 9.9456457073425542e-01, 3.6745069641528058e-17,
       -3.0603304105471010e-33, 1.0397872280487526e-49),
  qd_real( 9.9424044945318790e-01, 4.4129423472462673e-17,
       -3.0107231708238066e-33, 7.4201582906861892e-50),
  qd_real( 9.9390697000235606e-01, -1.8964849471123746e-17,
       -1.5980853777937752e-35, -8.5374807150597082e-52),
  qd_real( 9.9356413552059530e-01, 2.9752309927797428e-17,
       -4.5066707331134233e-34, -3.3548191633805036e-50),
  qd_real( 9.9321194923479450e-01, 3.3096906261272262e-17,
       1.5592811973249567e-33, 1.4373977733253592e-50),
  qd_real( 9.9285041445986510e-01, -1.4094517733693302e-17,
       -1.1954558131616916e-33, 1.8761873742867983e-50),
  qd_real( 9.9247953459870997e-01, 3.1093055095428906e-17,
       -1.8379594757818019e-33, -3.9885758559381314e-51),
  qd_real( 9.9209931314219180e-01, -3.9431926149588778e-17,
       -6.2758062911047230e-34, -1.2960929559212390e-50),
  qd_real( 9.9170975366909953e-01, -2.3372891311883661e-18,
       2.7073298824968591e-35, -1.2569459441802872e-51),
  qd_real( 9.9131085984611544e-01, -2.5192111583372105e-17,
       -1.2852471567380887e-33, 5.2385212584310483e-50),
  qd_real( 9.9090263542778001e-01, 1.5394565094566704e-17,
       -1.0799984133184567e-33, 2.7451115960133595e-51),
  qd_real( 9.9048508425645709e-01, -5.5411437553780867e-17,
       -1.4614017210753585e-33, -3.8339374397387620e-50),
  qd_real( 9.9005821026229712e-01, -1.7055485906233963e-17,
       1.3454939685758777e-33, 7.3117589137300036e-50),
  qd_real( 9.8962201746320089e-01, -5.2398217968132530e-17,
       1.3463189211456219e-33, 5.8021640554894872e-50),
  qd_real( 9.8917650996478101e-01, -4.0987309937047111e-17,
       -4.4857560552048437e-34, -3.9414504502871125e-50),
  qd_real( 9.8872169196032378e-01, -1.0976227206656125e-17,
       3.2311342577653764e-34, 9.6367946583575041e-51),
  qd_real( 9.8825756773074946e-01, 2.7030607784372632e-17,
       7.7514866488601377e-35, 2.1019644956864938e-51),
  qd_real( 9.8778414164457218e-01, -2.3600693397159021e-17,
       -1.2323283769707861e-33, 3.0130900716803339e-50),
  qd_real( 9.8730141815785843e-01, -5.2332261255715652e-17,
       -2.7937644333152473e-33, 1.2074160567958408e-49),
  qd_real( 9.8680940181418553e-01, -5.0287214351061075e-17,
       -2.2681526238144461e-33, 4.4003694320169133e-50),
  qd_real( 9.8630809724459867e-01, -2.1520877103013341e-17,
       1.1866528054187716e-33, -7.8532199199813836e-50),
  qd_real( 9.8579750916756748e-01, -5.1439452979953012e-17,
       2.6276439309996725e-33, 7.5423552783286347e-50),
  qd_real( 9.8527764238894122e-01, 2.3155637027900207e-17,
       -7.5275971545764833e-34, 1.0582231660456094e-50),
  qd_real( 9.8474850180190421e-01, 1.0548144061829957e-17,
       2.8786145266267306e-34, -3.6782210081466112e-51),
  qd_real( 9.8421009238692903e-01, 4.7983922627050691e-17,
       2.2597419645070588e-34, 1.7573875814863400e-50),
  qd_real( 9.8366241921173025e-01, 1.9864948201635255e-17,
       -1.0743046281211033e-35, 1.7975662796558100e-52),
  qd_real( 9.8310548743121629e-01, 4.2170007522888628e-17,
       8.2396265656440904e-34, -8.0803700139096561e-50),
  qd_real( 9.8253930228744124e-01, 1.5149580813777224e-17,
       -4.1802771422186237e-34, -2.2150174326226160e-50),
  qd_real( 9.8196386910955524e-01, 2.1108443711513084e-17,
       -1.5253013442896054e-33, -6.8388082079337969e-50),
  qd_real( 9.8137919331375456e-01, 1.3428163260355633e-17,
       -6.5294290469962986e-34, 2.7965412287456268e-51),
  qd_real( 9.8078528040323043e-01, 1.8546939997825006e-17,
       -1.0696564445530757e-33, 6.6668174475264961e-50),
  qd_real( 9.8018213596811743e-01, -3.6801786963856159e-17,
       6.3245171387992842e-34, 1.8600176137175971e-50),
  qd_real( 9.7956976568544052e-01, 1.5573991584990420e-17,
       -1.3401066029782990e-33, -1.7263702199862149e-50),
  qd_real( 9.7894817531906220e-01, -2.3817727961148053e-18,
       -1.0694750370381661e-34, -8.2293047196087462e-51),
  qd_real( 9.7831737071962765e-01, -2.1623082233344895e-17,
       1.0970403012028032e-33, 7.7091923099369339e-50),
  qd_real( 9.7767735782450993e-01, 5.0514136167059628e-17,
       -1.3254751701428788e-33, 7.0161254312124538e-50),
  qd_real( 9.7702814265775439e-01, -4.3353875751555997e-17,
       5.4948839831535478e-34, -9.2755263105377306e-51),
  qd_real( 9.7636973133002114e-01, 9.3093931526213780e-18,
       -4.1184949155685665e-34, -3.1913926031393690e-50),
  qd_real( 9.7570213003852857e-01, -2.5572556081259686e-17,
       -9.3174244508942223e-34, -8.3675863211646863e-51),
  qd_real( 9.7502534506699412e-01, 2.6642660651899135e-17,
       1.7819392739353853e-34, -3.3159625385648947e-51),
  qd_real( 9.7433938278557586e-01, 2.3041221476151512e-18,
       1.0758686005031430e-34, 5.1074116432809478e-51),
  qd_real( 9.7364424965081198e-01, -5.1729808691005871e-17,
       -1.5508473005989887e-33, -1.6505125917675401e-49),
  qd_real( 9.7293995220556018e-01, -3.1311211122281800e-17,
       -2.6874087789006141e-33, -2.1652434818822145e-51),
  qd_real( 9.7222649707893627e-01, 3.6461169785938221e-17,
       3.0309636883883133e-33, -1.2702716907967306e-51),
  qd_real( 9.7150389098625178e-01, -7.9865421122289046e-18,
       -4.3628417211263380e-34, 3.4307517798759352e-51),
  qd_real( 9.7077214072895035e-01, -4.7992163325114922e-17,
       3.0347528910975783e-33, 8.5989199506479701e-50),
  qd_real( 9.7003125319454397e-01, 1.8365300348428844e-17,
       -1.4311097571944918e-33, 8.5846781998740697e-51),
  qd_real( 9.6928123535654853e-01, -4.5663660261927896e-17,
       9.6147526917239387e-34, 8.1267605207871330e-51),
  qd_real( 9.6852209427441727e-01, 4.9475074918244771e-17,
       2.8558738351911241e-33, 6.2948422316507461e-50),
  qd_real( 9.6775383709347551e-01, -4.5512132825515820e-17,
       -1.4127617988719093e-33, -8.4620609089704578e-50),
  qd_real( 9.6697647104485207e-01, 3.8496228837337864e-17,
       -5.3881631542745647e-34, -3.5221863171458959e-50),
  qd_real( 9.6619000344541250e-01, 5.1298840401665493e-17,
       1.4564075904769808e-34, 1.0095973971377432e-50),
  qd_real( 9.6539444169768940e-01, -2.3745389918392156e-17,
       5.9221515590053862e-34, -3.8811192556231094e-50),
  qd_real( 9.6458979328981276e-01, -3.4189470735959786e-17,
       2.2982074155463522e-33, -4.5128791045607634e-50),
  qd_real( 9.6377606579543984e-01, 2.6463950561220029e-17,
       -2.9073234590199323e-36, -1.2938328629395601e-52),
  qd_real( 9.6295326687368388e-01, 8.9341960404313634e-18,
       -3.9071244661020126e-34, 1.6212091116847394e-50),
  qd_real( 9.6212140426904158e-01, 1.5236770453846305e-17,
       -1.3050173525597142e-33, 7.9016122394092666e-50),
  qd_real( 9.6128048581132064e-01, 2.0933955216674039e-18,
       1.0768607469015692e-34, -5.9453639304361774e-51),
  qd_real( 9.6043051941556579e-01, 2.4653904815317185e-17,
       -1.3792169410906322e-33, -4.7726598378506903e-51),
  qd_real( 9.5957151308198452e-01, 1.1000640085000957e-17,
       -4.2036030828223975e-34, 4.0023704842606573e-51),
  qd_real( 9.5870347489587160e-01, -4.3685014392372053e-17,
       2.2001800662729131e-33, -1.0553721324358075e-49),
  qd_real( 9.5782641302753291e-01, -1.7696710075371263e-17,
       1.9164034110382190e-34, 8.1489235071754813e-51),
  qd_real( 9.5694033573220882e-01, 4.0553869861875701e-17,
       -1.7147013364302149e-33, 2.5736745295329455e-50),
  qd_real( 9.5604525134999641e-01, 3.7705045279589067e-17,
       1.9678699997347571e-33, 8.5093177731230180e-50),
  qd_real( 9.5514116830577067e-01, 5.0088652955014668e-17,
       -2.6983181838059211e-33, 1.0102323575596493e-49),
  qd_real( 9.5422809510910567e-01, -3.7545901690626874e-17,
       1.4951619241257764e-33, -8.2717333151394973e-50),
  qd_real( 9.5330604035419386e-01, -2.5190738779919934e-17,
       -1.4272239821134379e-33, -4.6717286809283155e-50),
  qd_real( 9.5237501271976588e-01, -2.0269300462299272e-17,
       -1.0635956887246246e-33, -3.5514537666487619e-50),
  qd_real( 9.5143502096900834e-01, 3.1350584123266695e-17,
       -2.4824833452737813e-33, 9.5450335525380613e-51),
  qd_real( 9.5048607394948170e-01, 1.9410097562630436e-17,
       -8.1559393949816789e-34, -1.0501209720164562e-50),
  qd_real( 9.4952818059303667e-01, -7.5544151928043298e-18,
       -5.1260245024046686e-34, 1.8093643389040406e-50),
  qd_real( 9.4856134991573027e-01, 2.0668262262333232e-17,
       -5.9440730243667306e-34, 1.4268853111554300e-50),
  qd_real( 9.4758559101774109e-01, 4.3417993852125991e-17,
       -2.7728667889840373e-34, 5.5709160196519968e-51),
  qd_real( 9.4660091308328353e-01, 3.5056800210680730e-17,
       9.8578536940318117e-34, 6.6035911064585197e-50),
  qd_real( 9.4560732538052128e-01, 4.6019102478523738e-17,
       -6.2534384769452059e-34, 1.5758941215779961e-50),
  qd_real( 9.4460483726148026e-01, 8.8100545476641165e-18,
       5.2291695602757842e-34, -3.3487256018407123e-50),
  qd_real( 9.4359345816196039e-01, -2.4093127844404214e-17,
       1.0283279856803939e-34, -2.3398232614531355e-51),
  qd_real( 9.4257319760144687e-01, 1.3235564806436886e-17,
       -5.7048262885386911e-35, 3.9947050442753744e-51),
  qd_real( 9.4154406518302081e-01, -2.7896379547698341e-17,
       1.6273236356733898e-33, -5.3075944708471203e-51),
  qd_real( 9.4050607059326830e-01, 2.8610421567116268e-17,
       2.9261501147538827e-33, -2.6849867690896925e-50),
  qd_real( 9.3945922360218992e-01, -7.0152867943098655e-18,
       -5.6395693818011210e-34, 3.5568142678987651e-50),
  qd_real( 9.3840353406310806e-01, 5.4242545044795490e-17,
       -1.9039966607859759e-33, -1.5627792988341215e-49),
  qd_real( 9.3733901191257496e-01, -3.6570926284362776e-17,
       -1.1902940071273247e-33, -1.1215082331583223e-50),
  qd_real( 9.3626566717027826e-01, -1.3013766145497654e-17,
       5.2229870061990595e-34, -3.3972777075634108e-51),
  qd_real( 9.3518350993894761e-01, -3.2609395302485065e-17,
       -8.1813015218875245e-34, 5.5642140024928139e-50),
  qd_real( 9.3409255040425887e-01, 4.4662824360767511e-17,
       -2.5903243047396916e-33, 8.1505209004343043e-50),
  qd_real( 9.3299279883473885e-01, 4.2041415555384355e-17,
       9.0285896495521276e-34, 5.3019984977661259e-50),
  qd_real( 9.3188426558166815e-01, -4.0785944377318095e-17,
       1.7631450298754169e-33, 2.5776403305507453e-50),
  qd_real( 9.3076696107898371e-01, 1.9703775102838329e-17,
       6.5657908718278205e-34, -1.9480347966259524e-51),
  qd_real( 9.2964089584318121e-01, 5.1282530016864107e-17,
       2.3719739891916261e-34, -1.7230065426917127e-50),
  qd_real( 9.2850608047321559e-01, -2.3306639848485943e-17,
       -7.7799084333208503e-34, -5.8597558009300305e-50),
  qd_real( 9.2736252565040111e-01, -2.7677111692155437e-17,
       2.2110293450199576e-34, 2.0349190819680613e-50),
  qd_real( 9.2621024213831138e-01, -3.7303754586099054e-17,
       2.0464457809993405e-33, 1.3831799631231817e-49),
  qd_real( 9.2504924078267758e-01, 6.0529447412576159e-18,
       -8.8256517760278541e-35, 1.8285462122388328e-51),
  qd_real( 9.2387953251128674e-01, 1.7645047084336677e-17,
       -5.0442537321586818e-34, -4.0478677716823890e-50),
  qd_real( 9.2270112833387852e-01, 5.2963798918539814e-17,
       -5.7135699628876685e-34, 3.0163671797219087e-50),
  qd_real( 9.2151403934204190e-01, 4.1639843390684644e-17,
       1.1891485604702356e-33, 2.0862437594380324e-50),
  qd_real( 9.2031827670911059e-01, -2.7806888779036837e-17,
       2.7011013677071274e-33, 1.1998578792455499e-49),
  qd_real( 9.1911385169005777e-01, -2.6496484622344718e-17,
       6.5403604763461920e-34, -2.8997180201186078e-50),
  qd_real( 9.1790077562139050e-01, -3.9074579680849515e-17,
       2.3004636541490264e-33, 3.9851762744443107e-50),
  qd_real( 9.1667905992104270e-01, -4.1733978698287568e-17,
       1.2094444804381172e-33, 4.9356916826097816e-50),
  qd_real( 9.1544871608826783e-01, -1.3591056692900894e-17,
       5.9923027475594735e-34, 2.1403295925962879e-50),
  qd_real( 9.1420975570353069e-01, -3.6316182527814423e-17,
       -1.9438819777122554e-33, 2.8340679287728316e-50),
  qd_real( 9.1296219042839821e-01, -4.7932505228039469e-17,
       -1.7753551889428638e-33, 4.0607782903868160e-51),
  qd_real( 9.1170603200542988e-01, -2.6913273175034130e-17,
       -5.1928101916162528e-35, 1.1338175936090630e-51),
  qd_real( 9.1044129225806725e-01, -5.0433041673313820e-17,
       1.0938746257404305e-33, 9.5378272084170731e-51),
  qd_real( 9.0916798309052238e-01, -3.6878564091359894e-18,
       2.9951330310507693e-34, -1.2225666136919926e-50),
  qd_real( 9.0788611648766626e-01, -4.9459964301225840e-17,
       -1.6599682707075313e-33, -5.1925202712634716e-50),
  qd_real( 9.0659570451491533e-01, 3.0506718955442023e-17,
       -1.4478836557141204e-33, 1.8906373784448725e-50),
  qd_real( 9.0529675931811882e-01, -4.1153099826889901e-17,
       2.9859368705184223e-33, 5.1145293917439211e-50),
  qd_real( 9.0398929312344334e-01, -6.6097544687484308e-18,
       1.2728013034680357e-34, -4.3026097234014823e-51),
  qd_real( 9.0267331823725883e-01, -1.9250787033961483e-17,
       1.3242128993244527e-33, -5.2971030688703665e-50),
  qd_real( 9.0134884704602203e-01, -1.3524789367698682e-17,
       6.3605353115880091e-34, 3.6227400654573828e-50),
  qd_real( 9.0001589201616028e-01, -5.0639618050802273e-17,
       1.0783525384031576e-33, 2.8130016326515111e-50),
  qd_real( 8.9867446569395382e-01, 2.6316906461033013e-17,
       3.7003137047796840e-35, -2.3447719900465938e-51),
  qd_real( 8.9732458070541832e-01, -3.6396283314867290e-17,
       -2.3611649895474815e-33, 1.1837247047900082e-49),
  qd_real( 8.9596624975618511e-01, 4.9025099114811813e-17,
       -1.9440489814795326e-33, -1.7070486667767033e-49),
  qd_real( 8.9459948563138270e-01, -1.7516226396814919e-17,
       -1.3200670047246923e-33, -1.5953009884324695e-50),
  qd_real( 8.9322430119551532e-01, -4.1161239151908913e-18,
       2.5380253805715999e-34, 4.2849455510516192e-51),
  qd_real( 8.9184070939234272e-01, 4.6690228137124547e-18,
       1.6150254286841982e-34, -3.9617448820725012e-51),
  qd_real( 8.9044872324475788e-01, 1.1781931459051803e-17,
       -1.3346142209571930e-34, -9.4982373530733431e-51),
  qd_real( 8.8904835585466457e-01, -1.1164514966766675e-17,
       -3.4797636107798736e-34, -1.5605079997040631e-50),
  qd_real( 8.8763962040285393e-01, 1.2805091918587960e-17,
       3.9948742059584459e-35, 3.8940716325338136e-51),
  qd_real( 8.8622253014888064e-01, -6.7307369600274315e-18,
       1.2385593432917413e-34, 2.0364014759133320e-51),
  qd_real( 8.8479709843093779e-01, -9.4331469628972690e-18,
       -5.7106541478701439e-34, 1.8260134111907397e-50),
  qd_real( 8.8336333866573158e-01, 1.5822643380255127e-17,
       -7.8921320007588250e-34, -1.4782321016179836e-50),
  qd_real( 8.8192126434835505e-01, -1.9843248405890562e-17,
       -7.0412114007673834e-34, -1.0636770169389104e-50),
  qd_real( 8.8047088905216075e-01, 1.6311096602996350e-17,
       -5.7541360594724172e-34, -4.0128611862170021e-50),
  qd_real( 8.7901222642863353e-01, -4.7356837291118011e-17,
       1.4388771297975192e-33, -2.9085554304479134e-50),
  qd_real( 8.7754529020726124e-01, 5.0113311846499550e-17,
       2.8382769008739543e-34, 1.5550640393164140e-50),
  qd_real( 8.7607009419540660e-01, 5.8729024235147677e-18,
       2.7941144391738458e-34, -1.8536073846509828e-50),
  qd_real( 8.7458665227817611e-01, -5.7216617730397065e-19,
       -2.9705811503689596e-35, 8.7389593969796752e-52),
  qd_real( 8.7309497841829009e-01, 7.8424672990129903e-18,
       -4.8685015839797165e-34, -2.2815570587477527e-50),
  qd_real( 8.7159508665595109e-01, -5.5272998038551050e-17,
       -2.2104090204984907e-33, -9.7749763187643172e-50),
  qd_real( 8.7008699110871146e-01, -4.1888510868549968e-17,
       7.0900185861878415e-34, 3.7600251115157260e-50),
  qd_real( 8.6857070597134090e-01, 2.7192781689782903e-19,
       -1.6710140396932428e-35, -1.2625514734637969e-51),
  qd_real( 8.6704624551569265e-01, 3.0267859550930567e-18,
       -1.1559438782171572e-34, -5.3580556397808012e-52),
  qd_real( 8.6551362409056909e-01, -6.3723113549628899e-18,
       2.3725520321746832e-34, 1.5911880348395175e-50),
  qd_real( 8.6397285612158670e-01, 4.1486355957361607e-17,
       2.2709976932210266e-33, -8.1228385659479984e-50),
  qd_real( 8.6242395611104050e-01, 3.7008992527383130e-17,
       5.2128411542701573e-34, 2.6945600081026861e-50),
  qd_real( 8.6086693863776731e-01, -3.0050048898573656e-17,
       -8.8706183090892111e-34, 1.5005320558097301e-50),
  qd_real( 8.5930181835700836e-01, 4.2435655816850687e-17,
       7.6181814059912025e-34, -3.9592127850658708e-50),
  qd_real( 8.5772861000027212e-01, -4.8183447936336620e-17,
       -1.1044130517687532e-33, -8.7400233444645562e-50),
  qd_real( 8.5614732837519447e-01, 9.1806925616606261e-18,
       5.6328649785951470e-34, 2.3326646113217378e-51),
  qd_real( 8.5455798836540053e-01, -1.2991124236396092e-17,
       1.2893407722948080e-34, -3.6506925747583053e-52),
  qd_real( 8.5296060493036363e-01, 2.7152984251981370e-17,
       7.4336483283120719e-34, 4.2162417622350668e-50),
  qd_real( 8.5135519310526520e-01, -5.3279874446016209e-17,
       2.2281156380919942e-33, -4.0281886404138477e-50),
  qd_real( 8.4974176800085244e-01, 5.1812347659974015e-17,
       3.0810626087331275e-33, -2.5931308201994965e-50),
  qd_real( 8.4812034480329723e-01, 1.8762563415239981e-17,
       1.4048773307919617e-33, -2.4915221509958691e-50),
  qd_real( 8.4649093877405213e-01, -4.7969419958569345e-17,
       -2.7518267097886703e-33, -7.3518959727313350e-50),
  qd_real( 8.4485356524970712e-01, -4.3631360296879637e-17,
       -2.0307726853367547e-33, 4.3097229819851761e-50),
  qd_real( 8.4320823964184544e-01, 9.6536707005959077e-19,
       2.8995142431556364e-36, 9.6715076811480284e-53),
  qd_real( 8.4155497743689844e-01, -3.4095465391321557e-17,
       -8.4130208607579595e-34, -4.9447283960568686e-50),
  qd_real( 8.3989379419599952e-01, -1.6673694881511411e-17,
       -1.4759184141750289e-33, -7.5795098161914058e-50),
  qd_real( 8.3822470555483808e-01, -3.5560085052855026e-17,
       1.1689791577022643e-33, -5.8627347359723411e-50),
  qd_real( 8.3654772722351201e-01, -2.0899059027066533e-17,
       -9.8104097821002585e-35, -3.1609177868229853e-51),
  qd_real( 8.3486287498638001e-01, 4.6048430609159657e-17,
       -5.1827423265239912e-34, -7.0505343435504109e-51),
  qd_real( 8.3317016470191319e-01, 1.3275129507229764e-18,
       4.8589164115370863e-35, 4.5422281300506859e-51),
  qd_real( 8.3146961230254524e-01, 1.4073856984728024e-18,
       4.6951315383980830e-35, 5.1431906049905658e-51),
  qd_real( 8.2976123379452305e-01, -2.9349109376485597e-18,
       1.1496917934149818e-34, 3.5186665544980233e-51),
  qd_real( 8.2804504525775580e-01, -4.4196593225871532e-17,
       2.7967864855211251e-33, 1.0030777287393502e-49),
  qd_real( 8.2632106284566353e-01, -5.3957485453612902e-17,
       6.8976896130138550e-34, 3.8106164274199196e-50),
  qd_real( 8.2458930278502529e-01, -2.6512360488868275e-17,
       1.6916964350914386e-34, 6.7693974813562649e-51),
  qd_real( 8.2284978137582632e-01, 1.5193019034505495e-17,
       9.6890547246521685e-34, 5.6994562923653264e-50),
  qd_real( 8.2110251499110465e-01, 3.0715131609697682e-17,
       -1.7037168325855879e-33, -1.1149862443283853e-49),
  qd_real( 8.1934752007679701e-01, -4.8200736995191133e-17,
       -1.5574489646672781e-35, -9.5647853614522216e-53),
  qd_real( 8.1758481315158371e-01, -1.4883149812426772e-17,
       -7.8273262771298917e-34, 4.1332149161031594e-50),
  qd_real( 8.1581441080673378e-01, 8.2652693782130871e-18,
       -2.3028778135179471e-34, 1.5102071387249843e-50),
  qd_real( 8.1403632970594841e-01, -5.2127351877042624e-17,
       -1.9047670611316360e-33, -1.6937269585941507e-49),
  qd_real( 8.1225058658520388e-01, 3.1054545609214803e-17,
       2.2649541922707251e-34, -7.4221684154649405e-51),
  qd_real( 8.1045719825259477e-01, 2.3520367349840499e-17,
       -7.7530070904846341e-34, -7.2792616357197140e-50),
  qd_real( 8.0865618158817498e-01, 9.3251597879721674e-18,
       -7.1823301933068394e-34, 2.3925440846132106e-50),
  qd_real( 8.0684755354379922e-01, 4.9220603766095546e-17,
       2.9796016899903487e-33, 1.5220754223615788e-49),
  qd_real( 8.0503133114296355e-01, 5.1368289568212149e-17,
       6.3082807402256524e-34, 7.3277646085129827e-51),
  qd_real( 8.0320753148064494e-01, -3.3060609804814910e-17,
       -1.2242726252420433e-33, 2.8413673268630117e-50),
  qd_real( 8.0137617172314024e-01, -2.0958013413495834e-17,
       -4.3798162198006931e-34, 2.0235690497752515e-50),
  qd_real( 7.9953726910790501e-01, 2.0356723822005431e-17,
       -9.7448513696896360e-34, 5.3608109599696008e-52),
  qd_real( 7.9769084094339116e-01, -4.6730759884788944e-17,
       2.3075897077191757e-33, 3.1605567774640253e-51),
  qd_real( 7.9583690460888357e-01, -3.0062724851910721e-17,
       -2.2496210832042235e-33, -6.5881774117183040e-50),
  qd_real( 7.9397547755433717e-01, -7.4194631759921416e-18,
       2.4124341304631069e-34, -4.9956808616244972e-51),
  qd_real( 7.9210657730021239e-01, -3.7087850202326467e-17,
       -1.4874457267228264e-33, 2.9323097289153505e-50),
  qd_real( 7.9023022143731003e-01, 2.3056905954954492e-17,
       1.4481080533260193e-33, -7.6725237057203488e-50),
  qd_real( 7.8834642762660623e-01, 3.4396993154059708e-17,
       1.7710623746737170e-33, 1.7084159098417402e-49),
  qd_real( 7.8645521359908577e-01, -9.7841429939305265e-18,
       3.3906063272445472e-34, 5.7269505320382577e-51),
  qd_real( 7.8455659715557524e-01, -8.5627965423173476e-18,
       -2.1106834459001849e-34, -1.6890322182469603e-50),
  qd_real( 7.8265059616657573e-01, 9.0745866975808825e-18,
       6.7623847404278666e-34, -1.7173237731987271e-50),
  qd_real( 7.8073722857209449e-01, -9.9198782066678806e-18,
       -2.1265794012162715e-36, 3.0772165598957647e-54),
  qd_real( 7.7881651238147598e-01, -2.4891385579973807e-17,
       6.7665497024807980e-35, -6.5218594281701332e-52),
  qd_real( 7.7688846567323244e-01, 7.7418602570672864e-18,
       -5.9986517872157897e-34, 3.0566548232958972e-50),
  qd_real( 7.7495310659487393e-01, -5.2209083189826433e-17,
       -9.6653593393686612e-34, 3.7027750076562569e-50),
  qd_real( 7.7301045336273699e-01, -3.2565907033649772e-17,
       1.3860807251523929e-33, -3.9971329917586022e-50),
  qd_real( 7.7106052426181382e-01, -4.4558442347769265e-17,
       -2.9863565614083783e-33, -6.8795262083596236e-50),
  qd_real( 7.6910333764557959e-01, 5.1546455184564817e-17,
       2.6142829553524292e-33, -1.6199023632773298e-49),
  qd_real( 7.6713891193582040e-01, -1.8885903683750782e-17,
       -1.3659359331495433e-33, -2.2538834962921934e-50),
  qd_real( 7.6516726562245896e-01, -3.2707225612534598e-17,
       1.1177117747079528e-33, -3.7005182280175715e-50),
  qd_real( 7.6318841726338127e-01, 2.6314748416750748e-18,
       1.4048039063095910e-34, 8.9601886626630321e-52),
  qd_real( 7.6120238548426178e-01, 3.5315510881690551e-17,
       1.2833566381864357e-33, 8.6221435180890613e-50),
  qd_real( 7.5920918897838807e-01, -3.8558842175523123e-17,
       2.9720241208332759e-34, -1.2521388928220163e-50),
  qd_real( 7.5720884650648457e-01, -1.9909098777335502e-17,
       3.9409283266158482e-34, 2.0744254207802976e-50),
  qd_real( 7.5520137689653655e-01, -1.9402238001823017e-17,
       -3.7756206444727573e-34, -2.1212242308178287e-50),
  qd_real( 7.5318679904361252e-01, -3.7937789838736540e-17,
       -6.7009539920231559e-34, -6.7128562115050214e-51),
  qd_real( 7.5116513190968637e-01, 4.3499761158645868e-17,
       2.5227718971102212e-33, -6.5969709212757102e-50),
  qd_real( 7.4913639452345937e-01, -4.4729078447011889e-17,
       -2.4206025249983768e-33, 1.1336681351116422e-49),
  qd_real( 7.4710060598018013e-01, 1.1874824875965430e-17,
       2.1992523849833518e-34, 1.1025018564644483e-50),
  qd_real( 7.4505778544146595e-01, 1.5078686911877863e-17,
       8.0898987212942471e-34, 8.2677958765323532e-50),
  qd_real( 7.4300795213512172e-01, -2.5144629669719265e-17,
       7.1128989512526157e-34, 3.0181629077821220e-50),
  qd_real( 7.4095112535495911e-01, -1.4708616952297345e-17,
       -4.9550433827142032e-34, 3.1434132533735671e-50),
  qd_real( 7.3888732446061511e-01, 3.4324874808225091e-17,
       -1.3706639444717610e-33, -3.3520827530718938e-51),
  qd_real( 7.3681656887736990e-01, -2.8932468101656295e-17,
       -3.4649887126202378e-34, -1.8484474476291476e-50),
  qd_real( 7.3473887809596350e-01, -3.4507595976263941e-17,
       -2.3718000676666409e-33, -3.9696090387165402e-50),
  qd_real( 7.3265427167241282e-01, 1.8918673481573520e-17,
       -1.5123719544119886e-33, -9.7922152011625728e-51),
  qd_real( 7.3056276922782759e-01, -2.9689959904476928e-17,
       -1.1276871244239744e-33, -3.0531520961539007e-50),
  qd_real( 7.2846439044822520e-01, 1.1924642323370718e-19,
       5.9001892316611011e-36, 1.2178089069502704e-52),
  qd_real( 7.2635915508434601e-01, -3.1917502443460542e-17,
       7.7047912412039396e-34, 4.1455880160182123e-50),
  qd_real( 7.2424708295146689e-01, 2.9198471334403004e-17,
       2.3027324968739464e-33, -1.2928820533892183e-51),
  qd_real( 7.2212819392921535e-01, -2.3871262053452047e-17,
       1.0636125432862273e-33, -4.4598638837802517e-50),
  qd_real( 7.2000250796138165e-01, -2.5689658854462333e-17,
       -9.1492566948567925e-34, 4.4403780801267786e-50),
  qd_real( 7.1787004505573171e-01, 2.7006476062511453e-17,
       -2.2854956580215348e-34, 9.1726903890287867e-51),
  qd_real( 7.1573082528381871e-01, -5.1581018476410262e-17,
       -1.3736271349300259e-34, -1.2734611344111297e-50),
  qd_real( 7.1358486878079364e-01, -4.2342504403133584e-17,
       -4.2690366101617268e-34, -2.6352370883066522e-50),
  qd_real( 7.1143219574521643e-01, 7.9643298613856813e-18,
       2.9488239510721469e-34, 1.6985236437666356e-50),
  qd_real( 7.0927282643886569e-01, -3.7597359110245730e-17,
       1.0613125954645119e-34, 8.9465480185486032e-51),
  qd_real( 7.0710678118654757e-01, -4.8336466567264567e-17,
       2.0693376543497068e-33, 2.4677734957341755e-50)
};
#endif

/* Computes sin(a) and cos(a) using Taylor series.
   Assumes |a| <= pi/2048.                           */
void sincos_taylor(qd_real a,
                          inout qd_real sin_a, inout qd_real cos_a) {
  const float thresh = 0.5 * qd_eps * abs(to_float(a));
  qd_real p, s, t, x;

  if (is_zero(a)) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  x = -sqr(a);
  s = a;
  p = a;
  int i = 0;
  do {
    p *= x;
    t = p * inv_fact[i];
    s += t;
    i += 2;
  } while (i < n_inv_fact && abs(to_float(t)) > thresh);

  sin_a = s;
  cos_a = sqrt(1.0 - sqr(s));
}

qd_real sin_taylor(qd_real a) {
  const float thresh = 0.5 * qd_eps * abs(to_float(a));
  qd_real p, s, t, x;

  if (is_zero(a)) {
    return 0.0;
  }

  x = -sqr(a);
  s = a;
  p = a;
  int i = 0;
  do {
    p *= x;
    t = p * inv_fact[i];
    s += t;
    i += 2;
  } while (i < n_inv_fact && abs(to_float(t)) > thresh);

  return s;
}

qd_real cos_taylor(qd_real a) {
  const float thresh = 0.5 * qd_eps;
  qd_real p, s, t, x;

  if (is_zero(a)) {
    return 1.0;
  }

  x = -sqr(a);
  s = 1.0 + mul_pwr2(x, 0.5);
  p = x;
  int i = 1;
  do {
    p *= x;
    t = p * inv_fact[i];
    s += t;
    i += 2;
  } while (i < n_inv_fact && abs(to_float(t)) > thresh);

  return s;
}

qd_real sin(qd_real a) {

  /* Strategy.  To compute sin(x), we choose integers a, b so that

       x = s + a * (pi/2) + b * (pi/1024)

     and |s| <= pi/2048.  Using a precomputed table of
     sin(k pi / 1024) and cos(k pi / 1024), we can compute
     sin(x) from sin(s) and cos(s).  This greatly increases the
     convergence of the sine Taylor series.                          */

  if (is_zero(a)) {
    return 0.0;
  }

  // approximately reduce modulo 2*pi
  qd_real z = nint(a / qd_2pi);
  qd_real r = a - qd_2pi * z;

  // approximately reduce modulo pi/2 and then modulo pi/1024
  float q = floor(r.x[0] / qd_pi2[0] + 0.5);
  qd_real t = r - qd_pi2 * q;
  int j = int(q);
  q = floor(t.x[0] / _pi1024[0] + 0.5);
  t -= _pi1024 * q;
  int k = int(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    return qd_nan;
  }

  if (abs_k > 256) {
    return qd_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return sin_taylor(t);
      case 1:
        return cos_taylor(t);
      case -1:
        return -cos_taylor(t);
      default:
        return -sin_taylor(t);
    }
  }

  qd_real sin_t, cos_t;
  qd_real u = cos_table[abs_k-1];
  qd_real v = sin_table[abs_k-1];
  sincos_taylor(t, sin_t, cos_t);

  if (j == 0) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else {
      r = - u * cos_t - v * sin_t;
    }
  } else {
    if (k > 0) {
      r = - u * sin_t - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  }

  return r;
}

qd_real cos(qd_real a) {

  if (is_zero(a)) {
    return 1.0;
  }

  // approximately reduce modulo 2*pi
  qd_real z = nint(a / qd_2pi);
  qd_real r = a - qd_2pi * z;

  // approximately reduce modulo pi/2 and then modulo pi/1024
  float q = floor(r.x[0] / qd_pi2.x[0] + 0.5);
  qd_real t = r - qd_pi2 * q;
  int j = int(q);
  q = floor(t.x[0] / _pi1024.x[0] + 0.5);
  t -= _pi1024 * q;
  int k = int(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    return qd_nan;
  }

  if (abs_k > 256) {
    return qd_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return cos_taylor(t);
      case 1:
        return -sin_taylor(t);
      case -1:
        return sin_taylor(t);
      default:
        return -cos_taylor(t);
    }
  }

  qd_real sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);

  qd_real u = cos_table[abs_k-1];
  qd_real v = sin_table[abs_k-1];

  if (j == 0) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = - u * sin_t - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else {
      r = - u * cos_t - v * sin_t;
    }
  }

  return r;
}

void sincos(qd_real a, inout qd_real sin_a, inout qd_real cos_a) {

  if (is_zero(a)) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  // approximately reduce by 2*pi
  qd_real z = nint(a / qd_2pi);
  qd_real t = a - qd_2pi * z;

  // approximately reduce by pi/2 and then by pi/1024.
  float q = floor(t.x[0] / qd_pi2.x[0] + 0.5);
  t -= qd_pi2 * q;
  int j = int(q);
  q = floor(t.x[0] / _pi1024.x[0] + 0.5);
  t -= _pi1024 * q;
  int k = int(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    cos_a = sin_a = qd_nan;
    return;
  }

  if (abs_k > 256) {
    cos_a = sin_a = qd_nan;
    return;
  }

  qd_real sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);

  if (k == 0) {
    if (j == 0) {
      sin_a = sin_t;
      cos_a = cos_t;
    } else if (j == 1) {
      sin_a = cos_t;
      cos_a = -sin_t;
    } else if (j == -1) {
      sin_a = -cos_t;
      cos_a = sin_t;
    } else {
      sin_a = -sin_t;
      cos_a = -cos_t;
    }
    return;
  }

  qd_real u = cos_table[abs_k-1];
  qd_real v = sin_table[abs_k-1];

  if (j == 0) {
    if (k > 0) {
      sin_a = u * sin_t + v * cos_t;
      cos_a = u * cos_t - v * sin_t;
    } else {
      sin_a = u * sin_t - v * cos_t;
      cos_a = u * cos_t + v * sin_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      cos_a = - u * sin_t - v * cos_t;
      sin_a = u * cos_t - v * sin_t;
    } else {
      cos_a = v * cos_t - u * sin_t;
      sin_a = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      cos_a = u * sin_t + v * cos_t;
      sin_a =  v * sin_t - u * cos_t;
    } else {
      cos_a = u * sin_t - v * cos_t;
      sin_a = - u * cos_t - v * sin_t;
    }
  } else {
    if (k > 0) {
      sin_a = - u * sin_t - v * cos_t;
      cos_a = v * sin_t - u * cos_t;
    } else {
      sin_a = v * cos_t - u * sin_t;
      cos_a = - u * cos_t - v * sin_t;
    }
  }
}

qd_real atan(qd_real a) {
  return atan(a, qd_real(1.0));
}

qd_real atan(qd_real y, qd_real x) {
  /* Strategy: Instead of using Taylor series to compute
     arctan, we instead use Newton's iteration to solve
     the equation

        sin(z) = y/r    or    cos(z) = x/r

     where r = sqrt(x^2 + y^2).
     The iteration is given by

        z' = z + (y - sin(z)) / cos(z)          (for equation 1)
        z' = z - (x - cos(z)) / sin(z)          (for equation 2)

     Here, x and y are normalized so that x^2 + y^2 = 1.
     If |x| > |y|, then first iteration is used since the
     denominator is larger.  Otherwise, the second is used.
  */

  if (is_zero(x)) {

    if (is_zero(y)) {
      /* Both x and y is zero. */
      return qd_nan;
    }

    return (is_positive(y)) ? qd_pi2 : -qd_pi2;
  } else if (is_zero(y)) {
    return (is_positive(x)) ? qd_0 : qd_pi;
  }

  if (x == y) {
    return (is_positive(y)) ? qd_pi4 : -qd_3pi4;
  }

  if (x == -y) {
    return (is_positive(y)) ? qd_3pi4 : -qd_pi4;
  }

  qd_real r = sqrt(sqr(x) + sqr(y));
  qd_real xx = x / r;
  qd_real yy = y / r;

  /* Compute float precision approximation to atan. */
  qd_real z = atan(to_float(y), to_float(x));
  qd_real sin_z, cos_z;

  if (abs(xx.x[0]) > abs(yy.x[0])) {
    /* Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)  */
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
  } else {
    /* Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)  */
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
  }

  return z;
}


qd_real drem(qd_real a, qd_real b) {
  qd_real n = nint(a/b);
  return (a - n * b);
}

qd_real divrem(qd_real a, qd_real b, inout qd_real r) {
  qd_real n = nint(a/b);
  r = a - n * b;
  return n;
}

qd_real tan(qd_real a) {
  qd_real s, c;
  sincos(a, s, c);
  return s/c;
}

qd_real asin(qd_real a) {
  qd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    return qd_nan;
  }

  if (is_one(abs_a)) {
    return (is_positive(a)) ? qd_pi2 : -qd_pi2;
  }

  return atan(a, sqrt(1.0 - sqr(a)));
}

qd_real acos(qd_real a) {
  qd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    return qd_nan;
  }

  if (is_one(abs_a)) {
    return (is_positive(a)) ? qd_0 : qd_pi;
  }

  return atan(sqrt(1.0 - sqr(a)), a);
}

qd_real sinh(qd_real a) {
  if (is_zero(a)) {
    return 0.0;
  }

  if (abs(a) > 0.05) {
    qd_real ea = exp(a);
    return mul_pwr2(ea - inv(ea), 0.5);
  }

  /* Since a is small, using the above formula gives
     a lot of cancellation.   So use Taylor series. */
  qd_real s = a;
  qd_real t = a;
  qd_real r = sqr(t);
  float m = 1.0;
  float thresh = abs(to_float(a) * qd_eps);

  do {
    m += 2.0;
    t *= r;
    t /= (m-1) * m;

    s += t;
  } while (abs(t) > thresh);

  return s;
}

qd_real cosh(qd_real a) {
  if (is_zero(a)) {
    return 1.0;
  }

  qd_real ea = exp(a);
  return mul_pwr2(ea + inv(ea), 0.5);
}

qd_real tanh(qd_real a) {
  if (is_zero(a)) {
    return 0.0;
  }

  if (abs(to_float(a)) > 0.05) {
    qd_real ea = exp(a);
    qd_real inv_ea = inv(ea);
    return (ea - inv_ea) / (ea + inv_ea);
  } else {
    qd_real s, c;
    s = sinh(a);
    c = sqrt(1.0 + sqr(s));
    return s / c;
  }
}

void sincosh(qd_real a, inout qd_real s, inout qd_real c) {
  if (abs(to_float(a)) <= 0.05) {
    s = sinh(a);
    c = sqrt(1.0 + sqr(s));
  } else {
    qd_real ea = exp(a);
    qd_real inv_ea = inv(ea);
    s = mul_pwr2(ea - inv_ea, 0.5);
    c = mul_pwr2(ea + inv_ea, 0.5);
  }
}

qd_real asinh(qd_real a) {
  return log(a + sqrt(sqr(a) + 1.0));
}

qd_real acosh(qd_real a) {
  if (a < 1.0) {
    return qd_nan;
  }

  return log(a + sqrt(sqr(a) - 1.0));
}

qd_real atanh(qd_real a) {
  if (abs(a) >= 1.0) {
    return qd_nan;
  }

  return mul_pwr2(log((1.0 + a) / (1.0 - a)), 0.5);
}

qd_real fmod(qd_real a, qd_real b) {
  qd_real n = aint(a / b);
  return (a - b * n);
}

///=====================================================================
/// END qd-2.3.22+dfsg.1
///=====================================================================

#endif

const float nan = 0.0 / 0.0;
const float pi = 3.141592653;

float hypot1(float x, float y) { return sqrt(x * x + y * y); }
float hypot2(float x, float y) { return x * x + y * y; }

float srgb2lrgb(float s)
{
  if (s <= 0.04045)
    return s / 12.92;
  return pow((s + 0.055) / 1.055, 2.4);
}

vec3 srgb2lrgb(vec3 s)
{
  return vec3(srgb2lrgb(s.x), srgb2lrgb(s.y), srgb2lrgb(s.z));
}

float lrgb2srgb(float l)
{
  if (l <= 0.0031308)
    return l * 12.92;
  return 1.055  * pow(l, 1.0 / 2.4) - 0.055;
}

vec3 lrgb2srgb(vec3 s)
{
  return vec3(lrgb2srgb(s.x), lrgb2srgb(s.y), lrgb2srgb(s.z));
}

// <http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl>
vec3 hsv2rgb(vec3 c)
{
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  vec3 rgb = c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
  if (KFP_sRGB)
  {
    rgb = srgb2lrgb(rgb);
  }
  return rgb;
}

vec3 KF_Palette(float ix)
{
  // map [0..1) to [0..1024)
  ix -= floor(ix);
  ix *= 1024.0;
  // to match KF's regular colouring, need to sin-interpolate c0, c1
  // to get neighbouring colours in 1024-palette
  // and then linear-interpolate those (if smoothing is desired)
  // c0, c1 are neighbouring colours in n-palette with interpolant cf
  int m_nParts = textureSize(KFP_Palette, 0);
  vec3 c0, c1;
  {
    int i = int(floor(ix));
    float temp = float(i) * float(m_nParts) / 1024.0;
    int p = int(temp);
    int pn = (p + 1) % m_nParts;
    temp -= float(p);
    temp = sin((temp - 0.5) * pi) / 2.0 + 0.5;
    c0 = mix
      ( texelFetch(KFP_Palette, p, 0).rgb
      , texelFetch(KFP_Palette, pn, 0).rgb
      , temp
      );
  }
  if (KFP_Smooth)
  {
    int i = (int(floor(ix)) + 1) % 1024;
    float temp = float(i) * float(m_nParts) / 1024.0;
    int p = int(temp);
    int pn = (p + 1) % m_nParts;
    temp -= float(p);
    temp = sin((temp - 0.5) * pi) / 2.0 + 0.5;
    c1 = mix
      ( texelFetch(KFP_Palette, p, 0).rgb
      , texelFetch(KFP_Palette, pn, 0).rgb
      , temp
      );
    ix -= floor(ix);
    return mix(c0, c1, ix);
  }
  else
  {
    return c0;
  }
}

uint burtle_hash(uint a)
{
  a = (a+0x7ed55d16u) + (a<<12);
  a = (a^0xc761c23cu) ^ (a>>19);
  a = (a+0x165667b1u) + (a<<5);
  a = (a+0xd3a2646cu) ^ (a<<9);
  a = (a+0xfd7046c5u) + (a<<3);
  a = (a^0xb55a4f09u) ^ (a>>16);
  return a;
}

// uniform in [0,1)
float dither(int x, int y, int c)
{
  return float(burtle_hash(uint(x) + burtle_hash(uint(y) + burtle_hash(uint(c))))) / exp2(32.0);
}

vec2 GetPixelOffset(ivec2 ix)
{
  float x = 0.0, y = 0.0;
  int c = int(KFP_JitterSeed);
  if (c != 0)
  {
    float s = KFP_JitterScale;
    float u = dither(ix.x, ix.y, 2 * c + 0);
    float v = dither(ix.x, ix.y, 2 * c + 1);
    switch (KFP_JitterShape)
    {
      default:
      case 0: // uniform
        {
          x = s * (u - 0.5);
          y = s * (v - 0.5);
        }
        break;
      case 1: // Gaussian
        {
          // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
          float r = 0.0 < u && u < 1.0 ? sqrt(-2.0 * log(u)) : 0.0;
          float t = 2.0 * 3.141592653589793 * v;
          s *= 0.5;
          x = s * r * cos(t);
          y = s * r * sin(t);
        }
        break;
    }
  }
  else
  {
    x = 0.0;
    y = 0.0;
  }
  return vec2(x, y);
}

float49 cbrt(float49 a) { return nroot(a, 3); }
float49 wrap(float49 a) { return sub(a, floor(a)); }

float49 float49_(uint a)
{
  float hi = float(a);
  float lo = float(a - uint(hi));
  return float49_(hi, lo);
}
float49 float49_(uint a, uint b) { return add(ldexp(float49_(a), 32), float49_(b)); }
float49 float49_(uint a, uint b, float c) { return add(float49_(a, b), c); }

ivec2 Internal_TexCoord = Internal_TilePadding.yx + ivec2(Internal_TileSize.y - 1 - 2 * Internal_TilePadding.y - int(gl_FragCoord.y), int(gl_FragCoord.x));
ivec2 Internal_PixelCoord = Internal_TileOrigin.xy + ivec2(int(gl_FragCoord.x), Internal_TileSize.y - 1 - 2 * Internal_TilePadding.y - int(gl_FragCoord.y));

uint    getN1 (ivec2 offset) { return texelFetch(Internal_N1,  Internal_TexCoord + offset.yx, 0).r; }
uint    getN0 (ivec2 offset) { return texelFetch(Internal_N0,  Internal_TexCoord + offset.yx, 0).r; }
float   getNF (ivec2 offset) { return texelFetch(Internal_NF,  Internal_TexCoord + offset.yx, 0).r; }
float   getT  (ivec2 offset) { return texelFetch(Internal_T,   Internal_TexCoord + offset.yx, 0).r; }
float   getDEX(ivec2 offset) { return texelFetch(Internal_DEX, Internal_TexCoord + offset.yx, 0).r; }
float   getDEY(ivec2 offset) { return texelFetch(Internal_DEY, Internal_TexCoord + offset.yx, 0).r; }
uint    getN1 (void) { return getN1 (ivec2(0, 0)); }
uint    getN0 (void) { return getN0 (ivec2(0, 0)); }
float   getNF (void) { return getNF (ivec2(0, 0)); }
float   getT  (void) { return getT  (ivec2(0, 0)); }
float   getDEX(void) { return getDEX(ivec2(0, 0)); }
float   getDEY(void) { return getDEY(ivec2(0, 0)); }
bool    getGlitch  (ivec2 offset) { return getNF(offset) < 0.0; }
bool    getInterior(ivec2 offset) { return uvec2(getN0(offset), getN1(offset)) == KFP_Iterations; }
float49 getN  (ivec2 offset) { return float49_(getN1(offset), getN0(offset), 1.0 - getNF(offset)); }
vec2    getDE (ivec2 offset) { return vec2(getDEX(offset), getDEY(offset)); }
bool    getGlitch  (void) { return getGlitch(ivec2(0, 0)); }
bool    getInterior(void) { return getInterior(ivec2(0, 0)); }
float49 getN  (void) { return getN (ivec2(0, 0)); }
vec2    getDE (void) { return getDE(ivec2(0, 0)); }

bool inImage(ivec2 offset)
{
  ivec2 pixel = Internal_PixelCoord + offset;
  return 0 <= pixel.x && pixel.x < KFP_ImageSize.x &&
         0 <= pixel.y && pixel.y < KFP_ImageSize.y;
}

vec2 getJitter(ivec2 offset)
{
  if (inImage(offset))
  {
    return GetPixelOffset(Internal_PixelCoord + offset);
  }
  else
  {
    return vec2(0.0, 0.0); // FIXME reflect at image boundaries
  }
}

vec2 getJitter(void)
{
  return getJitter(ivec2(0, 0));
}

ivec2 getCoord(void)
{
  return Internal_TileOrigin.xy + ivec2(int(gl_FragCoord.x), int(gl_FragCoord.y));
}

float49 getN3x3(inout mat3 p, inout mat3 px, inout mat3 py)
{
  // load 3x3 stencil around the pixel
  float49 N = getN();
  for (int dj = -1; dj <= 1; ++dj)
  {
    for (int di = -1; di <= 1; ++di)
    {
      ivec2 offset = ivec2(di, dj);
      if (inImage(offset))
      {
        p[di + 1][dj + 1] = sub(getN(offset), N).x[0];
        vec2 delta = getJitter(offset);
        px[di + 1][dj + 1] = float(di) + delta.x;
        py[di + 1][dj + 1] = float(dj) + delta.y;
      }
    }
  }
  // reflect at image boundaries if necessary
  // this will break (result is infinite or NaN) for image size of 1 pixel
  p[1][1] *= 2.0;
  px[1][1] *= 2.0;
  py[1][1] *= 2.0;
  if (isnan(p[0][0])) {p[0][0] = p[1][1] - p[2][2];px[0][0] = px[1][1] - px[2][2];py[0][0] = py[1][1] - py[2][2];}
  if (isnan(p[0][1])) {p[0][1] = p[1][1] - p[2][1];px[0][1] = px[1][1] - px[2][1];py[0][1] = py[1][1] - py[2][1];}
  if (isnan(p[0][2])) {p[0][2] = p[1][1] - p[2][0];px[0][2] = px[1][1] - px[2][0];py[0][2] = py[1][1] - py[2][0];}
  if (isnan(p[1][0])) {p[1][0] = p[1][1] - p[1][2];px[1][0] = px[1][1] - px[1][2];py[1][0] = py[1][1] - py[1][2];}
  if (isnan(p[1][2])) {p[1][2] = p[1][1] - p[1][0];px[1][2] = px[1][1] - px[1][0];py[1][2] = py[1][1] - py[1][0];}
  if (isnan(p[2][0])) {p[2][0] = p[1][1] - p[0][2];px[2][0] = px[1][1] - px[0][2];py[2][0] = py[1][1] - py[0][2];}
  if (isnan(p[2][1])) {p[2][1] = p[1][1] - p[0][1];px[2][1] = px[1][1] - px[0][1];py[2][1] = py[1][1] - py[0][1];}
  if (isnan(p[2][2])) {p[2][2] = p[1][1] - p[0][0];px[2][2] = px[1][1] - px[0][0];py[2][2] = py[1][1] - py[0][0];}
  p[1][1] *= 0.5;
  px[1][1] *= 0.5;
  py[1][1] *= 0.5;
  return N;
}

float KF_Traditional(mat3 p, mat3 px, mat3 py)
{
  // traditional method reverse engineered from original code
  float gx = (p[0][1] - p[1][1]) * 1.414 / hypot1(px[0][1] - px[1][1], py[0][1] - py[1][1]);
  float gy = (p[1][0] - p[1][1]) * 1.414 / hypot1(px[1][0] - px[1][1], py[1][0] - py[1][1]);
  float gu = (p[0][0] - p[1][1]) * 1.414 / hypot1(px[0][0] - px[1][1], py[0][0] - py[1][1]);
  float gv = (p[0][2] - p[1][1]) * 1.414 / hypot1(px[0][2] - px[1][1], py[0][2] - py[1][1]);
  float g = abs(gx) + abs(gy) + abs(gu) + abs(gv);
  return 1.0 / g;
}

float KF_Traditional(void)
{
  mat3 p = mat3(nan, nan, nan, nan, nan, nan, nan, nan, nan);
  mat3 px = mat3(0.0);
  mat3 py = mat3(0.0);
  getN3x3(p, px, py);
  return KF_Traditional(p, px, py);
}

float KF_Forward3x3(mat3 p, mat3 px, mat3 py)
{
  // forward differencing in 8 directions from the target point
  float gx0 = sqr(p[0][1] - p[1][1]) / hypot2(px[0][1] - px[1][1], py[0][1] - py[1][1]);
  float gx2 = sqr(p[2][1] - p[1][1]) / hypot2(px[2][1] - px[1][1], py[2][1] - py[1][1]);
  float gy0 = sqr(p[1][0] - p[1][1]) / hypot2(px[1][0] - px[1][1], py[1][0] - py[1][1]);
  float gy2 = sqr(p[1][2] - p[1][1]) / hypot2(px[1][2] - px[1][1], py[1][2] - py[1][1]);
  float gu0 = sqr(p[0][0] - p[1][1]) / hypot2(px[0][0] - px[1][1], py[0][0] - py[1][1]);
  float gu2 = sqr(p[2][2] - p[1][1]) / hypot2(px[2][2] - px[1][1], py[2][2] - py[1][1]);
  float gv0 = sqr(p[2][0] - p[1][1]) / hypot2(px[2][0] - px[1][1], py[2][0] - py[1][1]);
  float gv2 = sqr(p[0][2] - p[1][1]) / hypot2(px[0][2] - px[1][1], py[0][2] - py[1][1]);
  float g = sqrt(0.25 * (gx0 + gx2 + gy0 + gy2 + gu0 + gu2 + gv0 + gv2));
  return 1.0 / (g * 2.8284271247461903);
}

float KF_Forward3x3(void)
{
  mat3 p = mat3(nan, nan, nan, nan, nan, nan, nan, nan, nan);
  mat3 px = mat3(0.0);
  mat3 py = mat3(0.0);
  getN3x3(p, px, py);
  return KF_Forward3x3(p, px, py);
}

float KF_Central3x3(mat3 p, mat3 px, mat3 py)
{
  // gerrit's central difference formula
  float gx = sqr(p[2][1] - p[0][1]) / hypot2(px[2][1] - px[0][1], py[2][1] - py[0][1]);
  float gy = sqr(p[1][2] - p[1][0]) / hypot2(px[1][2] - px[1][0], py[1][2] - py[1][0]);
  float g1 = sqr(p[2][2] - p[0][0]) / hypot2(px[2][2] - px[0][0], py[2][2] - py[0][0]);
  float g2 = sqr(p[0][2] - p[2][0]) / hypot2(px[0][2] - px[2][0], py[0][2] - py[2][0]);
  float g = sqrt(0.5 * (gx + gy + g1 + g2));
  return 1.0 / (g * 2.8284271247461903);
}

float KF_Central3x3(void)
{
  mat3 p = mat3(nan, nan, nan, nan, nan, nan, nan, nan, nan);
  mat3 px = mat3(0.0);
  mat3 py = mat3(0.0);
  getN3x3(p, px, py);
  return KF_Central3x3(p, px, py);
}

float KF_Diagonal2x2(mat3 p, mat3 px, mat3 py)
{
  // forward differencing in 2 diagonals of a 2x2 substencil
  if (KFP_JitterSeed == 0u)
  {
    float gu = sqr(p[0][0] - p[1][1]) / hypot2(px[0][0] - px[1][1], py[0][0] - py[1][1]);
    float gv = sqr(p[0][1] - p[1][0]) / hypot2(px[0][1] - px[1][0], py[0][1] - py[1][0]);
    float g = sqrt(gu + gv);
    return 1.0 / (g * 2.8284271247461903);
  }
  else
  {
    // with displacement correction by gerrit
    float nux = px[0][0] - px[1][1];
    float nuy = py[0][0] - py[1][1];
    float nvx = px[1][0] - px[0][1];
    float nvy = py[1][0] - py[0][1];
    float nu = hypot1(nux, nuy);
    float nv = hypot1(nvx, nvy);
    nux /= nu;
    nuy /= nu;
    nvx /= nv;
    nvy /= nv;
    float u = (p[0][0] - p[1][1]) / nu;
    float v = (p[1][0] - p[0][1]) / nv;
    float dotnunv = nux * nvx + nuy * nvy;
    float crossnunv = nux * nvy - nuy * nvx;
    float g = sqrt((u * u + v * v - 2 * u * v * dotnunv) / sqr(crossnunv));
    return 1.0 / (g * 2.8284271247461903);
  }
}

float KF_Diagonal2x2(void)
{
  mat3 p = mat3(nan, nan, nan, nan, nan, nan, nan, nan, nan);
  mat3 px = mat3(0.0);
  mat3 py = mat3(0.0);
  getN3x3(p, px, py);
  return KF_Diagonal2x2(p, px, py);
}

#if 0
            case Differences_LeastSquares3x3:
            {
              float gx = 0;
              float gy = 0;
              // compute_gradient_3x3(p, px, py, gx, gy); FIXME
              float g = hypot1(gx, gy);
              iter = float49_(g * 2.8284271247461903);
            }
            break;
            case Differences_LeastSquares2x2:
            {
              float gx = 0;
              float gy = 0;
              // compute_gradient_2x2(p, px, py, gx, gy); FIXME
              float g = hypot1(gx, gy);
              iter = float49_(g * 2.8284271247461903);
            }
            break;
#endif

float KF_Laplacian3x3(mat3 p, mat3 px, mat3 py)
{
  float L = 0.0;
  L +=  1.0 * p[0][0];
  L +=  4.0 * p[0][1];
  L +=  1.0 * p[0][2];
  L +=  4.0 * p[1][0];
  L -= 20.0 * p[1][1];
  L +=  4.0 * p[1][2];
  L +=  1.0 * p[2][0];
  L +=  4.0 * p[2][1];
  L +=  1.0 * p[2][2];
  L /=  6.0;
#define INV_LOG_2 1.4426950408889634
  float g = sqrt(abs(L * INV_LOG_2));
#undef INV_LOG_2
  return 1.0 / (g * 2.8284271247461903);
}

float KF_Laplacian3x3(void)
{
  mat3 p = mat3(nan, nan, nan, nan, nan, nan, nan, nan, nan);
  mat3 px = mat3(0.0);
  mat3 py = mat3(0.0);
  getN3x3(p, px, py);
  return KF_Laplacian3x3(p, px, py);
}

float KF_Analytic(ivec2 offset)
{
  return length(getDE(offset));
}

float KF_Analytic(void)
{
  return KF_Analytic(ivec2(0, 0));
}

float KF_DE(int method)
{
  switch (method)
  {
    default:
    case Differences_Traditional: return KF_Traditional();
    case Differences_Forward3x3: return KF_Forward3x3();
    case Differences_Central3x3: return KF_Central3x3();
    case Differences_Diagonal2x2: return KF_Diagonal2x2();
#if 0
    case Differences_LeastSquares2x2: return KF_LeastSquares2x2();
    case Differences_LeastSquares3x3: return KF_LeastSquares3x3();
#endif
    case Differences_Laplacian3x3: return KF_Laplacian3x3();
    case Differences_Analytic: return KF_Analytic();
  }
  return 0.0;
}

vec3 KF_InfiniteWaves(bool Smooth, float49 iter)
{
  float nH = 0.0, nS = 0.0, nB = 0.0;
  int nDR = 0, nDG = 0, nDB = 0;
  for (int i = 0; i < KFP_MultiWavesCount; i++){
    float nPeriod = float(KFP_MultiWaves[i].x);
    float g;
    if (Smooth)
      g = sin(2.0 * pi * wrap(div(iter, nPeriod * 2.0)).x[0]) / 2.0 + 0.5;
    else
      g = sin(2.0 * pi * wrap(div(floor(iter), nPeriod * 2.0)).x[0]) / 2.0 + 0.5;
    if (nPeriod < 0)
      g = -nPeriod / 100.0;
    int type = KFP_MultiWaves[i].z;
    if (type == 0){
      nH += g;
      nDR++;
    }
    if (type == 1){
      nS += g;
      nDG++;
    }
    if (type == 2){
      nB += g;
      nDB++;
    }
  }
  if (nDR > 0)
    nH /= float(nDR);
  if (nDG > 0)
    nS /= float(nDG);
  if (nDB > 0)
    nB /= float(nDB);
  return hsv2rgb(vec3(nH, nS, nB));
}

vec2 KF_TextureWarp(float TexturePower, float TextureRatio, vec2 SlopeDir)
{
  mat3 p = mat3(nan, nan, nan, nan, nan, nan, nan, nan, nan);
  mat3 px = mat3(0.0);
  mat3 py = mat3(0.0);
  getN3x3(p, px, py);
  float nImgOffs = TexturePower / 64.0;
  float diffx = p[0][1] - p[1][1];
  float diffy = p[1][0] - p[1][1];
  float diff = dot(vec2(diffx, diffy), SlopeDir);
  diff  = 1.0 + diff;
  diffx = 1.0 + diffx;
  diffy = 1.0 + diffy;
  diff  = pow(diff,  TexturePower);
  diffx = pow(diffx, TexturePower);
  diffy = pow(diffy, TexturePower);
  float sx = 1.0;
  float  sy = 1.0;
  if (diff  <= 1.0) { diff  = 1.0 / diff; }
  if (diffx <= 1.0) { diffx = 1.0 / diffx; sx = -sx; }
  if (diffy <= 1.0) { diffy = 1.0 / diffy; sy = -sy; }
  diff  = (atan(diff)  - pi / 4.0) / (pi / 4.0);
  diffx = (atan(diffx) - pi / 4.0) / (pi / 4.0);
  diffy = (atan(diffy) - pi / 4.0) / (pi / 4.0);
  diff  *= TextureRatio / 100.0;
  diffx *= TextureRatio / 100.0;
  diffy *= TextureRatio / 100.0;
  float dx = nImgOffs + sx * TexturePower * diffx;
  float dy = nImgOffs - sy * TexturePower * diffy;
  return vec2(dx, dy);
}

vec4 KF_Slopes(bool analytic, vec2 SlopeDir, float Power, float Ratio)
{
  Power *= float(KFP_ImageSize.x) / 640.0;
  Ratio /= 100.0;
  vec2 vdiff;
  if (analytic)
  {
    vec2 DE = getDE();
    vdiff = vec2(1.0, -1.0) * DE / dot(DE, DE);
  }
  else
  {
    float49 N = getN();
    vdiff.x = -sub(getN(ivec2(1, 0)), N).x[0];
    vdiff.y = sub(getN(ivec2(0, -1)), N).x[0];
  }
  float diff = dot(vdiff, SlopeDir);
  diff *= Power;
  if (diff >= 0.0)
  {
    diff = atan(diff) / (pi / 2.0);
    diff = diff * Ratio;
    return vec4(vec3(0.0), diff);
  }
  else
  {
    diff = -diff;
    diff = atan(diff) / (pi / 2.0);
    diff = diff * Ratio;
    return vec4(vec3(1.0), diff);
  }
}

float49 KF_InverseTransition(float49 iter)
{
  float49 iter_ = floor(iter);
  float offs = sub(iter, iter_).x[0];
  iter = add(iter_, 1.0 - offs);
  return iter;
}

float49 KF_IterTransform(float49 iter0)
{
  float49 iter = iter0;
  switch (KFP_ColorMethod)
  {
    default:
    case ColorMethod_Standard: iter = iter; break;
    case ColorMethod_SquareRoot: iter = sqrt(max(0.0, iter)); break;
    case ColorMethod_CubicRoot: iter = cbrt(max(0.0, iter)); break;
    case ColorMethod_Logarithm: iter = log(max(1.0, iter)); break;
    case ColorMethod_LogLog: iter = log(add(1.0, log(add(1.0, iter)))); break;
    case ColorMethod_ATan: iter = atan(iter); break;
    case ColorMethod_FourthRoot: iter = sqrt(sqrt(max(0.0, iter))); break;
    case ColorMethod_Stretched:
    {
      float49 imin = float49_(KFP_IterationsMin[1], KFP_IterationsMin[0]);
      float49 imax = float49_(KFP_IterationsMax[1], KFP_IterationsMax[0]);
      iter = mul(1024.0, div(sub(iter, imin), sub(imax, imin))); break;
    }
    case ColorMethod_DistanceLinear:
    case ColorMethod_DEPlusStandard:
    case ColorMethod_DistanceLog:
    case ColorMethod_DistanceSqrt:
    {
      iter = float49_(1.0 / KF_DE(KFP_Differences));
      // post differencing transfer functions
      iter = mul(iter, float(KFP_ImageSize.x) / 640.0);
      if (KFP_ColorMethod == ColorMethod_DistanceSqrt || KFP_ColorMethod == ColorMethod_DEPlusStandard)
        iter = sqrt(max(0.0, iter));
      else if (KFP_ColorMethod == ColorMethod_DistanceLog)
        iter = log(max(1.0, add(1.0, iter)));
      if(gt(iter, 1024.0))
        iter = float49_(1024.0);
      if(KFP_ColorMethod == ColorMethod_DEPlusStandard && gt(iter, KFP_IterDiv))
        iter = iter0;
      break;
    }
  }
  iter = div(iter, KFP_IterDiv);
  iter = add(iter, KFP_ColorOffset);
  if (KFP_InverseTransition)
  {
    iter = KF_InverseTransition(iter);
  }
  return iter;
}

vec3 KF_Colour(void)
{
  vec3 s = vec3(0.0);
  if (! KFP_ShowGlitches && getGlitch())
  {
    discard;
  }
  if (getInterior())
  {
    s = KFP_InteriorColor;
  }
  else
  {
    float49 iter = float49_(getN1(), getN0(), KFP_Flat ? 0.0 : 1.0 - getNF());
    iter = KF_IterTransform(iter);
    if (KFP_MultiWavesEnabled)
    {
      vec3 nRGB = KF_InfiniteWaves(KFP_Smooth, iter);
      if (KFP_MultiWavesBlend)
      {
        iter = add(iter, KFP_PhaseColorStrength / 100.0 * 1024.0 * getT());
        if (KFP_InverseTransition)
        {
          iter = KF_InverseTransition(iter);
        }
        vec3 nRGB2 = KF_Palette(wrap(div(iter, 1024.0)).x[0]);
        s = mix(nRGB, nRGB2, 0.5);
      }
      else
      {
        s = nRGB;
      }
    }
    else
    {
      iter = add(iter, KFP_PhaseColorStrength / 100.0 * 1024.0 * getT());
      if (KFP_InverseTransition)
      {
        iter = KF_InverseTransition(iter);
      }
      s = KF_Palette(wrap(div(iter, 1024.0)).x[0]);
    }
  }
  if (KFP_TextureEnabled)
  {
    vec2 tc = getCoord() + KF_TextureWarp(KFP_TexturePower, KFP_TextureRatio, KFP_SlopeDir);
    tc /= vec2(KFP_ImageSize.xy);
    s = mix(s, texture(KFP_Texture, tc).rgb, KFP_TextureMerge);
  }
  if (KFP_Slopes)
  {
    vec4 slope = KF_Slopes(KFP_Differences == Differences_Analytic, KFP_SlopeDir, KFP_SlopePower, KFP_SlopeRatio);
    s = mix(s, slope.rgb, slope.a);
  }
  if (getInterior() && ! KFP_TextureEnabled)
  {
    s = KFP_InteriorColor;
  }
  return s;
}

void main(void)
{
  Internal_One = Internal_Zero + KFP_ImageSize.x / KFP_ImageSize.x;
  Internal_Colour = vec4(clamp(colour(), vec3(0.0), vec3(65504.0)), 1.0);
}

#line 0 1
