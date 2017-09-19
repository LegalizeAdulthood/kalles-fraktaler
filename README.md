Kalles Fraktaler 2 + GMP
========================

As the orginal upstream author Karl Runmo says:

> Want to create DEEP Mandelbrot fractals 100 times faster than the commercial
> programs, for FREE? One hour or one minute? Three months or one day?
> Try Kalles Fraktaler!

I (Claude Heiland-Allen) forked the code and swapped out the custom arbitrary
precision floating point code for the highly optimized GMP library, making it
even faster.

Cross-compiled to Windows from Linux MINGW64, using GMP.  Now with many other
enhancements (mostly speed optimisations and bugfixes).

Original upstream version:

- <http://www.chillheimer.de/kallesfraktaler/>

This version:

- <https://mathr.co.uk/kf/kf.html>

Feedback:

- <https://fractalforums.org/kalles-fraktaler/> new forum (still in beta)
- <http://www.fractalforums.com/kalles-fraktaler/> legacy forum
- <mailto:claude@mathr.co.uk?subject=Kalles%20Fraktaler%202> personal mail


Known Bugs
----------

- "no newton.kfr" blank image on load and newton-raphson zoom fails with bad
  period detected (reported by Kalles Fraktaler)
- newton-raphson zooming to minibrot doesn't increase maxiters enough sometimes
- opencl support is very broken, proof of concept only
- may be difficult to build the source at the moment (out of date instructions)


Differences From Upstream 2.11.1
--------------------------------

### Incompatible Changes

- **In version `kf-2.12.1` and above**, DE colouring method #5 is once again
  backwards compatible with upstream `2.11.1`.  Parameter files made with
  `2.11.1+gmp.DATE` versions should be modified to use Distance (Square Root)
  colouring method #8.

- **In version `kf-2.11.1+gmp.20170822` only**, DE colouring method #5 used log
  instead of sqrt for a more perceptually linear effect.  In later versions,
  this log scaling is achieved with a new colouring method #7, while the DE
  colouring method #5 reverts to sqrt as before.  The new colouring method ID
  allows old `2.11.1+gmp.DATE`
  parameter files to be loaded into current versions and display as
  intended.  Any parameter files saved with the new Distance (Logarithm)
  colouring method will not display as intended in older versions.  Parameter
  files using Distance colouring method saved with this particular version
  should be modified to use Distance (Logarithm) in the latest version.

### Other Changes

- Makefile build system using MINGW to cross-compile to Windows from Linux
- uses GMP for arbitrary precision floating point instead of custom code
- uses Boost wrapper around GMP floats for higher-level coding
- used JPEG library downloaded as necessary at build time, instead of bundled
- long double support built into EXE (no separate DLL needed)
- virtually unlimited precision (memory needed for precise numbers is an issue)
- threaded calculations reimplemented with barriers to avoid WINE slowdown
- workaround for WINE issue artificially limiting image size (up to 2GiB now)
- bugfix: inflection performance issue (was converting number types needlessly)
- bugfix: cross-hair resource issue (reported and fixed by Kalles Fraktaler)
- miscellaneous code cleanups (-fpermissive fixes, const fixes, delete[] fixes,
  64bit compatibility paranoia)
- formula inner loops generated at compile time from high level specification
  XML using XSLT and a preprocessor implemented in Haskell
- optimized some reference calculations by floating temporaries out of loops
- optimized Newton-Raphson zooming by using lower-level GMP calls
- very experimental and broken OpenCL using CLEW (still disabled at build time)


Change Log
----------

- **kf-2.12.1** (2017-09-19)

    - simplified version numbering;
    - built for 64bit (as before) and 32bit (new);
    - documentation improvements;
    - fix division by zero assertion failure in File -> Examine zoom sequence;
    - fix crash in File -> Examine zoom sequence with only 1 image file;
    - adjust distance colour modes for backwards compatibility;

- **kf-2.11.1+gmp.20170913**

    - revert incompatible de log vs sqrt colouring change, instead add a new
      Distance (Logarithm) colouring method #7;
    - documentation improvements;
    - limit maximum series approximation terms to 60 to try to fix overskipping
      with large images

- **kf-2.11.1+gmp.20170822**

    - bugfix preprocessor for abs() formulas
    - de colouring with log instead of sqrt

- **kf-2.11.1+gmp.20170820**

    - bugfix preprocessor for diffabs() formulas

- **kf-2.11.1+gmp.20170714**

    - disabled OpenCL (be more compatible)

- **kf-2.11.1+gmp.20170713**

    - optimized Newton-Raphson zooming (3x faster in one test)

- **kf-2.11.1+gmp.20170711**

    - workaround for WINE issue artificially limiting image size (now bitmaps
      up to 2GiB can be created on all platforms)

- **kf-2.11.1+gmp.20170710**

    - optimized formulas (reference calculation for quadratic Mandlebrot is much
      faster due to lower-level calls to gmp)
    - very experimental opencl support (mostly broken)
    - bugfixes (fix hang loading deep zoom locations, fix newton size in new
      view radius calculation, more complete library credits in documentation)
    - prune dead code (incomplete jpeg library deleted from source, complete
      version downloaded at build time as needed, delete rudimentary openmp
      support, delete non-performant barrier variant, delete slower-than-gmp
      mpfr support, delete custom floating point support)

- **kf-2.11.1+gmp.20170703**

    - formulas now generated at compile time from formula definition XML using
      XSL stylesheet
    - used fixed format floats instead of scientific
    - try to hide command prompt window on Windows

- **kf-2.11.1+gmp.20170508**

    - restored threaded reference calculations (reimplemented with barrier()
      semantics to avoid single-threaded WINE SetEvent() rendezvous)

- **kf-2.11.1+gmp.20170504**

    - removed threaded reference calculations (too much overhead)
    - miscellaneous code cleanups (no need for -fpermissive, const fixes,
      delete[] fixes, 64bit compatibility paranoia)

- **kf-2.11.1+gmp.20170406**

    - fixed precision bugs (easy deep zoom, interactive failure)
    - fixed performance bug with inflections
    - fixed cross-hair resource bug
    - added WINDRES argument to build system
    - added more info to about dialog
    - include source code with release

- **kf-2.11.1+gmp.20170330.1**

    - fixes a crasher bug in the previous version

- **kf-2.11.1+gmp.20170330**

    - unlimited precision
    - separate compilation

- **kf-2.11.1+gmp.20170313**

    - long double compiled into exe (no dll)

- **kf-2.11.1+gmp.20170307**

    - kf-2.11.1 + gmp

- **kf-2.9.3+gmp.20170307**

    - kf-2.9.3 + gmp


TODO
----

- building: document the current system requirements
- user interface: batch mode
- user interface: PNG image export (JPEG is 8bit YUV which means colour gamut
  and precision is lost, even before lossy compression artifacts...)
- user interface: scripting interface
- calculations: implement scaled long double for e4900 to e9800
- calculations: optimize series approximation and probe point stuff
- calculations: work on OpenCL some more (try to get it working)
- preprocessor: float out temporaries from reference iterations
- preprocessor: flatten complex numbers to separate real and imaginary parts
- preprocessor: automatically parallelize reference iterations
- colouring: assume sRGB display and gamma-correct downscaling
- colouring: load/save palette to/from image (PNG required)
- colouring: rework entirely (now: 1024 colours with mandatory interpolation)
- colouring: implement Pauldelbrot's multiwave colouring


Getting The Code
----------------

I distribute EXEs bundled together with the corresponding source code.

The latest source code is available from my git repository:

    git clone https://code.mathr.co.uk/kalles-fraktaler-2.git
    cd kalles-fraktaler-2
    git checkout master       # for Karl's original upstream
    git checkout claude       # for MINGW build system and bug fixes
    git checkout claude-gmp   # for the GMP fork
    git checkout formulas     # for current development
    git tag -l                # list available release tags


Building On Linux
-----------------

(note: these instructions are out of date)

Build instructions for cross-compiling from GNU/Linux require about 3.5GB of
disk space and good internet download speed (or patience). About 410MB of
downloads after the chroot debootstrap step. If you have recent Debian you can
skip the chroot step and install natively.

0. Setup Debian Stretch chroot:

        mkdir ./vm
        sudo debootstrap stretch ./vm/
        sudo mount proc ./vm/proc -t proc
        sudo mount sysfs ./vm/sys -t sysfs
        sudo cp /etc/hosts ./vm/etc/hosts
        sudo chroot ./vm /bin/bash
        cd

1. Install dependencies (inside the chroot if you made one):

        apt-get install \
          build-essential \
          cabal-install \
          ghc \
          git \
          libghc-parsec3-dev \
          libtool \
          lzip \
          m4 \
          mingw-w64 \
          p7zip \
          wine64 \
          wine-binfmt \
          xsltproc \
          zip

2. Prepare non-root build user:

        adduser build
        # enter and confirm password
        su - build
        export CPPFLAGS=-D__USE_MINGW_ANSI_STDIO
        mkdir -p ~/win64/src
        # mkdir -p ~/win32/src

3. Download sources:

    Download the latest Boost (which is at time of writing is 1.65.1) and
    latest GMP (currently version 6.1.2) and clone kf git sources:

        cd ~/win64/src
        # cd ~/win32/src
        wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.7z
        wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.lz
        git clone https://code.mathr.co.uk/kalles-fraktaler-2.git

    Internet access is no longer required after this step.

4. Build GMP

        cd ~/win64/src
        tar xf gmp-6.1.2.tar.lz
        cd gmp-6.1.2
        ./configure --host=x86_64-w64-mingw32 --prefix=$HOME/win64
        # ./configure --host=i686-w64-mingw32 --prefix=$HOME/win32
        make -j 8
        make install
        make check

5. Prepare Boost headers

        cd ~/win64/src
        # cd ~/win32/src
        7zr x boost*.7z
        cd ~/win64/include
        # cd ~/win32/include
        ln -s ../src/boost*/boost/

6. Finally, build Kalles Fraktaler 2 + GMP

        cd ~/win64/src
        cd kalles-fraktaler-2
        git checkout formulas
        make -j 8 SYSTEM=64  # or SYSTEM=32 for 32bit version FIXME incomplete
        ./kf.exe  # test to see if it works

7. To cut a release bundle, use the script

        export VERSION=2.whatever
        git tag -s kf-${VERSION}
        ./release.sh ${VERSION}


Building on Windows
-------------------

(note: these instructions are out of date)

Build instructions for compiling on Windows (thanks to knighty!):

0. Remove any old msys2.

1. Downloaded latest version of msys2 (msys2-x86_64-20161025.exe).
  This is the 64 bit version. msys2-i686-20161025.exe is the 32 bit version.

2. After running it, it installs msys2. At the end the msys2 shell is launched.

3. In the msys2 shell, invoke pacman:

        pacman -Syuu

    This have to be done until is says there is nothing to do anymore.

4. Close the msys2 shell:

        exit

5. Reopen msys2 shell (from startup menu).

6. Install mingw/gcc 64 bit:

        pacman -S mingw-w64-x86_64-toolchain

    one can also install 32 bit version by:

        pacman -S mingw-w64-i686-toolchain

7. Install Boost

        pacman -S mingw-w64-x86_64-boost

    from msys shell

8. Close msys2 shell then open "msys2 mingw 64 bit" shell (in order to have all
  the environment variables properly set)

9. Change directory to the kalles fraktaler sources (where `Makefile` resides).

10. Compile

        mingw32-make WINDRES=windres

    (if this doesn't work edit the Makefile to replace the line

        WINDRES ?= x86_64-w64-mingw32-windres

    to

        WINDRES ?= windres

    and run `mingw32-make` without arguments)

11. Execute it this way from (msys2 mingw 64 bit) command line:

        ./fraktal_sft64    # for the claude branch
        ./kf.exe           # for the claude-gmp branch

    because it is linked dynamically to some libraries. In order to execute it
    from the explorer one needs to copy `libgmp-10.dll` and
    `libwinpthread-1.dll` from `msys64/mingw64/bin` next to the generated
    executable.


Configuration (`COMPILE_FLAGS` in `Makefile`)
---------------------------------------------

- add `-DKF_THREADED_REFERENCE_EVENT` to use original threaded reference
  calculations (too much overhead in WINE to make it worthwhile, except at very
  deep zooms)
- add `-DKF_THREADED_REFERENCE_BARRIER` to use barrier() threaded reference
  (acceptable overhead in WINE, CPU affinity is adjusted with zoom depth,
  enabled by default)


Legal
-----

- Copyright (c) 2013-2017 Karl Runmo, (c) 2017 Claude Heiland-Allen
- this software is based in part on the work of the Independent JPEG Group
- the GMP library is used under the conditions of the GNU Lesser General Public
  License version 3 and the GNU General Public License version 2
- the Boost library is used under the Boost Software License Version 1.0
- the CLEW library is used under the Boost Software License Version 1.0

NOTE: the binaries are statically linked with GMP, which is under dual LGPLv3 /
GPLv2 license. If you redistribute the binaries you must also be prepared to
distribute the source corresponding to those binaries to anyone you distribute
the binary to. To make this easier for you, the more recent zips include the
source too (though you'll also need to get the Boost and GMP sources). And of
course insert here the usual legal disclaimers about NO WARRANTY OF ANY KIND.
