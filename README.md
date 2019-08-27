# Kalles Fraktaler 2 +

As the orginal upstream author Karl Runmo says:

> Want to create DEEP Mandelbrot fractals 100 times faster than the commercial
> programs, for FREE? One hour or one minute? Three months or one day?
> Try Kalles Fraktaler!

I (Claude Heiland-Allen) forked the code and swapped out the custom arbitrary
precision floating point code for highly optimized libraries, making it
even faster.  Cross-compiled to Windows from Linux MINGW64.  Now with many other
enhancements.

Original upstream version:

- <http://www.chillheimer.de/kallesfraktaler/>

This version:

- <https://mathr.co.uk/kf/kf.html>

Feedback:

- <https://fractalforums.org/kalles-fraktaler/15> current forum
- <http://www.fractalforums.com/kalles-fraktaler/> legacy forum (read only)
- <mailto:claude@mathr.co.uk?subject=Kalles%20Fraktaler%202> personal mail


## Quick Start

- Download the latest zip from the website:
  <https://mathr.co.uk/kf/kf.html#download>

- Unzip it wherever you want, no installation required.

- Launch `kf.64.exe` for 64-bit (to be preferred), `kf.32.exe` for 32-bit.

- Start exploring!

  - Use the mouse scroll wheel to zoom.

  - Choose different fractal formulas in the iterations dialog (Ctrl-I).

  - Choose different colours in the colouring dialog (Ctrl-C).

  - Zoom deeper more quickly with Newton-Raphson zooming (Ctrl-D).


## Known Bugs

- Windows Defender sometimes reports that KF contains malware
  (false positive, I hope - I check with virustotal.com and all seems fine)
  (possibly due to screenshots for crosshair window, some code also moves
  the mouse cursor)
- out of memory conditions cause crashes (for example, if bitmap creation
  fails - also need to check huge sizes) (reported by gerrit)
- resizing window during examine zoom sequence auto solve glitches leads to
  corruption of the zoom sequence data
- "stop autosolve" during examine zoom sequence fails and corrupts zoom sequence
- minimizing window during zoom sequence rendering corrupts image (saves blank
  image or repeated frame) (reported by gerrit and CFJH)
- zoom size 1.2 sequence JPEGs have "boxes" while KFB is ok (reported by saka)
- resizing window by dragging frame corner in WINE on Linux sometimes crashes
- with "reuse reference", corrupt image at transition between number types
  (eg e600) (reported by CFJH) - workaround is to render in segments or force
  the number type higher ("use long double always", "use floatexp always")
- crash when zooming too quickly near interior black regions (reported by
  Foxxie) "usually near the elephant valley area or seahorse valley area of
  minibrots, happens worse the faster you zoom, usually if you try to zoom at
  the skinniest part very quickly"
- on special locations kf renders endless references and comes to no end
  (reported by CFJH)
- analytic DE broken with some power 3 Mandelbrot locations (reported by gerrit)
  (workaround is to disable series approximation)
- "resume zoom sequence" assumes "zoom size" is an integer
- "resume zoom sequence" re-uses last set zoom count limit
- "examine zoom sequence" doesn't save corrected PNG images during glitch solve
- speckles when rendering zoom out sequence
- black regions when rendering zoom out sequence (maximum iterations are reduced
  too much before spirals appear in next frame) (reported by gerrit)
- there is still a race conditions in guessing (doesn't wait for previous
  progressive interlacing pass to be 100% done before the next one starts)
- newton-raphson zooming to minibrot doesn't increase maxiters enough sometimes
- newton-raphson zoom preset depths are bad for formulas with power other than 2
- burning ship series approximation stops at first fold (typically 1 period of
  a central miniship), could potentially subdivide the region (and shift the
  series, for the folded part) and carry on?
- burning ship series approximation probe points might miss some folds near
  the edges of the image (need to increase probe point density?)
- bad combinations of skew, distance estimation, and series approximation
- scaled (long) double yr,yi can underflow to 0, eventually causing derivatives
  to be too small and de overflows to infinity -> blank screen: workaround is to
  force long double or floatexp as appropriate
- auto skew (escape) button doesn't work well with some formulas (eg SimonBrot)
- NR zoom doesn't work well in skewed locations
- nanomb1/2 OrderM, OrderN can only be changed by hand-editing .kfs Settings files
- nanomb2 RadiusScale can only be changed by hand-editing .kfs Settings files
- nanomb1/2 number type fixed to floatexp (long double or double would be faster)
- nanomb1/2 number type is not rescaled (only matters for long double / double)
- nanomb2 currently disables glitch detection and correction
- nanomb1/2 reference calculations are not multithreaded (single core only)
- nanomb1/2 reference calculations are using slow Boost C++ wrapper for MPFR
- kf-tile.exe doesn't support skew yet
- status bar reference count doesn't reset when zooming before it is "Done"
- help button in file browser does nothing
- opencl support is very broken, proof of concept only
- may be difficult to build the source at the moment
  (out of date instructions for Windows, dependency on 'et', ...)


## Differences From Upstream 2.11.1

### Incompatible Changes

- **In versions `2.14.7` and above**, the new "LogLog" and "ATan" iteration
  count transfer functions (ColorMethod 9 and 10) are not available in
  earlier versions.

- **In versions `2.14.5` and above**, the new "Low bailout" (SmoothMethod: 2)
  is not available in earlier versions.

- **In versions `2.14.5` and above**, the two new General Abs Quadratic formulas
  are not available in earlier versions.

- **In versions `2.14.4` and above**, the new "Mothbrot" formulas are not
  available in earlier versions.

- **In versions `2.14.3` and above**, colouring of zoomed out views is fixed,
  and may be different from earlier versions.

- **In versions `2.14.1` and above**, there are two new formulas which are
  not available in earlier versions.

- **In versions `2.13.1` through `2.13.9`**, derivatives are calculated by
  default.  In versions `2.14.1` and above, you need to enable the derivatives
  calculations (in the Iterations dialog, or via a settings file) if you
  want to use analytic distance estimation colouring.  Derivatives are not
  available in earlier versions.

- **In versions `2.13.1` through `2.13.4` inclusive** the interior is
  white.  In `2.13.5` and above it is user-settable, defaulting to black,
  Parameter files made with earlier `2.13` versions should be modified to
  explicitly set a white interior colour.  In earlier versions (including
  upstream `2.11.1`) the interior is black with no way of changing it.

- **In version `kf-2.12.10` only** the jitter is "Gaussian" with no way of
  changing the shape.  In `2.12.11` and above, the shape can be changed, and
  the default is now "uniform".  Earlier versions (including upstream `2.11.1`)
  do not support jitter at all.  To get the same results as `2.12.10` you must
  enable the Gaussian jitter checkbox and set the jitter scale box to 1.

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

- multiple finite difference methods for distance colouring
  (only the default `Traditional` is available in `2.11.1`).

- analytic DE colouring with derivatives (not available in `2.12.x` or earlier).


### Other Changes

- Makefile build system using MINGW to cross-compile to Windows from Linux
- uses GMP for arbitrary precision floating point instead of custom code
- uses Boost wrapper around GMP floats for higher-level coding
- use installed JPEG library, instead of bundled sources
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
- XML preprocessor optimizes more reference calculations in the same way
- optimized Newton-Raphson zooming by using lower-level GMP calls
- very experimental and broken OpenCL using CLEW (still disabled at build time)
- save images to PNG format as well as JPEG
- colouring uses floating point internally (fewer quantisation steps)
- dithering at end of colouring to improve perceptual quality (reduced banding)
- "glitch low tolerance" checkbox that can be enabled to detect glitches more
  aggressively (disabled by default, enable it if you get undetected glitches)
- updated program icon with transparent background and large version
- parameter data is saved as comment in image files (both PNG and JPEG)
- preferences (rendering settings not related to location) save and load
  (.kfs files and PNG/JPEG comments too)
- command line arguments to load settings and/or location
- command line arguments to render and save PNG/JPEG/KFB Map before quiting
- see change log below for more...


## Change Log

- **kf-2.14.7** (????-??-??)

    - feature: save as EXR (combined image with iteration data)
    - feature: LogLog and ATan color methods (iteration count transfer
      functions) (ported from Mandel Machine, suggested by Dinkydau)
    - feature: option to reduce thread count by 1 (suggested by saka)
    - bugfix: set thread background mode for improved responsiveness
      (reported by saka)
    - bugfix: remove annoying warning dialogs from store zoom out sequence
      (reported by saka)
    - bugfix: size status bar sections proportionally to window width
      (fixes accessibility issues with non-standard system font size)
      (reported by fractal)
    - bugfix: make Octave/Matlab scripts into function files for broader
      compatibility (usage example: `stratify("tile",4)`) (bug reported
      by Chronicler1701, fix suggested by gerrit).
    - new dependency: mingw-std-threads git
    - new dependency: ilmbase 2.3.0
    - new dependency: openexr 2.3.0
    - upgrade to boost 1.71.0
    - upgrade to gsl 2.6

- **kf-2.14.6.1** (2019-07-19)

    - bugfix: load palette after map in command line rendering so that the
      colours from the palette are actually used (reported by blob)

- **kf-2.14.6** (2019-07-18)

    - new feature: high quality image downscaling with anti-aliasing using
      the Pixman library
    - new feature: ball method for period detection in Newton-Raphson zoom
    - new feature: command line KFB map colouring
    - bugfix: image is refreshed from first calculated pixels, instead of 1%
      (reported by Fluoroantimonic_Acid)
    - bugfix: don't read Iteration Divisor from KFB map files (for historical
      reasons this is an integer in the file format, instead of the floating
      point value that it can take now in the rest of the program)
    - upgrade to boost 1.70.0
    - upgrade to mpfr 4.0.2p1
    - upgrade to png 1.6.37
    - upgrade to glm 0.9.9.5
    - upgrade to mingw64 g++ 8.3-win32 20190406 and rebuild everything

- **kf-2.14.5** (2019-03-20)

    - new feature: command line zoom sequence rendering with `-z NFRAMES` flag
    - new feature: low bailout (SmoothMethod: 2) (suggested by Dinkydau)
      <https://fractalforums.org/f/28/t/2594/msg12801#msg12801>
    - new feature: General Abs Quadratic formulas (inspired by Alef and gerrit)
      <https://fractalforums.org/f/28/t/2653/msg13242#msg13242>
    - bugfix: thin black line in exterior between iteration bands (reported by CFJH)
    - bugfix: malformed location no longer crashes (reported by Foxxie)
    - upgrade to mpfr 4.0.2
    - upgrade to glm 0.9.9.4

- **kf-2.14.4** (2019-01-18)

    - new feature: nanomb1 algorithm for power 2 Mandelbrot (experimental)
      (originally by knighty) (not quite as experimental as nanomb2, suggested
      by gerrit, requires exact period to be set by Newton or location dialog)
    - new feature: make interior checking for nanomb optional (new menu item)
    - new feature: RadiusScale setting for nanomb2 (only settable in .kfs so
      far, default value 0.1 may be safe but slow)
    - new feature: Mothbrot family of formulas (suggested by simon.snake)
      (not optimized yet, preliminary version for testing)
    - bugfix: incremental drawing looks much better with slopes and/or nde
      (there is a crosshair visible while rendering, which disappears when
      the image is finished, not sure why this happens)
    - upgrade to libpng 1.6.36
    - upgrade to libtiff 4.0.10
    - upgrade to libmpfr 4.0.1-p14
    - upgrade to boost 1.69.0

- **kf-2.14.3** (2018-10-29)

    - new feature: nanomb2 algorithm for power 2 Mandelbrot (experimental)
      (originally by knighty)
    - new feature: display Newton period (and set limit for nanomb2 algorithm)
      in Location dialog
    - new feature: Ctrl-Shift-W shortcut to set image size (suggested by gerrit)
    - bugfix: distinguish dialog titles for set window size and set image size
    - bugfix: fix corrupt images when zooming out from the default view
    - new feature: dialogs to en/disable derivatives when opening file with/out
      analytic DE; automatically enabling derivatives if needed in CLI mode
      (reported by gerrit)
    - default zoom size changed from 4 to 2 (press Ctrl-4 to reset to 4, or load
      a settings file with your preference) (suggested by gerrit, to avoid the
      dreaded "zoom size is not 2" warning when saving zoom out sequence).
    - merged changes from 2.13.11

- **kf-2.13.11** (2018-10-29)

    - bugfix: make preprocessor compatible with ghc-8.6 (MonadFail)
    - bugfix: fix translation with reuse reference enabled
      (reported by Dinkydau)
    - bugfix: don't add additional references if autosolve is disabled
    - upgrade to libpng 1.6.35

- **kf-2.14.2** (2018-10-23)

    - bugfix: initial window size was off by a few pixels
    - new feature: undo/redo (menu and Ctrl-Z, Ctrl-Y) (suggested by TwinDragon)
    - new feature: added copy/paste to menu (Ctrl-X, Ctrl-V worked already)
    - new feature: Ctrl-Shift-K shortcut to reset skew
    - new feature: Ctrl-Shift-T shortcut to reset rotation
    - bugfix: minor improvements to documentation
    - bugfix: fix series approximation derivative initialization
    - internal: refactor series approximation for Mandelbrot power > 2
    - merged changes from 2.13.10

- **kf-2.13.10** (2018-10-23)

    - bugfix: corrupt image at transition between number types (eg near e600)
      (reported by CFJH)
    - bugfix: changing "threads per CPU" during rendering crashes
      (reported by CFJH) (the menu is now disabled during rendering)
    - bugfix: set approx terms to 3 for other than power 2 Mandelbrot
    - bugfix: fix floatexp toString for negative values and 0
    - internal: update to et-kf API version 5

- **kf-2.14.1** (2018-09-27)

    - new feature: series approximation for Burning Ship power 2, using
      case analysis to know how it gets reflected at each iteration
    - two new formulas: General Quadratic Minus, General Quadratic Plus;
      "RedShiftRider factor A" sets coefficients (suggested by gerrit)
    - new feature: "Fast" preset activated by default
    - new feature: GUI menus reorganized to make more sense and updated
      the manual to reflect this
    - new feature: "check for updates" option in the file menu (requires
      internet access)
    - new feature: 1/4 and 1/2 threads per CPU settings (minimum 1 thread)
    - new feature: optional use of d/dZ in autoskew escape (suggested by
      gerrit)

- **kf-2.13.9** (2018-09-06)

    - new feature: auto skew without miniset: new button in Newton-Raphson
      zooming dialog, "Auto Skew (Escape)", that skews at the center pixels
      (algorithm suggested by gerrit)
    - bugfix: set window title on parameter drag and drop (reported by gerrit)
    - bugfix: reading PNG metadata works even if it is moved after the image
      data chunks and has a miscapitalized "Comment" keyword
    - bugfix: fix Mandelbar derivative calculations for ADE (reported by gerrit)
    - bugfix: fix 4th False Quasi formulas Newton-Raphson zooming

- **kf-2.13.8** (2018-08-28)

    - new feature: auto skew (via Newton-Raphson zooming dialog)
    - new feature: show/hide crosshair window (suggested by gerrit and others)
    - new feature: quality presets ("fast" but inaccurate, "best" but slow)
    - new feature: page up / page down keyboard shortcuts to zoom in / out
    - new feature: drag-and-drop parameter files from the file manager to
      the main window to open them
    - new feature: drag-and-drop palette files from the file manager to
      the color dialog window to open them
    - new feature: CLI prints total remaining pixels (suggested by gerrit)
    - bugfix: suppress error dialogs when loading metadata from TIFF
    - bugfix: remove obsolete .ini stuff that was overriding default settings
    - bugfix: enabled "no reuse center" by default
      (without it zoom out sequence sometimes glitches)
    - major documentation improvements
    - new dependency: GLM 0.9.9.0

- **kf-2.13.7** (2018-08-14)

    - export to uncompressed TIFF images (faster than PNG, but larger)
      (suggested by gerrit)
    - nicer default palette
    - upgrade to Boost 1.68.0
    - upgrade to MPFR 4.0.1-p13

- **kf-2.13.6** (2018-08-02)

    - fix for "a single dot appears in the middle of a mini"
      (reported by gerrit)
    - fix for "zooming out saves Zoom 0 in KFR parameter file"
      (reported by CFJH)
    - "open map" function (file menu)
      workflow: open KFR, ESC to abort, open KFB, adjust colours, save image
      (ignoring the warning about a KFB map existing with the same name)

- **kf-2.13.5** (2018-06-26)

    - colouring is parallelized for speedup when multiple cores are available
    - only sort pixel indices in interactive mode (speeds up batch mode)
    - palette loading fixed to load all colouring data
    - interior colour can be set now, defaults to black
    - upgrade to GSL 2.5

- **kf-2.13.4** (2018-06-21)

    - derivative computations can be switched off in the iterations
      dialog (increases rendering speed, breaks analytic DE colouring)

- **kf-2.12.13.1** (2018-06-04)

    - fix severe performance bug introduced in previous version
    - upgrade to Boost 1.67.0
    - upgrade to MPFR 4.0.1-p6

- **kf-2.12.13** (2018-05-15)

    - fix crasher bug relating to series approximation memory allocation
      (reported by gerrit)

- **kf-2.13.3** (2018-05-08)

    - Newton-Raphson zooming precision loss bug fix (reported by gerrit)
    - Newton-Raphson zooming for Mandelbrot powers 6 7 8 9 10
      (reported by gerrit)
    - Newton-Raphson zooming size estimate fixed for more formulas (now
      correctly uses smallest power instead of largest, for renormalization
      behaviour near zero)
    - Mandelbrot power 6 7 8 9 10 analytic distance estimation bug fixes
      (reported by gerrit)

- **kf-2.13.2** (2018-05-01)

    - Newton-Raphson zooming for all formulas (using code generated by 'et')
      (no progress reporting yet, nor convergence detection for early exit)
    - R^2 matrix derivatives for all formulas (fixes analytic DE with skew)
    - fix image corruption when enabling guessing with analytic DE
    - fix assertion failed popup (race condition between parameter changes and
      colouring during rendering) (reported by gerrit)

- **kf-2.13.1** (2018-03-30)

    - added derivative calculations
    - known bug: calculated even if not needed for colouring
    - known bug: some derivatives are C approximations instead of R^2 matrix
    - analytic DE colouring mode
    - known bug: DE data not saved to KFB map files
    - interior is white instead of black
    - known bug: interior colour should be customizable
    - series approximation uses Horner's rule to evaluate polynomials
      (minor optimisation)
    - default settings are best-quality by default

- **kf-2.12.12** (2018-03-30)

    - Free Software under GNU Affero General Public License, version 3 or greater
    - least squares numerical differences (experimental; subject to change)
    - Laplacian numerical differences (suggested by gerrit)
    - (possibly stratified) tiled settings/parameters generator program
    - stratified tiling reassembler GNU Octave script

- **kf-2.12.11** (2018-03-12)

    - uniform jitter by default (suggested by gerrit)
    - Gaussian jitter available as an option
    - jitter scale factor setting (suggested by gerrit)
    - log verbosity flag for command line rendering
    - upgrade to libjpeg6b2 <https://jpegclub.org/support/>
    - fix bug where Ctrl-S would overwrite the wrong file (reported by gerrit)
    - fix bug where lines would appear aligned with references (reported by
      gerrit and CFJH)
    - fix "color method and infinite waves are immediately applied" even when
      colouring is disabled (reported by gerrit)

- **kf-2.12.10** (2018-02-22)

    - upgrade to MPFR 4.0.1
    - removed some build instructions from README (see the prepare.sh script)
      (it was too annoying to have to keep updating everything in two places)
    - command line non-interactive mode works without opening a window
      (suggested by gerrit)
    - pixel jitter (dithering of position): set non-zero seed in iterations
      dialog to enable it, different seeds give slightly different images which
      can be stacked later for anti-aliasing as an alternative to rendering
      large and downscaling later; helps reduce Moiré grid artifacts at the
      cost of noise
    - "Enable" toggle in colouring dialog: colouring takes long for large
      images, disable the toggle if you want to make many changes at once, then
      click "Apply" or re-enable to see the results (discussed with gerrit)

- **kf-2.12.9.1** (2018-01-24)

    - fixed broken "approx low tolerance" checkbox (thanks to Dinkydau)
    - fixed NR-zoom dialog to remember custom zoom factor between reopenings
      (reported by gerrit)
    - fixed initialization of some menu item state on load

- **kf-2.12.9** (2018-01-24)

    - improved Newton-Raphson zooming dialog, now zooms to a factor between
      current zoom and minibrot size (suggested by Dinkydau)
    - implemented "scaled long double" iterations for Mandelbrot power 2 and 3,
      which might speed up some zooms between e4900 and e9800
    - fixed "the maxiterations suddenly resets to some very low value"
      (reported by gerrit)
    - fixed "Newton iterates 100 times before failing", now fails fast if
      the C goes outside the target (reported by gerrit and Dinkydau, fix
      suggested by knighty)
    - fixed "Newton zooms to blank image", reported by gerrit, fix was to
      remove a check on the derivative being huge that was falsely reporting
      success
    - fixed "manual ApproxTerms gets reset to AutoApproxTerms value"
      (reported by gerrit)
    - barrier implementation yields if there are fewer CPUs than threads
      (fixes very slow Newton-Raphson zooming and reference calculations
      on dual-core laptop) (this is the behaviour of 2.12.5 and earlier, but
      conditional on number of cores)
    - upgrade to libpng 1.6.34
    - added shell script containing some of the build instructions from
      the documentation

- **kf-2.12.8** (2017-12-22)

    - fixed the fix for hang in normalisation (was generating bad images)
      (reported by gerrit with bug7.png)
    - restore ignoring isolated glitches
    - fixed DE spots bug (somehow the reference was being reset to the center
      of the screen sometimes without its pixel position being updated)
      (reported by gerrit)
    - fixed bugs with guessing and glitch status (only interpolate when the
      glitch status of both neighbours is the same) (caused lines when dragging
      the view, probably other badness too)
    - added logic to prevent too-huge image size being entered via the GUI
      (settings files are still not checked...)
    - clarified store zoom count user interface (reported by CFJH)
    - fixed store zoom out for deep zooms (bug reported by gerrit) (the issue
      was exceptions thrown by std::stod(), remember to catch them)
    - refactor auto iterations (issue reported by gerrit remains unresolved)
    - fixed hardcoded count in glitch correction
    - use Taylor intervals in ball-period method (code copied from knighty)

- **kf-2.12.7** (2017-12-07)

    - renabled guessing conditional on menu option (was disabled in 2.12.4 as
      enabling it made some random speckles, possibly due to a race condition)
      (requested by Kalles Fraktaler and Fractal universe)
    - added "threads per cpu core" setting
    - compiles clean with -Wwrite-strings
    - even lower resolution preview for more intensive locations (suggested by
      Foxxie) (implemented with Adam7-style interlacing with circular sorting)
    - copy (Ctrl-X) and paste (Ctrl-V) parameters from the system clipboard
    - fixed hang crash bug when normalizing smooth iteration values
    - use interval arithmetic ball-period method instead of box-period (speeds
      up Newton-Raphson zooming a bit) (idea from knighty and gerrit)

- **kf-2.12.6** (2017-11-24)

    - fix central differences (reported by gerrit)
    - fix insufficient precision in Zoom: saved in .kfr (reported by CFJH)
    - option to render zoom out sequence without saving KFB maps (suggested by
      CFJH)
    - option to stop rendering zoom out sequence after a certain number of
      frames (suggested by CFJH)
    - resume zoom sequence works without KFB maps saved every frame (still needs
      a "last.kfb", this is saved automatically when needed)
    - barrier no longer yields (fixes priority inversion on heavily loaded
      systems) (reported by gerrit)
    - switch from GMP mpf_t to MPFR (fixes some blank images on load, also some
      Newton-Raphson zoom failures - bug involved incorrect normalization)
      (reported by Kalles Fraktaler and gerrit)

- **kf-2.12.5** (2017-11-02)

    - preferences (rendering settings not related to location) save and load
      (.kfs files and PNG/JPEG comments too)
    - command line arguments to load settings and/or location
    - command line arguments to render and save PNG/JPEG/KFB Map before quiting
    - auto-added new references recalculate only all glitched pixels (in earlier
      versions it would recalculate all pixels with same integer iteration
      count, which may or may not have been glitched, and may have missed some
      glitches)
    - glitch correction now uses glitch flag instead of just iteration count
      (this ensures the reference is added in a really glitched pixel, so at
      least one pixel will be fixed by each reference, ensuring termination
      with a finite number of references)
    - fix bugs with references when calculating their own pixels
    - fix off-by-one when references fix their own pixels (reported by gerrit)
    - single pixel glitches are no longer fixed by copying neighbour
    - fixed glitch at image boundary correction
    - fixed memory leak in glitch correction
    - fixed horizontal line corruption in examine zoom sequence glitch
      correction (reported by Fractal universe)
    - fixed "Mandelbar Celtic" formula (reported by Kalles Fraktaler)
    - fixed "Mandelbar" formula (reported by Foxxie)
    - fixed "Burning Ship Power 4" formula (reported by Foxxie)
    - fixed complex formulas reference precision problems (reported by Foxxie)
    - fix for auto-iterations (now respects GUI) (reported by Foxxie)
    - fix for crash selecting invalid power (reported by Foxxie)
    - fix for gigantic zoom value bug in Newton zooming (reported by Foxxie and
      gerrit)
    - added "no reuse center" option to prevent rectangle pasting (suggested by
      quaz0r)
    - fix for signed integer overflow reports negative percentage in status bar
      when the number of iterations is large (reported by Foxxie and another)
    - fix for confusing PNG save options dialog called "JPEG properties"
    - major code refactoring into multiple files for ease of maintenance
    - delete no-longer-used single-threaded Newton-Raphson zooming code

- **kf-2.12.4** (2017-10-06)

    - "glitch low tolerance" checkbox that can be enabled to detect glitches
      more aggressively (enabling it allows "Olbaid-ST-023.kfr" to render
      correctly, but taking 16x longer than with it disabled) (incorrect
      render reported by Kalles Fraktaler)
    - updated program icon with transparent background and large version,
      and use it for child windows too
    - long double off-by-one bug fixed (incorrect render reported by CFJH)
    - floatexp implementation bug fixes (0.0 was implemented incorrectly)
      (caused a series approximation underskipping failure that was reported
      by Kalles Fraktaler)
    - fix distance colouring grid artifacts with small zoom size by disabling
      "reuse center" when zoom size is not an integer (reported by gerrit)
    - fix some iteration band edge artifacts (with external postprocessing of
      kfb map files) by increasing the smooth bailout radius from 100 to 10000
      (reported by gerrit)
    - ensure added references fix their corresponding pixel (suggested by
      Pauldelbrot) (prevents loop in auto solve glitches whereby a reference
      was repeatedly being added at the same location without progress being
      made)
    - disabled "guessing" (was causing occasional randomly bright single pixels
      at low zoom levels, possibly a race condition?)
    - build against an installed libjpeg instead of each time after clean
    - fix broken complex formula reference calculations
    - fix compilation warnings (now almost clean with -Wall -Wextra)
    - compile as C++17 (should also still work as C++11 and C++14)
    - delete bitrotten code paths for SetEvent()-based multithreading
    - delete unused thread affinity setting code
    - parameter data is saved as comment in image files (both PNG and JPEG)
    - can load parameters from image file comments (both PNG and JPEG)
    - delete obsolete VS build system
    - fix "infinite waves" colouring (reported by gerrit)

- **kf-2.12.3** (2017-09-25)

    - multiple finite difference methods for distance colouring
      (suggested by gerrit)
    - fix bug in examine zoom sequence auto solve glitches (first frame only)
      (reported by Dinkydau and Fractal universe)
    - raise limit for maximum number of references from 199 to 10000
      (default still 69, you can change it in iterations dialog)
    - settable number of references per pass for examine zoom sequence auto
      solve glitches (default still 10)

- **kf-2.12.2** (2017-09-20)

    - PNG image saving support using libpng and zlib;
    - JPEG default quality to 100 (was 99);
    - colouring uses floating point internally to reduce quantisation steps;
    - dithering at end of colouring to improve perceptual quality;
    - `formula.cpp` included in source zip so GHC is not needed unless changing
      formula code;
    - optimized `diffabs()` code: one test Burning Ship location is 7.5% faster;
    - preprocessor optimizes reference calculations by floating temporary
      variable (re)allocations out of the inner loops: one test Burning Ship
      location is 30% faster;

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


## TODO

### User Interface

- crosshair cursor with more contrast (suggested by CFJH)
- show box-region for nr-zoom before clicking (suggested by
  Foxxie)
- adjust the size of the box via slider or like shift_scroll
  wheel or something like that? (suggested by Foxxie for nr-zoom, could also
  be useful for ctrl-left-click zoom)
- undo history for calculation data (suggested by TwinDragon)
- online help within program (suggested by TwinDragon)
- save image now function (without waiting for calculations)
- command line: print total runtime (suggested by gerrit)
- log window for diagnostics/debugging

### Calculations

- increase ref count limit without restarting from scratch
- increase maxiters limit without restarting from scratch
- optimize series approximation and probe point stuff
- work on OpenCL some more (try to get it working)
- calculate series approximation in parallel with reference
- refine minibrot using interior distance estimates
- refine minibrot using boundary shrinking (calculate edges only)
- lowest-connected-bailout radius (suggested by Dinkydau)
- use minimum |z| pixels for new references (suggested by quaz0r)
- option to build with int64_t iteration counts (audit -Wconversion)

### Newton-Raphson Zooming

- add zoom depth offset in addition to zoom depth factor to (eg) zoom
  closer to minibrot or embedded Julia
- zoom to Misiurewicz points (custom zoom factor, manual preperiod
  selection) (suggested by gerrit)
- properly debug huge zoom values from size estimate
- make it work better in hard-skewed locations (need to skew the box
  period coordinates?)
- add disk save/resume functionality for long-running calculations in
  case of system snafu

### NanoMB

- automatic number type selection
- automatic number type promotion on overflow
- re-scaled number types
- SSA escape radius factor control (around 0.1 to 1.0)
- automatic period detection

### Preprocessor

- flatten complex numbers to separate real and imaginary parts
- common subexpression elimination (share results, might be
  especially useful for large powers of complex numbers)
- automatically parallelize reference iterations

### Colouring

- high bit depth image export (eg 16bit PNG) (suggested by Dinkydau)
- assume sRGB display and gamma-correct downscaling
  (pixman currently supports 32bit sRGB+A, but no 24bit sRGB without alpha)
- load/save palette to/from image
- rework entirely (now: 1024 colours with mandatory interpolation)
- implement Pauldelbrot's multiwave colouring
- colour cycling (suggested by blob)
- show .kfp filename in colour dialog (suggested by gerrit)
- more flexible colouring with lighting/layers/etc (suggested by Fraktalist)
- import gradients from various formats (.ugr/.gradient, .map, .xml (flam3))
  inspired by Mandel Meute + padleywood's gradient convert tool:
  <https://fractalforums.org/f/11/t/2934>


## Getting The Code

I distribute EXEs bundled together with the corresponding source code.

The latest source code is available from my git repository:

    git clone https://code.mathr.co.uk/kalles-fraktaler-2.git
    cd kalles-fraktaler-2
    git checkout master       # for Karl's original upstream
    git checkout claude       # for MINGW build system and bug fixes
    git checkout kf-2.12      # for old stable bugfixes
    git checkout kf-2.13      # for current stable development
    git checkout kf-2.14      # for current experimental development
    git tag -l                # list available release tags

You also need `et` to generate the formula code for Newton-Raphson zooming:

    git clone https://code.mathr.co.uk/et.git

This is outside the scope of this document at the moment, easiest is to
download the release bundle and copy the `formula/generated` folder from
the included src zip.  Get in touch for help if you want to regenerate
these sources yourself.


## Building On Linux

Build instructions for cross-compiling from GNU/Linux require about 10GB of
disk space and good internet download speed (or patience).  About 750MB of
downloads including the chroot debootstrap step.  If you have recent Debian
you can skip the chroot step and install natively.

0. Setup Debian Buster chroot:

        mkdir ./vm
        sudo debootstrap buster ./vm/
        sudo mount proc ./vm/proc -t proc
        sudo mount sysfs ./vm/sys -t sysfs
        sudo cp /etc/hosts ./vm/etc/hosts
        sudo chroot ./vm /bin/bash
        cd

1. Install dependencies (inside the chroot if you made one):

        dpkg --add-architecture i386
        apt-get update
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
          pkg-config \
          wget \
          wine32 \
          wine64 \
          wine-binfmt \
          xsltproc \
          zip
        apt-get install \
          pandoc \
          texlive-fonts-recommended \
          texlive-latex-recommended   # for PDF manual
        mount binfmt_misc /proc/sys/fs/binfmt_misc -t binfmt_misc
        update-binfmts --import /usr/share/binfmts/wine

    For Ubuntu replace "wine32 wine64 wine-binfmt" with "wine" (but see note
    about build failures with some versions).

2. Prepare non-root build user:

        adduser build
        # enter and confirm password
        su - build
        mkdir -p ~/win64/src
        mkdir -p ~/win32/src

3. Download Kalles Fraktaler 2 + sources:

        cd ~/win64/src
        git clone https://code.mathr.co.uk/kalles-fraktaler-2.git
        cd kalles-fraktaler-2
        git checkout kf-2.14
        cd ..
        cp -avit ~/win32/src kalles-fraktaler-2

4. Download and build and install 3rd party library sources (inspect the script
if you want to be sure it isn't doing anything dodgy, or to copy/paste parts if
necessary), the script builds both 64bit and 32bit variants:

        cd ~/win64/src/kalles-fraktaler-2
        bash ./prepare.sh

5. Download the previous version and copy the `et`-generated formulas from it:

        cd ~/win64/src/kalles-fraktaler-2
        wget -c https://mathr.co.uk/kf/kf-2.14.6.1.zip
        unzip kf-*.zip
        cd kf-*/
        unzip kf-*src.zip
        cd kf-*-src/
        cp -avit ../../formula/generated formula/generated/*.c

6. Finally, build Kalles Fraktaler 2 + (64bit and 32bit):

        cd ~/win64/src/kalles-fraktaler-2
        ./release.sh $(git describe)

Note: build fails on Ubuntu 16.04.3 LTS (xenial):

    $ make
    x86_64-w64-mingw32-g++ -mfpmath=sse -xc++ -Wno-write-strings -pipe -MMD -g -O3 -ffast-math -I/home/claude/win64/include -DKF_THREADED_REFERENCE_BARRIER -o fraktal_sft/CDecNumber.o -c fraktal_sft/CDecNumber.cpp
    In file included from fraktal_sft/CDecNumber.cpp:1:0:
    fraktal_sft/CDecNumber.h:5:76: error: ‘decNumber’ was not declared in this scope
     typedef boost::multiprecision::number<boost::multiprecision::gmp_float<0>> decNumber;
                                                                                ^
    fraktal_sft/CDecNumber.h:5:76: error: template argument 1 is invalid
    fraktal_sft/CDecNumber.h:5:62: error: template argument 1 is invalid
     typedef boost::multiprecision::number<boost::multiprecision::gmp_float<0>> decNumber;
                                                                  ^
    fraktal_sft/CDecNumber.h:5:62: error: template argument 2 is invalid
    fraktal_sft/CDecNumber.h:5:32: warning: ‘typedef’ was ignored in this declaration
     typedef boost::multiprecision::number<boost::multiprecision::gmp_float<0>> decNumber;
                                    ^
    ...
    $ x86_64-w64-mingw32-g++ --version
    x86_64-w64-mingw32-g++ (GCC) 5.3.1 20160211
    Copyright (C) 2015 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

The working Debian Stretch has:

    $ x86_64-w64-mingw32-g++ --version
    x86_64-w64-mingw32-g++ (GCC) 6.3.0 20170516
    Copyright (C) 2016 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

The working Debian Buster has:

    $ x86_64-w64-mingw32-g++ --version
    x86_64-w64-mingw32-g++ (GCC) 8.3-win32 20190406
    Copyright (C) 2018 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


## Building on Windows

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


## Legal

Kalles Fraktaler 2
Copyright (C) 2013-2017 Karl Runmo
Copyright (C) 2017-2019 Claude Heiland-Allen

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

- this software is based in part on the work of the Independent JPEG Group
  <https://jpegclub.org/reference/libjpeg-license/>
- the TIFF library is used under the libtiff license:
  <https://gitlab.com/libtiff/libtiff/blob/master/COPYRIGHT>
- the PNG library is used under the libpng license
  <https://libpng.org/pub/png/src/libpng-LICENSE.txt>
- the ZLIB library is used under the zlib license
  <https://zlib.net/zlib_license.html>
- the GMP library is used under the conditions of the GNU Lesser General Public
  License version 3 and the GNU General Public License version 2
  <https://www.gnu.org/licenses/lgpl-3.0.en.html>
  <https://www.gnu.org/licenses/gpl-2.0.html>
- the MPFR library is used under the conditions of the GNU Lesser General Public
  License version 3
  <https://www.gnu.org/licenses/lgpl-3.0.en.html>
- the GSL library is used under the conditions of the GNU General Public License
  <https://www.gnu.org/licenses/gpl.html>
- the PIXMAN library is used under the conditions of the MIT License
  <https://cgit.freedesktop.org/pixman/tree/COPYING>
- the MINGW-STD-THREADS library is used under the conditions of the Simplified
  BSD License
  <https://github.com/meganz/mingw-std-threads/blob/master/LICENSE>
- the ILMBASE library is used under the conditions of the Modified BSD License
  <https://www.openexr.com/license.html>
- the OPENEXR library is used under the conditions of the Modified BSD License
  <https://www.openexr.com/license.html>
- the GLM library is used under the conditions of the MIT License
  <https://glm.g-truc.net/copying.txt>
- the BOOST library is used under the Boost Software License Version 1.0
  <https://www.boost.org/LICENSE_1_0.txt>
- the CLEW library is used under the Boost Software License Version 1.0
  <https://www.boost.org/LICENSE_1_0.txt>

**NOTE**:
If you redistribute the binaries or provide access to the binaries as a
service, you must also be prepared to distribute the source corresponding
to those binaries.  To make this easier for you, the more recent zips
include the source too (though you'll also need to get the third party
library sources).


## Acknowledgements

Thanks to:

- K.I.Martin for applying Perturbation and Series Approximation on the
  Mandelbrot set and generously sharing the theory and Java source code!
- Pauldelbrot for finding the reliable glitch detection method
- Botond Kósa and knighty for the extensions of Series Approximation
- laser blaster for the Burning ship formula
- stardust4ever for other fractal types
- claude for the Newton-Raphson method
- gerrit for the distance colouring differencing variations
- Dinkydau, Fractal universe, CFJH, Foxxie and others for reporting bugs
- Chillheimer for hosting my program

Claude also thanks Karl for releasing the source code and assigning a Free
Software license.


# User Manual

## File

  - **Open**

    Opens the current location from a parameter file (*.kfr)
    You can also load metadata from images saved by KF.
    You can also drag-and-drop files from the file manager on to the
    main window to open them as parameters.

  - **Save**

    Saves the current location in the current parameter file (*.kfr)

  - **Save as**

      - **KFR** Saves the current location in a new parameter file (*.kfr)

      - **PNG**
        Saves the current location in a PNG image file (*.png).
        The location and settings are saved in the file metadata.

      - **JPEG**
        Saves the current location in a JPEG image file (*.jpg).
        The location and settings are saved in the file metadata.

      - **TIFF**
        Saves the current location in a TIFF image file (*.tif).
        The location and settings are saved in the file metadata.

      - **EXR**
        Saves the current location in a EXR image file (*.exr).
        The iteration data is also saved in the file.
        The location and settings are saved in the file metadata.

      - **KFB**
        Saves the iteration data in a map file (*.kfb). This file can be used by
        the KeyFramMovie program.  **The location and settings are not saved.**

    Saved image dimensions can be smaller than the calculated image size set in
    the settings.  The image will be downscaled according to the shrink quality
    setting.  If the saved image dimensions are larger, the image will be
    recalculated at the new size.
    
    Note that for EXR, the saved image dimensions are only for the "preview image",
    the full calculated image size is used for the main image and the iteration
    data.

  - **Open map**

    Load the iteration data from map file (*.kfb).

    Note: you must set the aspect ratio of the window to match the KFB data
    before opening the map file.  If you have also saved images, you can do
    that by loading the image as a settings file.

    Note: if you want to continue zooming from the location, you must load
    the KFR file before opening the map file.  You can also load a saved image
    as a location file.  Location information is *not* stored in KFB files.

    Note: there is a historical accident whereby the iteration divider is saved
    in the KFB as an integer, losing any fractional part and sometimes resetting
    to 1 on load.  A workaround is to load the palette from a KFP file (which is
    just a .kfr renamed to .kfp) after opening the map file.

  - **Store zoom-out images**

    Zoom out automatically with the selected Zoom size and store JPEG/PNG/TIFF
    images and map file (*.kfb) for each zoom out. The zoom out stops when the
    depth is lower than 1. The resulting files can be used by the KeyFramMovie
    program to create a zoom-in animation.

  - **Examine Zoom sequence**

    Make sure you store the end location as a kfr file in the same directory
    as you store the zoom sequence frames. This function allows you to examine
    the frames one by one and add references to remove eventual visible glitch
    blobs, or choose another pixel as the main reference.

  - **Resume Zoom sequence**

    Make sure you store the end location as a kfr file in the same directory as
    you store the zoom sequence frames. This function allows you to resume and
    continue the zoom out sequnce, if it got interrupted.

  - **Open Settings**

    Opens rendering settings from a settings file (*.kfs)
    You can also load metadata from images saved by KF.

  - **Save Settings**

    Saves the current rendering settings to a settings file (*.kfs)

  - **Check for update**

    Retrieves information from KF's homepage on the internet, to let you
    know if a newer version is available.

  - **Exit**

    Exit this program

## Fraktal

  - **Presets**

    Set groups of settings to suggested preset values.

    - **Fast** accuracy may be compromised but it's fast for browsing.  Sets
      ignore isolated glitch neighbourhood to 4, enables guessing, disables
      low tolerance for glitches and approximation, and disables derivatives
      computation unless analytic DE colouring is currently in use.

    - **Best** highest quality settings for important images, but slow.  Sets
      ignore isolated glitch neighbourhood to 0 (disabled), disables guessing,
      enables low tolerance for glitches and approximation, and enables jitter
      with a default seed of 1.

  - **Location...**

    Displays the Location dialog where the coordinates for this location is
    displayed and can be edited.

  - **Iterations...**

    Displays the Iterations dialog where the maximum iteration number for this
    location is displayed and can be edited.

    The smooth color transition method is also set here, and the power on the
    Mandelbrot function.

    The fractal types is also set here - Mandelbrot, Burning Ship, Buffalo or
    Celtic.

    This dialog also displays

    - **Min:** The minimum iteration count for a pixel in this location

    - **Max:** The maximum iteration count for a pixel in this location

    - **Appr:** The number of iterations given by Series approximation

    - **Calculations:** The number of calculations performed and also the number of
      calculations per second is shown if this dialog is displayed while the
      image is rendered

  - **Set colors...**

    Displays the Number of colors dialog where the colors can be edited.

  - **Refresh**

    Render the current location

  - **Cancel rendering**

    Cancel the current rendering

  - **Reset**

    Set the location to the start point

  - **Undo**

    Go back to previous location(s).

  - **Redo**

    Go forward in the undo history.

  - **Copy**

    Save the current location to the system clipboard.

  - **Paste**

    Set the current location from the system clipboard.


## View

  - **Zoom size**

    Set the level of zoom, left mouse click to zoom in, right to zoom out

  - **Animate soom**

    Turns animation on or off when zooming

  - **Set window size**

    Set the size of the display window.

  - **Set image size**

    Set the size of the internal image size. If this is larger than the window
    size, an anti-alias effect is achieved

  - **Arbitrary size**

    (I don't know what this does?)

  - **Rotate**

    Activate rotation, drag to rotate the image

  - **Reset rotation**

    Clear any rotation

  - **Skew**

    Opens the Skew dialog which allows to "un-skew" locations that are skewed.

    Automatic unskew is available via the Newton-Raphson zooming dialog.

  - **Set Ratio**

    Enables changing the ratio between height and width of the background image
    in order to enable stretching locations. Combinated with rotation, an
    almost infinite skewing ability is enabled, useful when exploring the
    hidden treasures of the new Fractals!  (This is another name for skew.)

  - **Reset Ratio**

    Reset skew ratio to default

  - **Skew animation**

    Activates or deactivas skew animation. If activated, a popup allows you to
    specify end skew parameters and number of frames. The fractal will be
    rendered frame by frame, and can be combined with frame by frame rendering
    in KeyFrameMovieMaker or MMY3D

  - **Shrink quality**

    Set quality of image scaling.  "Fast" and "Default" are better suited
    for exploring, but "Best" looks much better at the cost of slower colouring
    (it does properly filtered anti-aliasing).

    You shouldn't need to adjust this unless you set the image size larger than
    the window size, or intend saving image files at smaller resolutions than the
    image size.

    Note: this feature is still a bit experimental and subject to change in
    future versions.

## Navigation

  - **Newton-Raphson zooming**

    When activated, a dialog will be displayed, which allows you to select if
    the zoom should jump directly to the minibrot, or to 3/4 zooms to the
    minibrot, where the current pattern is doubled, etc.

    The zoom level of the current pattern is set when opening the dialog, and
    can be changed with the capture button (which gets the current zoom level
    from the image view).

    Click on the fractal to specify the start point of the search of the
    minibrot. The current zoom size is used to set the boundaries of search
    around the selected point

    Notice that it can take a long time to calculate the position of
    deep minibrots. However, that should be still much faster than
    zooming to the minibrot manually by selecting the center of the pattern
    in the view, or with the automatic search of minibrot that is also using
    the pattern center.

    When "auto skew (newton)" is enabled before activating, the view will be
    skewed to make features near the minibrot approximately circular.

    When "auto skew (escape)" is activated, the view will be skewed to make
    features in the current view approximately circular (without zooming).

    "Auto skew escape" has an option to take into account d/dz derivatives as
    well as d/dc, some fractals work better with it enabled, some without.
    Experiment and report back!

  - **Find Minibrot**

    Starts an automatic zoom-in in the image's pattern center, until a Minibrot
    is found or if it fails to find the center.

    It's probably better to use Newton-Raphson zooming if possible.

  - **Non exact find Minibrot**

    Makes the Find Minibrot function fail every 20 zoom-in, in order to gain
    depth automatically without ending up in a Minibrot

    Newton-Raphson zooming may be a more useful option.

  - **Center cursor**

    Center the cursor to image's pattern center

  - **Find center of glitch (Color)**

    Centers the mouse pointer over the glitch blob found, if any

  - **Find highest iteration**

    Centers the mouse pointer over the pixel with the highest iteration

  - **Show Inflection**

    Activate or deactivate display of Inflection

  - **Show crosshair window**

    Display a small window that magnifies the area around the mouse cursor.
    Perfect for precisely picking particular pixels for zooming etc.

## Advanced

  - **Auto solve glitches**

    Turns the Auto solve glitches function on or off

  - **Solve glitch with near pixel method**

    Instead of re-render all pixels with the same iteration count value(color)
    only the connected pixels are re-rendered. On some locations other areas in
    the same view have the exact same iteration count values. These pixels may
    be correctly rendered and may be incorrect if re-rendered with another
    reference

  - **Ignore isolated small glitches**

    When enabled, ignores single-pixel glitches by interpolating their value
    from neighbouring pixels.  If the image size is very large, there may be
    a very large number of these tiny glitches, whose incorrect rendering may
    be invisible to the eye, and whose correct rendering may take forever.

  - **Set main reference**

    Let you click the image and select the main reference for the whole image.
    This can be useful when glitches appears on top of minibrots when the
    reference is outside this minibrot. The glitch pattern disappears from the
    minibrot if the main reference is selected inside the minibrot.

  - **Add reference (Color)**

    Add a reference and re-calculates the pixels with the same iteration count
    as the reference. This is useful if the Auto solve glitches function fails
    to find and solve glitches in the image

  - **Reuse reference**

    Do not re-calculate the reference for further zooming. This can be useful
    when during automatic zoom-out and to test different reference points, but
    must not be used together with the Auto solve glitches function active.

    Note: reuse reference cannot be used for zoom sequences in which the number
    type used for calculations changes (which happens near 1e9864 1e4932 1e616
    1e308 for power 2 formulas).  To avoid corrupt zoom out images when reuse
    reference is enabled, also enable "Use floatexp always" (deeper than 1e4932)
    or "Use long double always" (deeper than 1e308).  This may slow down
    calculations.  Alternatively, render the zoom out sequence in several
    segments, one for each number type (floatexp, scaled long double, long
    double, scaled double, double).

  - **No reuse center**

    Don't paste the previous image in the middle when zooming out.  Disabling
    this (ie, do reuse center) can be faster but can also lead to bad images.

  - **Show iterations**

    Displays the image black-and-white with the pixels with the highest
    iteration as white and the pixels with the lowest iteration as black

  - **Show smooth transition colors**

    Displays the image black-and-white representing the smoothing coefficient

  - **Show glitches**

    When activated, glitches are displayed with a solid color

  - **Mirror**

    mirrors the image around the x-axis. Can be used on the deeper half of a
    zoom sequence to a minibrot - but not too close to the minibrot and too
    close to the half...

  - **No series approximation**

    Turns the Series approximation function on or off.

  - **Use long double always**

    Use always the 80-bit long double hardware data type. This can solve some
    type of glitches

  - **Use floatexp always**

    Use always the double mantissa/integer exponent data type. This probably
    only make the render slower

  - **Use NanoMB1 (experimental)**

    For power 2 Mandelbrot only.

    Use knighty's experimental NanoMB1 algorithm for bivariate super-
    series-approximation.  Calculations are done with floatexp always.
    A regular series approximation plus perturbation pass follows for
    glitch correction.

    It is required to set the period in the Location dialog before
    enabling NanoMB1. Using Newton zoom sets the period automatically.

    It is recommended to enable Reuse Reference after recalculating with
    NanoMB1.  A new reference is calculated automatically after each Newton
    zooming, whatever the setting of Reuse Reference.

    Whether NanoMB1 is faster or not depends heavily on the location: views
    close to minis should be significantly faster than the regular 'fast'
    preset.

  - **Use NanoMB2 (experimental)**

    For power 2 Mandelbrot only.

    Use knighty's experimental NanoMB2 algorithm for bivariate super-
    series-approximation.  Calculations are done with floatexp always,
    glitch detection and correction is disabled (and may be unnecessary?).

    It is recommended to set the period limit in the Location dialog before
    enabling NanoMB2, otherwise it uses the maximum iteration count which
    may take significantly longer.  Using Newton zoom sets the period limit
    automatically.

    It is recommended to enable Reuse Reference after recalculating with
    NanoMB2.  A new reference is calculated automatically after each Newton
    zooming, whatever the setting of Reuse Reference.

    Whether NanoMB2 is faster or not depends heavily on the location: views
    close to minis should be significantly faster than the regular 'fast'
    preset.

    The RadiusScale field in the .kfs settings file controls the scaling of
    the escape radius calculated for each minibrot in the chain.  Increasing
    it may be faster but lead to visible distortion, decreasing it may help
    the distortion but slow things down.  No GUI for this control yet.

  - **Interior checking (NanoMB only)**

    For power 2 Mandelbrot rendered with NanoMB1 or NanoMB2 only.

    Use an interior checking algorithm, which may or may not speed up
    per-pixel calculations in locations with large interior regions
    visible.

    Experimental.  Correctness is to be evaluated.  Subject to change.

  - **Use auto iterations**

    Turns automatic iteration control on or off. This is on per default.

  - **Use guessing**

    Enable interpolation of neighbouring pixel data when the iteration count
    is the same.  This speeds up rendering of interior regions, but some
    colouring can lead to visible artifacts in the exterior.

  - **Threads per CPU**

    The default of 1 is sensible if KF is the only thing you're running.
    Increase the count if you want KF to peg your system, decrease it if
    you want to do other things at the same time.

    The last menu item reduces the calculated number of threads by 1,
    which can help improve system responsiveness in some configurations.

    There will always be at least one thread.

    These options only affect perturbation rendering and image colouring,
    not Newton-Raphson zooming or reference calculations.

    Note: changing this during rendering could lead to crashes, so the menu
    is disabled during rendering to prevent that.

## About

At the very top right:

  - **?**

    Open about dialog, with version information and credits.

    This also functions as a lock mechanism, preventing accidental zooming
    while a long render is taking place.

## Iterations dialog

  - Number of iterations. Increase this if the interior is "blobby".

  - Minimum iteration count achieved in the image (display only).

  - Maximum iteration count achieved in the image (display only).

  - Series approximation iteration count (display only).

  - Smooth method

      - High bailout: large escape radius gives a smoother appearance.

      - Bailout=2: small escape radius can help finding features.

  - Value of the power `p` in the fractal formula (not used in every
    formula)

  - Fractal type (formula name).  See below for details.

  - Maximum number of secondary reference points for automatic glitch
    correction.  There is a hard limit of 10000, which is also the
    default.

  - Glitch low tolerance.  When checked, glitches are more likely to be
    detected.  Disabling it can lead to bad images, but is faster.

  - Series approximation low tolerance.  When checked, series approximation
    is stricter.  Disabling it can lead to bad images, but is faster.

  - Automatic approximation terms based on number of pixels.

  - Approximation terms: number of terms used for series approximation.

  - Calculations per second (display only).

  - Real and Imag checkboxes: use these parts of `z` when considering
    bailout past the escape radius.

  - Seed R and I number boxes: start iterating `z` from these coordinates
    (default `0 + 0 i`, for best semantics it should be a critical point
    of the iteration formula, where it's `d/dz` derivative is zero).

  - Factor `a` R and I number boxes: set the complex number `a` (denoted
    `f = d + e i` in the formula list below) for TheRedshiftRider formulas.

  - Jitter seed: non-zero enables jitter with a pseudo-random-number generator
    seed value.

  - Jitter scale in units of a pixel (1 pixel is sensible in most cases).

  - Gaussian jitter is probably best left disabled (uniform jitter looks
    better).

  - Derivatives calculation can be enabled (if needed for analytic DE
    colouring) or disabled (speeds up rendering).

### Formulas

Notation:

    i            -- imaginary unit (square root of -1)
    p            -- integer power between 2 and 5 (10 for Mandelbrot)
    c = a + i b  -- pixel coordinates (parameter plane)
    z = x + i y  -- iteration variable
    w = u + i v  -- temporary variable for two-stage formulas
    f = d + i e  -- constant 'a' for TheRedshiftRider formulas
    l m ...      -- juxtaposition is multiplication
    ^            -- raise an expression to a positive integer power
    |.|          -- surrounding a real-valued expression: absolute value

Formulas:

  -   Mandelbrot

          z := z^p + c

  -   Burning Ship

          z := (|x| + i |y|)^p + c

  -   Buffalo

          w := z^p
          z := (|u| + i |v|) + c

  -   Celtic

          w := z^p
          z := (|u| + i v) + c

  -   Mandelbar

          z := (x - i y)^p + c

  -   Mandelbar Celtic

          w := (x - i y)^2
          z := (|u| + i v) + c

  -   Perpendicular Mandelbrot

          z := (|x| - i y)^2 + c

  -   Perpendicular Burning Ship

          z := (x - i |y|)^2 + c

  -   Perpendicular Celtic

          w := (|x| - i y)^2
          z := (|u| + i v) + c

  -   Perpendicular Buffalo

          w := (x - i |y|)^2
          z := (|u| + i v) + c

  -   Cubic Quasi Burning Ship

          z := (|x| (x^2 - 3 y^2) - i |y (3 x^2 - y^2)|) + c

  -   Cubic Partial BS Real

          z := (|x| (x^2 - 3 y^2) + i y (3 x^2 - y^2)) + c

  -   Cubic Partial BS Imag

          z := (x (x^2 - 3 y^2) + i |y| (3 x^2 - y^2)) + c

  -   Cubic Flying Squirrel (Buffalo Imag)

          z := (x (x^2 - 3 y^2) + i |y (3 x^2 - y^2)|) + c

  -   Cubic Quasi Perpendicular

          z := (|x| (x^2 - 3 y^2) - i y |3 x^2 - y^2|) + c

  -   4th Burning Ship Partial Imag

          z := (x + i |y|)^4 + c

  -   4th Burning Ship Partial Real

          z := (|x| + i y)^4 + c

  -   4th Burning Ship Partial Real Mbar

          z := (|x| - i y)^4 + c

  -   4th Celtic Burning Ship Partial Imag

          w := (x + i |y|)^4
          z := (|u| + i v) + c

  -   4th Celtic Burning Ship Partial Real

          w := (|x| + i y)^4
          z := (|u| + i v) + c

  -   4th Celtic Burning Ship Partial Real Mbar

          w := (|x| - i |y|)^4
          z := (|u| + i v) + c

  -   4th Buffalo Partial Imag

          w := z^4
          z := (u + i |v|) + c

  -   4th Celtic Mbar

          w := (x - i y)^4
          z := (|u| + i v) + c

  -   4th False Quasi Perpendicular

          z := ((x^4 + y^4 - 6 x^2 y^2) - i 4 x y |x^2 - y^2|) + c

  -   4th False Quasi Heart

          z := ((x^4 + y^4 - 6 x^2 y^2) + i 4 x y |x^2 - y^2|) + c

  -   4th Celtic False Quasi Perpendicular

          z := (|x^4 + y^4 - 6 x^2 y^2| - i 4 x y |x^2 - y^2|) + c

  -   4th Celtic False Quasi Heart

          z := (|x^4 + y^4 - 6 x^2 y^2| + i 4 x y |x^2 - y^2|) + c

  -   5th Burning Ship Partial

          z := (|x| + i y)^5 + c

  -   5th Burning Ship Partial Mbar

          z := (|x| - i y)^5 + c

  -   5th Celtic Mbar

          w := (x - i y)^5
          z := (|u| + i v) + c

  -   5th Quasi Burning Ship (BS/Buffalo Hybrid)

          w := (|x| + i y)^5
          z := (u - i |v|) + c

  -   5th Quasi Perpendicular

          z := (|x| (x^4 + 5 y^4 - 10 x^2 y^2) - i y (|5 x^4 + y^4 - 10 x^2 y^2|)) + c

  -   5th Quasi Heart

          z := (|x| (x^4 + 5 y^4 - 10 x^2 y^2) + i y (|5 x^4 + y^4 - 10 x^2 y^2|)) + c

  -   SimonBrot 4th

          z := z^2 (|x| + i |y|)^2 + c

  -   4th Imag Quasi Perpendicular / Heart

          z := ((x^4 + y^4 - 6 x^2 y^2) + i 4 x |y (x^2 - y^2)|) + c

  -   4th Real Quasi Perpendicular

          z := ((x^4 + y^4 - 6 x^2 y^2) - i 4 y |x (x^2 - y^2)|) + c

  -   4th Real Quasi Heart

          z := ((x^4 + y^4 - 6 x^2 y^2) + i 4 y |x (x^2 - y^2)|) + c

  -   4th Celtic Imag Quasi Perpendicular / Heart

          z := (|x^4 + y^4 - 6 x^2 y^2| + i 4 x |y (x^2 - y^2)|) + c

  -   4th Celtic Real Quasi Perpendicular

          z := (|x^4 + y^4 - 6 x^2 y^2| - i 4 y |x (x^2 - y^2)|) + c

  -   4th Celtic Real Quasi Heart

          z := (|x^4 + y^4 - 6 x^2 y^2| + i 4 y |x (x^2 - y^2)|) + c

  -   SimonBrot 6th

          z := z^3 (|x| + i |y|)^3 + c

  -   HPDZ Buffalo

          z := (((x^2 - y^2) - |x|) + i (|2xy| - |y|)) + c

  -   TheRedshiftRider 1: `a*z^2+z^3+c`

          z := (f z^2 + z^3) + c

  -   TheRedshiftRider 2: `a*z^2-z^3+c`

          z := (f z^2 - z^3) + c

  -   TheRedshiftRider 3: `2*z^2-z^3+c`

          z := (2 z^2 - z^3) + c

  -   TheRedshiftRider 4: `a*z^2+z^4+c`

          z := (f z^2 + z^4) + c

  -   TheRedshiftRider 5: `a*z^2-z^4+c`

          z := (f z^2 - z^4) + c

  -   TheRedshiftRider 6: `a*z^2+z^5+c`

          z := (f z^2 + z^5) + c

  -   TheRedshiftRider 7: `a*z^2-z^5+c`

          z := (f z^2 - z^5) + c

  -   TheRedshiftRider 8: `a*z^2+z^6+c`

          z := (f z^2 + z^6) + c

  -   TheRedshiftRider 9: `a*z^2-z^6+c`

          z := (f z^2 - z^6) + c

  -   SimonBrot2 4th

          w := z^2
          z := w (|u| + i |v|) + c

  -   General Quadratic Minus

          z := ((x^2 - y^2) + i (2 d x y + e x^2)) + c

  -   General Quadratic Plus

          z := ((x^2 + y^2) + i (2 d x y + e x^2)) + c

  -   Mothbrot 2nd 1x1

          z := z^1 (|x| + i |y|)^1 + c

  -   Mothbrot 3rd 1x2

          z := z^1 (|x| + i |y|)^2 + c

  -   Mothbrot 3rd 2x1

          z := z^2 (|x| + i |y|)^1 + c

  -   Mothbrot 4th 1x3

          z := z^1 (|x| + i |y|)^3 + c

  -   Mothbrot 4th 2x2 (aka SimonBrot 4th)

          z := z^2 (|x| + i |y|)^2 + c

  -   Mothbrot 4th 3x1

          z := z^3 (|x| + i |y|)^1 + c

  -   Mothbrot 5th 1x4

          z := z^1 (|x| + i |y|)^4 + c

  -   Mothbrot 5th 2x3

          z := z^2 (|x| + i |y|)^3 + c

  -   Mothbrot 5th 3x2

          z := z^3 (|x| + i |y|)^2 + c

  -   Mothbrot 5th 4x1

          z := z^4 (|x| + i |y|)^1 + c

  -   Mothbrot 6th 1x5

          z := z^1 (|x| + i |y|)^5 + c

  -   Mothbrot 6th 2x4 (Simon's Mothbrot)

          z := z^2 (|x| + i |y|)^4 + c

  -   Mothbrot 6th 3x3 (aka SimonBrot 6th)

          z := z^3 (|x| + i |y|)^3 + c

  -   Mothbrot 6th 4x2

          z := z^4 (|x| + i |y|)^2 + c

  -   Mothbrot 6th 5x1

          z := z^5 (|x| + i |y|)^1 + c

  -   General Abs Quadratic Minus

          z := ((x^2 - y^2) + i (2 d |x y| + e x^2)) + c

  -   General Abs Quadratic Plus

          z := ((x^2 + y^2) + i (2 d |x y| + e x^2)) + c

A machine-readable version of this formula list is found in the 'et' repository:
<https://code.mathr.co.uk/et/blob/kf:/kf/formulas.et>  This is used by 'et' when
generating formula code (for Newton-Raphson zooming, etc).


## Number of colors dialog

  - **Number of key colors**

    Set the number of key colors between 1 and 1024.

  - **Divide iteration**

    Divide each iteration number with this value, for dense images this value
    can be greater than 1.  For DE, values less than 1 can be useful.

  - **Color offset**

    Offset the colors in the palette

  - **Random**

    Fill the palette with random colors made from the Seed value. The Seed
    button select a seed value randomly.

  - **More contrast**

    Move RGB values closer to max or min

  - **Less contrast**

    Move RGB values closer to the middle

  - **Show slopes**

    Enable slope encoding for 3D effect.

    First value is the magnification of the slopes. The start value of 100 is
    suitable for the unzoomed view. Deep views requires a couple of magnitudes
    higher value.

    The second value is the percentage with which the slope encoding is applied
    on the coloring. 100 is max, however flat areas will still have the palette
    color visible.

  - **Save palette**

    Save the current palette in KFP (*.kfp) file

  - **Open palette**

    Load palette from a KFP (*.kfp) file

  - **Expand double**

    Double the number of key colors without changing the palette. This allows
    finer control of individual colors without changing the palette for other
    colors

  - **Expand all**

    Increase the number of key color to maximum 1024 without changing the
    palette

  - **Double**

    Double the key colors by repeating them

  - **Merge Colors**

    Allows a selected color to be merged to every specied key color

  - **Show index**

    Capture the mouse, hover the mouse over the fractal image and the
    corresponding color in the list will be highlighted. Click and the color
    selection dialog will be displayed for the active color

  - **Smooth color transition**

    Makes the transitions of colors smooth

  - **Inverse smooth color transition**

    Inverse the smooth color transition which makes edges more visible

  - **Unnamed dropdown box**

    Specifies handling of the iteration count values prior to coloring

  - **Palette waves**

    The palette can be filled from sine waves applied on Red, Green, Blue and
    Black-and-white. Each input box specifies the number of periods applied on
    the number of key colors in the palette. If the input box is left empty, no
    wave of this color is applied. At right of each input box the "P"-button
    makes the number you entered prime, since different prime numbers probably
    give more variation. The last input box specifies the waves offset.

    The button "Generate" applies the waves on the palette, the "Seed" button
    fills the fields with random values

  - **Infinite waves**

    Waves can be applied on Hue, Saturation and Brightness rather than RGB
    values. The Period value specifies the length of the period (not the
    number of periods as for the Palette waves). Periods with prime numbers
    should be able to produce an infinite number unique colors

    A negative value on Hue, Saturation or Brightness makes a flat percentage
    value to be applied on all iterations.


## Command Line Usage

    kf.exe [options]
        -o, --load-map      [FILE.kfb]  load map file
        -c, --load-palette  [FILE.kfp]  load palette file
        -l, --load-location [FILE.kfr]  load location file
        -s, --load-settings [FILE.kfs]  load settings file
        -x, --save-exr      [FILE.exr]  save EXR
        -t, --save-tif      [FILE.tif]  save TIFF
        -p, --save-png      [FILE.png]  save PNG
        -j, --save-jpg      [FILE.jpg]  save JPEG
        -m, --save-map      [FILE.kfb]  save KFB
        -z, --zoom-out      [NFRAMES]   zoom sequence
        --log (debug|status|info|warn|error)
                                        logging verbosity
        -v, -V, --version               show version
        -h, -H, -?, --help              show this help

Locations and settings can also be image files with embedded comments.

If any of the save options are give, KF switches to a non-interactive mode - it
will render the image and save to all specified types before quitting.  No GUI.

A typical workflow would be to start KF without arguments, set up the window
size (eg 640x360), image size (eg 3840x2160), glitch low tolerance flag, etc,
then save the settings to a .kfs file, before quitting.

Then launch KF from the command line telling it to load the settings file you
just saved, plus the location file you want to render, and where to save the
output images to.  Then wait patiently.  You can write a script that renders
multiple locations in succession, either in batch .BAT on Windows, or in Shell
.sh on *nix with WINE.

Note that you might have to double up backslashes within quoted paths (if they
contain spaces).  Maybe forward slashes work instead, but you do need quotes
(either single `''` or double `""`, in a matching pair around the whole argument)
if there are spaces.  Your shell might also do funky stuff with slashes and
quotes, so you might need to double up or quadruple up the backslashes and
quotes.  Easiest to avoid spaces and keep your files in the current working
directory...  Example:

    kf.exe -s settings.kfs -l location.kfr -p out.png -j out.jpg -m out.kfb

Use `--log info` to disable the status updates, use `--log warn` to output only
important messages.  The default is `--log status`.

With `-z` you can specify how many frames to render, or `-1` to zoom all the
way out.  In zoom out mode the PNG, JPEG, and KFB filenames should contain a
printf flag for an integer, for example `image-%08d.png` will have 8 decimal
digits padded with leading 0. This is filled by the frame number, which always
starts from 0.  Zooming is by the zoom size in the settings file.

New in 2.14.6 is standalone KFB map colouring support with the `-o`/`--load-map`
flag:

    kf.exe -o map.kfb -c palette.kfp -p out.png
