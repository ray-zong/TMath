#ifndef SOFAMHDEF
#define SOFAMHDEF

/*
**  - - - - - - - -
**   s o f a m . h
**  - - - - - - - -
**
**  Macros used by SOFA library.
**
**  This file is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Please note that the constants defined below are to be used only in
**  the context of the SOFA software, and have no other official IAU
**  status.  In addition, self consistency is not guaranteed.
**
**  This revision:   2013 August 27
**
**  SOFA release 2016-05-03
**
**  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
*/

/* Star-independent astrometry parameters */
typedef struct {
   double pmt;        /* PM time interval (SSB, Julian years) */
   double eb[3];      /* SSB to observer (vector, au) */
   double eh[3];      /* Sun to observer (unit vector) */
   double em;         /* distance from Sun to observer (au) */
   double v[3];       /* barycentric observer velocity (vector, c) */
   double bm1;        /* sqrt(1-|v|^2): reciprocal of Lorenz factor */
   double bpn[3][3];  /* bias-precession-nutation matrix */
   double along;      /* longitude + s' + dERA(DUT) (radians) */
   double phi;        /* geodetic latitude (radians) */
   double xpl;        /* polar motion xp wrt local meridian (radians) */
   double ypl;        /* polar motion yp wrt local meridian (radians) */
   double sphi;       /* sine of geodetic latitude */
   double cphi;       /* cosine of geodetic latitude */
   double diurab;     /* magnitude of diurnal aberration vector */
   double eral;       /* "local" Earth rotation angle (radians) */
   double refa;       /* refraction constant A (radians) */
   double refb;       /* refraction constant B (radians) */
} iauASTROM;
/* (Vectors eb, eh, em and v are all with respect to BCRS axes.) */

/* Body parameters for light deflection */
typedef struct {
   double bm;         /* mass of the body (solar masses) */
   double dl;         /* deflection limiter (radians^2/2) */
   double pv[2][3];   /* barycentric PV of the body (au, au/day) */
} iauLDBODY;

/* Pi */
#define DPI (3.141592653589793238462643)

/* 2Pi */
#define D2PI (6.283185307179586476925287)

/* Radians to degrees */
#define DR2D (57.29577951308232087679815)

/* Degrees to radians */
#define DD2R (1.745329251994329576923691e-2)

/* Radians to arcseconds */
#define DR2AS (206264.8062470963551564734)

/* Arcseconds to radians */
#define DAS2R (4.848136811095359935899141e-6)

/* Seconds of time to radians */
#define DS2R (7.272205216643039903848712e-5)

/* Arcseconds in a full circle */
#define TURNAS (1296000.0)

/* Milliarcseconds to radians */
#define DMAS2R (DAS2R / 1e3)

/* Length of tropical year B1900 (days) */
#define DTY (365.242198781)

/* Seconds per day. */
#define DAYSEC (86400.0)

/* Days per Julian year */
#define DJY (365.25)

/* Days per Julian century */
#define DJC (36525.0)

/* Days per Julian millennium */
#define DJM (365250.0)

/* Reference epoch (J2000.0), Julian Date */
#define DJ00 (2451545.0)

/* Julian Date of Modified Julian Date zero */
#define DJM0 (2400000.5)

/* Reference epoch (J2000.0), Modified Julian Date */
#define DJM00 (51544.5)

/* 1977 Jan 1.0 as MJD */
#define DJM77 (43144.0)

/* TT minus TAI (s) */
#define TTMTAI (32.184)

/* Astronomical unit (m) */
#define DAU (149597870e3)

/* Speed of light (m/s) */
#define CMPS 299792458.0

/* Light time for 1 au (s) */
#define AULT 499.004782

/* Speed of light (AU per day) */
#define DC (DAYSEC / AULT)

/* L_G = 1 - d(TT)/d(TCG) */
#define ELG (6.969290134e-10)

/* L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
#define ELB (1.550519768e-8)
#define TDB0 (-6.55e-5)

/* Schwarzschild radius of the Sun (au) */
/* = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11 */
#define SRS 1.97412574336e-8

/* dint(A) - truncate to nearest whole number towards zero (double) */
#define dint(A) ((A)<0.0?ceil(A):floor(A))

/* dnint(A) - round to nearest whole number (double) */
#define dnint(A) ((A)<0.0?ceil((A)-0.5):floor((A)+0.5))

/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-fabs(A):fabs(A))

/* max(A,B) - larger (most +ve) of two numbers (generic) */
#define gmax(A,B) (((A)>(B))?(A):(B))

/* min(A,B) - smaller (least +ve) of two numbers (generic) */
#define gmin(A,B) (((A)<(B))?(A):(B))

/* Reference ellipsoids */
#define WGS84 1
#define GRS80 2
#define WGS72 3

#endif
