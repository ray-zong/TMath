#ifndef TSOFA_H
#define TSOFA_H

/*
**  Prototype function declarations for SOFA library.
**
**  This file is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  This revision:   2016 March 10
**
**  SOFA release 2016-05-03
*/
#include <cmath>

#include "TSofam.h"
#include "TGlobal.h"

TMATH_BEGIN_NAMESPACE

class TMATH_EXPORT TSofa
{
public:
    /* Astronomy/Calendars 历元*/
    static inline int iauCal2jd(int iy, int im, int id, double *djm0, double *djm);
    static inline double iauEpb(double dj1, double dj2);
    static inline void iauEpb2jd(double epb, double *djm0, double *djm);
    static inline double iauEpj(double dj1, double dj2);
    static inline void iauEpj2jd(double epj, double *djm0, double *djm);
    static inline int iauJd2cal(double dj1, double dj2,
        int *iy, int *im, int *id, double *fd);
    static inline int iauJdcalf(int ndp, double dj1, double dj2, int iymdf[4]);

    /* Astronomy/Astrometry 天体测量*/
    static inline void iauAb(double pnat[3], double v[3], double s, double bm1,
        double ppr[3]);
    static inline void iauApcg(double date1, double date2,
        double ebpv[2][3], double ehp[3],
        iauASTROM *astrom);
    static inline void iauApcg13(double date1, double date2, iauASTROM *astrom);
    static inline void iauApci(double date1, double date2,
        double ebpv[2][3], double ehp[3],
        double x, double y, double s,
        iauASTROM *astrom);
    static inline void iauApci13(double date1, double date2,
        iauASTROM *astrom, double *eo);
    static inline void iauApco(double date1, double date2,
        double ebpv[2][3], double ehp[3],
        double x, double y, double s, double theta,
        double elong, double phi, double hm,
        double xp, double yp, double sp,
        double refa, double refb,
        iauASTROM *astrom);
    static inline int iauApco13(double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        iauASTROM *astrom, double *eo);
    static inline void iauApcs(double date1, double date2, double pv[2][3],
        double ebpv[2][3], double ehp[3],
        iauASTROM *astrom);
    static inline void iauApcs13(double date1, double date2, double pv[2][3],
        iauASTROM *astrom);
    static inline void iauAper(double theta, iauASTROM *astrom);
    static inline void iauAper13(double ut11, double ut12, iauASTROM *astrom);
    static inline void iauApio(double sp, double theta,
        double elong, double phi, double hm, double xp, double yp,
        double refa, double refb,
        iauASTROM *astrom);
    static inline int iauApio13(double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        iauASTROM *astrom);
    static inline void iauAtci13(double rc, double dc,
        double pr, double pd, double px, double rv,
        double date1, double date2,
        double *ri, double *di, double *eo);
    static inline void iauAtciq(double rc, double dc, double pr, double pd,
        double px, double rv, iauASTROM *astrom,
        double *ri, double *di);
    static inline void iauAtciqn(double rc, double dc, double pr, double pd,
        double px, double rv, iauASTROM *astrom,
        int n, iauLDBODY b[], double *ri, double *di);
    static inline void iauAtciqz(double rc, double dc, iauASTROM *astrom,
        double *ri, double *di);
    static inline int iauAtco13(double rc, double dc,
        double pr, double pd, double px, double rv,
        double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        double *aob, double *zob, double *hob,
        double *dob, double *rob, double *eo);
    static inline void iauAtic13(double ri, double di,
        double date1, double date2,
        double *rc, double *dc, double *eo);
    static inline void iauAticq(double ri, double di, iauASTROM *astrom,
        double *rc, double *dc);
    static inline void iauAticqn(double ri, double di, iauASTROM *astrom,
        int n, iauLDBODY b[], double *rc, double *dc);
    static inline int iauAtio13(double ri, double di,
        double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        double *aob, double *zob, double *hob,
        double *dob, double *rob);
    static inline void iauAtioq(double ri, double di, iauASTROM *astrom,
        double *aob, double *zob,
        double *hob, double *dob, double *rob);
    static inline int iauAtoc13(const char *type, double ob1, double ob2,
        double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        double *rc, double *dc);
    static inline int iauAtoi13(const char *type, double ob1, double ob2,
        double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        double *ri, double *di);
    static inline void iauAtoiq(const char *type,
        double ob1, double ob2, iauASTROM *astrom,
        double *ri, double *di);
    static inline void iauLd(double bm, double p[3], double q[3], double e[3],
        double em, double dlim, double p1[3]);
    static inline void iauLdn(int n, iauLDBODY b[], double ob[3], double sc[3],
        double sn[3]);
    static inline void iauLdsun(double p[3], double e[3], double em, double p1[3]);
    static inline void iauPmpx(double rc, double dc, double pr, double pd,
        double px, double rv, double pmt, double pob[3],
        double pco[3]);
    static inline int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
        double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b,
        double *ra2, double *dec2, double *pmr2, double *pmd2,
        double *px2, double *rv2);
    static inline void iauPvtob(double elong, double phi, double height, double xp,
        double yp, double sp, double theta, double pv[2][3]);
    static inline void iauRefco(double phpa, double tc, double rh, double wl,
        double *refa, double *refb);

    /* Astronomy/Ephemerides 星历表*/
    static inline int iauEpv00(double date1, double date2,
        double pvh[2][3], double pvb[2][3]);
    static inline int iauPlan94(double date1, double date2, int np, double pv[2][3]);

    /* Astronomy/FundamentalArgs 基本*/
    static inline double iauFad03(double t);
    static inline double iauFae03(double t);
    static inline double iauFaf03(double t);
    static inline double iauFaju03(double t);
    static inline double iauFal03(double t);
    static inline double iauFalp03(double t);
    static inline double iauFama03(double t);
    static inline double iauFame03(double t);
    static inline double iauFane03(double t);
    static inline double iauFaom03(double t);
    static inline double iauFapa03(double t);
    static inline double iauFasa03(double t);
    static inline double iauFaur03(double t);
    static inline double iauFave03(double t);

    /* Astronomy/PrecNutPolar */
    static inline void iauBi00(double *dpsibi, double *depsbi, double *dra);
    static inline void iauBp00(double date1, double date2,
        double rb[3][3], double rp[3][3], double rbp[3][3]);
    static inline void iauBp06(double date1, double date2,
        double rb[3][3], double rp[3][3], double rbp[3][3]);
    static inline void iauBpn2xy(double rbpn[3][3], double *x, double *y);
    static inline void iauC2i00a(double date1, double date2, double rc2i[3][3]);
    static inline void iauC2i00b(double date1, double date2, double rc2i[3][3]);
    static inline void iauC2i06a(double date1, double date2, double rc2i[3][3]);
    static inline void iauC2ibpn(double date1, double date2, double rbpn[3][3],
        double rc2i[3][3]);
    static inline void iauC2ixy(double date1, double date2, double x, double y,
        double rc2i[3][3]);
    static inline void iauC2ixys(double x, double y, double s, double rc2i[3][3]);
    static inline void iauC2t00a(double tta, double ttb, double uta, double utb,
        double xp, double yp, double rc2t[3][3]);
    static inline void iauC2t00b(double tta, double ttb, double uta, double utb,
        double xp, double yp, double rc2t[3][3]);
    static inline void iauC2t06a(double tta, double ttb, double uta, double utb,
        double xp, double yp, double rc2t[3][3]);
    static inline void iauC2tcio(double rc2i[3][3], double era, double rpom[3][3],
        double rc2t[3][3]);
    static inline void iauC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
        double rc2t[3][3]);
    static inline void iauC2tpe(double tta, double ttb, double uta, double utb,
        double dpsi, double deps, double xp, double yp,
        double rc2t[3][3]);
    static inline void iauC2txy(double tta, double ttb, double uta, double utb,
        double x, double y, double xp, double yp,
        double rc2t[3][3]);
    static inline double iauEo06a(double date1, double date2);
    static inline double iauEors(double rnpb[3][3], double s);
    static inline void iauFw2m(double gamb, double phib, double psi, double eps,
        double r[3][3]);
    static inline void iauFw2xy(double gamb, double phib, double psi, double eps,
        double *x, double *y);
    static inline void iauLtp(double epj, double rp[3][3]);
    static inline void iauLtpb(double epj, double rpb[3][3]);
    static inline void iauLtpecl(double epj, double vec[3]);
    static inline void iauLtpequ(double epj, double veq[3]);
    static inline void iauNum00a(double date1, double date2, double rmatn[3][3]);
    static inline void iauNum00b(double date1, double date2, double rmatn[3][3]);
    static inline void iauNum06a(double date1, double date2, double rmatn[3][3]);
    static inline void iauNumat(double epsa, double dpsi, double deps, double rmatn[3][3]);
    static inline void iauNut00a(double date1, double date2, double *dpsi, double *deps);
    static inline void iauNut00b(double date1, double date2, double *dpsi, double *deps);
    static inline void iauNut06a(double date1, double date2, double *dpsi, double *deps);
    static inline void iauNut80(double date1, double date2, double *dpsi, double *deps);
    static inline void iauNutm80(double date1, double date2, double rmatn[3][3]);
    static inline double iauObl06(double date1, double date2);
    static inline double iauObl80(double date1, double date2);
    static inline void iauP06e(double date1, double date2,
        double *eps0, double *psia, double *oma, double *bpa,
        double *bqa, double *pia, double *bpia,
        double *epsa, double *chia, double *za, double *zetaa,
        double *thetaa, double *pa,
        double *gam, double *phi, double *psi);
    static inline void iauPb06(double date1, double date2,
        double *bzeta, double *bz, double *btheta);
    static inline void iauPfw06(double date1, double date2,
        double *gamb, double *phib, double *psib, double *epsa);
    static inline void iauPmat00(double date1, double date2, double rbp[3][3]);
    static inline void iauPmat06(double date1, double date2, double rbp[3][3]);
    static inline void iauPmat76(double date1, double date2, double rmatp[3][3]);
    static inline void iauPn00(double date1, double date2, double dpsi, double deps,
        double *epsa,
        double rb[3][3], double rp[3][3], double rbp[3][3],
        double rn[3][3], double rbpn[3][3]);
    static inline void iauPn00a(double date1, double date2,
        double *dpsi, double *deps, double *epsa,
        double rb[3][3], double rp[3][3], double rbp[3][3],
        double rn[3][3], double rbpn[3][3]);
    static inline void iauPn00b(double date1, double date2,
        double *dpsi, double *deps, double *epsa,
        double rb[3][3], double rp[3][3], double rbp[3][3],
        double rn[3][3], double rbpn[3][3]);
    static inline void iauPn06(double date1, double date2, double dpsi, double deps,
        double *epsa,
        double rb[3][3], double rp[3][3], double rbp[3][3],
        double rn[3][3], double rbpn[3][3]);
    static inline void iauPn06a(double date1, double date2,
        double *dpsi, double *deps, double *epsa,
        double rb[3][3], double rp[3][3], double rbp[3][3],
        double rn[3][3], double rbpn[3][3]);
    static inline void iauPnm00a(double date1, double date2, double rbpn[3][3]);
    static inline void iauPnm00b(double date1, double date2, double rbpn[3][3]);
    static inline void iauPnm06a(double date1, double date2, double rnpb[3][3]);
    static inline void iauPnm80(double date1, double date2, double rmatpn[3][3]);
    static inline void iauPom00(double xp, double yp, double sp, double rpom[3][3]);
    static inline void iauPr00(double date1, double date2,
        double *dpsipr, double *depspr);
    static inline void iauPrec76(double date01, double date02,
        double date11, double date12,
        double *zeta, double *z, double *theta);
    static inline double iauS00(double date1, double date2, double x, double y);
    static inline double iauS00a(double date1, double date2);
    static inline double iauS00b(double date1, double date2);
    static inline double iauS06(double date1, double date2, double x, double y);
    static inline double iauS06a(double date1, double date2);
    static inline double iauSp00(double date1, double date2);
    static inline void iauXy06(double date1, double date2, double *x, double *y);
    static inline void iauXys00a(double date1, double date2,
        double *x, double *y, double *s);
    static inline void iauXys00b(double date1, double date2,
        double *x, double *y, double *s);
    static inline void iauXys06a(double date1, double date2,
        double *x, double *y, double *s);

    /* Astronomy/RotationAndTime */
    static inline double iauEe00(double date1, double date2, double epsa, double dpsi);
    static inline double iauEe00a(double date1, double date2);
    static inline double iauEe00b(double date1, double date2);
    static inline double iauEe06a(double date1, double date2);
    static inline double iauEect00(double date1, double date2);
    static inline double iauEqeq94(double date1, double date2);
    static inline double iauEra00(double dj1, double dj2);
    static inline double iauGmst00(double uta, double utb, double tta, double ttb);
    static inline double iauGmst06(double uta, double utb, double tta, double ttb);
    static inline double iauGmst82(double dj1, double dj2);
    static inline double iauGst00a(double uta, double utb, double tta, double ttb);
    static inline double iauGst00b(double uta, double utb);
    static inline double iauGst06(double uta, double utb, double tta, double ttb,
        double rnpb[3][3]);
    static inline double iauGst06a(double uta, double utb, double tta, double ttb);
    static inline double iauGst94(double uta, double utb);

    /* Astronomy/SpaceMotion */
    static inline int iauPvstar(double pv[2][3], double *ra, double *dec,
        double *pmr, double *pmd, double *px, double *rv);
    static inline int iauStarpv(double ra, double dec,
        double pmr, double pmd, double px, double rv,
        double pv[2][3]);

    /* Astronomy/StarCatalogs */
    static inline void iauFk52h(double r5, double d5,
        double dr5, double dd5, double px5, double rv5,
        double *rh, double *dh,
        double *drh, double *ddh, double *pxh, double *rvh);
    static inline void iauFk5hip(double r5h[3][3], double s5h[3]);
    static inline void iauFk5hz(double r5, double d5, double date1, double date2,
        double *rh, double *dh);
    static inline void iauH2fk5(double rh, double dh,
        double drh, double ddh, double pxh, double rvh,
        double *r5, double *d5,
        double *dr5, double *dd5, double *px5, double *rv5);
    static inline void iauHfk5z(double rh, double dh, double date1, double date2,
        double *r5, double *d5, double *dr5, double *dd5);
    static inline int iauStarpm(double ra1, double dec1,
        double pmr1, double pmd1, double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b,
        double *ra2, double *dec2,
        double *pmr2, double *pmd2, double *px2, double *rv2);

    /* Astronomy/EclipticCoordinates 黄道坐标*/
    static inline void iauEceq06(double date1, double date2, double dl, double db,
        double *dr, double *dd);
    static inline void iauEcm06(double date1, double date2, double rm[3][3]);
    static inline void iauEqec06(double date1, double date2, double dr, double dd,
        double *dl, double *db);
    static inline void iauLteceq(double epj, double dl, double db, double *dr, double *dd);
    static inline void iauLtecm(double epj, double rm[3][3]);
    static inline void iauLteqec(double epj, double dr, double dd, double *dl, double *db);

    /* Astronomy/GalacticCoordinates 银河坐标*/
    static inline void iauG2icrs(double dl, double db, double *dr, double *dd);
    static inline void iauIcrs2g(double dr, double dd, double *dl, double *db);

    /* Astronomy/GeodeticGeocentric 大地地心*/
    static inline int iauEform(int n, double *a, double *f);
    static inline int iauGc2gd(int n, double xyz[3],
        double *elong, double *phi, double *height);
    static inline int iauGc2gde(double a, double f, double xyz[3],
        double *elong, double *phi, double *height);
    static inline int iauGd2gc(int n, double elong, double phi, double height,
        double xyz[3]);
    static inline int iauGd2gce(double a, double f,
        double elong, double phi, double height, double xyz[3]);

    /* Astronomy/Timescales 时间尺度*/
    static inline int iauD2dtf(const char *scale, int ndp, double d1, double d2,
        int *iy, int *im, int *id, int ihmsf[4]);
    static inline int iauDat(int iy, int im, int id, double fd, double *deltat);
    static inline double iauDtdb(double date1, double date2,
        double ut, double elong, double u, double v);
    static inline int iauDtf2d(const char *scale, int iy, int im, int id,
        int ihr, int imn, double sec, double *d1, double *d2);
    static inline int iauTaitt(double tai1, double tai2, double *tt1, double *tt2);
    static inline int iauTaiut1(double tai1, double tai2, double dta,
        double *ut11, double *ut12);
    static inline int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2);
    static inline int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2);
    static inline int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2);
    static inline int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2);
    static inline int iauTdbtt(double tdb1, double tdb2, double dtr,
        double *tt1, double *tt2);
    static inline int iauTttai(double tt1, double tt2, double *tai1, double *tai2);
    static inline int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2);
    static inline int iauTttdb(double tt1, double tt2, double dtr,
        double *tdb1, double *tdb2);
    static inline int iauTtut1(double tt1, double tt2, double dt,
        double *ut11, double *ut12);
    static inline int iauUt1tai(double ut11, double ut12, double dta,
        double *tai1, double *tai2);
    static inline int iauUt1tt(double ut11, double ut12, double dt,
        double *tt1, double *tt2);
    static inline int iauUt1utc(double ut11, double ut12, double dut1,
        double *utc1, double *utc2);
    static inline int iauUtctai(double utc1, double utc2, double *tai1, double *tai2);
    static inline int iauUtcut1(double utc1, double utc2, double dut1,
        double *ut11, double *ut12);

    /* VectorMatrix/AngleOps 角度运算*/
    static inline void iauA2af(int ndp, double angle, char *sign, int idmsf[4]);
    static inline void iauA2tf(int ndp, double angle, char *sign, int ihmsf[4]);
    static inline int iauAf2a(char s, int ideg, int iamin, double asec, double *rad);
    static inline double iauAnp(double a);
    static inline double iauAnpm(double a);
    static inline void iauD2tf(int ndp, double days, char *sign, int ihmsf[4]);
    static inline int iauTf2a(char s, int ihour, int imin, double sec, double *rad);
    static inline int iauTf2d(char s, int ihour, int imin, double sec, double *days);

    /* VectorMatrix/BuildRotations 旋转*/
    static inline void iauRx(double phi, double r[3][3]);
    static inline void iauRy(double theta, double r[3][3]);
    static inline void iauRz(double psi, double r[3][3]);

    /* VectorMatrix/CopyExtendExtract */
    static inline void iauCp(double p[3], double c[3]);
    static inline void iauCpv(double pv[2][3], double c[2][3]);
    static inline void iauCr(double r[3][3], double c[3][3]);
    static inline void iauP2pv(double p[3], double pv[2][3]);
    static inline void iauPv2p(double pv[2][3], double p[3]);

    /* VectorMatrix/Initialization 初始化*/
    static inline void iauIr(double r[3][3]);
    static inline void iauZp(double p[3]);
    static inline void iauZpv(double pv[2][3]);
    static inline void iauZr(double r[3][3]);

    /* VectorMatrix/MatrixOps 矩阵运算*/
    static inline void iauRxr(double a[3][3], double b[3][3], double atb[3][3]);
    static inline void iauTr(double r[3][3], double rt[3][3]);

    /* VectorMatrix/MatrixVectorProducts 矩阵向量乘积*/
    static inline void iauRxp(double r[3][3], double p[3], double rp[3]);
    static inline void iauRxpv(double r[3][3], double pv[2][3], double rpv[2][3]);
    static inline void iauTrxp(double r[3][3], double p[3], double trp[3]);
    static inline void iauTrxpv(double r[3][3], double pv[2][3], double trpv[2][3]);

    /* VectorMatrix/RotationVectors 向量旋转*/
    static inline void iauRm2v(double r[3][3], double w[3]);
    static inline void iauRv2m(double w[3], double r[3][3]);

    /* VectorMatrix/SeparationAndAngle */
    static inline double iauPap(double a[3], double b[3]);
    static inline double iauPas(double al, double ap, double bl, double bp);
    static inline double iauSepp(double a[3], double b[3]);
    static inline double iauSeps(double al, double ap, double bl, double bp);

    /* VectorMatrix/SphericalCartesian 空间笛卡尔坐标系*/
    static inline void iauC2s(double p[3], double *theta, double *phi);
    static inline void iauP2s(double p[3], double *theta, double *phi, double *r);
    static inline void iauPv2s(double pv[2][3],
        double *theta, double *phi, double *r,
        double *td, double *pd, double *rd);
    static inline void iauS2c(double theta, double phi, double c[3]);
    static inline void iauS2p(double theta, double phi, double r, double p[3]);
    static inline void iauS2pv(double theta, double phi, double r,
        double td, double pd, double rd,
        double pv[2][3]);

    /* VectorMatrix/VectorOps 向量运算*/
    static inline double iauPdp(double a[3], double b[3]);
    static inline double iauPm(double p[3]);
    static inline void iauPmp(double a[3], double b[3], double amb[3]);
    static inline void iauPn(double p[3], double *r, double u[3]);
    static inline void iauPpp(double a[3], double b[3], double apb[3]);
    static inline void iauPpsp(double a[3], double s, double b[3], double apsb[3]);
    static inline void iauPvdpv(double a[2][3], double b[2][3], double adb[2]);
    static inline void iauPvm(double pv[2][3], double *r, double *s);
    static inline void iauPvmpv(double a[2][3], double b[2][3], double amb[2][3]);
    static inline void iauPvppv(double a[2][3], double b[2][3], double apb[2][3]);
    static inline void iauPvu(double dt, double pv[2][3], double upv[2][3]);
    static inline void iauPvup(double dt, double pv[2][3], double p[3]);
    static inline void iauPvxpv(double a[2][3], double b[2][3], double axb[2][3]);
    static inline void iauPxp(double a[3], double b[3], double axb[3]);
    static inline void iauS2xpv(double s1, double s2, double pv[2][3], double spv[2][3]);
    static inline void iauSxp(double s, double p[3], double sp[3]);
    static inline void iauSxpv(double s, double pv[2][3], double spv[2][3]);
};

TMATH_END_NAMESPACE

#endif
