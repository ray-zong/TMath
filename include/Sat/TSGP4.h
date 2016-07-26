#ifndef TSGP4_H
#define TSGP4_H

#include "TGlobal.h"
#include "TOrbitalElements.h"
#include "TECI.h"
#include "TTle.h"
#include "TDateTime.h"


TMATH_BEGIN_NAMESPACE
/**
 * @brief The simplified perturbations model 4 propagater.
 */

class TMATH_EXPORT TSGP4
{
public:
    TSGP4(const TTle &tle);
    virtual ~TSGP4();

    void setTle(const TTle &tle);
    //∑÷÷”//
    TECI findPosition(double tsince) const;
    TECI findPosition(const TDateTime &dateTime) const;

private:
    struct CommonConstants
    {
        double cosio;
        double sinio;
        double eta;
        double t2cof;
        double a3ovk2;
        double x1mth2;
        double x3thm1;
        double x7thm1;
        double aycof;
        double xlcof;
        double xnodcf;
        double c1;
        double c4;
        double omgdot; // secular rate of omega (rad/s)
        double xnodot; // secular rate of xnode (rad/s)
        double xmdot;  // secular rate of xmo   (rad/s)
    };

    struct NearSpaceConstants
    {
        double c5;
        double omgcof;
        double xmcof;
        double delmo;
        double sinmo;
        double d2;
        double d3;
        double d4;
        double t3cof;
        double t4cof;
        double t5cof;
    };

    struct DeepSpaceConstants
    {
        double gsto;
        double zmol;
        double zmos;

        /*
         * lunar / solar constants for epoch
         * applied during DeepSpaceSecular()
         */
        double sse;
        double ssi;
        double ssl;
        double ssg;
        double ssh;
        /*
         * lunar / solar constants
         * used during DeepSpaceCalculateLunarSolarTerms()
         */
        double se2;
        double si2;
        double sl2;
        double sgh2;
        double sh2;
        double se3;
        double si3;
        double sl3;
        double sgh3;
        double sh3;
        double sl4;
        double sgh4;
        double ee2;
        double e3;
        double xi2;
        double xi3;
        double xl2;
        double xl3;
        double xl4;
        double xgh2;
        double xgh3;
        double xgh4;
        double xh2;
        double xh3;
        /*
         * used during DeepSpaceCalcDotTerms()
         */
        double d2201;
        double d2211;
        double d3210;
        double d3222;
        double d4410;
        double d4422;
        double d5220;
        double d5232;
        double d5421;
        double d5433;
        double del1;
        double del2;
        double del3;
        /*
         * whether the deep space orbit is
         * geopotential resonance for 12 hour orbits
         */
        bool resonance_flag;
        /*
         * whether the deep space orbit is
         * 24h synchronous resonance
         */
        bool synchronous_flag;
    };

    struct IntegratorValues
    {
        double xndot;
        double xnddt;
        double xldot;
    };

    struct IntegratorConstants
    {
        /*
         * integrator constants
         */
        double xfact;
        double xlamo;

        /*
         * integrator values for epoch
         */
        struct IntegratorValues values_0;
    };

    struct IntegratorParams
    {
        /*
         * integrator values
         */
        double xli;
        double xni;
        double atime;
        /*
         * itegrator values for current d_atime_
         */
        struct IntegratorValues values_t;
    };
    
    void initialise();
    TECI findPositionSDP4(const double tsince) const;
    TECI findPositionSGP4(double tsince) const;
    TECI calculateFinalPositionVelocity(const double tsince, const double e, const double a,
            const double omega, const double xl, const double xnode, 
            const double xincl, const double xlcof, const double aycof, 
            const double x3thm1, const double x1mth2, const double x7thm1, 
            const double cosio, const double sinio) const;
    /**
     * Deep space initialisation
     */
    void deepSpaceInitialise(
            const double eosq, const double sinio, const double cosio,
            const double betao, const double theta2, const double betao2,
            const double xmdot, const double omgdot, const double xnodot);
    /*
     * Calculate lunar / solar terms
     */
    void deepSpaceCalculateLunarSolarTerms(
            const double tsince, double& pe, double& pinc,
            double& pl, double& pgh, double& ph) const;
    /**
     * Calculate lunar / solar periodics and apply
     */
    void deepSpacePeriodics(
            const double tsince, double& em, double& xinc,
            double& omgasm, double& xnodes, double& xll) const;
    /**
     * Deep space secular effects
     */
    void deepSpaceSecular(
            const double tsince,
            double& xll,
            double& omgasm,
            double& xnodes,
            double& em,
            double& xinc,
            double& xn) const;
    /**
     * Calculate dot terms
     * @param[in,out] values the integrator values
     */
    void deepSpaceCalcDotTerms(struct IntegratorValues& values) const;
    /**
     * Deep space integrator for time period of delt
     */
    void deepSpaceIntegrator(
            const double delt,
            const double step2,
            const struct IntegratorValues& values) const;
    void reset();

    /*
     * the constants used
     */
    struct CommonConstants common_consts_;
    struct NearSpaceConstants nearspace_consts_;
    struct DeepSpaceConstants deepspace_consts_;
    struct IntegratorConstants integrator_consts_;
    mutable struct IntegratorParams integrator_params_;

    /*
     * the orbit data
     */
    TOrbitalElements m_orbitalElements;

    /*
     * flags
     */
    bool m_use_simple_model;
    bool m_use_deep_space;

    static const struct TSGP4::CommonConstants Empty_CommonConstants;
    static const struct TSGP4::NearSpaceConstants Empty_NearSpaceConstants;
    static const struct TSGP4::DeepSpaceConstants Empty_DeepSpaceConstants;
    static const struct TSGP4::IntegratorConstants Empty_IntegratorConstants;
    static const struct TSGP4::IntegratorParams Empty_IntegratorParams;
};

TMATH_END_NAMESPACE

#endif