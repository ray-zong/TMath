#ifndef TORBITALELEMENTS_H
#define TORBITALELEMENTS_H

#include "TGlobal.h"
#include "TDateTime.h"
#include "TTle.h"

TMATH_BEGIN_NAMESPACE
/**
 * @brief The extracted orbital elements used by the SGP4 propagator.
 */
class TMATH_EXPORT TOrbitalElements
{
public:
    TOrbitalElements(const TTle& tle);

    /*
     * XMO
     */
    double meanAnomoly() const;

    /*
     * XNODEO
     */
    double ascendingNode() const;

    /*
     * OMEGAO
     */
    double argumentPerigee() const;

    /*
     * EO
     */
    double eccentricity() const;

    /*
     * XINCL
     */
    double inclination() const;

    /*
     * XNO
     */
    double meanMotion() const;

    /*
     * BSTAR
     */
    double bStar() const;

    /*
     * AODP
     */
    double recoveredSemiMajorAxis() const;

    /*
     * XNODP
     */
    double recoveredMeanMotion() const;

    /*
     * PERIGE
     */
    double perigee() const;

    /*
     * Period in minutes
     */
    double period() const;

    /*
     * EPOCH
     */
    TDateTime epoch() const;

private:
    double m_mean_anomoly;
    double m_ascending_node;
    double m_argument_perigee;
    double m_eccentricity;
    double m_inclination;
    double m_mean_motion;
    double m_bstar;
    double m_recovered_semi_major_axis;
    double m_recovered_mean_motion;
    double m_perigee;
    double m_period;
    TDateTime m_epoch;
};

TMATH_END_NAMESPACE

#endif
