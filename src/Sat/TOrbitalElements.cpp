#include "TOrbitalElements.h"

#include "TConstants.h"

TMATH_BEGIN_NAMESPACE

TOrbitalElements::TOrbitalElements(const TTle& tle)
{
    /*
     * extract and format tle data
     */
    m_mean_anomoly = tle.meanAnomaly(false);
    m_ascending_node = tle.rightAscendingNode(false);
    m_argument_perigee = tle.argumentPerigee(false);
    m_eccentricity = tle.eccentricity();
    m_inclination = tle.inclination(false);
    m_mean_motion = tle.meanMotion() * kTWOPI / kMINUTES_PER_DAY;
    m_bstar = tle.bStar();
    m_epoch = tle.epoch();

    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    const double a1 = pow(kXKE / meanMotion(), kTWOTHIRD);
    const double cosio = cos(inclination());
    const double theta2 = cosio * cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double eosq = eccentricity() * eccentricity();
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);
    const double temp = (1.5 * kCK2) * x3thm1 / (betao * betao2);
    const double del1 = temp / (a1 * a1);
    const double a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)));
    const double del0 = temp / (a0 * a0);

    m_recovered_mean_motion = meanMotion() / (1.0 + del0);
    /*
     * alternative way to calculate
     * doesnt affect final results
     * recovered_semi_major_axis_ = pow(XKE / RecoveredMeanMotion(), TWOTHIRD);
     */
    m_recovered_semi_major_axis = a0 / (1.0 - del0);

    /*
     * find perigee and period
     */
    m_perigee = (recoveredSemiMajorAxis() * (1.0 - eccentricity()) - kAE) * kXKMPER;
    m_period = kTWOPI / recoveredMeanMotion();
}

double TOrbitalElements::meanAnomoly() const
{
    return m_mean_anomoly;
}

double TOrbitalElements::ascendingNode() const
{
    return m_ascending_node;
}

double TOrbitalElements::argumentPerigee() const
{
    return m_argument_perigee;
}

double TOrbitalElements::eccentricity() const
{
    return m_eccentricity;
}

double TOrbitalElements::inclination() const
{
    return m_inclination;
}

double TOrbitalElements::meanMotion() const
{
    return m_mean_motion;
}

double TOrbitalElements::bStar() const
{
    return m_bstar;
}

double TOrbitalElements::recoveredSemiMajorAxis() const
{
    return m_recovered_semi_major_axis;
}

double TOrbitalElements::recoveredMeanMotion() const
{
    return m_recovered_mean_motion;
}

double TOrbitalElements::perigee() const
{
    return m_perigee;
}

double TOrbitalElements::period() const
{
    return m_period;
}

TDateTime TOrbitalElements::epoch() const
{
    return m_epoch;
}

TMATH_END_NAMESPACE
