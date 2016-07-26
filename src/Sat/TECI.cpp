#include "TECI.h"

#include "TUtilities.h"
#include "TConstants.h"

TMATH_BEGIN_NAMESPACE

TECI::TECI( const TDateTime& dt, const double latitude, const double longitude, const double altitude )
{
    toECI(dt, TCoordGeodetic(latitude, longitude, altitude));
}

TECI::TECI( const TDateTime& dt, const TCoordGeodetic& geo )
{
    toECI(dt, geo);
}

TECI::TECI( const TDateTime &dt, const TVector<double, 3> &position )
    : m_dt(dt),
    m_position(position)
{
}

TECI::TECI( const TDateTime &dt, const TVector<double, 3> &position, const TVector<double, 3> &velocity )
    : m_dt(dt),
    m_position(position),
    m_velocity(velocity)
{
}

bool TECI::operator==( const TDateTime& dt ) const
{
    return m_dt == dt;
}

bool TECI::operator!=( const TDateTime& dt ) const
{
    return m_dt != dt;
}

void TECI::update( const TDateTime& dt, const TCoordGeodetic& geo )
{
    toECI(dt, geo);
}

void TECI::setPosition(const TVector<double, 3> &pos)
{
    m_position = pos;
}


void TECI::setVelocity(const TVector<double, 3> &vel)
{
    m_velocity = vel;
}

//TVector<double, 3> TECI::position() const
//{
//    return m_position;
//}

//TVector<double, 3> TECI::velocity() const
//{
//    return m_velocity;
//}

TDateTime TECI::getDateTime() const
{
    return m_dt;
}

TCoordGeodetic TECI::toGeodetic() const
{
    const double theta = AcTan(m_position.y(), m_position.x());

    const double lon = WrapNegPosPI(theta
        - m_dt.toGreenwichSiderealTime());

    const double r = sqrt((m_position.x() * m_position.x())
        + (m_position.y() * m_position.y()));

    static const double e2 = kF * (2.0 - kF);

    double lat = AcTan(m_position.z(), r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;

    do
    {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = AcTan(m_position.z() + kXKMPER * c * e2 * sinphi, r);
        cnt++;
    }
    while (fabs(lat - phi) >= 1e-10 && cnt < 10);

    const double alt = r / cos(lat) - kXKMPER * c;

    return TCoordGeodetic(lat, lon, alt, true);
}

void TECI::toECI( const TDateTime& dt, const TCoordGeodetic& geo )
{
    /*
     * set date
     */
    m_dt = dt;

    static const double mfactor = kTWOPI * (kOMEGA_E / kSECONDS_PER_DAY);

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    const double theta = m_dt.toLocalMeanSiderealTime(geo.longitude());

    /*
     * take into account earth flattening
     */
    const double c = 1.0
        / sqrt(1.0 + kF * (kF - 2.0) * pow(sin(geo.latitude()), 2.0));
    const double s = pow(1.0 - kF, 2.0) * c;
    const double achcp = (kXKMPER * c + geo.altitude()) * cos(geo.latitude());

    /*
     * X position in km
     * Y position in km
     * Z position in km
     * W magnitude in km
     */
    m_position[0] = achcp * cos(theta);
    m_position[1] = achcp * sin(theta);
    m_position[2] = (kXKMPER * s + geo.altitude()) * sin(geo.latitude());

    /*
     * X velocity in km/s
     * Y velocity in km/s
     * Z velocity in km/s
     * W magnitude in km/s
     */
    m_velocity[0] = -mfactor * m_position.y();
    m_velocity[1] = mfactor * m_position.x();
    m_velocity[2] = 0.0;
}

TMATH_END_NAMESPACE
