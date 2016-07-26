#include "TCoordGeodetic.h"

#include "TUtilities.h"

TMATH_BEGIN_NAMESPACE

TCoordGeodetic::TCoordGeodetic()
    : m_latitude(0.0),
    m_longitude(0.0),
    m_altitude(0.0)
{
}

TCoordGeodetic::TCoordGeodetic( double lat, double lon, double alt, bool is_radians /*= false*/ )
{
    if (is_radians)
    {
        m_latitude = lat;
        m_longitude = lon;
    }
    else
    {
        m_latitude = DegreesToRadians(lat);
        m_longitude = DegreesToRadians(lon);
    }
    m_altitude = alt;
}

TCoordGeodetic::TCoordGeodetic( const TCoordGeodetic& geo )
{
    m_latitude = geo.m_latitude;
    m_longitude = geo.m_longitude;
    m_altitude = geo.m_altitude;
}

TCoordGeodetic& TCoordGeodetic::operator=( const TCoordGeodetic& geo )
{
    if (this != &geo)
    {
        m_latitude = geo.m_latitude;
        m_longitude = geo.m_longitude;
        m_altitude = geo.m_altitude;
    }
    return *this;
}

bool TCoordGeodetic::operator==( const TCoordGeodetic& geo ) const
{
    return isEqual(geo);
}

bool TCoordGeodetic::operator!=( const TCoordGeodetic& geo ) const
{
    return !isEqual(geo);
}

std::string TCoordGeodetic::toString() const
{
    std::stringstream ss;
    ss << std::right << std::fixed << std::setprecision(3);
    ss << "Lat: " << std::setw(7) << RadiansToDegrees(m_latitude);
    ss << ", Lon: " << std::setw(7) << RadiansToDegrees(m_longitude);
    ss << ", Alt: " << std::setw(9) << m_altitude;
    return ss.str();
}

double TCoordGeodetic::latitude() const
{
    return m_latitude;
}

double TCoordGeodetic::longitude() const
{
    return m_longitude;
}

double TCoordGeodetic::altitude() const
{
    return m_altitude;
}

bool TCoordGeodetic::isEqual( const TCoordGeodetic& geo ) const
{
    bool equal = false;
    if (m_latitude == geo.m_latitude &&
        m_longitude == geo.m_longitude &&
        m_altitude == geo.m_altitude)
    {
        equal = false;
    }
    return equal;
}

TMATH_END_NAMESPACE
