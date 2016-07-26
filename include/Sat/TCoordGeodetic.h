#ifndef TCOORDGEODETIC_H_
#define TCOORDGEODETIC_H_

#include <string>
#include <sstream>
#include <iomanip>

#include "TGlobal.h"

TMATH_BEGIN_NAMESPACE

/**
 * @brief Stores a geodetic location (latitude, longitude, altitude).
 *
 * Internally the values are stored in radians and kilometres.
 */
class TMATH_EXPORT TCoordGeodetic
{
public:
    /**
     * Default constructor
     */
    TCoordGeodetic();

    /**
     * Constructor
     * @param[in] lat the latitude (degrees by default)
     * @param[in] lon the longitude (degrees by default)
     * @param[in] alt the altitude in kilometers
     * @param[in] is_radians whether the latitude/longitude is in radians
     */
    TCoordGeodetic(double lat, double lon, double alt,
            bool is_radians = false);

    /**
     * Copy constructor
     * @param[in] geo object to copy from
     */
    TCoordGeodetic(const TCoordGeodetic& geo);

    /**
     * Assignment operator
     * @param[in] geo object to copy from
     */
    TCoordGeodetic& operator=(const TCoordGeodetic& geo);

    /**
     * Equality operator
     * @param[in] geo the object to compare with
     * @returns whether the object is equal
     */
    bool operator==(const TCoordGeodetic& geo) const;

    /**
     * Inequality operator
     * @param[in] geo the object to compare with
     * @returns whether the object is not equal
     */
    bool operator!=(const TCoordGeodetic& geo) const;

    /**
     * Dump this object to a string
     * @returns string
     */
    std::string toString() const;

    double latitude() const;

    double longitude() const;

    double altitude() const;

private:
    bool isEqual(const TCoordGeodetic& geo) const;

    /** latitude in radians (-PI >= latitude < PI) */
    double m_latitude;
    /** latitude in radians (-PI/2 >= latitude <= PI/2) */
    double m_longitude;
    /** altitude in kilometers */
    double m_altitude;
};

/**
 * Dump a Coordgeodetic to a stream
 * @param[in,out] strm stream to output to
 * @param[in] g the CoordGeodetic to print
 */
inline std::ostream& operator<<(std::ostream& strm, const TCoordGeodetic& g)
{
    return strm << g.toString();
}

TMATH_END_NAMESPACE

#endif
