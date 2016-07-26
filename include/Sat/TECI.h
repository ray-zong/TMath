#ifndef TECI_H
#define TECI_H

#include "TGlobal.h"
#include "TCoordGeodetic.h"
#include "TVector.h"
#include "TDateTime.h"

TMATH_BEGIN_NAMESPACE

template class TMATH_EXPORT TVector<double, 4>;

/**
 * @brief Stores an Earth-centered inertial position for a particular time.
 */

class TMATH_EXPORT TECI
{
public:
    /**
     * @param[in] dt the date to be used for this position
     * @param[in] latitude the latitude in degrees
     * @param[in] longitude the longitude in degrees
     * @param[in] altitude the altitude in kilometers
     */
    TECI(const TDateTime& dt,
            const double latitude,
            const double longitude,
            const double altitude);

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] geo the position
     */
    TECI(const TDateTime& dt, const TCoordGeodetic& geo);

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] position the position
     */
    TECI(const TDateTime &dt, const TVector<double, 3> &position);

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] position the position
     * @param[in] velocity the velocity
     */
    TECI(const TDateTime &dt, const TVector<double, 3> &position, const TVector<double, 3> &velocity);

    /**
     * Equality operator
     * @param dt the date to compare
     * @returns true if the object matches
     */
    bool operator==(const TDateTime& dt) const;

    /**
     * Inequality operator
     * @param dt the date to compare
     * @returns true if the object doesn't match
     */
    bool operator!=(const TDateTime& dt) const;

    /**
     * Update this object with a new date and geodetic position
     * @param dt new date
     * @param geo new geodetic position
     */
    void update(const TDateTime& dt, const TCoordGeodetic& geo);

    /**
     * @brief 设置位置
     */
    void setPosition(const TVector<double, 3> &pos);

    /**
     * @brief 设置速度
     */
    void setVelocity(const TVector<double, 3> &vel);

    /**
     * @returns the position
     */
    TVector<double, 3> position() const{return m_position;}

    /**
     * @returns the velocity
     */
    TVector<double, 3> velocity() const{return m_velocity;}

    /**
     * @returns the date
     */
    TDateTime getDateTime() const;

    /**
     * @returns the position in geodetic form
     */
    TCoordGeodetic toGeodetic() const;

private:
    /**
     * Converts a DateTime and Geodetic position to Eci coordinates
     * @param[in] dt the date
     * @param[in] geo the geodetic position
     */
    void toECI(const TDateTime& dt, const TCoordGeodetic& geo);

    TDateTime m_dt;
    TVector<double, 3> m_position;
    TVector<double, 3> m_velocity;
};

TMATH_END_NAMESPACE


#endif
