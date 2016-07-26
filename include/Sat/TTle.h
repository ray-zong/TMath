#ifndef TTLE_H
#define TTLE_H

#include "TGlobal.h"
#include "TDateTime.h"
#include "TTleException.h"

TMATH_BEGIN_NAMESPACE


/**
 * @brief Processes a two-line element set used to convey OrbitalElements.
 *
 * Used to extract the various raw fields from a two-line element set.
 */
class TMATH_EXPORT TTle
{
public:
    /**
     * @details Initialise given the two lines of a tle
     * @param[in] line_one Tle line one
     * @param[in] line_two Tle line two
     */
    TTle(const std::string& line_one,
            const std::string& line_two);

    /**
     * @details Initialise given the satellite name and the two lines of a tle
     * @param[in] name Satellite name
     * @param[in] line_one Tle line one
     * @param[in] line_two Tle line two
     */
    TTle(const std::string& name,
            const std::string& line_one,
            const std::string& line_two);

    /**
     * Copy constructor
     * @param[in] tle Tle object to copy from
     */
    TTle(const TTle& tle);

    /**
     * Get the satellite name
     * @returns the satellite name
     */
    std::string name() const;

    /**
     * Get the first line of the tle
     * @returns the first line of the tle
     */
    std::string line1() const;

    /**
     * Get the second line of the tle
     * @returns the second line of the tle
     */
    std::string line2() const;

    /**
     * Get the norad number
     * @returns the norad number
     */
    unsigned int noradNumber() const;

    /**
     * Get the international designator
     * @returns the international designator
     */
    std::string intDesignator() const;

    /**
     * Get the tle epoch
     * @returns the tle epoch
     */
    TDateTime epoch() const;

    /**
     * Get the first time derivative of the mean motion divided by two
     * @returns the first time derivative of the mean motion divided by two
     */
    double meanMotionDt2() const;

    /**
     * Get the second time derivative of mean motion divided by six
     * @returns the second time derivative of mean motion divided by six
     */
    double meanMotionDdt6() const;

    /**
     * Get the BSTAR drag term
     * @returns the BSTAR drag term
     */
    double bStar() const;

    /**
     * Get the inclination
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the inclination
     */
    double inclination(bool in_degrees) const;

    /**
     * Get the right ascension of the ascending node
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the right ascension of the ascending node
     */
    double rightAscendingNode(const bool in_degrees) const;

    /**
     * Get the eccentricity
     * @returns the eccentricity
     */
    double eccentricity() const;

    /**
     * Get the argument of perigee
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the argument of perigee
     */
    double argumentPerigee(const bool in_degrees) const;

    /**
     * Get the mean anomaly
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the mean anomaly
     */
    double meanAnomaly(const bool in_degrees) const;

    /**
     * Get the mean motion
     * @returns the mean motion (revolutions per day)
     */
    double meanMotion() const;

    /**
     * Get the orbit number
     * @returns the orbit number
     */
    unsigned int orbitNumber() const;

    /**
     * Get the expected tle line length
     * @returns the tle line length
     */
    static unsigned int lineLength();
    
    /**
     * Dump this object to a string
     * @returns string
     */
    std::string toString() const;

private:
    void initialize();
    static bool isValidLineLength(const std::string& str);
    void extractInteger(const std::string& str, unsigned int& val);
    void extractDouble(const std::string& str, int point_pos, double& val);
    void extractExponential(const std::string& str, double& val);

private:
    std::string m_name;
    std::string m_line_one;
    std::string m_line_two;

    std::string m_int_designator;
    TDateTime m_epoch;
    double m_mean_motion_dt2;
    double m_mean_motion_ddt6;
    double m_bstar;
    double m_inclination;
    double m_right_ascending_node;
    double m_eccentricity;
    double m_argument_perigee;
    double m_mean_anomaly;
    double m_mean_motion;
    unsigned int m_norad_number;
    unsigned int m_orbit_number;

    static const unsigned int TLE_LEN_LINE_DATA = 69;
    static const unsigned int TLE_LEN_LINE_NAME = 22;
};


inline std::ostream& operator<<(std::ostream& strm, const TTle& t)
{
    return strm << t.toString();
}

TMATH_END_NAMESPACE

#endif
