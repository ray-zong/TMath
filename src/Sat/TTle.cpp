#include "TTle.h"

#include <locale> 

#include "TUtilities.h"

TMATH_BEGIN_NAMESPACE

static const unsigned int TLE1_COL_NORADNUM = 2;
static const unsigned int TLE1_LEN_NORADNUM = 5;
static const unsigned int TLE1_COL_INTLDESC_A = 9;
static const unsigned int TLE1_LEN_INTLDESC_A = 2;
//  static const unsigned int TLE1_COL_INTLDESC_B = 11;
static const unsigned int TLE1_LEN_INTLDESC_B = 3;
//  static const unsigned int TLE1_COL_INTLDESC_C = 14;
static const unsigned int TLE1_LEN_INTLDESC_C = 3;
static const unsigned int TLE1_COL_EPOCH_A = 18;
static const unsigned int TLE1_LEN_EPOCH_A = 2;
static const unsigned int TLE1_COL_EPOCH_B = 20;
static const unsigned int TLE1_LEN_EPOCH_B = 12;
static const unsigned int TLE1_COL_MEANMOTIONDT2 = 33;
static const unsigned int TLE1_LEN_MEANMOTIONDT2 = 10;
static const unsigned int TLE1_COL_MEANMOTIONDDT6 = 44;
static const unsigned int TLE1_LEN_MEANMOTIONDDT6 = 8;
static const unsigned int TLE1_COL_BSTAR = 53;
static const unsigned int TLE1_LEN_BSTAR = 8;
//  static const unsigned int TLE1_COL_EPHEMTYPE = 62;
//  static const unsigned int TLE1_LEN_EPHEMTYPE = 1;
//  static const unsigned int TLE1_COL_ELNUM = 64;
//  static const unsigned int TLE1_LEN_ELNUM = 4;

static const unsigned int TLE2_COL_NORADNUM = 2;
static const unsigned int TLE2_LEN_NORADNUM = 5;
static const unsigned int TLE2_COL_INCLINATION = 8;
static const unsigned int TLE2_LEN_INCLINATION = 8;
static const unsigned int TLE2_COL_RAASCENDNODE = 17;
static const unsigned int TLE2_LEN_RAASCENDNODE = 8;
static const unsigned int TLE2_COL_ECCENTRICITY = 26;
static const unsigned int TLE2_LEN_ECCENTRICITY = 7;
static const unsigned int TLE2_COL_ARGPERIGEE = 34;
static const unsigned int TLE2_LEN_ARGPERIGEE = 8;
static const unsigned int TLE2_COL_MEANANOMALY = 43;
static const unsigned int TLE2_LEN_MEANANOMALY = 8;
static const unsigned int TLE2_COL_MEANMOTION = 52;
static const unsigned int TLE2_LEN_MEANMOTION = 11;
static const unsigned int TLE2_COL_REVATEPOCH = 63;
static const unsigned int TLE2_LEN_REVATEPOCH = 5;

TTle::TTle( const std::string& line_one, const std::string& line_two )
    : m_line_one(line_one)
    , m_line_two(line_two)
{
    initialize();
}

TTle::TTle( const std::string& name, const std::string& line_one, const std::string& line_two )
    : m_name(name)
    , m_line_one(line_one)
    , m_line_two(line_two)
{
    initialize();
}

TTle::TTle( const TTle& tle )
{
    m_name = tle.m_name;
    m_line_one = tle.m_line_one;
    m_line_two = tle.m_line_two;

    m_norad_number = tle.m_norad_number;
    m_int_designator = tle.m_int_designator;
    m_epoch = tle.m_epoch;
    m_mean_motion_dt2 = tle.m_mean_motion_dt2;
    m_mean_motion_ddt6 = tle.m_mean_motion_ddt6;
    m_bstar = tle.m_bstar;
    m_inclination = tle.m_inclination;
    m_right_ascending_node = tle.m_right_ascending_node;
    m_eccentricity = tle.m_eccentricity;
    m_argument_perigee = tle.m_argument_perigee;
    m_mean_anomaly = tle.m_mean_anomaly;
    m_mean_motion = tle.m_mean_motion;
    m_orbit_number = tle.m_orbit_number;
}

std::string TTle::name() const
{
    return m_name;
}

std::string TTle::line1() const
{
    return m_line_one;
}

std::string TTle::line2() const
{
    return m_line_two;
}

unsigned int TTle::noradNumber() const
{
    return m_norad_number;
}

std::string TTle::intDesignator() const
{
    return m_int_designator;
}

TDateTime TTle::epoch() const
{
    return m_epoch;
}

double TTle::meanMotionDt2() const
{
    return m_mean_motion_dt2;
}

double TTle::meanMotionDdt6() const
{
    return m_mean_motion_ddt6;
}

double TTle::bStar() const
{
    return m_bstar;
}

double TTle::inclination( bool in_degrees ) const
{
    if (in_degrees)
    {
        return m_inclination;
    }
    else
    {
        return DegreesToRadians(m_inclination);
    }
}

double TTle::rightAscendingNode( const bool in_degrees ) const
{
    if (in_degrees)
    {
        return m_right_ascending_node;
    }
    else
    {
        return DegreesToRadians(m_right_ascending_node);
    }
}

double TTle::eccentricity() const
{
    return m_eccentricity;
}

double TTle::argumentPerigee( const bool in_degrees ) const
{
    if (in_degrees)
    {
        return m_argument_perigee;
    }
    else
    {
        return DegreesToRadians(m_argument_perigee);
    }
}

double TTle::meanAnomaly( const bool in_degrees ) const
{
    if (in_degrees)
    {
        return m_mean_anomaly;
    }
    else
    {
        return DegreesToRadians(m_mean_anomaly);
    }
}

double TTle::meanMotion() const
{
    return m_mean_motion;
}

unsigned int TTle::orbitNumber() const
{
    return m_orbit_number;
}

unsigned int TTle::lineLength()
{
    return TLE_LEN_LINE_DATA;
}

std::string TTle::toString() const
{
    std::stringstream ss;
    ss << std::right << std::fixed;
    ss << "Norad Number:         " << noradNumber() << std::endl;
    ss << "Int. Designator:      " << intDesignator() << std::endl;
    ss << "Epoch:                " << epoch() << std::endl;
    ss << "Orbit Number:         " << orbitNumber() << std::endl;
    ss << std::setprecision(8);
    ss << "Mean Motion Dt2:      ";
    ss << std::setw(12) << meanMotionDt2() << std::endl;
    ss << "Mean Motion Ddt6:     ";
    ss << std::setw(12) << meanMotionDdt6() << std::endl;
    ss << "Eccentricity:         ";
    ss << std::setw(12) << eccentricity() << std::endl;
    ss << "BStar:                ";
    ss << std::setw(12) << bStar() << std::endl;
    ss << "Inclination:          ";
    ss << std::setw(12) << inclination(true) << std::endl;
    ss << "Right Ascending Node: ";
    ss << std::setw(12) << rightAscendingNode(true) << std::endl;
    ss << "Argument Perigee:     ";
    ss << std::setw(12) << argumentPerigee(true) << std::endl;
    ss << "Mean Anomaly:         ";
    ss << std::setw(12) << meanAnomaly(true) << std::endl;
    ss << "Mean Motion:          ";
    ss << std::setw(12) << meanMotion() << std::endl;
    return ss.str();
}

/**
 * Initialise the tle object.
 * @exception TleException
 */
void TTle::initialize()
{
    if (!isValidLineLength(m_line_one))
    {
        throw TTleException("Invalid length for line one");
    }

    if (!isValidLineLength(m_line_two))
    {
        throw TTleException("Invalid length for line two");
    }

    if (m_line_one[0] != '1')
    {
        throw TTleException("Invalid line beginning for line one");
    }
        
    if (m_line_two[0] != '2')
    {
        throw TTleException("Invalid line beginning for line two");
    }

    unsigned int sat_number_1;
    unsigned int sat_number_2;

    extractInteger(m_line_one.substr(TLE1_COL_NORADNUM,
                TLE1_LEN_NORADNUM), sat_number_1);
    extractInteger(m_line_two.substr(TLE2_COL_NORADNUM,
                TLE2_LEN_NORADNUM), sat_number_2);

    if (sat_number_1 != sat_number_2)
    {
        throw TTleException("Satellite numbers do not match");
    }

    m_norad_number = sat_number_1;

    if (m_name.empty())
    {
        m_name = m_line_one.substr(TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM);
    }

    m_int_designator = m_line_one.substr(TLE1_COL_INTLDESC_A,
            TLE1_LEN_INTLDESC_A + TLE1_LEN_INTLDESC_B + TLE1_LEN_INTLDESC_C);

    unsigned int year = 0;
    double day = 0.0;

    extractInteger(m_line_one.substr(TLE1_COL_EPOCH_A,
                TLE1_LEN_EPOCH_A), year);
    extractDouble(m_line_one.substr(TLE1_COL_EPOCH_B,
                TLE1_LEN_EPOCH_B), 4, day);
    extractDouble(m_line_one.substr(TLE1_COL_MEANMOTIONDT2,
                TLE1_LEN_MEANMOTIONDT2), 2, m_mean_motion_dt2);
    extractExponential(m_line_one.substr(TLE1_COL_MEANMOTIONDDT6,
                TLE1_LEN_MEANMOTIONDDT6), m_mean_motion_ddt6);
    extractExponential(m_line_one.substr(TLE1_COL_BSTAR,
                TLE1_LEN_BSTAR), m_bstar);

    /*
     * line 2
     */
    extractDouble(m_line_two.substr(TLE2_COL_INCLINATION,
                TLE2_LEN_INCLINATION), 4, m_inclination);
    extractDouble(m_line_two.substr(TLE2_COL_RAASCENDNODE,
                TLE2_LEN_RAASCENDNODE), 4, m_right_ascending_node);
    extractDouble(m_line_two.substr(TLE2_COL_ECCENTRICITY,
                TLE2_LEN_ECCENTRICITY), -1, m_eccentricity);
    extractDouble(m_line_two.substr(TLE2_COL_ARGPERIGEE,
                TLE2_LEN_ARGPERIGEE), 4, m_argument_perigee);
    extractDouble(m_line_two.substr(TLE2_COL_MEANANOMALY,
                TLE2_LEN_MEANANOMALY), 4, m_mean_anomaly);
    extractDouble(m_line_two.substr(TLE2_COL_MEANMOTION,
                TLE2_LEN_MEANMOTION), 3, m_mean_motion);
    extractInteger(m_line_two.substr(TLE2_COL_REVATEPOCH,
                TLE2_LEN_REVATEPOCH), m_orbit_number);
    
    if (year < 57)
        year += 2000;
    else
        year += 1900;

    m_epoch = TDateTime(year, day);
}

/**
 * Check 
 * @param str The string to check
 * @returns Whether true of the string has a valid length
 */
bool TTle::isValidLineLength(const std::string& str)
{
    return str.length() == lineLength() ? true : false;
}

/**
 * Convert a string containing an integer
 * @param[in] str The string to convert
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void TTle::extractInteger(const std::string& str, unsigned int& val)
{
    bool found_digit = false;
    unsigned int temp = 0;

    for (std::string::const_iterator i = str.begin(); i != str.end(); ++i)
    {
        if (isdigit(*i))
        {
            found_digit = true;
            temp = (temp * 10) + (static_cast<unsigned char>(*i) - '0');
        }
        else if (found_digit)
        {
            throw TTleException("Unexpected non digit");
        }
        else if (*i != ' ')
        {
            throw TTleException("Invalid character");
        }
    }

    if (!found_digit)
    {
        val = 0;
    }
    else
    {
        val = temp;
    }
}

/**
 * Convert a string containing an double
 * @param[in] str The string to convert
 * @param[in] point_pos The position of the decimal point. (-1 if none)
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void TTle::extractDouble(const std::string& str, int point_pos, double& val)
{
    std::string temp;
    bool found_digit = false;

    for (std::string::const_iterator i = str.begin(); i != str.end(); ++i)
    {
        /*
         * integer part
         */
        if (point_pos >= 0 && i < str.begin() + point_pos - 1)
        {
            bool done = false;

            if (i == str.begin())
            {
                if(*i == '-' || *i == '+')
                {
                    /*
                     * first character could be signed
                     */
                    temp += *i;
                    done = true;
                }
            }

            if (!done)
            {
                if (isdigit(*i))
                {
                    found_digit = true;
                    temp += *i;
                }
                else if (found_digit)
                {
                    throw TTleException("Unexpected non digit");
                }
                else if (*i != ' ')
                {
                    throw TTleException("Invalid character");
                }
            }
        }
        /*
         * decimal point
         */
        else if (point_pos >= 0 && i == str.begin() + point_pos - 1)
        {
            if (temp.length() == 0)
            {
                /*
                 * integer part is blank, so add a '0'
                 */
                temp += '0';
            }

            if (*i == '.')
            {
                /*
                 * decimal point found
                 */
                temp += *i;
            }
            else
            {
                throw TTleException("Failed to find decimal point");
            }
        }
        /*
         * fraction part
         */
        else
        {
            if (i == str.begin() && point_pos == -1)
            {
                /*
                 * no decimal point expected, add 0. beginning
                 */
                temp += '0';
                temp += '.';
            }
            
            /*
             * should be a digit
             */
            if (isdigit(*i))
            {
                temp += *i;
            }
            else
            {
                throw TTleException("Invalid digit");
            }
        }
    }

    if (!FromString<double>(temp, val))
    {
        throw TTleException("Failed to convert value to double");
    }
}

/**
 * Convert a string containing an exponential
 * @param[in] str The string to convert
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void TTle::extractExponential(const std::string& str, double& val)
{
    std::string temp;

    for (std::string::const_iterator i = str.begin(); i != str.end(); ++i)
    {
        if (i == str.begin())
        {
            if (*i == '-' || *i == '+' || *i == ' ')
            {
                if (*i == '-')
                {
                    temp += *i;
                }
                temp += '0';
                temp += '.';
            }
            else
            {
                throw TTleException("Invalid sign");
            }
        }
        else if (i == str.end() - 2)
        {
            if (*i == '-' || *i == '+')
            {
                temp += 'e';
                temp += *i;
            }
            else
            {
                throw TTleException("Invalid exponential sign");
            }
        }
        else
        {
            if (isdigit(*i))
            {
                temp += *i;
            }
            else
            {
                throw TTleException("Invalid digit");
            }
        }
    }

    if (!FromString<double>(temp, val))
    {
        throw TTleException("Failed to convert value to double");
    }
}

TMATH_END_NAMESPACE

