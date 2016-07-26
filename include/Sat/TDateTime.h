#ifndef TDATETIME_H
#define TDATETIME_H

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdint>

#include "TGlobal.h"
#include "TTimeSpan.h"

TMATH_BEGIN_NAMESPACE

static const int kDaysInMonth[2][13] = {
        //  1   2   3   4   5   6   7   8   9   10  11  12
        {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
        {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};
static const int kCumulDaysInMonth[2][13] = {
    //  1  2   3   4   5    6    7    8    9    10   11   12
    {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
    {0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
};

/**
 * @brief Represents an instance in time.
 */

class TMATH_EXPORT TDateTime
{
public:
    /**
     * Default contructor
     * Initialise to 0001/01/01 00:00:00.000000
     */
    TDateTime();

    /**
     * Constructor
     * @param[in] ticks raw tick value
     */
    TDateTime(int64_t ticks);

    /**
     * Constructor
     * @param[in] year the year
     * @param[in] doy the day of the year
     */
    TDateTime(unsigned int year, double doy);

    /**
     * Constructor
     * @param[in] year the year
     * @param[in] month the month
     * @param[in] day the day
     */
    TDateTime(int year, int month, int day);

    /**
     * Constructor
     * @param[in] year the year
     * @param[in] month the month
     * @param[in] day the day
     * @param[in] hour the hour
     * @param[in] minute the minute
     * @param[in] second the second
     */
    TDateTime(int year, int month, int day, 
        int hour, int minute, int second);

    /**
     * Constructor
     * @param[in] year the year
     * @param[in] month the month
     * @param[in] day the day
     * @param[in] hour the hour
     * @param[in] minute the minute
     * @param[in] second the second
     * @param[in] microsecond the microsecond
     */
    void initialise(int year,
            int month,
            int day,
            int hour,
            int minute,
            int second,
            int microsecond);

    /**
     * Find whether a year is a leap year
     * @param[in] year the year to check
     * @returns whether the year is a leap year
     */
    static bool isLeapYear(int year);

    /**
     * Checks whether the given year is valid
     * @param[in] year the year to check
     * @returns whether the year is valid
     */
    static bool isValidYear(int year);

    /**
     * Check whether the year/month is valid
     * @param[in] year the year to check
     * @param[in] month the month to check
     * @returns whether the year/month is valid
     */
    static bool isValidYearMonth(int year, int month);

    /**
     * Check whether the year/month/day is valid
     * @param[in] year the year to check
     * @param[in] month the month to check
     * @param[in] day the day to check
     * @returns whether the year/month/day is valid
     */
    static bool isValidYearMonthDay(int year, int month, int day);

    /**
     * Find the number of days in a month given the year/month
     * @param[in] year the year
     * @param[in] month the month
     * @returns the days in the given month
     */
    static int daysInMonth(int year, int month);

    /**
     * Find the day of the year given the year/month/day
     * @param[in] year the year
     * @param[in] month the month
     * @param[in] day the day
     * @returns the day of the year
     */
    int dayOfYear(int year, int month, int day) const;

    /**
     *
     */
    double absoluteDays(unsigned int year, double doy) const;

    int absoluteDays(int year, int month, int day) const;

    TTimeSpan timeOfDay() const;

    int dayOfWeek() const;

    bool equals(const TDateTime& dt) const;

    int compare(const TDateTime& dt) const;

    TDateTime addYears(const int years) const;

    TDateTime addMonths(const int months) const;

    /**
     * Add a TimeSpan to this DateTime
     * @param[in] t the TimeSpan to add
     * @returns a DateTime which has the given TimeSpan added
     */
    TDateTime add(const TTimeSpan& t) const;

    TDateTime addDays(const double days) const;

    TDateTime addHours(const double hours) const;

    TDateTime addMinutes(const double minutes) const;

    TDateTime addSeconds(const double seconds) const;

    TDateTime addMicroseconds(const double microseconds) const;

    TDateTime addTicks(int64_t ticks) const;

    /**
     * Get the number of ticks
     * @returns the number of ticks
     */
    int64_t ticks() const;

    void fromTicks(int& year, int& month, int& day) const;

    int year() const;

    int month() const;

    int day() const;

    /**
     * Hour component
     * @returns the hour component
     */
    int hour() const;

    /**
     * Minute component
     * @returns the minute component
     */
    int minute() const;

    /**
     * Second component
     * @returns the Second component
     */
    int second() const;

    /**
     * Microsecond component
     * @returns the microsecond component
     */
    int microsecond() const;

    /**
     * Convert to a julian date
     * @returns the julian date
     */
    double toJulian() const;

    /**
     * Convert to greenwich sidereal time
     * @returns the greenwich sidereal time
     */
    double toGreenwichSiderealTime() const;

    /**
     * Convert to local mean sidereal time (GMST plus the observer's longitude)
     * @param[in] lon observers longitude
     * @returns the local mean sidereal time
     */
    double toLocalMeanSiderealTime(const double lon) const;

    std::string toString() const;

    private:
        int64_t m_encoded;
};

inline std::ostream& operator<<(std::ostream& strm, const TDateTime& dt)
{
    return strm << dt.toString();
}

inline TDateTime operator+(const TDateTime& dt, TTimeSpan ts)
{
    return TDateTime(dt.ticks() + ts.ticks());
}

inline TDateTime operator-(const TDateTime& dt, const TTimeSpan& ts)
{
    return TDateTime(dt.ticks() - ts.ticks());
}

inline TTimeSpan operator-(const TDateTime& dt1, const TDateTime& dt2)
{
    return TTimeSpan(dt1.ticks() - dt2.ticks());
}

inline bool operator==(const TDateTime& dt1, const TDateTime& dt2)
{
    return dt1.equals(dt2);
}

inline bool operator>(const TDateTime& dt1, const TDateTime& dt2)
{
    return (dt1.compare(dt2) > 0);
}

inline bool operator>=(const TDateTime& dt1, const TDateTime& dt2)
{
    return (dt1.compare(dt2) >= 0);
}

inline bool operator!=(const TDateTime& dt1, const TDateTime& dt2)
{
    return !dt1.equals(dt2);
}

inline bool operator<(const TDateTime& dt1, const TDateTime& dt2)
{
    return (dt1.compare(dt2) < 0);
}

inline bool operator<=(const TDateTime& dt1, const TDateTime& dt2)
{
    return (dt1.compare(dt2) <= 0);
}

TMATH_END_NAMESPACE

#endif
