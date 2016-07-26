#include "TDateTime.h"

#include <ctime>

#include "TUtilities.h"

TMATH_BEGIN_NAMESPACE

TDateTime::TDateTime()
{
    initialise(1, 1, 1, 0, 0, 0, 0);
}

TDateTime::TDateTime( int64_t ticks )
    : m_encoded(ticks)
{
}

TDateTime::TDateTime( unsigned int year, double doy )
{
    m_encoded = TTimeSpan(
        static_cast<int64_t>(absoluteDays(year, doy) * TicksPerDay)).ticks();
}

TDateTime::TDateTime( int year, int month, int day )
{
    initialise(year, month, day, 0, 0, 0, 0);
}

TDateTime::TDateTime( int year, int month, int day, int hour, int minute, int second )
{
    initialise(year, month, day, hour, minute, second, 0);
}

void TDateTime::initialise( int year, int month, int day, int hour, int minute, int second, int microsecond )
{
    if (!isValidYearMonthDay(year, month, day) ||
        hour < 0 || hour > 23 ||
        minute < 0 || minute > 59 ||
        second < 0 || second > 59 ||
        microsecond < 0 || microsecond > 999999)
    {
        throw 1;
    }
    m_encoded = TTimeSpan(
        absoluteDays(year, month, day),
        hour,
        minute,
        second,
        microsecond).ticks();
}

bool TDateTime::isLeapYear( int year )
{
    if (!isValidYear(year))
    {
        throw 1;
    }

    return (((year % 4) == 0 && (year % 100) != 0) || (year % 400) == 0);
}

bool TDateTime::isValidYear( int year )
{
    bool valid = true;
    if (year < 1 || year > 9999)
    {
        valid = false;
    }
    return valid;
}

bool TDateTime::isValidYearMonth( int year, int month )
{
    bool valid = true;
    if (isValidYear(year))
    {
        if (month < 1 || month > 12)
        {
            valid = false;
        }
    }
    else 
    {
        valid = false;
    }
    return valid;
}

bool TDateTime::isValidYearMonthDay( int year, int month, int day )
{
    bool valid = true;
    if (isValidYearMonth(year, month))
    {
        if (day < 1 || day > daysInMonth(year, month))
        {
            valid = false;
        }
    }
    else
    {
        valid = false;
    }
    return valid;
}

int TDateTime::daysInMonth( int year, int month )
{
    if (!isValidYearMonth(year, month))
    {
        throw 1;
    }

    const int* daysInMonthPtr;

    if (isLeapYear(year))
    {
        daysInMonthPtr = kDaysInMonth[1];
    }
    else
    {
        daysInMonthPtr = kDaysInMonth[0];
    }

    return daysInMonthPtr[month];
}

int TDateTime::dayOfYear( int year, int month, int day ) const
{
    if (!isValidYearMonthDay(year, month, day))
    {
        throw 1;
    }

    int daysThisYear = day;

    if (isLeapYear(year))
    {
        daysThisYear += kCumulDaysInMonth[1][month];
    }
    else
    {
        daysThisYear += kCumulDaysInMonth[0][month];
    }

    return daysThisYear;
}

double TDateTime::absoluteDays( unsigned int year, double doy ) const
{
    int64_t previousYear = year - 1;

        /*
         * + days in previous years ignoring leap days
         * + Julian leap days before this year
         * - minus prior century years
         * + plus prior years divisible by 400 days
         */
        int64_t daysSoFar = 365 * previousYear
            + previousYear / 4LL
            - previousYear / 100LL
            + previousYear / 400LL;

        return static_cast<double>(daysSoFar) + doy - 1.0;
}

int TDateTime::absoluteDays( int year, int month, int day ) const
{
    int previousYear = year - 1;

        /*
         * days this year (0 - ...)
         * + days in previous years ignoring leap days
         * + Julian leap days before this year
         * - minus prior century years
         * + plus prior years divisible by 400 days
         */
        int result = dayOfYear(year, month, day) - 1
            + 365 * previousYear
            + previousYear / 4
            - previousYear / 100
            + previousYear / 400;

        return result;
}

TTimeSpan TDateTime::timeOfDay() const
{
    return TTimeSpan(ticks() % TicksPerDay);
}

int TDateTime::dayOfWeek() const
{
    /*
         * The fixed day 1 (January 1, 1 Gregorian) is Monday.
         * 0 Sunday
         * 1 Monday
         * 2 Tuesday
         * 3 Wednesday
         * 4 Thursday
         * 5 Friday
         * 6 Saturday
         */
        return static_cast<int>(((m_encoded / TicksPerDay) + 1LL) % 7LL);
}

bool TDateTime::equals( const TDateTime& dt ) const
{
    return (m_encoded == dt.m_encoded);
}

int TDateTime::compare( const TDateTime& dt ) const
{
    int ret = 0;

    if (m_encoded < dt.m_encoded)
    {
        return -1;
    }
    else if (m_encoded > dt.m_encoded)
    {
        return 1;
    }

    return ret;
}

TDateTime TDateTime::addYears( const int years ) const
{
    return addMonths(years * 12);
}

TDateTime TDateTime::addMonths( const int months ) const
{
    int year;
    int month;
    int day;
    fromTicks(year, month, day);
    month += months % 12;
    year += months / 12;

    if (month < 1)
    {
        month += 12;
        --year;
    }
    else if (month > 12)
    {
        month -= 12;
        ++year;
    }

    int maxday = daysInMonth(year, month);
    day = std::min(day, maxday);

    return TDateTime(year, month, day).add(timeOfDay());
}

TDateTime TDateTime::add( const TTimeSpan& t ) const
{
    return addTicks(t.ticks());
}

TDateTime TDateTime::addDays( const double days ) const
{
    return addMicroseconds(days * 86400000000.0);
}

TDateTime TDateTime::addHours( const double hours ) const
{
    return addMicroseconds(hours * 3600000000.0);
}

TDateTime TDateTime::addMinutes( const double minutes ) const
{
    return addMicroseconds(minutes * 60000000.0);
}

TDateTime TDateTime::addSeconds( const double seconds ) const
{
    return addMicroseconds(seconds * 1000000.0);
}

TDateTime TDateTime::addMicroseconds( const double microseconds ) const
{
    int64_t ticks = static_cast<int64_t>(microseconds * TicksPerMicrosecond);
    return addTicks(ticks);
}

TDateTime TDateTime::addTicks( int64_t ticks ) const
{
    return TDateTime(m_encoded + ticks);
}

std::int64_t TDateTime::ticks() const
{
    return m_encoded;
}

void TDateTime::fromTicks( int& year, int& month, int& day ) const
{
    int totalDays = static_cast<int>(m_encoded / TicksPerDay);
        
        /*
         * number of 400 year cycles
         */
        int num400 = totalDays / 146097;
        totalDays -= num400 * 146097;
        /*
         * number of 100 year cycles
         */
        int num100 = totalDays / 36524;
        if (num100 == 4)
        {
            /*
             * last day of the last leap century
             */
            num100 = 3;
        }
        totalDays -= num100 * 36524;
        /*
         * number of 4 year cycles
         */
        int num4 = totalDays / 1461;
        totalDays -= num4 * 1461;
        /*
         * number of years
         */
        int num1 = totalDays / 365;
        if (num1 == 4)
        {
            /*
             * last day of the last leap olympiad
             */
            num1 = 3;
        }
        totalDays -= num1 * 365;

        /*
         * find year
         */
        year = (num400 * 400) + (num100 * 100) + (num4 * 4) + num1 + 1;
        
        /*
         * convert day of year to month/day
         */
        const int* daysInMonthPtr;
        if (isLeapYear(year))
        {
            daysInMonthPtr = kDaysInMonth[1];
        }
        else
        {
            daysInMonthPtr = kDaysInMonth[0];
        }

        month = 1;
        while (totalDays >= daysInMonthPtr[month] && month <= 12)
        {
            totalDays -= daysInMonthPtr[month++];
        }

        day = totalDays + 1;
}

int TDateTime::year() const
{
    int year;
    int month;
    int day;
    fromTicks(year, month, day);
    return year;
}

int TDateTime::month() const
{
    int year;
    int month;
    int day;
    fromTicks(year, month, day);
    return month;
}

int TDateTime::day() const
{
    int year;
    int month;
    int day;
    fromTicks(year, month, day);
    return day;
}

int TDateTime::hour() const
{
    return static_cast<int>(m_encoded % TicksPerDay / TicksPerHour);
}

int TDateTime::minute() const
{
    return static_cast<int>(m_encoded % TicksPerHour / TicksPerMinute);
}

int TDateTime::second() const
{
    return static_cast<int>(m_encoded % TicksPerMinute / TicksPerSecond);
}

int TDateTime::microsecond() const
{
    return static_cast<int>(m_encoded % TicksPerSecond / TicksPerMicrosecond);
}

double TDateTime::toJulian() const
{
    TTimeSpan ts = TTimeSpan(ticks());
    return ts.totalDays() + 1721425.5;
}

double TDateTime::toGreenwichSiderealTime() const
{
    // t = Julian centuries from 2000 Jan. 1 12h UT1
    const double t = (toJulian() - 2451545.0) / 36525.0;

    // Rotation angle in arcseconds
    double theta = 67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * t
        + 0.093104 * t * t
        - 0.0000062 * t * t * t;

    // 360.0 / 86400.0 = 1.0 / 240.0
    return WrapTwoPI(DegreesToRadians(theta / 240.0));
}

double TDateTime::toLocalMeanSiderealTime( const double lon ) const
{
    return WrapTwoPI(toGreenwichSiderealTime() + lon);
}

std::string TDateTime::toString() const
{
    std::stringstream ss;
    int year;
    int month;
    int day;
    fromTicks(year, month, day);
    ss << std::right << std::setfill('0');
    ss << std::setw(4) << year << "-";
    ss << std::setw(2) << month << "-";
    ss << std::setw(2) << day << " ";
    ss << std::setw(2) << hour() << ":";
    ss << std::setw(2) << minute() << ":";
    ss << std::setw(2) << second() << ".";
    ss << std::setw(6) << microsecond() << " UTC";
    return ss.str();
}

TMATH_END_NAMESPACE
