#ifndef TTIMESPAN_H_
#define TTIMESPAN_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdint>

#include "../TGlobal.h"

TMATH_BEGIN_NAMESPACE

static const int64_t TicksPerDay =  86400000000LL;
static const int64_t TicksPerHour =  3600000000LL;
static const int64_t TicksPerMinute =  60000000LL;
static const int64_t TicksPerSecond =   1000000LL;
static const int64_t TicksPerMillisecond = 1000LL;
static const int64_t TicksPerMicrosecond =    1LL;

static const int64_t UnixEpoch = 62135596800000000LL;

static const int64_t MaxValueTicks = 315537897599999999LL;

// 1582-Oct-15
static const int64_t GregorianStart = 49916304000000000LL;


/**
* @brief Represents a time interval.
*
* Represents a time interval (duration/elapsed) that is measured as a positive
* or negative number of days, hours, minutes, seconds, and fractions
* of a second.
*/
class TMATH_EXPORT TTimeSpan
{
public:
    TTimeSpan(int64_t ticks);

    TTimeSpan(int hours, int minutes, int seconds);

    TTimeSpan(int days, int hours, int minutes, int seconds);

    TTimeSpan(int days, int hours, int minutes, int seconds, int microseconds);

    TTimeSpan add(const TTimeSpan& ts) const;

    TTimeSpan subtract(const TTimeSpan& ts) const;

    int compare(const TTimeSpan& ts) const;

    bool equals(const TTimeSpan& ts) const;

    int days() const;

    int hours() const;

    int minutes() const;

    int seconds() const;

    int milliseconds() const;

    int microseconds() const;

    int64_t ticks() const;

    double totalDays() const;

    double totalHours() const;

    double totalMinutes() const;

    double totalSeconds() const;

    double totalMilliseconds() const;

    double totalMicroseconds() const;

    std::string toString() const;

private:
    int64_t m_ticks;

    void calculateTicks(int days,
        int hours,
        int minutes,
        int seconds,
        int microseconds);
};

inline std::ostream& operator<<(std::ostream& strm, const TTimeSpan& t)
{
    return strm << t.toString();
}

inline TTimeSpan operator+(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return ts1.add(ts2);
}

inline TTimeSpan operator-(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return ts1.subtract(ts2);
}

inline bool operator==(const TTimeSpan& ts1, TTimeSpan& ts2)
{
    return ts1.equals(ts2);
}

inline bool operator>(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return (ts1.compare(ts2) > 0);
}

inline bool operator>=(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return (ts1.compare(ts2) >= 0);
}

inline bool operator!=(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return !ts1.equals(ts2);
}

inline bool operator<(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return (ts1.compare(ts2) < 0);
}

inline bool operator<=(const TTimeSpan& ts1, const TTimeSpan& ts2)
{
    return (ts1.compare(ts2) <= 0);
}

TMATH_END_NAMESPACE

#endif
