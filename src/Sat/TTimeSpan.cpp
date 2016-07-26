#include "TTimeSpan.h"

TMATH_BEGIN_NAMESPACE

TTimeSpan::TTimeSpan( int64_t ticks )
    : m_ticks(ticks)
{

}

TTimeSpan::TTimeSpan( int hours, int minutes, int seconds )
{
    calculateTicks(0, hours, minutes, seconds, 0);
}

TTimeSpan::TTimeSpan( int days, int hours, int minutes, int seconds )
{
    calculateTicks(days, hours, minutes, seconds, 0);
}

TTimeSpan::TTimeSpan( int days, int hours, int minutes, int seconds, int microseconds )
{
    calculateTicks(days, hours, minutes, seconds, microseconds);
}

TTimeSpan TTimeSpan::add( const TTimeSpan& ts ) const
{
    return TTimeSpan(m_ticks + ts.m_ticks);
}

TTimeSpan TTimeSpan::subtract( const TTimeSpan& ts ) const
{
    return TTimeSpan(m_ticks - ts.m_ticks);
}

int TTimeSpan::compare( const TTimeSpan& ts ) const
{
    int ret = 0;

    if (m_ticks < ts.m_ticks)
    {
        ret = -1;
    }
    if (m_ticks > ts.m_ticks)
    {
        ret = 1;
    }
    return ret;
}

bool TTimeSpan::equals( const TTimeSpan& ts ) const
{
    return m_ticks == ts.m_ticks;
}

int TTimeSpan::days() const
{
    return static_cast<int>(m_ticks / TicksPerDay);
}

int TTimeSpan::hours() const
{
    return static_cast<int>(m_ticks % TicksPerDay / TicksPerHour);
}

int TTimeSpan::minutes() const
{
    return static_cast<int>(m_ticks % TicksPerHour / TicksPerMinute);
}

int TTimeSpan::seconds() const
{
    return static_cast<int>(m_ticks % TicksPerMinute / TicksPerSecond);
}

int TTimeSpan::milliseconds() const
{
    return static_cast<int>(m_ticks % TicksPerSecond / TicksPerMillisecond);
}

int TTimeSpan::microseconds() const
{
    return static_cast<int>(m_ticks % TicksPerSecond / TicksPerMicrosecond);
}

std::int64_t TTimeSpan::ticks() const
{
    return m_ticks;
}

double TTimeSpan::totalDays() const
{
    return static_cast<double>(m_ticks) / TicksPerDay;
}

double TTimeSpan::totalHours() const
{
    return static_cast<double>(m_ticks) / TicksPerHour;
}

double TTimeSpan::totalMinutes() const
{
    return static_cast<double>(m_ticks) / TicksPerMinute;
}

double TTimeSpan::totalSeconds() const
{
    return static_cast<double>(m_ticks) / TicksPerSecond;
}

double TTimeSpan::totalMilliseconds() const
{
    return static_cast<double>(m_ticks) / TicksPerMillisecond;
}

double TTimeSpan::totalMicroseconds() const
{
    return static_cast<double>(m_ticks) / TicksPerMicrosecond;
}

std::string TTimeSpan::toString() const
{
    std::stringstream ss;

    ss << std::right << std::setfill('0');

    if (m_ticks < 0)
    {
        ss << '-';
    }

    if (days() != 0)
    {
        ss << std::setw(2) << std::abs(days()) << '.';
    }

    ss << std::setw(2) << std::abs(hours()) << ':';
    ss << std::setw(2) << std::abs(minutes()) << ':';
    ss << std::setw(2) << std::abs(seconds());

    if (microseconds() != 0)
    {
        ss << '.' << std::setw(6) << std::abs(microseconds());
    }

    return ss.str();
}

void TTimeSpan::calculateTicks( int days, int hours, int minutes, int seconds, int microseconds )
{
    m_ticks = days * TicksPerDay +
        (hours * 3600LL + minutes * 60LL + seconds) * TicksPerSecond + 
        microseconds * TicksPerMicrosecond;
}


TMATH_END_NAMESPACE
