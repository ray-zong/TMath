#include "TUtilities.h"

#include <algorithm>
#include <functional>
#include <locale>

#include "TConstants.h"

TMATH_BEGIN_NAMESPACE


double Mod( const double x, const double y )
{
    if (y == 0.0)
    {
        return x;
    }

    return x - y * floor(x / y);
}

double WrapNegPosPI(const double a)
{
    return Mod(a + kPI, kTWOPI) - kPI;
}

double WrapTwoPI(const double a)
{
    return Mod(a, kTWOPI);
}

double WrapNegPos180(const double a)
{
    return Mod(a + 180.0, 360.0) - 180.0;
}

double Wrap360(const double a)
{
    return Mod(a, 360.0);
}

double DegreesToRadians(const double degrees)
{
    return degrees * kPI / 180.0;
}

double RadiansToDegrees(const double radians)
{
    return radians * 180.0 / kPI;
}

double AcTan(const double sinx, const double cosx)
{
    if (cosx == 0.0)
    {
        if (sinx > 0.0)
        {
            return kPI / 2.0;
        }
        else
        {
            return 3.0 * kPI / 2.0;
        }
    }
    else
    {
        if (cosx > 0.0)
        {
            return atan(sinx / cosx);
        }
        else
        {
            return kPI + atan(sinx / cosx);
        }
    }
}

struct IsDigit: std::unary_function<char, bool>
{
    bool operator()(char c) const
    {
        return std::isdigit(c, std::locale::classic()) == 0;
    }
};

void TrimLeft( std::string &s )
{
    s.erase(s.begin(),
        std::find_if(s.begin(), s.end(), std::not1(IsDigit())));
}

void TrimRight( std::string &s )
{
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(IsDigit())).base(),
        s.end());
}

void Trim(std::string &s )
{
    TrimLeft(s);
    TrimRight(s);
}

TMATH_END_NAMESPACE
