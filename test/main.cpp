#include <iostream>

#include "TVector.h"
#include "TMatrix.h"
#include "TTle.h"
#include "TSGP4.h"
#include "TECI.h"
#include "TSatelliteException.h"
#include <QDebug>
#include <QString>
#include <fstream>
#include "TDateTime.h"
#include "TKepler.h"

using TMath::TVector;
using TMath::TMatrix;
using TMath::TTle;
using TMath::TSGP4;
using TMath::TECI;
using TMath::TSatelliteException;

using namespace TMath;

int main()
{

    TVector<double, 6> elements;
    //六根数 半长轴,离心率,轨道倾角,升交点赤经,近地点幅角,平近点角
    elements[0] = 6.67814e6;
    elements[1] = 1.44389e-15;
    elements[2] = 0.497419;
    elements[3] = 1.93691e-17;
    elements[4] = 0.0982197;
    elements[5] = 6.25438;

    double time = 0;
    TVector<double, 3> rVector;
    TVector<double, 3> vVector;
    TKepler::elemToCoor(elements, time, rVector, vVector);

    rVector = TVector<double, 3>(6662058.85140649,407041.030561599,221005.247530123);
    vVector = TVector<double, 3>(-535.827862619766,6773.17939681257,3677.53635817301);
    TKepler::coorToElem(rVector, vVector, elements);

    TKepler::elemToCoor(elements, time, rVector, vVector);

    return 0;
}
