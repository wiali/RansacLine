#include "stdafx.h"
#include"fitline.h"
#include <math.h>

const double eps = 1e-6;
//https://www.cnblogs.com/paiandlu/p/7843236.html
int FitLine::FitLine2D(const Point2D* pPoints, int nCount, float* pWeights, Line& line)
{
    float   sumx = 0.0;                      /* sum of x     */
    float   sumx2 = 0.0;                     /* sum of x**2  */
    float   sumxy = 0.0;                     /* sum of x * y */
    float   sumy = 0.0;                      /* sum of y     */
    float   sumy2 = 0.0;                     /* sum of y**2  */

    for (int i = 0; i < nCount; i++)
    {
        sumx += (pWeights ? pWeights[i] : 1) * pPoints[i].x;
        sumx2 += (pWeights ? pWeights[i] : 1) * sqr(pPoints[i].x);
        sumxy += (pWeights ? pWeights[i] : 1) * pPoints[i].x * pPoints[i].y;
        sumy += (pWeights ? pWeights[i] : 1) * pPoints[i].y;
        sumy2 += (pWeights ? pWeights[i] : 1) * sqr(pPoints[i].y);
    }

    float denom = (nCount * sumx2 - sqr(sumx));
    if (denom == 0)
    {
        // singular matrix. can't solve the problem.
        line.a = 0;
        line.b = 0;
        line.r = 0;
        return 1;
    }

    line.a = (nCount * sumxy - sumx * sumy) / denom;;
    line.b = (sumy * sumx2 - sumx * sumxy) / denom;
    /* compute correlation coeff */
    line.r = (sumxy - sumx * sumy / nCount) /
             sqrt((sumx2 - sqr(sumx) / nCount) *
                  (sumy2 - sqr(sumy) / nCount));

    return 0;
}

//http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
double FitLine::CalcDist2D( Point2D* points, int count, Line& line, float* dist)
{
    double sum_dist = 0.;

    for(int j = 0; j < count; j++ )
    {
        float x = points[j].x;
        float y = points[j].y;

        dist[j] = (float) fabs( -line.a * x + y - line.b ) / sqrt(1 + sqr(line.a));
        sum_dist += dist[j];
    }

    return sum_dist;
}

int FitLine::FitLine2D(Point2D* points, int count, Line& line)
{
    //first no weight to least square a initial line.
    FitLine2D(points, count, NULL, line);
    //least square with weight.
    float* dist = new float[count];
    float* weight = new float[count];

    for(int i = 0; i < 5; i++)
    {
        CalcDist2D(points, count, line, dist);
        WeightL1( dist, count, weight);
        FitLine2D(points, count, weight, line);
    }
    delete [] dist;

    return 0;
}

void FitLine::WeightL1( float* d, int count, float* w )
{
    int i;

    for( i = 0; i < count; i++ )
    {
        double t = fabs( (double) d[i] );
        w[i] = (float)(1. / max(t, eps));
    }
}

void FitLine::WeightL12( float* d, int count, float* w )
{
    int i;

    for( i = 0; i < count; i++ )
    {
        w[i] = 1.0f / (float) sqrt( 1 + (double) (d[i] * d[i] * 0.5) );
    }
}


void FitLine::WeightHuber( float* d, int count, float* w, float _c )
{
    int i;
    const float c = _c <= 0 ? 1.345f : _c;

    for( i = 0; i < count; i++ )
    {
        if( d[i] < c )
            w[i] = 1.0f;
        else
            w[i] = c / d[i];
    }
}


void FitLine::WeightFair( float* d, int count, float* w, float _c )
{
    int i;
    const float c = _c == 0 ? 1 / 1.3998f : 1 / _c;

    for( i = 0; i < count; i++ )
    {
        w[i] = 1 / (1 + d[i] * c);
    }
}

void FitLine::WeightWelsch( float* d, int count, float* w, float _c )
{
    int i;
    const float c = _c == 0 ? 1 / 2.9846f : 1 / _c;

    for( i = 0; i < count; i++ )
    {
        w[i] = (float) exp( -d[i] * d[i] * c * c );
    }
}



