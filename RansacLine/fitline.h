#ifndef FITLINE_H
#define FITLINE_H

#include"type.h"
#include<cmath>

class FitLine
{
public:
    static int FitLine2D(const Point2D* pPoints, int nCount, float* pWeights, Line& line);
    static double CalcDist2D(Point2D* points, int count, Line& line, float* dist);
    static int FitLine2D(Point2D* points, int count, Line& line);

private:
    inline static float sqr(float x)
    {
        return x*x;
    }

    inline static double max(double a, double b)
    {
        return a > b ? a : b;
    }

    static void WeightL1(float *d, int count, float *w);
    static void WeightL12(float *d, int count, float *w);
    static void WeightHuber(float *d, int count, float *w, float _c);
    static void WeightFair(float *d, int count, float *w, float _c);
    static void WeightWelsch(float *d, int count, float *w, float _c);
};

#endif
