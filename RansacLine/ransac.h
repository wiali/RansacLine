#ifndef _RANSAC_H_
#define _RANSAC_H_

#include<cstdlib>
#include"fitline.h"
#include<iostream>
#include<ctime>
#include<cstring>
#include<cstdio>

class Ransac
{
public:
    static float consensus(Point2D* points, size_t Cnt, Line& line);
    static float consensus(Point2D* points, size_t Cnt, Line& line, int numForEstimate, float successProbability, float maxOutliersPercentage);
};

#endif
