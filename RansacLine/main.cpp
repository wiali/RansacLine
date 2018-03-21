#include "stdafx.h"
#include"fitline.h"
#include"ransac.h"
#include<cstdlib>
#include<cstdio>
#include<ctime>

const int COUNT = 1000;
const int inlierCnt = 400;
const int outlierCnt = COUNT - inlierCnt;
Point2D points[COUNT];

int main()
{
	srand((unsigned int)time(0));

	float a = rand()%100/5.0;
	float b = rand()%100/5.0;
   	printf("Line: (%f, %f)\n", a, b);

   	//inliers
    for(int i = 0 ; i < inlierCnt ; i++)
    {
        points[i].x = i;
        points[i].y = i * a + b + rand()%100/1024.0;
    }
    
    //outliers
    for(int i = inlierCnt; i < COUNT; i++)
    {
		points[i].x = 20 + rand()%100;
		points[i].y = 30 + rand()%100;
    }
 	 
    Line line;
 	//least square fit with only inliers
    FitLine::FitLine2D(points, inlierCnt, NULL, line);

    printf("least square fit: a: %f  b: %f\n", line.a, line.b);
   	
	//Ransac fit with inlier and outlier; note: the result maybe wrong
	Ransac::consensus(points, COUNT, line);
	printf("ransac fit(including outliers): a: %f  b: %f\n", line.a, line.b);
	
	int numForEstimate = 5;
	float successProbability = 0.999f;
	float maxOutliersPercentage = (float)outlierCnt/COUNT; 
	Ransac::consensus(points, COUNT, line, numForEstimate, successProbability, maxOutliersPercentage);
	printf("ransac fit(including outliers): a: %f  b: %f\n", line.a, line.b);
	
    getchar();

	return 0;
}
