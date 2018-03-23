#include "stdafx.h"
#include"ransac.h"

using namespace std;

//return the percentage of inlines
float Ransac::consensus(Point2D* points, size_t Cnt, Line& line)
{
    int numDataObjects = Cnt;
    int numForEstimate = Cnt * 0.1;
    int maxVoteCnt = 0;
    float inliersPercentage = 0.0;

    int ransac_times = 500;
    int* Chosen = new int[numDataObjects];

    Point2D* subPoints = new Point2D[numForEstimate];
    int pointCnt = 0;
    int voteCnt = 0;
    for(int i = 0; i < ransac_times; i++)
    {
        //randomly select data for exact model fit ('numForEstimate' objects).
        memset(Chosen, 0, numDataObjects * sizeof(int));
        int maxIndex = numDataObjects - 1;
        for(int j = 0; j < numForEstimate; j++)
        {
            int selectedIndex = rand() % numDataObjects;
            Chosen[selectedIndex] = 1;
        }

        pointCnt = 0;
        for(int k = 0; k < numDataObjects; k++)
        {
            if(Chosen[k])
            {
                subPoints[pointCnt].x = points[k].x;
                subPoints[pointCnt].y = points[k].y;
                pointCnt++;
            }
        }

        Line tempLine;
        FitLine::FitLine2D(subPoints, pointCnt, tempLine);

        voteCnt = 0;
        for(int k = 0; k < Cnt; k++)
        {
            if(abs(points[k].y - tempLine.a * points[k].x - tempLine.b) < 2)
            {
                voteCnt++;
            }
        }

        if(voteCnt > maxVoteCnt)
        {
            maxVoteCnt = voteCnt;
            inliersPercentage = (float)maxVoteCnt / Cnt;
            //			printf("a: %f\tb%f\tpercent: %f\n", a, b, inliersPercentage);
            line = tempLine;

        }

        //		if(inliersPercentage > 0.2)
        //		{
        //			return inliersPercentage;
        //		}
    }
    return inliersPercentage;

}

//http://blog.csdn.net/gggg_ggg/article/details/45694901
float Ransac::consensus(Point2D* points, size_t Cnt, Line& line, int numForEstimate,
                        float successProbability, float maxOutliersPercentage)
{

    //1 - p = (1 - w^n)^k
    //p =
    //float outlierPercentage = maxOutliersPercentage;//???
    float numerator = log(1.0f - successProbability);
    float denominator = log(1.0f - pow(1.0 - maxOutliersPercentage, numForEstimate));

    int ransac_times = (int)(numerator / denominator + 0.5);

    printf("ransac_iterator_times: %d\n", ransac_times);
    int numDataObjects = Cnt;
    //int numForEstimate = Cnt*0.1;
    int maxVoteCnt = 0;
    float inliersPercentage = 0.0;

    int* Chosen = new int[numDataObjects];

    Point2D* subPoints = new Point2D[numForEstimate];
    int pointCnt = 0;
    int voteCnt = 0;
    for(int i = 0; i < ransac_times; i++)
    {
        //randomly select data for exact model fit ('numForEstimate' objects).
        memset(Chosen, 0, numDataObjects * sizeof(int));
        int maxIndex = numDataObjects - 1;
        for(int j = 0; j < numForEstimate; j++)
        {
            int selectedIndex = rand() % numDataObjects;
            Chosen[selectedIndex] = 1;
        }
        //fitting
        pointCnt = 0;
        for(int k = 0; k < numDataObjects; k++)
        {
            if(Chosen[k])
            {
                subPoints[pointCnt].x = points[k].x;
                subPoints[pointCnt].y = points[k].y;
                pointCnt++;
            }
        }

        Line tempLine;
        FitLine::FitLine2D(subPoints, pointCnt, tempLine);

        //get the best line
        voteCnt = 0;
        for(int k = 0; k < Cnt; k++)
        {
            //if on the line or near the line
            if(abs(points[k].y - tempLine.a * points[k].x - tempLine.b) < 2)
            {
                voteCnt++;
            }
        }

        if(voteCnt > maxVoteCnt)
        {
            maxVoteCnt = voteCnt;
            inliersPercentage = (float)maxVoteCnt / Cnt;
            //			printf("a: %f\tb%f\tpercent: %f\n", a, b, inliersPercentage);
            line = tempLine;

        }
        //if inliers percentage is high then we can take this as best solve directly.
        //		if(inliersPercentage > 0.2)
        //		{
        //			return inliersPercentage;
        //		}
    }
    return inliersPercentage;
}




