#include "stdafx.h"
#include"fitline.h"
#include"ransac.h"
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include <math.h>
#include <vector>

#include "KdTree.h"
#include "PointCloud.h"
#include "levmarq_test.h"

const int COUNT = 1000;
const int inlierCnt = 400;
const int outlierCnt = COUNT - inlierCnt;
Point2D points[COUNT];

using namespace std;

//-----------------------------------------------------------------------------------------------------------
void initialData()
{
    srand((unsigned int)time(0));

    float a = rand() % 100 / 5.0;
    float b = rand() % 100 / 5.0;
    printf("Line: (%f, %f)\n", a, b);

    //inliers
    for (int i = 0; i < inlierCnt; i++)
    {
        points[i].x = i;
        points[i].y = i * a + b + rand() % 100 / 1024.0;
    }

    //outliers
    for (int i = inlierCnt; i < COUNT; i++)
    {
        points[i].x = 20 + rand() % 100;
        points[i].y = 30 + rand() % 100;
    }
}

//-----------------------------------------------------------------------------------------------------------
void runFPFH()
{
    double rN = 50;
    double rH = 25;
    double binsize = 0.1;
    #define M_PI  3.14159265358979323846   // pi
    double minsize = -M_PI;
    double maxsize = M_PI;

    PointCloud cloud;

    for (int i = 0; i < COUNT; ++i)
    {
        float z = 1;
        if (i % 100 == 0)
            z = rand() % 100;

        Point* pts = new Point(i, points[i].x, points[i].y, z);
        cloud.addPoint(pts);
    }

    //build kd-tree from cloud
    KdTreeNode root;
    root.buildTree(cloud.getPoints(), 0);

    //compute normals for all points with radius rN
    int count = 0;
    std::cout << "computing normals..." << std::endl;
    for (unsigned int i = 0; i < cloud.getSize(); i++) {
        count++;
        if (count % 1000 == 0) {
            std::cout << count << "..." << std::endl;
        }
        std::vector<Point*> neighbors;
        root.fixedRadiusSearch(cloud.getPoints()[i], rN, neighbors);
        cloud.getPoints()[i]->computeNormal(neighbors);
    }

    //compute SPFHs for all points with radius rH
    count = 0;
    std::cout << "computing SPFHs..." << std::endl;
    for (unsigned int i = 0; i < cloud.getSize(); i++) {
        count++;
        if (count % 1000 == 0) {
            std::cout << count << "..." << std::endl;
        }
        std::vector<Point*> neighbors;
        root.fixedRadiusSearch(cloud.getPoints()[i], rH, neighbors);

        SPFH* s = new SPFH(minsize, maxsize, binsize);
        for (unsigned int j = 0; j < neighbors.size(); j++) {
            if (cloud.getPoints()[i] == neighbors[j]) {
                continue;
            }
            Features* f = new Features(cloud.getPoints()[i], neighbors[j]);
            s->addToHistogram(f);
            delete f;
        }
        cloud.getPoints()[i]->setSimplePointFeatureHistogram(s);
    }

    //compute SPFHs for all points with radius rH
    count = 0;
    std::cout << "computing FPFHs..." << std::endl;
    for (unsigned int i = 0; i < cloud.getSize(); i++) {
        count++;
        if (count % 1000 == 0) {
            std::cout << count << "..." << std::endl;
        }
        std::vector<Point*> neighbors;
        root.fixedRadiusSearch(cloud.getPoints()[i], rH, neighbors);
        FPFH* f = new FPFH(minsize, maxsize, binsize);
        f->addToHistogram(
            cloud.getPoints()[i]->getSimplePointFeatureHistogram(), 1, 1);

        for (unsigned int j = 0; j < neighbors.size(); j++) {
            if (cloud.getPoints()[i] == neighbors[j]) {
                continue;
            }
            f->addToHistogram(neighbors[j]->getSimplePointFeatureHistogram(),
                neighbors.size(),
                cloud.getPoints()[i]->euclideanDistance(neighbors[j]));
        }
        cloud.getPoints()[i]->setFastPointFeatureHistogram(f);

    }
}

//-----------------------------------------------------------------------------------------------------------
void runRansac()
{
    Line line;
    //least square fit with only inliers
    FitLine::FitLine2D(points, inlierCnt, NULL, line);

    printf("least square fit: a: %f  b: %f\n", line.a, line.b);

    //Ransac fit with inlier and outlier; note: the result maybe wrong
    Ransac::consensus(points, COUNT, line);
    printf("ransac fit(including outliers): a: %f  b: %f\n", line.a, line.b);

    int numForEstimate = 5;
    float successProbability = 0.999f;
    float maxOutliersPercentage = (float)outlierCnt / COUNT;
    Ransac::consensus(points, COUNT, line, numForEstimate, successProbability, maxOutliersPercentage);
    printf("ransac fit(including outliers): a: %f  b: %f\n", line.a, line.b);

}

//-----------------------------------------------------------------------------------------------------------
int main()
{
    LevMarq_Test::run();

    initialData();

    runFPFH();

    runRansac();

    getchar();

	return 0;
}
