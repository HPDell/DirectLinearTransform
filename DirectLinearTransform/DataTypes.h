#pragma once

#include <vector>

using namespace std;
struct ImagePoint
{
    double x;
    double y;
};

typedef ImagePoint ImageGcp;

struct GroundPoint
{
    double X;
    double Y;
    double Z;
};

typedef GroundPoint GroundGcp;

struct GCP
{
    size_t id;
    ImageGcp image;
    GroundGcp ground;
};

struct GCPpp : public GCP {
    double residual;
};

struct ExteriorElements
{
    double Xs;
    double Ys;
    double Zs;
    double phi;
    double omega;
    double kappa;
};

struct InneriorElements
{
    double fx;
    double fy;
    double x0;
    double y0;
    double ds;
    double db;
};

struct SingleTiePoint
{
    size_t imgId;
    ImagePoint ip;
};

struct TiePoint
{
    vector<SingleTiePoint> ties;
    GroundPoint gp;
};


struct ImagePointResidual
{
    size_t id;
    double vx;
    double vy;
    double vxy;
};

struct ResectionReport
{
    InneriorElements inn_rep;
    ExteriorElements ext_rep;
    vector<double> l;
    vector<double> calib_params;
    double m0;
    vector<ImagePointResidual> image_point_residual;
};

struct GroundPointResidual
{
    size_t id;
    double vX;
    double vY;
    double vZ;
};

struct IntersectionReport
{
    double mXY;
    double mZ;
    vector<GroundPointResidual> ground_point_residual;
};