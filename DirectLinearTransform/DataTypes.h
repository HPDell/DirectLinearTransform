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