#pragma once

struct ImageGcp
{
    double x;
    double y;
};

struct GroundGcp
{
    double X;
    double Y;
    double Z;
};

struct GCP
{
    size_t id;
    ImageGcp image;
    GroundGcp ground;
};

struct ExteriorElements
{
    double Xs;
    double Ys;
    double Zs;
    double varphi;
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