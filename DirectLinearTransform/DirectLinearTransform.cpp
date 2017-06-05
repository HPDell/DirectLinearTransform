#include "stdafx.h"
#include "DirectLinearTransform.h"
#include <Eigen>
#include <time.h>

#ifdef _DEBUG
#include <iostream>
#include <iomanip>
#endif // _DEBUG


#define random(x) (rand() % x)

using namespace Eigen;

CDirectLinearTransform::CDirectLinearTransform()
    : m_nCalibNums(0)
{
    ext.Xs = 0.0; ext.Ys = 0.0; ext.Zs = 0.0;
    ext.phi = 0.0; ext.omega = 0.0; ext.kappa = 0.0;
    inn.x0 = 0.0; inn.y0 = 0.0; 
    inn.fx = 0.0; inn.fy = 0.0;
    inn.ds = 0.0; inn.db = 0.0;
    k = nullptr;
    memset(l, 0, 11 * sizeof(double));
}

CDirectLinearTransform::CDirectLinearTransform(size_t nCalibNums)
    : m_nCalibNums(nCalibNums)
{
    ext.Xs = 0.0; ext.Ys = 0.0; ext.Zs = 0.0;
    ext.phi = 0.0; ext.omega = 0.0; ext.kappa = 0.0;
    k = new double[m_nCalibNums];
    memset(k, 0, m_nCalibNums * sizeof(double));
    memset(l, 0, 11 * sizeof(double));
}


CDirectLinearTransform::~CDirectLinearTransform()
{
    vGcps.clear();
    delete[] k;
    k = nullptr;
}

void CDirectLinearTransform::Resection()
{
    // 入口处检查
    if (vGcps.size() < 6)
    {
        return;
    }

    // l系数近似值解算
    Resection_l_Approx();
    // 内方位元素近似值解算
    CalcInneriorElements();
    // 设置迭代条件
    double fx0 = -DBL_MAX;
    double fx1 = inn.fx;
    while (abs(fx1 - fx0) > 0.01)
    {
        fx0 = fx1;
        // 求解l系数精确值
        Resection_l_Exact();
        // 内方位元素
        CalcInneriorElements();
        fx1 = inn.fx;
    }
    CalcExteriorElements();
}

void CDirectLinearTransform::Intersection(CDirectLinearTransform& right)
{
    // 修正像点坐标
    this->RectifyImagePoint();
    right.RectifyImagePoint();
    // 前方交会
    size_t nImagePointNum = vTiePoints.size();
    for (size_t i = 0; i < nImagePointNum; i++)
    {
        // 计算近似值
        ImagePoint& lp = vTiePoints.at(i).image;
        ImagePoint& rp = right.vTiePoints.at(i).image;
        GroundPoint& gp = vTiePoints.at(i).ground;
        double* ll = l;
        double* rl = right.l;
        Intersection_SpaceApprox(lp, rp, ll, rl, gp);
        // 精确值
        Intersection_SpaceExact(lp, rp, ll, rl, gp);
    }

}

void CDirectLinearTransform::Intersection(CDirectLinearTransform * dlt, size_t dlt_num)
{
}

double CDirectLinearTransform::ResectionPrecision(double* pointResidual)
{
    double m0 = 0.0;
    size_t nPointNum = vGcps.size();
    for (size_t i = 0; i < nPointNum; i++)
    {
        GCP& iGcp = vGcps.at(i);
        double err_x = 0.0, err_y = 0.0;
        ImagePoint imgPointCalc = CalcImageResidual(iGcp);
        err_x = imgPointCalc.x - iGcp.image.x/* / (l[8] * iGcp.ground.X + l[9] * iGcp.ground.Y + l[10] * iGcp.ground.Z + 1)*/;
        err_y = imgPointCalc.y - iGcp.image.y/* / (l[8] * iGcp.ground.X + l[9] * iGcp.ground.Y + l[10] * iGcp.ground.Z + 1)*/;
        double residual = err_x*err_x + err_y*err_y;
        m0 += residual;
        if (pointResidual != nullptr) pointResidual[i] = err_x;
    }
    m0 = sqrt(m0 / (2 * nPointNum - 11));
    return m0;
}

void CDirectLinearTransform::CalcInneriorElements()
{
    double gamma3t = l[8] * l[8] + l[9] * l[9] + l[10] * l[10];
    // 像主点位置
    inn.x0 = -(l[0] * l[8] + l[1] * l[9] + l[2] * l[10]) / gamma3t;
    inn.y0 = -(l[4] * l[8] + l[5] * l[9] + l[6] * l[10]) / gamma3t;
    // A B C
    double A = (l[0] * l[0] + l[1] * l[1] + l[2] * l[2]) / gamma3t - inn.x0 * inn.x0;
    double B = (l[4] * l[4] + l[5] * l[5] + l[6] * l[6]) / gamma3t - inn.y0 * inn.y0;
    double C = (l[0] * l[4] + l[1] * l[5] + l[2] * l[6]) / gamma3t - inn.x0 * inn.y0;
    inn.fx = sqrt((A*B - C*C) / B);
    inn.fy = sqrt((A*B - C*C) / A);
    inn.ds = sqrt(A / B) - 1;
    inn.db = sqrt((C*C) / (A*B));
    inn.db = (C < 0) ? asin(inn.db) : -asin(inn.db);
}

void CDirectLinearTransform::CalcExteriorElements()
{
    // 外放为线元素
    Matrix3d matCoeff;
    matCoeff << l[0], l[1], l[2], l[4], l[5], l[6], l[8], l[9], l[10];
    Vector3d matConst;
    matConst << l[3], l[7], -1;
    Vector3d matExtLine = matCoeff.lu().solve(matConst);
    ext.Xs = matExtLine[0];
    ext.Ys = matExtLine[1];
    ext.Zs = matExtLine[2];
    // 外方位角元素
    ext.phi = atan(l[8] / l[10]);
    //double b1, b2, b3, a3, c3;
    double b3;
    //a3 = l[8] / sqrt(l[8] * l[8] + l[9] * l[9] + l[10] * l[10]);
    b3 = l[9] / sqrt(l[8] * l[8] + l[9] * l[9] + l[10] * l[10]);
    //c3 = l[10] / sqrt(l[8] * l[8] + l[9] * l[9] + l[10] * l[10]);
    ext.omega = asin(-b3);
    //ext.kappa = 
}

void CDirectLinearTransform::Resection_l_Approx()
{
    srand((unsigned)time(NULL));
    size_t nGcpNum = vGcps.size();
    GCP* pGcp = vGcps.data();
    // 随机抽取6个点
    // 构建矩阵
    MatrixXd matCoeff(12, 11);
    VectorXd vecConst(12);
    for (size_t i = 0; i < 6; i++)
    {
        size_t iPoint = random(nGcpNum);
        GCP& iGcp = *(pGcp + i);
        // 系数矩阵
        matCoeff(i*2, 0) = iGcp.ground.X;  matCoeff(i*2 + 1, 0) = 0;
        matCoeff(i*2, 1) = iGcp.ground.Y;  matCoeff(i*2 + 1, 1) = 0;
        matCoeff(i*2, 2) = iGcp.ground.Z;  matCoeff(i*2 + 1, 2) = 0;
        matCoeff(i*2, 3) = 1;  matCoeff(i*2 + 1, 3) = 0;
        matCoeff(i*2, 4) = 0;  matCoeff(i*2 + 1, 4) = iGcp.ground.X;
        matCoeff(i*2, 5) = 0;  matCoeff(i*2 + 1, 5) = iGcp.ground.Y;
        matCoeff(i*2, 6) = 0;  matCoeff(i*2 + 1, 6) = iGcp.ground.Z;
        matCoeff(i*2, 7) = 0;  matCoeff(i*2 + 1, 7) = 1;
        matCoeff(i*2, 8) = iGcp.image.x * iGcp.ground.X;  matCoeff(i*2 + 1, 8) = iGcp.image.y * iGcp.ground.X;
        matCoeff(i*2, 9) = iGcp.image.x * iGcp.ground.Y;  matCoeff(i*2 + 1, 9) = iGcp.image.y * iGcp.ground.Y;
        matCoeff(i*2, 10) = iGcp.image.x * iGcp.ground.Z; matCoeff(i*2 + 1, 10) = iGcp.image.y * iGcp.ground.Z;
        // 常数项矩阵
        vecConst(i*2) = iGcp.image.x;
        vecConst(i*2 + 1) = iGcp.image.y;
    }
    //cout << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(6) << setw(8) << matCoeff;
    VectorXd vecL = (matCoeff.transpose() * matCoeff).lu().solve(matCoeff.transpose() * vecConst);
    for (size_t i = 0; i < 11; i++) l[i] = vecL[i];
}

void CDirectLinearTransform::Resection_l_Exact()
{
    size_t nGcpNums = vGcps.size();
    GCP* pGcp = vGcps.data();

    // 逐点法化
    MatrixXd MTM = MatrixXd::Zero(11 + m_nCalibNums, 11 + m_nCalibNums);
    VectorXd MTW = VectorXd::Zero(11 + m_nCalibNums);
    MatrixXd iM(2, 11 + m_nCalibNums);
    VectorXd iW(2);
    for (size_t i = 0; i < nGcpNums; i++)
    {
        GCP& iGcp = *(pGcp + i);
        double A = l[8] * iGcp.ground.X + l[9] * iGcp.ground.Y + l[10] * iGcp.ground.Z;
        double Dx = iGcp.image.x - inn.x0, Dy = iGcp.image.y - inn.y0;
        double r2 = Dx*Dx + Dy*Dy;
        // 系数矩阵
        iM(0, 0) = -iGcp.ground.X / A; iM(1, 0) = 0;
        iM(0, 1) = -iGcp.ground.Y / A; iM(1, 1) = 0;
        iM(0, 2) = -iGcp.ground.Z / A; iM(1, 2) = 0;
        iM(0, 3) = -1 / A;             iM(1, 3) = 0;
        iM(0, 4) = 0; iM(1, 4) = -iGcp.ground.X / A;
        iM(0, 5) = 0; iM(1, 5) = -iGcp.ground.Y / A;
        iM(0, 6) = 0; iM(1, 6) = -iGcp.ground.Z / A;
        iM(0, 7) = 0; iM(1, 7) = -1 / A;
        iM(0, 8) = -iGcp.image.x * iGcp.ground.X / A;  iM(1, 8) = -iGcp.image.y * iGcp.ground.X / A;
        iM(0, 9) = -iGcp.image.x * iGcp.ground.Y / A;  iM(1, 9) = -iGcp.image.y * iGcp.ground.Y / A;
        iM(0, 10) = -iGcp.image.x * iGcp.ground.Z / A; iM(1, 10) = -iGcp.image.y * iGcp.ground.Z / A;
        for (size_t ik = 0; ik < m_nCalibNums; ik++)
        {
            iM(0, 11 + ik) = -Dx * pow(r2, ik + 1);
            iM(1, 11 + ik) = -Dy * pow(r2, ik + 1);
        }
        // 常数项向量
        iW(0) = iGcp.image.x / A;
        iW(1) = iGcp.image.y / A;
        // 法化
        MTM += iM.transpose() * iM;
        MTW += iM.transpose() * iW;
#ifdef SHOW_MATRIXS
        // 输出
        cout << i << "M:" << endl;
        cout << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(3) << iM << endl;
        cout << i << "W:" << endl;
        cout << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(3) << iW << endl;
        cout << i << "MTM:" << endl;
        cout << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(3) << MTM << endl;
        cout << i << "MTW:" << endl;
        cout << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(3) << MTW << endl;
#endif // _DEBUG
    }
    VectorXd vecL = MTM.lu().solve(MTW);
#ifdef SHOW_MATRIXS
    cout << MTM.determinant() << endl;
    cout << "MTM" << endl << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(3) << MTM << endl;
    cout << "MTW" << endl << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(3) << MTW << endl;
    cout << "L" << endl << setiosflags(ios::right) << setiosflags(ios::fixed) << setprecision(6) << vecL << endl;
#endif // _DEBUG
    for (size_t i = 0; i < 11; i++) l[i] = vecL(i);
    for (size_t i = 0; i < m_nCalibNums; i++) k[i] = vecL(11 + i);
}

ImagePoint CDirectLinearTransform::CalcImageResidual(GCP & gcp)
{
    GroundPoint& gp = gcp.ground;
    ImagePoint& ip = gcp.image;
    double x = (l[0] * gp.X + l[1] * gp.Y + l[2] * gp.Z + l[3]) / -(l[8] * gp.X + l[9] * gp.Y + l[10] * gp.Z + 1);
    double y = (l[4] * gp.X + l[5] * gp.Y + l[6] * gp.Z + l[7]) / -(l[8] * gp.X + l[9] * gp.Y + l[10] * gp.Z + 1);
    ImagePoint imagePoint = { x, y };
    //ImagePoint& ip = gcp.image;
    double r2 = (ip.x - inn.x0)*(ip.x - inn.x0) + (ip.y - inn.y0)*(ip.y - inn.y0);
    for (size_t i = 0; i < m_nCalibNums; i++)
    {
        imagePoint.x -= k[i] * (ip.x - inn.x0) * pow(r2, i + 1);
        imagePoint.y -= k[i] * (ip.y - inn.y0) * pow(r2, i + 1);
    }
    return imagePoint;
}

void CDirectLinearTransform::RectifyImagePoint()
{
    size_t nImagePointNum = vTiePoints.size();
    for (size_t i = 0; i < nImagePointNum; i++)
    {
        ImagePoint& iIp = vTiePoints.at(i).image;
        double dx = 0.0, dy = 0.0;
        double Dx = iIp.x - inn.x0, Dy = iIp.y - inn.y0;
        double r2 = (Dx * Dx) + (Dy * Dy);
        for (size_t ik = 0; ik < m_nCalibNums; ik++)
        {
            dx += Dx * pow(r2, ik + 1) * k[ik];
            dy += Dy * pow(r2, ik + 1) * k[ik];
        }
        iIp.x += dx;
        iIp.y += dy;
    }
}

void CDirectLinearTransform::Intersection_SpaceApprox(ImagePoint& lp, ImagePoint& rp, double* ll, double* rl, GroundPoint& gp)
{
    MatrixXd matCoeff(4, 3);
    VectorXd vecConst(4);
    matCoeff(0, 0) = ll[0] + lp.x * ll[8]; matCoeff(0, 1) = ll[1] + lp.x * ll[9]; matCoeff(0, 2) = ll[2] + lp.x * ll[10];
    matCoeff(1, 0) = ll[4] + lp.x * ll[8]; matCoeff(1, 1) = ll[5] + lp.x * ll[9]; matCoeff(1, 2) = ll[6] + lp.x * ll[10];
    matCoeff(2, 0) = rl[0] + rp.x * rl[8]; matCoeff(2, 1) = rl[1] + rp.x * rl[9]; matCoeff(2, 2) = rl[2] + rp.x * rl[10];
    matCoeff(3, 0) = rl[4] + rp.x * rl[8]; matCoeff(3, 1) = rl[5] + rp.x * rl[9]; matCoeff(3, 2) = rl[6] + rp.x * rl[10];
    vecConst(0) = -ll[3] * lp.x; vecConst(1) = -ll[7] * lp.y;
    vecConst(2) = -rl[3] * rp.x; vecConst(3) = -rl[7] * rp.y;
    VectorXd matX = (matCoeff.transpose() * matCoeff).lu().solve(matCoeff.transpose() * vecConst);
    gp.X = matX[0]; gp.Y = matX[1]; gp.Z = matX[2];
}

void CDirectLinearTransform::Intersection_SpaceExact(ImagePoint& lp, ImagePoint& rp, double* ll, double* rl, GroundPoint& gp)
{
    //double X = gp.X, X0 = -DBL_MAX;
    double delta = DBL_MAX;
    MatrixXd matCoeff(4, 3);
    Vector4d vecConst(4);
    VectorXd matX;
    while (abs(delta) > 1.0e-4)
    {
        double X0 = gp.X, Y0 = gp.Y, Z0 = gp.Z;
        double Al = ll[8] * gp.Y + ll[9] * gp.Y + ll[10] * gp.Z + 1;
        double Ar = rl[8] * gp.Y + rl[9] * gp.Y + rl[10] * gp.Z + 1;
        matCoeff(0, 0) = -(ll[0] + lp.x * ll[8]) / Al; matCoeff(0, 1) = -(ll[1] + lp.x * ll[9]) / Al; matCoeff(0, 2) = -(ll[2] + lp.x * ll[10]) / Al;
        matCoeff(1, 0) = -(ll[4] + lp.x * ll[8]) / Al; matCoeff(1, 1) = -(ll[5] + lp.x * ll[9]) / Al; matCoeff(1, 2) = -(ll[6] + lp.x * ll[10]) / Al;
        matCoeff(2, 0) = -(rl[0] + rp.x * rl[8]) / Ar; matCoeff(2, 1) = -(rl[1] + rp.x * rl[9]) / Ar; matCoeff(2, 2) = -(rl[2] + rp.x * rl[10]) / Ar;
        matCoeff(3, 0) = -(rl[4] + rp.x * rl[8]) / Ar; matCoeff(3, 1) = -(rl[5] + rp.x * rl[9]) / Ar; matCoeff(3, 2) = -(rl[6] + rp.x * rl[10]) / Ar;
        vecConst(0) = (-ll[3] * lp.x) / Al; vecConst(1) = (-ll[7] * lp.y) / Al;
        vecConst(2) = (-rl[3] * rp.x) / Ar; vecConst(3) = (-rl[7] * rp.y) / Ar;
        matX = (matCoeff.transpose() * matCoeff).lu().solve(matCoeff.transpose() * vecConst);
        gp.X = matX[0]; gp.Y = matX[1]; gp.Z = matX[2];
        delta = abs(max(gp.X - X0, max(gp.Y - Y0, gp.Z - Z0)));
    }
}

