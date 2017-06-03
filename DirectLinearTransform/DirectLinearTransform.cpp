#include "stdafx.h"
#include "DirectLinearTransform.h"
#include <Eigen>

#define random(x) (rand() % x)

using namespace Eigen;

CDirectLinearTransform::CDirectLinearTransform()
    : m_nCalibNums(0)
{
    ext.Xs = 0.0; ext.Ys = 0.0; ext.Zs = 0.0;
    ext.varphi = 0.0; ext.omega = 0.0; ext.kappa = 0.0;
    k = nullptr;
    memset(l, 0, 11 * sizeof(double));
}

CDirectLinearTransform::CDirectLinearTransform(size_t nCalibNums)
    : m_nCalibNums(nCalibNums)
{
    ext.Xs = 0.0; ext.Ys = 0.0; ext.Zs = 0.0;
    ext.varphi = 0.0; ext.omega = 0.0; ext.kappa = 0.0;
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
    // ��ڴ����
    if (vGcps.size() < 6)
    {
        return;
    }

    // lϵ������ֵ����
    Resection_l_Approx();
    // �ڷ�λԪ�ؽ���ֵ����
    CalcInneriorElements();
    // ���õ�������
    double fx0 = -DBL_MAX;
    double fx1 = inn.fx;
    while (abs(fx1 - fx0) < 0.01)
    {
        fx0 = fx1;
        // ���lϵ����ȷֵ
        Resection_l_Precise();
        // �ڷ�λԪ��
        CalcInneriorElements();
        fx1 = inn.fx;
    }
}

void CDirectLinearTransform::CalcInneriorElements()
{
    double gamma3t = l[8] * l[8] + l[9] * l[9] + l[10] * l[10];
    // ������λ��
    inn.x0 = -(l[0] * l[8] + l[1] * l[9] + l[2] * l[10]) * gamma3t;
    inn.y0 = -(l[4] * l[8] + l[5] * l[9] + l[6] * l[10]) * gamma3t;
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

void CDirectLinearTransform::Resection_l_Approx()
{
    size_t nGcpNum = vGcps.size();
    GCP* pGcp = vGcps.data();
    // �����ȡ6����
    // ��������
    Matrix<double, 12, 11> matCoeff;
    VectorXd vecConst(12);
    for (size_t i = 0; i < 6; i++)
    {
        size_t iPoint = random(nGcpNum);
        GCP& iGcp = *(pGcp + i);
        // ϵ������
        matCoeff(i, 0) = iGcp.ground.X;  matCoeff(i + 1, 0) = 0;
        matCoeff(i, 1) = iGcp.ground.Y;  matCoeff(i + 1, 1) = 0;
        matCoeff(i, 2) = iGcp.ground.Z;  matCoeff(i + 1, 2) = 0;
        matCoeff(i, 3) = 1;  matCoeff(i + 1, 3) = 0;
        matCoeff(i, 4) = 0;  matCoeff(i + 1, 4) = iGcp.ground.X;
        matCoeff(i, 5) = 0;  matCoeff(i + 1, 5) = iGcp.ground.Y;
        matCoeff(i, 6) = 0;  matCoeff(i + 1, 6) = iGcp.ground.Z;
        matCoeff(i, 7) = 0;  matCoeff(i + 1, 7) = 1;
        matCoeff(i, 8) = iGcp.image.x * iGcp.ground.X;  matCoeff(i + 1, 8) = iGcp.image.y * iGcp.ground.X;
        matCoeff(i, 9) = iGcp.image.x * iGcp.ground.Y;  matCoeff(i + 1, 9) = iGcp.image.y * iGcp.ground.Y;
        matCoeff(i, 10) = iGcp.image.x * iGcp.ground.Z; matCoeff(i + 1, 10) = iGcp.image.y * iGcp.ground.Z;
        // ���������
        vecConst(i) = iGcp.image.x;
        vecConst(i + 1) = iGcp.image.y;
    }
    VectorXd vecL(11);
    vecL = matCoeff.lu().solve(vecConst);
    for (size_t i = 0; i < 11; i++)
    {
        l[i] = vecL[i];
    }
}

void CDirectLinearTransform::Resection_l_Precise()
{
    size_t nGcpNums = vGcps.size();
    GCP* pGcp = vGcps.data();

    // ��㷨��
    MatrixXd MTM(11 + m_nCalibNums, 11 + m_nCalibNums);
    VectorXd MTW(11 + m_nCalibNums);
    MatrixXd iM(2, 11 + m_nCalibNums);
    VectorXd iW(2);
    for (size_t i = 0; i < nGcpNums; i++)
    {
        GCP& iGcp = *(pGcp + i);
        double A = l[8] * iGcp.ground.X + l[9] * iGcp.ground.Y + l[10] * iGcp.ground.Z;
        double Dx = iGcp.image.x - inn.x0, Dy = iGcp.image.y - inn.y0;
        double r2 = Dx*Dx + Dy*Dy;
        // ϵ������
        iM(0, 0) = iGcp.ground.X / A; iM(1, 0) = 0;
        iM(0, 1) = iGcp.ground.Y / A; iM(1, 1) = 0;
        iM(0, 2) = iGcp.ground.Z / A; iM(1, 2) = 0;
        iM(0, 3) = 1 / A; iM(1, 3) = 0;
        iM(0, 4) = 0; iM(1, 4) = iGcp.ground.X;
        iM(0, 5) = 0; iM(1, 5) = iGcp.ground.Y;
        iM(0, 6) = 0; iM(1, 6) = iGcp.ground.Z;
        iM(0, 8) = iGcp.image.x * iGcp.ground.X / A;  iM(1, 8) = iGcp.image.y * iGcp.ground.X / A;
        iM(0, 9) = iGcp.image.x * iGcp.ground.Y / A;  iM(1, 9) = iGcp.image.y * iGcp.ground.Y / A;
        iM(0, 10) = iGcp.image.x * iGcp.ground.Z / A; iM(1, 10) = iGcp.image.y * iGcp.ground.Z / A;
        for (size_t ik = 0; ik < m_nCalibNums; ik++)
        {
            iM(i, 11 + ik) = Dx * pow(r2, ik + 1);
        }
        // ����������
        iW(0) = iGcp.image.x / A;
        iW(1) = iGcp.image.y / A;
        // ����
        MTM += iM.transpose() * iM;
        MTW += iM.transpose() * iW;
    }
    VectorXd vecL(11);
    vecL = MTM.lu().solve(MTW);
    for (size_t i = 0; i < 11; i++) l[i] = vecL(i);
}

