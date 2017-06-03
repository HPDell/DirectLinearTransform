#pragma once

#include <vector>

#include "DataTypes.h"

using namespace std;
class CDirectLinearTransform
{
private:
    vector<GCP> vGcps;
    // 11����
    ExteriorElements ext;
    InneriorElements inn;
    // �������
    size_t m_nCalibNums;
    double* k;
    // 11��lϵ��
    double l[11];
    //double A;

public:
    CDirectLinearTransform();
    CDirectLinearTransform(size_t nCalibNums);
    ~CDirectLinearTransform();

    /// <summary>
    /// ����һ��GCP
    /// </summary>
    /// <param name="id">���</param>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <param name="X"></param>
    /// <param name="Y"></param>
    /// <param name="Z"></param>
    /// <created>HuYG,2017/6/3</created>
    void AddGcp(size_t id, double x, double y, double X, double Y, double Z)
    {
        GCP gcp = { id, {x, y}, {X,Y,Z} };
        vGcps.push_back(gcp);
    }

    void Resection();

private:
    /// <summary>
    /// ����ڷ�λԪ��
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void CalcInneriorElements();
    /// <summary>
    /// �ռ�󷽽���
    /// ���lϵ����ֵ
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void Resection_l_Approx();
    /// <summary>
    /// �ռ�󷽽���
    /// ����lϵ����ȷֵ
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void Resection_l_Precise();
};

