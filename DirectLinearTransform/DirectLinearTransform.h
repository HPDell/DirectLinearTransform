#pragma once

#include <vector>

#include "DataTypes.h"

using namespace std;
class CDirectLinearTransform
{
private:
    vector<GCP> vGcps;
    // 11参数
    ExteriorElements ext;
    InneriorElements inn;
    // 畸变参数
    size_t m_nCalibNums;
    double* k;
    // 11个l系数
    double l[11];
    //double A;

public:
    CDirectLinearTransform();
    CDirectLinearTransform(size_t nCalibNums);
    ~CDirectLinearTransform();

    /// <summary>
    /// 增加一个GCP
    /// </summary>
    /// <param name="id">编号</param>
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
    /// 求解内方位元素
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void CalcInneriorElements();
    /// <summary>
    /// 空间后方交会
    /// 求解l系数初值
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void Resection_l_Approx();
    /// <summary>
    /// 空间后方交会
    /// 计算l系数精确值
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void Resection_l_Precise();
};

