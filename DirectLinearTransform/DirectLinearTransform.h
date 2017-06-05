#pragma once

#include <vector>

#include "DataTypes.h"

using namespace std;
class CDirectLinearTransform
{
private:
    vector<GCP> vGcps;
    vector<GCP> vTiePoints;
    // 11参数
    ExteriorElements ext;
    InneriorElements inn;
    // 畸变参数
    size_t m_nCalibNums;
    double* k;
    // 11个l系数
protected:
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
    /// <param name="x">像点x坐标</param>
    /// <param name="y">像点y坐标</param>
    /// <param name="X">物方点X坐标</param>
    /// <param name="Y">物方点Y坐标</param>
    /// <param name="Z">物方点Z坐标</param>
    /// <created>HuYG,2017/6/3</created>
    void AddGcp(size_t id, double x, double y, double X, double Y, double Z)
    {
        GCP gcp = { id, {x, y}, {X,Y,Z} };
        vGcps.push_back(gcp);
    }
    /// <summary>
    /// 增加一个像点
    /// </summary>
    /// <param name="id">编号</param>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <param name="X"></param>
    /// <param name="Y"></param>
    /// <param name="Z"></param>
    /// <created>HuYG,2017/6/3</created>
    void AddTiePoint(size_t id, double x, double y, double X, double Y, double Z)
    {
        GCP gcp = { id,{ x, y },{ X,Y,Z } };
        vGcps.push_back(gcp);
    }
    /// <summary>
    /// 解算l系数
    /// </summary>
    /// <created>HuYG,2017/6/4</created>
    void Resection();
    /// <summary>
    /// l系数计算精度
    /// </summary>
    /// <param name="pointResidual"></param>
    /// <returns></returns>
    /// <created>HuYG,2017/6/4</created>
    double ResectionPrecision(double* pointResidual = nullptr);
    /// <summary>
    /// 根据l系数计算像点坐标
    /// </summary>
    /// <param name="gcp">地面点坐标</param>
    /// <returns>像点坐标</returns>
    /// <created>HuYG,2017/6/4</created>
    ImagePoint CalcImagePoint(GroundPoint& gp)
    {
        ImagePoint imagePoint = { 
            inn.x0 - (l[0] * gp.X + l[1] * gp.Y + l[2] * gp.Z + l[3]) / (l[8] * gp.X + l[9] * gp.Y + l[10] * gp.Z + 1),
            inn.y0 - (l[4] * gp.X + l[5] * gp.Y + l[6] * gp.Z + l[7]) / (l[8] * gp.X + l[9] * gp.Y + l[10] * gp.Z + 1)
        };
        double x = inn.x0 - (l[0] * gp.X + l[1] * gp.Y + l[2] * gp.Z + l[3]) / (l[8] * gp.X + l[9] * gp.Y + l[10] * gp.Z + 1);
        double y = inn.y0 - (l[4] * gp.X + l[5] * gp.Y + l[6] * gp.Z + l[7]) / (l[8] * gp.X + l[9] * gp.Y + l[10] * gp.Z + 1);
        return imagePoint;
    }
    void Intersection(CDirectLinearTransform& right);
    static void Intersection(CDirectLinearTransform* dlt, size_t dlt_num);
    /// <summary>
    /// 获取l系数
    /// </summary>
    /// <returns>l系数指针</returns>
    /// <created>HuYG,2017/6/4</created>
    double* GetL()
    {
        return l;
    }
    /// <summary>
    /// 获取畸变参数
    /// </summary>
    /// <returns>基本参数指针</returns>
    /// <created>HuYG,2017/6/4</created>
    double* GetCalibParams()
    {
        return k;
    }
    ExteriorElements GetExteriorElements()
    {
        return ext;
    }
    InneriorElements GetInneriorElements()
    {
        return inn;
    }

private:
    /// <summary>
    /// 求解内方位元素
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void CalcInneriorElements();
    /// <summary>
    /// 计算外方位元素
    /// </summary>
    /// <created>HuYG,2017/6/4</created>
    void CalcExteriorElements();
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
    void Resection_l_Exact();
    /// <summary>
    /// 计算像点残差
    /// </summary>
    /// <param name="gcp"></param>
    /// <returns></returns>
    /// <created>HuYG,2017/6/4</created>
    ImagePoint CalcImageResidual(GCP& gcp);
    /// <summary>
    /// 纠正像点坐标
    /// </summary>
    /// <created>HuYG,2017/6/4</created>
    void RectifyImagePoint();
    void Intersection_SpaceApprox(ImagePoint& lp, ImagePoint& rp, double* ll, double* rl, GroundPoint& gp);
    void Intersection_SpaceExact(ImagePoint& lp, ImagePoint& rp, double* ll, double* rl, GroundPoint& gp);
};

