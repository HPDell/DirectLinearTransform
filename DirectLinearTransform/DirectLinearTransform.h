#pragma once

#include <vector>

#include "DataTypes.h"

using namespace std;
class CDirectLinearTransform
{
private:
    vector<GCP> vGcps;
    vector<GCP> vTiePoints;
    // 11����
    ExteriorElements ext;
    InneriorElements inn;
    // �������
    size_t m_nCalibNums;
    double* k;
    // 11��lϵ��
protected:
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
    /// <param name="x">���x����</param>
    /// <param name="y">���y����</param>
    /// <param name="X">�﷽��X����</param>
    /// <param name="Y">�﷽��Y����</param>
    /// <param name="Z">�﷽��Z����</param>
    /// <created>HuYG,2017/6/3</created>
    void AddGcp(size_t id, double x, double y, double X, double Y, double Z)
    {
        GCP gcp = { id, {x, y}, {X,Y,Z} };
        vGcps.push_back(gcp);
    }
    /// <summary>
    /// ����һ�����
    /// </summary>
    /// <param name="id">���</param>
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
    /// ����lϵ��
    /// </summary>
    /// <created>HuYG,2017/6/4</created>
    void Resection();
    /// <summary>
    /// lϵ�����㾫��
    /// </summary>
    /// <param name="pointResidual"></param>
    /// <returns></returns>
    /// <created>HuYG,2017/6/4</created>
    double ResectionPrecision(double* pointResidual = nullptr);
    /// <summary>
    /// ����lϵ�������������
    /// </summary>
    /// <param name="gcp">���������</param>
    /// <returns>�������</returns>
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
    /// ��ȡlϵ��
    /// </summary>
    /// <returns>lϵ��ָ��</returns>
    /// <created>HuYG,2017/6/4</created>
    double* GetL()
    {
        return l;
    }
    /// <summary>
    /// ��ȡ�������
    /// </summary>
    /// <returns>��������ָ��</returns>
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
    /// ����ڷ�λԪ��
    /// </summary>
    /// <created>HuYG,2017/6/3</created>
    void CalcInneriorElements();
    /// <summary>
    /// �����ⷽλԪ��
    /// </summary>
    /// <created>HuYG,2017/6/4</created>
    void CalcExteriorElements();
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
    void Resection_l_Exact();
    /// <summary>
    /// �������в�
    /// </summary>
    /// <param name="gcp"></param>
    /// <returns></returns>
    /// <created>HuYG,2017/6/4</created>
    ImagePoint CalcImageResidual(GCP& gcp);
    /// <summary>
    /// �����������
    /// </summary>
    /// <created>HuYG,2017/6/4</created>
    void RectifyImagePoint();
    void Intersection_SpaceApprox(ImagePoint& lp, ImagePoint& rp, double* ll, double* rl, GroundPoint& gp);
    void Intersection_SpaceExact(ImagePoint& lp, ImagePoint& rp, double* ll, double* rl, GroundPoint& gp);
};

