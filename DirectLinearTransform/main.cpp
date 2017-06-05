// DirectLinearTransform.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "DirectLinearTransform.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#define PIXELSIZE 0.0051966

//#define RESECTION_REPORT
#define INTERSECTION
#define INTERSECTION_REPORT

using namespace std;

int main(int argc, char** argv)
{
    // 左片求解L系数
    CDirectLinearTransform dlt_L(2);
    ifstream fin_gcp_L(argv[1]);
    if (!fin_gcp_L.is_open())
    {
        cout << "打开控制点文件失败" << endl;
        return -1;
    }
    int nGcpNum_L = 0; // 控制点数量
    fin_gcp_L >> nGcpNum_L;
    vector<size_t> vPointId_L;
    for (size_t i = 0; i < nGcpNum_L && !fin_gcp_L.eof(); i++)
    {
        size_t gcp_id;
        double gcp_coord[5];
        fin_gcp_L >> gcp_id >> gcp_coord[0] >> gcp_coord[1] >> gcp_coord[2] >> gcp_coord[3] >> gcp_coord[4];
        vPointId_L.push_back(gcp_id);
        dlt_L.AddGcp(gcp_id, gcp_coord[0], gcp_coord[1], gcp_coord[2], gcp_coord[3], gcp_coord[4]);
    }
    fin_gcp_L.close();

    dlt_L.Resection();

#ifdef RESECTION_REPORT
    // 平差报告结构
    ResectionReport report = dlt_L.GetResectionReport();
    // 输出结果
    cout << "直接线性变换l系数解算报告" << endl;

    cout << endl << "内方位元素" << endl;
    InneriorElements& inn = report.inn_rep;
    cout << setw(3) << left << "fx" << right << fixed << setw(14) << setprecision(3) << inn.fx * PIXELSIZE << endl;
    cout << setw(3) << left << "fy" << right << fixed << setw(14) << setprecision(3) << inn.fy * PIXELSIZE << endl;
    cout << setw(3) << left << "x0" << right << fixed << setw(14) << setprecision(3) << inn.x0 << endl;
    cout << setw(3) << left << "y0" << right << fixed << setw(14) << setprecision(3) << inn.y0 << endl;
    cout << setw(3) << left << "ds" << right << scientific << setw(14) << setprecision(3) << inn.ds << endl;
    cout << setw(3) << left << "db" << right << scientific << setw(14) << setprecision(3) << inn.db << endl;
    cout << endl << "外方位元素" << endl;
    ExteriorElements& ext = report.ext_rep;
    cout << setw(6) << left << "Xs" << right << fixed << setw(10) << setprecision(3) << ext.Xs << endl;
    cout << setw(6) << left << "Ys" << right << fixed << setw(10) << setprecision(3) << ext.Ys << endl;
    cout << setw(6) << left << "Zs" << right << fixed << setw(10) << setprecision(3) << ext.Zs << endl;
    cout << setw(6) << left << "phi" << right << fixed << setw(10) << setprecision(3) << ext.phi << endl;
    cout << setw(6) << left << "omega" << right << fixed << setw(10) << setprecision(3) << ext.omega << endl;
    //cout << setw(6) << left << "kappa" << right << fixed << setw(8) << setprecision(3) << ext.kappa << endl;
    cout << endl << "l系数" << endl;
    vector<double>& l = report.l;
    for (size_t i = 0; i < 11; i++)
    {
        cout << "l" << setw(2) << left << (i + 1) << setw(14) << right << fixed << setprecision(6) << l[i] << endl;
    }
    cout << endl << "畸变参数" << endl;
    vector<double>& k = report.calib_params;
    for (size_t i = 0; i < 2; i++)
    {
        cout << "k" << setw(2) << left << (i + 1) << setw(16) << right << scientific << setprecision(6) << k[i] << endl;
    }
    //double* pResiduals = new double[nGcpNum_L];
    vector<ImagePointResidual>& residual = report.image_point_residual;
    cout << endl << "像点中误差: " << setprecision(6) << fixed << report.m0 << endl;
    for (size_t i = 0; i < nGcpNum_L; i++)
    {
        cout << setw(4) << left << residual[i].id << "号点：" << fixed << setprecision(3)
            << left << setw(4) << "x:" << right << setw(6) << residual[i].vx << " 像素"
            << left << setw(4) << "y:" << right << setw(6) << residual[i].vy << " 像素"
            << left << setw(4) << "xy:" << right << setw(6) << residual[i].vxy << " 像素" << endl;
    }
    //system("pause");
#endif // RESECTION_REPORT

#ifdef INTERSECTION
    // 右片求解l系数
    CDirectLinearTransform dlt_R(2);
    ifstream fin_gcp_R(argv[2]);
    if (!fin_gcp_R.is_open())
    {
        cout << "打开控制点文件失败" << endl;
        return -1;
    }
    int nGcpNum_R = 0; // 控制点数量
    fin_gcp_R >> nGcpNum_R;
    vector<size_t> vPointId_R;
    for (size_t i = 0; i < nGcpNum_R && !fin_gcp_R.eof(); i++)
    {
        size_t gcp_id;
        double gcp_coord[5];
        fin_gcp_R >> gcp_id >> gcp_coord[0] >> gcp_coord[1] >> gcp_coord[2] >> gcp_coord[3] >> gcp_coord[4];
        vPointId_R.push_back(gcp_id);
        dlt_R.AddGcp(gcp_id, gcp_coord[0], gcp_coord[1], gcp_coord[2], gcp_coord[3], gcp_coord[4]);
    }
    fin_gcp_R.close();

    dlt_R.Resection();

    // 读取物方点坐标
    ifstream fin_tie(argv[3]);
    if (!fin_tie.is_open())
    {
        cout << "打开连接点文件失败" << endl;
        return -1;
    }
    int nTiePoint = 0;
    fin_tie >> nTiePoint;
    vector<size_t> vTieId;
    for (size_t i = 0; i < nTiePoint && !fin_tie.eof(); i++)
    {
        size_t tie_id;
        double tie_coord[7];
        fin_tie >> tie_id >> tie_coord[0] >> tie_coord[1] >> tie_coord[2] >> tie_coord[3] >> tie_coord[4] >> tie_coord[5] >> tie_coord[6];
        vPointId_R.push_back(tie_id);
        //dlt_R.AddGcp(tie_id, tie_coord[0], tie_coord[1], tie_coord[2], tie_coord[3], tie_coord[4]);
        dlt_L.AddTiePoint(tie_id, tie_coord[0], tie_coord[1], 0.0, 0.0, 0.0);
        dlt_R.AddTiePoint(tie_id, tie_coord[2], tie_coord[3], tie_coord[4], tie_coord[5], tie_coord[6]);
    }

    // 求解物方坐标精确值
    dlt_L.Intersection(dlt_R);
    //system("pause");
#endif // INTERSECTION

#ifdef INTERSECTION_REPORT
    IntersectionReport intersection = dlt_L.GetIntersectionReport(dlt_R);
    cout << "直接线性变换空间前方交会报告" << endl;
    cout << endl << "检查点" << setw(2) << "XY" << "中误差：" 
        << fixed << right << setw(10) << setprecision(3) << intersection.mXY << endl;
    cout << endl << "检查点" << setw(2) << "Z" << "中误差：" 
        << fixed << right << setw(10) << setprecision(3) << intersection.mZ << endl;
    cout << endl << "检查点像点残差" << endl;
    for each (GroundPointResidual ground in intersection.ground_point_residual)
    {
        cout << setw(4) << left << ground.id << "号点：" << fixed << setprecision(3)
            << right << setw(4) << "X:" << right << setw(10) << ground.vX << " mm"
            << right << setw(4) << "Y:" << right << setw(10) << ground.vY << " mm"
            << right << setw(4) << "Z:" << right << setw(10) << ground.vZ << " mm" << endl;
    }
    //system("pause");
#endif // INTERSECTION_REPORT

    return 0;
}

