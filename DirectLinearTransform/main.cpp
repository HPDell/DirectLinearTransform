// DirectLinearTransform.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "DirectLinearTransform.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#define PIXELSIZE 0.0051966

#define INTERSECTION    

using namespace std;

int main(int argc, char** argv)
{
    // ��Ƭ���Lϵ��
    CDirectLinearTransform dlt_L(2);
    ifstream fin_gcp_L(argv[1]);
    if (!fin_gcp_L.is_open())
    {
        cout << "�򿪿��Ƶ��ļ�ʧ��" << endl;
        return -1;
    }
    int nGcpNum_L = 0; // ���Ƶ�����
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
    // lϵ��
    double* l = dlt_L.GetL();
    double* k = dlt_L.GetCalibParams();
    // ������
    cout << "ֱ�����Ա任lϵ�����㱨��" << endl;

    cout << endl << "�ڷ�λԪ��" << endl;
    InneriorElements inn = dlt_L.GetInneriorElements();
    cout << setw(3) << left << "fx" << right << fixed << setw(14) << setprecision(3) << inn.fx * PIXELSIZE << endl;
    cout << setw(3) << left << "fy" << right << fixed << setw(14) << setprecision(3) << inn.fy * PIXELSIZE << endl;
    cout << setw(3) << left << "x0" << right << fixed << setw(14) << setprecision(3) << inn.x0 << endl;
    cout << setw(3) << left << "y0" << right << fixed << setw(14) << setprecision(3) << inn.y0 << endl;
    cout << setw(3) << left << "ds" << right << scientific << setw(14) << setprecision(3) << inn.ds << endl;
    cout << setw(3) << left << "db" << right << scientific << setw(14) << setprecision(3) << inn.db << endl;
    cout << endl << "�ⷽλԪ��" << endl;
    ExteriorElements ext = dlt_L.GetExteriorElements();
    cout << setw(6) << left << "Xs" << right << fixed << setw(10) << setprecision(3) << ext.Xs << endl;
    cout << setw(6) << left << "Ys" << right << fixed << setw(10) << setprecision(3) << ext.Ys << endl;
    cout << setw(6) << left << "Zs" << right << fixed << setw(10) << setprecision(3) << ext.Zs << endl;
    cout << setw(6) << left << "phi" << right << fixed << setw(10) << setprecision(3) << ext.phi << endl;
    cout << setw(6) << left << "omega" << right << fixed << setw(10) << setprecision(3) << ext.omega << endl;
    //cout << setw(6) << left << "kappa" << right << fixed << setw(8) << setprecision(3) << ext.kappa << endl;
    cout << endl << "lϵ��" << endl;
    for (size_t i = 0; i < 11; i++)
    {
        cout << "l" << setw(2) << left << (i + 1) << setw(14) << right << fixed << setprecision(6) << l[i] << endl;
    }
    cout << endl << "�������" << endl;
    for (size_t i = 0; i < 2; i++)
    {
        cout << "k" << setw(2) << left << (i + 1) << setw(12) << right << scientific << setprecision(3) << k[i] << endl;
    }
    double* pResiduals = new double[nGcpNum_L];
    cout << endl << "��������" << setprecision(6) << fixed << dlt_L.ResectionPrecision(pResiduals) << endl;
    for (size_t i = 0; i < nGcpNum_L; i++)
    {
        cout << setw(4) << left << vPointId_L[i] << "�ŵ㣺"
            << right << setw(8) << fixed << setprecision(3) << pResiduals[i] << " ����" << endl;
    }
#endif // RESECTION_REPORT

#ifdef INTERSECTION
    // ��Ƭ���lϵ��
    CDirectLinearTransform dlt_R(2);
    ifstream fin_gcp_R(argv[2]);
    if (!fin_gcp_R.is_open())
    {
        cout << "�򿪿��Ƶ��ļ�ʧ��" << endl;
        return -1;
    }
    int nGcpNum_R = 0; // ���Ƶ�����
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

    // ��ȡ�﷽������
    ifstream fin_tie(argv[3]);
    if (!fin_tie.is_open())
    {

    }

    // ����﷽���꾫ȷֵ
    dlt_L.Intersection(dlt_R);
#endif // INTERSECTION

    system("pause");
    return 0;
}

