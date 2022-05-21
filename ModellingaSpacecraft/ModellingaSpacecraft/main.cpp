#include "mainwindow.h"

#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

double R0 = 6371;

class Coordinates{
    public:
    double x, y, z, vx, vy, vz;
    Coordinates(){
        x = 0;
        y = 0;
        z = 0;
        vx = 0;
        vy = 0;
        vz = 0;
    }
    Coordinates(double x0, double y0, double z0, double vx0, double vy0, double vz0){
        x = x0;
        y = y0;
        z = z0;
        vx = vx0;
        vy = vy0;
        vz = vz0;
    }
    void print(){
        printf("x = %f\t", x);
        printf("y = %f\t", y);
        printf("z = %f\t", z);
        printf("vx = %f\t", vx);
        printf("vy = %f\t", vy);
        printf("vz = %f\t", vz);
    }
};
Coordinates operator +(Coordinates a, Coordinates b)
{
    Coordinates c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    c.vx = a.vx + b.vx;
    c.vy = a.vy + b.vy;
    c.vz = a.vz + b.vz;
    return c;
}
Coordinates operator -(Coordinates a, Coordinates b)
{
    Coordinates c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    c.vx = a.vx - b.vx;
    c.vy = a.vy - b.vy;
    c.vz = a.vz - b.vz;
    return c;
}
Coordinates operator *(Coordinates a, double b)
{
    Coordinates c;
    c.x = b*a.x;
    c.y = b*a.y;
    c.z = b*a.z;
    c.vx = b*a.vx;
    c.vy = b*a.vy;
    c.vz = b*a.vz;
    return c;
}
Coordinates operator *(double a, Coordinates b)
{
    Coordinates c;
    c.x = a*b.x;
    c.y = a*b.y;
    c.z = a*b.z;
    c.vx = a*b.vx;
    c.vy = a*b.vy;
    c.vz = a*b.vz;
    return c;
}


Coordinates F(Coordinates Sc, double t)
{
    Coordinates res;
    double mu = 398600.4415;
    res.x = Sc.vx;
    res.y = Sc.vy;
    res.z = Sc.vz;
    res.vx = (-1)*mu*Sc.x/pow(sqrt(Sc.x*Sc.x + Sc.y*Sc.y + Sc.z*Sc.z), 3);
    res.vy = (-1)*mu*Sc.y/pow(sqrt(Sc.x*Sc.x + Sc.y*Sc.y + Sc.z*Sc.z), 3);
    res.vz = (-1)*mu*Sc.z/pow(sqrt(Sc.x*Sc.x + Sc.y*Sc.y + Sc.z*Sc.z), 3);
    return res;
}

double norm(vector<double> a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

vector<Coordinates> RK4(double x0, double y0, double z0, double vx0, double vy0, double vz0, double t0, double dt, double t)
{
    vector<Coordinates> XYZ;
    Coordinates Sc(x0, y0, z0, vx0, vy0, vz0);
    XYZ.push_back(Sc);
    Coordinates k1, k2, k3, k4;
    for(int i = t0; i<t; i++)
    {
        k1 = F(Sc, t0)*dt;
        k2 = F(Sc + 0.5*k1, t0 + 0.5*dt)*dt;
        k3 = F(Sc + 0.5*k2, t0 + 0.5*dt)*dt;
        k4 = F(Sc + k3, t0 + dt)*dt;
        Sc = Sc + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        XYZ.push_back(Sc);
    }
    return XYZ;
}

vector<double> Apoghei(vector<Coordinates> a)
{
    double dist = 0;
    double A_h = 0;
    double x, y, z, vx, vy, vz, v;
    vector<double> r_a;
    for(int i = 0; i<a.size(); i++)
    {
        dist = sqrt(a[i].x*a[i].x + a[i].y*a[i].y + a[i].z*a[i].z);
        if (A_h <= dist)
        {
            A_h = dist;
            x = a[i].x;
            y = a[i].y;
            z = a[i].z;
            vx = a[i].vx;
            vy = a[i].vy;
            vz = a[i].vz;
        }
    }
    r_a.push_back(x);
    r_a.push_back(y);
    r_a.push_back(z);
    return r_a;
}

vector<double> Perighei(vector<Coordinates> a)
{
    double dist = 0;
    double P_h = 1000000;
    double x, y, z, vx, vy, vz, v;
    vector<double> r_p;
    for(int i = 0; i<a.size(); i++)
    {
        dist = sqrt(a[i].x*a[i].x + a[i].y*a[i].y + a[i].z*a[i].z);
        if (P_h >= dist)
        {
            P_h = dist;
            x = a[i].x;
            y = a[i].y;
            z = a[i].z;
            vx = a[i].vx;
            vy = a[i].vy;
            vz = a[i].vz;
        }
    }
    r_p.push_back(x);
    r_p.push_back(y);
    r_p.push_back(z);
    return r_p;
}

void Orbit_parameters(vector<Coordinates> s)
{
    double a = (norm(Apoghei(s)) - R0 + norm(Perighei(s)) - R0)/2;
    double e = (norm(Apoghei(s)) - norm(Perighei(s)))/(norm(Apoghei(s)) - R0 + norm(Perighei(s)) - R0);
    double x1 = s[0].x;
    double x2 = s[1000].x;
    double x3 = s[2000].x;
    double y1 = s[0].y;
    double y2 = s[1000].y;
    double y3 = s[2000].y;
    double z1 = s[0].z;
    double z2 = s[1000].z;
    double z3 = s[2000].z;
    double nx = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
    double ny = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1);
    double nz = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
    vector<double> n;
    n.push_back(nx);
    n.push_back(ny);
    n.push_back(nz);
    nx = nx/(norm(n));
    ny = ny/(norm(n));
    nz = nz/(norm(n));
    n.push_back(nx);
    n.push_back(ny);
    n.push_back(nz);
    vector<double> h;
    double tt1;
    tt1 = s[0].y*s[0].vz - s[0].z*s[0].vy;
    h.push_back(tt1);
    tt1 = s[0].z*s[0].vx - s[0].x*s[0].vz;
    h.push_back(tt1);
    tt1 = s[0].x*s[0].vy - s[0].y*s[0].vx;
    h.push_back(tt1);
    vector<double> nn;
    nn.push_back(-h[1]);
    nn.push_back(h[0]);
    nn.push_back(0);
    double Omega;
    if(norm(nn) == 0)
        Omega = 0;
    else
    {
    if(nn[1] >= 0)
        Omega = acos(nn[0]/norm(nn));
        else
            Omega = 8*atan(1) - acos(nn[0]/norm(nn));
    }
    double i = acos(nz);
    double px = ny;
    double py = -nx;
    double pz = 0;
    double p_foc = a*(1 - e*e);
    vector<double> p;
    p.push_back(px);
    p.push_back(py);
    p.push_back(pz);
    vector<double> r_p = Perighei(s);
    vector<double> r_a = Apoghei(s);
    r_p[0] = r_p[0];
    r_p[1] = r_p[1];
    r_p[2] = r_p[2];
    double omega = acos(p[0]*r_p[0]/norm(r_p) + p[1]*r_p[1]/norm(r_p) + p[2]*r_p[2]/norm(r_p));
    cout << endl;
    cout << "PARAMETERS OF THE ORBIT:" << endl;
    cout << "------------------------" << endl;
    cout << "a = " << a << endl;
    cout << "p = " << p_foc << endl;
    cout << "e = " << e << endl;
    cout << endl;
    cout << "PERIGEE" << endl;
    cout << "Height: " << norm(r_p)-R0 << "  km" << endl;
    cout << "Coordinates: x = " << r_p[0] << ";  y = " << r_p[1] << ";  z = " << r_p[2] << endl;
    cout << endl;
    cout << "APOGEE" << endl;
    cout << "Height: " << norm(r_a)-R0 << "  km" << endl;
    cout << "Coordinates: x = " << r_a[0] << ";  y = " << r_a[1] << ";  z = " << r_a[2] << endl;
    cout << endl;
    cout << "Inclination i = " << i << " rad" << endl;
    cout << "Longitude of the ascending node: Omega = " << Omega << " rad" << endl;
    cout << "Pericentr argument:  omega = " << omega << " rad" << endl;
    cout << "Normal: nx = " << nx << "    ny = " << ny << "    nz = " << nz << endl;
    cout << endl;
}

void recording_data2(vector<Coordinates> D)
{
    ofstream f("Orbit_Data_new.txt");
    for(int i=0; i<D.size(); i++)
    {
        f << D[i].x << " " << D[i].y << " " << D[i].z << " " << D[i].vx <<  " " << D[i].vy << " " << D[i].vz << endl;
    }
}

void IncreasingApoghei(vector<Coordinates> Sc, double H)
{
    double u;
    double mu = 398600.4415;
    double r1 = norm(Perighei(Sc));
    double r2 = norm(Apoghei(Sc));
    u = sqrt(2*mu/r1)*(sqrt((r2 + H)/(r1 + r2 + H)) - sqrt(r2/(r1 + r2)));
    cout << endl;
    cout << "To increase the height of the apoghei you need to increse velocity in  u = " << u << " km/s  times" << endl;
    cout << endl;
}

void recording_data(vector<Coordinates> D)
{
    ofstream f("Orbit_Data.txt");
    for(int i=0; i<D.size(); i++)
    {
        f << D[i].x << " " << D[i].y << " " << D[i].z << " " << D[i].vx <<  " " << D[i].vy << " " << D[i].vz << endl;
    }
}


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QTranslator translator;
    const QStringList uiLanguages = QLocale::system().uiLanguages();
    for (const QString &locale : uiLanguages) {
        const QString baseName = "ModellingaSpacecraft_" + QLocale(locale).name();
        if (translator.load(":/i18n/" + baseName)) {
            a.installTranslator(&translator);
            break;
        }
    }
    MainWindow w;
    w.show();

    double t0 = 0;
    double dt = 1;
    double t = 20000;
    double x0 = R0 + 2000;
    double y0 = 0;
    double z0 = 0;
    double vx0 = 0;
    double vy0 = 8;
    double vz0 = 0;
    vector<Coordinates> Sc, Sc2;
    Sc = RK4(x0, y0, z0, vx0, vy0, vz0, t0, dt, t);
    vector<double> P = Perighei(Sc);
    vector<double> A = Apoghei(Sc);
    Orbit_parameters(Sc);
    double H = 5000;
    IncreasingApoghei(Sc, H);
    recording_data(Sc);
    double u = 0.31362;
    Sc2 = Sc = RK4(x0, y0, z0, vx0, vy0 + u, vz0, t0, dt, t);
    recording_data2(Sc2);
    return a.exec();
}
