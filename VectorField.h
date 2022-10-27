#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

struct Vector
{
    double vx,vy,vz;

    void operator+=(const Vector& v)
    {
        this->vx += v.vx;
        this->vy += v.vy;
        this->vz += v.vz;
    }
};

Vector operator+(const Vector& v1,const Vector& v2)
{
    return Vector{v1.vx+v2.vx,v1.vy+v2.vy,v1.vz+v2.vz};
}

Vector operator-(const Vector& v1,const Vector& v2)
{
    return Vector{v1.vx-v2.vx,v1.vy-v2.vy,v1.vz-v2.vz};
}

Vector operator*(const double& d,const Vector& v)
{
    return Vector{d*v.vx,d*v.vy,d*v.vz};
}

struct VectorField
{
    double xmin,ymin,zmin,xmax,ymax,zmax;
    double step;
    int ximax,yimax,zimax;

    vector<vector<vector<Vector>>> field;

    VectorField(double xmin,double xmax,double ymin, double ymax, double zmin, double zmax, double step)
    : xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax),zmin(zmin),zmax(zmax),step(step)
    {        
        for (int i = 0; i <= round((xmax-xmin)/step); i++)
        {
            vector<vector<Vector>> v1;
            for (int j = 0; j <= round((ymax-ymin)/step); j++)
            {
                vector<Vector> v2;
                for (int k = 0; k <= round((zmax-zmin)/step); k++)
                {
                    v2.push_back({0,0,0});
                }
                v1.push_back(v2);
            }
            field.push_back(v1);
        }

        this->GetVectorIndexes(xmax,ymax,zmax,ximax,yimax,zimax);
    }

    void GetVectorIndexes(double x, double y, double z, int& xi, int& yi, int& zi)
    {
        if(x<xmin||x>xmax||y<ymin||y>ymax||z<zmin||z>zmax)
            cerr << "Cannot read field out of bounds.\n";
        
        xi = round((x-xmin)/step);
        yi = round((y-ymin)/step);
        zi = round((z-zmin)/step);
    }

    Vector* GetVector(double x, double y, double z)
    {
        int xi,yi,zi;
        this->GetVectorIndexes(x,y,z,xi,yi,zi);

        return &(field[xi][yi][zi]);
    }

    void LoadField(const char* filename)
    {
        std::ifstream inf {filename};
        cout << "\nLoading field from " << filename << "\n";

        int lines_read = 0;
        int lines_processed = 0;
        int lines_expected = round((xmax-xmin)/step+1)*round((ymax-ymin)/step+1)*round((zmax-zmin)/step+1);

        while (inf)
        {
            std::string X,Y,Z,VX,VY,VZ;
            inf >> X; inf >> Y; inf >> Z; inf >> VX; inf >> VY; inf >> VZ;
            lines_read++;
            try
            {                
                if (X != "")
                {
                    double x,y,z,vx,vy,vz;
                    x = stod(X); y = stod(Y); z = stod(Z);
                    vx = stod(VX); vy = stod(VY); vz = stod(VZ);

                    *(this->GetVector(x,y,z)) = Vector{vx,vy,vz};
                    lines_processed++;
                }
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                cout << "X: " << X << " Y: " << Y << " Z: " << Z << " VX: " << VX << " VY: " << VY << " VZ: " << VZ << "\n";
            }
        }

        cout << "Lines read: " << lines_read << " processed: " << lines_processed << " expected: " << lines_expected << "\n\n";
    }

    Vector GetField(double x, double y, double z)
    {
        int xi,yi,zi;
        this->GetVectorIndexes(x,y,z,xi,yi,zi);
        
        int xi2,yi2,zi2;
        if(x-(xmin+step*xi)<0) xi2 = xi-1; else xi2 = xi + 1; if(xi2>ximax) xi2=xi;
        if(y-(ymin+step*yi)<0) yi2 = yi-1; else yi2 = yi + 1; if(yi2>yimax) yi2=yi;
        if(z-(zmin+step*zi)<0) zi2 = zi-1; else zi2 = zi + 1; if(zi2>zimax) zi2=zi;

        double dx = abs((x-(xmin+step*xi))/step);
        double dy = abs((y-(ymin+step*yi))/step);
        double dz = abs((z-(zmin+step*zi))/step);

        Vector& c000 = field[xi][yi][zi];
        Vector& c001 = field[xi][yi][zi2];
        Vector& c010 = field[xi][yi2][zi];
        Vector& c011 = field[xi][yi2][zi2];
        Vector& c100 = field[xi2][yi][zi];
        Vector& c101 = field[xi2][yi][zi2];
        Vector& c110 = field[xi2][yi2][zi];
        Vector& c111 = field[xi2][yi2][zi2];

        Vector c00 = (1-dx)*c000+dx*c100;
        Vector c01 = (1-dx)*c001+dx*c101;
        Vector c10 = (1-dx)*c010+dx*c110;
        Vector c11 = (1-dx)*c011+dx*c111;

        Vector c0 = (1-dy)*c00+dy*c10;
        Vector c1 = (1-dy)*c01+dy*c11;

        return (1-dz)*c0+dz*c1;
    }
};


struct SensorData
{
    double x1,y1,z1,t1;
    double x1dev,y1dev,z1dev,t1dev;

    void operator+=(const SensorData& s)
    {
        this->x1 += s.x1;
        this->y1 += s.y1;
        this->z1 += s.z1;
        this->t1 += s.t1;
    }   
};

SensorData operator+(const SensorData& s1,const SensorData& s2)
{
    return SensorData{s1.x1+s2.x1,s1.y1+s2.y1,s1.z1+s2.z1,s1.t1+s2.t1,0,0,0,0};
}

SensorData operator-(const SensorData& s1,const SensorData& s2)
{
    return SensorData{s1.x1-s2.x1,s1.y1-s2.y1,s1.z1-s2.z1,s1.t1-s2.t1,0,0,0,0};
}

SensorData operator*(const double& d,const SensorData& s)
{
    return SensorData{d*s.x1,d*s.y1,d*s.z1,d*s.t1,0,0,0,0};
}

template<typename T>
struct Field
{
    double xmin,ymin,zmin,xmax,ymax,zmax;
    double step;
    T def;

    int ximax,yimax,zimax;
    vector<vector<vector<T>>> field;

    void SetDefault(T def) {this->def = def;}
    void InitField(double xmin,double xmax,double ymin, double ymax, double zmin, double zmax, double step)
    {
        this->xmin = xmin;
        this->xmax = xmax;
        this->ymin = ymin;
        this->ymax = ymax;
        this->zmin = zmin;
        this->zmax = zmax;
        this->step = step;

        for (int i = 0; i <= round((xmax-xmin)/step); i++)
        {
            vector<vector<T>> v1;
            for (int j = 0; j <= round((ymax-ymin)/step); j++)
            {
                vector<T> v2;
                for (int k = 0; k <= round((zmax-zmin)/step); k++)
                {
                    v2.push_back(def);
                }
                v1.push_back(v2);
            }
            field.push_back(v1);
        }

        this->GetPointIndexes(xmax,ymax,zmax,ximax,yimax,zimax);
    }
    
    void GetPointIndexes(double x, double y, double z, int& xi, int& yi, int& zi)
    {
        if(x<xmin||x>xmax||y<ymin||y>ymax||z<zmin||z>zmax)
            cerr << "Cannot read field out of bounds.\n";

        if (x < xmin) xi = 0;
        else if (x > xmax) xi = ximax;
        if (y < ymin) yi = 0;
        else if (y > ymax) yi = yimax;
        if (z < zmin) zi = 0;
        else if (z > zmax) zi = zimax;        
        
        xi = round((x-xmin)/step);
        yi = round((y-ymin)/step);
        zi = round((z-zmin)/step);
    }

    T* GetPoint(double x, double y, double z)
    {
        int xi,yi,zi;
        this->GetPointIndexes(x,y,z,xi,yi,zi);

        return &(field[xi][yi][zi]);
    }

    void SetPoint(double x, double y, double z, T new_value)
    {
        *(this->GetPoint(x,y,z)) = new_value;
    }

    T GetField(double x, double y, double z)
    {
        int xi,yi,zi;
        this->GetPointIndexes(x,y,z,xi,yi,zi);
        
        int xi2,yi2,zi2;
        if(x-(xmin+step*xi)<0) xi2 = xi-1; else xi2 = xi + 1; if(xi2>ximax) xi2=xi;
        if(y-(ymin+step*yi)<0) yi2 = yi-1; else yi2 = yi + 1; if(yi2>yimax) yi2=yi;
        if(z-(zmin+step*zi)<0) zi2 = zi-1; else zi2 = zi + 1; if(zi2>zimax) zi2=zi;

        double dx = abs((x-(xmin+step*xi))/step);
        double dy = abs((y-(ymin+step*yi))/step);
        double dz = abs((z-(zmin+step*zi))/step);

        T& c000 = field[xi][yi][zi];
        T& c001 = field[xi][yi][zi2];
        T& c010 = field[xi][yi2][zi];
        T& c011 = field[xi][yi2][zi2];
        T& c100 = field[xi2][yi][zi];
        T& c101 = field[xi2][yi][zi2];
        T& c110 = field[xi2][yi2][zi];
        T& c111 = field[xi2][yi2][zi2];

        T c00 = (1-dx)*c000+dx*c100;
        T c01 = (1-dx)*c001+dx*c101;
        T c10 = (1-dx)*c010+dx*c110;
        T c11 = (1-dx)*c011+dx*c111;

        T c0 = (1-dy)*c00+dy*c10;
        T c1 = (1-dy)*c01+dy*c11;

        return (1-dz)*c0+dz*c1;
    }
};