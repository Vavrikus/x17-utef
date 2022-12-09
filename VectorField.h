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

    double SqMagnitude() const {return vx*vx+vy*vy+vz*vz;}
    double Magnitude() const {return sqrt(this->SqMagnitude());}
    double Angle(const Vector& other) const
    {
        return acos((vx*other.vx+vy*other.vy+vz*other.vz)/(this->Magnitude()*other.Magnitude()));
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

        if(xi2 < 0) xi2 = 0; if(yi2 < 0) yi2 = 0; if(zi2 < 0) zi2 = 0;

        // trilinear interpolation
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

    T Invert(double,double,double){return nullptr;} // Only for SensorData
};

template<>
SensorData Field<SensorData>::Invert(double x1, double z1, double t1)
{
    // Find 8 closest points using binary search (assuming ordering)
        int ximin = 0;               // Starting minimal search x index
        int ximax = this->ximax;     // Starting maximal search x index
        int yimin = this->yimax;     // Starting minimal search y index
        int yimax = 0;               // Starting maximal search y index
        int zimin = 0;               // Starting minimal search z index
        int zimax = this->zimax;     // Starting maximal search z index

        int  ximid = (ximin+ximax)/2; // x midpoint variable
        int  yimid = (yimin+yimax)/2; // y midpoint variable
        int  zimid = (zimin+zimax)/2; // z midpoint variable
        bool ximid_is_max,yimid_is_max,zimid_is_max; // Do midpoint variables contain higher value?

        // z index search
        while ((zimin+1<zimax))
        {
            zimid = (zimin+zimax)/2;

            double z1low  = field[ximid][yimid][zimin].z1;
            double z1high = field[ximid][yimid][zimax].z1;
            double z1mid  = field[ximid][yimid][zimid].z1;

            if(z1mid < z1low || z1mid > z1high) cerr << "WARNING: Ordering mistake in binary search (z).\n";
            if(z1 <= z1mid) {zimax = zimid; zimid_is_max = true;}
            if(z1 >  z1mid) {zimin = zimid; zimid_is_max = false;}
        }

        // Some part of the field is initialized with zeros and isn't therefore ordered
        if(zimid != 0)
        {
            while(field[ximin][yimid][zimid].x1 == 0) ximin++;
            while(field[ximax][yimid][zimid].x1 == 0) ximax--;
        }

        // x index search
        while ((ximin+1<ximax))
        {
            ximid = (ximin+ximax)/2;

            double x1low  = field[ximin][yimid][zimid].x1;
            double x1high = field[ximax][yimid][zimid].x1;
            double x1mid  = field[ximid][yimid][zimid].x1;

            if(x1mid < x1low || x1mid > x1high) cerr << "WARNING: Ordering mistake in binary search (x).\n";
            if(x1 <= x1mid) {ximax = ximid; ximid_is_max = true;}
            if(x1 >  x1mid) {ximin = ximid; ximid_is_max = false;}
        }

        // y index search
        while ((yimin>yimax+1))
        {
            yimid = (yimin+yimax)/2;

            double t1low  = field[ximid][yimin][zimid].t1;
            double t1high = field[ximid][yimax][zimid].t1;
            double t1mid  = field[ximid][yimid][zimid].t1;

            if(t1mid < t1low || t1mid > t1high) cerr << "WARNING: Ordering mistake in binary search (y).\n";
            if(t1 <= t1mid) {yimax = yimid; yimid_is_max = true;}
            if(t1 >  t1mid) {yimin = yimid; yimid_is_max = false;}
        }

        // Set up actual min/max variables for the cube
        if(ximid_is_max) {ximax = ximid; ximin = ximax-1;}
        else             {ximin = ximid; ximax = ximin+1;}
        if(yimid_is_max) {yimax = yimid; yimin = yimax+1;}
        else             {yimin = yimid; yimax = yimin-1;}
        if(zimid_is_max) {zimax = zimid; zimin = zimax-1;}
        else             {zimin = zimid; zimax = zimin+1;}

        // Sanity check
            cout << "Cube for (x1,z1,t1) = (" << x1 << "," << z1 << "," << t1 << "): \n";
            cout << "000: [" << ximin << "][" << yimin << "][" << zimin << "], (x,z,t) = (" << field[ximin][yimin][zimin].x1 << "," << field[ximin][yimin][zimin].z1 << "," << field[ximin][yimin][zimin].t1 << ")\n";
            cout << "001: [" << ximin << "][" << yimin << "][" << zimax << "], (x,z,t) = (" << field[ximin][yimin][zimax].x1 << "," << field[ximin][yimin][zimax].z1 << "," << field[ximin][yimin][zimax].t1 << ")\n";
            cout << "010: [" << ximin << "][" << yimax << "][" << zimin << "], (x,z,t) = (" << field[ximin][yimax][zimin].x1 << "," << field[ximin][yimax][zimin].z1 << "," << field[ximin][yimax][zimin].t1 << ")\n";
            cout << "011: [" << ximin << "][" << yimax << "][" << zimax << "], (x,z,t) = (" << field[ximin][yimax][zimax].x1 << "," << field[ximin][yimax][zimax].z1 << "," << field[ximin][yimax][zimax].t1 << ")\n";
            cout << "100: [" << ximax << "][" << yimin << "][" << zimin << "], (x,z,t) = (" << field[ximax][yimin][zimin].x1 << "," << field[ximax][yimin][zimin].z1 << "," << field[ximax][yimin][zimin].t1 << ")\n";
            cout << "101: [" << ximax << "][" << yimin << "][" << zimax << "], (x,z,t) = (" << field[ximax][yimin][zimax].x1 << "," << field[ximax][yimin][zimax].z1 << "," << field[ximax][yimin][zimax].t1 << ")\n";
            cout << "110: [" << ximax << "][" << yimax << "][" << zimin << "], (x,z,t) = (" << field[ximax][yimax][zimin].x1 << "," << field[ximax][yimax][zimin].z1 << "," << field[ximax][yimax][zimin].t1 << ")\n";
            cout << "111: [" << ximax << "][" << yimax << "][" << zimax << "], (x,z,t) = (" << field[ximax][yimax][zimax].x1 << "," << field[ximax][yimax][zimax].z1 << "," << field[ximax][yimax][zimax].t1 << ")\n\n";

            if(field[ximin][yimin][zimin].x1 > x1 || field[ximin][yimin][zimax].x1 > x1 ||
               field[ximin][yimax][zimin].x1 > x1 || field[ximin][yimax][zimax].x1 > x1)
                    cerr << "ERROR: Minimal x bound not minimal.\n";
            if(field[ximax][yimin][zimin].x1 < x1 || field[ximax][yimin][zimax].x1 < x1 ||
               field[ximax][yimax][zimin].x1 < x1 || field[ximax][yimax][zimax].x1 < x1)
                    cerr << "ERROR: Maximal x bound not maximal.\n";

            if(field[ximin][yimin][zimin].t1 > t1 || field[ximin][yimin][zimax].t1 > t1 ||
               field[ximax][yimin][zimin].t1 > t1 || field[ximax][yimin][zimax].t1 > t1)
                    cerr << "ERROR: Minimal y bound not minimal.\n";
            if(field[ximin][yimax][zimin].t1 < t1 || field[ximin][yimax][zimax].t1 < t1 ||
               field[ximin][yimax][zimin].t1 < t1 || field[ximin][yimax][zimax].t1 < t1)
                    cerr << "ERROR: Maximal y bound not maximal.\n";

            if(field[ximin][yimin][zimin].z1 > z1 || field[ximax][yimin][zimin].z1 > z1 ||
               field[ximin][yimax][zimin].z1 > z1 || field[ximax][yimax][zimin].z1 > z1)
                    cerr << "ERROR: Minimal z bound not minimal.\n";
            if(field[ximin][yimin][zimax].z1 < z1 || field[ximax][yimin][zimax].z1 < z1 ||
               field[ximin][yimax][zimax].z1 < z1 || field[ximax][yimax][zimax].z1 < z1)
                    cerr << "ERROR: Maximal z bound not maximal.\n";



    // Solve system of 8 linear equations to get the interpolating polynomial

    return {};
}

// estimates difference between two points with given values in map
double Offset(SensorData s, double x1, double z1, double t1)
{
    constexpr double tfact = 0.00327; // time is measured at different scale, it needs weight
    return sqrt(pow(x1-s.x1,2)+pow(z1-s.z1,2)+pow(tfact*(t1-s.t1),2));
}

SensorData RecoPoint(Field<SensorData>* map, double x1, double z1, double t1, double max_err)
{
    // start looking at the same position
    double x = x1;
    double y = (map->ymax+map->ymin)/2;
    double z = z1;
    double step = map->step/10;

    double offset;      // metric of distance between points
    int iterations = 0; // number of iterations should not exceed 100
    double damp = 0.005;  // damping coefficient

    // loop for offset minimization
    do
    {
        // calculate offset gradient
        SensorData xa = map->GetField(x+step,y,z);
        SensorData xb = map->GetField(x-step,y,z);
        SensorData ya = map->GetField(x,y+step,z);
        SensorData yb = map->GetField(x,y-step,z);
        SensorData za = map->GetField(x,y,z+step);
        SensorData zb = map->GetField(x,y,z-step);

        double oxa = Offset(xa,x1,z1,t1);
        double oxb = Offset(xb,x1,z1,t1);
        double oya = Offset(ya,x1,z1,t1);
        double oyb = Offset(yb,x1,z1,t1);
        double oza = Offset(za,x1,z1,t1);
        double ozb = Offset(zb,x1,z1,t1);

        double gradx = (oxa-oxb)/(2*step);
        double grady = (oya-oyb)/(2*step);
        double gradz = (oza-ozb)/(2*step);

        //adjust current guess by minus gradient
        x -= damp*gradx; y -= damp*grady; z -= damp*gradz;

        //check bounds
        if (x < map->xmin) x = map->xmin;
        if (x > map->xmax) x = map->xmax;
        if (y < map->ymin) y = map->ymin;
        if (y > map->ymax) y = map->ymax;
        if (z < map->zmin) z = map->zmin;
        if (z > map->zmax) z = map->zmax;
        //cout << "gradx: " << x << " grady: " << y << " gradz: " << z << "\n";

        // calculate values at current position
        SensorData cur = map->GetField(x,y,z);
        offset = Offset(cur,x1,z1,t1);

        //make sure step isn't too high
        if(offset < 10*step) step /= 10;

        iterations++;
        // cout << "iter: " << iterations << "\n";        
        if (iterations == 10000) cout << "1000 iterations.\n";
    }
    while ((offset > max_err) && (iterations < 1000));

    return {x,y,z,0};
}