#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stack>
#include <tuple>
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
    void Normalize();
    double Angle(const Vector& other) const
    {
        return acos((vx*other.vx+vy*other.vy+vz*other.vz)/(this->Magnitude()*other.Magnitude()));
    }
    double SqDist(const Vector& other) const {return pow(vx-other.vx,2)+pow(vy-other.vy,2)+pow(vz-other.vz,2);}
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

Vector operator/(const Vector& v,const double& d)
{
    return Vector{v.vx/d,v.vy/d,v.vz/d};
}

double operator*(const Vector& v1, const Vector& v2)
{
    return v1.vx*v2.vx+v1.vy*v2.vy+v1.vz*v2.vz;
}

void Vector::Normalize() {*this = (*this)/(this->Magnitude());}

double LineSqDist(const Vector& origin, Vector orientation, double max_param, const Vector& point)
{
    orientation.Normalize();
    double t_close = orientation*(point-origin);
    if (t_close > max_param) t_close = max_param;
    if (t_close < 0)      t_close = 0;

    Vector line_vector = origin+t_close*orientation-point;

    return line_vector.SqMagnitude();
}

template <int M, int N>
struct Matrix
{
    double elements[M*N]; // rows before columns

    Matrix() = default;

    Matrix(const double (&arr)[M*N])
    {
        copy(std::begin(arr),std::end(arr),std::begin(elements));
        // cout << "NEW MATRIX:\n"; this->Print();
    }

    void operator+=(const Matrix<M,N>& A)
    {
        for (int r = 0; r < M; r++)
        {
            for (int c = 0; c < N; c++)
            {
                this->at(r,c) += A.at(r,c);
            }        
        }
    }

    void operator*=(const double& d)
    {
        for (int r = 0; r < M; r++)
        {
            for (int c = 0; c < N; c++)
            {
                this->at(r,c) *= d;
            }        
        }
    }
 
    const double& at(int row, int column) const
    {
        if(row > -1 && row < M && column > -1 && column < N) return elements[row*N+column];
        else cerr << "ERROR: Invalid matrix element (" << row << "," << column << ") of " << M << "x" << N << " matrix.\n";
        return elements[0];
    }

    double& at(int row, int column)
    {
        if(row > -1 && row < M && column > -1 && column < N) return elements[row*N+column];
        else cerr << "ERROR: Invalid matrix element (" << row << "," << column << ") of " << M << "x" << N << " matrix.\n";
        return elements[0];
    }

    vector<double> GetColumn(int c)
    {
        if(c > -1 && c < N)
        {
            vector<double> column;
            for (int r = 0; r < M; r++) column.push_back(this->at(r,c));
            return column;
        }
        else cerr << "ERROR: Invalid matrix column index.\n";
        return vector<double>();
    }

    void Print()
    {
        for (int r = 0; r < M; r++)
        {
            for (int c = 0; c < N; c++)
            {
                cout << this->at(r,c) << " ";
            }
            cout << "\n";
        }        
    }

    void SwitchRows(int row1, int row2)
    {
        if (row1 > -1 && row1 < M && row2 > -1 && row2 < M)
        {
            for (int c = 0; c < N; c++)
            {
                double temp = this->at(row1,c);
                this->at(row1,c) = this->at(row2,c);
                this->at(row2,c) = temp;
            }            
        }

        else cerr << "ERROR: Invalid matrix row number.\n";        
    }

    void Reduce()
    {
        for (int c = 0; c < M; c++)
        {
            // Find row with non zero c-th element
            bool not_zero = false;
            for (int r = c; r < M; r++)
            {
                if (this->at(r,c) != 0)
                {
                    not_zero = true;
                    if (r != c) 
                        this->SwitchRows(r,c);
                    break;
                }
            }            
            if (!not_zero) cerr << "WARNING: Singular matrix.\n";

            // Normalize selected row
            double first = this->at(c,c);
            if(first == 0) continue;
            for (int c2 = c; c2 < N; c2++) this->at(c,c2) /= first;            

            // Subtracting rows
            for (int r = 0; r < M; r++)
            {
                if (r != c)
                {
                    double factor = this->at(r,c);
                    for (int c2 = c; c2 < N; c2++)
                    {
                        this->at(r,c2) -= factor*this->at(c,c2);
                    }
                }
            }
            // cout << "REDUCTION OF MATRIX IN PROGRESS:\n"; this->Print();
        }

        // cout << "REDUCED MATRIX:\n"; this->Print();
    }
};

template <int M, int N>
Matrix<M,N> operator+(const Matrix<M,N>& A, const Matrix<M,N>& B)
{
    Matrix<M,N> result;
    for (int r = 0; r < M; r++)
    {
        for (int c = 0; c < N; c++)
        {
            result.at(r,c) = A.at(r,c)+B.at(r,c);
        }        
    }
    
    return result;
}

template <int M, int N>
Matrix<M,N> operator*(const double& d, const Matrix<M,N>& A)
{
    Matrix<M,N> result;
    for (int r = 0; r < M; r++) for (int c = 0; c < N; c++) result.at(r,c) = d*A.at(r,c);
    
    return result;
}

template <int M, int N>
Matrix<M,N> operator/(const Matrix<M,N>& A,const double& d)
{
    Matrix<M,N> result;
    for (int r = 0; r < M; r++) for (int c = 0; c < N; c++) result.at(r,c) = A.at(r,c)/d;
    
    return result;
}

template<int M,int N,int P>
Matrix<M,P> operator*(const Matrix<M,N>& A, const Matrix<N,P>& B)
{
    Matrix<M,P> result;
    for (int r = 0; r < M; r++)
    {
        for (int c = 0; c < P; c++)
        {
            double sum = 0;
            for (int i = 0; i < N; i++) sum += A.at(r,i)*B.at(i,c);
            result.at(r,c) = sum;
        }        
    }
    
    return result;
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
                    // changing coordinate system from magnetic simulation (x,y,z) --> (y,z,x)
                        // x = stod(Z); y = stod(X); z = stod(Y);
                        // vx = stod(VZ); vy = stod(VX); vz = stod(VY);
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

    Vector GetField(const Vector& vec)
    {
        return GetField(vec.vx,vec.vy,vec.vz);
    }
};


struct SensorData
{
    double x1,y1,z1,t1;
    double x1dev,y1dev,z1dev,t1dev;
    double xycorr,xzcorr,xtcorr,yzcorr,ytcorr,ztcorr;

    void operator+=(const SensorData& s)
    {
        this->x1 += s.x1;
        this->y1 += s.y1;
        this->z1 += s.z1;
        this->t1 += s.t1;
    }

    double operator[](int i)
    {
        switch (i)
        {
        case 0:
            return this->x1;
            break;
        case 1:
            return this->y1;
            break;
        case 2:
            return this->z1;
            break;
        case 3:
            return this->t1;
            break;
        
        default:
            cerr << "ERROR: SensorData[] index out of range.\n";
            return -1;
            break;
        }
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

void PrintMap(const Field<SensorData>* map)
{
    for (int i = map->yimax; i > -1; i--)
    {
        cout << "LAYER " << i << "\n";
        cout << "=================================================================\n";

        for (int j = map->zimax; j > -1; j--)
        {
            for (int k = 0; k <= map->ximax; k++)
            {
                cout << map->field[k][i][j].t1 << " ";
            }
            cout << "\n";
        }

        cout << "\n\n";
    }
}

/// @brief Templated function for cube corners selection
/// @tparam T Any number type
/// @param indexes Array containing xmin,xmax,ymin,ymax,zmin,zmax in this order
/// @param corner  Index of the corner (0-7 <----> xyz: 000 --> 111) 
/// @return Tuple with x,y,z of corner in this order
template<typename T>
array<T,3> CubeCorner(const T (&indexes)[6], int corner)
{
    T x,y,z;

    if((int)corner/4 == 0)     x = indexes[0];
    else                       x = indexes[1];

    if(((int)corner/2)%2 == 0) y = indexes[2];
    else                       y = indexes[3];

    if(corner % 2 == 0)        z = indexes[4];
    else                       z = indexes[5];

    return {x,y,z};
}

vector<vector<double>> GetInterpolCoef(Field<SensorData>& map, const int (&indexes)[6])
{
    double bounds[6]; // {xmin, xmax, ymin, ymax, zmin, zmax}

    bounds[0] = indexes[0]*map.step+map.xmin;
    bounds[1] = indexes[1]*map.step+map.xmin;

    bounds[2] = indexes[2]*map.step+map.ymin;
    bounds[3] = indexes[3]*map.step+map.ymin;

    bounds[4] = indexes[4]*map.step+map.zmin;
    bounds[5] = indexes[5]*map.step+map.zmin;

    const int coor_ids[] = {0,1,3}; // x, y and t coordinate indexes in SensorData (operator[])
    vector<vector<double>> xytvalues;
    vector<vector<double>> output;

    for (int id : coor_ids)
    {
        vector<double> coor_values;
        for (int i = 0; i < 8; i++)
        {
            const auto [xi,yi,zi] = CubeCorner(indexes,i);
            coor_values.push_back(map.field[xi][yi][zi][id]);
        }
        xytvalues.push_back(coor_values);
    }

    for (int i = 0; i < 3; i++)
    {
        double arr[8*9];
        for (int j = 0; j < 8; j++)
        {
            arr[9*j]   = 1;
            arr[9*j+1] = xytvalues[0][j];
            arr[9*j+2] = xytvalues[1][j];
            arr[9*j+3] = xytvalues[2][j];
            arr[9*j+4] = xytvalues[0][j]*xytvalues[1][j];
            arr[9*j+5] = xytvalues[0][j]*xytvalues[2][j];
            arr[9*j+6] = xytvalues[1][j]*xytvalues[2][j];
            arr[9*j+7] = xytvalues[0][j]*xytvalues[1][j]*xytvalues[2][j];

            const auto corner = CubeCorner(bounds,j);
            arr[9*j+8] = corner[i];
        }
        
        Matrix<8,9> matrix = Matrix<8,9>(arr);
        matrix.Reduce();
        output.push_back(matrix.GetColumn(8));
    }

    return output;
}  

template<>
SensorData Field<SensorData>::Invert(double x1, double y1, double t1)
{
    // Find 8 closest points using binary search (assuming ordering)
        int ximin = 0;               // Starting minimal search x index
        int ximax = this->ximax;     // Starting maximal search x index
        int yimin = 0;               // Starting minimal search y index
        int yimax = this->yimax;     // Starting maximal search y index
        int zimin = this->zimax;     // Starting minimal search z index (inverted - time is maximal for z minimal)
        int zimax = 0;               // Starting maximal search z index (inverted - time is maximal for z minimal)

        int  ximid = (ximin+ximax)/2; // x midpoint variable
        int  yimid = (yimin+yimax)/2; // y midpoint variable
        int  zimid = (zimin+zimax)/2; // z midpoint variable
        bool ximid_is_max,yimid_is_max,zimid_is_max; // Do midpoint variables contain higher value?

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

        // Some part of the field is initialized with zeros and isn't therefore ordered
        if(ximid != 0)
        {
            while(field[ximid][yimin][zimid].y1 == 0) yimin++;
            while(field[ximid][yimax][zimid].y1 == 0) yimax--;
        }

        // y index search
        while ((yimin+1<yimax))
        {
            yimid = (yimin+yimax)/2;

            double y1low  = field[ximid][yimin][zimid].y1;
            double y1high = field[ximid][yimax][zimid].y1;
            double y1mid  = field[ximid][yimid][zimid].y1;

            if(y1mid < y1low || y1mid > y1high) cerr << "WARNING: Ordering mistake in binary search (y).\n";
            if(y1 <= y1mid) {yimax = yimid; yimid_is_max = true;}
            if(y1 >  y1mid) {yimin = yimid; yimid_is_max = false;}
        }

        // z index search
        while ((zimin>zimax+1))
        {
            zimid = (zimin+zimax)/2;

            double t1low  = field[ximid][yimid][zimin].t1;
            double t1high = field[ximid][yimid][zimax].t1;
            double t1mid  = field[ximid][yimid][zimid].t1;

            if(t1mid < t1low || t1mid > t1high) cerr << "WARNING: Ordering mistake in binary search (z).\n";
            if(t1 <= t1mid) {zimax = zimid; zimid_is_max = true;}
            if(t1 >  t1mid) {zimin = zimid; zimid_is_max = false;}
        }

        // Set up actual min/max variables for the cube
        if(ximid_is_max) {ximax = ximid; ximin = ximax-1;}
        else             {ximin = ximid; ximax = ximin+1;}
        if(yimid_is_max) {yimax = yimid; yimin = yimax-1;}
        else             {yimin = yimid; yimax = yimin+1;}
        if(zimid_is_max) {zimax = zimid; zimin = zimax+1;}
        else             {zimin = zimid; zimax = zimin-1;}

        // Sanity check
            // cout << "Cube for (x1,y1,t1) = (" << x1 << "," << y1 << "," << t1 << "): \n";
            // cout << "000: [" << ximin << "][" << yimin << "][" << zimin << "], (x,y,t) = (" << field[ximin][yimin][zimin].x1 << "," << field[ximin][yimin][zimin].y1 << "," << field[ximin][yimin][zimin].t1 << ")\n";
            // cout << "001: [" << ximax << "][" << yimin << "][" << zimin << "], (x,y,t) = (" << field[ximax][yimin][zimin].x1 << "," << field[ximax][yimin][zimin].y1 << "," << field[ximax][yimin][zimin].t1 << ")\n";
            // cout << "010: [" << ximin << "][" << yimin << "][" << zimax << "], (x,y,t) = (" << field[ximin][yimin][zimax].x1 << "," << field[ximin][yimin][zimax].y1 << "," << field[ximin][yimin][zimax].t1 << ")\n";
            // cout << "011: [" << ximax << "][" << yimin << "][" << zimax << "], (x,y,t) = (" << field[ximax][yimin][zimax].x1 << "," << field[ximax][yimin][zimax].y1 << "," << field[ximax][yimin][zimax].t1 << ")\n";
            // cout << "100: [" << ximin << "][" << yimax << "][" << zimin << "], (x,y,t) = (" << field[ximin][yimax][zimin].x1 << "," << field[ximin][yimax][zimin].y1 << "," << field[ximin][yimax][zimin].t1 << ")\n";
            // cout << "101: [" << ximax << "][" << yimax << "][" << zimin << "], (x,y,t) = (" << field[ximax][yimax][zimin].x1 << "," << field[ximax][yimax][zimin].y1 << "," << field[ximax][yimax][zimin].t1 << ")\n";
            // cout << "110: [" << ximin << "][" << yimax << "][" << zimax << "], (x,y,t) = (" << field[ximin][yimax][zimax].x1 << "," << field[ximin][yimax][zimax].y1 << "," << field[ximin][yimax][zimax].t1 << ")\n";
            // cout << "111: [" << ximax << "][" << yimax << "][" << zimax << "], (x,y,t) = (" << field[ximax][yimax][zimax].x1 << "," << field[ximax][yimax][zimax].y1 << "," << field[ximax][yimax][zimax].t1 << ")\n\n";

            // if(field[ximin][yimin][zimin].x1 > x1 || field[ximin][yimin][zimax].x1 > x1 ||
            //    field[ximin][yimax][zimin].x1 > x1 || field[ximin][yimax][zimax].x1 > x1)
            //         cerr << "ERROR: Minimal x bound not minimal.\n";
            // if(field[ximax][yimin][zimin].x1 < x1 || field[ximax][yimin][zimax].x1 < x1 ||
            //    field[ximax][yimax][zimin].x1 < x1 || field[ximax][yimax][zimax].x1 < x1)
            //         cerr << "ERROR: Maximal x bound not maximal.\n";

            // if(field[ximin][yimin][zimin].y1 > y1 || field[ximin][yimin][zimax].y1 > y1 ||
            //    field[ximax][yimin][zimin].y1 > y1 || field[ximax][yimin][zimax].y1 > y1)
            //         cerr << "ERROR: Minimal y bound not minimal.\n";
            // if(field[ximin][yimax][zimin].y1 < y1 || field[ximin][yimax][zimax].y1 < y1 ||
            //    field[ximin][yimax][zimin].y1 < y1 || field[ximin][yimax][zimax].y1 < y1)
            //         cerr << "ERROR: Maximal y bound not maximal.\n";

            // if(field[ximin][yimin][zimin].t1 > t1 || field[ximax][yimin][zimin].t1 > t1 ||
            //    field[ximin][yimax][zimin].t1 > t1 || field[ximax][yimax][zimin].t1 > t1)
            //         cerr << "ERROR: Minimal z bound not minimal.\n";
            // if(field[ximin][yimin][zimax].t1 < t1 || field[ximax][yimin][zimax].t1 < t1 ||
            //    field[ximin][yimax][zimax].t1 < t1 || field[ximax][yimax][zimax].t1 < t1)
            //         cerr << "ERROR: Maximal z bound not maximal.\n";



    vector<vector<double>> interpol = GetInterpolCoef(*this,{ximin,ximax,yimin,yimax,zimin,zimax});

    double xout = interpol[0][0]+interpol[0][1]*x1+interpol[0][2]*y1+interpol[0][3]*t1+interpol[0][4]*x1*y1+interpol[0][5]*x1*t1+interpol[0][6]*y1*t1+interpol[0][7]*x1*y1*t1;
    double yout = interpol[1][0]+interpol[1][1]*x1+interpol[1][2]*y1+interpol[1][3]*t1+interpol[1][4]*x1*y1+interpol[1][5]*x1*t1+interpol[1][6]*y1*t1+interpol[1][7]*x1*y1*t1;
    double zout = interpol[2][0]+interpol[2][1]*x1+interpol[2][2]*y1+interpol[2][3]*t1+interpol[2][4]*x1*y1+interpol[2][5]*x1*t1+interpol[2][6]*y1*t1+interpol[2][7]*x1*y1*t1;

    return {xout,yout,zout};
}

// estimates difference between two points with given values in map
double Offset(SensorData s, double x1, double y1, double t1)
{
    constexpr double tfact = 0.00327; // time is measured at different scale, it needs weight
    return sqrt(pow(x1-s.x1,2)+pow(y1-s.y1,2)+pow(tfact*(t1-s.t1),2));
}

SensorData RecoPoint(Field<SensorData>* map, double x1, double y1, double t1, double max_err)
{
    // start looking at the same position
    double x = x1;
    double y = y1;
    double z = (map->zmax+map->zmin)/2;
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

        double oxa = Offset(xa,x1,y1,t1);
        double oxb = Offset(xb,x1,y1,t1);
        double oya = Offset(ya,x1,y1,t1);
        double oyb = Offset(yb,x1,y1,t1);
        double oza = Offset(za,x1,y1,t1);
        double ozb = Offset(zb,x1,y1,t1);

        double gradx = (oxa-oxb)/(2*step);
        double grady = (oya-oyb)/(2*step);
        double gradz = (oza-ozb)/(2*step);

        // adjust current guess by minus gradient
        x -= damp*gradx; y -= damp*grady; z -= damp*gradz;

        // check bounds
        if (x < map->xmin) x = map->xmin;
        if (x > map->xmax) x = map->xmax;
        if (y < map->ymin) y = map->ymin;
        if (y > map->ymax) y = map->ymax;
        if (z < map->zmin) z = map->zmin;
        if (z > map->zmax) z = map->zmax;
        //cout << "gradx: " << x << " grady: " << y << " gradz: " << z << "\n";

        // calculate values at current position
        SensorData cur = map->GetField(x,y,z);
        offset = Offset(cur,x1,y1,t1);

        // make sure step isn't too high
        if(offset < 10*step) step /= 10;

        iterations++;
        // cout << "iter: " << iterations << "\n";        
        if (iterations == 10000) cout << "1000 iterations.\n";
    }
    while ((offset > max_err) && (iterations < 1000));

    return {x,y,z,0};
}