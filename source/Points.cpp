// C++ dependencies
#include <iostream>
#include <vector>

// X17 dependencies
#include "Field.h"
#include "Points.h"
#include "X17Utilities.h"

namespace X17
{
    //// Public MicroPoint methods.

    void MicroPoint::MakeTTreeBranches(TTree* tree)
    {
        tree->Branch("x0",&start.point.x);
        tree->Branch("y0",&start.point.y);
        tree->Branch("z0",&start.point.z);
        tree->Branch("t0",&start.t);
        tree->Branch("e0",&e0);
        tree->Branch("x1",&end.point.x);
        tree->Branch("y1",&end.point.y);
        tree->Branch("z1",&end.point.z);
        tree->Branch("t1",&end.t);
        tree->Branch("e1",&e1);
    }

    void MicroPoint::SetTChainBranches(TChain* chain, bool old_data)
    {
        if (old_data)
        {
            std::cout << "MAKE SURE TO CHANGE THE COORDINATE ASSIGNMENT IF USING NEW SIMULATIONS!!!!!!!\n";
            chain->SetBranchAddress("z0",&start.point.x);
            chain->SetBranchAddress("x0",&start.point.y);
            chain->SetBranchAddress("y0",&start.point.z);
            chain->SetBranchAddress("t0",&start.t);
            chain->SetBranchAddress("e0",&e0);
            chain->SetBranchAddress("z1",&end.point.x);
            chain->SetBranchAddress("x1",&end.point.y);
            chain->SetBranchAddress("y1",&end.point.z);
            chain->SetBranchAddress("t1",&end.t);
            chain->SetBranchAddress("e1",&e1);
        }

        else
        {
            chain->SetBranchAddress("x0",&start.point.x);
            chain->SetBranchAddress("y0",&start.point.y);
            chain->SetBranchAddress("z0",&start.point.z);
            chain->SetBranchAddress("t0",&start.t);
            chain->SetBranchAddress("e0",&e0);
            chain->SetBranchAddress("x1",&end.point.x);
            chain->SetBranchAddress("y1",&end.point.y);
            chain->SetBranchAddress("z1",&end.point.z);
            chain->SetBranchAddress("t1",&end.t);
            chain->SetBranchAddress("e1",&e1);
        }
    }





    //// Public MapPoint methods.

    double MapPoint::operator[](int i) const
    {
        switch (i)
        {
        case 0:
            return this->point.point.x;
            break;
        case 1:
            return this->point.point.y;
            break;
        case 2:
            return this->point.point.z;
            break;
        case 3:
            return this->point.t;
            break;
        
        default:
            throw std::out_of_range("MapPoint[] index out of range.\n");
            break;
        }
    }





    //// Functions related to classes from Points.h.

    std::vector<TMarker3DBox*> GetDataMarkers(const std::vector<RecoPoint>& data, double zbin_size)
    {
        std::vector<TMarker3DBox*> markers;
        constexpr double max_size = 0.75;

        // Find maximal count.
        int max_count = 0;
        for (RecoPoint p : data) if (p.count > max_count) max_count = p.count;

        // Create markers.
        for (RecoPoint p : data)
        {
            using namespace constants;

            double rel_size = max_size * p.count / max_count;
            double xlen = rel_size * pad_width  / 2.0;
            double ylen = rel_size * pad_height / 2.0;
            double zlen = rel_size * zbin_size  / 2.0;

            markers.push_back(new TMarker3DBox(p.x(),p.y(),p.z(),xlen,ylen,zlen,0,0));
        }
        
        return markers;
    }
} // namespace X17