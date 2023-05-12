// X17 dependencies
#include "Points.h"

namespace X17
{
    void MicroPoint::MakeTTreeBranches(TTree* tree)
    {
        tree->Branch("x0",&x0);
        tree->Branch("y0",&y0);
        tree->Branch("z0",&z0);
        tree->Branch("t0",&t0);
        tree->Branch("e0",&e0);
        tree->Branch("x1",&x1);
        tree->Branch("y1",&y1);
        tree->Branch("z1",&z1);
        tree->Branch("t1",&t1);
        tree->Branch("e1",&e1);
    }

    void MicroPoint::SetTChainBranches(TChain* chain, bool old_data)
    {
        if (old_data)
        {
            std::cout << "MAKE SURE TO CHANGE THE COORDINATE ASSIGNMENT IF USING NEW SIMULATIONS!!!!!!!\n";
            chain->SetBranchAddress("z0",&x0);
            chain->SetBranchAddress("x0",&y0);
            chain->SetBranchAddress("y0",&z0);
            chain->SetBranchAddress("t0",&t0);
            chain->SetBranchAddress("e0",&e0);
            chain->SetBranchAddress("z1",&x1);
            chain->SetBranchAddress("x1",&y1);
            chain->SetBranchAddress("y1",&z1);
            chain->SetBranchAddress("t1",&t1);
            chain->SetBranchAddress("e1",&e1);
        }

        else
        {
            chain->SetBranchAddress("x0",&x0);
            chain->SetBranchAddress("y0",&y0);
            chain->SetBranchAddress("z0",&z0);
            chain->SetBranchAddress("t0",&t0);
            chain->SetBranchAddress("e0",&e0);
            chain->SetBranchAddress("x1",&x1);
            chain->SetBranchAddress("y1",&y1);
            chain->SetBranchAddress("z1",&z1);
            chain->SetBranchAddress("t1",&t1);
            chain->SetBranchAddress("e1",&e1);
        }
    }

    double MapPoint::operator[](int i) const
    {
        switch (i)
        {
        case 0:
            return this->point.x;
            break;
        case 1:
            return this->point.y;
            break;
        case 2:
            return this->point.z;
            break;
        case 3:
            return this->point.t;
            break;
        
        default:
            throw std::out_of_range("MapPoint[] index out of range.\n");
            break;
        }
    }

    std::vector<TMarker3DBox*> GetDataMarkers(std::vector<RecoPoint> data, double zbin_size)
    {
        std::vector<TMarker3DBox*> markers;
        constexpr double max_size = 0.75;

        // find maximal count
        int max_count = 0;
        for (RecoPoint p : data) if (p.count > max_count) max_count = p.count;

        // create markers
        for (RecoPoint p : data)
        {
            using namespace constants;

            double rel_size = max_size * p.count / max_count;
            double xlen = rel_size * pad_width  / 2.0;
            double ylen = rel_size * pad_height / 2.0;
            double zlen = rel_size * zbin_size  / 2.0;

            markers.push_back(new TMarker3DBox(p.x,p.y,p.z,xlen,ylen,zlen,0,0));
        }
        
        return markers;
    }
} // namespace X17