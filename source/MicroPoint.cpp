// X17 dependencies
#include "MicroPoint.h"

ClassImp(X17::MicroPoint)

namespace X17
{   
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

    void MicroPoint::SetTTreeBranches(TTree* tree)
    {
        tree->SetBranchAddress("start.point.x",&start.point.x);
        tree->SetBranchAddress("start.point.y",&start.point.y);
        tree->SetBranchAddress("start.point.z",&start.point.z);
        tree->SetBranchAddress("start.t",&start.t);
        tree->SetBranchAddress("end.point.x",&end.point.x);
        tree->SetBranchAddress("end.point.y",&end.point.y);
        tree->SetBranchAddress("end.point.z",&end.point.z);
        tree->SetBranchAddress("end.t",&end.t);
        tree->SetBranchAddress("e0",&e0);
        tree->SetBranchAddress("e1",&e1);
    }
}