#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TApplication.h>
#include <sys/stat.h>
#include <unistd.h> // for sleep()
#include <iostream>

long getFileMTime(const char* filename)
{
    struct stat st;
    if (stat(filename, &st) == 0) return st.st_mtime;
    return -1;
}

void watchAndDisplay(const char* filename, const char* canvasName)
{
    long lastMTime = getFileMTime(filename);
    TCanvas* canvas = nullptr;
    bool firstTime = true;

    while (true)
    {
        long newMTime = getFileMTime(filename);
        if (firstTime)
        {
            firstTime = false;
            newMTime = -1;
        }
        
        if (newMTime != lastMTime)
        {
            lastMTime = newMTime;

            TFile* f = TFile::Open(filename);
            if (!f || f->IsZombie()) 
            {
                std::cerr << "Failed to open file\n";
                delete f;
                lastMTime = -1;
                sleep(1);
                continue;
            }

            canvas = (TCanvas*)f->Get(canvasName);
            if (canvas)
            {
                canvas->Draw();
                canvas->Update();
            } 
            else 
            {
                std::cerr << "Canvas not found in file\n";
                lastMTime = -1;
            }

            f->Close();
            delete f;
        }
        
        gSystem->ProcessEvents(); // keep GUI responsive
        sleep(1); // check every second
        // break;
    }
}

void cview_update()
{
    std::string filename, canvasName;
    std::cout << "Enter filename: ";
    std::getline(std::cin, filename);
    std::cout << "Enter canvas name: ";
    std::getline(std::cin, canvasName);
    watchAndDisplay(filename.c_str(), canvasName.c_str());
}