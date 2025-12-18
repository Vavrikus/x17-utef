#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TApplication.h>
#include <sys/stat.h>
#include <sys/select.h>
#include <unistd.h>
#include <iostream>
#include <string>

// --- Helper Functions ---

// Check for console input without stopping the program
bool isInputAvailable()
{
    struct timeval tv;
    fd_set fds;
    tv.tv_sec = 0;
    tv.tv_usec = 0; 
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    select(STDIN_FILENO + 1, &fds, NULL, NULL, &tv);
    return FD_ISSET(STDIN_FILENO, &fds);
}

long getFileMTime(const char* filename)
{
    struct stat st;
    if (stat(filename, &st) == 0) return st.st_mtime;
    return -1;
}

// --- Main Watcher Logic ---

void watchAndDisplay(const char* filename, const char* canvasName)
{
    long lastMTime = getFileMTime(filename);
    TCanvas* displayedCanvas = nullptr; // We will own this copy
    bool firstTime = true;

    std::cout << ">>> Watching '" << canvasName << "'. Type 'stop' + Enter to change canvas.\n";

    while (true)
    {
        // 1. Check for User Command
        if (isInputAvailable()) 
        {
            std::string input;
            std::getline(std::cin, input);
            if (input == "stop") 
            {
                std::cout << "Stopping watch...\n";
                break;
            }
        }

        // 2. Check for File Update
        long newMTime = getFileMTime(filename);
        if (firstTime) { firstTime = false; newMTime = -1; }
        
        if (newMTime != lastMTime)
        {
            lastMTime = newMTime;

            TFile* f = TFile::Open(filename);
            if (!f || f->IsZombie()) 
            {
                // File might be in the middle of writing, wait a bit
                if (f) delete f;
                lastMTime = -1;
                usleep(500000); 
                continue;
            }

            TCanvas* fileCanvas = (TCanvas*)f->Get(canvasName);
            if (fileCanvas)
            {
                // CRITICAL: Delete old canvas to close window/free memory
                if (displayedCanvas) delete displayedCanvas;

                // CRITICAL: Clone the canvas so it survives after file close
                displayedCanvas = (TCanvas*)fileCanvas->Clone();
                displayedCanvas->Draw();
            } 
            else 
            {
                std::cerr << "Canvas '" << canvasName << "' not found.\n";
            }

            f->Close();
            delete f;
        }
        
        gSystem->ProcessEvents(); 
        usleep(100000); // 0.1 seconds
    }

    // Cleanup when changing canvas
    if (displayedCanvas) delete displayedCanvas;
}

// --- Menu Logic ---

void cview_update()
{
    std::string filename, canvasName;

    while (true) // OUTER LOOP: Select File
    {
        std::cout << "\n=== Select File ===\n";
        std::cout << "Enter filename (or 'exit'): ";
        
        if (!std::cin.good()) { std::cin.clear(); std::cin.ignore(1000, '\n'); }
        std::getline(std::cin, filename);

        if (filename == "exit" || filename.empty()) break;

        // Quick check if file exists
        TFile* fTest = TFile::Open(filename.c_str());
        if (!fTest || fTest->IsZombie()) {
            std::cerr << "Error: Cannot open file '" << filename << "'\n";
            if (fTest) delete fTest;
            continue;
        }
        delete fTest;

        while (true) // INNER LOOP: Select Canvas
        {
            // List contents to assist selection
            TFile* f = TFile::Open(filename.c_str());
            std::cout << "\n--- Contents of " << filename << " ---\n";
            f->ls();
            f->Close();
            delete f;

            std::cout << "Enter canvas name (or 'back' to change file, 'exit' to quit): ";
            std::getline(std::cin, canvasName);

            if (canvasName == "back") break; // Break to Outer Loop
            if (canvasName == "exit") return; // Exit Program entirely
            if (canvasName.empty()) continue;

            // Start watching
            watchAndDisplay(filename.c_str(), canvasName.c_str());
        }
    }
    
    std::cout << "Exiting application.\n";
}