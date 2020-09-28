#ifndef _GETFILTER_H_
#define _GETFILTER_H_

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TROOT.h"
#include <fstream>
#include <iostream>
#include <TTree.h>
#include <TSystem.h>
#include <sstream>
#include "TVirtualFFT.h"
#include "math.h"
#include "TCanvas.h"

bool GetFilter(int, char**);

#endif
