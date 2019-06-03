// EnergyPlus, Copyright (c) 1996-2019, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without the U.S. Department of Energy's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// EnergyPlus::SolarShading Unit Tests

// Google Test Headers
#include <gtest/gtest.h>
#include <chrono>

// EnergyPlus Headers
#include <EnergyPlus/DataBSDFWindow.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataGlobals.hh>
#include <EnergyPlus/DataHeatBalance.hh>
#include <EnergyPlus/DataShadowingCombinations.hh>
#include <EnergyPlus/DataSurfaces.hh>
#include <EnergyPlus/DataSystemVariables.hh>
#include <EnergyPlus/DataVectorTypes.hh>
#include <EnergyPlus/HeatBalanceManager.hh>
#include <EnergyPlus/ScheduleManager.hh>
#include <EnergyPlus/SimulationManager.hh>
#include <EnergyPlus/SizingManager.hh>
#include <EnergyPlus/SolarShading.hh>
#include <EnergyPlus/SurfaceGeometry.hh>
#include <EnergyPlus/UtilityRoutines.hh>

#include "Fixtures/EnergyPlusFixture.hh"

using namespace EnergyPlus;
using namespace EnergyPlus::SolarShading;
using namespace ObjexxFCL;

#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

TEST_F(EnergyPlusFixture, PolygonTest_TestTest)
{
    HCX.dimension(100, 100, 0.0);
    HCY.dimension(100, 100, 0.0);
    
    HCA.dimension(100, 16, 0.0);
    HCB.dimension(100, 16, 0.0);
    HCC.dimension(100, 16, 0.0);
    XTEMP.dimension(100, 0.0);
    YTEMP.dimension(100, 0.0);
    ATEMP.dimension(100, 0.0);
    BTEMP.dimension(100, 0.0);
    CTEMP.dimension(100, 0.0);
    XTEMP1.dimension(100, 0.0);
    YTEMP1.dimension(100, 0.0);
    int NS1 = 1;
    int NS2 = 2;
    int NV1 = 4;
    int NV2 = 4;
    int NV3 = 0;
    //populate arrays
    HCX[0] = 600000;
    HCY[0] = 220000;

    HCX[1] = 700000;
    HCY[1] = 220000;

    HCX[2] = 700000;
    HCY[2] = 180000;

    HCX[3] = 600000;
    HCY[3] = 180000;

    //clipping
    HCX[16] = 650000;   
    HCY[16] = 200000;

    HCX[17] = 890000;
    HCY[17] = 200000;

    HCX[18] = 890000;
    HCY[18] = 150000;

    HCX[19] = 650000;
    HCY[19] = 150000;

    //Populate ABC arrays
    
    for (int i = 0; i < 4; i++) {
        HCA[i] = HCY[i] - HCY[(i+1)%4];
        HCB[i] = HCX[(i+1)%4] - HCX[i];
        HCC[i] = HCX[i] * HCY[(i+1)%4] - HCY[i]*HCX[(i+1)%4];
    }

    for (int i = 16; i < 20; i++) {
        HCA[i] = HCY[i] - HCY[16 + ((i+1)%4)];
        HCB[i] = HCX[16 + ((i+1)%4)] - HCX[i];
        HCC[i] = HCX[i] * HCY[16 + ((i+1)%4)] - HCY[i]*HCX[16 + ((i+1)%4)];
    }



    int iterations = 10000;
    int t_i = 100;
    int k = 0;
    for (int i = 0; i < (iterations/4)*4; i++) { //heat up cache
        CLIPPOLY_baseline(NS1, NS2, NV1, NV2, NV3);
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    double baseline_ClipPoly = 0;
    for (int m = 0; m < t_i; m++) {
        start = std::chrono::high_resolution_clock::now(); 
        for (int i = 0; i < (iterations/4)*4; i++) {
            CLIPPOLY_baseline(NS1, NS2, NV1, NV2, NV3);
        }
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        baseline_ClipPoly += duration.count();
    }
    double attempt_ClipPoly = 0;
    for (int m = 0; m < t_i; m++) {
        start = std::chrono::high_resolution_clock::now(); 
        for (int i = 0; i < (iterations/4)*4; i++) {
            CLIPPOLY(NS1, NS2, NV1, NV2, NV3);
        }
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        attempt_ClipPoly += duration.count();
    }
    double ratio_ClipPoly = baseline_ClipPoly/attempt_ClipPoly;
    //EXPECT_TRUE(false) << "\nSpeedup: " << ratio_ClipPoly << "\n";
}
