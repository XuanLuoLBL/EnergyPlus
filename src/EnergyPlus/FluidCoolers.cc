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

// C++ Headers
#include <cassert>
#include <cmath>

// ObjexxFCL Headers
#include <ObjexxFCL/Array.functions.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/gio.hh>
#include <ObjexxFCL/string.functions.hh>

// EnergyPlus Headers
#include <BranchNodeConnections.hh>
#include <CurveManager.hh>
#include <DataBranchAirLoopPlant.hh>
#include <DataEnvironment.hh>
#include <DataHVACGlobals.hh>
#include <DataIPShortCuts.hh>
#include <DataLoopNode.hh>
#include <DataPlant.hh>
#include <DataPrecisionGlobals.hh>
#include <DataSizing.hh>
#include <FluidCoolers.hh>
#include <FluidProperties.hh>
#include <General.hh>
#include <GlobalNames.hh>
#include <InputProcessing/InputProcessor.hh>
#include <NodeInputManager.hh>
#include <OutAirNodeManager.hh>
#include <OutputProcessor.hh>
#include <OutputReportPredefined.hh>
#include <PlantUtilities.hh>
#include <Psychrometrics.hh>
#include <ReportSizingManager.hh>
#include <ScheduleManager.hh>
#include <UtilityRoutines.hh>

namespace EnergyPlus {

namespace FluidCoolers {

    // Module containing the routines dealing with the objects FluidCooler:SingleSpeed and
    // FluidCooler:TwoSpeed

    // MODULE INFORMATION:
    //       AUTHOR         Chandan Sharma
    //       DATE WRITTEN   August 2008
    //       MODIFIED       April 2010, Chandan Sharma, FSEC
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS MODULE:
    // Model the performance of fluid coolers

    // REFERENCES:
    // Based on cooling tower by Shirey, Raustad: Dec 2000; Shirey, Sept 2002

    // MODULE PARAMETER DEFINITIONS:
    std::string const cFluidCooler_SingleSpeed("FluidCooler:SingleSpeed");
    std::string const cFluidCooler_TwoSpeed("FluidCooler:TwoSpeed");
    static std::string const BlankString;

    bool GetFluidCoolerInputFlag(true);
    int NumSimpleFluidCoolers(0); // Number of similar fluid coolers

    Array1D_bool CheckEquipName;

    // Object Data
    Array1D<FluidCoolerspecs> SimpleFluidCooler; // dimension to number of machines
    std::unordered_map<std::string, std::string> UniqueSimpleFluidCoolerNames;

    void SimFluidCoolers(std::string &FluidCoolerType,
                         std::string &FluidCoolerName,
                         int &CompIndex,
                         bool &RunFlag,
                         bool const InitLoopEquip,
                         Real64 &MaxCap,
                         Real64 &MinCap,
                         Real64 &OptCap)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // Main fluid cooler driver subroutine. Gets called from PlantCondLoopSupplySideManager

        // METHODOLOGY EMPLOYED:
        // After being called by PlantCondLoopSupplySideManager, this subroutine
        // calls GetFluidCoolerInput to get all fluid cooler input info (one time only),
        // then calls the appropriate subroutine to calculate fluid cooler performance,
        // update records (node info) and writes output report info.

        // REFERENCES:
        // Based on SimTowers subroutine by Fred Buhl, May 2002; Richard Raustad, FSEC, Feb 2005

        // LOCAL VARIABLE DECLARATIONS:
        int FluidCoolerNum;

        // GET INPUT
        if (GetFluidCoolerInputFlag) {
            GetFluidCoolerInput();
            GetFluidCoolerInputFlag = false;
        }
        // INITIALIZE
        // Find the correct Equipment
        if (CompIndex == 0) {
            FluidCoolerNum = UtilityRoutines::FindItemInList(FluidCoolerName, SimpleFluidCooler);
            if (FluidCoolerNum == 0) {
                ShowFatalError("SimFluidCoolers: Unit not found = " + FluidCoolerName);
            }
            CompIndex = FluidCoolerNum;
        } else {
            FluidCoolerNum = CompIndex;
            if (FluidCoolerNum > NumSimpleFluidCoolers || FluidCoolerNum < 1) {
                ShowFatalError("SimFluidCoolers:  Invalid CompIndex passed = " + General::TrimSigDigits(FluidCoolerNum) +
                               ", Number of Units = " + General::TrimSigDigits(NumSimpleFluidCoolers) + ", Entered Unit name = " + FluidCoolerName);
            }
            if (CheckEquipName(FluidCoolerNum)) {
                if (FluidCoolerName != SimpleFluidCooler(FluidCoolerNum).Name) {
                    ShowFatalError("SimFluidCoolers: Invalid CompIndex passed = " + General::TrimSigDigits(FluidCoolerNum) + ", Unit name = " +
                                   FluidCoolerName + ", stored Unit Name for that index = " + SimpleFluidCooler(FluidCoolerNum).Name);
                }
                CheckEquipName(FluidCoolerNum) = false;
            }
        }

        // CALCULATE
        {
            auto const SELECT_CASE_var(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num);

            if (SELECT_CASE_var == DataPlant::TypeOf_FluidCooler_SingleSpd) {
                if (InitLoopEquip) {
                    InitFluidCooler(FluidCoolerNum, RunFlag);
                    SizeFluidCooler(FluidCoolerNum);
                    MinCap = 0.0;
                    MaxCap = SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                    OptCap = SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                    return;
                }
                InitFluidCooler(FluidCoolerNum, RunFlag);
                SingleSpeedFluidCooler(FluidCoolerNum);
                UpdateFluidCooler(FluidCoolerNum);
                ReportFluidCooler(RunFlag, FluidCoolerNum);

            } else if (SELECT_CASE_var == DataPlant::TypeOf_FluidCooler_TwoSpd) {
                if (InitLoopEquip) {
                    InitFluidCooler(FluidCoolerNum, RunFlag);
                    SizeFluidCooler(FluidCoolerNum);
                    MinCap = 0.0; // signifies non-load based model (i.e. forward
                    MaxCap = SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                    OptCap = SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                    return;
                }
                InitFluidCooler(FluidCoolerNum, RunFlag);
                TwoSpeedFluidCooler(FluidCoolerNum);
                UpdateFluidCooler(FluidCoolerNum);
                ReportFluidCooler(RunFlag, FluidCoolerNum);

            } else {
                ShowFatalError("SimFluidCoolers: Invalid Fluid Cooler Type Requested = " + FluidCoolerType);
            }
        } // TypeOfEquip
    }

    void GetFluidCoolerInput()
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Chandan Sharma
        //       DATE WRITTEN:    August 2008
        //       MODIFIED         Chandan Sharma, FSEC, April 2010
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        // Obtains input data for fluid coolers and stores it in SimpleFluidCooler data structure.

        // METHODOLOGY EMPLOYED:
        // Uses "Get" routines to read in the data.

        // REFERENCES:
        // Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002;

        // Using/Aliasing
        using namespace DataIPShortCuts; // Data for field names, blank numerics

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int FluidCoolerNum;               // Fluid cooler number, reference counter for SimpleFluidCooler data array
        int NumSingleSpeedFluidCoolers;   // Total number of single-speed Fluid Coolers
        int SingleSpeedFluidCoolerNumber; // Specific single-speed fluid cooler of interest
        int NumTwoSpeedFluidCoolers;      // Number of two-speed Fluid Coolers
        int TwoSpeedFluidCoolerNumber;    // Specific two-speed fluid cooler of interest
        int NumAlphas;                    // Number of elements in the alpha array
        int NumNums;                      // Number of elements in the numeric array
        int IOStat;                       // IO Status when calling get input subroutine
        bool ErrorsFound(false);   // Logical flag set .TRUE. if errors found while getting input data
        Array1D<Real64> NumArray(16);     // Numeric input data array
        Array1D_string AlphArray(5);      // Character string input data array

        //! LKL - still more renaming stuff to go.

        // Get number of all Fluid Coolers specified in the input data file (idf)
        NumSingleSpeedFluidCoolers = inputProcessor->getNumObjectsFound("FluidCooler:SingleSpeed");
        NumTwoSpeedFluidCoolers = inputProcessor->getNumObjectsFound("FluidCooler:TwoSpeed");
        NumSimpleFluidCoolers = NumSingleSpeedFluidCoolers + NumTwoSpeedFluidCoolers;

        if (NumSimpleFluidCoolers <= 0)
            ShowFatalError("No fluid cooler objects found in input, however, a branch object has specified a fluid cooler. Search the input for "
                           "fluid cooler to determine the cause for this error.");

        // See if load distribution manager has already gotten the input
        if (allocated(SimpleFluidCooler)) return;
        GetFluidCoolerInputFlag = false;

        // Allocate data structures to hold fluid cooler input data, report data and fluid cooler inlet conditions
        SimpleFluidCooler.allocate(NumSimpleFluidCoolers);
        UniqueSimpleFluidCoolerNames.reserve(NumSimpleFluidCoolers);
        CheckEquipName.dimension(NumSimpleFluidCoolers, true);

        // Load data structures with fluid cooler input data
        cCurrentModuleObject = cFluidCooler_SingleSpeed;
        for (SingleSpeedFluidCoolerNumber = 1; SingleSpeedFluidCoolerNumber <= NumSingleSpeedFluidCoolers; ++SingleSpeedFluidCoolerNumber) {
            FluidCoolerNum = SingleSpeedFluidCoolerNumber;
            inputProcessor->getObjectItem(cCurrentModuleObject,
                                          SingleSpeedFluidCoolerNumber,
                                          AlphArray,
                                          NumAlphas,
                                          NumArray,
                                          NumNums,
                                          IOStat,
                                          lNumericFieldBlanks,
                                          lAlphaFieldBlanks,
                                          cAlphaFieldNames,
                                          cNumericFieldNames);
            GlobalNames::VerifyUniqueInterObjectName(
                UniqueSimpleFluidCoolerNames, AlphArray(1), cCurrentModuleObject, cAlphaFieldNames(1), ErrorsFound);

            SimpleFluidCooler(FluidCoolerNum).Name = AlphArray(1);
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerType = cCurrentModuleObject;
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num = DataPlant::TypeOf_FluidCooler_SingleSpd;
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerMassFlowRateMultiplier = 2.5;
            SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum = NodeInputManager::GetOnlySingleNode(
                AlphArray(2), ErrorsFound, cCurrentModuleObject, AlphArray(1), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Inlet, 1, DataLoopNode::ObjectIsNotParent);
            SimpleFluidCooler(FluidCoolerNum).WaterOutletNodeNum = NodeInputManager::GetOnlySingleNode(
                AlphArray(3), ErrorsFound, cCurrentModuleObject, AlphArray(1), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Outlet, 1, DataLoopNode::ObjectIsNotParent);
            BranchNodeConnections::TestCompSet(cCurrentModuleObject, AlphArray(1), AlphArray(2), AlphArray(3), "Chilled Water Nodes");
            SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA = NumArray(1);
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity = NumArray(2);
            SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp = NumArray(3);
            SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp = NumArray(4);
            SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp = NumArray(5);
            SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate = NumArray(6);
            if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRateWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate = NumArray(7);
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRateWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower = NumArray(8);
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPowerWasAutoSized = true;
            }

            //   outdoor air inlet node
            if (AlphArray(5).empty()) {
                SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum = 0;
            } else {
                SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum = NodeInputManager::GetOnlySingleNode(AlphArray(5),
                                                                                             ErrorsFound,
                                                                                             cCurrentModuleObject,
                                                                                             SimpleFluidCooler(FluidCoolerNum).Name,
                                                                                             DataLoopNode::NodeType_Air,
                                                                                             DataLoopNode::NodeConnectionType_OutsideAirReference,
                                                                                             1,
                                                                                             DataLoopNode::ObjectIsNotParent);
                if (!OutAirNodeManager::CheckOutAirNodeNumber(SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum)) {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name + "\" " + cAlphaFieldNames(5) + "= \"" +
                                    AlphArray(5) + "\" not valid.");
                    ShowContinueError("...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.");
                    ErrorsFound = true;
                }
            }

            ErrorsFound = ErrorsFound || TestFluidCoolerSingleSpeedInputForDesign(
                                             cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames, FluidCoolerNum);

        } // End Single-Speed fluid cooler Loop

        cCurrentModuleObject = cFluidCooler_TwoSpeed;
        for (TwoSpeedFluidCoolerNumber = 1; TwoSpeedFluidCoolerNumber <= NumTwoSpeedFluidCoolers; ++TwoSpeedFluidCoolerNumber) {
            FluidCoolerNum = NumSingleSpeedFluidCoolers + TwoSpeedFluidCoolerNumber;
            inputProcessor->getObjectItem(cCurrentModuleObject,
                                          TwoSpeedFluidCoolerNumber,
                                          AlphArray,
                                          NumAlphas,
                                          NumArray,
                                          NumNums,
                                          IOStat,
                                          lNumericFieldBlanks,
                                          lAlphaFieldBlanks,
                                          cAlphaFieldNames,
                                          cNumericFieldNames);
            GlobalNames::VerifyUniqueInterObjectName(
                UniqueSimpleFluidCoolerNames, AlphArray(1), cCurrentModuleObject, cAlphaFieldNames(1), ErrorsFound);

            SimpleFluidCooler(FluidCoolerNum).Name = AlphArray(1);
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerType = cCurrentModuleObject;
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num = DataPlant::TypeOf_FluidCooler_TwoSpd;
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerMassFlowRateMultiplier = 2.5;
            SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum = NodeInputManager::GetOnlySingleNode(
                AlphArray(2), ErrorsFound, cCurrentModuleObject, AlphArray(1), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Inlet, 1, DataLoopNode::ObjectIsNotParent);
            SimpleFluidCooler(FluidCoolerNum).WaterOutletNodeNum = NodeInputManager::GetOnlySingleNode(
                AlphArray(3), ErrorsFound, cCurrentModuleObject, AlphArray(1), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Outlet, 1, DataLoopNode::ObjectIsNotParent);
            BranchNodeConnections::TestCompSet(cCurrentModuleObject, AlphArray(1), AlphArray(2), AlphArray(3), "Chilled Water Nodes");

            SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA = NumArray(1);
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUAWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA = NumArray(2);
            if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUAWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUASizingFactor = NumArray(3);
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity = NumArray(4);
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap = NumArray(5);
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCapWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCapSizingFactor = NumArray(6);
            SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp = NumArray(7);
            SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp = NumArray(8);
            SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp = NumArray(9);
            SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate = NumArray(10);
            if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRateWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate = NumArray(11);
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRateWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower = NumArray(12);
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPowerWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate = NumArray(13);
            if (SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRateWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRateSizingFactor = NumArray(14);
            SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower = NumArray(15);
            if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower == DataSizing::AutoSize) {
                SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPowerWasAutoSized = true;
            }
            SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPowerSizingFactor = NumArray(16);

            //   outdoor air inlet node
            if (AlphArray(5).empty()) {
                SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum = 0;
            } else {
                SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum = NodeInputManager::GetOnlySingleNode(AlphArray(5),
                                                                                             ErrorsFound,
                                                                                             cCurrentModuleObject,
                                                                                             SimpleFluidCooler(FluidCoolerNum).Name,
                                                                                             DataLoopNode::NodeType_Air,
                                                                                             DataLoopNode::NodeConnectionType_OutsideAirReference,
                                                                                             1,
                                                                                             DataLoopNode::ObjectIsNotParent);
                if (!OutAirNodeManager::CheckOutAirNodeNumber(SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum)) {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name + "\" " + cAlphaFieldNames(5) + "= \"" +
                                    AlphArray(5) + "\" not valid.");
                    ShowContinueError("...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.");
                    ErrorsFound = true;
                }
            }

            ErrorsFound = ErrorsFound || TestFluidCoolerTwoSpeedInputForDesign(
                                             cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames, FluidCoolerNum);
        }

        if (ErrorsFound) {
            ShowFatalError("Errors found in getting fluid cooler input.");
        }

        // Set up output variables, CurrentModuleObject='FluidCooler:SingleSpeed'
        for (FluidCoolerNum = 1; FluidCoolerNum <= NumSingleSpeedFluidCoolers; ++FluidCoolerNum) {
            SetupOutputVariable("Cooling Tower Inlet Temperature",
                                OutputProcessor::Unit::C,
                                SimpleFluidCooler(FluidCoolerNum).InletWaterTemp,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Outlet Temperature",
                                OutputProcessor::Unit::C,
                                SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Mass Flow Rate",
                                OutputProcessor::Unit::kg_s,
                                SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Heat Transfer Rate",
                                OutputProcessor::Unit::W,
                                SimpleFluidCooler(FluidCoolerNum).Qactual,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Fan Electric Power",
                                OutputProcessor::Unit::W,
                                SimpleFluidCooler(FluidCoolerNum).FanPower,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Fan Electric Energy",
                                OutputProcessor::Unit::J,
                                SimpleFluidCooler(FluidCoolerNum).FanEnergy,
                                "System",
                                "Sum",
                                SimpleFluidCooler(FluidCoolerNum).Name,
                                _,
                                "Electric",
                                "HeatRejection",
                                _,
                                "Plant");
        }

        // CurrentModuleObject='FluidCooler:TwoSpeed'
        for (FluidCoolerNum = NumSingleSpeedFluidCoolers + 1; FluidCoolerNum <= NumSingleSpeedFluidCoolers + NumTwoSpeedFluidCoolers;
             ++FluidCoolerNum) {
            SetupOutputVariable("Cooling Tower Inlet Temperature",
                                OutputProcessor::Unit::C,
                                SimpleFluidCooler(FluidCoolerNum).InletWaterTemp,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Outlet Temperature",
                                OutputProcessor::Unit::C,
                                SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Mass Flow Rate",
                                OutputProcessor::Unit::kg_s,
                                SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Heat Transfer Rate",
                                OutputProcessor::Unit::W,
                                SimpleFluidCooler(FluidCoolerNum).Qactual,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Fan Electric Power",
                                OutputProcessor::Unit::W,
                                SimpleFluidCooler(FluidCoolerNum).FanPower,
                                "System",
                                "Average",
                                SimpleFluidCooler(FluidCoolerNum).Name);
            SetupOutputVariable("Cooling Tower Fan Electric Energy",
                                OutputProcessor::Unit::J,
                                SimpleFluidCooler(FluidCoolerNum).FanEnergy,
                                "System",
                                "Sum",
                                SimpleFluidCooler(FluidCoolerNum).Name,
                                _,
                                "Electric",
                                "HeatRejection",
                                _,
                                "Plant");
        }
    }

    bool TestFluidCoolerSingleSpeedInputForDesign(std::string const &cCurrentModuleObject,
                                                  Array1D<std::string> const &AlphArray,
                                                  Array1D<std::string> const &cNumericFieldNames,
                                                  Array1D<std::string> const &cAlphaFieldNames,
                                                  int const &FluidCoolerNum)
    {
        // FUNCTION INFORMATION:
        //       AUTHOR:          Chandan Sharma
        //       DATE WRITTEN:    August 2008
        //       MODIFIED         Chandan Sharma, FSEC, April 2010
        //       RE-ENGINEERED    Jason Glazer, GARD Analytics, February 2015, refactor into a separate function

        // PURPOSE OF THIS FUNCTION:
        // Separate the testing of inputs related to design so that it could be called from the unit tests

        // REFERENCES:
        // Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002;

        // FUNCTION LOCAL VARIABLE DECLARATIONS:
        bool ErrorsFound = false;

        //   Design entering water temperature, design entering air temperature and design entering air
        //   wetbulb temperature must be specified for the both the performance input methods
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp <= 0.0) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(3) +
                            "\", entered value <= 0.0, but must be > 0 ");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp <= 0.0) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(4) +
                            "\", entered value <= 0.0, but must be > 0 ");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp <= 0.0) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(5) +
                            "\", entered value <= 0.0, but must be > 0 ");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp) {
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\"," + cNumericFieldNames(3) + " must be greater than " +
                            cNumericFieldNames(4) + '.');
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp) {
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\"," + cNumericFieldNames(4) + " must be greater than " +
                            cNumericFieldNames(5) + '.');
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate <= 0.0 && SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate != DataSizing::AutoSize) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(7) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate <= 0.0 && !SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRateWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(6) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower <= 0.0 && SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower != DataSizing::AutoSize) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(8) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }

        //   Check various inputs for both the performance input methods
        if (UtilityRoutines::SameString(AlphArray(4), "UFactorTimesAreaAndDesignWaterFlowRate")) {
            SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num = PerfInputMethod::U_FACTOR;
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA <= 0.0 &&
                SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA != DataSizing::AutoSize) {
                ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(1) +
                                "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
                ErrorsFound = true;
            }
        } else if (UtilityRoutines::SameString(AlphArray(4), "NominalCapacity")) {
            SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num = PerfInputMethod::NOMINAL_CAPACITY;
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity <= 0.0) {
                ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(2) +
                                "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA != 0.0) {
                if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA > 0.0) {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Nominal fluid cooler capacity and design fluid cooler UA have been specified.");
                } else {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Nominal fluid cooler capacity has been specified and design fluid cooler UA is being autosized.");
                }
                ShowContinueError(
                    "Design fluid cooler UA field must be left blank when nominal fluid cooler capacity performance input method is used.");
                ErrorsFound = true;
            }
        } else { // Fluid cooler performance input method is not specified as a valid "choice"
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\", invalid " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
            ShowContinueError(R"(... must be "UFactorTimesAreaAndDesignWaterFlowRate" or "NominalCapacity".)");
            ErrorsFound = true;
        }
        return ErrorsFound;
    }

    bool TestFluidCoolerTwoSpeedInputForDesign(std::string const &cCurrentModuleObject,
                                               Array1D<std::string> const &AlphArray,
                                               Array1D<std::string> const &cNumericFieldNames,
                                               Array1D<std::string> const &cAlphaFieldNames,
                                               int const &FluidCoolerNum)
    {
        // FUNCTION INFORMATION:
        //       AUTHOR:          Chandan Sharma
        //       DATE WRITTEN:    August 2008
        //       MODIFIED         Chandan Sharma, FSEC, April 2010
        //       RE-ENGINEERED    Jason Glazer, GARD Analytics, February 2015, refactor into a separate function

        // PURPOSE OF THIS FUNCTION:
        // Separate the testing of inputs related to design so that it could be called from the unit tests

        // REFERENCES:
        // Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002;

        // FUNCTION LOCAL VARIABLE DECLARATIONS:
        bool ErrorsFound = false;

        //   Design entering water temperature, design entering air temperature and design entering air
        //   wetbulb temperature must be specified for the both the performance input methods
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp <= 0.0) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(7) +
                            "\", entered value <= 0.0, but must be > 0 ");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp <= 0.0) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(8) +
                            "\", entered value <= 0.0, but must be > 0 ");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp <= 0.0) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(9) +
                            "\", entered value <= 0.0, but must be > 0 ");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", " + cNumericFieldNames(7) + " must be greater than " +
                            cNumericFieldNames(8) + '.');
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", " + cNumericFieldNames(8) + " must be greater than " +
                            cNumericFieldNames(9) + '.');
            ErrorsFound = true;
        }

        //   Check various inputs for both the performance input methods
        if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate <= 0.0 && !SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRateWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(10) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + "= \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate <= 0.0 && !SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRateWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(11) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + "= \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate <= 0.0 && !SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRateWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(13) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + "= \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        //   High speed air flow rate must be greater than low speed air flow rate.
        //   Can't tell yet if autosized, check later in InitFluidCooler.
        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate <= SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate &&
            !SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRateWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                            "\". Fluid cooler air flow rate at low fan speed must be less than the air flow rate at high fan speed.");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower <= 0.0 && !SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPowerWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(12) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower <= 0.0 && !SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPowerWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(15) +
                            "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
            ErrorsFound = true;
        }
        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower <= SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower &&
            !SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPowerWasAutoSized) {
            ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                            "\". Fluid cooler low speed fan power must be less than high speed fan power.");
            ErrorsFound = true;
        }

        if (UtilityRoutines::SameString(AlphArray(4), "UFactorTimesAreaAndDesignWaterFlowRate")) {
            SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num = PerfInputMethod::U_FACTOR;
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA <= 0.0 &&
                !SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUAWasAutoSized) {
                ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(1) +
                                "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA <= 0.0 &&
                !SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUAWasAutoSized) {
                ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(2) +
                                "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + " = \"" + AlphArray(4) + "\".");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA <= SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA &&
                !SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUAWasAutoSized) {
                ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                "\". Fluid cooler UA at low fan speed must be less than the fluid cooler UA at high fan speed.");
                ErrorsFound = true;
            }
        } else if (UtilityRoutines::SameString(AlphArray(4), "NominalCapacity")) {
            SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num = PerfInputMethod::NOMINAL_CAPACITY;
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity <= 0.0) {
                ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(4) +
                                "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + "= \"" + AlphArray(4) + "\".");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap <= 0.0 &&
                !SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCapWasAutoSized) {
                ShowSevereError(cCurrentModuleObject + " = \"" + AlphArray(1) + "\", invalid data for \"" + cNumericFieldNames(5) +
                                "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames(4) + "= \"" + AlphArray(4) + "\".");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA != 0.0) {
                if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA > 0.0) {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Nominal capacity input method and fluid cooler UA at high fan speed have been specified.");
                } else {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Nominal capacity input method has been specified and fluid cooler UA at high fan speed is being autosized.");
                }
                ShowContinueError(
                    "Fluid cooler UA at high fan speed must be left blank when nominal fluid cooler capacity performance input method is used.");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA != 0.0) {
                if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA > 0.0) {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Nominal capacity input method and fluid cooler UA at low fan speed have been specified.");
                } else {
                    ShowSevereError(cCurrentModuleObject + "= \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Nominal capacity input method has been specified and fluid cooler UA at low fan speed is being autosized.");
                }
                ShowContinueError(
                    "Fluid cooler UA at low fan speed must be left blank when nominal fluid cooler capacity performance input method is used.");
                ErrorsFound = true;
            }
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap >= SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity) {
                ShowSevereError(cCurrentModuleObject + " = \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                "\". Low-speed nominal capacity must be less than the high-speed nominal capacity.");
                ErrorsFound = true;
            }
        } else { // Fluid cooler performance input method is not specified as a valid "choice"
            ShowSevereError(cCurrentModuleObject + "= \"" + AlphArray(1) + "\", invalid " + cAlphaFieldNames(4) + "= \"" + AlphArray(4) + "\".");
            ShowContinueError(R"(... must be "UFactorTimesAreaAndDesignWaterFlowRate" or "NominalCapacity".)");
            ErrorsFound = true;
        }
        return ErrorsFound;
    }

    void InitFluidCooler(int const FluidCoolerNum,     // Number of the current fluid cooler being simulated
                         bool const EP_UNUSED(RunFlag) // TRUE if fluid cooler is ON
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for initializations of the fluid cooler components and for
        // final checking of fluid cooler inputs (post autosizing)

        // METHODOLOGY EMPLOYED:
        // Uses the status flags to trigger initializations.

        // REFERENCES:
        // Based on InitTower subroutine by Don Shirey Sept/Oct 2002, F Buhl Oct 2002

        // SUBROUTINE PARAMETER DEFINITIONS:
        static std::string const RoutineName("InitFluidCooler");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        bool ErrorsFound(false); // Flag if input data errors are found
        int TypeOf_Num(0);
        Real64 rho; // local density of fluid

        if (SimpleFluidCooler(FluidCoolerNum).oneTimeInit) {

            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_SingleSpd) {
                TypeOf_Num = DataPlant::TypeOf_FluidCooler_SingleSpd;
            } else if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_TwoSpd) {
                TypeOf_Num = DataPlant::TypeOf_FluidCooler_TwoSpd;
            } else {
                assert(false);
            }

            // Locate the tower on the plant loops for later usage
            PlantUtilities::ScanPlantLoopsForObject(SimpleFluidCooler(FluidCoolerNum).Name,
                                    TypeOf_Num,
                                    SimpleFluidCooler(FluidCoolerNum).LoopNum,
                                    SimpleFluidCooler(FluidCoolerNum).LoopSideNum,
                                    SimpleFluidCooler(FluidCoolerNum).BranchNum,
                                    SimpleFluidCooler(FluidCoolerNum).CompNum,
                                    ErrorsFound,
                                    _,
                                    _,
                                    _,
                                    _,
                                    _);

            if (ErrorsFound) {
                ShowFatalError("InitFluidCooler: Program terminated due to previous condition(s).");
            }

            SimpleFluidCooler(FluidCoolerNum).oneTimeInit = false;
        }

        // Begin environment initializations
        if (SimpleFluidCooler(FluidCoolerNum).beginEnvrnInit && DataGlobals::BeginEnvrnFlag && (DataPlant::PlantFirstSizesOkayToFinalize)) {

            rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                   DataGlobals::InitConvTemp,
                                                    DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                   RoutineName);
            SimpleFluidCooler(FluidCoolerNum).DesWaterMassFlowRate = SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate * rho;
            PlantUtilities::InitComponentNodes(0.0,
                               SimpleFluidCooler(FluidCoolerNum).DesWaterMassFlowRate,
                               SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum,
                               SimpleFluidCooler(FluidCoolerNum).WaterOutletNodeNum,
                               SimpleFluidCooler(FluidCoolerNum).LoopNum,
                               SimpleFluidCooler(FluidCoolerNum).LoopSideNum,
                               SimpleFluidCooler(FluidCoolerNum).BranchNum,
                               SimpleFluidCooler(FluidCoolerNum).CompNum);
            SimpleFluidCooler(FluidCoolerNum).beginEnvrnInit = false;
        }

        if (!DataGlobals::BeginEnvrnFlag) {
            SimpleFluidCooler(FluidCoolerNum).beginEnvrnInit = true;
        }

        // Each time initializations
        SimpleFluidCooler(FluidCoolerNum).WaterTemp = DataLoopNode::Node(SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum).Temp;

        if (SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum != 0) {
            SimpleFluidCooler(FluidCoolerNum).AirTemp = DataLoopNode::Node(SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum).Temp;
            SimpleFluidCooler(FluidCoolerNum).AirHumRat = DataLoopNode::Node(SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum).HumRat;
            SimpleFluidCooler(FluidCoolerNum).AirPress = DataLoopNode::Node(SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum).Press;
            SimpleFluidCooler(FluidCoolerNum).AirWetBulb = DataLoopNode::Node(SimpleFluidCooler(FluidCoolerNum).OutdoorAirInletNodeNum).OutAirWetBulb;
        } else {
            SimpleFluidCooler(FluidCoolerNum).AirTemp = DataEnvironment::OutDryBulbTemp;
            SimpleFluidCooler(FluidCoolerNum).AirHumRat = DataEnvironment::OutHumRat;
            SimpleFluidCooler(FluidCoolerNum).AirPress = DataEnvironment::OutBaroPress;
            SimpleFluidCooler(FluidCoolerNum).AirWetBulb = DataEnvironment::OutWetBulbTemp;
        }

        SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate = PlantUtilities::RegulateCondenserCompFlowReqOp(SimpleFluidCooler(FluidCoolerNum).LoopNum,
                                                           SimpleFluidCooler(FluidCoolerNum).LoopSideNum,
                                                           SimpleFluidCooler(FluidCoolerNum).BranchNum,
                                                           SimpleFluidCooler(FluidCoolerNum).CompNum,
                                                           SimpleFluidCooler(FluidCoolerNum).DesWaterMassFlowRate *
                                                               SimpleFluidCooler(FluidCoolerNum).FluidCoolerMassFlowRateMultiplier);

        PlantUtilities::SetComponentFlowRate(SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate,
                             SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum,
                             SimpleFluidCooler(FluidCoolerNum).WaterOutletNodeNum,
                             SimpleFluidCooler(FluidCoolerNum).LoopNum,
                             SimpleFluidCooler(FluidCoolerNum).LoopSideNum,
                             SimpleFluidCooler(FluidCoolerNum).BranchNum,
                             SimpleFluidCooler(FluidCoolerNum).CompNum);
    }

    void SizeFluidCooler(int const FluidCoolerNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       April 2010, Chandan Sharma, FSEC
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for sizing fluid cooler Components for which capacities and flow rates
        // have not been specified in the input. This subroutine also calculates fluid cooler UA if the user
        // has specified fluid cooler performance via the "Nominal Capacity" method.

        // METHODOLOGY EMPLOYED:
        // Obtains condenser flow rate from the plant sizing array. If fluid cooler performance is specified
        // via the "Nominal Capacity" method, the water flow rate is directly proportional to capacity.

        // REFERENCES:
        // Based on SizeTower by Don Shirey, Sept/Oct 2002; Richard Raustad, Feb 2005

        // SUBROUTINE PARAMETER DEFINITIONS:
        int const MaxIte(500);    // Maximum number of iterations
        Real64 const Acc(0.0001); // Accuracy of result
        static std::string const CalledFrom("SizeFluidCooler");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int PltSizCondNum(0);           // Plant Sizing index for condenser loop
        int SolFla;                     // Flag of solver
        Real64 DesFluidCoolerLoad(0.0); // Design fluid cooler load [W]
        Real64 UA0;                     // Lower bound for UA [W/C]
        Real64 UA1;                     // Upper bound for UA [W/C]
        Real64 UA;                      // Calculated UA value
        Real64 OutWaterTempAtUA0;       // Water outlet temperature at UA0
        Real64 OutWaterTempAtUA1;       // Water outlet temperature at UA1
        Array1D<Real64> Par(5);         // Parameter array need for RegulaFalsi routine
        std::string equipName;
        Real64 Cp;                            // local specific heat for fluid
        Real64 rho;                           // local density for fluid
        Real64 tmpDesignWaterFlowRate;        // local temporary for water volume flow rate
        Real64 tmpHighSpeedFanPower;          // local temporary for high speed fan power
        Real64 tmpHighSpeedAirFlowRate;       // local temporary for high speed air flow rate
        Real64 tmpHighSpeedEvapFluidCoolerUA; // local temporary for high speed cooler UA
        bool ErrorsFound;

        tmpDesignWaterFlowRate = SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate;
        tmpHighSpeedAirFlowRate = SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate;
        // Find the appropriate Plant Sizing object
        PltSizCondNum = DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).PlantSizNum;

        if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRateWasAutoSized) {
            if (PltSizCondNum > 0) {
                if (DataSizing::PlantSizData(PltSizCondNum).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow) {
                    tmpDesignWaterFlowRate = DataSizing::PlantSizData(PltSizCondNum).DesVolFlowRate;
                    if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate = tmpDesignWaterFlowRate;
                } else {
                    tmpDesignWaterFlowRate = 0.0;
                    if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate = tmpDesignWaterFlowRate;
                }
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Design Water Flow Rate [m3/s]",
                                           SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Design Water Flow Rate [m3/s]",
                                           SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate);
                    }
                }
            } else {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    ShowSevereError("Autosizing error for fluid cooler object = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    ShowFatalError("Autosizing of fluid cooler condenser flow rate requires a loop Sizing:Plant object.");
                }
            }
            // This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
            // temperature is less than design inlet air dry bulb temperature
            if (DataSizing::PlantSizData(PltSizCondNum).ExitTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp && DataPlant::PlantFirstSizesOkayToFinalize) {
                ShowSevereError("Error when autosizing the UA value for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name + '.');
                ShowContinueError("Design Loop Exit Temperature (" + General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2) +
                                  " C) must be greater than design entering air dry-bulb temperature (" +
                                          General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp, 2) +
                                  " C) when autosizing the fluid cooler UA.");
                ShowContinueError("It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the Fluid Cooler "
                                  "design approach temperature (e.g., 4 C).");
                ShowContinueError("If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be "
                                  "> design inlet air dry-bulb temp if autosizing the Fluid Cooler.");
                ShowFatalError("Review and revise design input values as appropriate.");
            }
        }

        PlantUtilities::RegisterPlantCompDesignFlow(SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum, tmpDesignWaterFlowRate);

        if (SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num == PerfInputMethod::U_FACTOR &&
            SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUAWasAutoSized) {
            if (PltSizCondNum > 0) {
                rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                       DataGlobals::InitConvTemp,
                                                        DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                       CalledFrom);
                Cp = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                                            DataSizing::PlantSizData(PltSizCondNum).ExitTemp,
                                                            DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                           CalledFrom);
                DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * DataSizing::PlantSizData(PltSizCondNum).DeltaT;
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity = DesFluidCoolerLoad;
            } else {
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity = 0.0;
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPowerWasAutoSized) {
            // We assume the nominal fan power is 0.0105 times the design load
            if (SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num == PerfInputMethod::NOMINAL_CAPACITY) {
                tmpHighSpeedFanPower = 0.0105 * SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower = tmpHighSpeedFanPower;
            } else {
                if (DesFluidCoolerLoad > 0.0) {
                    tmpHighSpeedFanPower = 0.0105 * DesFluidCoolerLoad;
                    if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower = tmpHighSpeedFanPower;
                } else if (PltSizCondNum > 0) {
                    if (DataSizing::PlantSizData(PltSizCondNum).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow) {
                        // This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
                        // temperature is less than design inlet air dry bulb temperature
                        if (DataSizing::PlantSizData(PltSizCondNum).ExitTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp &&
                                DataPlant::PlantFirstSizesOkayToFinalize) {
                            ShowSevereError("Error when autosizing the UA value for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name + '.');
                            ShowContinueError("Design Loop Exit Temperature (" + General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2) +
                                              " C) must be greater than design entering air dry-bulb temperature (" +
                                                      General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp, 2) +
                                              " C) when autosizing the fluid cooler UA.");
                            ShowContinueError("It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the "
                                              "Fluid Cooler design approach temperature (e.g., 4 C).");
                            ShowContinueError("If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design "
                                              "Setpoint must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler.");
                            ShowFatalError("Review and revise design input values as appropriate.");
                        }
                        rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                               DataGlobals::InitConvTemp,
                                                                DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                               CalledFrom);
                        Cp = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                                                    DataSizing::PlantSizData(PltSizCondNum).ExitTemp,
                                                                    DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                                   CalledFrom);
                        DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * DataSizing::PlantSizData(PltSizCondNum).DeltaT;
                        tmpHighSpeedFanPower = 0.0105 * DesFluidCoolerLoad;
                        if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower = tmpHighSpeedFanPower;
                    } else {
                        tmpHighSpeedFanPower = 0.0;
                        if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower = tmpHighSpeedFanPower;
                    }
                } else {
                    if (DataPlant::PlantFirstSizesOkayToFinalize) {
                        ShowSevereError("Autosizing of fluid cooler fan power requires a loop Sizing:Plant object.");
                        ShowFatalError(" Occurs in fluid cooler object = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    }
                }
            }
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_SingleSpd) {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Fan Power at Design Air Flow Rate [W]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Fan Power at Design Air Flow Rate [W]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower);
                    }
                }
            } else if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_TwoSpd) {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Fan Power at High Fan Speed [W]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Fan Power at High Fan Speed [W]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower);
                    }
                }
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRateWasAutoSized) {
            if (SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num == PerfInputMethod::NOMINAL_CAPACITY) {
                tmpHighSpeedAirFlowRate =
                    SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity /
                    (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp - SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp) * 4.0;
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
            } else {
                if (DesFluidCoolerLoad > 0.0) {
                    tmpHighSpeedAirFlowRate =
                        DesFluidCoolerLoad /
                        (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp - SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp) * 4.0;
                    if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
                } else if (PltSizCondNum > 0) {
                    if (DataSizing::PlantSizData(PltSizCondNum).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow) {
                        // This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
                        // temperature is less than design inlet air dry bulb temperature
                        if (DataSizing::PlantSizData(PltSizCondNum).ExitTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp &&
                                DataPlant::PlantFirstSizesOkayToFinalize) {
                            ShowSevereError("Error when autosizing the UA value for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name + '.');
                            ShowContinueError("Design Loop Exit Temperature (" + General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2) +
                                              " C) must be greater than design entering air dry-bulb temperature (" +
                                                      General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp, 2) +
                                              " C) when autosizing the fluid cooler UA.");
                            ShowContinueError("It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the "
                                              "Fluid Cooler design approach temperature (e.g., 4 C).");
                            ShowContinueError("If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design "
                                              "Setpoint must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler.");
                            ShowFatalError("Review and revise design input values as appropriate.");
                        }
                        rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                               DataGlobals::InitConvTemp,
                                                                DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                               CalledFrom);
                        Cp = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                                                    DataSizing::PlantSizData(PltSizCondNum).ExitTemp,
                                                                    DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                                   CalledFrom);
                        DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * DataSizing::PlantSizData(PltSizCondNum).DeltaT;
                        tmpHighSpeedAirFlowRate =
                            DesFluidCoolerLoad /
                            (SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp - SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp) *
                            4.0;
                        if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
                    } else {
                        tmpHighSpeedAirFlowRate = 0.0;
                        if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
                    }
                } else {
                    if (DataPlant::PlantFirstSizesOkayToFinalize) {
                        ShowSevereError("Autosizing of fluid cooler air flow rate requires a loop Sizing:Plant object");
                        ShowFatalError(" Occurs in fluid cooler object = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    }
                }
            }
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_SingleSpd) {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Design Air Flow Rate [m3/s]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Design Air Flow Rate [m3/s]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate);
                    }
                }
            } else if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType == "FluidCooler:TwoSpeed") {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Air Flow Rate at High Fan Speed [m3/s]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Air Flow Rate at High Fan Speed [m3/s]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate);
                    }
                }
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUAWasAutoSized && DataPlant::PlantFirstSizesOkayToFinalize) {
            if (PltSizCondNum > 0) {
                if (DataSizing::PlantSizData(PltSizCondNum).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow) {
                    // This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
                    // temperature is less than design inlet air dry bulb temperature
                    if (DataSizing::PlantSizData(PltSizCondNum).ExitTemp <= SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp &&
                            DataPlant::PlantFirstSizesOkayToFinalize) {
                        ShowSevereError("Error when autosizing the UA value for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name + '.');
                        ShowContinueError("Design Loop Exit Temperature (" + General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2) +
                                          " C) must be greater than design entering air dry-bulb temperature (" +
                                                  General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp, 2) +
                                          " C) when autosizing the fluid cooler UA.");
                        ShowContinueError("It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the Fluid "
                                          "Cooler design approach temperature (e.g., 4 C).");
                        ShowContinueError("If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint "
                                          "must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler.");
                        ShowFatalError("Review and revise design input values as appropriate.");
                    }
                    rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                           DataGlobals::InitConvTemp,
                                                            DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                           CalledFrom);
                    Cp = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                                                DataSizing::PlantSizData(PltSizCondNum).ExitTemp,
                                                                DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                               CalledFrom);
                    DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * DataSizing::PlantSizData(PltSizCondNum).DeltaT;
                    Par(1) = DesFluidCoolerLoad;
                    Par(2) = double(FluidCoolerNum);
                    Par(3) = rho * tmpDesignWaterFlowRate; // design water mass flow rate
                    Par(4) = tmpHighSpeedAirFlowRate;      // design air volume flow rate
                    Par(5) = Cp;
                    UA0 = 0.0001 * DesFluidCoolerLoad; // Assume deltaT = 10000K (limit)
                    UA1 = DesFluidCoolerLoad;          // Assume deltaT = 1K
                    SimpleFluidCooler(FluidCoolerNum).WaterTemp = DataSizing::PlantSizData(PltSizCondNum).ExitTemp + DataSizing::PlantSizData(PltSizCondNum).DeltaT;
                    SimpleFluidCooler(FluidCoolerNum).AirTemp = SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp;
                    SimpleFluidCooler(FluidCoolerNum).AirWetBulb = SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp;
                    SimpleFluidCooler(FluidCoolerNum).AirPress = DataEnvironment::StdBaroPress;
                    SimpleFluidCooler(FluidCoolerNum).AirHumRat = Psychrometrics::PsyWFnTdbTwbPb(SimpleFluidCooler(FluidCoolerNum).AirTemp,
                                                                                      SimpleFluidCooler(FluidCoolerNum).AirWetBulb,
                                                                                      SimpleFluidCooler(FluidCoolerNum).AirPress,
                                                                                      CalledFrom);
                    General::SolveRoot(Acc, MaxIte, SolFla, UA, SimpleFluidCoolerUAResidual, UA0, UA1, Par);
                    if (SolFla == -1) {
                        ShowWarningError("Iteration limit exceeded in calculating fluid cooler UA.");
                        ShowContinueError("Autosizing of fluid cooler UA failed for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name);
                        ShowContinueError("The final UA value =" + General::RoundSigDigits(UA, 2) + " W/K, and the simulation continues...");
                    } else if (SolFla == -2) {
                        SimSimpleFluidCooler(int(Par(2)), Par(3), Par(4), UA0, OutWaterTempAtUA0);
                        SimSimpleFluidCooler(int(Par(2)), Par(3), Par(4), UA1, OutWaterTempAtUA1);
                        ShowSevereError(CalledFrom + ": The combination of design input values did not allow the calculation of a ");
                        ShowContinueError("reasonable UA value. Review and revise design input values as appropriate. Specifying hard");
                        ShowContinueError(R"(sizes for some "autosizable" fields while autosizing other "autosizable" fields may be )");
                        ShowContinueError("contributing to this problem.");
                        ShowContinueError("This model iterates on UA to find the heat transfer required to provide the design outlet ");
                        ShowContinueError("water temperature. Initially, the outlet water temperatures at high and low UA values are ");
                        ShowContinueError("calculated. The Design Exit Water Temperature should be between the outlet water ");
                        ShowContinueError("temperatures calculated at high and low UA values. If the Design Exit Water Temperature is ");
                        ShowContinueError("out of this range, the solution will not converge and UA will not be calculated. ");
                        ShowContinueError("The possible solutions could be to manually input adjusted water and/or air flow rates based ");
                        ShowContinueError("on the autosized values shown below or to adjust design fluid cooler air inlet dry-bulb temperature.");
                        ShowContinueError("Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp).");
                        ShowContinueError("Inputs to the fluid cooler object:");
                        ShowContinueError("Design Fluid Cooler Load [W]                       = " + General::RoundSigDigits(Par(1), 2));
                        ShowContinueError("Design Fluid Cooler Water Volume Flow Rate [m3/s]  = " +
                                                  General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate, 6));
                        ShowContinueError("Design Fluid Cooler Air Volume Flow Rate [m3/s]    = " + General::RoundSigDigits(Par(4), 2));
                        ShowContinueError("Design Fluid Cooler Air Inlet Dry-bulb Temp [C]    = " +
                                                  General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).AirTemp, 2));
                        ShowContinueError("Inputs to the plant sizing object:");
                        ShowContinueError("Design Exit Water Temp [C]                         = " +
                                                  General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2));
                        ShowContinueError("Loop Design Temperature Difference [C]             = " +
                                                  General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).DeltaT, 2));
                        ShowContinueError("Design Fluid Cooler Water Inlet Temp [C]           = " +
                                                  General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).WaterTemp, 2));
                        ShowContinueError("Calculated water outlet temp at low UA [C] (UA = " + General::RoundSigDigits(UA0, 2) +
                                          " W/K) = " + General::RoundSigDigits(OutWaterTempAtUA0, 2));
                        ShowContinueError("Calculated water outlet temp at high UA [C](UA = " + General::RoundSigDigits(UA1, 2) +
                                          " W/K) = " + General::RoundSigDigits(OutWaterTempAtUA1, 2));
                        ShowFatalError("Autosizing of Fluid Cooler UA failed for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    }
                    tmpHighSpeedEvapFluidCoolerUA = UA;
                    if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA = tmpHighSpeedEvapFluidCoolerUA;
                    SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity = DesFluidCoolerLoad;
                } else {
                    tmpHighSpeedEvapFluidCoolerUA = 0.0;
                    if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA = tmpHighSpeedEvapFluidCoolerUA;
                }
                if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_SingleSpd) {
                    if (DataPlant::PlantFirstSizesOkayToFinalize) {
                        if (DataPlant::PlantFinalSizesOkayToReport) {
                            ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                               SimpleFluidCooler(FluidCoolerNum).Name,
                                               "U-factor Times Area Value at Design Air Flow Rate [W/K]",
                                               SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                        }
                        if (DataPlant::PlantFirstSizesOkayToReport) {
                            ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                               SimpleFluidCooler(FluidCoolerNum).Name,
                                               "Initial U-factor Times Area Value at Design Air Flow Rate [W/K]",
                                               SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                        }
                    }
                } else if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_TwoSpd) {
                    if (DataPlant::PlantFirstSizesOkayToFinalize) {
                        if (DataPlant::PlantFinalSizesOkayToReport) {
                            ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                               SimpleFluidCooler(FluidCoolerNum).Name,
                                               "U-factor Times Area Value at High Fan Speed [W/K]",
                                               SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                        }
                        if (DataPlant::PlantFirstSizesOkayToReport) {
                            ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                               SimpleFluidCooler(FluidCoolerNum).Name,
                                               "Initial U-factor Times Area Value at High Fan Speed [W/K]",
                                               SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                        }
                    }
                }
            } else {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    ShowSevereError("Autosizing error for fluid cooler object = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    ShowFatalError("Autosizing of fluid cooler UA requires a loop Sizing:Plant object.");
                }
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num == PerfInputMethod::NOMINAL_CAPACITY) {
            if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate >= DataHVACGlobals::SmallWaterVolFlow) {
                rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                       DataGlobals::InitConvTemp,
                                                        DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                       CalledFrom);
                Cp = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                           SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp,
                                                            DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                           CalledFrom);
                DesFluidCoolerLoad = SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                Par(1) = DesFluidCoolerLoad;
                Par(2) = double(FluidCoolerNum);
                Par(3) = rho * tmpDesignWaterFlowRate; // design water mass flow rate
                Par(4) = tmpHighSpeedAirFlowRate;      // design air volume flow rate
                Par(5) = Cp;
                UA0 = 0.0001 * DesFluidCoolerLoad; // Assume deltaT = 10000K (limit)
                UA1 = DesFluidCoolerLoad;          // Assume deltaT = 1K
                SimpleFluidCooler(FluidCoolerNum).WaterTemp =
                    SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp; // design inlet water temperature
                SimpleFluidCooler(FluidCoolerNum).AirTemp =
                    SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp; // design inlet air dry-bulb temp
                SimpleFluidCooler(FluidCoolerNum).AirWetBulb =
                    SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp; // design inlet air wet-bulb temp
                SimpleFluidCooler(FluidCoolerNum).AirPress = DataEnvironment::StdBaroPress;
                SimpleFluidCooler(FluidCoolerNum).AirHumRat = Psychrometrics::PsyWFnTdbTwbPb(SimpleFluidCooler(FluidCoolerNum).AirTemp,
                                                                                  SimpleFluidCooler(FluidCoolerNum).AirWetBulb,
                                                                                  SimpleFluidCooler(FluidCoolerNum).AirPress);
                General::SolveRoot(Acc, MaxIte, SolFla, UA, SimpleFluidCoolerUAResidual, UA0, UA1, Par);
                if (SolFla == -1) {
                    ShowWarningError("Iteration limit exceeded in calculating fluid cooler UA.");
                    ShowContinueError("Autosizing of fluid cooler UA failed for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    ShowContinueError("The final UA value =" + General::RoundSigDigits(UA, 2) + " W/K, and the simulation continues...");
                } else if (SolFla == -2) {
                    SimSimpleFluidCooler(int(Par(2)), Par(3), Par(4), UA0, OutWaterTempAtUA0);
                    SimSimpleFluidCooler(int(Par(2)), Par(3), Par(4), UA1, OutWaterTempAtUA1);
                    ShowSevereError(CalledFrom + ": The combination of design input values did not allow the calculation of a ");
                    ShowContinueError("reasonable UA value. Review and revise design input values as appropriate. Specifying hard");
                    ShowContinueError(R"(sizes for some "autosizable" fields while autosizing other "autosizable" fields may be )");
                    ShowContinueError("contributing to this problem.");
                    ShowContinueError("This model iterates on UA to find the heat transfer required to provide the design outlet ");
                    ShowContinueError("water temperature. Initially, the outlet water temperatures at high and low UA values are ");
                    ShowContinueError("calculated. The Design Exit Water Temperature should be between the outlet water ");
                    ShowContinueError("temperatures calculated at high and low UA values. If the Design Exit Water Temperature is ");
                    ShowContinueError("out of this range, the solution will not converge and UA will not be calculated. ");
                    ShowContinueError("The possible solutions could be to manually input adjusted water and/or air flow rates based ");
                    ShowContinueError("on the autosized values shown below or to adjust design fluid cooler air inlet dry-bulb temperature.");
                    ShowContinueError("Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp).");
                    ShowContinueError("Inputs to the fluid cooler object:");
                    ShowContinueError("Design Fluid Cooler Load [W]                       = " + General::RoundSigDigits(Par(1), 2));
                    ShowContinueError("Design Fluid Cooler Water Volume Flow Rate [m3/s]  = " +
                                              General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate, 6));
                    ShowContinueError("Design Fluid Cooler Air Volume Flow Rate [m3/s]    = " + General::RoundSigDigits(Par(4), 2));
                    ShowContinueError("Design Fluid Cooler Air Inlet Dry-bulb Temp [C]    = " +
                                              General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).AirTemp, 2));
                    ShowContinueError("Inputs to the plant sizing object:");
                    ShowContinueError("Design Exit Water Temp [C]                         = " +
                                              General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2));
                    ShowContinueError("Loop Design Temperature Difference [C]             = " +
                                              General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).DeltaT, 2));
                    ShowContinueError("Design Fluid Cooler Water Inlet Temp [C]           = " +
                                              General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).WaterTemp, 2));
                    ShowContinueError("Calculated water outlet temp at low UA [C] (UA = " + General::RoundSigDigits(UA0, 2) +
                                      " W/K) = " + General::RoundSigDigits(OutWaterTempAtUA0, 2));
                    ShowContinueError("Calculated water outlet temp at high UA [C] (UA = " + General::RoundSigDigits(UA1, 2) +
                                      " W/K) = " + General::RoundSigDigits(OutWaterTempAtUA1, 2));
                    ShowFatalError("Autosizing of Fluid Cooler UA failed for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name);
                }
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA = UA;
            } else {
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA = 0.0;
            }
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_SingleSpd) {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Fluid cooler UA value at design air flow rate based on nominal capacity input [W/K]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Fluid cooler UA value at design air flow rate based on nominal capacity input [W/K]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                    }
                }
            } else if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_TwoSpd) {
                if (DataPlant::PlantFirstSizesOkayToFinalize) {
                    if (DataPlant::PlantFinalSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Fluid cooler UA value at high fan speed based on nominal capacity input [W/K]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                    }
                    if (DataPlant::PlantFirstSizesOkayToReport) {
                        ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                           SimpleFluidCooler(FluidCoolerNum).Name,
                                           "Initial Fluid cooler UA value at high fan speed based on nominal capacity input [W/K]",
                                           SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA);
                    }
                }
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRateWasAutoSized && DataPlant::PlantFirstSizesOkayToFinalize) {
            SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate =
                SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRateSizingFactor * SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate;
            if (DataPlant::PlantFinalSizesOkayToReport) {
                ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                   SimpleFluidCooler(FluidCoolerNum).Name,
                                   "Air Flow Rate at Low Fan Speed [m3/s]",
                                   SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate);
            }
            if (DataPlant::PlantFirstSizesOkayToReport) {
                ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                   SimpleFluidCooler(FluidCoolerNum).Name,
                                   "Initial Air Flow Rate at Low Fan Speed [m3/s]",
                                   SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate);
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPowerWasAutoSized && DataPlant::PlantFirstSizesOkayToFinalize) {
            SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower =
                SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPowerSizingFactor * SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower;
            if (DataPlant::PlantFinalSizesOkayToReport) {
                ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                   SimpleFluidCooler(FluidCoolerNum).Name,
                                   "Fan Power at Low Fan Speed [W]",
                                   SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower);
            }
            if (DataPlant::PlantFirstSizesOkayToReport) {
                ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                   SimpleFluidCooler(FluidCoolerNum).Name,
                                   "Initial Fan Power at Low Fan Speed [W]",
                                   SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower);
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUAWasAutoSized && DataPlant::PlantFirstSizesOkayToFinalize) {
            SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA =
                SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUASizingFactor * SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA;
            if (DataPlant::PlantFinalSizesOkayToReport) {
                ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                   SimpleFluidCooler(FluidCoolerNum).Name,
                                   "U-factor Times Area Value at Low Fan Speed [W/K]",
                                   SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA);
            }
            if (DataPlant::PlantFirstSizesOkayToReport) {
                ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                   SimpleFluidCooler(FluidCoolerNum).Name,
                                   "Initial U-factor Times Area Value at Low Fan Speed [W/K]",
                                   SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA);
            }
        }

        if (SimpleFluidCooler(FluidCoolerNum).PerformanceInputMethod_Num == PerfInputMethod::NOMINAL_CAPACITY &&
            SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_TwoSpd) {
            if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCapWasAutoSized && DataPlant::PlantFirstSizesOkayToFinalize) {
                SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap =
                    SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCapSizingFactor *
                    SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity;
                if (DataPlant::PlantFinalSizesOkayToReport) {
                    ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                       SimpleFluidCooler(FluidCoolerNum).Name,
                                       "Low Fan Speed Nominal Capacity [W]",
                                       SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap);
                }
                if (DataPlant::PlantFirstSizesOkayToReport) {
                    ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                       SimpleFluidCooler(FluidCoolerNum).Name,
                                       "Initial Low Fan Speed Nominal Capacity [W]",
                                       SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap);
                }
            }

            if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate >= DataHVACGlobals::SmallWaterVolFlow &&
                SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap > 0.0) {
                rho = FluidProperties::GetDensityGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                       DataGlobals::InitConvTemp,
                                                        DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                       CalledFrom);
                Cp = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                           SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp,
                                                            DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                           CalledFrom);
                DesFluidCoolerLoad = SimpleFluidCooler(FluidCoolerNum).FluidCoolerLowSpeedNomCap;
                Par(1) = DesFluidCoolerLoad;
                Par(2) = double(FluidCoolerNum);
                Par(3) = rho * tmpDesignWaterFlowRate;                          // design water mass flow rate
                Par(4) = SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate; // Air volume flow rate at low fan speed
                Par(5) = Cp;
                UA0 = 0.0001 * DesFluidCoolerLoad; // Assume deltaT = 10000K (limit)
                UA1 = DesFluidCoolerLoad;          // Assume deltaT = 1K
                SimpleFluidCooler(FluidCoolerNum).WaterTemp =
                    SimpleFluidCooler(FluidCoolerNum).DesignEnteringWaterTemp; // design inlet water temperature
                SimpleFluidCooler(FluidCoolerNum).AirTemp =
                    SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirTemp; // design inlet air dry-bulb temp
                SimpleFluidCooler(FluidCoolerNum).AirWetBulb =
                    SimpleFluidCooler(FluidCoolerNum).DesignEnteringAirWetBulbTemp; // design inlet air wet-bulb temp
                SimpleFluidCooler(FluidCoolerNum).AirPress = DataEnvironment::StdBaroPress;
                SimpleFluidCooler(FluidCoolerNum).AirHumRat = Psychrometrics::PsyWFnTdbTwbPb(SimpleFluidCooler(FluidCoolerNum).AirTemp,
                                                                                  SimpleFluidCooler(FluidCoolerNum).AirWetBulb,
                                                                                  SimpleFluidCooler(FluidCoolerNum).AirPress,
                                                                                  CalledFrom);
                General::SolveRoot(Acc, MaxIte, SolFla, UA, SimpleFluidCoolerUAResidual, UA0, UA1, Par);
                if (SolFla == -1) {
                    ShowWarningError("Iteration limit exceeded in calculating fluid cooler UA.");
                    ShowContinueError("Autosizing of fluid cooler UA failed for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name);
                    ShowContinueError("The final UA value at low fan speed =" + General::RoundSigDigits(UA, 2) + " W/C, and the simulation continues...");
                } else if (SolFla == -2) {
                    SimSimpleFluidCooler(int(Par(2)), Par(3), Par(4), UA0, OutWaterTempAtUA0);
                    SimSimpleFluidCooler(int(Par(2)), Par(3), Par(4), UA1, OutWaterTempAtUA1);
                    ShowSevereError(CalledFrom + ": The combination of design input values did not allow the calculation of a ");
                    ShowContinueError("reasonable low-speed UA value. Review and revise design input values as appropriate. ");
                    ShowContinueError(R"(Specifying hard sizes for some "autosizable" fields while autosizing other "autosizable" )");
                    ShowContinueError("fields may be contributing to this problem.");
                    ShowContinueError("This model iterates on UA to find the heat transfer required to provide the design outlet ");
                    ShowContinueError("water temperature. Initially, the outlet water temperatures at high and low UA values are ");
                    ShowContinueError("calculated. The Design Exit Water Temperature should be between the outlet water ");
                    ShowContinueError("temperatures calculated at high and low UA values. If the Design Exit Water Temperature is ");
                    ShowContinueError("out of this range, the solution will not converge and UA will not be calculated. ");
                    ShowContinueError("The possible solutions could be to manually input adjusted water and/or air flow rates based ");
                    ShowContinueError("on the autosized values shown below or to adjust design fluid cooler air inlet dry-bulb temperature.");
                    ShowContinueError("Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp).");
                    ShowContinueError("Inputs to the fluid cooler object:");
                    ShowContinueError("Design Fluid Cooler Load [W]                         = " + General::RoundSigDigits(Par(1), 2));
                    ShowContinueError("Design Fluid Cooler Water Volume Flow Rate [m3/s]    = " +
                                              General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate, 6));
                    ShowContinueError("Design Fluid Cooler Air Volume Flow Rate [m3/s]      = " + General::RoundSigDigits(Par(4), 2));
                    ShowContinueError("Design Fluid Cooler Air Inlet Dry-bulb Temp [C]      = " +
                                              General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).AirTemp, 2));
                    ShowContinueError("Inputs to the plant sizing object:");
                    ShowContinueError("Design Exit Water Temp [C]                           = " +
                                              General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).ExitTemp, 2));
                    ShowContinueError("Loop Design Temperature Difference [C]               = " +
                                              General::RoundSigDigits(DataSizing::PlantSizData(PltSizCondNum).DeltaT, 2));
                    ShowContinueError("Design Fluid Cooler Water Inlet Temp [C]             = " +
                                              General::RoundSigDigits(SimpleFluidCooler(FluidCoolerNum).WaterTemp, 2));
                    ShowContinueError("Calculated water outlet temp at low UA [C](UA = " + General::RoundSigDigits(UA0, 2) +
                                      " W/C) = " + General::RoundSigDigits(OutWaterTempAtUA0, 2));
                    ShowContinueError("Calculated water outlet temp at high UA [C](UA = " + General::RoundSigDigits(UA1, 2) +
                                      " W/C) = " + General::RoundSigDigits(OutWaterTempAtUA1, 2));
                    ShowFatalError("Autosizing of Fluid Cooler UA failed for fluid cooler = " + SimpleFluidCooler(FluidCoolerNum).Name);
                }
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA = UA;
            } else {
                if (DataPlant::PlantFirstSizesOkayToFinalize) SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA = 0.0;
            }
            if (DataPlant::PlantFirstSizesOkayToFinalize) {
                if (DataPlant::PlantFinalSizesOkayToReport) {
                    ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                       SimpleFluidCooler(FluidCoolerNum).Name,
                                       "U-factor Times Area Value at Low Fan Speed [W/C]",
                                       SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA);
                }
                if (DataPlant::PlantFirstSizesOkayToReport) {
                    ReportSizingManager::ReportSizingOutput(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType,
                                       SimpleFluidCooler(FluidCoolerNum).Name,
                                       "Initial U-factor Times Area Value at Low Fan Speed [W/C]",
                                       SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA);
                }
            }
        }

        ErrorsFound = false;

        if (DataPlant::PlantFinalSizesOkayToReport) {
            // create predefined report
            equipName = SimpleFluidCooler(FluidCoolerNum).Name;
            OutputReportPredefined::PreDefTableEntry(OutputReportPredefined::pdchMechType, equipName, SimpleFluidCooler(FluidCoolerNum).FluidCoolerType);
            OutputReportPredefined::PreDefTableEntry(OutputReportPredefined::pdchMechNomCap, equipName, SimpleFluidCooler(FluidCoolerNum).FluidCoolerNominalCapacity);
        }

        if (SimpleFluidCooler(FluidCoolerNum).FluidCoolerType_Num == DataPlant::TypeOf_FluidCooler_TwoSpd && DataPlant::PlantFirstSizesOkayToFinalize) {
            if (SimpleFluidCooler(FluidCoolerNum).DesignWaterFlowRate > 0.0) {
                if (SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate <= SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate) {
                    ShowSevereError("FluidCooler:TwoSpeed  \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Low speed air flow rate must be less than high speed air flow rate.");
                    ErrorsFound = true;
                }
                if (SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA <= SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA) {
                    ShowSevereError("FluidCooler:TwoSpeed  \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                    "\". Fluid cooler UA at low fan speed must be less than the fluid cooler UA at high fan speed.");
                    ErrorsFound = true;
                }
            }
        }

        if (ErrorsFound) {
            ShowFatalError("SizeFluidCooler: Program terminated due to previous condition(s).");
        }
    }

    void SingleSpeedFluidCooler(int &FluidCoolerNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       Dec. 2008. BG. added RunFlag logic per original methodology
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // To simulate the operation of a single-speed fan fluid cooler.

        // METHODOLOGY EMPLOYED:
        // The fluid cooler is modeled using effectiveness-NTU relationships for
        // cross flow heat exchangers (both stream unmixed)based on cooling tower model.
        // The subroutine calculates the period of time required to meet a
        // leaving water temperature setpoint. It assumes that part-load
        // operation represents a linear interpolation of two steady-state regimes.
        // Cyclic losses are neglected. The period of time required to meet the
        // leaving water temperature setpoint is used to determine the required
        // fan power and energy.
        // A RunFlag is passed by the upper level manager to indicate the ON/OFF status,
        // or schedule, of the fluid cooler. If the fluid cooler is OFF, outlet water
        // temperature and flow rate are passed through the model from inlet node to
        // outlet node without intervention. Reports are also updated with fan power
        // and energy being zero.
        // When the RunFlag indicates an ON condition for thefluid cooler, the
        // mass flow rate and water temperature are read from the inlet node of the
        // fluid cooler (water-side). The outdoor air dry-bulb temperature is used
        // as the entering condition to thefluid cooler (air-side).Thefluid cooler
        // fan is turned on and design parameters are used to calculate the leaving
        // water temperature.If the calculated leaving water temperature is below the setpoint,
        // a fan run-time fraction is calculated and used to determine fan power. The leaving
        // water temperature setpoint is placed on the outlet node. If the calculated
        // leaving water temperature is at or above the setpoint, the calculated
        // leaving water temperature is placed on the outlet node and the fan runs at
        // full power. Water mass flow rate is passed from inlet node to outlet node
        // with no intervention.

        // REFERENCES:
        // ASHRAE HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculation. 1999.
        // Based on SingleSpeedTower subroutine by Dan Fisher ,Sept 1998.

        // SUBROUTINE PARAMETER DEFINITIONS:
        static std::string const RoutineName("SingleSpeedFluidCooler");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 AirFlowRate;
        Real64 UAdesign; // UA value at design conditions (entered by user or calculated)
        Real64 OutletWaterTempOFF;
        Real64 FanModeFrac;
        Real64 FanPowerOn;
        Real64 CpWater;
        Real64 TempSetPoint = 0.0;
        int LoopNum;
        int LoopSideNum;

        // set inlet and outlet nodes
        auto & waterInletNode = SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum;
        SimpleFluidCooler(FluidCoolerNum).Qactual = 0.0;
        FanModeFrac = 0.0;
        SimpleFluidCooler(FluidCoolerNum).FanPower = 0.0;
        SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = DataLoopNode::Node(waterInletNode).Temp;
        LoopNum = SimpleFluidCooler(FluidCoolerNum).LoopNum;
        LoopSideNum = SimpleFluidCooler(FluidCoolerNum).LoopSideNum;
        {
            auto const SELECT_CASE_var(DataPlant::PlantLoop(LoopNum).LoopDemandCalcScheme);
            if (SELECT_CASE_var == DataPlant::SingleSetPoint) {
                TempSetPoint = DataPlant::PlantLoop(LoopNum).LoopSide(LoopSideNum).TempSetPoint;
            } else if (SELECT_CASE_var == DataPlant::DualSetPointDeadBand) {
                TempSetPoint = DataPlant::PlantLoop(LoopNum).LoopSide(LoopSideNum).TempSetPointHi;
            }
        }

        //   MassFlowTol is a parameter to indicate a no flow condition
        if (SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate <= DataBranchAirLoopPlant::MassFlowTolerance) return;

        if (SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp < TempSetPoint) { // already there don't need to run the cooler
            return;
        }

        //   Initialize local variables
        OutletWaterTempOFF = DataLoopNode::Node(waterInletNode).Temp;
        SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = OutletWaterTempOFF;

        UAdesign = SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA;
        AirFlowRate = SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate;
        FanPowerOn = SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower;

        SimSimpleFluidCooler(FluidCoolerNum, SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate, AirFlowRate, UAdesign, SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp);

        if (SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp <= TempSetPoint) {
            //   Setpoint was met with pump ON and fan ON, calculate run-time fraction or just wasn't needed at all
            if (SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp != OutletWaterTempOFF) { // don't divide by zero
                FanModeFrac = (TempSetPoint - OutletWaterTempOFF) / (SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp - OutletWaterTempOFF);
            }
            SimpleFluidCooler(FluidCoolerNum).FanPower = max(FanModeFrac * FanPowerOn, 0.0); // BG change
            SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = TempSetPoint;
        } else {
            //    Setpoint was not met, fluid cooler ran at full capacity
            SimpleFluidCooler(FluidCoolerNum).FanPower = FanPowerOn;
        }
        CpWater = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                                         DataLoopNode::Node(waterInletNode).Temp,
                                                         DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                        RoutineName);
        SimpleFluidCooler(FluidCoolerNum).Qactual = SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate * CpWater * (DataLoopNode::Node(waterInletNode).Temp - SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp);
    }

    void TwoSpeedFluidCooler(int &FluidCoolerNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       Dec. 2008. BG. added RunFlag logic per original methodology
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // To simulate the operation of a fluid cooler with a two-speed fan.

        // METHODOLOGY EMPLOYED:
        // The fluid cooler is modeled using effectiveness-NTU relationships for
        // cross flow heat exchangers (both stream unmixed)based on cooling tower model.
        // The subroutine calculates the period of time required to meet a
        // leaving water temperature setpoint. It assumes that part-load
        // operation represents a linear interpolation of two steady-state regimes
        // (high-speed fan operation and low-speed fan operation).
        // Cyclic losses are neglected. The period of time required to meet the
        // leaving water temperature setpoint is used to determine the required
        // fan power and energy.
        // A RunFlag is passed by the upper level manager to indicate the ON/OFF status,
        // or schedule, of the fluid cooler. If the fluid cooler is OFF, outlet water
        // temperature and flow rate are passed through the model from inlet node to
        // outlet node without intervention.Reports are also updated with fan power
        // and fan energy being zero.
        // When the RunFlag indicates an ON condition for the fluid cooler, the
        // mass flow rate and water temperature are read from the inlet node of the
        // fluid cooler (water-side). The outdoor air dry-bulb temperature is used
        // as the entering condition to the fluid cooler (air-side). Input deck
        // parameters are read for the low fan speed and a leaving water temperature
        // is calculated.
        // If the calculated leaving water temperature is below the setpoint,
        // a fan run-time fraction (FanModeFrac) is calculated and used to determine fan power.
        // The leaving water temperature setpoint is placed on the outlet node.
        // If the calculated leaving water temperature is at or above
        // the setpoint, the fluid cooler fan is turned on 'high speed' and the routine is
        // repeated. If the calculated leaving water temperature is below the setpoint,
        // a fan run-time fraction is calculated for the second stage fan and fan power
        // is calculated as FanModeFrac*HighSpeedFanPower+(1-FanModeFrac)*LowSpeedFanPower.
        // If the calculated leaving water temperature is above the leaving water temp.
        // setpoint, the calculated leaving water temperature is placed on the outlet
        // node and the fan runs at full power (High Speed Fan Power). Water mass flow
        // rate is passed from inlet node to outlet node with no intervention.

        // REFERENCES:
        // ASHRAE HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculation. 1999.
        // Based on TwoSpeedTower by Dan Fisher ,Sept. 1998.

        // SUBROUTINE PARAMETER DEFINITIONS:
        static std::string const RoutineName("TwoSpeedFluidCooler");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 AirFlowRate;
        Real64 UAdesign; // UA value at design conditions (entered by user) [W/C]
        Real64 OutletWaterTempOFF;
        Real64 OutletWaterTemp1stStage;
        Real64 OutletWaterTemp2ndStage;
        Real64 FanModeFrac;
        Real64 FanPowerLow;
        Real64 FanPowerHigh;
        Real64 CpWater;
        Real64 TempSetPoint = 0.0;
        int LoopNum;
        int LoopSideNum;

        auto & waterInletNode = SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum;
        SimpleFluidCooler(FluidCoolerNum).Qactual = 0.0;
        SimpleFluidCooler(FluidCoolerNum).FanPower = 0.0;
        SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = DataLoopNode::Node(waterInletNode).Temp;
        LoopNum = SimpleFluidCooler(FluidCoolerNum).LoopNum;
        LoopSideNum = SimpleFluidCooler(FluidCoolerNum).LoopSideNum;
        {
            auto const SELECT_CASE_var(DataPlant::PlantLoop(LoopNum).LoopDemandCalcScheme);
            if (SELECT_CASE_var == DataPlant::SingleSetPoint) {
                TempSetPoint = DataPlant::PlantLoop(LoopNum).LoopSide(LoopSideNum).TempSetPoint;
            } else if (SELECT_CASE_var == DataPlant::DualSetPointDeadBand) {
                TempSetPoint = DataPlant::PlantLoop(LoopNum).LoopSide(LoopSideNum).TempSetPointHi;
            }
        }

        // MassFlowTol is a parameter to indicate a no flow condition
        if (SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate <= DataBranchAirLoopPlant::MassFlowTolerance || DataPlant::PlantLoop(LoopNum).LoopSide(LoopSideNum).FlowLock == 0) return;

        // set local variable for fluid cooler
        SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate = DataLoopNode::Node(waterInletNode).MassFlowRate;
        OutletWaterTempOFF = DataLoopNode::Node(waterInletNode).Temp;
        OutletWaterTemp1stStage = OutletWaterTempOFF;
        OutletWaterTemp2ndStage = OutletWaterTempOFF;
        FanModeFrac = 0.0;

        if (OutletWaterTempOFF < TempSetPoint) { // already there don't need to run the cooler
            return;
        }

        UAdesign = SimpleFluidCooler(FluidCoolerNum).LowSpeedFluidCoolerUA;
        AirFlowRate = SimpleFluidCooler(FluidCoolerNum).LowSpeedAirFlowRate;
        FanPowerLow = SimpleFluidCooler(FluidCoolerNum).LowSpeedFanPower;

        SimSimpleFluidCooler(FluidCoolerNum, SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp1stStage);

        if (OutletWaterTemp1stStage <= TempSetPoint) {
            // Setpoint was met with pump ON and fan ON 1st stage, calculate fan mode fraction
            if (OutletWaterTemp1stStage != OutletWaterTempOFF) { // don't divide by zero
                FanModeFrac = (TempSetPoint - OutletWaterTempOFF) / (OutletWaterTemp1stStage - OutletWaterTempOFF);
            }
            SimpleFluidCooler(FluidCoolerNum).FanPower = FanModeFrac * FanPowerLow;
            SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = TempSetPoint;
            SimpleFluidCooler(FluidCoolerNum).Qactual *= FanModeFrac;
        } else {
            // Setpoint was not met, turn on fluid cooler 2nd stage fan
            UAdesign = SimpleFluidCooler(FluidCoolerNum).HighSpeedFluidCoolerUA;
            AirFlowRate = SimpleFluidCooler(FluidCoolerNum).HighSpeedAirFlowRate;
            FanPowerHigh = SimpleFluidCooler(FluidCoolerNum).HighSpeedFanPower;

            SimSimpleFluidCooler(FluidCoolerNum, SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp2ndStage);

            if ((OutletWaterTemp2ndStage <= TempSetPoint) && UAdesign > 0.0) {
                // Setpoint was met with pump ON and fan ON 2nd stage, calculate fan mode fraction
                FanModeFrac = (TempSetPoint - OutletWaterTemp1stStage) / (OutletWaterTemp2ndStage - OutletWaterTemp1stStage);
                SimpleFluidCooler(FluidCoolerNum).FanPower = max((FanModeFrac * FanPowerHigh) + (1.0 - FanModeFrac) * FanPowerLow, 0.0);
                SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = TempSetPoint;
            } else {
                // Setpoint was not met, fluid cooler ran at full capacity
                SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = OutletWaterTemp2ndStage;
                SimpleFluidCooler(FluidCoolerNum).FanPower = FanPowerHigh;
            }
        }
        CpWater = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                                         DataLoopNode::Node(waterInletNode).Temp,
                                                         DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                        RoutineName);
        SimpleFluidCooler(FluidCoolerNum).Qactual = SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate * CpWater * (DataLoopNode::Node(waterInletNode).Temp - SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp);
    }

    void SimSimpleFluidCooler(
        int const FluidCoolerNum, Real64 const _WaterMassFlowRate, Real64 const AirFlowRate, Real64 const UAdesign, Real64 &_OutletWaterTemp)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       April 2010, Chandan Sharma, FSEC
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // See purpose for Single Speed or Two Speed Fluid Cooler model

        // METHODOLOGY EMPLOYED:
        // See methodology for Single Speed or Two Speed Fluid Cooler model

        // Locals
        Real64 _InletWaterTemp; // Water inlet temperature
        Real64 _Qactual;        // Actual heat transfer rate between fluid cooler water and air [W]

        // SUBROUTINE PARAMETER DEFINITIONS:
        static std::string const RoutineName("SimSimpleFluidCooler");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 MdotCpWater;      // Water mass flow rate times the heat capacity [W/K]
        Real64 InletAirTemp;     // Dry-bulb temperature of air entering the fluid cooler [C]
        Real64 CpWater;          // Heat capacity of water [J/kg/K]
        Real64 CpAir;            // Heat capacity of air [J/kg/K]
        Real64 AirDensity;       // Density of air [kg/m3]
        Real64 AirMassFlowRate;  // Mass flow rate of air [kg/s]
        Real64 effectiveness;    // Effectiveness of the heat exchanger [-]
        Real64 AirCapacity;      // MdotCp of air through the fluid cooler
        Real64 CapacityRatioMin; // Minimum capacity of airside and waterside
        Real64 CapacityRatioMax; // Maximum capacity of airside and waterside
        Real64 CapacityRatio;    // Ratio of minimum to maximum capacity
        Real64 NumTransferUnits; // Number of transfer Units [NTU]
        Real64 ETA;              // initialize some local variables
        Real64 A;                // initialize some local variables

        // set local fluid cooler inlet and outlet temperature variables
        _InletWaterTemp = SimpleFluidCooler(FluidCoolerNum).WaterTemp;
        _OutletWaterTemp = _InletWaterTemp;
        InletAirTemp = SimpleFluidCooler(FluidCoolerNum).AirTemp;

        if (UAdesign == 0.0) return;

        // set water and air properties
        AirDensity =
                Psychrometrics::PsyRhoAirFnPbTdbW(SimpleFluidCooler(FluidCoolerNum).AirPress, InletAirTemp, SimpleFluidCooler(FluidCoolerNum).AirHumRat);
        AirMassFlowRate = AirFlowRate * AirDensity;
        CpAir = Psychrometrics::PsyCpAirFnWTdb(SimpleFluidCooler(FluidCoolerNum).AirHumRat, InletAirTemp);
        CpWater = FluidProperties::GetSpecificHeatGlycol(DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidName,
                                        _InletWaterTemp,
                                                         DataPlant::PlantLoop(SimpleFluidCooler(FluidCoolerNum).LoopNum).FluidIndex,
                                        RoutineName);

        // Calcluate mass flow rates
        MdotCpWater = _WaterMassFlowRate * CpWater;
        AirCapacity = AirMassFlowRate * CpAir;

        // calculate the minimum to maximum capacity ratios of airside and waterside
        CapacityRatioMin = min(AirCapacity, MdotCpWater);
        CapacityRatioMax = max(AirCapacity, MdotCpWater);
        CapacityRatio = CapacityRatioMin / CapacityRatioMax;

        // Calculate number of transfer units (NTU)
        NumTransferUnits = UAdesign / CapacityRatioMin;
        ETA = std::pow(NumTransferUnits, 0.22);
        A = CapacityRatio * NumTransferUnits / ETA;
        effectiveness = 1.0 - std::exp((std::exp(-A) - 1.0) / (CapacityRatio / ETA));

        // calculate water to air heat transfer
        _Qactual = effectiveness * CapacityRatioMin * (_InletWaterTemp - InletAirTemp);

        if (_Qactual >= 0.0) {
            _OutletWaterTemp = _InletWaterTemp - _Qactual / MdotCpWater;
        } else {
            _OutletWaterTemp = _InletWaterTemp;
        }
    }

    Real64 SimpleFluidCoolerUAResidual(Real64 const UA,          // UA of fluid cooler
                                       Array1<Real64> const &Par // par(1) = design fluid cooler load [W]
    )
    {

        // FUNCTION INFORMATION:
        //       AUTHOR         Chandan Sharma
        //       DATE WRITTEN   August 2008
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS FUNCTION:
        // Calculates residual function (Design fluid cooler load - fluid cooler Output) / Design fluid cooler load.
        // Fluid cooler output depends on the UA which is being varied to zero the residual.

        // METHODOLOGY EMPLOYED:
        // Puts UA into the fluid cooler data structure, calls SimSimpleFluidCooler, and calculates
        // the residual as defined above.

        // REFERENCES:
        // Based on SimpleTowerUAResidual by Fred Buhl, May 2002

        // par(2) = Fluid cooler number
        // par(3) = design water mass flow rate [kg/s]
        // par(4) = design air volume flow rate [m3/s]
        // par(5) = water specific heat [J/(kg*C)]

        // FUNCTION LOCAL VARIABLE DECLARATIONS:
        int FluidCoolerIndex; // index of this fluid cooler
        Real64 OutWaterTemp;  // outlet water temperature [C]
        Real64 Output;        // Fluid cooler  output [W]

        FluidCoolerIndex = int(Par(2));
        SimSimpleFluidCooler(FluidCoolerIndex, Par(3), Par(4), UA, OutWaterTemp);
        Output = Par(5) * Par(3) * (SimpleFluidCooler(FluidCoolerIndex).WaterTemp - OutWaterTemp);
        Real64 Residuum = (Par(1) - Output) / Par(1);
        return Residuum;
    }

    void UpdateFluidCooler(int const FluidCoolerNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Chandan Sharma
        //       DATE WRITTEN:    August 2008
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for passing results to the outlet water node.

        // SUBROUTINE PARAMETER DEFINITIONS:
        static ObjexxFCL::gio::Fmt LowTempFmt("(' ',F6.2)");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        std::string CharErrOut;
        std::string CharLowOutletTemp;
        int LoopNum;
        int LoopSideNum;
        Real64 LoopMinTemp;

        auto & waterOutletNode = SimpleFluidCooler(FluidCoolerNum).WaterOutletNodeNum;
        DataLoopNode::Node(waterOutletNode).Temp = SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp;

        LoopNum = SimpleFluidCooler(FluidCoolerNum).LoopNum;
        LoopSideNum = SimpleFluidCooler(FluidCoolerNum).LoopSideNum;
        if (DataPlant::PlantLoop(LoopNum).LoopSide(LoopSideNum).FlowLock == 0 || DataGlobals::WarmupFlag) return;

        // Check flow rate through fluid cooler and compare to design flow rate, show warning if greater than Design * Mulitplier
        if (DataLoopNode::Node(waterOutletNode).MassFlowRate >
            SimpleFluidCooler(FluidCoolerNum).DesWaterMassFlowRate * SimpleFluidCooler(FluidCoolerNum).FluidCoolerMassFlowRateMultiplier) {
            ++SimpleFluidCooler(FluidCoolerNum).HighMassFlowErrorCount;
            if (SimpleFluidCooler(FluidCoolerNum).HighMassFlowErrorCount < 2) {
                ShowWarningError(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType + " \"" + SimpleFluidCooler(FluidCoolerNum).Name + "\"");
                ShowContinueError(" Condenser Loop Mass Flow Rate is much greater than the fluid coolers design mass flow rate.");
                ShowContinueError(" Condenser Loop Mass Flow Rate = " + General::TrimSigDigits(DataLoopNode::Node(waterOutletNode).MassFlowRate, 6));
                ShowContinueError(" Fluid Cooler Design Mass Flow Rate   = " +
                                          General::TrimSigDigits(SimpleFluidCooler(FluidCoolerNum).DesWaterMassFlowRate, 6));
                ShowContinueErrorTimeStamp("");
            } else {
                ShowRecurringWarningErrorAtEnd(
                    SimpleFluidCooler(FluidCoolerNum).FluidCoolerType + " \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                        "\"  Condenser Loop Mass Flow Rate is much greater than the fluid coolers design mass flow rate error continues...",
                    SimpleFluidCooler(FluidCoolerNum).HighMassFlowErrorIndex,
                    DataLoopNode::Node(waterOutletNode).MassFlowRate,
                    DataLoopNode::Node(waterOutletNode).MassFlowRate);
            }
        }

        // Check if OutletWaterTemp is below the minimum condenser loop temp and warn user
        LoopMinTemp = DataPlant::PlantLoop(LoopNum).MinTemp;
        if (SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp < LoopMinTemp && SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate > 0.0) {
            ++SimpleFluidCooler(FluidCoolerNum).OutletWaterTempErrorCount;
            ObjexxFCL::gio::write(CharLowOutletTemp, LowTempFmt) << LoopMinTemp;
            ObjexxFCL::gio::write(CharErrOut, LowTempFmt) << SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp;
            strip(CharErrOut);
            if (SimpleFluidCooler(FluidCoolerNum).OutletWaterTempErrorCount < 2) {
                ShowWarningError(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType + " \"" + SimpleFluidCooler(FluidCoolerNum).Name + "\"");
                ShowContinueError(" Fluid cooler water outlet temperature (" + CharErrOut +
                                  " C) is below the specified minimum condenser loop temp of " + stripped(CharLowOutletTemp) + " C");
                ShowContinueErrorTimeStamp("");
            } else {
                ShowRecurringWarningErrorAtEnd(
                    SimpleFluidCooler(FluidCoolerNum).FluidCoolerType + " \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                        "\" Fluid cooler water outlet temperature is below the specified minimum condenser loop temp error continues...",
                    SimpleFluidCooler(FluidCoolerNum).OutletWaterTempErrorIndex,
                    SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp,
                    SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp);
            }
        }

        // Check if water mass flow rate is small (e.g. no flow) and warn user
        if (SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate > 0.0 && SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate <= DataBranchAirLoopPlant::MassFlowTolerance) {
            ++SimpleFluidCooler(FluidCoolerNum).SmallWaterMassFlowErrorCount;
            if (SimpleFluidCooler(FluidCoolerNum).SmallWaterMassFlowErrorCount < 2) {
                ShowWarningError(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType + " \"" + SimpleFluidCooler(FluidCoolerNum).Name + "\"");
                ShowContinueError(" Fluid cooler water mass flow rate near zero.");
                ShowContinueErrorTimeStamp("");
                ShowContinueError("Actual Mass flow = " + General::TrimSigDigits(SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate, 2));
            } else {
                ShowRecurringWarningErrorAtEnd(SimpleFluidCooler(FluidCoolerNum).FluidCoolerType + " \"" + SimpleFluidCooler(FluidCoolerNum).Name +
                                                   "\" Fluid cooler water mass flow rate near zero error continues...",
                                               SimpleFluidCooler(FluidCoolerNum).SmallWaterMassFlowErrorIndex,
                                               SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate,
                                               SimpleFluidCooler(FluidCoolerNum).WaterMassFlowRate);
            }
        }
    }

    void ReportFluidCooler(bool const RunFlag, int const FluidCoolerNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Chandan Sharma
        //       DATE WRITTEN:    August 2008
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine updates the report variables for the fluid cooler.

        Real64 ReportingConstant = DataHVACGlobals::TimeStepSys * DataGlobals::SecInHour;
        auto & waterInletNode = SimpleFluidCooler(FluidCoolerNum).WaterInletNodeNum;
        if (!RunFlag) {
            SimpleFluidCooler(FluidCoolerNum).InletWaterTemp = DataLoopNode::Node(waterInletNode).Temp;
            SimpleFluidCooler(FluidCoolerNum).OutletWaterTemp = DataLoopNode::Node(waterInletNode).Temp;
            SimpleFluidCooler(FluidCoolerNum).Qactual = 0.0;
            SimpleFluidCooler(FluidCoolerNum).FanPower = 0.0;
            SimpleFluidCooler(FluidCoolerNum).FanEnergy = 0.0;
        } else {
            SimpleFluidCooler(FluidCoolerNum).InletWaterTemp = DataLoopNode::Node(waterInletNode).Temp;
            SimpleFluidCooler(FluidCoolerNum).FanEnergy = SimpleFluidCooler(FluidCoolerNum).FanPower * ReportingConstant;
        }
    }

    void clear_state()
    {
        UniqueSimpleFluidCoolerNames.clear();
        GetFluidCoolerInputFlag = true;
    }

} // namespace FluidCoolers

} // namespace EnergyPlus
