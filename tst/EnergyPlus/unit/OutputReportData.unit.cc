// EnergyPlus, Copyright (c) 1996-2016, The Board of Trustees of the University of Illinois and
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights
// reserved.
//
// If you have questions about your rights to use or distribute this software, please contact
// Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
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
//     similar designation, without Lawrence Berkeley National Laboratory's prior written consent.
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
//
// You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
// features, functionality or performance of the source code ("Enhancements") to anyone; however,
// if you choose to make your Enhancements available either publicly, or directly to Lawrence
// Berkeley National Laboratory, without imposing a separate written license agreement for such
// Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free
// perpetual license to install, use, modify, prepare derivative works, incorporate into other
// computer software, distribute, and sublicense such enhancements or derivative works thereof,
// in binary and source code form.

// EnergyPlus::OutputReportTabular Unit Tests

// Google Test Headers
#include <gtest/gtest.h>
// ObjexxFCL Headers
#include <ObjexxFCL/Array1D.hh>
// EnergyPlus Headers
#include "Fixtures/EnergyPlusFixture.hh"
#include <EnergyPlus/OutputReportData.hh>
#include <EnergyPlus/UtilityRoutines.hh>
#include <EnergyPlus/OutputProcessor.hh>
#include <EnergyPlus/DataOutputs.hh>

using namespace EnergyPlus;
using namespace ObjexxFCL;
using namespace OutputProcessor;
using namespace DataOutputs;
TEST_F( EnergyPlusFixture, OutputReportData_AnnualFieldSetConstructor )
{
	std::string varNameTest = "TestReport";
	AnnualFieldSet::AggregationKind kindOfAggregationTest = AnnualFieldSet::AggregationKind::sumOrAvg;
	int numDigitsShownTest = 3;
	AnnualFieldSet fldStTest = AnnualFieldSet( varNameTest, kindOfAggregationTest, numDigitsShownTest );
	EXPECT_EQ( fldStTest.m_variMeter, varNameTest );
	EXPECT_EQ( fldStTest.m_aggregate, kindOfAggregationTest );
	EXPECT_EQ( fldStTest.m_showDigits, numDigitsShownTest );
}

TEST_F( EnergyPlusFixture, OutputReportData_getVariableKeys )
{
	std::string varNameTest = "TestReport";
	AnnualFieldSet::AggregationKind kindOfAggregationTest = AnnualFieldSet::AggregationKind::sumOrAvg;
	int numDigitsShownTest = 3;
	AnnualFieldSet fldStTest = AnnualFieldSet( varNameTest, kindOfAggregationTest, numDigitsShownTest );

	Real64 extLitPow;
	Real64 extLitUse;

	SetupOutputVariable( "Exterior Lights Electric Energy [J]", extLitUse, "Zone", "Sum", "Lite1", _, "Electricity", "Exterior Lights", "General" );
	SetupOutputVariable( "Exterior Lights Electric Energy [J]", extLitUse, "Zone", "Sum", "Lite2", _, "Electricity", "Exterior Lights", "General" );
	SetupOutputVariable( "Exterior Lights Electric Energy [J]", extLitUse, "Zone", "Sum", "Lite3", _, "Electricity", "Exterior Lights", "General" );
	SetupOutputVariable( "Exterior Lights Electric Power [W]", extLitPow, "Zone", "Average", "Lite1" );
	SetupOutputVariable( "Exterior Lights Electric Power [W]", extLitPow, "Zone", "Average", "Lite2" );
	SetupOutputVariable( "Exterior Lights Electric Power [W]", extLitPow, "Zone", "Average", "Lite3" );

	int keyCount = 0;
	int typeVar = 0;
	int avgSumVar = 0;
	int stepTypeVar = 0;
	std::string unitsVar = "";

	fldStTest.m_variMeter = "EXTERIOR LIGHTS ELECTRIC ENERGY";
	keyCount = fldStTest.getVariableKeyCountandTypeFromFldSt( typeVar, avgSumVar, stepTypeVar, unitsVar );
	EXPECT_EQ( keyCount, 3 );

	fldStTest.getVariableKeysFromFldSt( typeVar, keyCount, fldStTest.m_namesOfKeys, fldStTest.m_indexesForKeyVar );

	EXPECT_EQ( fldStTest.m_namesOfKeys[0], "LITE1" );
	EXPECT_EQ( fldStTest.m_namesOfKeys[1], "LITE2" );
	EXPECT_EQ( fldStTest.m_namesOfKeys[2], "LITE3" );
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "Outside Air Inlet Node,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "Relief Air Outlet Node,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "(Relief|Outside) Air (Outlet|Inlet) Node,", "System Node Temperature,", "timestep;",
			                                                 " Output:Variable,", "Mixed Air Node,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "(Mixed|Single) Air Node,", "System Node Temperature,", "timestep;",
			                                                 " Output:Variable,", "*,", "Unitary System Compressor Part Load Ratio,", "timestep;",
			                                                 " Output:Variable,", ".*,", "Zone Air System Sensible Heating Rate,", "timestep;",
			                                                 " Output:Variable,", "SALESFLOOR OUTLET NODE,", "System Node Temperature,", "timestep;",
			                                                 " Output:Variable,", "BackRoom(.*),", "System Node Temperature,", "timestep;",
			                                                 " Output:Variable,", "(.*)N(ode|ODE),", "System Node Humidity Ratio,", "timestep;",
	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 5);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 10);

	OnTheList = FindItemInVariableList( "Outside Air Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "OUTSIDE AIR INLET NODE", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Mixed Air Node", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Inlet Node", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Inlet Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Any Node Here", "Zone Air System Sensible Heating Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Salesfloor Outlet Node", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "AnySalesfloor Outlet Node", "System Node Temperature" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "AnyOutside Air Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom OUTLET NODE", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Inlet NODE", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Node", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom OUTLET NODE", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Inlet NODE", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Any Node", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_Plus )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "(.+)Inlet(.+),", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "(.+)Inlet,", "System Node Humidity Ratio,", "timestep;",
			                                                 " Output:Variable,", "(.+)Node,", "Zone Air System Sensible Heating Rate,", "timestep;",
			                                                 " Output:Variable,", "Outside Air (.+) Node,", "Unitary System Compressor Part Load Ratio,", "timestep;",
			                                                 " Output:Variable,", "Outside Air .+ Node,", "Unitary System Load Ratio,", "timestep;",
			                                                 " Output:Variable,", ".+,", "System Node Temperature,", "timestep;",
	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 6);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 6);

	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor INLET Node", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Inlet Node", "System Node Humidity Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Any Node", "Zone Air System Sensible Heating Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Inlet Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Outlet Node", "Unitary System Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Any Node", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "", "System Node Temperature" );
	EXPECT_NE(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_Star )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "(.*)Inlet(.*),", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "(.*)Inlet,", "System Node Humidity Ratio,", "timestep;",
			                                                 " Output:Variable,", "(.*)Node,", "Zone Air System Sensible Heating Rate,", "timestep;",
			                                                 " Output:Variable,", "Outside Air(.*) Node,", "Unitary System Compressor Part Load Ratio,", "timestep;",
			                                                 " Output:Variable,", "Outside Air.* Node,", "Unitary System Load Ratio,", "timestep;",
			                                                 " Output:Variable,", ".*,", "System Node Temperature,", "timestep;",
			                                                 " Output:Variable,", "*,", "Refrigeration Compressor Rack Electric Power,", "timestep;",
	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 7);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 7);

	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor INLET Node", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Inlet Node", "System Node Humidity Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Any Inlet", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom Any Node", "Zone Air System Sensible Heating Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Node", "Zone Air System Sensible Heating Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "NODE", "Zone Air System Sensible Heating Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Inlet Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Outlet Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Outlet Node", "Unitary System Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside Air Node", "Unitary System Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outside AirNode", "Unitary System Load Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Any Node", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "", "System Node Temperature" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Any Node", "Refrigeration Compressor Rack Electric Power" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "", "Refrigeration Compressor Rack Electric Power" );
	EXPECT_EQ(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_Pipe )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "SalesFloor I(nlet|NLET) Node,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "SalesFloor O(utlet|UTLET) Node,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "System (Inlet|Outlet) Node,", "Unitary System Compressor Part Load Ratio,", "timestep;",
			                                                 " Output:Variable,", "(BackRoom|BACKROOM|SALESFLOOR|SalesFloor) (Outlet|OUTLET) (NODE|Node),", "System Node Humidity Ratio,", "timestep;",
	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 3);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 4);

	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor INLET Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor Outlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor OUTLET Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "System Inlet Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "System Outlet Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "System Another Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom OUTLET NODE", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SALESFLOOR OUTLET NODE", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor Outlet Node", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "BACKROOM Outlet Node", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_Brackets )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "([A-Za-z] ?)+,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "[A-Za-z0-9_]+,", "System Node Humidity Ratio,", "timestep;",
			                                                 " Output:Variable,", "[A-Z]{4},", "Unitary System Compressor Part Load Ratio,", "timestep;",
			                                                 " Output:Variable,", "[A-Za-z]{5,6},", "Zone Air System Sensible Heating Rate,", "timestep;",
	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 4);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 4);

	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom OUTLET NODE", "System Node Humidity Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "BackRoom_NODE1", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "NODE", "Unitary System Compressor Part Load Ratio" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Node", "Unitary System Compressor Part Load Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "NOD", "Unitary System Compressor Part Load Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet", "Zone Air System Sensible Heating Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outlet", "Zone Air System Sensible Heating Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Any Node", "Zone Air System Sensible Heating Rate" );
	EXPECT_NE(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_SpecChars )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "\\w,", "System Node Mass Flow Rate,", "timestep;",

	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	//bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 1);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 1);

//	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
//	EXPECT_NE(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_Carrot )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "^Inlet(.*)Node,", "System Node Mass Flow Rate,", "timestep;",
			                                                 " Output:Variable,", "[^0-9]+,", "System Node Humidity Ratio,", "timestep;",

	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 2);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 2);

	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet System Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor1", "System Node Humidity Ratio" );
	EXPECT_NE(true, OnTheList);
	OnTheList = FindItemInVariableList( "SalesFloor", "System Node Humidity Ratio" );
	EXPECT_EQ(true, OnTheList);
}

TEST_F( EnergyPlusFixture, OutputReportData_Regex_Dollar )
{
	std::string const idf_objects = delimited_string({
			                                                 "Version,8.6;",
			                                                 " Output:Variable,", "(.*)Node$,", "System Node Mass Flow Rate,", "timestep;",
	                                                 });
	ASSERT_TRUE( process_idf( idf_objects ) );
	bool OnTheList;

	EXPECT_EQ( (int)OutputVariablesNames.size(), 1);
	EXPECT_EQ( (int)OutputVariablesForSimulation.size(), 1);

	OnTheList = FindItemInVariableList( "SalesFloor Inlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Outlet Node", "System Node Mass Flow Rate" );
	EXPECT_EQ(true, OnTheList);
	OnTheList = FindItemInVariableList( "Inlet Node1 ", "System Node Mass Flow Rate" );
	EXPECT_NE(true, OnTheList);
}