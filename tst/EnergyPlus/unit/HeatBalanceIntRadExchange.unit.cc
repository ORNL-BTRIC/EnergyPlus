// EnergyPlus, Copyright (c) 1996-2018, The Board of Trustees of the University of Illinois,
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

// EnergyPlus::HeatBalFiniteDiffManager Unit Tests

// Google Test Headers
#include <gtest/gtest.h>

// EnergyPlus Headers
#include "Fixtures/EnergyPlusFixture.hh"
#include <EnergyPlus/DataGlobals.hh>
#include <EnergyPlus/DataHeatBalSurface.hh>
#include <EnergyPlus/DataHeatBalance.hh>
#include <EnergyPlus/DataSurfaces.hh>
#include <EnergyPlus/HeatBalanceIntRadExchange.hh>

using namespace EnergyPlus::HeatBalanceIntRadExchange;

namespace EnergyPlus {

	TEST_F( EnergyPlusFixture, HeatBalanceIntRadExchange_FixViewFactorsTest)
	{

		int N; // NUMBER OF SURFACES
		Array1D< Real64 > A; // AREA VECTOR- ASSUMED,BE N ELEMENTS LONG
		Array2D< Real64 > F; // APPROXIMATE DIRECT VIEW FACTOR MATRIX (N X N)
		int ZoneNum; // Zone number being fixed
		Real64 OriginalCheckValue; // check of SUM(F) - N
		Real64 FixedCheckValue; // check after fixed of SUM(F) - N
		Real64 FinalCheckValue; // the one to go with
		int NumIterations; // number of iterations to fixed
		Real64 RowSum; // RowSum of Fixed

		N = 3;

		A.allocate( N );
		F.allocate( N, N );

		A( 1 ) = 1.0;
		A( 2 ) = 1.0;
		A( 3 ) = 1.0;

		F( 1, 1 ) = 0.0;
		F( 1, 2 ) = 0.5;
		F( 1, 3 ) = 0.5;
		F( 2, 1 ) = 0.5;
		F( 2, 2 ) = 0.0;
		F( 2, 3 ) = 0.5;
		F( 3, 1 ) = 0.5;
		F( 3, 2 ) = 0.5;
		F( 3, 3 ) = 0.0;

		ZoneNum = 1;

		DataHeatBalance::Zone.allocate( ZoneNum );
		DataHeatBalance::Zone( ZoneNum ).Name = "Test";

		FixViewFactors( N, A, F, ZoneNum, OriginalCheckValue, FixedCheckValue, FinalCheckValue, NumIterations, RowSum );

		std::string const error_string = delimited_string( {
			"   ** Warning ** Surfaces in Zone=\"Test\" do not define an enclosure.",
			"   **   ~~~   ** Number of surfaces <= 3, view factors are set to force reciprocity but may not fulfill completeness.",
			"   **   ~~~   ** Reciprocity means that radiant exchange between two surfaces will match and not lead to an energy loss.",
			"   **   ~~~   ** Completeness means that all of the view factors between a surface and the other surfaces in a zone add up to unity.",
			"   **   ~~~   ** So, when there are three or less surfaces in a zone, EnergyPlus will make sure there are no losses of energy but",
			"   **   ~~~   ** it will not exchange the full amount of radiation with the rest of the zone as it would if there was a completed enclosure.",
		} );

		EXPECT_TRUE( compare_err_stream( error_string, true ) );

		// Tests for correction of view factors based on GitHub Issue #5772

		A( 1 ) = 20.0;
		A( 2 ) = 180.0;
		A( 3 ) = 180.0;
		F( 1, 1 ) = 0.0;
		F( 1, 2 ) = 0.5;
		F( 1, 3 ) = 0.5;
		F( 2, 1 ) = 0.1;
		F( 2, 2 ) = 0.0;
		F( 2, 3 ) = 0.9;
		F( 3, 1 ) = 0.1;
		F( 3, 2 ) = 0.9;
		F( 3, 3 ) = 0.0;

		FixViewFactors( N, A, F, ZoneNum, OriginalCheckValue, FixedCheckValue, FinalCheckValue, NumIterations, RowSum );
		EXPECT_NEAR( F( 1, 2 ), 0.07986, 0.001 );
		EXPECT_NEAR( F( 2, 1 ), 0.71875, 0.001 );
		EXPECT_NEAR( F( 3, 2 ), 0.28125, 0.001 );

		A( 1 ) = 100.0;
		A( 2 ) = 100.0;
		A( 3 ) = 200.0;
		F( 1, 1 ) = 0.0;
		F( 1, 2 ) = 1.0/3.0;
		F( 1, 3 ) = 2.0/3.0;
		F( 2, 1 ) = 1.0/3.0;
		F( 2, 2 ) = 0.0;
		F( 2, 3 ) = 2.0/3.0;
		F( 3, 1 ) = 0.5;
		F( 3, 2 ) = 0.5;
		F( 3, 3 ) = 0.0;

		FixViewFactors( N, A, F, ZoneNum, OriginalCheckValue, FixedCheckValue, FinalCheckValue, NumIterations, RowSum );
		EXPECT_NEAR( F( 1, 2 ), 0.181818, 0.001 );
		EXPECT_NEAR( F( 2, 3 ), 0.25, 0.001 );
		EXPECT_NEAR( F( 3, 2 ), 0.5, 0.001 );

		A( 1 ) = 100.0;
		A( 2 ) = 150.0;
		A( 3 ) = 200.0;
		F( 1, 1 ) = 0.0;
		F( 1, 2 ) = 150.0/350.0;
		F( 1, 3 ) = 200.0/350.0;
		F( 2, 1 ) = 1.0/3.0;
		F( 2, 2 ) = 0.0;
		F( 2, 3 ) = 2.0/3.0;
		F( 3, 1 ) = 0.4;
		F( 3, 2 ) = 0.6;
		F( 3, 3 ) = 0.0;

		FixViewFactors( N, A, F, ZoneNum, OriginalCheckValue, FixedCheckValue, FinalCheckValue, NumIterations, RowSum );
		EXPECT_NEAR( F( 1, 2 ), 0.21466, 0.001 );
		EXPECT_NEAR( F( 1, 3 ), 0.25445, 0.001 );
		EXPECT_NEAR( F( 2, 1 ), 0.32199, 0.001 );
		EXPECT_NEAR( F( 2, 3 ), 0.36832, 0.001 );
		EXPECT_NEAR( F( 3, 1 ), 0.50890, 0.001 );
		EXPECT_NEAR( F( 3, 2 ), 0.49110, 0.001 );

		A.deallocate();
		F.deallocate();

	}

	TEST_F( EnergyPlusFixture, HeatBalanceIntRadExchange_UpdateMovableInsulationFlagTest)
	{

		bool DidMIChange;
		int SurfNum;

		DataHeatBalance::Construct.allocate( 1 );
		DataHeatBalance::Material.allocate( 1 );
		DataSurfaces::Surface.allocate( 1 );

		SurfNum = 1;
		DataSurfaces::Surface( 1 ).MaterialMovInsulInt = 1;
		DataSurfaces::Surface( 1 ).MovInsulIntPresent = false;
		DataSurfaces::Surface( 1 ).MovInsulIntPresentPrevTS = false;
		DataSurfaces::Surface( 1 ).Construction = 1;
		DataSurfaces::Surface( 1 ).MaterialMovInsulInt = 1;
		DataHeatBalance::Construct( 1 ).InsideAbsorpThermal = 0.9;
		DataHeatBalance::Material( 1 ).AbsorpThermal = 0.5;
		DataHeatBalance::Material( 1 ).Resistance = 1.25;
		DataSurfaces::Surface( 1 ).SchedMovInsulInt = -1;
		DataHeatBalance::Material( 1 ).AbsorpSolar = 0.25;

		// Test 1: Movable insulation present but wasn't in previous time step, also movable insulation emissivity different than base construction
		//         This should result in a true value from the algorithm which will cause interior radiant exchange matrices to be recalculated
		HeatBalanceIntRadExchange::UpdateMovableInsulationFlag( DidMIChange, SurfNum );
		EXPECT_TRUE( DidMIChange );

		// Test 2: Movable insulation present and was also present in previous time step.  This should result in a false value since nothing has changed.
		DataSurfaces::Surface( 1 ).MovInsulIntPresentPrevTS = true;
		HeatBalanceIntRadExchange::UpdateMovableInsulationFlag( DidMIChange, SurfNum );
		EXPECT_TRUE( !DidMIChange );

		// Test 2: Movable insulation present but wasn't in previous time step.  However, the emissivity of the movable insulation and that of the
  		// 		   construction are the same so nothing has actually changed.  This should result in a false value.
		DataSurfaces::Surface( 1 ).MovInsulIntPresentPrevTS = false;
		DataHeatBalance::Material( 1 ).AbsorpThermal = DataHeatBalance::Construct( 1 ).InsideAbsorpThermal;
		HeatBalanceIntRadExchange::UpdateMovableInsulationFlag( DidMIChange, SurfNum );
		EXPECT_TRUE( !DidMIChange );

	}

//	TEST_F( EnergyPlusFixture, HeatBalanceIntRadExchange_EigenScriptF) {
//
//
//		 std::vector<double> A = {73.200000000000003, 12.558216433873083, 12.558216433873085, 55.440000000000005,
//						 99.160000000000011, 99.160000000000011};
//
//		std::vector< std::vector<double> > F = {
//				{0, 0.19960798646887035, 0.19960798646887032, 0.22557917636479394, 0.28092364369066486, 0.28092364369066486},
//				{0.03424481278696289, 0, 0.027873485297038016, 0.031783681166627727, 0.040047073848088373, 0.040047073848088373},
//				{0.03424481278696289, 0.027873485297038019, 0, 0.031783681166627734, 0.040047073848088373, 0.040047073848088373},
//				{0.17084849095169641, 0.14031349858925751, 0.14031349858925751, 0, 0.19882295509714726, 0.19882295509714726},
//				{0.38055175557877496, 0.3162127252453894, 0.3162127252453894, 0.35561479486711978, 0, 0.43986318709772698},
//				{0.38055175557877496, 0.3162127252453894, 0.3162127252453894, 0.35561479486711978, 0.43986318709772698, 0}
//		};
//
//		auto const N = 6;
//
//		std::vector<double> Emiss = {
//				0.90000000000000002, 0.90000000000000002, 0.90000000000000002, 0.90000000000000002,
//				0.90000000000000002, 0.65000000000000002
//		};
//
////		std::vector<std::vector<double>> scriptF(6, std::vector<double>(6));
//
//		Eigen::Matrix<Real64, Eigen::Dynamic, Eigen::Dynamic> scriptF(6, 6);
//
//		CalcScriptF(6, A, F, Emiss, scriptF);
//
//		auto const b = 5.6697e-8;
//		scriptF *= b;
//
////		std::vector<std::vector<double>> expectedScriptF = {
////				{0.048532695662563359, 0.034715888073639685, 0.034715888073639692, 0.17083130738775773, 0.36988200893234513, 0.24168020164255244},
////				{0.20235381515929896, 0.0064654109430229987, 0.028980177532077409, 0.14360013674321756, 0.31376147881716809, 0.20501115392738298},
////				{0.20235381515929904, 0.028980177532077402, 0.0064654109430259963, 0.14360013674321759, 0.3137614788171682, 0.20501115392738301},
////				{0.22555648810937706, 0.032528167336852383, 0.03252816733685239, 0.033628496809217652, 0.34841186109414951, 0.2276516478510118},
////				{0.27304722724735447, 0.039736633315833832, 0.039736633315833832, 0.19479582068434503, 0.078265650059336514, 0.27412878490770132},
////				{0.17840853933274353, 0.025963840685539984, 0.025963840685539988, 0.12727921900827044, 0.27412878490770137, 0.018066779673117081}
////		};
//
//		std::vector<std::vector<double>> expectedScriptF = {
//				{0.048532695662563359 * b, 0.034715888073639685 * b, 0.034715888073639692 * b, 0.17083130738775773 * b, 0.36988200893234513 * b, 0.24168020164255244 * b},
//				{0.20235381515929896 * b, 0.0064654109430229987 * b, 0.028980177532077409 * b, 0.14360013674321756 * b, 0.31376147881716809 * b, 0.20501115392738298 * b},
//				{0.20235381515929904 * b, 0.028980177532077402 * b, 0.0064654109430259963 * b, 0.14360013674321759 * b, 0.3137614788171682 * b, 0.20501115392738301 * b},
//				{0.22555648810937706 * b, 0.032528167336852383 * b, 0.03252816733685239 * b, 0.033628496809217652 * b, 0.34841186109414951 * b, 0.2276516478510118 * b},
//				{0.27304722724735447 * b, 0.039736633315833832 * b, 0.039736633315833832 * b, 0.19479582068434503 * b, 0.078265650059336514 * b, 0.27412878490770132 * b},
//				{0.17840853933274353 * b, 0.025963840685539984 * b, 0.025963840685539988 * b, 0.12727921900827044 * b, 0.27412878490770137 * b, 0.018066779673117081 * b}
//		};
//
//		std::vector<std::vector<double>> actualScriptF = {
//				{scriptF(0, 0), scriptF(0, 1), scriptF(0, 2), scriptF(0, 3), scriptF(0, 4), scriptF(0, 5)},
//				{scriptF(1, 0), scriptF(1, 1), scriptF(1, 2), scriptF(1, 3), scriptF(1, 4), scriptF(1, 5)},
//				{scriptF(2, 0), scriptF(2, 1), scriptF(2, 2), scriptF(2, 3), scriptF(2, 4), scriptF(2, 5)},
//				{scriptF(3, 0), scriptF(3, 1), scriptF(3, 2), scriptF(3, 3), scriptF(3, 4), scriptF(3, 5)},
//				{scriptF(4, 0), scriptF(4, 1), scriptF(4, 2), scriptF(4, 3), scriptF(4, 4), scriptF(4, 5)},
//				{scriptF(5, 0), scriptF(5, 1), scriptF(5, 2), scriptF(5, 3), scriptF(5, 4), scriptF(5, 5)}
//		};
//
//		compare_containers(expectedScriptF, actualScriptF);
//	}

}
