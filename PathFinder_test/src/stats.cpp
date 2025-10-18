/**
 * SmartCGMS - continuous glucose monitoring and controlling framework
 * https://diabetes.zcu.cz/
 *
 * Copyright (c) since 2018 University of West Bohemia.
 *
 * Contact:
 * diabetes@mail.kiv.zcu.cz
 * Medical Informatics, Department of Computer Science and Engineering
 * Faculty of Applied Sciences, University of West Bohemia
 * Univerzitni 8, 301 00 Pilsen
 * Czech Republic
 *
 *
 * Purpose of this software:
 * This software is intended to demonstrate work of the diabetes.zcu.cz research
 * group to other scientists, to complement our published papers. It is strictly
 * prohibited to use this software for diagnosis or treatment of any medical condition,
 * without obtaining all required approvals from respective regulatory bodies.
 *
 * Especially, a diabetic patient is warned that unauthorized use of this software
 * may result into severe injure, including death.
 *
 *
 * Licensing terms:
 * Unless required by applicable law or agreed to in writing, software
 * distributed under these license terms is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *
 * a) For non-profit, academic research, this software is available under the
 *      GPLv3 license.
 * b) For any other use, especially commercial use, you must contact us and
 *       obtain specific terms and conditions for the use of the software.
 * c) When publishing work with results obtained using this software, you agree to cite the following paper:
 *       Tomas Koutny and Martin Ubl, "Parallel software architecture for the next generation of glucose
 *       monitoring", Procedia Computer Science, Volume 141C, pp. 279-286, 2018
 */

#include "stats.h"

#include <algorithm>


void CStats::Calculate_Stats() {
	mStats.avg = std::numeric_limits<double>::quiet_NaN();
	mStats.stddev = std::numeric_limits<double>::quiet_NaN();
	mStats.min = std::numeric_limits<double>::quiet_NaN();
	mStats.q25 = std::numeric_limits<double>::quiet_NaN();
	mStats.med = std::numeric_limits<double>::quiet_NaN();
	mStats.q75 = std::numeric_limits<double>::quiet_NaN();
	mStats.max = std::numeric_limits<double>::quiet_NaN();


	if (empty()) return;
		
	

	std::sort(begin(), end());	//needed to calculate quartiles
	mStats.min = operator[](0);
	mStats.max = operator[](size() - 1);
	if (mStats.min != mStats.max) {
		{ //quartiles		
			const size_t N4 = size() / 4;
			const size_t N2 = size() / 2;

			mStats.q25 = operator[](N4);
			mStats.med = operator[](N2);
			mStats.q75 = operator[](N4 + N2);
		}

		double N = static_cast<double>(size());

		{//avg
			mStats.avg = 0.0;
			for (const auto &value : *this) {
				mStats.avg += value;
			}
			mStats.avg /= N;
		}

		{ //stddev
			mStats.stddev = 0.0;
			for (const auto &value : *this) {
				const double dev = mStats.avg - value;
				mStats.stddev += dev * dev;
			}

			//first, try Unbiased estimation of standard deviation
			if (N > 1.5) N -= 1.5;
			else if (N > 1.0) N -= 1.0;	//if not, try to fall back to Bessel's Correction at least

			mStats.stddev = sqrt(mStats.stddev / N);
		}
	} else {
		//on the paper, this branch would be unecessary
		//in reality, there can be an error due to the limited number of bits, which represent the floating-point number
		//hence, completely deterministic methods such as NewUOA and Pathfinder could look like non-deterministic or erroneously implemented
		mStats.q25 = mStats.med = mStats.q75 = mStats.avg = mStats.min;
		mStats.stddev = 0.0;
	}
}

const TStats& CStats::Get_Stats() const {
	return mStats;
}