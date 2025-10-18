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

#pragma once

#include <scgms/iface/SolverIface.h>


#include <atomic>
#include <vector>
#include <algorithm>
#include <string>
#include <mutex>

extern const double pi; 
extern const double E;

#undef max
#undef min


class CSolution : public std::vector<double> {
protected:
public:
	CSolution();
	void setConstant(const double new_value, const size_t new_size = std::numeric_limits<size_t>::max());		
};


class CCommon_Problem {
protected:
	const size_t mProblem_Size;
	double mUpper_Bound_1D = std::numeric_limits<double>::quiet_NaN();	
	double mOptimum_Fitness = std::numeric_limits<double>::quiet_NaN(); //cached value to save some calls from Check_Objective_Call
	
	CSolution mAnalytic_Optimum, mOptimum_Parameters, mShift, mParams001;
	std::atomic<bool> m001_reached{ false }, mOptimum_Reached{ false };
	std::mutex mStat_Guard;
protected:
	//double mShift = -4.0;	//e.g., to eliminate too good results which occurs by a chance when the algorithm constructs
								//initial estimate as average of the bounds and accidentally hits the optimum
	std::atomic<uint64_t> mObjective_Calls{ 0 };	//number of callings to the objective functions
	uint64_t mLeast_Objective_Call = std::numeric_limits<uint64_t>::max() ;		//number of objective call, when the difference between optimum got to zero
	uint64_t mLeast_Objective_Call_001 = std::numeric_limits<uint64_t>::max();		//number of objective call, when the difference between optimum got to zero

protected:
	inline void Check_Objective_Call(const double fitness, const double *solution);
	virtual std::unique_ptr<CCommon_Problem> Clone_Internal() = 0;
	virtual const char* Get_Name_Internal() = 0;	
public:
	CCommon_Problem(const size_t problem_size);
	virtual ~CCommon_Problem() {};

	virtual void Init_Optimum();

	std::unique_ptr<CCommon_Problem> Clone();
		
	void get_bounds(CSolution &lower, CSolution &upper);
	void get_optimum(CSolution &parameters, double &fitness);
	void randomize_shift();
	void reset_counters();

	void Get_Objective_Calls(double &total, double &least, double &least001, CSolution &params001);
	virtual double Calculate_Fitness(const double *solution) = 0;
	virtual bool Can_Be_Solved();	//tests whether the problem can be solved for the given problem size
	std::string Get_Name();
	size_t Problem_Size();
};


using TProblem_Collection = std::vector<std::unique_ptr<CCommon_Problem>>;

TProblem_Collection Create_Problem_Collection(const size_t problem_size);

BOOL IfaceCalling Fitness_Wrapper(const void* data, const size_t count, const double* solution, double* const fitnesss);
