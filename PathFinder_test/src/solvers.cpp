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

#include "solvers.h"

#include <scgms/rtl/scgmsLib.h>
#include <scgms/rtl/SolverLib.h>

#include <iostream>
#include <vector>
#include <chrono>
#include <set>
#include <numeric>
#include <map>

namespace diagnostic {
	namespace mt_metade {	//mersenne twister initialized with linear random generator
		constexpr GUID id = { 0x1b21b62f, 0x7c6c, 0x4027,{ 0x89, 0xbc, 0x68, 0x7d, 0x8b, 0xd3, 0x2b, 0x3c } };	// {1B21B62F-7C6C-4027-89BC-687D8BD32B3C}
	}

	namespace halton_metade {
		constexpr GUID id = { 0x1274b08, 0xf721, 0x42bc, { 0xa5, 0x62, 0x5, 0x56, 0x71, 0x4c, 0x56, 0x85 } };
	}

	namespace rnd_metade {		//std::random_device, should be cryprographicallly secure depending on the implementation
		constexpr GUID id = { 0x2332f9a7, 0x39a2, 0x4fd6, { 0x90, 0xd5, 0x90, 0xb8, 0x85, 0x20, 0x18, 0x69 } };
	}

	namespace pathfinder {
		constexpr GUID id_fast = { 0x787223e7, 0x6363, 0x41d0, { 0xb1, 0x3d, 0x93, 0xb5, 0xd4, 0xd9, 0x4b, 0xbd } };
		constexpr GUID id_spiral = { 0xbd3baa39, 0xc447, 0x436f, { 0xbf, 0x6c, 0xdf, 0x21, 0x38, 0xd5, 0x31, 0xd4 } };
		constexpr GUID id_landscape = { 0xc919978d, 0x1b6f, 0x4b34, { 0xb3, 0xe5, 0xf8, 0x38, 0xaf, 0xe, 0x8e, 0xfb } };

	}

	//The nlopt, pagmo and ppr solvers require a special library, which does not ship with SmartCGMS release by default. You may request it.

	namespace nlopt {
		constexpr GUID newuoa_id = { 0xfd4f3f19, 0xcd6b, 0x4598,{ 0x86, 0x32, 0x40, 0x84, 0x7a, 0xad, 0x9f, 0x5 } };	// {FD4F3F19-CD6B-4598-8632-40847AAD9F05}
		constexpr GUID bobyqa_id = { 0x46559ef3, 0xa11c, 0x4387,{ 0x8f, 0x96, 0x47, 0x2d, 0x26, 0x71, 0xa, 0x1e } };	// {46559EF3-A11C-4387-8F96-472D26710A1E}
		constexpr GUID simplex_id = { 0x5f19b7de, 0x2a16, 0x4e2b, { 0x80, 0x9, 0x25, 0x1a, 0x2e, 0x9e, 0x1, 0xf0 } };
		constexpr GUID subplex_id = { 0x8342205, 0x2014, 0x4709, { 0x85, 0x76, 0x72, 0x95, 0x96, 0xc, 0x63, 0x78 } };
		constexpr GUID praxis_id = { 0x4b29deea, 0x72b8, 0x42e3, { 0xab, 0x56, 0x5, 0x39, 0xe4, 0xda, 0xab, 0xa0 } };

	}

	namespace pagmo {
		constexpr GUID pso_id = { 0x48b7e2f6, 0xa915, 0x4b63, { 0xb0, 0xf7, 0x18, 0x3e, 0xc0, 0x9b, 0x2, 0x5d } };
		constexpr GUID sade_id = { 0x90f4d682, 0xd8e, 0x4126, { 0x95, 0xc4, 0x9d, 0xb6, 0x19, 0x63, 0x3a, 0xa8 } };
		constexpr GUID de1220_id = { 0xb5ca8160, 0xf646, 0x4d53, { 0x8d, 0x66, 0xd9, 0xa8, 0xab, 0x25, 0x19, 0xbc } };
		constexpr GUID abc_id = { 0x1663854d, 0xed2e, 0x4493, { 0xa6, 0xe4, 0x46, 0xcc, 0x3f, 0xec, 0x98, 0x34 } };
		constexpr GUID cmaes_id = { 0x4e44d9f0, 0xd5d2, 0x430a, { 0x99, 0x8, 0x98, 0x13, 0xb0, 0x7c, 0x77, 0xc6 } };
		constexpr GUID xnes_id = { 0x100f6539, 0x63bd, 0x4a64, { 0x82, 0xe9, 0x86, 0x70, 0xb8, 0x0, 0xca, 0x11 } };
		constexpr GUID gpso_id = { 0x1f24727f, 0xe423, 0x4e2a, { 0xb8, 0xfe, 0xd3, 0x1, 0x64, 0xf7, 0x5a, 0x29 } };
		constexpr GUID ihs_id = { 0xc3522bcd, 0x9ddc, 0x4655, { 0xaf, 0x7f, 0xe, 0x24, 0xdf, 0x66, 0x2d, 0xa0 } };
	}


	namespace ppr {
		constexpr GUID spo_id = { 0xf8461a8e, 0x365, 0x4b90, { 0x8c, 0x99, 0xad, 0x54, 0x3a, 0x0, 0xa7, 0x8c } };
	}

	namespace sequential_brute_force_scan {
		constexpr GUID id = { 0x33d92b0, 0xb49c, 0x45d1, { 0x95, 0x7f, 0x57, 0x68, 0x2d, 0x56, 0xab, 0xd2 } };  // {033D92B0-B49C-45D1-957F-57682D56ABD2}
	}

	namespace pso {
		constexpr GUID id = { 0x93fff43a, 0x50e8, 0x4c7b, { 0x82, 0xf0, 0xf2, 0x90, 0xdf, 0xf2, 0x8, 0x9c } };	// {93FFF43A-50E8-4C7B-82F0-F290DFF2089C}
	}

	namespace rumoropt {
		constexpr GUID id = { 0x6bad021f, 0x6f68, 0x4246, { 0xa4, 0xe6, 0x7b, 0x19, 0x50, 0xca, 0x71, 0xcb } };	// {6BAD021F-6F68-4246-A4E6-7B1950CA71CB}
	}


	constexpr bool debugging = true;

	constexpr bool run_deterministic_solver_once = true;
	std::set<GUID> deterministic_solvers = {nlopt::newuoa_id, nlopt::bobyqa_id, nlopt::simplex_id, nlopt::subplex_id, 
										    mt_metade::id, pathfinder::id_fast,  pathfinder::id_spiral, pathfinder::id_landscape };
		//deterministic solvers, which are known to always give the very same result

	std::set<GUID> no_population_solvers = { nlopt::newuoa_id, nlopt::bobyqa_id, nlopt::simplex_id, nlopt::subplex_id, nlopt::praxis_id };

	std::map<GUID, uint64_t> faulty_solvers = { {nlopt::bobyqa_id, 100} };

	std::set<GUID> allowed_solvers = {nlopt::bobyqa_id, nlopt::newuoa_id, nlopt::praxis_id, nlopt::simplex_id, nlopt::subplex_id,
									  pagmo::abc_id, pagmo::cmaes_id, pagmo::sade_id, pagmo::de1220_id, pagmo::pso_id, pagmo::xnes_id, pagmo::ihs_id, pagmo::gpso_id,
									  halton_metade::id, mt_metade::id, rnd_metade::id, ppr::spo_id, 
									  pathfinder::id_fast, pathfinder::id_spiral, pathfinder::id_landscape,
									  //sequential_brute_force_scan::id, // disable for preliminary analysis (for full test, this should be enabled as a reference algorithm)
									  pso::id, rumoropt::id
									};

		//some solvers fail on the problem size=> we need to check it to a avoid forcefull cancellation of a long computation
}


void Run_Solver(const scgms::TSolver_Descriptor &desc, CCommon_Problem * working_problem, const size_t max_generations, const size_t population_size, TSolver_Result &result) {

		
	CSolution lower_bound, upper_bound;
	std::unique_ptr<CSolution> optimum = std::make_unique<CSolution>();	//when used for the semestral project, some people might have noticed that the optimum 
											//vector sits right after the upper_bound
	double optimum_fitness;
	working_problem->get_bounds(lower_bound, upper_bound);
	working_problem->get_optimum(*optimum, optimum_fitness);		

	bool failed = false;

	CSolution local_parameters;
	local_parameters.setConstant(std::numeric_limits<double>::quiet_NaN(), lower_bound.size());

	//using TObjective_Function = BOOL(IfaceCalling*)(const void* data, const size_t count, const double* solution, double* const fitness);

	solver::TSolver_Setup solver_setup{ lower_bound.size(), 1,
							lower_bound.data(), upper_bound.data(),
							nullptr, 0,			//no hints
							local_parameters.data(),
							working_problem, Fitness_Wrapper, nullptr,
							max_generations, population_size, std::numeric_limits<double>::min() };



	solver::TSolver_Progress solver_progress{ 0 };

	std::chrono::high_resolution_clock::time_point Solve_Start_Time = std::chrono::high_resolution_clock::now();
	try {
		if (solver::Solve_Generic(desc.id, solver_setup, solver_progress) != S_OK)
			failed = true;
	}
	catch (...) { failed = true; }
		
	std::chrono::high_resolution_clock::time_point Solve_Stop_Time = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> secs_duration = Solve_Stop_Time - Solve_Start_Time;
	result.seconds.push_back(secs_duration.count()*0.001);

	double total_calls, least_call, least_call_001;
	CSolution params_001;
	working_problem->Get_Objective_Calls(total_calls, least_call, least_call_001, params_001);
	result.total_objective_calls.push_back(total_calls);
	result.least_objective_call.push_back(least_call);
	result.least_objective_call_001.push_back(least_call_001);

	const double local_fitness = working_problem->Calculate_Fitness(local_parameters.data());
	if (isnan(local_fitness)) failed = true;

	result.optimum_fitness.push_back(optimum_fitness);
	result.fitness.push_back(local_fitness);
	result.fitness_error.push_back(fabs(local_fitness - optimum_fitness));
		
	for (size_t i = 0; i < local_parameters.size(); i++) {
		result.optimum[i].push_back((*optimum)[i]);

		result.parameters[i].push_back(local_parameters[i]);
		result.abs_parameter_error.push_back(fabs(local_parameters[i] - (*optimum)[i]));			
	}

	for (size_t i = 0; i < params_001.size(); i++) {
		result.abs_parameter_error_001.push_back(fabs(params_001[i] - (*optimum)[i]));
	}

	if (failed) result.fail_count++;

}

std::vector<TSolver_Result> Run_Solvers(size_t repetitions, CCommon_Problem *problem, bool randomize_optimum) {

	std::vector<TSolver_Result> results;

	

	//check, whether the problem can be actually solved
	if (!problem->Can_Be_Solved()) return results; //likely, the problem cannot be solved for this particular problem size, thus causing some algorithms to fail or run forever, such as Pagmo::ABC

	size_t Max_Generations = 100'000;
	std::vector<size_t> population_size = { 7, 15, 25, 40, 60, 100 };

	if (diagnostic::debugging) {
		//Max_Generations = 100'000;
		population_size = { 100 };
	}
	
	const auto solvers = scgms::get_solver_descriptor_list();
	auto working_problem = problem->Clone();


	for (const size_t current_population_size : population_size) {

		std::map<GUID, TSolver_Result> working_results;
		auto run_solver = [&working_problem, Max_Generations, current_population_size, randomize_optimum](const scgms::TSolver_Descriptor& solver, TSolver_Result &result) {
			if (diagnostic::debugging) {


				if (diagnostic::allowed_solvers.find(solver.id) == diagnostic::allowed_solvers.end()) return;

				//if (wcsstr(solver.description, L"Pathfinder") == nullptr) return;
				//if (wcscmp(solver.description, L"Pathfinder")) return;
			}

			
			working_problem->reset_counters();
			Run_Solver(solver, working_problem.get(), Max_Generations, current_population_size, result);
		};


		//initialize the results
		for (const auto& solver : solvers) {

			TSolver_Result result;
			result.optimum.resize(problem->Problem_Size());
			result.parameters.resize(problem->Problem_Size());
			result.name = solver.description;
			if (current_population_size > 0) {
				result.name += L"_";
				result.name += std::to_wstring(current_population_size);
			}
			result.fail_count = 0;
			
			working_results[solver.id] = result;
		}



		std::wcout << L"Executing repetition... " << std::endl;
		std::wcout << L"Repetition; Optimum parameter ";
		for (size_t j = 0; j < working_problem->Problem_Size(); j++)
			std::wcout << j << "; ";
		std::wcout << std::endl;

		for (size_t repetition = 0; repetition < repetitions; repetition++) {			
			if (randomize_optimum) working_problem->randomize_shift();	//we need to ensure that in each iteration each solver has exactly the same problem

			//print optimum parameters:			
			std::wcout << repetition << L"; ";
			CSolution optimum_params;
			double optimum_fitness;
			working_problem->get_optimum(optimum_params, optimum_fitness);
			for (size_t j = 0; j < working_problem->Problem_Size(); j++)
				std::wcout << optimum_params[j] << "; ";
			std::wcout << std::endl << std::flush;
			
			
			for (const auto& solver : solvers) {
				TSolver_Result &result = working_results[solver.id];

				bool faulty = solver.specialized;	//skip specilazed solvers as well
				//check if it is faulty solver
				if (!faulty) {
					const auto fs = diagnostic::faulty_solvers.find(solver.id);
					if (fs != diagnostic::faulty_solvers.end()) {
						if (problem->Problem_Size() >= fs->second) {
							result.fail_count++;
							faulty = true;
							break;
						}

					}
				}

				if (!faulty) run_solver(solver, result);
			}
		}
		std::wcout << std::endl;

		//put the results to the overall results
		for (auto& result : working_results)
			results.push_back(std::move(result.second));
	}
	


	return results;
}

void Evaluate_Solvers(CCommon_Problem *problem, const size_t repetitions, const size_t problem_ordinal_number, bool randomize_optimum) {

	const size_t problem_size = problem->Problem_Size();

	//write prolog
	{
		double fitness;
		CSolution optimum, ub;
		problem->get_bounds(optimum, ub);
		problem->get_optimum(optimum, fitness);
		std::cout << "--=== Evaluating solvers on " << problem->Get_Name() << " with problem size = " << problem_size << "... ===--" << std::endl;
		std::cout << "Problem ordinal number: " << problem_ordinal_number << std::endl;
		std::cout.precision(3);
		std::cout << "upper_bound[0]=" << ub[0] << ", optimum_analytical_parameter[0]=" << optimum[0] << ", optimum_fitness=" << fitness << std::endl << std::endl;
	}

	//1. run and collect results
	std::vector<TSolver_Result> results = Run_Solvers(repetitions, problem, randomize_optimum);

	if (!results.empty()) {
		//print the fitness and its optimium as avg +- stdev
		CStats global_optimum_fitness;
		std::vector<CStats> global_optimum{ problem_size };
		


		//2. evaluate the results 
		for (auto &stats : results) {
			for (auto &param : stats.parameters)
				param.Calculate_Stats();

			stats.abs_parameter_error.Calculate_Stats();
			stats.fitness.Calculate_Stats();

			stats.fitness_error.Calculate_Stats();
			stats.abs_parameter_error.Calculate_Stats();
			stats.abs_parameter_error_001.Calculate_Stats();

			stats.total_objective_calls.Calculate_Stats();
			stats.least_objective_call.Calculate_Stats();
			stats.least_objective_call_001.Calculate_Stats();

			stats.seconds.Calculate_Stats();

			if (randomize_optimum) {
				global_optimum_fitness.insert(global_optimum_fitness.end(), stats.optimum_fitness.begin(), stats.optimum_fitness.end());
				for (size_t i = 0; i < problem_size; i++) {
					global_optimum[i].insert(global_optimum[i].begin(), stats.optimum[i].begin(), stats.optimum[i].end());
				}
			}
		}

		if (randomize_optimum) {
			global_optimum_fitness.Calculate_Stats();
			for (auto& opt : global_optimum)
				opt.Calculate_Stats();
			std::cout << "randomized optimum_fitness=" << global_optimum_fitness.Get_Stats().avg << " +/- " << global_optimum_fitness.Get_Stats().stddev << std::endl;
			std::cout << std::endl;
		}



		//find the best algorithm by fitness error, avg param error and then by time
		std::sort(results.begin(), results.end(), [](const TSolver_Result &a, TSolver_Result &b) {
			int fails = a.fail_count - b.fail_count;
			if (fails != 0) return fails < 0; // a.fail_count <= b.fail_count;

			double diff = a.abs_parameter_error.Get_Stats().avg - b.abs_parameter_error.Get_Stats().avg;
			if (diff != 0.0) return diff < 0.0; //a.abs_parameter_error.Get_Stats().avg < b.abs_parameter_error.Get_Stats().avg;

			diff = a.fitness_error.Get_Stats().avg - b.fitness_error.Get_Stats().avg;
			if (diff != 0.0) return diff < 0.0; //a.fitness_error.Get_Stats().avg < b.fitness_error.Get_Stats().avg;

			diff = a.abs_parameter_error_001.Get_Stats().avg - b.abs_parameter_error_001.Get_Stats().avg;
			if (diff != 0.0) return diff < 0.0;

			diff = a.least_objective_call_001.Get_Stats().avg - b.least_objective_call_001.Get_Stats().avg;
			if (diff != 0.0) return diff < 0.0;

			return wcscmp(a.name.c_str(), b.name.c_str()) < 0;

			//return a.seconds.Get_Stats().avg < b.seconds.Get_Stats().avg; - comparing the seconds is no longer relevant due to the massively parallel execution
		});


		std::vector<CStats> parameters;
		CStats fitness;

		CStats fitness_error;
		CStats abs_parameter_error;

		CStats total_objective_calls;
		CStats least_objective_call;
		CStats least_objective_call_001;

		CStats seconds;


		//3. print the results
		std::cout << std::endl;
		std::string title_line = "general;;;;;;;;";
		std::string header_line = "solver; fails; param_err; fitness_err; least_calls; lc_001; param_err001; ";

		auto add_params = [&title_line, &header_line, problem_size](const char* title) {
			title_line += title;
			title_line += ";;;;;;;;;";
			header_line += "; fitness; fitness_err; param_err; time; total calls; least calls; lc_001; pe_001; ";
			for (size_t i = 0; i < problem_size; i++) {
				title_line += "; ";
				header_line += std::to_string(i);
				header_line += "; ";
			}
		};

		add_params("Average");
		add_params("Standard Deviation");
		add_params("Minimum");
		add_params("Q25");
		add_params("Median");
		add_params("Q75");
		add_params("Maximum");

		std::cout << title_line << std::endl;
		std::cout << header_line << std::endl;


		auto write_marker = [problem_size](auto getter, const TSolver_Result &result) {
			std::cout.precision(std::numeric_limits< double >::max_digits10); 
			std::cout << std::scientific;
			std::cout << getter(result.fitness) << "; ";
			std::cout << getter(result.fitness_error) << "; ";
			std::cout << getter(result.abs_parameter_error) << "; ";

			std::cout.precision(3);
			std::cout << getter(result.seconds) << "; ";

			std::cout.precision(std::numeric_limits< double >::max_digits10);
			std::cout << std::scientific;
			std::cout << getter(result.total_objective_calls) << "; ";
			std::cout << getter(result.least_objective_call) << "; ";
			std::cout << getter(result.least_objective_call_001) << "; ";
			std::cout << getter(result.abs_parameter_error_001);


			for (size_t i = 0; i < problem_size; i++)
				std::cout << "; " << getter(result.parameters[i]);

			std::cout << ";;";
		};

		auto write_avg_stddev = [](auto getter) {
			const auto &stats = getter().Get_Stats();
			std::cout << stats.avg << " ± " << stats.stddev << "; ";
		};


		for (size_t i = 0; i < results.size(); i++) {
			const auto &result = results[i];
			std::wcout << result.name << "; " << result.fail_count << "; ";
			std::cout.precision(3);
			std::cout << std::scientific;

			write_avg_stddev([&result]() {return result.abs_parameter_error; });
			write_avg_stddev([&result]() {return result.fitness_error; });
			write_avg_stddev([&result]() {return result.least_objective_call; });
			write_avg_stddev([&result]() {return result.least_objective_call_001; });
			write_avg_stddev([&result]() {return result.abs_parameter_error_001; });

			std::cout << ";";
			write_marker([](const CStats &stats) {return stats.Get_Stats().avg; }, result);
			write_marker([](const CStats &stats) {return stats.Get_Stats().stddev; }, result);
			write_marker([](const CStats &stats) {return stats.Get_Stats().min; }, result);
			write_marker([](const CStats &stats) {return stats.Get_Stats().q25; }, result);
			write_marker([](const CStats &stats) {return stats.Get_Stats().med; }, result);
			write_marker([](const CStats &stats) {return stats.Get_Stats().q75; }, result);
			write_marker([](const CStats &stats) {return stats.Get_Stats().max; }, result);
			std::cout << std::endl;
		}
	} else 
	  std::cout << "This problem cannot be solved with the chosen problem size.";

	std::cout << std::endl << "--=== " << problem->Get_Name() << " evaluation completed. ===--" << std::endl << std::endl;
}