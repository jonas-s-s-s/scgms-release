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

#include <iostream>

#ifdef WIN32
	#include <Windows.h>
#endif

int __cdecl main(int argc, char* argv[]) {
#ifdef WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	std::cout << "Welcome to the test of the solvers against the Pathfinder." << std::endl << std::endl;

	size_t problem_size = 3;
	size_t repetitions = 1;

	if (argc > 1 && isdigit(argv[1][0])) {
		problem_size = std::atoi(argv[1]);
	}
	else
		std::cout << "Usage: problem_size [repetitions] [problem_ordinal_number] [-randomize]" << std::endl << std::endl;

	if (argc > 2 && isdigit(argv[2][0])) {
		repetitions = std::atoi(argv[2]);
	}


	bool randomize_optimum = false;
	for (size_t i = 1; i < argc; i++) {

		if (strcmp(argv[i], "-randomize") == 0) {

			randomize_optimum = true;
			std::cout << "Will randomize optimum solutions." << std::endl;
			break;
		}
	}

	const auto problems = Create_Problem_Collection(problem_size);
	size_t low_problem_number = 0;
	size_t high_problem_number = problems.size()-1;
	if (argc > 3 && isdigit(argv[3][0])) {

		const size_t user_problem_number = std::atoi(argv[3]);
		if ((user_problem_number >= low_problem_number) && (user_problem_number <= high_problem_number))
			low_problem_number = high_problem_number = user_problem_number;
		else
			std::cout << "problem_ordinal_number out of bounds, ignoring it..." << std::endl;
	}
	
	for (size_t problem_number = low_problem_number; problem_number <= high_problem_number; problem_number++) {
		Evaluate_Solvers(problems[problem_number].get(), repetitions, problem_number, randomize_optimum);
	}

	return 0;
}