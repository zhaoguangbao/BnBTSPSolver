//
// Created by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//

/**
 * @file main.cpp
 *
 * @brief execution of the programm
 */
#include <iostream>
#include <cstring>
#include <ctime>
#include "../header/tsp.hpp"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "No parameters were given. Please give an --instance ./instance.tsp as an program argument" << std::endl;
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "--instance") != 0) {
        std::cerr << "First argument should be an instance of TSP. Execute like ./program --instance ./dir_to_instance.tsp [--solution ./dir_to_output.opt.tour]";
        return EXIT_FAILURE;
    }
    std::string file = argv[2];

    std::clock_t begin = clock();

    TSP::Instance<double,double> myTSP(file);
    myTSP.compute_optimal_tour();
    std::cout << myTSP.length() << std::endl;
    std::clock_t end = clock();
    if(argc > 4){
        myTSP.print_optimal_tour(argv[4]);
    }

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cerr << "reading, initializing and computing the optimal tour took " << elapsed_secs << " s." << std::endl;
    return EXIT_SUCCESS;
}

//Hello