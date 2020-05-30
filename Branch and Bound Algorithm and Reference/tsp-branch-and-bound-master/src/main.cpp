#include <iostream>

#include "include/Solver.hpp"
#include "include/Parser.hpp"

void runAndPrintResults(Parser& parser)
{
    Solver solver{parser.getCitiesMatrix()};
    solver.run();
    solver.printSolution();
}

int main(int argc, char** argv)
{
    char choice;
    std::string filename;
    do
    {
        std::cout << "1. Test using full matrix." << std::endl;
        std::cout << "2. Test using EUC_2D list." << std::endl;
        std::cout << "3. Test using Lower Diagonal Row matrix." << std::endl;
        std::cout << "0. Exit" << std::endl;
        std::cout << std::endl;
        std::cout.flush();
        std::cin >> choice;
        switch (choice)
        {
        case '1':
        {
            try
            {
                std::cout << "Filename: \n";
                std::cin >> filename;

                Parser parser(filename);
                parser.loadCitiesMatrix();

                runAndPrintResults(parser);
            }
            catch (const std::runtime_error& ex)
            {
                std::cout << ex.what() << std::endl;
            }
            break;
        }
        case '2':
        {
            try
            {
                std::cout << "Filename: \n";
                std::cin >> filename;

                Parser parser(filename);
                parser.loadCitiesList();

                runAndPrintResults(parser);
            }
            catch (const std::runtime_error& ex)
            {
                std::cout << ex.what() << std::endl;
            }
            break;
        }
        case '3':
        {
            try
            {
                std::cout << "Filename: \n";
                std::cin >> filename;

                Parser parser(filename);
                parser.loadLowerDiagonalRow();

                runAndPrintResults(parser);
            }
            catch (const std::runtime_error& ex)
            {
                std::cout << ex.what() << std::endl;
            }
            break;
        }
        default:
            break;
        }
    } while (choice != '0');

    return 0;
}
