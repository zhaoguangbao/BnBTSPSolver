#include <string>
#include <ctime>
#include "Parser.hpp"
#include "bnbtspsolver.hpp"

using namespace std;

int main(int argc, char** argv)
{
    if(argc<2)
        return 0;
    const string filename = string(argv[1]);

    std::clock_t begin = clock();

    Parser<uint32_t> parser(filename);
    parser.loadCitiesMatrix(Parser<uint32_t>::convertType::tp_uint_t);
    //parser.loadAtspCitiesMatrix(Parser<uint32_t>::convertType::tp_uint_t);
    auto D=parser.getCitiesMatrix();
    size_t n=D.size();
    BnBTSPSol::branch_and_bound(n, D);

    std::clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "reading, initializing and computing the optimal tour took " << elapsed_secs << " s." << std::endl;
    getchar();
    return 0;
}
