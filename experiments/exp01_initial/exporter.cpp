#include <iostream>
#include <string>
#include <vector>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
#include <vcp/ldbase.hpp>

typedef kv::dd AppData;
typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 1500 > > DataType;
typedef AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;
typedef vcp::mats< AppData > POLICY;

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <uh_filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    vcp::matrix< AppData, POLICY > uh;
    vcp::load(uh, filename.c_str());

    vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
    Approximate_Generator.setting(40, 3, 1, 1, 50);
    Approximate_Generator.setting_list();
    Approximate_Generator.setting_uh(uh);

    vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);

    // Grafics is probably a matrix of size 100 x 2 (or similar) containing x and u(x).
    // Let's print it in a CSV format.
    for (int i = 0; i < Grafics.rowsize(); ++i) {
        for (int j = 0; j < Grafics.columnsize(); ++j) {
            std::cout << double(Grafics(i, j));
            if (j < Grafics.columnsize() - 1) std::cout << ",";
        }
        std::cout << std::endl;
    }

    return 0;
}
