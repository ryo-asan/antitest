#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <filesystem>

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
#include <vcp/vcp_timer.hpp>

typedef kv::dd AppData;
typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 10000 > > DataType;
typedef AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;
typedef vcp::mats< AppData > POLICY;

namespace fs = std::filesystem;

int main() {
    std::cout.precision(17);

    // 共通のパラメータ設定
    int p = 3;
    int Dimension = 1;
    int Number_of_variables = 1;
    double D = 0.06;
    double beta = 1.25;
    double gamma = 1.0;
    double a = 2.5;
    double b = 0.1;

    int initial_order = 40;
    int target_order = 80;
    int step_order = 10;

    std::string in_dir = "data";
    std::string out_dir = "80data";

    if (!fs::exists(in_dir)) {
        std::cerr << "Error: Directory '" << in_dir << "' does not exist." << std::endl;
        return 1;
    }

    if (!fs::exists(out_dir)) {
        fs::create_directory(out_dir);
    }

    // 重複を排除するためのマップ: key=(eps, alpha, rounded_H10) -> (filename_base)
    std::map<std::tuple<double, double, double>, std::string> unique_tasks;

    for (const auto& entry : fs::directory_iterator(in_dir)) {
        std::string filename = entry.path().filename().string();
        if (filename.find("_uh.matrix_kvdd") != std::string::npos) {
            std::string base_name = filename.substr(0, filename.find("_uh.matrix_kvdd"));
            
            // パラメータをパース
            // 形式: eps_0.0350_alpha_0.25_trial_1_H10_3.50103
            std::stringstream ss(base_name);
            std::string token;
            std::vector<std::string> parts;
            while (std::getline(ss, token, '_')) {
                parts.push_back(token);
            }

            if (parts.size() >= 8) {
                double eps = std::stod(parts[1]);
                double alpha = std::stod(parts[3]);
                double h10 = std::stod(parts[7]);

                if (h10 < 1e-3) continue; // 自明解を除外

                double rounded_h10 = std::round(h10 * 100.0) / 100.0;
                auto key = std::make_tuple(eps, alpha, rounded_h10);

                if (unique_tasks.find(key) == unique_tasks.end()) {
                    unique_tasks[key] = base_name;
                }
            }
        }
    }

    if (unique_tasks.empty()) {
        std::cout << "対象となる非自明な近似解データが見つかりませんでした。" << std::endl;
        return 0;
    }

    std::cout << "抽出された固有の非自明解数: " << unique_tasks.size() << std::endl;
    std::cout << "次数の拡張を開始します (" << initial_order << " -> " << target_order << ")...\n" << std::endl;

    int task_idx = 0;
    for (const auto& pair : unique_tasks) {
        task_idx++;
        double epsilon = std::get<0>(pair.first);
        double alpha = std::get<1>(pair.first);
        std::string base_name = pair.second;

        std::cout << "[Task " << task_idx << "/" << unique_tasks.size() << "] eps=" << epsilon << ", alpha=" << alpha << std::endl;
        std::cout << "Target File: " << base_name << std::endl;

        vcp::matrix< AppData, POLICY > uh, vh, wh;
        vcp::load(uh, (in_dir + "/" + base_name + "_uh").c_str());
        vcp::load(vh, (in_dir + "/" + base_name + "_vh").c_str());
        vcp::load(wh, (in_dir + "/" + base_name + "_wh").c_str());

        bool success = true;
        for (int current_order = initial_order + step_order; current_order <= target_order; current_order += step_order) {
            std::cout << "  -> Updating to order " << current_order << "..." << std::endl;

            vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
            Approximate_Generator.setting(current_order, p, Dimension, Number_of_variables, 50);
            Approximate_Generator.setting_list();
            
            vcp::matrix< int > list_uh = Approximate_Generator.output_list();
            uh.resize(list_uh.rowsize(), Number_of_variables);

            vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
            vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();
            vcp::matrix< AppData, POLICY > uhphi, uh2phi, uh3phi, uhphiphi, uh2phiphi;
            vcp::matrix< AppData, POLICY > DF, F, syuusei, check;

            AppData cc;
            int iter_count = 0;
            bool converged = false;
            while(1){
                Approximate_Generator.setting_uh(uh);
                uhphi = Approximate_Generator.uhphi(1);
                uh2phi = Approximate_Generator.uhphi(2);
                uh3phi = Approximate_Generator.uhphi(3);
                uhphiphi = Approximate_Generator.uhphiphi(1);
                uh2phiphi = Approximate_Generator.uhphiphi(2);

                using std::pow;
                DF = DL - (-alpha * L + 2 * beta * uhphiphi - 3 * gamma * uh2phiphi - epsilon * (a * L * lss(DL + L, L) + b * L * lss(D * D * DL + L, L))) / pow(epsilon, 2);
                F = DL * uh - (-alpha * uhphi + beta * uh2phi - gamma * uh3phi - epsilon * (a * L * lss(DL + L, uhphi) + b * L * lss(D * D * DL + L, uhphi))) / pow(epsilon, 2);  
                
                syuusei = lss(DF, F);
                uh = uh - syuusei;
                check = max(abs(syuusei));
                cc = check(0);
                iter_count++;
                
                if (cc < pow(2.0, -30)) {
                    vh = lss((DL + L), L) * uh;
                    wh = lss((D * D * DL + L), L) * uh;
                    converged = true;
                    break;
                }
                if (iter_count > 100) {
                    std::cout << "     Newton Method Failed at order " << current_order << std::endl;
                    success = false;
                    break;
                }
            }
            if (!success) break;
        }

        if (success) {
            vcp::save(uh, (out_dir + "/" + base_name + "_order80_uh").c_str());
            vcp::save(vh, (out_dir + "/" + base_name + "_order80_vh").c_str());
            vcp::save(wh, (out_dir + "/" + base_name + "_order80_wh").c_str());
            std::cout << "  => Successfully saved to " << out_dir << "/" << base_name << "_order80_*" << std::endl;
        }
    }

    return 0;
}
