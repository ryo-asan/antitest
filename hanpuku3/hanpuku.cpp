// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.

#include <iostream>
#include <cmath>
#include <algorithm> // std::max, std::abs
#include <sstream>   // std::stringstream
#include <iomanip>   // std::fixed, std::setprecision
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
#include <vcp/vcp_timer.hpp>

// =========================================================================
// データ型の定義
// =========================================================================
typedef kv::dd AppData;
typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 1500 > > DataType; // 高精度計算用
typedef AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;
typedef vcp::mats< AppData > POLICY;
typedef vcp::pidblas VPOLICY;

int main(void){

    std::cout.precision(17);
    
    // 変数・行列の初期化
    vcp::matrix< AppData, POLICY > uh;
    vcp::matrix< AppData, POLICY > vh;
    vcp::matrix< AppData, POLICY > wh;
    
    vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
    
    // 共通のパラメータ
    int uh_Order_legendre = 40;
    int p = 3;
    int Dimension = 1;
    int Number_of_variables = 1;

    double D = 0.06;
    double beta = 1.25;
    double gamma = 1.0;
    
    // 近似解生成器のセットアップ
    Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 50);
    Approximate_Generator.setting_list();
    vcp::matrix< int > list_uh = Approximate_Generator.output_list();

    // =========================================================================
    // 【ループによるパラメータ探索】
    // veri.cppを参考に、複数のパラメータ（例: alpha, epsilon, a, b）を変化させて近似解を求める。
    // =========================================================================
    
    std::vector<double> epsilons;
    for(double e = 0.060; e >= 0.0459; e -= 0.001) {
        epsilons.push_back(e);
    }
    std::vector<double> alphas = {0.25};
    double a = 2.5;
    double b = 0.1;

    for (double alpha : alphas) {
        for (double epsilon : epsilons) {
            std::cout << "\n=======================================================" << std::endl;
            std::cout << "パラメータ: alpha = " << alpha << ", epsilon = " << epsilon 
                      << ", a = " << a << ", b = " << b << std::endl;
            std::cout << "=======================================================" << std::endl;

            // 同じパラメータで複数回（例：100回）試行し、異なる非自明解を探す
            for (int trial = 1; trial <= 100; ++trial) {
                std::cout << "\n--- Trial " << trial << " ---" << std::endl;
                // =========================================================================
                // 【初期近似解の設定】
                // ゼロ解への吸い込みを防ぐため、毎回乱数で初期化して新しい非自明解を探索します。
                // =========================================================================
                uh.rand(list_uh.rowsize(), Number_of_variables);
                vh.zeros(list_uh.rowsize(), Number_of_variables);
                wh.zeros(list_uh.rowsize(), Number_of_variables);
                
                // 今回の乱数初期値を後で比較・再現できるように保存する
                std::stringstream ss_init;
                ss_init << "data/eps_" << std::fixed << std::setprecision(4) << epsilon 
                        << "_alpha_" << std::setprecision(2) << alpha
                        << "_trial_" << trial << "_init_uh";
                vcp::save(uh, ss_init.str().c_str());

                AppData H10_norm;
                AppData L2_norm;
                bool converged = false;

                vcp::time.tic(); // 計算時間の計測開始

                {
                    // 内積行列の作成
                    vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
                    vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

                    vcp::matrix< AppData, POLICY > uhphi;
                    vcp::matrix< AppData, POLICY > uh2phi;
                    vcp::matrix< AppData, POLICY > uh3phi;
                    vcp::matrix< AppData, POLICY > uhphiphi;
                    vcp::matrix< AppData, POLICY > uh2phiphi;
                    vcp::matrix< AppData, POLICY > DF;
                    vcp::matrix< AppData, POLICY > F;
                    vcp::matrix< AppData, POLICY > syuusei;
                    vcp::matrix< AppData, POLICY > check;

                    AppData cc;
                    int iter_count = 0;
                    
                    // ニュートン法の反復ループ
                    while(1){
                        Approximate_Generator.setting_uh(uh);
                        
                        uhphi = Approximate_Generator.uhphi(1);
                        uh2phi = Approximate_Generator.uhphi(2);
                        uh3phi = Approximate_Generator.uhphi(3);
                        uhphiphi = Approximate_Generator.uhphiphi(1);
                        uh2phiphi = Approximate_Generator.uhphiphi(2);

                        using std::pow;

                        // ヤコビ行列 (DF) と残差ベクトル (F) の構築
                        DF = DL - (-alpha * L + 2 * beta * uhphiphi - 3 * gamma * uh2phiphi - epsilon * (a * L * lss(DL + L, L) + b * L * lss(D * D * DL + L, L))) / pow(epsilon, 2);
                        F = DL * uh - (-alpha * uhphi + beta * uh2phi - gamma * uh3phi - epsilon * (a * L * lss(DL + L, uhphi) + b * L * lss(D * D * DL + L, uhphi))) / pow(epsilon, 2);  
                        
                        syuusei = lss(DF, F);
                        uh = uh - syuusei;
                        
                        check = max(abs(syuusei));
                        cc = check(0);
                        iter_count++;
                        std::cout << "Iteration " << iter_count << " Correction max: " << cc << std::endl;
                        
                        // 収束判定
                        if (cc < pow(2.0, -30)) {
                            vh = lss((DL + L), L) * uh;
                            wh = lss((D * D * DL + L), L) * uh;
                            Approximate_Generator.setting_uh(uh);

                            H10_norm = sqrt(transpose(uh) * DL * uh)(0);
                            L2_norm  = sqrt(transpose(uh) * L * uh)(0);

                            std::cout << "Newton Method Converged." << std::endl;
                            converged = true;
                            break;
                        }
                        if (iter_count > 100) {
                            std::cout << "Newton Method Failed (Max iteration reached)." << std::endl;
                            break;
                        }
                    }
                }

                vcp::time.toc(); // 計算時間の計測終了

                if (converged) {
                    std::cout << "|| uh ||_H10  : " << H10_norm << std::endl;
                    std::cout << "|| uh ||_L2   : " << L2_norm << std::endl;

                    // =========================================================================
                    // 近似解の保存 (hanpukuディレクトリからの相対パスで data/ に保存)
                    // どのパラメータ・初期値による結果かがわかるようにファイル名を設定
                    // =========================================================================
                    std::stringstream ss_uh, ss_vh, ss_wh;
                    ss_uh << "data/eps_" << std::fixed << std::setprecision(4) << epsilon 
                          << "_alpha_" << std::setprecision(2) << alpha
                          << "_trial_" << trial
                          << "_H10_" << std::setprecision(5) << (double)H10_norm << "_uh";
                    ss_vh << "data/eps_" << std::fixed << std::setprecision(4) << epsilon 
                          << "_alpha_" << std::setprecision(2) << alpha
                          << "_trial_" << trial
                          << "_H10_" << std::setprecision(5) << (double)H10_norm << "_vh";
                    ss_wh << "data/eps_" << std::fixed << std::setprecision(4) << epsilon 
                          << "_alpha_" << std::setprecision(2) << alpha
                          << "_trial_" << trial
                          << "_H10_" << std::setprecision(5) << (double)H10_norm << "_wh";

                    // 保存処理
                    vcp::save(uh, ss_uh.str().c_str());
                    vcp::save(vh, ss_vh.str().c_str());
                    vcp::save(wh, ss_wh.str().c_str());
                    
                    std::cout << "Saved converged solution to: " << ss_uh.str() << std::endl;
                }
            } // trialループの終わり

        }
    }

    Approximate_Generator.clear();
    return 0;
}
