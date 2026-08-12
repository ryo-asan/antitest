// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.

#include <iostream>
#include <cmath>
#include <algorithm> // std::max, std::abs を使用するために追加
#include <sstream>   // std::stringstream を使用するために追加
#include <iomanip>   // std::fixed, std::setprecision を使用するために追加

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
// 【詳細コメント】データ型の定義
// 計算の精度やアルゴリズムに応じたデータ型を定義します。
// 近似解の計算では丸め誤差を抑えるため kv::dd（倍精度拡張）を使用します。
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

    // 【詳細コメント】出力精度を高めに設定
    std::cout.precision(17);
    
    // =========================================================================
    // 【詳細コメント】変数・行列の初期化
    // uh, vh, wh はそれぞれ対象とする方程式の未知関数に対応する係数行列です。
    // =========================================================================
    vcp::matrix< AppData, POLICY > uh;
    vcp::matrix< AppData, POLICY > vh;
    vcp::matrix< AppData, POLICY > wh;
    
    // 【詳細コメント】Legendre多項式に基づく基底関数の生成器
    vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
    
    // 【詳細コメント】パラメータの定義
    // 方程式の性質を決定する各種パラメータ（次元、多項式の次数、係数など）を設定します。
    int uh_Order_legendre = 40; // ルジャンドル多項式の展開次数
    int p = 3;                  // 非線形項の次数などに対応
    int Dimension = 1;          // 空間次元
    int Number_of_variables = 1;// 変数の数

    double D = 0.06;
    double alpha = 0.25;
    double beta = 1.25;
    double gamma = 1.0;
    double epsilon = 0.045000;
    double a = 2.5;
    double b = 0.1;

    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "gamma = " << gamma << std::endl;
    std::cout << "epsilon = " << epsilon << std::endl;
    
    // =========================================================================
    // 【詳細コメント】近似解生成器のセットアップ
    // 展開次数や次元を指定し、使用する基底のリストを作成します。
    // =========================================================================
    Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 50);
    Approximate_Generator.setting_list();
    vcp::matrix< int > list_uh = Approximate_Generator.output_list();

    // =========================================================================
    // 【詳細コメント】初期近似解の設定 (ones)
    // ファイルから読み込む代わりに、すべての要素を 1 に初期化した行列を用います。
    // =========================================================================
    uh.ones(list_uh.rowsize(), Number_of_variables);
    vh.ones(list_uh.rowsize(), Number_of_variables);
    wh.ones(list_uh.rowsize(), Number_of_variables);

    // 各種ノルムを保持する変数の宣言
    AppData H10_norm;
    AppData L2_norm;

    vcp::time.tic(); // 計算時間の計測開始

    {
        // =========================================================================
        // 【詳細コメント】内積行列の作成
        // DL: (nabla phi_i, nabla phi_j)_L2 （剛性行列に相当）
        // L:  (phi_i, phi_j)_L2 （質量行列に相当）
        // =========================================================================
        vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
        vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

        // ニュートン法で用いる各種行列の宣言
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
        
        // =========================================================================
        // 【詳細コメント】ニュートン法の反復ループ
        // 残差が十分に小さくなるまで、解の更新（syuuseiの減算）を繰り返します。
        // =========================================================================
        while(1){
            // 現在の近似解 uh をジェネレータにセットし、非線形項の計算を準備
            Approximate_Generator.setting_uh(uh);
            
            // 各種高次項・微分項の計算
            uhphi = Approximate_Generator.uhphi(1);
            uh2phi = Approximate_Generator.uhphi(2);
            uh3phi = Approximate_Generator.uhphi(3);
            uhphiphi = Approximate_Generator.uhphiphi(1);
            uh2phiphi = Approximate_Generator.uhphiphi(2);

            using std::pow;

            // =========================================================================
            // 【詳細コメント】ヤコビ行列 (DF) と残差ベクトル (F) の構築
            // 対象となる微分方程式系に基づいてDFおよびFを組み立てます。
            // =========================================================================
            DF = DL - (-alpha * L + 2 * beta * uhphiphi - 3 * gamma * uh2phiphi - epsilon * (a * L * lss(DL + L, L) + b * L * lss(D * D * DL + L, L))) / pow(epsilon, 2);
            F = DL * uh - (-alpha * uhphi + beta * uh2phi - gamma * uh3phi - epsilon * (a * L * lss(DL + L, uhphi) + b * L * lss(D * D * DL + L, uhphi))) / pow(epsilon, 2);  
            
            // 線形方程式 DF * syuusei = F を解いて修正量 (syuusei) を求める
            syuusei = lss(DF, F);
            
            // 近似解の更新
            uh = uh - syuusei;
            
            // 修正量の最大絶対値（無限大ノルム）をチェック
            check = max(abs(syuusei));
            cc = check(0);
            std::cout << "Correction max: " << cc << std::endl;
            
            // =========================================================================
            // 【詳細コメント】収束判定
            // 修正量が十分に小さくなった場合（ここでは 2^-30 を閾値とする）、
            // 収束とみなしてループを抜け、必要なノルムを計算します。
            // =========================================================================
            if (cc < pow(2.0, -30)) {
                // 従属する変数 vh, wh を更新
                vh = lss((DL + L), L) * uh;
                wh = lss((D * D * DL + L), L) * uh;
                Approximate_Generator.setting_uh(uh);

                // H10ノルムとL2ノルムの計算
                H10_norm = sqrt(transpose(uh) * DL * uh)(0);
                L2_norm  = sqrt(transpose(uh) * L * uh)(0);

                std::cout << "Newton Method Converged.\n" << std::endl;
                break;
            }
        }
    }

    vcp::time.toc(); // 計算時間の計測終了

    // =========================================================================
    // 【詳細コメント】近似解の最大値・最小値の計算
    // 指定した区間（ここでは各次元 [0, 0.5]）での関数 uh の大域的な最小値・最大値を求めます。
    // =========================================================================
    std::cout << "\nCalculate the maximum and minimum value" << std::endl;
    std::vector< kv::interval< double > > x;
    x.resize(Dimension);
    for (int d = 0; d < Dimension; d++) {
        x[d] = kv::interval< double >(0, 0.5);
    }
    
    std::vector< double > uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
    std::vector< double > uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

    for (int i = 0; i < Number_of_variables; i++) {
        std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
    }

    // =========================================================================
    // 【詳細コメント】各種ノルムの出力
    // L-infinityノルムは先ほど求めた最大値・最小値の絶対値の大きい方を採用します。
    // =========================================================================
    double Linf_norm = std::max(std::abs(uh_min[0]), std::abs(uh_max[0]));
    std::cout << "|| uh ||_H10  : " << H10_norm << std::endl;
    std::cout << "|| uh ||_L2   : " << L2_norm << std::endl;
    std::cout << "|| uh ||_Linf : " << Linf_norm << std::endl;

    // =========================================================================
    // 【詳細コメント】計算された近似解の保存
    // 動的にファイル名を生成し、vcp::save を用いて uh, vh, wh を保存します。
    // H10_norm(0).upper() を double にキャストしてファイル名に使用します。
    // =========================================================================
    std::stringstream ss_uh, ss_vh, ss_wh;
    ss_uh << "value/eps_" << std::fixed << std::setprecision(6) << epsilon 
          << "_H10_" << std::setprecision(5) << (double)H10_norm << "_uh";
    ss_vh << "value/eps_" << std::fixed << std::setprecision(6) << epsilon 
          << "_H10_" << std::setprecision(5) << (double)H10_norm << "_vh";
    ss_wh << "value/eps_" << std::fixed << std::setprecision(6) << epsilon 
          << "_H10_" << std::setprecision(5) << (double)H10_norm << "_wh";

    vcp::save(uh, ss_uh.str().c_str());
    vcp::save(vh, ss_vh.str().c_str());
    vcp::save(wh, ss_wh.str().c_str());
    
    std::cout << "\nSaved approximated solutions to `value/` directory." << std::endl;

    // メモリの解放と後始末
    Approximate_Generator.clear();

    return 0;
}
