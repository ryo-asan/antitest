// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

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
typedef kv::interval< kv::mpfr< 1500 > > DataType;
typedef AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;
typedef vcp::mats< AppData > POLICY;
typedef vcp::pidblas VPOLICY;

// タスク管理用の構造体
struct VerificationTask {
    int order;
    double epsilon;
    std::string uh_file;
    std::string vh_file;
    std::string wh_file;
};

// txtファイルからタスクを読み込む関数
std::vector<VerificationTask> read_tasks_from_file(const std::string& filename) {
    std::vector<VerificationTask> tasks;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return tasks;
    }

    while (std::getline(file, line)) {
        if (line.find("{") != std::string::npos && line.find("}") != std::string::npos && line.find("std::vector") == std::string::npos) {
            for (char& c : line) {
                if (c == '{' || c == '}' || c == ',' || c == '"') {
                    c = ' ';
                }
            }
            std::stringstream ss(line);
            VerificationTask task;
            if (ss >> task.order >> task.epsilon >> task.uh_file >> task.vh_file >> task.wh_file) {
                tasks.push_back(task);
            }
        }
    }
    return tasks;
}

// 文字列置換用のヘルパー関数
std::string replaceString(std::string subject, const std::string& search, const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

int main(void){
    std::cout.precision(17);
    
    // unique_tasks_list_2.txt からタスクを読み込む
    std::vector<VerificationTask> tasks = read_tasks_from_file("unique_tasks_list.txt");
    
    if (tasks.empty()) {
        std::cerr << "タスクが見つかりませんでした。プログラムを終了します。" << std::endl;
        return 1;
    }

    // 共通のパラメータ設定
    int p = 3;
    int Dimension = 1;
    int Number_of_variables = 1;
    double D = 0.06;
    double alpha = 0.25;
    double beta = 1.25;
    double gamma = 1.0;
    double a = 2.5;
    double b = 0.1;

    // 次数の設定
    int initial_order = 40;
    int target_order = 80;
    int step_order = 10;

    // =========================================================================
    // Phase 1: 近似解の探索フェーズ (全タスクに対して次数を上げて保存)
    // =========================================================================
    std::cout << "##################################################" << std::endl;
    std::cout << " Phase 1: Approximation (Order " << initial_order << " -> " << target_order << ")" << std::endl;
    std::cout << "##################################################\n" << std::endl;

    for (size_t task_idx = 0; task_idx < tasks.size(); ++task_idx) {
        const auto& task = tasks[task_idx];
        double epsilon = task.epsilon;

        std::cout << "[Phase 1] Task " << task_idx + 1 << "/" << tasks.size() 
                  << " : Epsilon = " << epsilon << " | File = " << task.uh_file << std::endl;

        vcp::matrix< AppData, POLICY > uh, vh, wh;
        
        // 初期の次数40のデータをロード
        vcp::load(uh, ("data/" + task.uh_file + "_uh").c_str());
        vcp::load(vh, ("data/" + task.vh_file + "_vh").c_str());
        vcp::load(wh, ("data/" + task.wh_file + "_wh").c_str());

        // ゼロ解への飛び移りを防ぐため、次数を段階的に上げる
        for (int current_order = initial_order + step_order; current_order <= target_order; current_order += step_order) {
            
            vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
            Approximate_Generator.setting(current_order, p, Dimension, Number_of_variables, 50);
            Approximate_Generator.setting_list(); // フルリストで統一
            
            vcp::matrix< int > list_uh = Approximate_Generator.output_list();

            // 新しい次数のサイズに合わせてリサイズ（拡張分はゼロ埋めされる）
            uh.resize(list_uh.rowsize(), Number_of_variables);
            vh.resize(list_uh.rowsize(), Number_of_variables);
            wh.resize(list_uh.rowsize(), Number_of_variables);

            // ニュートン法による近似解の計算
            vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
            vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();
            vcp::matrix< AppData, POLICY > uhphi, uh2phi, uh3phi, uhphiphi, uh2phiphi;
            vcp::matrix< AppData, POLICY > DF, F, syuusei, check;

            AppData cc;
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
                
                if (cc < pow(2.0, -30)) {
                    vh = lss((DL + L), L) * uh;
                    wh = lss((D * D * DL + L), L) * uh;
                    break;
                }
            }
        } // 段階的な次数の引き上げ終了

        // 目標の次数(80)に到達したら保存
        std::string new_uh_file = replaceString(task.uh_file, "order40", "order80");
        std::string new_vh_file = replaceString(task.vh_file, "order40", "order80");
        std::string new_wh_file = replaceString(task.wh_file, "order40", "order80");

        vcp::save(uh, ("data_80/" + new_uh_file + "_uh").c_str());
        vcp::save(vh, ("data_80/" + new_vh_file + "_vh").c_str());
        vcp::save(wh, ("data_80/" + new_wh_file + "_wh").c_str());
        
        std::cout << "  -> Saved as: data_80/" << new_uh_file << std::endl;
    }

    // =========================================================================
    // Phase 2: 精度保証フェーズ (保存された次数80の解をロードして検証)
    // =========================================================================
    std::cout << "\n##################################################" << std::endl;
    std::cout << " Phase 2: Rigorous Verification (Order " << target_order << ")" << std::endl;
    std::cout << "##################################################\n" << std::endl;

    for (size_t task_idx = 0; task_idx < tasks.size(); ++task_idx) {
        const auto& task = tasks[task_idx];
        double epsilon = task.epsilon;
        
        std::string new_uh_file = replaceString(task.uh_file, "order40", "order80");
        std::string new_vh_file = replaceString(task.vh_file, "order40", "order80");
        std::string new_wh_file = replaceString(task.wh_file, "order40", "order80");

        std::cout << "\n==================================================" << std::endl;
        std::cout << "[Phase 2] Task " << task_idx + 1 << "/" << tasks.size() << " : Epsilon = " << epsilon << std::endl;
        std::cout << "File = data_80/" << new_uh_file << std::endl;
        std::cout << "==================================================" << std::endl;

        vcp::matrix< AppData, POLICY > uh, vh, wh;
        
        // Phase 1で作成した次数80のデータをロード
        vcp::load(uh, ("data_80/" + new_uh_file + "_uh").c_str());
        vcp::load(vh, ("data_80/" + new_vh_file + "_vh").c_str());
        vcp::load(wh, ("data_80/" + new_wh_file + "_wh").c_str());

        vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
        Approximate_Generator.setting(target_order, p, Dimension, Number_of_variables, 50);
        Approximate_Generator.setting_list();
        vcp::matrix< int > list_uh = Approximate_Generator.output_list();
        Approximate_Generator.setting_uh(uh);

        vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
        vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();
        AppData H10_norm = sqrt(transpose(uh) * DL * uh)(0);
        AppData L2_norm  = sqrt(transpose(uh) * L * uh)(0);

        std::vector< kv::interval< double > > x(Dimension);
        for (int d = 0; d < Dimension; d++) {
            x[d] = kv::interval< double >(0, 0.5);
        }
        std::vector< double > uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
        std::vector< double > uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

        double Linf_norm = std::max(std::abs(uh_min[0]), std::abs(uh_max[0]));
        std::cout << "|| uh ||_H10  : " << H10_norm << std::endl;
        std::cout << "|| uh ||_L2   : " << L2_norm << std::endl;
        std::cout << "|| uh ||_Linf : " << Linf_norm << std::endl;

        Approximate_Generator.clear();

        // ---------------------------------------------------------------------
        // Eigenvalue of tilde{F}
        // ---------------------------------------------------------------------
        vcp::time.tic();
        VData K = VData(0);
        VData CN, Cs3, Cs4, CpFtilde;
        VData uh_infsup;
        VData uh_Ftilde_norm;
        uh_infsup.lower() = uh_min[0];
        uh_infsup.upper() = uh_max[0];
        
        {
            vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
            Verification_Generator.setting(target_order, p, Dimension, Number_of_variables, 1, target_order);
            Verification_Generator.setting_list();
            
            vcp::matrix< VData, VPOLICY > uhi;
            vcp::convert(uh, uhi);
            Verification_Generator.setting_uh(uhi, list_uh, 1);

            vcp::matrix< VData, VPOLICY > DL_v = Verification_Generator.dphidphi();
            vcp::matrix< VData, VPOLICY > uh2phiphi = Verification_Generator.uhphiphi(2);
            vcp::matrix< VData, VPOLICY > L_v = Verification_Generator.phiphi();
            
            vcp::matrix< VData, VPOLICY > sigma_tmp =(VData(alpha)*L_v + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L_v*lss(DL_v+L_v,L_v)+b*L_v*lss(D*D*DL_v+L_v,L_v)))/pow(epsilon,2);
            VData sigma = max(max(sigma_tmp))(0);
            
            vcp::matrix< VData, VPOLICY > G = DL_v + (VData(alpha)*L_v + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L_v*lss(DL_v+L_v,L_v)+b*L_v*lss(D*D*DL_v+L_v,L_v)))/pow(epsilon,2)+sigma;
            DL_v.clear();
            uh_Ftilde_norm = (transpose(uhi)*G*uhi)(0);
        }

        {
            vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
            Verification_Generator.setting(target_order, p, Dimension, Number_of_variables, 1, target_order);
            Verification_Generator.setting_list();
            
            vcp::matrix< VData, VPOLICY > uhi;
            vcp::convert(uh, uhi);
            Verification_Generator.setting_uh(uhi, list_uh, 1);

            vcp::matrix< VData, VPOLICY > DL_v = Verification_Generator.dphidphi();
            vcp::matrix< VData, VPOLICY > uhphiphi = Verification_Generator.uhphiphi(1);
            vcp::matrix< VData, VPOLICY > uh2phiphi = Verification_Generator.uhphiphi(2);
            vcp::matrix< VData, VPOLICY > L_v = Verification_Generator.phiphi();

            vcp::matrix< VData, VPOLICY > sigma_tmp =(VData(alpha)*L_v + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L_v*lss(DL_v+L_v,L_v)+b*L_v*lss(D*D*DL_v+L_v,L_v)))/pow(epsilon,2);
            VData sigma = max(max(sigma_tmp))(0);
            
            vcp::matrix< VData, VPOLICY > fdtildeuh_LL2 = (alpha*L_v + 3*gamma*uh2phiphi +epsilon*(a*L_v*lss(DL_v+L_v,L_v)+b*L_v*lss(D*D*DL_v+L_v,L_v)))/pow(epsilon,2)+sigma;
            vcp::matrix< VData, VPOLICY > G = DL_v + (VData(alpha)*L_v + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L_v*lss(DL_v+L_v,L_v)+b*L_v*lss(D*D*DL_v+L_v,L_v)))/pow(epsilon,2)+sigma;
            DL_v.clear();
            uh2phiphi.clear();
            
            vcp::matrix< VData, VPOLICY > E;
            compsym(G);
            eigsymge(L_v, G, E);
            L_v.clear();
            E = 1/diag(E);
            
            CN = Verification_Generator.Ritz_projection_error< VData >();
            Cs3 = Verification_Generator.Sobolev_constant< VData >(3);
            Cs4 = Verification_Generator.Sobolev_constant< VData >(4);
            
            VData CFtilde = CN*(1 + 1/min(E)(0)* max(max(fdtildeuh_LL2))(0));
            VData lambda_Ftilde = min(E)(0)/(1+pow(CFtilde,2)*min(E)(0));

            CpFtilde = 1/sqrt(lambda_Ftilde);
            CpFtilde.lower() = CpFtilde.upper();
            
            eigsymge( uhphiphi *VData(2.0)*beta/(VData(epsilon)*VData(epsilon)),G, E);
            
            G.clear();
            uhphiphi.clear();
            E = diag(E);
            
            vcp::matrix< VData, VPOLICY > O;
            int n = 0;
            int place = 0;
            for(int i=0 ;i < E.rowsize();i++){
                if(E(i).lower() > VData(1.0).upper()){
                    place=i;
                    n=1;
                    break;
                }
            }
            if( n == 1){
                O.zeros(2,1);
                O(0).lower()=E(place-1).lower();
                O(0).upper()=E(place-1).upper();
                O(1).lower()=E(place).lower();
                O(1).upper()=E(place).upper();
            }
            else {
                O.zeros(1,1);
                E = max(abs( E )); 
                O(0).lower()=E(0).lower();
                O(0).upper()=E(0).upper();
            }
            
            O = 1/O;
            VData CFtilde2 = sqrt(abs(2*VData(beta)/pow(epsilon,2)*uh_infsup+sigma))*CFtilde;
            CFtilde2.lower() = CFtilde2.upper();
            
            for (int i = 0; i < O.rowsize(); i++){
                O(i).lower() = (O(i)/(1+pow(CFtilde2,2)*O(i))).lower();
            }
            for (int i = 0; i < O.rowsize(); i++){
                if ( O(i).upper() < (VData(2.0)*beta/pow(epsilon,2)).lower() ){
                    O(i).lower() = O(i).upper();
                }
                else if ( O(i).lower() > (VData(2.0)*beta/pow(epsilon,2)).upper() ){
                    O(i).upper() = O(i).lower();
                }
            }
            for (int i = 0; i < O.rowsize(); i++ ){
                O(i) = 1/(1-1/O(i));
            }
            O = max(abs( O ));
            K = O(0).upper();
        }

        // ---------------------------------------------------------------------
        // Calculate Residual Norm || Laplace(uh) - f(uh) ||_L2
        // ---------------------------------------------------------------------
        VData Res = VData(0);
        {
            int Local_Number_of_variables = 2;
            vcp::Legendre_Bases_Generator< DataType, VResData, VResPOLICY > Verification_Generator;
            Verification_Generator.setting(target_order, p, Dimension, Local_Number_of_variables, 2);
            Verification_Generator.setting_list();
            
            vcp::matrix< VResData, VResPOLICY > uhi, vhi, whi, zhi, zhi2, zhi3;
            vcp::convert(uh, uhi);
            vcp::convert(vh, vhi);
            vcp::convert(wh, whi);

            zhi = horzcat(uhi, vhi);
            Verification_Generator.setting_uh(zhi);

            VResData uh2 = Verification_Generator.integral_uh(2,0);
            VResData uh3 = Verification_Generator.integral_uh(3,0);        
            VResData uh4 = Verification_Generator.integral_uh(4,0);
            VResData uh5 = Verification_Generator.integral_uh(5,0);
            VResData uh6 = Verification_Generator.integral_uh(6,0);
            VResData vh2 = Verification_Generator.integral_uh(0, 2);
            VResData uhvh = Verification_Generator.integral_uh(1, 1); 
            VResData uh2vh = Verification_Generator.integral_uh(2, 1); 
            VResData uh3vh = Verification_Generator.integral_uh(3, 1); 
                                           
            VResData LuhLuh = Verification_Generator.integral_LuhLuh(0);
            VResData Luh_uh1 = Verification_Generator.integral_Luhuh(0, 1, 0);
            VResData Luh_uh2 = Verification_Generator.integral_Luhuh(0, 2, 0);
            VResData Luh_uh3 = Verification_Generator.integral_Luhuh(0, 3, 0); 
            VResData Luh_vh = Verification_Generator.integral_Luhuh(0, 0, 1); 
            VResData Lvh_Lvh = Verification_Generator.integral_LuhLuh(1); 
            VResData Lvh_uh = Verification_Generator.integral_Luhuh(1, 1, 0); 
            VResData Lvh_vh = Verification_Generator.integral_Luhuh(1, 0, 1); 

            VResData firstuv = LuhLuh;
            VResData seconduv = -VResData(2)*(VResData(alpha)*Luh_uh1-VResData(beta)*Luh_uh2+VResData(gamma)*Luh_uh3+VResData(epsilon)*VResData(a)*Luh_vh)/pow(epsilon,2);
            VResData thirduv = 2*alpha*a*uhvh/pow(epsilon,3)+pow(alpha,2)*uh2/pow(epsilon,4)-2*a*beta*uh2vh/pow(epsilon,3)+2*a*gamma*uh3vh/pow(epsilon,3)-2*alpha*beta*uh3/pow(epsilon,4)+(2*alpha*gamma+pow(beta,2))*uh4/pow(epsilon,4)-2*beta*gamma*uh5/pow(epsilon,4)+pow(gamma,2)*uh6/pow(epsilon,4);
            
            VData CRonetmp = abs(a/epsilon)*CN*sqrt(1+pow(CN,2));
            VResData CRone;
            vcp::convert(CRonetmp,  CRone);

            VResData fourth = uh2 + 2*Lvh_uh - 2*uhvh + Lvh_Lvh - 2*Lvh_vh + vh2;

            zhi2 = horzcat(uhi, whi);
            Verification_Generator.setting_uh(zhi2);

            VResData wh2 = Verification_Generator.integral_uh(0, 2);
            VResData uhwh = Verification_Generator.integral_uh(1, 1); 
            VResData uh2wh = Verification_Generator.integral_uh(2, 1); 
            VResData uh3wh = Verification_Generator.integral_uh(3, 1); 

            VResData Luh_wh = Verification_Generator.integral_Luhuh(0, 0, 1); 
            VResData Lwh_Lwh = Verification_Generator.integral_LuhLuh(1); 
            VResData Lwh_uh = Verification_Generator.integral_Luhuh(1, 1, 0); 
            VResData Lwh_wh = Verification_Generator.integral_Luhuh(1, 0, 1); 

            VResData seconduw = -2*(VResData(b)*Luh_wh)/epsilon;
            VResData thirduw = (2*alpha*b*uhwh-2*b*beta*uh2wh+2*b*gamma*uh3wh)/pow(epsilon,3);
            
            VData CRtwotmp = abs(b/epsilon)*CN*sqrt(1+pow(CN,2)/pow(D,2));
            VResData CRtwo;
            vcp::convert(CRtwotmp,  CRtwo);
            
            VResData fifth =  uh2/pow(D,4) + 2*Lwh_uh/pow(D,2) - 2*uhwh/pow(D,4) + Lwh_Lwh - 2*Lwh_wh/pow(D,2) + wh2/pow(D,4);
            
            zhi3 = horzcat(vhi, whi);
            Verification_Generator.setting_uh(zhi3);
            VResData vh2_vw = Verification_Generator.integral_uh(2, 0);
            VResData wh2_vw = Verification_Generator.integral_uh(0, 2);
            VResData vhwh = Verification_Generator.integral_uh(1, 1); 

            VResData thirdvw =(pow(a,2)*vh2_vw + 2*a*b*vhwh + pow(b,2)*wh2_vw)/pow(epsilon,2);

            VResData first = firstuv;
            VResData second = seconduv + seconduw;
            VResData third = thirduv + thirduw + thirdvw;

            {
                using std::sqrt;
                using std::abs;
                vcp::convert(sqrt(abs(first + second + third)) + CRone*sqrt(abs(fourth)) + CRtwo*sqrt(abs(fifth)), Res);
                std::cout << "Residual Norm : || Laplace(uh) - f(uh) ||_L2 <= " << Res << std::endl;
            }
            Res = CpFtilde * Res;
            std::cout << "Residual Norm : || F(uh) ||_(H-1) <= " << Res << std::endl;
        }

        // ---------------------------------------------------------------------
        // G & Kantorovich
        // ---------------------------------------------------------------------
        VData G = VData(0);
        {
            G = 2*beta/pow(epsilon,2)*pow(Cs3,3) + 3*gamma/pow(epsilon,2)*pow(Cs4,4)*(2*uh_Ftilde_norm + 4*K*Res); 
        }

        VData Check = pow(K,2)*Res*G;
        std::cout << "K^2RG " << Check << std::endl;
        if (Check.upper() <= 0.5 ){
            std::cout << "Verification Succeed!" << std::endl;
            VData rho = (1 - sqrt(1 - 2*Check))/(K*G);
            std::cout << "|| u*-uh || <= " << rho << std::endl;
        }
        else {
            std::cout << "Verification failed..." << std::endl;
        }

        vcp::time.toc();

    } // Phase 2 タスクループ終了

    return 0;
}