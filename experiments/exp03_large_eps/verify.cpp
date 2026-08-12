#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
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
typedef kv::interval< kv::mpfr< 1500 > > DataType;
typedef AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;
typedef vcp::mats< AppData > POLICY;
typedef vcp::pidblas VPOLICY;

namespace fs = std::filesystem;

int main() {
    std::cout.precision(17);

    int p = 3;
    int Dimension = 1;
    int Number_of_variables = 1;
    double D = 0.06;
    double alpha = 0.25;
    double beta = 1.25;
    double gamma = 1.0;
    double a = 2.5;
    double b = 0.1;
    int target_order = 80;

    std::string in_dir = "80data";
    std::string out_base_dir = "results";

    if (!fs::exists(in_dir)) {
        std::cerr << "Error: Directory '" << in_dir << "' does not exist." << std::endl;
        return 1;
    }

    if (!fs::exists(out_base_dir)) {
        fs::create_directory(out_base_dir);
    }

    for (const auto& entry : fs::directory_iterator(in_dir)) {
        std::string filename = entry.path().filename().string();
        if (filename.find("_uh.matrix_kvdd") != std::string::npos) {
            std::string base_name = filename.substr(0, filename.find("_uh.matrix_kvdd"));
            
            // パラメータをパース
            std::stringstream ss(base_name);
            std::string token;
            std::vector<std::string> parts;
            while (std::getline(ss, token, '_')) {
                parts.push_back(token);
            }
            if (parts.size() < 8) continue;

            double epsilon = std::stod(parts[1]);
            double current_alpha = std::stod(parts[3]);
            alpha = current_alpha; // Update alpha based on file

            // Create result directory
            std::string result_dir = out_base_dir + "/" + base_name;
            if (!fs::exists(result_dir)) {
                fs::create_directory(result_dir);
            }

            // Output to text file
            std::string log_file = result_dir + "/verification.txt";
            std::ofstream out(log_file);
            out.precision(17);

            out << "==================================================" << std::endl;
            out << "Verification for: " << base_name << std::endl;
            out << "Epsilon = " << epsilon << ", Alpha = " << alpha << std::endl;
            out << "==================================================" << std::endl;

            vcp::matrix< AppData, POLICY > uh, vh, wh;
            vcp::load(uh, (in_dir + "/" + base_name + "_uh").c_str());
            vcp::load(vh, (in_dir + "/" + base_name + "_vh").c_str());
            vcp::load(wh, (in_dir + "/" + base_name + "_wh").c_str());

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
                x[d] = kv::interval< double >(0, 1.0);
            }
            std::vector< double > uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
            std::vector< double > uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

            double Linf_norm = std::max(std::abs(uh_min[0]), std::abs(uh_max[0]));
            out << "|| uh ||_H10  : " << H10_norm << std::endl;
            out << "|| uh ||_L2   : " << L2_norm << std::endl;
            out << "|| uh ||_Linf : " << Linf_norm << std::endl;
            
            // 出力用データ作成
            vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
            std::ofstream csv_out(result_dir + "/graph_data.csv");
            for (int i = 0; i < Grafics.rowsize(); ++i) {
                for (int j = 0; j < Grafics.columnsize(); ++j) {
                    csv_out << double(Grafics(i, j));
                    if (j < Grafics.columnsize() - 1) csv_out << ",";
                }
                csv_out << std::endl;
            }
            csv_out.close();

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
                    out << "Residual Norm : || Laplace(uh) - f(uh) ||_L2 <= " << Res << std::endl;
                }
                Res = CpFtilde * Res;
                out << "Residual Norm : || F(uh) ||_(H-1) <= " << Res << std::endl;
            }

            // ---------------------------------------------------------------------
            // G & Kantorovich
            // ---------------------------------------------------------------------
            VData G = VData(0);
            {
                G = 2*beta/pow(epsilon,2)*pow(Cs3,3) + 3*gamma/pow(epsilon,2)*pow(Cs4,4)*(2*uh_Ftilde_norm + 4*K*Res); 
            }

            VData Check = pow(K,2)*Res*G;
            out << "K^2RG " << Check << std::endl;
            if (Check.upper() <= 0.5 ){
                out << "Verification Succeed!" << std::endl;
                VData rho = (1 - sqrt(1 - 2*Check))/(K*G);
                out << "|| u*-uh || <= " << rho << std::endl;
            }
            else {
                out << "Verification failed..." << std::endl;
            }
            
            vcp::time.toc();
            out.close();
            std::cout << "Verified and logged: " << result_dir << std::endl;
        }
    }

    return 0;
}
