// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and / or other materials provided with the distribution.
// * Neither the name of the Kouta Sekine nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL KOUTA SEKINE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <sstream>      // std::stringstream のために必要
#include <iomanip>      // std::setw, std::setfill, std::fixed, std::setprecision のために必要
#include <string>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>

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

// typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;
// typedef kv::dd AppData_ddtodouble;

typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 1500 > > DataType;
typedef  AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;

//typedef vcp::pdblas POLICY;
typedef vcp::mats< AppData > POLICY;

typedef vcp::pidblas VPOLICY;



int main(void){

    // std::cout.precision(17);
    
    

    vcp::matrix< AppData, POLICY > uh;
    vcp::matrix< AppData, POLICY > vh;
    vcp::matrix< AppData, POLICY > wh;
    vcp::matrix< AppData, POLICY > uh_load;
    vcp::matrix< AppData, POLICY > vh_load;
    vcp::matrix< AppData, POLICY > wh_load;
    vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
    
    int Order_legendre = 40;
    int uh_Order_legendre = 40;
    int p = 3;
    int Dimension = 1;
    int Number_of_variables = 1;

    int count = 0;

    double D=0.06;
    double alpha=0.25;
    double beta=1.25;
    double gamma=1;
    double epsilon=0.025;
    double H10=0;
    double L2=0;
    double a=2.5;
    double b=0.1;
    // double N=1.29061;
    double N=1.289;
    // double ch =1;
    // double cha =1;
    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "gamma = " << gamma << std::endl;
    //std::cout << "delta = " << delta << std::endl;
    std::cout << "epsilon = " << epsilon << std::endl;
    
    // Setting of Approximate_Generator
    
    Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 50);
    // Setting the list of Approximate_Generator
    //Approximate_Generator.setting_list();
    Approximate_Generator.setting_evenlist();

    // output the list => list_uh
    vcp::matrix< int > list_uh = Approximate_Generator.output_list();
    
        
    // setting initialization value of uh
    vcp::load(uh, "value/test2_uh_eps_0.025000_count_01_uh");
    vcp::load(vh, "value/test2_uh_eps_0.025000_count_01_vh");
    vcp::load(wh, "value/test2_uh_eps_0.025000_count_01_wh");

        {
        // Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
        vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
        std::cout << "DL" << std::endl;
        std::cout << DL << std::endl;
        
        N = (double)sqrt(transpose(uh)*DL*uh)(0);
        std::cout << "Initial N set to: " << N << std::endl;
        // Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
        vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();
        std::cout << "L" << std::endl;
        std::cout << L << std::endl;
        vcp::matrix< AppData, POLICY > uhphi;
        vcp::matrix< AppData, POLICY > uh2phi;//追加
        vcp::matrix< AppData, POLICY > uh3phi;
        vcp::matrix< AppData, POLICY > uhphiphi;//追加
        vcp::matrix< AppData, POLICY > uh2phiphi;
        vcp::matrix< AppData, POLICY > DF;
        vcp::matrix< AppData, POLICY > DF_u;
        vcp::matrix< AppData, POLICY > F;
        vcp::matrix< AppData, POLICY > F_U;
        vcp::matrix< AppData, POLICY > F_S;        
        vcp::matrix< AppData, POLICY > F_epsolon;
        vcp::matrix< AppData, POLICY > G_u;
        vcp::matrix< AppData, POLICY > G_epsilon;
        vcp::matrix< AppData, POLICY > Y_U;
        vcp::matrix< AppData, POLICY > Y_S;        
        vcp::matrix< AppData, POLICY > Y;        
        vcp::matrix< AppData, POLICY > syuusei_extend;
        vcp::matrix< AppData, POLICY > syuusei_u;
        vcp::matrix< AppData, POLICY > syuusei_epsilon;
        vcp::matrix< AppData, POLICY > check;
        vcp::matrix< AppData, POLICY > A;
        

        double dN = -0.002;
        int max_steps = 50;
        
        for (int step = 1; step <= max_steps; step++) {
            vcp::matrix< AppData, POLICY > uh_backup = uh;
            double epsilon_backup = epsilon;
            double N_backup = N;
            
            N += dN;
            std::cout << "=== Step " << step << " / N = " << N << " ===" << std::endl;

            AppData cc;
            int iter = 0;
            bool converged = false;
            while(1){
                iter++;
                Approximate_Generator.setting_uh(uh);
                uhphi = Approximate_Generator.uhphi(1);
                uh2phi = Approximate_Generator.uhphi(2);
                uh3phi = Approximate_Generator.uhphi(3);
                uhphiphi = Approximate_Generator.uhphiphi(1);
                uh2phiphi = Approximate_Generator.uhphiphi(2);

                using std::pow;
                
                DF_u = DL - (-alpha *L +2*beta * uhphiphi -3 * gamma * uh2phiphi-epsilon*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)) )/pow(epsilon,2);
                F_epsolon = (2*(-alpha *uhphi +beta * uh2phi -gamma * uh3phi)-epsilon*(a*L*lss(DL+L,uhphi)+b*L*lss(D*D*DL+L,uhphi)) )/(epsilon*epsilon*epsilon);
                G_u=transpose(uh)*DL/sqrt(transpose(uh)*DL*uh)(0);
                G_epsilon.zeros(1,1);
                
                Y_U=horzcat(DF_u,F_epsolon);
                Y_S=horzcat(G_u,G_epsilon);
                Y=vercat(Y_U,Y_S);
                
                F_U = DL * uh - (-alpha *uhphi +beta * uh2phi -gamma * uh3phi-epsilon*(a*L*lss(DL+L,uhphi)+b*L*lss(D*D*DL+L,uhphi)) )/pow(epsilon,2);
                F_S =sqrt(transpose(uh)*DL*uh)-N;
                F=vercat(F_U,F_S);
                
                std::cout << "Y" << std::endl;
                std::cout << Y << std::endl;
                std::cout << "F" << std::endl;
                std::cout << F << std::endl;
                
                syuusei_extend = lss(Y, F);
                
                std::cout << "syuusei_extend" << std::endl;
                std::cout << syuusei_extend << std::endl;
                
                syuusei_u.zeros(uh.rowsize(),1);
                for(int i = 0; i < syuusei_extend.rowsize()-1; i++){
                    syuusei_u(i) = syuusei_extend(i);
                }
                
                uh = uh - syuusei_u;
                epsilon = epsilon - (double)syuusei_extend(syuusei_extend.rowsize()-1);
                
                check = max(abs(syuusei_u));
                cc = check(0);
                
                if (cc < pow(2.0,-30)) {
                    converged = true;
                    break;
                }
                if (iter > 50 || std::isnan((double)cc) || std::isinf((double)cc)) {
                    break; // failed
                }
            } // end of while(1) Newton loop
            
            if (converged) {
                double diff_u = (double)max(abs(uh - uh_backup))(0);
                double diff_eps = std::abs(epsilon - epsilon_backup);
                
                if (diff_u > 0.3 || diff_eps > 0.05) {
                    std::cout << "Jump detected! diff_u = " << diff_u << ", diff_eps = " << diff_eps << std::endl;
                    uh = uh_backup;
                    epsilon = epsilon_backup;
                    N = N_backup;
                    dN /= 2.0;
                    std::cout << "Halving dN to " << dN << std::endl;
                    step--; // Retry
                    if (std::abs(dN) < 1e-8) {
                        std::cout << "dN too small. Stopping." << std::endl;
                        break;
                    }
                } else {
                    vh=lss((DL+L),L)*uh;
                    wh=lss((D*D*DL+L),L)*uh;
                    Approximate_Generator.setting_uh(uh);
                    H10=(double)sqrt(transpose(uh)*DL*uh)(0);
                    L2=(double)sqrt(transpose(uh)*L*uh)(0);
                    
                    std::stringstream ss_uh, ss_vh, ss_wh;
                    ss_uh << std::fixed << std::setprecision(6) << "value/eps_" << epsilon << "_H10_" << H10 << "_uh";
                    ss_vh << std::fixed << std::setprecision(6) << "value/eps_" << epsilon << "_H10_" << H10 << "_vh";
                    ss_wh << std::fixed << std::setprecision(6) << "value/eps_" << epsilon << "_H10_" << H10 << "_wh";
                    
                    vcp::save(uh, ss_uh.str().c_str());
                    vcp::save(vh, ss_vh.str().c_str());
                    vcp::save(wh, ss_wh.str().c_str());
                    
                    std::ofstream csv_file("output_txtdata/bifurcation_data.csv", std::ios::app);
                    if (csv_file.is_open()) {
                        csv_file << N << "," << epsilon << "," << H10 << "\n";
                        csv_file.close();
                    }
                    
                    std::cout << "Step " << step << " successful. N = " << N << ", epsilon = " << epsilon << ", H10 = " << H10 << std::endl;
                }
            } else {
                std::cout << "Newton failed. Halving dN." << std::endl;
                uh = uh_backup;
                epsilon = epsilon_backup;
                N = N_backup;
                dN /= 2.0;
                step--; // Retry
                if (std::abs(dN) < 1e-8) {
                    std::cout << "dN too small. Stopping." << std::endl;
                    break;
                }
            }
        }
    }
    //(-DL+r*L),I.eye(L.rowsize())
    // vcp::time.toc();
    // vcp::save(uh, "test_vernot_1029/test1_uh");
    // vcp::save(vh, "test_vernot_1029/test1_vh");
    // vcp::save(wh, "test_vernot_1029/test1_wh");
    
    // uh data for Grafics
    vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
    // std::cout << Grafics << std::endl;//実験する時にコメント外す
    // std::cout << uh << std::endl;//実験する時にコメント外す

    // minimal and maximum value of approximate solution uh

    // std::cout << "\nCalculate the maximum and minimum value" << std::endl;
    std::vector< kv::interval< double > > x;

    x.resize(Dimension);
    for (int d = 0; d < Dimension; d++) {
        x[d] = kv::interval< double >(0, 0.5);
    }
    std::vector< double > uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
    std::vector< double > uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

    //std::vector< double > vh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
    //std::vector< double > vh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

    //vcp::save(uh_min, "Data_Nagumo/list_uh_min");
    //vcp::save(uh_max, "Data_Nagumo/list_uh_max");

    // for (int i = 0; i < Number_of_variables; i++) {
    //     std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
    //     std::cout << "vh in [" << vh_min[i] << ", " << vh_max[i] << "]" << std::endl;
    // }
        // vcp::save(uh, "txtdata_norm/eps=0.255_uh");
        // vcp::save(vh, "txtdata_norm/eps=0.255_vh");
        // vcp::save(wh, "txtdata_norm/eps=0.255_wh");
    
    
        
        // std::cout << "anser set start" << std::endl;
        // std::cout << "count:"<< count << std::endl;
        std::cout << "|| uh ||_H10" << std::endl;
        std::cout <<  H10 << std::endl;
        std::cout << "|| uh ||_L2" << std::endl;
        std::cout <<  L2 << std::endl;
        std::cout << "epsilon:"  <<  std::endl;
        std::cout << epsilon << std::endl;
        std::cout << "Linf:" << std::endl;
        std::cout << uh_max[Number_of_variables-1] << std::endl;
        std::cout << "uh_folm" << std::endl;
        std::cout << Grafics << std::endl;//実験する時にコメント外す
        std::cout << "anser set end" << std::endl;
    Approximate_Generator.clear();
    vcp::time.toc();

    
    return 0;
}