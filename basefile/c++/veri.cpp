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


#include <iostream>
#include <cmath>
#include <algorithm> // std::max, std::abs を使用するために追加

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
typedef kv::interval< kv::mpfr< 1500 > > DataType;//高め調整
typedef  AppData ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;

//typedef vcp::pdblas POLICY;
typedef vcp::mats< AppData > POLICY;

typedef vcp::pidblas VPOLICY;



int main(void){

    std::cout.precision(17);
    
    

    vcp::matrix< AppData, POLICY > uh;
    vcp::matrix< AppData, POLICY > vh;
    vcp::matrix< AppData, POLICY > wh;
    vcp::matrix< AppData, POLICY > uh_load;
    vcp::matrix< AppData, POLICY > vh_load;
    vcp::matrix< AppData, POLICY > wh_load;
    vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
    // vcp::matrix< AppData_ddtodouble , POLICY > uh_before;
    // vcp::matrix< AppData_ddtodouble , POLICY > vh_before;
    // vcp::matrix< AppData_ddtodouble , POLICY > wh_before;
    vcp::matrix< kv::dd , vcp::mats< kv::dd > > uh_before;
    vcp::matrix< kv::dd , vcp::mats< kv::dd > > vh_before;
    vcp::matrix< kv::dd , vcp::mats< kv::dd > > wh_before;
    int Order_legendre = 40;
    int uh_Order_legendre = 40;
    int p = 3;
    int Dimension = 1;
    int Number_of_variables = 1;

    double D=0.06;
    double alpha=0.25;
    double beta=1.25;
    double gamma=1.0;
    double epsilon=0.045000;
    
    double a=2.5;
    double b=0.1;

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
    Approximate_Generator.setting_list();
    // Approximate_Generator.setting_evenlist();

    // output the list => list_uh
    vcp::matrix< int > list_uh = Approximate_Generator.output_list();

    // setting initialization value of uh
     
        //  uh.zeros(list_uh.rowsize(), Number_of_variables);
        //  vh.zeros(list_uh.rowsize(), Number_of_variables);
        //  wh.zeros(list_uh.rowsize(), Number_of_variables);
        //  uh.ones(list_uh.rowsize(), Number_of_variables);
        //  vh.ones(list_uh.rowsize(), Number_of_variables);
        //  wh.ones(list_uh.rowsize(), Number_of_variables);
        //  vcp::load(uh, "remember/saver/N=1000_uh");
        //  vcp::load(vh, "remember/saver/N=1000_uh");
        //  vcp::load(wh, "remember/saver/N=1000_uh");
        //vcp::save(uh, "Data_Nagumo/uh_rand");
        // std::cout << "r = " << r << std::endl;
        // std::cout << "delta = " << delta << std::endl;
        // std::cout << "ch = " << ch << std::endl;
        // std::cout << "count = " << count << std::endl;
        // vcp::load(uh, "remember/saver/finishuh");
        // vcp::load(uh, "remember/saver/finishvh");
        // vcp::load(uh, "remember/saver/finishwh");
        //vcp::load(uh, "remember/saver/uh");
         vcp::load(uh, "order40_eps0.045000_H10_1.03966_count05_uh");
         vcp::load(vh, "order40_eps0.045000_H10_1.03966_count05_vh");
         vcp::load(wh, "order40_eps0.045000_H10_1.03966_count05_wh");
        //  vcp::load(uh_before, "remember/saver/uh_3fun/uh_0.02");
        //  vcp::load(vh_before, "remember/saver/uh_3fun/vh_0.02");
        //  vcp::load(wh_before, "remember/saver/uh_3fun/wh_0.02");
        // vcp::load(vh, "remember/saver/uh_3fun/vh_0.02");
        // vcp::load(wh, "remember/saver/uh_3fun/wh_0.02");
        // vcp::convert(uh_before,uh);
        // vcp::convert(vh_before,vh);
        // vcp::convert(wh_before,wh);
    //uh = uh+0.1*ch;
    
    //uh(0) = uh(0)+10;
    //uh(1) = 5;
    //uh(0) = 0.1*ch;
    //uh =uh*ch;
    //uh =uh*ch;
    //vcp::load(uh, "Data_Nagumo/uh_test2_base40_0.0835");
    //vcp::load(vh, "Data_Nagumo/vh_test2_base40_0.0835");
    //vcp::save(uh, "Data_Nagumo/uh_rand");
    // if (uh.rowsize() <= uh_load.rowsize()){
    //         std::cout << "count = "  << std::endl;
    //             uh = uh_load.submatrix({0, uh.rowsize() - 1}, {0, uh.columnsize() - 1});
    //     }
    //     else {
    //         std::cout << "check else start"  << std::endl;
    //             vcp::matrix<AppData, POLICY> zero;
    //             zero.zeros(uh.rowsize() - uh_load.rowsize(), uh_load.columnsize());
    //             std::cout << "check else "  << std::endl;
    //             uh = vercat(uh_load, zero);
    //             std::cout << "check else fin "  << std::endl;
    //     }
    //     if (vh.rowsize() <= vh_load.rowsize()){
    //             vh = vh_load.submatrix({0, vh.rowsize() - 1}, {0, vh.columnsize() - 1});
    //     }
    //     else {
    //         std::cout << "check else start"  << std::endl;
    //             vcp::matrix<AppData, POLICY> zero;
    //             zero.zeros(vh.rowsize() - vh_load.rowsize(), vh_load.columnsize());
    //             std::cout << "check else "  << std::endl;
    //             vh = vercat(vh_load, zero);
    //             std::cout << "check else fin "  << std::endl;
    //     }
    //     if (wh.rowsize() <= wh_load.rowsize()){
    //             wh = wh_load.submatrix({0, wh.rowsize() - 1}, {0, wh.columnsize() - 1});
    //     }
    //     else {
    //             vcp::matrix<AppData, POLICY> zero;
    //             zero.zeros(wh.rowsize() - wh_load.rowsize(), wh_load.columnsize());
    //             wh = vercat(wh_load, zero);
    //     }
    // uh.resize(list_uh.rowsize(), Number_of_variables);
    // vh.resize(list_uh.rowsize(), Number_of_variables);
    // wh.resize(list_uh.rowsize(), Number_of_variables);

        // ノルム保存用の変数を追加
        AppData H10_norm ;
        AppData L2_norm ;

        {
        // Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
        vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
        // Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
        vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

        vcp::matrix< AppData, POLICY > uhphi;
        vcp::matrix< AppData, POLICY > uh2phi;//追加
        vcp::matrix< AppData, POLICY > uh3phi;
        vcp::matrix< AppData, POLICY > uhphiphi;//追加
        vcp::matrix< AppData, POLICY > uh2phiphi;
        vcp::matrix< AppData, POLICY > DF;
        vcp::matrix< AppData, POLICY > F;
        vcp::matrix< AppData, POLICY > syuusei;
        vcp::matrix< AppData, POLICY > check;
        vcp::matrix< AppData, POLICY > A;

        {
            AppData cc;
            while(1){
                Approximate_Generator.setting_uh(uh);
                uhphi = Approximate_Generator.uhphi(1);
                uh2phi = Approximate_Generator.uhphi(2);//追加
                uh3phi = Approximate_Generator.uhphi(3);
                uhphiphi = Approximate_Generator.uhphiphi(1);//追加
                uh2phiphi = Approximate_Generator.uhphiphi(2);

                using std::pow;

                /*
                DF = DL - LAMBDA * ( -a*L + 2.0*(1.0+a)*uhphiphi  - 3.0*uh2phiphi);
                F = DL * uh - LAMBDA *( -a*uhphi + (1.0+a)*uh2phi -uh3phi);
                */
                //std::cout << "Newton Method Real Start " << std::endl;
                DF = DL - (-alpha *L +2*beta * uhphiphi -3 * gamma * uh2phiphi-epsilon*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)) )/pow(epsilon,2);
                //std::cout << "DFOK " << std::endl;
                //F = transpose(DL * uh) - (-alpha *transpose(uhphi) +beta * transpose(uh2phi) -gamma * transpose(uh3phi)-delta *transpose(uhphi)*L*lss(-DL+r*L,L))/pow(epsilon,2);                  
                F = DL * uh - (-alpha *uhphi +beta * uh2phi -gamma * uh3phi-epsilon*(a*L*lss(DL+L,uhphi)+b*L*lss(D*D*DL+L,uhphi)) )/pow(epsilon,2);  
                //std::cout << "FOK " << std::endl;
                //syuusei = lss(DF, transpose(F));
                syuusei = lss(DF, F);
                uh = uh - syuusei;
                check = max(abs(syuusei));
                cc = check(0);
                std::cout << cc << std::endl;
                if (cc < pow(2.0,-30)) {
                    //A.eye((L*(DL+r*L)).columnsize());
                    std::cout << cc << std::endl;
                    vh=lss((DL+L),L)*uh;
                    wh=lss((D*D*DL+L),L)*uh;
                    Approximate_Generator.setting_uh(uh);

                    // --- 収束したタイミングで H10ノルム と L2ノルム を計算 ---
                    H10_norm = sqrt(transpose(uh) * DL * uh)(0);
                    L2_norm  = sqrt(transpose(uh) * L * uh)(0);

                    //std::cout << "Convergence \n" << std::endl;
                    break;
                }
            }
        }
    }
    //(-DL+r*L),I.eye(L.rowsize())
    vcp::time.toc();

    //vcp::save(list_uh, "Data_Nagumo/test_list_Base40");
    //vcp::save(uh, "Data_3Nagumo/uh_test2_base40_0.01");
    //vcp::save(vh, "Data_Nagumo/vh_test2_base40_0.0825");
    //vcp::save(uh, "syushi/rand");
    // vcp::save(uh, "remember/saver/uh_3fun/uh_0.025");
    // vcp::save(vh, "remember/saver/uh_3fun/vh_0.025");
    // vcp::save(wh, "remember/saver/uh_3fun/wh_0.025");
    // vcp::save(uh, "remember/saver/N=50_randuh");
    // vcp::save(vh, "remember/saver/N=50_randvh");
    // vcp::save(wh, "remember/saver/N=50_randwh");
    // vcp::save(uh, "asahi/saver/N=40_randuh");
    // uh data for Grafics
    vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
    std::cout << Grafics << std::endl;//実験する時にコメント外す
    std::cout << uh << std::endl;//実験する時にコメント外す

    // minimal and maximum value of approximate solution uh

    std::cout << "\nCalculate the maximum and minimum value" << std::endl;
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

    for (int i = 0; i < Number_of_variables; i++) {
        std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
        //std::cout << "vh in [" << vh_min[i] << ", " << vh_max[i] << "]" << std::endl;
    }

    // --- L_inf ノルムの計算と、3種類のノルムの出力 ---
    double Linf_norm = std::max(std::abs(uh_min[0]), std::abs(uh_max[0]));
    std::cout << "|| uh ||_H10  : " << H10_norm << std::endl;
    std::cout << "|| uh ||_L2   : " << L2_norm << std::endl;
    std::cout << "|| uh ||_Linf : " << Linf_norm << std::endl;


    // r=r+0.1;
    // if(uh_max[0] >=0.0000001){
    //     std::cout<< ch <<std::endl;
    // break;
    // }

    // // if(ch>10.0){
    // //     r =r +0.2;
    // //     ch=0;
    // //     if(r>2.0){
    // //         delta =delta +0.2;
    // //         r=0.2;
    // //     }
    // // }
    // // if(count == 30){
    // //     ch=ch+1;
    // //     count = 0;
    // //     }
    // // count =count +1;
    // if(a>=2.0){
    //     //ch+=1;
    // a=0.1;
    // if(b>=2.0){
    // b=0.1;
    
    // if(c>=2.0){
    //     c=0.1;
    //     alpha+=0.1;
    // }
    // c+=0.1;
    // }
    // b+=0.1;
    // }
    // a+=0.1;
    //}
    // ch+=1;


    Approximate_Generator.clear();
    vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/********************************** Eigenvalue of tilde{F} *************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
        vcp::time.tic();
    VData K = VData(0);
    VData CN, CNw, Cp, Cpw, Cs, Cs3, Cs4, Csw4, CpFtilde;
    VData uh_infsup;
    VData uh_Ftilde_norm;
    uh_infsup.lower() = uh_min[0];
    uh_infsup.upper() = uh_max[0];
{
    std::cout << "\nuh_Ftilde_norm" << std::endl;
    vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
    Verification_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
    Verification_Generator.setting_list();
    // Verification_Generator.setting_evenlist();
    vcp::matrix< VData, VPOLICY > uhi;
    vcp::convert(uh, uhi);
    // uh setting : Last Argument is list divide : full list => 1 , even list => 2 
    Verification_Generator.setting_uh(uhi, list_uh, 1);

    // Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
    vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();
    vcp::matrix< VData, VPOLICY > uhphiphi = Verification_Generator.uhphiphi(1);
    vcp::matrix< VData, VPOLICY > uh2phiphi = Verification_Generator.uhphiphi(2);
    vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();
//    Verification_Generator.clear();
    vcp::matrix< VData, VPOLICY > sigma_tmp =(VData(alpha)*L + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)))/pow(epsilon,2);
    VData sigma =max(max(sigma_tmp))(0);
    std::cout << "sigma=" << sigma<< std::endl;

    vcp::matrix< VData, VPOLICY > G = DL + (VData(alpha)*L + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)))/pow(epsilon,2)+sigma;
    DL.clear();
    uh_Ftilde_norm = (transpose(uhi)*G*uhi)(0);
//uh_Ftilde_normの計算
    std::cout << "uh_Ftilde_norm = (transpose(uhi)*G*uhi)(0); OK" << std::endl;
     
}
//残りの計算
{
    std::cout << "\nEigenvalue of tilde{F}" << std::endl;
    vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
    Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
    Verification_Generator.setting_list();
    // Verification_Generator.setting_evenlist();
    vcp::matrix< VData, VPOLICY > uhi;
    vcp::convert(uh, uhi);
    // uh setting : Last Argument is list divide : full list => 1 , even list => 2 
    Verification_Generator.setting_uh(uhi, list_uh, 1);

    // Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
    vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();
    vcp::matrix< VData, VPOLICY > uhphiphi = Verification_Generator.uhphiphi(1);
    vcp::matrix< VData, VPOLICY > uh2phiphi = Verification_Generator.uhphiphi(2);
    vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();

    vcp::matrix< VData, VPOLICY > sigma_tmp =(VData(alpha)*L + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)))/pow(epsilon,2);
    VData sigma =max(max(sigma_tmp))(0);
    std::cout << "sigma=" << sigma<< std::endl;
    //VData fdtildeuh_LL2 = (VData(alpha)*L + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)))/pow(epsilon,2);
     vcp::matrix< VData, VPOLICY > fdtildeuh_LL2 = (alpha*L + 3*gamma*uh2phiphi +epsilon*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)))/pow(epsilon,2)+sigma;//+sigma
//    Verification_Generator.clear();
    // for (int i = 0; i < uh2phiphi.rowsize(); i++){
    //     std::cout << "uh2phiphi("<< i<<") ="<< uh2phiphi(i)<< std::endl;  
    // }
    // for (int i = 0; i < uh.rowsize(); i++){
    //     std::cout << "uh2phiphi("<< i<<") ="<< uh2phiphi(i)<< std::endl;  
    // }
    //std::cout << uh << std::endl;
    
    std::cout << "tete" << std::endl;
    vcp::matrix< VData, VPOLICY > G = DL + (VData(alpha)*L + 3*VData(gamma)*uh2phiphi +VData(epsilon)*(a*L*lss(DL+L,L)+b*L*lss(D*D*DL+L,L)))/pow(epsilon,2)+sigma;//+sigma
    std::cout << "tete" << std::endl;
    //std::cout << G << std::endl;
    DL.clear();
    std::cout << "DL.clear(); OK" << std::endl;
//    uh_Ftilde_norm = (transpose(uhi)*G*uhi)(0);
//    std::cout << "uh_Ftilde_norm = (transpose(uhi)*G*uhi)(0); OK" << std::endl;
    uh2phiphi.clear();
    std::cout << "uh2phiphi.clear(); OK" << std::endl;
    vcp::matrix< VData, VPOLICY > E;
    
    compsym(G);
    std::cout << " compsym(G); OK" << std::endl;
    eigsymge(L, G, E);
    std::cout << " eigsymge(L, G, E); OK" << std::endl;
    L.clear();
    std::cout << " L.clear(); OK" << std::endl;
    E = 1/diag(E);
    std::cout << " E = 1/diag(E); OK" << std::endl;
    std::cout << "Minmum Eigenvalue of tilde{F}:" << std::endl;
    std::cout << min(E) << std::endl;
    
    CN = Verification_Generator.Ritz_projection_error< VData >();
    std::cout << "CN = " << CN << std::endl;

    Cs3 = Verification_Generator.Sobolev_constant< VData >(3);
    std::cout << "Cs3 = " << Cs3 << ", p =" << "3" << std::endl;

    Cs4 = Verification_Generator.Sobolev_constant< VData >(4);
    std::cout << "Cs4 = " << Cs4 << ", p =" << "4" << std::endl;
    std::cout << "fdtildeuh_LL2(0) = " << fdtildeuh_LL2(0) << std::endl;
    std::cout << "fdtildeuh_LL2.rowsize() = " << fdtildeuh_LL2.rowsize() << std::endl;
    std::cout << "fdtildeuh_LL2.columnsize() = " << fdtildeuh_LL2.columnsize() << std::endl;
    std::cout << "max(max(fdtildeuh_LL2))(0) = " << max(max(fdtildeuh_LL2))(0) << std::endl;
    std::cout << "max(max(fdtildeuh_LL2)).rowsize() = " << max(max(fdtildeuh_LL2)).rowsize()<< std::endl;
    std::cout << "max(max(fdtildeuh_LL2)).columnsize() = " << max(max(fdtildeuh_LL2)).columnsize()<< std::endl;
    VData CFtilde = CN*(1 + 1/min(E)(0)* max(max(fdtildeuh_LL2))(0));
    std::cout << "CFtilde = " << CFtilde << std::endl;

    VData lambda_Ftilde = min(E)(0)/(1+pow(CFtilde,2)*min(E)(0));
    std::cout << "lambda_Ftilde = " << lambda_Ftilde << std::endl;

    CpFtilde = 1/sqrt(lambda_Ftilde);
    CpFtilde.lower() = CpFtilde.upper();
    std::cout << "CpFtilde = " << CpFtilde << std::endl;
    //std::cout <<  uhphiphi << std::endl;
    //std::cout << G << std::endl;
    eigsymge( uhphiphi *VData(2.0)*beta/(VData(epsilon)*VData(epsilon)),G, E);
    std::cout << " eigsymge(uhphiphi, G, E); OK" << std::endl;
    //std::cout << E << std::endl;
    G.clear();
    uhphiphi.clear();
    E=diag(E);
    std::cout << E << std::endl;
    std::cout << " E=diag(E); OK" << std::endl;
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
    
    else if(n == 0){
        O.zeros(1,1);
        E = max(abs( E )); 
        O(0).lower()=E(0).lower();
        O(0).upper()=E(0).upper();
    }
    std::cout << " n" << n << std::endl;
    std::cout << " E to O OK" << std::endl;
    std::cout << O << std::endl;
    O =1/O;
    //std::cout << E << std::endl;
    std::cout << O << std::endl;
    //E = 1/diag(E);
    //std::cout << " E = 1/diag(E); OK" << std::endl;
    VData CFtilde2 = sqrt(abs(2*VData(beta)/pow(epsilon,2)*uh_infsup+sigma))*CFtilde;//+sigma
    CFtilde2.lower() = CFtilde2.upper();
    std::cout << "CFtilde2 = " << CFtilde2 << std::endl;
    // for (int i = 0; i < E.rowsize(); i++ ){
    //     std::cout << "E("<< i<<") ="<< E(i)<< std::endl;
    // }
    for (int i = 0; i < O.rowsize(); i++){
        O(i).lower() = (O(i)/(1+pow(CFtilde2,2)*O(i))).lower();
    }
    std::cout << "for1 OK " << O << std::endl;
    for (int i = 0; i < O.rowsize(); i++){
        if ( O(i).upper() < (VData(2.0)*beta/pow(epsilon,2)).lower() ){
            O(i).lower() = O(i).upper();
        }
        else if ( O(i).lower() > (VData(2.0)*beta/pow(epsilon,2)).upper() ){
            O(i).upper() = O(i).lower();
        }
        else {
            std::cout << "Verification failed..." << std::endl;
        }
    }
    std::cout << "for2 OK " << O << std::endl;
    // for (int i = 0; i < E.rowsize(); i++ ){
    //     E(i) = E(i)/(E(i)-VData(2.0)*beta/pow(epsilon,2));
    // }
    for (int i = 0; i < O.rowsize(); i++ ){
        O(i) = 1/(1-1/O(i));
    }
    std::cout << "for3 OK " << O << std::endl;
    std::cout << "O = " << O << std::endl;    
    O = max(abs( O ));
    std::cout << "O = " << O << std::endl;    

    K = O(0).upper();
    std::cout << "K = " << K << std::endl;
}

    vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/******************* Calculate Residual Norm || Laplace(uh) - f(uh) ||_L2 **********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
    vcp::time.tic();
    VData Res = VData(0);
    {
        Number_of_variables=2;//追加
        std::cout << "\nCalculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
        vcp::Legendre_Bases_Generator< DataType, VResData, VResPOLICY > Verification_Generator;
        Verification_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 2);
        // Setting the list of Verification_Generator 
        Verification_Generator.setting_list();
        // Verification_Generator.setting_evenlist();
        vcp::matrix< VResData, VResPOLICY > uhi;
        vcp::matrix< VResData, VResPOLICY > vhi;
        vcp::matrix< VResData, VResPOLICY > whi;
        vcp::matrix< VResData, VResPOLICY > zhi;
        vcp::matrix< VResData, VResPOLICY > zhi2;
        vcp::matrix< VResData, VResPOLICY > zhi3;

        // vcp::interval(uh, uhi);
        // vcp::interval(vh, vhi);
        // vcp::interval(wh, whi);
        vcp::convert(uh, uhi);
        vcp::convert(vh, vhi);
        vcp::convert(wh, whi);
/************************************************uhとvhでできる計算 ***********************************************/

        
        //vcp::convert(uh, uhi);
        //vcp::convert(vh, vhi);
        zhi=horzcat(uhi,vhi);
        //zhi=vercat(uhi,vhi);
        Verification_Generator.setting_uh(zhi);
        // || Laplace(uh) - f(uh) ||_L2 = sqrt( | (Laplace(uh), Laplace(uh))_L2 + 2(-Laplace(uh), f(uh))_L2 + (f(uh), f(uh))_L2 | )
        // std::cout<< "OK" << std::endl;;
        VResData uh2 = Verification_Generator.integral_uh(2,0);
        VResData uh3 = Verification_Generator.integral_uh(3,0);        
        VResData uh4 = Verification_Generator.integral_uh(4,0);
        VResData uh5 = Verification_Generator.integral_uh(5,0);
        VResData uh6 = Verification_Generator.integral_uh(6,0);
        VResData vh2 = Verification_Generator.integral_uh(0, 2);//v^2
        VResData uhvh = Verification_Generator.integral_uh(1, 1); // (uhphi, vhphi)
        VResData uh2vh = Verification_Generator.integral_uh(2, 1); // (uhphi, uhvhphi)
        VResData uh3vh = Verification_Generator.integral_uh(3, 1); // (uh^2phi, uhvhphi)
        std::cout << "uh2 = " << uh2 << std::endl;
        std::cout << "uh3 = " << uh3 << std::endl;
        std::cout << "uh4 = " << uh4 << std::endl;
        std::cout << "uh5 = " << uh5 << std::endl;
        std::cout << "uh6 = " << uh6 << std::endl; 
        std::cout << " vh2= " << vh2 << std::endl;
        std::cout << " uhvh= " << uhvh << std::endl;
        std::cout << " uh2vh= " << uh2vh << std::endl;
        std::cout << " uh3vh= " << uh3vh << std::endl;
                                       

        VResData LuhLuh = Verification_Generator.integral_LuhLuh(0);
        
        // std::cout<< "OK" << std::endl;
        VResData Luh_uh1 = Verification_Generator.integral_Luhuh(0, 1, 0);
        // std::cout<< "OK2" << std::endl;
        VResData Luh_uh2 = Verification_Generator.integral_Luhuh(0, 2, 0);
        // std::cout<< "OK3" << std::endl;
        VResData Luh_uh3 = Verification_Generator.integral_Luhuh(0, 3, 0); 
        // std::cout<< "OK4" << std::endl;
        VResData Luh_vh = Verification_Generator.integral_Luhuh(0, 0, 1); 
        // std::cout<< "OK5" << std::endl;
        VResData Lvh_Lvh = Verification_Generator.integral_LuhLuh(1); // (DL vh, DL vh)
        // std::cout<< "OK6" << std::endl;
        VResData Lvh_uh = Verification_Generator.integral_Luhuh(1, 1, 0); // (DL vh, vhphi)
        VResData Lvh_vh = Verification_Generator.integral_Luhuh(1, 0, 1); // (DL vh, vhphi) 
        std::cout << "LuhLuh = " << LuhLuh << std::endl;      
        std::cout << "Luh_uh1 = " << Luh_uh1 << std::endl;
        std::cout << "Luh_uh2 = " << Luh_uh2 << std::endl;
        std::cout << "Luh_uh3 = " << Luh_uh3 << std::endl;
        std::cout << " Luh_vh= " << Luh_vh << std::endl;
        std::cout << " Lvh_Lvh= " << Lvh_Lvh << std::endl;
        std::cout << " Lvh_uh= " << Lvh_uh << std::endl;
        std::cout << " Lvh_vh= " << Lvh_vh << std::endl;


        VResData firstuv = LuhLuh;
        VResData seconduv = -VResData(2)*(VResData(alpha)*Luh_uh1-VResData(beta)*Luh_uh2+VResData(gamma)*Luh_uh3+VResData(epsilon)*VResData(a)*Luh_vh)/pow(epsilon,2);
        // VResData thirduv = (pow(epsilon,2)*pow(a,2)*vh2-2*epsilon*a*uhvh+pow(alpha,2)*uh2-2*a*beta*epsilon*uh2vh+2*epsilon*gamma*a*uh3vh-2*alpha*beta*uh3+2*alpha*gamma*uh4+pow(beta,2)*uh4-2*beta*gamma*uh5+pow(gamma,2)*uh6)/pow(epsilon,4);
        VResData thirduv = 2*alpha*a*uhvh/pow(epsilon,3)+pow(alpha,2)*uh2/pow(epsilon,4)-2*a*beta*uh2vh/pow(epsilon,3)+2*a*gamma*uh3vh/pow(epsilon,3)-2*alpha*beta*uh3/pow(epsilon,4)+(2*alpha*gamma+pow(beta,2))*uh4/pow(epsilon,4)-2*beta*gamma*uh5/pow(epsilon,4)+pow(gamma,2)*uh6/pow(epsilon,4);
        VData CRonetmp = abs(a/epsilon)*CN*sqrt(1+pow(CN,2));
        VResData CRone;
        vcp::convert(CRonetmp,  CRone);

        // VResData fourth = Lvh_Lvh+2*Lvh_vh + vh2;
        VResData fourth = uh2 +2*Lvh_uh-2*uhvh+Lvh_Lvh-2*Lvh_vh+vh2;
        std::cout << "firstuv : " << firstuv << std::endl;
        std::cout << "seconduv : " << seconduv << std::endl;        
        std::cout << "thirduv : " << thirduv << std::endl;
        std::cout << "CRone : " << CRone << std::endl;

        /************************************************uhとwhでできる計算 ***********************************************/  
        zhi2=horzcat(uhi,whi);
        //zhi=vercat(uhi,vhi);
        Verification_Generator.setting_uh(zhi2);
        // || Laplace(uh) - f(uh) ||_L2 = sqrt( | (Laplace(uh), Laplace(uh))_L2 + 2(-Laplace(uh), f(uh))_L2 + (f(uh), f(uh))_L2 | )
        std::cout<< "OK" << std::endl;;
        VResData wh2 = Verification_Generator.integral_uh(0, 2);//v^2
        VResData uhwh = Verification_Generator.integral_uh(1, 1); // (uhphi, vhphi)
        VResData uh2wh = Verification_Generator.integral_uh(2, 1); // (uhphi, uhvhphi)
        VResData uh3wh = Verification_Generator.integral_uh(3, 1); // (uh^2phi, uhvhphi)
        std::cout << " vh2= " << uh6 << std::endl;
        std::cout << " uhwh= " << uhwh << std::endl;
        std::cout << " uh2wh= " << uh2wh << std::endl;
        std::cout << " uh3wh= " << uh3wh << std::endl;

        VResData Luh_wh = Verification_Generator.integral_Luhuh(0, 0, 1); 
        std::cout<< "OK5" << std::endl;
        VResData Lwh_Lwh = Verification_Generator.integral_LuhLuh(1); // (DL vh, DL vh)
        std::cout<< "OK6" << std::endl;
        VResData Lwh_uh = Verification_Generator.integral_Luhuh(1, 1, 0); // (DL vh, vhphi) 
        VResData Lwh_wh = Verification_Generator.integral_Luhuh(1, 0, 1); // (DL vh, vhphi) 

        VResData seconduw = -2*(VResData(b)*Luh_wh)/epsilon;
        VResData thirduw = (2*alpha*b*uhwh-2*b*beta*uh2wh+2*b*gamma*uh3wh)/pow(epsilon,3);
        VData CRtwotmp = abs(b/epsilon)*CN*sqrt(1+pow(CN,2)/pow(D,2));
        VResData CRtwo;
        vcp::convert(CRtwotmp,  CRtwo);
        // VResData fifth =  pow(D,4)*Lwh_Lwh+2*pow(D,2)*Lwh_wh + wh2;
        VResData fifth =  uh2/pow(D,4)+2*Lwh_uh/pow(D,2)-2*uhwh/pow(D,4)+Lwh_Lwh-2*Lwh_wh/pow(D,2)+wh2/pow(D,4);
        /************************************************vhとwhでできる計算 ***********************************************/
        zhi3=horzcat(vhi,whi);
        Verification_Generator.setting_uh(zhi3);
        VResData vh2_vw = Verification_Generator.integral_uh(2, 0);//v^2
        VResData wh2_vw = Verification_Generator.integral_uh(0, 2);//v^2
        VResData vhwh = Verification_Generator.integral_uh(1, 1); 

        VResData thirdvw =(pow(a,2)*vh2_vw+2*a*b*vhwh+pow(b,2)*wh2_vw)/pow(epsilon,2);

        //ここですべてまとめる。
        VResData first = firstuv;
        VResData second = seconduv +seconduw ;
        VResData third = thirduv + thirduw + thirdvw ;
        VResData firstothird=first+second+third;
        std::cout << " first = " << first << std::endl;
        std::cout << " second = " << second << std::endl;
        std::cout << " third = " << third << std::endl;
        std::cout << " fourth = " << fourth << std::endl;
        std::cout << " fifth = " << fifth << std::endl;
        std::cout << " CRone = " << CRone << std::endl;
        std::cout << " CRtwo = " << CRtwo << std::endl;
        std::cout << " first to third = " << sqrt(abs(first + second + third)) << std::endl;

        {
            using std::sqrt;
            using std::abs;
            vcp::convert(sqrt(abs(first + second + third))+ CRone*sqrt(abs(fourth))+CRtwo*sqrt(abs(fifth)), Res);
            std::cout << "Residual Norm : || Laplace(uh) - f(uh) ||_L2 <= " << Res << std::endl;
        }
        Res = CpFtilde * Res;
        std::cout << "Residual Norm : || F(uh) ||_(H-1) <= " << Res << std::endl;
    }
    vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/******************* G **********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
    vcp::time.tic();
    VData G = VData(0);
    {
        G = 2*beta/pow(epsilon,2)*pow(Cs3,3)+3*gamma/pow(epsilon,2)*pow(Cs4,4)*(2*uh_Ftilde_norm+4*K*Res); //+ 3*VData(LAMBDA)*pow(Cs4,4)*(2*uh_Ftilde_norm + 4*K*Res);
        std::cout << "G = " << G << std::endl;
    }


/////////////////////////////////////////////////////////////////////////////////////////////////
/******************* Kantorovich **********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
    VData Check = pow(K,2)*Res*G;
    std::cout << "K^2*delta*G = " << Check << std::endl;
    if (Check.upper() <= 0.5 ){
        std::cout << "Verification Succeed!" << std::endl;
        VData rho = (1 - sqrt(1 - 2*Check))/(K*G);
        std::cout << "|| u*-uh || <= " << rho << std::endl;
    }
    else {
        std::cout << "Verification failed..." << std::endl;
    }
    return 0;
}