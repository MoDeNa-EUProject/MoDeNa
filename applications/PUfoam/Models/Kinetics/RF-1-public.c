//
// State variables globals
//


enum statevar{n_Catalyst_1=0, n_CE_A0=1, n_CE_A1=2, 
n_CE_B=3, n_CE_B2=4, n_CE_I0=5, 
n_CE_I1=6, n_CE_I2=7, n_CE_PBA=8, 
n_CE_Breac=9, n_CE_Areac0=10, n_CE_Areac1=11, 
n_CE_Ireac0=12, n_CE_Ireac1=13, n_CE_Ireac2=14, 
n_Bulk=15, n_R_1=16, n_R_1_mass=17, 
n_R_1_temp=18, n_R_1_vol=19}; 
int DIM = 20;
double XINI[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
char* XNAMES[20] = {"Catalyst_1","CE_A0","CE_A1","CE_B","CE_B2","CE_I0","CE_I1","CE_I2","CE_PBA","CE_Breac","CE_Areac0","CE_Areac1","CE_Ireac0","CE_Ireac1","CE_Ireac2","Bulk","R_1","R_1_mass","R_1_temp","R_1_vol"};
double XTYPES[20] = {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,2,2,2,2};
// static double* xGlobal = 0;
static double tGlobal = 0;
double max(double x1, double x2) {return (x1 > x2) ? x1 : x2;}
double min(double x1, double x2) {return (x1 < x2) ? x1 : x2;}
double power(double x1, double x2) {return pow(x1,x2);}
//
// Parameter globals
//
enum params{n_kH_prim_OH_EO_Solvent=0, n_kS_prim_OH_EO_Solvent=1, n_kS_sec_OH_PO_Solvent=2, 
n_kH_sec_OH_PO_Solvent=3, n_kS_prim_OH_BDO_Solvent=4, n_kH_prim_OH_BDO_Solvent=5, 
n_kS_prim_OH_EO_Cat1=6, n_kS_prim_OH_EO_Cat2=7, n_kH_sec_OH_PO_Cat1=8, 
n_kS_prim_OH_EO_Cat3=9, n_kH_prim_OH_EO_Cat1=10, n_kS_sec_OH_PO_Cat1=11, 
n_kH_prim_OH_EO_Cat2=12, n_kH_prim_OH_EO_Cat3=13, n_kH_sec_OH_PO_Cat2=14, 
n_kH_sec_OH_PO_Cat3=15, n_kS_sec_OH_PO_Cat2=16, n_kS_sec_OH_PO_Cat3=17, 
n_kS_prim_OH_BDO_Cat1=18, n_kS_prim_OH_BDO_Cat2=19, n_kS_prim_OH_BDO_Cat3=20, 
n_kH_prim_OH_BDO_Cat1=21, n_kH_prim_OH_BDO_Cat2=22, n_kH_prim_OH_BDO_Cat3=23, 
n_kH_H2O_prim=24, n_kS_H2O_prim=25, n_kH_H2O_sec=26, 
n_kS_H2O_sec=27, n_kH_H2O_BDO=28, n_kS_H2O_BDO=29, 
n_kH_H2O_H2O=30, n_kS_H2O_H2O=31, n_kH_H2O_Cat1=32, 
n_kH_H2O_Cat1s=33, n_kH_H2O_Cat2=34, n_kH_H2O_Cat2s=35, 
n_kH_H2O_Cat3=36, n_kH_H2O_Cat3s=37, n_kS_H2O_Cat1=38, 
n_kS_H2O_Cat1s=39, n_kS_H2O_Cat2=40, n_kS_H2O_Cat2s=41, 
n_kS_H2O_Cat3=42, n_kS_H2O_Cat3s=43, n_f_cat=44, 
n_CS=45, n_CH=46, n_k_A1=47, 
n_k_A0=48, n_k_I1=49, n_k_I2=50, 
n_f_oh_a0=51, n_f_oh_a2=52, n_deltah_a0=53, 
n_deltah_a1=54, n_deltah_a2=55, n_deltah_b1=56, 
n_deltah_b2=57, n_f_oh_a1=58, n_f_oh_b=59, 
n_f_oh_b2=60, n_convcrit=61, n_t_mix=62, 
n_kzero=63, n_kone=64, n_ktwo=65, 
n_kthree=66, n_f_coeff_cross=67, n_f_coeff_end=68, 
n_f_coeff_cat=69, n_t_phase_transfer_A=70, n_kslow=71, 
n_kfast=72, n_kdummy=73, n_bubble_conc=74, 
n_bubble_diam=75, n_diffusion=76, n_sherwood=77, 
n_cp_A0=78, n_cp_I0=79, n_cp_I1=80, 
n_cp_I2=81, n_cp_B2=82, n_cp_PBA=83, 
n_cp_Cat=84, n_cp_A1=85}; 
//
int noConst = 86;
double pGlobal[86] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double pGlobalIso[86] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double pGlobalE[86] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double pGlobalF[86] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double pGlobalType[86] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double pControl[86] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//Mark here by a 1 which parameters should be controlled from outside
//These parameters may not be of type Arrhenius, see pGlobalType
//
// General globals
//
double TEND = 5000; // end time
double TOL = 0.001; // tolerance
double TINI = 0.0001; // initial stepsize
double RGas =1;
double Kelvin = 273.15;

//Prepared for later extension: 
int DIMLIB = 1;
char* LIBNAMES[1] = {"Name_1"};
//

void SetIniState(double* x){
x[n_Catalyst_1] = 6.73000e-02;
x[n_CE_A0] = 1.92250e+00;
x[n_CE_A1] = 2.26920e+00;
x[n_CE_B] = 0.00000e+00;
x[n_CE_B2] = 5.46200e-01;
x[n_CE_I0] = 2.19790e+00;
x[n_CE_I1] = 1.64000e+00;
x[n_CE_I2] = 1.71030e+00;
x[n_CE_PBA] = 0.00000e+00;
x[n_CE_Breac] = 0.00000e+00;
x[n_CE_Areac0] = 0.00000e+00;
x[n_CE_Areac1] = 0.00000e+00;
x[n_CE_Ireac0] = 0.00000e+00;
x[n_CE_Ireac1] = 0.00000e+00;
x[n_CE_Ireac2] = 0.00000e+00;
x[n_Bulk] = 4.45849e+00;
x[n_R_1] = 0.00000e+00;
x[n_R_1_mass] = 1.00000e+00;
x[n_R_1_temp] = 2.27000e+01;
x[n_R_1_vol] = 8.46382e-01;
}

void SetGlobalIniState(){
  SetIniState(XINI);
}

void SetMolweight(double* x){
x[n_Catalyst_1] = 0.00000e+00;
x[n_CE_A0] = 0.00000e+00;
x[n_CE_A1] = 0.00000e+00;
x[n_CE_B] = 0.00000e+00;
x[n_CE_B2] = 0.00000e+00;
x[n_CE_I0] = 0.00000e+00;
x[n_CE_I1] = 0.00000e+00;
x[n_CE_I2] = 0.00000e+00;
x[n_CE_PBA] = 0.00000e+00;
x[n_CE_Breac] = 0.00000e+00;
x[n_CE_Areac0] = 0.00000e+00;
x[n_CE_Areac1] = 0.00000e+00;
x[n_CE_Ireac0] = 0.00000e+00;
x[n_CE_Ireac1] = 0.00000e+00;
x[n_CE_Ireac2] = 0.00000e+00;
x[n_Bulk] = 2.65000e-01;
x[n_R_1] = 0.00000e+00;
x[n_R_1_mass] = 0.00000e+00;
x[n_R_1_temp] = 0.00000e+00;
x[n_R_1_vol] = 0.00000e+00;
}

void SetIniParam(double* p){
p[n_kH_prim_OH_EO_Solvent] = 0.00000e+00;
p[n_kS_prim_OH_EO_Solvent] = 0.00000e+00;
p[n_kS_sec_OH_PO_Solvent] = 0.00000e+00;
p[n_kH_sec_OH_PO_Solvent] = 0.00000e+00;
p[n_kS_prim_OH_BDO_Solvent] = 0.00000e+00;
p[n_kH_prim_OH_BDO_Solvent] = 0.00000e+00;
p[n_kS_prim_OH_EO_Cat1] = 0.00000e+00;
p[n_kS_prim_OH_EO_Cat2] = 0.00000e+00;
p[n_kH_sec_OH_PO_Cat1] = 0.00000e+00;
p[n_kS_prim_OH_EO_Cat3] = 0.00000e+00;
p[n_kH_prim_OH_EO_Cat1] = 0.00000e+00;
p[n_kS_sec_OH_PO_Cat1] = 0.00000e+00;
p[n_kH_prim_OH_EO_Cat2] = 0.00000e+00;
p[n_kH_prim_OH_EO_Cat3] = 0.00000e+00;
p[n_kH_sec_OH_PO_Cat2] = 0.00000e+00;
p[n_kH_sec_OH_PO_Cat3] = 0.00000e+00;
p[n_kS_sec_OH_PO_Cat2] = 0.00000e+00;
p[n_kS_sec_OH_PO_Cat3] = 0.00000e+00;
p[n_kS_prim_OH_BDO_Cat1] = 0.00000e+00;
p[n_kS_prim_OH_BDO_Cat2] = 0.00000e+00;
p[n_kS_prim_OH_BDO_Cat3] = 0.00000e+00;
p[n_kH_prim_OH_BDO_Cat1] = 0.00000e+00;
p[n_kH_prim_OH_BDO_Cat2] = 0.00000e+00;
p[n_kH_prim_OH_BDO_Cat3] = 0.00000e+00;
p[n_kH_H2O_prim] = 0.00000e+00;
p[n_kS_H2O_prim] = 0.00000e+00;
p[n_kH_H2O_sec] = 0.00000e+00;
p[n_kS_H2O_sec] = 0.00000e+00;
p[n_kH_H2O_BDO] = 0.00000e+00;
p[n_kS_H2O_BDO] = 0.00000e+00;
p[n_kH_H2O_H2O] = 0.00000e+00;
p[n_kS_H2O_H2O] = 0.00000e+00;
p[n_kH_H2O_Cat1] = 0.00000e+00;
p[n_kH_H2O_Cat1s] = 0.00000e+00;
p[n_kH_H2O_Cat2] = 0.00000e+00;
p[n_kH_H2O_Cat2s] = 0.00000e+00;
p[n_kH_H2O_Cat3] = 0.00000e+00;
p[n_kH_H2O_Cat3s] = 0.00000e+00;
p[n_kS_H2O_Cat1] = 0.00000e+00;
p[n_kS_H2O_Cat1s] = 0.00000e+00;
p[n_kS_H2O_Cat2] = 0.00000e+00;
p[n_kS_H2O_Cat2s] = 0.00000e+00;
p[n_kS_H2O_Cat3] = 0.00000e+00;
p[n_kS_H2O_Cat3s] = 0.00000e+00;
p[n_f_cat] = 1.00000e+00;
p[n_CS] = 0.00000e+00;
p[n_CH] = 1.00000e+00;
p[n_k_A1] = 3.00000e+00;
p[n_k_A0] = 3.00000e+00;
p[n_k_I1] = 3.00000e+00;
p[n_k_I2] = 4.00000e+00;
p[n_f_oh_a0] = 2.00000e+00;
p[n_f_oh_a2] = 2.00000e+00;
p[n_deltah_a0] = -8.30000e+01;
p[n_deltah_a1] = -8.30000e+01;
p[n_deltah_a2] = -8.00000e+01;
p[n_deltah_b1] = -4.70000e+01;
p[n_deltah_b2] = -1.27000e+02;
p[n_f_oh_a1] = 2.00000e+00;
p[n_f_oh_b] = 3.00000e+00;
p[n_f_oh_b2] = 4.00000e+00;
p[n_convcrit] = 9.99000e-01;
p[n_t_mix] = 1.00000e+00;
p[n_kzero] = 0.00000e+00;
p[n_kone] = 1.00000e+00;
p[n_ktwo] = 2.00000e+00;
p[n_kthree] = 3.00000e+00;
p[n_f_coeff_cross] = 0.00000e+00;
p[n_f_coeff_end] = 1.00000e+00;
p[n_f_coeff_cat] = 1.00000e+00;
p[n_t_phase_transfer_A] = 1.00000e+00;
p[n_kslow] = 1.00000e+00;
p[n_kfast] = 1.00000e+02;
p[n_kdummy] = 0.00000e+00;
p[n_bubble_conc] = 1.00000e+07;
p[n_bubble_diam] = 2.00000e-05;
p[n_diffusion] = 1.00000e-09;
p[n_sherwood] = 2.00000e+00;
p[n_cp_A0] = 3.08000e-01;
p[n_cp_I0] = 2.12500e-01;
p[n_cp_I1] = 2.28900e-01;
p[n_cp_I2] = 2.31000e-01;
p[n_cp_B2] = 7.52000e-02;
p[n_cp_PBA] = 2.00000e-01;
p[n_cp_Cat] = 2.00000e-01;
p[n_cp_A1] = 1.25400e-01;
}

void SetIniParamE(double* p){
p[n_kH_prim_OH_EO_Solvent] = 3.40000e+03;
p[n_kS_prim_OH_EO_Solvent] = 3.40000e+03;
p[n_kS_sec_OH_PO_Solvent] = 6.00000e+03;
p[n_kH_sec_OH_PO_Solvent] = 5.49500e+03;
p[n_kS_prim_OH_BDO_Solvent] = 4.80000e+03;
p[n_kH_prim_OH_BDO_Solvent] = 3.21000e+03;
p[n_kS_prim_OH_EO_Cat1] = 5.00000e+03;
p[n_kS_prim_OH_EO_Cat2] = 2.50000e+03;
p[n_kH_sec_OH_PO_Cat1] = 4.60000e+03;
p[n_kS_prim_OH_EO_Cat3] = 4.00000e+03;
p[n_kH_prim_OH_EO_Cat1] = 3.00000e+03;
p[n_kS_sec_OH_PO_Cat1] = 0.00000e+00;
p[n_kH_prim_OH_EO_Cat2] = 3.50000e+03;
p[n_kH_prim_OH_EO_Cat3] = 2.50000e+03;
p[n_kH_sec_OH_PO_Cat2] = 4.50000e+03;
p[n_kH_sec_OH_PO_Cat3] = 3.50000e+03;
p[n_kS_sec_OH_PO_Cat2] = 4.50000e+03;
p[n_kS_sec_OH_PO_Cat3] = 4.50000e+03;
p[n_kS_prim_OH_BDO_Cat1] = 3.20000e+03;
p[n_kS_prim_OH_BDO_Cat2] = 3.20000e+03;
p[n_kS_prim_OH_BDO_Cat3] = 2.50000e+03;
p[n_kH_prim_OH_BDO_Cat1] = 3.00000e+03;
p[n_kH_prim_OH_BDO_Cat2] = 2.50000e+03;
p[n_kH_prim_OH_BDO_Cat3] = 2.50000e+03;
p[n_kH_H2O_prim] = 2.50000e+03;
p[n_kS_H2O_prim] = 2.40000e+03;
p[n_kH_H2O_sec] = 2.40000e+03;
p[n_kS_H2O_sec] = 2.40000e+03;
p[n_kH_H2O_BDO] = 2.40000e+03;
p[n_kS_H2O_BDO] = 2.40000e+03;
p[n_kH_H2O_H2O] = 4.50000e+03;
p[n_kS_H2O_H2O] = 4.50000e+03;
p[n_kH_H2O_Cat1] = 5.80000e+03;
p[n_kH_H2O_Cat1s] = 5.50000e+03;
p[n_kH_H2O_Cat2] = 2.50000e+03;
p[n_kH_H2O_Cat2s] = 2.50000e+03;
p[n_kH_H2O_Cat3] = 2.50000e+03;
p[n_kH_H2O_Cat3s] = 2.50000e+03;
p[n_kS_H2O_Cat1] = 3.00000e+03;
p[n_kS_H2O_Cat1s] = 5.50000e+03;
p[n_kS_H2O_Cat2] = 2.50000e+03;
p[n_kS_H2O_Cat2s] = 2.50000e+03;
p[n_kS_H2O_Cat3] = 2.50000e+03;
p[n_kS_H2O_Cat3s] = 2.50000e+03;
p[n_f_cat] = 0.00000e+00;
p[n_CS] = 0.00000e+00;
p[n_CH] = 0.00000e+00;
p[n_k_A1] = 0.00000e+00;
p[n_k_A0] = 0.00000e+00;
p[n_k_I1] = 0.00000e+00;
p[n_k_I2] = 0.00000e+00;
p[n_f_oh_a0] = 0.00000e+00;
p[n_f_oh_a2] = 0.00000e+00;
p[n_deltah_a0] = 0.00000e+00;
p[n_deltah_a1] = 0.00000e+00;
p[n_deltah_a2] = 0.00000e+00;
p[n_deltah_b1] = 0.00000e+00;
p[n_deltah_b2] = 0.00000e+00;
p[n_f_oh_a1] = 0.00000e+00;
p[n_f_oh_b] = 0.00000e+00;
p[n_f_oh_b2] = 0.00000e+00;
p[n_convcrit] = 0.00000e+00;
p[n_t_mix] = 0.00000e+00;
p[n_kzero] = 0.00000e+00;
p[n_kone] = 0.00000e+00;
p[n_ktwo] = 0.00000e+00;
p[n_kthree] = 0.00000e+00;
p[n_f_coeff_cross] = 0.00000e+00;
p[n_f_coeff_end] = 0.00000e+00;
p[n_f_coeff_cat] = 0.00000e+00;
p[n_t_phase_transfer_A] = 0.00000e+00;
p[n_kslow] = 0.00000e+00;
p[n_kfast] = 0.00000e+00;
p[n_kdummy] = 0.00000e+00;
p[n_bubble_conc] = 0.00000e+00;
p[n_bubble_diam] = 0.00000e+00;
p[n_diffusion] = 0.00000e+00;
p[n_sherwood] = 0.00000e+00;
p[n_cp_A0] = 0.00000e+00;
p[n_cp_I0] = 0.00000e+00;
p[n_cp_I1] = 0.00000e+00;
p[n_cp_I2] = 0.00000e+00;
p[n_cp_B2] = 0.00000e+00;
p[n_cp_PBA] = 0.00000e+00;
p[n_cp_Cat] = 0.00000e+00;
p[n_cp_A1] = 0.00000e+00;
}

void SetIniParamF(double* p){
p[n_kH_prim_OH_EO_Solvent] = 1.30000e+02;
p[n_kS_prim_OH_EO_Solvent] = 1.30000e+02;
p[n_kS_sec_OH_PO_Solvent] = 2.00000e+03;
p[n_kH_sec_OH_PO_Solvent] = 1.70000e+04;
p[n_kS_prim_OH_BDO_Solvent] = 2.70000e+02;
p[n_kH_prim_OH_BDO_Solvent] = 2.10000e+02;
p[n_kS_prim_OH_EO_Cat1] = 6.80000e+02;
p[n_kS_prim_OH_EO_Cat2] = 7.00000e+02;
p[n_kH_sec_OH_PO_Cat1] = 3.00000e+02;
p[n_kS_prim_OH_EO_Cat3] = 7.00000e+02;
p[n_kH_prim_OH_EO_Cat1] = 5.00000e+02;
p[n_kS_sec_OH_PO_Cat1] = 4.00000e+02;
p[n_kH_prim_OH_EO_Cat2] = 8.00000e+02;
p[n_kH_prim_OH_EO_Cat3] = 5.00000e+02;
p[n_kH_sec_OH_PO_Cat2] = 7.50000e+03;
p[n_kH_sec_OH_PO_Cat3] = 8.00000e+03;
p[n_kS_sec_OH_PO_Cat2] = 7.00000e+03;
p[n_kS_sec_OH_PO_Cat3] = 8.00000e+03;
p[n_kS_prim_OH_BDO_Cat1] = 5.00000e+02;
p[n_kS_prim_OH_BDO_Cat2] = 8.00000e+02;
p[n_kS_prim_OH_BDO_Cat3] = 9.00000e+02;
p[n_kH_prim_OH_BDO_Cat1] = 9.00000e+02;
p[n_kH_prim_OH_BDO_Cat2] = 9.00000e+02;
p[n_kH_prim_OH_BDO_Cat3] = 7.00000e+02;
p[n_kH_H2O_prim] = 1.00000e+01;
p[n_kS_H2O_prim] = 1.00000e+01;
p[n_kH_H2O_sec] = 1.00000e+01;
p[n_kS_H2O_sec] = 1.00000e+01;
p[n_kH_H2O_BDO] = 1.00000e+01;
p[n_kS_H2O_BDO] = 1.00000e+01;
p[n_kH_H2O_H2O] = 1.00000e+01;
p[n_kS_H2O_H2O] = 1.00000e+01;
p[n_kH_H2O_Cat1] = 1.00000e+06;
p[n_kH_H2O_Cat1s] = 1.00000e+06;
p[n_kH_H2O_Cat2] = 5.00000e+01;
p[n_kH_H2O_Cat2s] = 5.00000e+01;
p[n_kH_H2O_Cat3] = 5.00000e+01;
p[n_kH_H2O_Cat3s] = 5.00000e+01;
p[n_kS_H2O_Cat1] = 1.00000e+02;
p[n_kS_H2O_Cat1s] = 1.00000e+06;
p[n_kS_H2O_Cat2] = 5.00000e+01;
p[n_kS_H2O_Cat2s] = 5.00000e+01;
p[n_kS_H2O_Cat3] = 5.00000e+01;
p[n_kS_H2O_Cat3s] = 5.00000e+01;
p[n_f_cat] = 0.00000e+00;
p[n_CS] = 0.00000e+00;
p[n_CH] = 0.00000e+00;
p[n_k_A1] = 0.00000e+00;
p[n_k_A0] = 0.00000e+00;
p[n_k_I1] = 0.00000e+00;
p[n_k_I2] = 0.00000e+00;
p[n_f_oh_a0] = 0.00000e+00;
p[n_f_oh_a2] = 0.00000e+00;
p[n_deltah_a0] = 0.00000e+00;
p[n_deltah_a1] = 0.00000e+00;
p[n_deltah_a2] = 0.00000e+00;
p[n_deltah_b1] = 0.00000e+00;
p[n_deltah_b2] = 0.00000e+00;
p[n_f_oh_a1] = 0.00000e+00;
p[n_f_oh_b] = 0.00000e+00;
p[n_f_oh_b2] = 0.00000e+00;
p[n_convcrit] = 0.00000e+00;
p[n_t_mix] = 0.00000e+00;
p[n_kzero] = 0.00000e+00;
p[n_kone] = 0.00000e+00;
p[n_ktwo] = 0.00000e+00;
p[n_kthree] = 0.00000e+00;
p[n_f_coeff_cross] = 0.00000e+00;
p[n_f_coeff_end] = 0.00000e+00;
p[n_f_coeff_cat] = 0.00000e+00;
p[n_t_phase_transfer_A] = 0.00000e+00;
p[n_kslow] = 0.00000e+00;
p[n_kfast] = 0.00000e+00;
p[n_kdummy] = 0.00000e+00;
p[n_bubble_conc] = 0.00000e+00;
p[n_bubble_diam] = 0.00000e+00;
p[n_diffusion] = 0.00000e+00;
p[n_sherwood] = 0.00000e+00;
p[n_cp_A0] = 0.00000e+00;
p[n_cp_I0] = 0.00000e+00;
p[n_cp_I1] = 0.00000e+00;
p[n_cp_I2] = 0.00000e+00;
p[n_cp_B2] = 0.00000e+00;
p[n_cp_PBA] = 0.00000e+00;
p[n_cp_Cat] = 0.00000e+00;
p[n_cp_A1] = 0.00000e+00;
}

void SetIniParamType(double* p){
p[n_kH_prim_OH_EO_Solvent] = 1;
p[n_kS_prim_OH_EO_Solvent] = 1;
p[n_kS_sec_OH_PO_Solvent] = 1;
p[n_kH_sec_OH_PO_Solvent] = 1;
p[n_kS_prim_OH_BDO_Solvent] = 1;
p[n_kH_prim_OH_BDO_Solvent] = 1;
p[n_kS_prim_OH_EO_Cat1] = 1;
p[n_kS_prim_OH_EO_Cat2] = 1;
p[n_kH_sec_OH_PO_Cat1] = 1;
p[n_kS_prim_OH_EO_Cat3] = 1;
p[n_kH_prim_OH_EO_Cat1] = 1;
p[n_kS_sec_OH_PO_Cat1] = 1;
p[n_kH_prim_OH_EO_Cat2] = 1;
p[n_kH_prim_OH_EO_Cat3] = 1;
p[n_kH_sec_OH_PO_Cat2] = 1;
p[n_kH_sec_OH_PO_Cat3] = 1;
p[n_kS_sec_OH_PO_Cat2] = 1;
p[n_kS_sec_OH_PO_Cat3] = 1;
p[n_kS_prim_OH_BDO_Cat1] = 1;
p[n_kS_prim_OH_BDO_Cat2] = 1;
p[n_kS_prim_OH_BDO_Cat3] = 1;
p[n_kH_prim_OH_BDO_Cat1] = 1;
p[n_kH_prim_OH_BDO_Cat2] = 1;
p[n_kH_prim_OH_BDO_Cat3] = 1;
p[n_kH_H2O_prim] = 1;
p[n_kS_H2O_prim] = 1;
p[n_kH_H2O_sec] = 1;
p[n_kS_H2O_sec] = 1;
p[n_kH_H2O_BDO] = 1;
p[n_kS_H2O_BDO] = 1;
p[n_kH_H2O_H2O] = 1;
p[n_kS_H2O_H2O] = 1;
p[n_kH_H2O_Cat1] = 1;
p[n_kH_H2O_Cat1s] = 1;
p[n_kH_H2O_Cat2] = 1;
p[n_kH_H2O_Cat2s] = 1;
p[n_kH_H2O_Cat3] = 1;
p[n_kH_H2O_Cat3s] = 1;
p[n_kS_H2O_Cat1] = 1;
p[n_kS_H2O_Cat1s] = 1;
p[n_kS_H2O_Cat2] = 1;
p[n_kS_H2O_Cat2s] = 1;
p[n_kS_H2O_Cat3] = 1;
p[n_kS_H2O_Cat3s] = 1;
p[n_f_cat] = 0;
p[n_CS] = 0;
p[n_CH] = 0;
p[n_k_A1] = 0;
p[n_k_A0] = 0;
p[n_k_I1] = 0;
p[n_k_I2] = 0;
p[n_f_oh_a0] = 0;
p[n_f_oh_a2] = 0;
p[n_deltah_a0] = 0;
p[n_deltah_a1] = 0;
p[n_deltah_a2] = 0;
p[n_deltah_b1] = 0;
p[n_deltah_b2] = 0;
p[n_f_oh_a1] = 0;
p[n_f_oh_b] = 0;
p[n_f_oh_b2] = 0;
p[n_convcrit] = 0;
p[n_t_mix] = 0;
p[n_kzero] = 0;
p[n_kone] = 0;
p[n_ktwo] = 0;
p[n_kthree] = 0;
p[n_f_coeff_cross] = 0;
p[n_f_coeff_end] = 0;
p[n_f_coeff_cat] = 0;
p[n_t_phase_transfer_A] = 0;
p[n_kslow] = 0;
p[n_kfast] = 0;
p[n_kdummy] = 0;
p[n_bubble_conc] = 0;
p[n_bubble_diam] = 0;
p[n_diffusion] = 0;
p[n_sherwood] = 0;
p[n_cp_A0] = 0;
p[n_cp_I0] = 0;
p[n_cp_I1] = 0;
p[n_cp_I2] = 0;
p[n_cp_B2] = 0;
p[n_cp_PBA] = 0;
p[n_cp_Cat] = 0;
p[n_cp_A1] = 0;
}

void SetFixedIniParam(double* p, double* pC, double Temp, double RG){
 for (int i=0; i< 86; i++){
   if (pC[i] == 0){
     p[i] = pGlobalIso[i];
     if (pGlobalType[i] == 1)
        p[i] = pGlobalF[i]*exp(-pGlobalE[i]/RG/(Temp + Kelvin));
  }
 }
}

void SetGlobalIniParam(){
  SetIniParam(pGlobalIso);
  SetIniParamE(pGlobalE);
  SetIniParamF(pGlobalF);
  SetIniParamType(pGlobalType);
//To be called inside the rhs function:
//  SetFixedIniParam(pGlobal, pControl, 100, 1);
}

//
//CiT functions
//
int GetReacIndex(int reactorIndex);
double exp(double expA, double expB);
double getkp(double* pArray, int paramIndex);
double gettemp(double* xArray, int reacIndex);
double getco(double* xArray, int stateIndex);
//################
int GetReacIndex(int reactorIndex)
{
double result1=0.0;
int no = n_R_1;
result1 = no;
return result1;
}
//################
double exp(double expA, double expB)
{
double result1=0.0;
result1 = pow(expA, expB);
return result1;
}
//################
double getkp(double* pArray, int paramIndex)
{
double result1=0.0;
int no = paramIndex;
if(no > noConst) return 0.0;
result1 = pArray[no];
return result1;
}
//################
double gettemp(double* xArray, int reacIndex)
{
double result1=0.0;
int no = GetReacIndex(reacIndex) + 2;
if(no >= DIM) return 0.0;
result1 = xArray[no];
return result1;
}
//################
double getco(double* xArray, int stateIndex)
{
double result1=0.0;
int no = stateIndex;
if(no >= DIM) return 0.0;
if (xArray != 0)
return xArray[no];
return result1;
}

//
//User scripts
//
void f_f_rate_coeff_a0(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2);
double f_f_rate_coefficient(double* xArray, double* pArray, double arg1, double arg2);
void f_f_rate_coeff_a1(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2);
void f_f_rate_coeff_b2(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2);
void f_f_rate_coeff_b1(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2);
//################Script##################################
void f_f_rate_coeff_a0(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2)
{
double f_oh_a = 0;
double deltah = 0;
double f_rate_coeff = 0;
double time_mix = 0;
double time = 0;
double result3 = 0;
double result4 = 0;

// This routine calculates the rate coefficient for the reaction of the end
// group A(0) with the isocyanate end group I
f_oh_a = getkp(pArray, n_f_oh_a0) ;
deltah = getkp(pArray, n_deltah_a0) ;
f_rate_coeff = f_f_rate_coefficient(xArray, pArray,f_oh_a,0) ;
time_mix = getkp(pArray, n_t_phase_transfer_A) ;
time = tGlobal ;
if (time < time_mix) 
{
f_rate_coeff = 0. ;
}
result1 = f_rate_coeff ;
result2 = 0 ;
result3 = deltah ;
result4 = 0 ;

}

//################Script##################################
double f_f_rate_coefficient(double* xArray, double* pArray, double arg1, double arg2)
{
double f_end_group = 0;
double f_oh_a0 = 0;
double f_oh_a1 = 0;
double f_cat = 0;
double fka_Cat_S1 = 0;
double fka_Cat_H1 = 0;
double fka_Cat_S2 = 0;
double fka_Cat_H2 = 0;
double fka_Cat_S3 = 0;
double fka_Cat_H3 = 0;
double fka_prim_S = 0;
double fka_sec_S = 0;
double fka_BDO_S = 0;
double fka_H2O_S = 0;
double fka_prim_H = 0;
double fka_sec_H = 0;
double fka_BDO_H = 0;
double fka_H2O_H = 0;
double cA0 = 0;
double cA1 = 0;
double cH2O = 0;
double cCat1 = 0;
double cprim = 0;
double csec = 0;
double cBDO = 0;
double fka_S = 0;
double fka_H = 0;
double CCA1 = 0;
double CCI1 = 0;
double CCI2 = 0;
double C0 = 0;
double CS = 0;
double CH = 0;
double factor = 0;
double coefficient = 0;
double result1 = 0;

//********************************************************************************************************
// This routine calculates the rate coefficient for the reaction of the various OH and water end
// groups with the isocyanate end group I
// Autor: Andreas Daiss, GME/M
// Status: 16/02/2015
//********************************************************************************************************
// Get type of end group for which reaction shall be calculated: 0=BDO,1=primary,2=secondary,3=water
f_end_group = arg1 ;
// Get the type of the OH-end groups of A and C (primary OH, secondary OH or a primary OH from BDO)
f_oh_a0 = getkp(pArray, n_f_oh_a0) ;
f_oh_a1 = getkp(pArray, n_f_oh_a1) ;
// Determination of rate coefficient as a function of catalyst type: 1=Amin, 2=Tin
// The two rate coefficients for each catalyst type refer to either a non cross-linked "soft" system (S) or a strongly cross-linked "hard" system (H)
f_cat = getkp(pArray, n_f_cat) ;
//********************************************************************************************
// BDO end group
//********************************************************************************************
if ( f_end_group == 0 ) 
{
fka_Cat_S1 = getkp(pArray, n_kS_prim_OH_BDO_Cat1) ;
fka_Cat_H1 = getkp(pArray, n_kH_prim_OH_BDO_Cat1) ;
fka_Cat_S2 = getkp(pArray, n_kS_prim_OH_BDO_Cat2) ;
fka_Cat_H2 = getkp(pArray, n_kH_prim_OH_BDO_Cat2) ;
fka_Cat_S3 = getkp(pArray, n_kS_prim_OH_BDO_Cat3) ;
fka_Cat_H3 = getkp(pArray, n_kH_prim_OH_BDO_Cat3) ;
// Get the rate coefficients for the OH end group and water catalyzed reaction. Again it is distinguished between a non (S) and a strongly (H) cross-linked system.
// Rate coefficients for low cross linker content (S)
fka_prim_S = getkp(pArray, n_kS_prim_OH_BDO_Solvent) ;
fka_sec_S = getkp(pArray, n_kS_prim_OH_BDO_Solvent) ;
fka_BDO_S = getkp(pArray, n_kS_prim_OH_BDO_Solvent) ;
fka_H2O_S = getkp(pArray, n_kS_prim_OH_BDO_Solvent) ;
// Rate coefficients for high cross linker content (H)
fka_prim_H = getkp(pArray, n_kH_prim_OH_BDO_Solvent) ;
fka_sec_H = getkp(pArray, n_kH_prim_OH_BDO_Solvent) ;
fka_BDO_H = getkp(pArray, n_kH_prim_OH_BDO_Solvent) ;
fka_H2O_H = getkp(pArray, n_kH_prim_OH_BDO_Solvent) ;
}
//********************************************************************************************
// Primary end group (non BDO!)
//********************************************************************************************
if ( f_end_group == 1 ) 
{
fka_Cat_S1 = getkp(pArray, n_kS_prim_OH_EO_Cat1) ;
fka_Cat_H1 = getkp(pArray, n_kH_prim_OH_EO_Cat1) ;
fka_Cat_S2 = getkp(pArray, n_kS_prim_OH_EO_Cat2) ;
fka_Cat_H2 = getkp(pArray, n_kH_prim_OH_EO_Cat2) ;
fka_Cat_S3 = getkp(pArray, n_kS_prim_OH_EO_Cat3) ;
fka_Cat_H3 = getkp(pArray, n_kH_prim_OH_EO_Cat3) ;
// Get the rate coefficients for the OH end group and water catalyzed reaction. Again it is distinguished between a non (S) and a strongly (H) cross-linked system.
// Rate coefficients for low cross linker content (S)
fka_prim_S = getkp(pArray, n_kS_prim_OH_EO_Solvent) ;
fka_sec_S = getkp(pArray, n_kS_prim_OH_EO_Solvent) ;
fka_BDO_S = getkp(pArray, n_kS_prim_OH_EO_Solvent) ;
fka_H2O_S = getkp(pArray, n_kS_prim_OH_EO_Solvent) ;
// Rate coefficients for high cross linker content (H)
fka_prim_H = getkp(pArray, n_kH_prim_OH_EO_Solvent) ;
fka_sec_H = getkp(pArray, n_kH_prim_OH_EO_Solvent) ;
fka_BDO_H = getkp(pArray, n_kH_prim_OH_EO_Solvent) ;
fka_H2O_H = getkp(pArray, n_kH_prim_OH_EO_Solvent) ;
}
//********************************************************************************************
// Secondary end group
//********************************************************************************************
if ( f_end_group == 2 ) 
{
fka_Cat_S1 = getkp(pArray, n_kS_sec_OH_PO_Cat1) ;
fka_Cat_H1 = getkp(pArray, n_kH_sec_OH_PO_Cat1) ;
fka_Cat_S2 = getkp(pArray, n_kS_sec_OH_PO_Cat2) ;
fka_Cat_H2 = getkp(pArray, n_kH_sec_OH_PO_Cat2) ;
fka_Cat_S3 = getkp(pArray, n_kS_sec_OH_PO_Cat3) ;
fka_Cat_H3 = getkp(pArray, n_kH_sec_OH_PO_Cat3) ;
;
// Get the rate coefficients for the OH end group and water catalyzed reaction. Again it is distinguished between a non (S) and a strongly (H) cross-linked system.
// Rate coefficients for low cross linker content (S)
fka_prim_S = getkp(pArray, n_kS_sec_OH_PO_Solvent) ;
fka_sec_S = getkp(pArray, n_kS_sec_OH_PO_Solvent) ;
fka_BDO_S = getkp(pArray, n_kS_sec_OH_PO_Solvent) ;
fka_H2O_S = getkp(pArray, n_kS_sec_OH_PO_Solvent) ;
// Rate coefficients for high cross linker content (H)
fka_prim_H = getkp(pArray, n_kH_sec_OH_PO_Solvent) ;
fka_sec_H = getkp(pArray, n_kH_sec_OH_PO_Solvent) ;
fka_BDO_H = getkp(pArray, n_kH_sec_OH_PO_Solvent) ;
fka_H2O_H = getkp(pArray, n_kH_sec_OH_PO_Solvent) ;
}
//********************************************************************************************
// Water end group
//********************************************************************************************
if ( f_end_group == 3 ) 
{
fka_Cat_S1 = getkp(pArray, n_kS_H2O_Cat1) ;
fka_Cat_H1 = getkp(pArray, n_kH_H2O_Cat1) ;
fka_Cat_S2 = getkp(pArray, n_kS_H2O_Cat2) ;
fka_Cat_H2 = getkp(pArray, n_kH_H2O_Cat2) ;
fka_Cat_S3 = getkp(pArray, n_kS_H2O_Cat3) ;
fka_Cat_H3 = getkp(pArray, n_kH_H2O_Cat3) ;
// Get the rate coefficients for the OH end group and water catalyzed reaction. Again it is distinguished between a non (S) and a strongly (H) cross-linked system.
// Rate coefficients for low cross linker content (S)
fka_prim_S = getkp(pArray, n_kS_H2O_prim) ;
fka_sec_S = getkp(pArray, n_kS_H2O_sec) ;
fka_BDO_S = getkp(pArray, n_kS_H2O_BDO) ;
fka_H2O_S = getkp(pArray, n_kS_H2O_H2O) ;
// Rate coefficients for high cross linker content (H)
fka_prim_H = getkp(pArray, n_kH_H2O_prim) ;
fka_sec_H = getkp(pArray, n_kH_H2O_sec) ;
fka_BDO_H = getkp(pArray, n_kH_H2O_BDO) ;
fka_H2O_H = getkp(pArray, n_kH_H2O_H2O) ;
}
if ( f_end_group == 4 ) 
{
fka_Cat_S1 = getkp(pArray, n_kS_H2O_Cat1s) ;
fka_Cat_H1 = getkp(pArray, n_kH_H2O_Cat1s) ;
fka_Cat_S2 = getkp(pArray, n_kS_H2O_Cat2s) ;
fka_Cat_H2 = getkp(pArray, n_kH_H2O_Cat2s) ;
fka_Cat_S3 = getkp(pArray, n_kS_H2O_Cat3s) ;
fka_Cat_H3 = getkp(pArray, n_kH_H2O_Cat3s) ;
}
// Get concentrations of non-reacted OH end groups for determination of solvent catalysis effect
cA0 = getco(xArray, n_CE_A0) ;
cA1 = getco(xArray, n_CE_A1) ;
cH2O = getco(xArray, n_CE_B2) ;
cCat1 = getco(xArray, n_Catalyst_1) ;
// Sum up separately the concentrations of of primary, secondary and BDO primary end groups.
cprim = 0 ;
csec = 0 ;
cBDO = 0 ;
if ( f_oh_a0 == 0 ) 
{
cBDO = cBDO + cA0 ;
}
else
{
if (f_oh_a0 == 1 ) 
{
cprim = cprim + cA0 ;
}
else
{
csec = csec + cA0 ;
}
}
if ( f_oh_a1 == 1 ) 
{
cBDO = cBDO + cA1 ;
}
else
{
if (f_oh_a1 == 1 ) 
{
cprim = cprim + cA1 ;
}
else
{
csec = csec + cA1 ;
}
}
// Summing of all rate coefficient (OH end group catatyzed as well as amine or tin catalyzed
//fka_H = cprim*fka_prim_H + csec*fka_sec_H + cBDO*fka_BDO_H + cH2O*fka_H2O_H + cCat1*fka_Cat_H1
fka_S = cprim * fka_prim_S + csec * fka_sec_S + cBDO * fka_BDO_S + cH2O * fka_H2O_S + cCat1 * fka_Cat_S1 ;
fka_H = cCat1 * fka_Cat_H1 + fka_sec_H ;
// Determination of rate coefficient as function of cross linker concentration
// Get concentration of potential cross links
CCA1 = (getco(xArray, n_CE_Areac1) + getco(xArray, n_CE_A1)) * (getkp(pArray, n_k_A1) - 2)/getkp(pArray, n_k_A1) ;
CCI1 = (getco(xArray, n_CE_Ireac1) + getco(xArray, n_CE_I1)) * (getkp(pArray, n_k_I1) - 2)/getkp(pArray, n_k_I1) ;
CCI2 = (getco(xArray, n_CE_Ireac2) + getco(xArray, n_CE_I2)) * (getkp(pArray, n_k_I2) - 2)/getkp(pArray, n_k_I2) ;
C0 = CCA1 + CCI1 + CCI2 ;
// Get upper and lower bound for cross link concentration. Lower bound: Below this concentration the system is treated as non cross-linked.
// Upper bound: Above this concentration the system is treated as "fully" cross-linked. In between the rate coefficient is a function of the potential cross link
// concentration
CS = getkp(pArray, n_CS) ;
CH = getkp(pArray, n_CH) ;
factor = (min(CH,max(C0,CS)) - CS)/(CH - CS) ;
//coefficient = factor*fka_H+(1-factor)*fka_S
coefficient = fka_H ;
result1 = coefficient ;

return result1;
}

//################Script##################################
void f_f_rate_coeff_a1(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2)
{
double f_oh_a = 0;
double deltah = 0;
double f_rate_coeff = 0;
double time_mix = 0;
double time = 0;
double result3 = 0;
double result4 = 0;

// This routine calculates the rate coefficient for the reaction of the end
// group A(1) with the isocyanate end group I
f_oh_a = getkp(pArray, n_f_oh_a1) ;
deltah = getkp(pArray, n_deltah_a1) ;
f_rate_coeff = f_f_rate_coefficient(xArray, pArray,f_oh_a,0) ;
time_mix = getkp(pArray, n_t_phase_transfer_A) ;
time = tGlobal ;
if (time < time_mix) 
{
f_rate_coeff = 0. ;
}
result1 = f_rate_coeff ;
result2 = 0 ;
result3 = deltah ;
result4 = 0 ;

}

//################Script##################################
void f_f_rate_coeff_b2(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2)
{
double f_oh_b2 = 0;
double deltah = 0;
double f_rate_coeff = 0;
double result3 = 0;
double result4 = 0;

// This routine calculates the rate coefficient for the reaction of the end
// group B with the isocyanate end group I
f_oh_b2 = getkp(pArray, n_f_oh_b2) ;
deltah = getkp(pArray, n_deltah_b2) ;
f_rate_coeff = f_f_rate_coefficient(xArray, pArray,f_oh_b2,0) ;
result1 = f_rate_coeff ;
result2 = 0 ;
result3 = deltah ;
result4 = 0 ;

}

//################Script##################################
void f_f_rate_coeff_b1(double* xArray, double* pArray, double arg1, double arg2, double& result1, double& result2)
{
double f_oh_b = 0;
double deltah = 0;
double f_rate_coeff = 0;
double time_mix = 0;
double time = 0;
double result3 = 0;
double result4 = 0;

// This routine calculates the rate coefficient for the reaction of the end
// group B with the isocyanate end group I
f_oh_b = getkp(pArray, n_f_oh_b) ;
deltah = getkp(pArray, n_deltah_b1) ;
f_rate_coeff = f_f_rate_coefficient(xArray, pArray,f_oh_b,0) ;
time_mix = getkp(pArray, n_t_mix) ;
time = tGlobal ;
if (time < time_mix) 
{
f_rate_coeff = 0.01 * f_rate_coeff ;
}
result1 = f_rate_coeff ;
result2 = 0 ;
result3 = deltah ;
result4 = 0 ;

}

bool F(double* x, double t, double* fx)
{
  tGlobal = t;
double* p = pGlobal;
double k1 = 0;
double k2 = 0;
// Assignments for low molecular species
	double x_Catalyst_1 = x[n_Catalyst_1];
	double x_CE_A0 = x[n_CE_A0];
	double x_CE_A1 = x[n_CE_A1];
	double x_CE_B = x[n_CE_B];
	double x_CE_B2 = x[n_CE_B2];
	double x_CE_I0 = x[n_CE_I0];
	double x_CE_I1 = x[n_CE_I1];
	double x_CE_I2 = x[n_CE_I2];
	double x_CE_PBA = x[n_CE_PBA];
	double x_CE_Breac = x[n_CE_Breac];
	double x_CE_Areac0 = x[n_CE_Areac0];
	double x_CE_Areac1 = x[n_CE_Areac1];
	double x_CE_Ireac0 = x[n_CE_Ireac0];
	double x_CE_Ireac1 = x[n_CE_Ireac1];
	double x_CE_Ireac2 = x[n_CE_Ireac2];
	double x_Bulk = x[n_Bulk];
// Assignments for polymer species
// Assignments for coefficients
	double RGas = 1; // Activation energies normalized
// as soon as RGas is defined you can evaluate the parameter:
double T = gettemp(x, 1);
SetFixedIniParam(p, pControl, T, RGas);
	double kzero = getkp(p, n_kzero);
	double kone = getkp(p, n_kone);
	double kdummy = getkp(p, n_kdummy);
	int ii;
	for (ii = 0; ii<20; ii++) fx[ii] = 0;
// Reaction steps
// Kinetic
// kone*CE_A0^kone+kone*CE_I0^kone<-->kone*CE_Areac0^kone+kone*CE_Ireac0^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a0, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a0(x, p, k1, k2, k1, k2);
	fx[n_CE_A0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Areac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
// Kinetic
// kone*CE_A1^kone+kone*CE_I0^kone<-->kone*CE_Areac1^kone+kone*CE_Ireac0^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a1(x, p, k1, k2, k1, k2);
	fx[n_CE_A1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Areac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
// Kinetic
// kone*CE_A0^kone+kone*CE_I1^kone<-->kone*CE_Areac0^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a0, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a0(x, p, k1, k2, k1, k2);
	fx[n_CE_A0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Areac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_A1^kone+kone*CE_I1^kone<-->kone*CE_Areac1^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a1(x, p, k1, k2, k1, k2);
	fx[n_CE_A1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Areac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_A0^kone+kone*CE_I2^kone<-->kone*CE_Areac0^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a0, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a0(x, p, k1, k2, k1, k2);
	fx[n_CE_A0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Areac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
// Kinetic
// kone*CE_A1^kone+kone*CE_I2^kone<-->kone*CE_Areac1^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a1(x, p, k1, k2, k1, k2);
	fx[n_CE_A1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Areac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
// Kinetic
// kone*CE_B^kone+kone*CE_I0^kone<-->kone*CE_Breac^kone+kone*CE_Ireac0^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b2, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b2(x, p, k1, k2, k1, k2);
	fx[n_CE_B] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
// Kinetic
// kone*CE_B^kone+kone*CE_I2^kone<-->kone*CE_Breac^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b2, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b2(x, p, k1, k2, k1, k2);
	fx[n_CE_B] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
// Kinetic
// kone*CE_B^kone+kone*CE_I1^kone<-->kone*CE_Breac^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b2, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b2(x, p, k1, k2, k1, k2);
	fx[n_CE_B] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_B2^kone+kone*CE_I0^kone<-->kone*CE_B^kone+kone*CE_Ireac0^kone+kone*CE_Breac^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b1(x, p, k1, k2, k1, k2);
	fx[n_CE_B2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_B] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
// Kinetic
// kone*CE_B2^kone+kone*CE_I1^kone<-->kone*CE_B^kone+kone*CE_Breac^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b1(x, p, k1, k2, k1, k2);
	fx[n_CE_B2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_B] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_B2^kone+kone*CE_I2^kone<-->kone*CE_B^kone+kone*CE_Breac^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b1(x, p, k1, k2, k1, k2);
	fx[n_CE_B2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_B] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
return true; 
}
bool A(double* x, double t, double* fx) // diagonal derivatives in vector fx
{
  tGlobal = t;
double* p = pGlobal;
double k1 = 0;
double k2 = 0;
// Assignments for low molecular species
	double x_Catalyst_1 = x[n_Catalyst_1];
	double x_CE_A0 = x[n_CE_A0];
	double x_CE_A1 = x[n_CE_A1];
	double x_CE_B = x[n_CE_B];
	double x_CE_B2 = x[n_CE_B2];
	double x_CE_I0 = x[n_CE_I0];
	double x_CE_I1 = x[n_CE_I1];
	double x_CE_I2 = x[n_CE_I2];
	double x_CE_PBA = x[n_CE_PBA];
	double x_CE_Breac = x[n_CE_Breac];
	double x_CE_Areac0 = x[n_CE_Areac0];
	double x_CE_Areac1 = x[n_CE_Areac1];
	double x_CE_Ireac0 = x[n_CE_Ireac0];
	double x_CE_Ireac1 = x[n_CE_Ireac1];
	double x_CE_Ireac2 = x[n_CE_Ireac2];
	double x_Bulk = x[n_Bulk];
// Assignments for polymer species
// Assignments for coefficients
	double RGas = 1; // Activation energies normalized
// as soon as RGas is defined you can evaluate the parameter:
double T = gettemp(x, 1);
SetFixedIniParam(p, pControl, T, RGas);
	double kzero = getkp(p, n_kzero);
	double kone = getkp(p, n_kone);
	double kdummy = getkp(p, n_kdummy);
	int ii;
	for (ii = 0; ii<20; ii++) fx[ii] = 0;
// Reaction steps
// Kinetic
// kone*CE_A0^kone+kone*CE_I0^kone<-->kone*CE_Areac0^kone+kone*CE_Ireac0^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a0, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a0(x, p, k1, k2, k1, k2);
	fx[n_CE_A0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Areac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac0, kone);
// Kinetic
// kone*CE_A1^kone+kone*CE_I0^kone<-->kone*CE_Areac1^kone+kone*CE_Ireac0^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a1(x, p, k1, k2, k1, k2);
	fx[n_CE_A1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Areac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac0, kone);
// Kinetic
// kone*CE_A0^kone+kone*CE_I1^kone<-->kone*CE_Areac0^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a0, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a0(x, p, k1, k2, k1, k2);
	fx[n_CE_A0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Areac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_A1^kone+kone*CE_I1^kone<-->kone*CE_Areac1^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a1(x, p, k1, k2, k1, k2);
	fx[n_CE_A1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Areac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_A0^kone+kone*CE_I2^kone<-->kone*CE_Areac0^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a0, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a0(x, p, k1, k2, k1, k2);
	fx[n_CE_A0] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Areac0] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_A0, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac0, kone) * pow(  x_CE_Ireac2, kone);
// Kinetic
// kone*CE_A1^kone+kone*CE_I2^kone<-->kone*CE_Areac1^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_a1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_a1(x, p, k1, k2, k1, k2);
	fx[n_CE_A1] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Areac1] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_A1, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Areac1, kone) * pow(  x_CE_Ireac2, kone);
// Kinetic
// kone*CE_B^kone+kone*CE_I0^kone<-->kone*CE_Breac^kone+kone*CE_Ireac0^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b2, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b2(x, p, k1, k2, k1, k2);
	fx[n_CE_B] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac0, kone);
// Kinetic
// kone*CE_B^kone+kone*CE_I2^kone<-->kone*CE_Breac^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b2, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b2(x, p, k1, k2, k1, k2);
	fx[n_CE_B] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
// Kinetic
// kone*CE_B^kone+kone*CE_I1^kone<-->kone*CE_Breac^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b2, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b2(x, p, k1, k2, k1, k2);
	fx[n_CE_B] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_B, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_B2^kone+kone*CE_I0^kone<-->kone*CE_B^kone+kone*CE_Ireac0^kone+kone*CE_Breac^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b1(x, p, k1, k2, k1, k2);
	fx[n_CE_B2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_I0] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_B] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_Ireac0] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I0, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Ireac0, kone) * pow(  x_CE_Breac, kone);
// Kinetic
// kone*CE_B2^kone+kone*CE_I1^kone<-->kone*CE_B^kone+kone*CE_Breac^kone+kone*CE_Ireac1^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b1(x, p, k1, k2, k1, k2);
	fx[n_CE_B2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_I1] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_B] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
	fx[n_CE_Ireac1] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I1, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac1, kone);
// Kinetic
// kone*CE_B2^kone+kone*CE_I2^kone<-->kone*CE_B^kone+kone*CE_Breac^kone+kone*CE_Ireac2^kone,  kdummy,  kzero
//	double EvalFunction(f_rate_coeff_b1, k1, k2)		<------
k1 = kone;
k2 = kone;
f_f_rate_coeff_b1(x, p, k1, k2, k1, k2);
	fx[n_CE_B2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_I2] += -k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) +k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_B] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Breac] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
	fx[n_CE_Ireac2] += +k1 * kone * pow(  x_CE_B2, kone) * pow(  x_CE_I2, kone) -k2 * kone * pow(  x_CE_B, kone) * pow(  x_CE_Breac, kone) * pow(  x_CE_Ireac2, kone);
return true; 
}
