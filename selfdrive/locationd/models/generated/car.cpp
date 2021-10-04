#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5370574731188966315) {
   out_5370574731188966315[0] = delta_x[0] + nom_x[0];
   out_5370574731188966315[1] = delta_x[1] + nom_x[1];
   out_5370574731188966315[2] = delta_x[2] + nom_x[2];
   out_5370574731188966315[3] = delta_x[3] + nom_x[3];
   out_5370574731188966315[4] = delta_x[4] + nom_x[4];
   out_5370574731188966315[5] = delta_x[5] + nom_x[5];
   out_5370574731188966315[6] = delta_x[6] + nom_x[6];
   out_5370574731188966315[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4122355017533011283) {
   out_4122355017533011283[0] = -nom_x[0] + true_x[0];
   out_4122355017533011283[1] = -nom_x[1] + true_x[1];
   out_4122355017533011283[2] = -nom_x[2] + true_x[2];
   out_4122355017533011283[3] = -nom_x[3] + true_x[3];
   out_4122355017533011283[4] = -nom_x[4] + true_x[4];
   out_4122355017533011283[5] = -nom_x[5] + true_x[5];
   out_4122355017533011283[6] = -nom_x[6] + true_x[6];
   out_4122355017533011283[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5400649318337851136) {
   out_5400649318337851136[0] = 1.0;
   out_5400649318337851136[1] = 0.0;
   out_5400649318337851136[2] = 0.0;
   out_5400649318337851136[3] = 0.0;
   out_5400649318337851136[4] = 0.0;
   out_5400649318337851136[5] = 0.0;
   out_5400649318337851136[6] = 0.0;
   out_5400649318337851136[7] = 0.0;
   out_5400649318337851136[8] = 0.0;
   out_5400649318337851136[9] = 1.0;
   out_5400649318337851136[10] = 0.0;
   out_5400649318337851136[11] = 0.0;
   out_5400649318337851136[12] = 0.0;
   out_5400649318337851136[13] = 0.0;
   out_5400649318337851136[14] = 0.0;
   out_5400649318337851136[15] = 0.0;
   out_5400649318337851136[16] = 0.0;
   out_5400649318337851136[17] = 0.0;
   out_5400649318337851136[18] = 1.0;
   out_5400649318337851136[19] = 0.0;
   out_5400649318337851136[20] = 0.0;
   out_5400649318337851136[21] = 0.0;
   out_5400649318337851136[22] = 0.0;
   out_5400649318337851136[23] = 0.0;
   out_5400649318337851136[24] = 0.0;
   out_5400649318337851136[25] = 0.0;
   out_5400649318337851136[26] = 0.0;
   out_5400649318337851136[27] = 1.0;
   out_5400649318337851136[28] = 0.0;
   out_5400649318337851136[29] = 0.0;
   out_5400649318337851136[30] = 0.0;
   out_5400649318337851136[31] = 0.0;
   out_5400649318337851136[32] = 0.0;
   out_5400649318337851136[33] = 0.0;
   out_5400649318337851136[34] = 0.0;
   out_5400649318337851136[35] = 0.0;
   out_5400649318337851136[36] = 1.0;
   out_5400649318337851136[37] = 0.0;
   out_5400649318337851136[38] = 0.0;
   out_5400649318337851136[39] = 0.0;
   out_5400649318337851136[40] = 0.0;
   out_5400649318337851136[41] = 0.0;
   out_5400649318337851136[42] = 0.0;
   out_5400649318337851136[43] = 0.0;
   out_5400649318337851136[44] = 0.0;
   out_5400649318337851136[45] = 1.0;
   out_5400649318337851136[46] = 0.0;
   out_5400649318337851136[47] = 0.0;
   out_5400649318337851136[48] = 0.0;
   out_5400649318337851136[49] = 0.0;
   out_5400649318337851136[50] = 0.0;
   out_5400649318337851136[51] = 0.0;
   out_5400649318337851136[52] = 0.0;
   out_5400649318337851136[53] = 0.0;
   out_5400649318337851136[54] = 1.0;
   out_5400649318337851136[55] = 0.0;
   out_5400649318337851136[56] = 0.0;
   out_5400649318337851136[57] = 0.0;
   out_5400649318337851136[58] = 0.0;
   out_5400649318337851136[59] = 0.0;
   out_5400649318337851136[60] = 0.0;
   out_5400649318337851136[61] = 0.0;
   out_5400649318337851136[62] = 0.0;
   out_5400649318337851136[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8379005402228353804) {
   out_8379005402228353804[0] = state[0];
   out_8379005402228353804[1] = state[1];
   out_8379005402228353804[2] = state[2];
   out_8379005402228353804[3] = state[3];
   out_8379005402228353804[4] = state[4];
   out_8379005402228353804[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8379005402228353804[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8379005402228353804[7] = state[7];
}
void F_fun(double *state, double dt, double *out_4099916354956946501) {
   out_4099916354956946501[0] = 1;
   out_4099916354956946501[1] = 0;
   out_4099916354956946501[2] = 0;
   out_4099916354956946501[3] = 0;
   out_4099916354956946501[4] = 0;
   out_4099916354956946501[5] = 0;
   out_4099916354956946501[6] = 0;
   out_4099916354956946501[7] = 0;
   out_4099916354956946501[8] = 0;
   out_4099916354956946501[9] = 1;
   out_4099916354956946501[10] = 0;
   out_4099916354956946501[11] = 0;
   out_4099916354956946501[12] = 0;
   out_4099916354956946501[13] = 0;
   out_4099916354956946501[14] = 0;
   out_4099916354956946501[15] = 0;
   out_4099916354956946501[16] = 0;
   out_4099916354956946501[17] = 0;
   out_4099916354956946501[18] = 1;
   out_4099916354956946501[19] = 0;
   out_4099916354956946501[20] = 0;
   out_4099916354956946501[21] = 0;
   out_4099916354956946501[22] = 0;
   out_4099916354956946501[23] = 0;
   out_4099916354956946501[24] = 0;
   out_4099916354956946501[25] = 0;
   out_4099916354956946501[26] = 0;
   out_4099916354956946501[27] = 1;
   out_4099916354956946501[28] = 0;
   out_4099916354956946501[29] = 0;
   out_4099916354956946501[30] = 0;
   out_4099916354956946501[31] = 0;
   out_4099916354956946501[32] = 0;
   out_4099916354956946501[33] = 0;
   out_4099916354956946501[34] = 0;
   out_4099916354956946501[35] = 0;
   out_4099916354956946501[36] = 1;
   out_4099916354956946501[37] = 0;
   out_4099916354956946501[38] = 0;
   out_4099916354956946501[39] = 0;
   out_4099916354956946501[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4099916354956946501[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4099916354956946501[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4099916354956946501[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4099916354956946501[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4099916354956946501[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4099916354956946501[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4099916354956946501[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4099916354956946501[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4099916354956946501[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4099916354956946501[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4099916354956946501[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4099916354956946501[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4099916354956946501[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4099916354956946501[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4099916354956946501[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4099916354956946501[56] = 0;
   out_4099916354956946501[57] = 0;
   out_4099916354956946501[58] = 0;
   out_4099916354956946501[59] = 0;
   out_4099916354956946501[60] = 0;
   out_4099916354956946501[61] = 0;
   out_4099916354956946501[62] = 0;
   out_4099916354956946501[63] = 1;
}
void h_25(double *state, double *unused, double *out_6190470910566514058) {
   out_6190470910566514058[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2914461629953627239) {
   out_2914461629953627239[0] = 0;
   out_2914461629953627239[1] = 0;
   out_2914461629953627239[2] = 0;
   out_2914461629953627239[3] = 0;
   out_2914461629953627239[4] = 0;
   out_2914461629953627239[5] = 0;
   out_2914461629953627239[6] = 1;
   out_2914461629953627239[7] = 0;
}
void h_24(double *state, double *unused, double *out_3299126626044070150) {
   out_3299126626044070150[0] = state[4];
   out_3299126626044070150[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4127584231619464799) {
   out_4127584231619464799[0] = 0;
   out_4127584231619464799[1] = 0;
   out_4127584231619464799[2] = 0;
   out_4127584231619464799[3] = 0;
   out_4127584231619464799[4] = 1;
   out_4127584231619464799[5] = 0;
   out_4127584231619464799[6] = 0;
   out_4127584231619464799[7] = 0;
   out_4127584231619464799[8] = 0;
   out_4127584231619464799[9] = 0;
   out_4127584231619464799[10] = 0;
   out_4127584231619464799[11] = 0;
   out_4127584231619464799[12] = 0;
   out_4127584231619464799[13] = 1;
   out_4127584231619464799[14] = 0;
   out_4127584231619464799[15] = 0;
}
void h_30(double *state, double *unused, double *out_5485437936792686622) {
   out_5485437936792686622[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6699437280995556041) {
   out_6699437280995556041[0] = 0;
   out_6699437280995556041[1] = 0;
   out_6699437280995556041[2] = 0;
   out_6699437280995556041[3] = 0;
   out_6699437280995556041[4] = 1;
   out_6699437280995556041[5] = 0;
   out_6699437280995556041[6] = 0;
   out_6699437280995556041[7] = 0;
}
void h_26(double *state, double *unused, double *out_6653181737868828359) {
   out_6653181737868828359[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8619054784892836512) {
   out_8619054784892836512[0] = 0;
   out_8619054784892836512[1] = 0;
   out_8619054784892836512[2] = 0;
   out_8619054784892836512[3] = 0;
   out_8619054784892836512[4] = 0;
   out_8619054784892836512[5] = 0;
   out_8619054784892836512[6] = 0;
   out_8619054784892836512[7] = 1;
}
void h_27(double *state, double *unused, double *out_5587155025890036541) {
   out_5587155025890036541[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7987019268832181353) {
   out_7987019268832181353[0] = 0;
   out_7987019268832181353[1] = 0;
   out_7987019268832181353[2] = 0;
   out_7987019268832181353[3] = 1;
   out_7987019268832181353[4] = 0;
   out_7987019268832181353[5] = 0;
   out_7987019268832181353[6] = 0;
   out_7987019268832181353[7] = 0;
}
void h_29(double *state, double *unused, double *out_1725471549461534889) {
   out_1725471549461534889[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7923508023585163622) {
   out_7923508023585163622[0] = 0;
   out_7923508023585163622[1] = 1;
   out_7923508023585163622[2] = 0;
   out_7923508023585163622[3] = 0;
   out_7923508023585163622[4] = 0;
   out_7923508023585163622[5] = 0;
   out_7923508023585163622[6] = 0;
   out_7923508023585163622[7] = 0;
}
void h_28(double *state, double *unused, double *out_1446770817468421996) {
   out_1446770817468421996[0] = state[5];
   out_1446770817468421996[1] = state[6];
}
void H_28(double *state, double *unused, double *out_2890248246617030067) {
   out_2890248246617030067[0] = 0;
   out_2890248246617030067[1] = 0;
   out_2890248246617030067[2] = 0;
   out_2890248246617030067[3] = 0;
   out_2890248246617030067[4] = 0;
   out_2890248246617030067[5] = 1;
   out_2890248246617030067[6] = 0;
   out_2890248246617030067[7] = 0;
   out_2890248246617030067[8] = 0;
   out_2890248246617030067[9] = 0;
   out_2890248246617030067[10] = 0;
   out_2890248246617030067[11] = 0;
   out_2890248246617030067[12] = 0;
   out_2890248246617030067[13] = 0;
   out_2890248246617030067[14] = 1;
   out_2890248246617030067[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5370574731188966315) {
  err_fun(nom_x, delta_x, out_5370574731188966315);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4122355017533011283) {
  inv_err_fun(nom_x, true_x, out_4122355017533011283);
}
void car_H_mod_fun(double *state, double *out_5400649318337851136) {
  H_mod_fun(state, out_5400649318337851136);
}
void car_f_fun(double *state, double dt, double *out_8379005402228353804) {
  f_fun(state,  dt, out_8379005402228353804);
}
void car_F_fun(double *state, double dt, double *out_4099916354956946501) {
  F_fun(state,  dt, out_4099916354956946501);
}
void car_h_25(double *state, double *unused, double *out_6190470910566514058) {
  h_25(state, unused, out_6190470910566514058);
}
void car_H_25(double *state, double *unused, double *out_2914461629953627239) {
  H_25(state, unused, out_2914461629953627239);
}
void car_h_24(double *state, double *unused, double *out_3299126626044070150) {
  h_24(state, unused, out_3299126626044070150);
}
void car_H_24(double *state, double *unused, double *out_4127584231619464799) {
  H_24(state, unused, out_4127584231619464799);
}
void car_h_30(double *state, double *unused, double *out_5485437936792686622) {
  h_30(state, unused, out_5485437936792686622);
}
void car_H_30(double *state, double *unused, double *out_6699437280995556041) {
  H_30(state, unused, out_6699437280995556041);
}
void car_h_26(double *state, double *unused, double *out_6653181737868828359) {
  h_26(state, unused, out_6653181737868828359);
}
void car_H_26(double *state, double *unused, double *out_8619054784892836512) {
  H_26(state, unused, out_8619054784892836512);
}
void car_h_27(double *state, double *unused, double *out_5587155025890036541) {
  h_27(state, unused, out_5587155025890036541);
}
void car_H_27(double *state, double *unused, double *out_7987019268832181353) {
  H_27(state, unused, out_7987019268832181353);
}
void car_h_29(double *state, double *unused, double *out_1725471549461534889) {
  h_29(state, unused, out_1725471549461534889);
}
void car_H_29(double *state, double *unused, double *out_7923508023585163622) {
  H_29(state, unused, out_7923508023585163622);
}
void car_h_28(double *state, double *unused, double *out_1446770817468421996) {
  h_28(state, unused, out_1446770817468421996);
}
void car_H_28(double *state, double *unused, double *out_2890248246617030067) {
  H_28(state, unused, out_2890248246617030067);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
