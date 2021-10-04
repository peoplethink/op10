#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5370574731188966315);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4122355017533011283);
void car_H_mod_fun(double *state, double *out_5400649318337851136);
void car_f_fun(double *state, double dt, double *out_8379005402228353804);
void car_F_fun(double *state, double dt, double *out_4099916354956946501);
void car_h_25(double *state, double *unused, double *out_6190470910566514058);
void car_H_25(double *state, double *unused, double *out_2914461629953627239);
void car_h_24(double *state, double *unused, double *out_3299126626044070150);
void car_H_24(double *state, double *unused, double *out_4127584231619464799);
void car_h_30(double *state, double *unused, double *out_5485437936792686622);
void car_H_30(double *state, double *unused, double *out_6699437280995556041);
void car_h_26(double *state, double *unused, double *out_6653181737868828359);
void car_H_26(double *state, double *unused, double *out_8619054784892836512);
void car_h_27(double *state, double *unused, double *out_5587155025890036541);
void car_H_27(double *state, double *unused, double *out_7987019268832181353);
void car_h_29(double *state, double *unused, double *out_1725471549461534889);
void car_H_29(double *state, double *unused, double *out_7923508023585163622);
void car_h_28(double *state, double *unused, double *out_1446770817468421996);
void car_H_28(double *state, double *unused, double *out_2890248246617030067);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}