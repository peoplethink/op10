#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_803611454938064350);
void live_err_fun(double *nom_x, double *delta_x, double *out_8121643576604858258);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_9110887307405050186);
void live_H_mod_fun(double *state, double *out_6190164851874934599);
void live_f_fun(double *state, double dt, double *out_5982242886378791203);
void live_F_fun(double *state, double dt, double *out_7119455599707917982);
void live_h_3(double *state, double *unused, double *out_2170545346926773647);
void live_H_3(double *state, double *unused, double *out_3947804948072223774);
void live_h_4(double *state, double *unused, double *out_1491348554285626960);
void live_H_4(double *state, double *unused, double *out_3964879641017559692);
void live_h_9(double *state, double *unused, double *out_6977411688525780381);
void live_H_9(double *state, double *unused, double *out_5403655191488153103);
void live_h_10(double *state, double *unused, double *out_3060866906772579977);
void live_H_10(double *state, double *unused, double *out_2377626970716605579);
void live_h_12(double *state, double *unused, double *out_7539190307902695261);
void live_H_12(double *state, double *unused, double *out_1858324129278350123);
void live_h_31(double *state, double *unused, double *out_295030898892930956);
void live_H_31(double *state, double *unused, double *out_135056521437401816);
void live_h_32(double *state, double *unused, double *out_7455682882139262716);
void live_H_32(double *state, double *unused, double *out_5803290872635112154);
void live_h_13(double *state, double *unused, double *out_9059686845966255663);
void live_H_13(double *state, double *unused, double *out_7519097344678673706);
void live_h_14(double *state, double *unused, double *out_6977411688525780381);
void live_H_14(double *state, double *unused, double *out_5403655191488153103);
void live_h_19(double *state, double *unused, double *out_6665703002444968454);
void live_H_19(double *state, double *unused, double *out_8610409031118967332);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}