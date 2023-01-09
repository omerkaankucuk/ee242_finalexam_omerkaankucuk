#ifndef SYSTEM_CONTROLLER_H
#define SYSTEM_CONTROLLER_H

#include "mbed.h"
#include "menu.h"
#include "Ccontrol.h"
#include "matrix.h"

typedef struct {
    bool flag;
    int type;
    float x1;
} controller_struct;

typedef struct {
	float *F; 
	float *G;	
	float *C; 
	float *D;	
	int no_states;   
	float *K; 
	float *L;	
	float k1;
    float *xest0; 	
} ss_controller_struct;

extern controller_struct system_controller;

class RLS_Controller
{
public:
    RLS_Controller(unsigned int na, unsigned int nb, double gain, unsigned int cnt_max);
    float obtain_controller_output(float xnew, float ynew, int cnt);
    float *_num;
    float *_den;
    unsigned int _na;
    unsigned int _nb;
private:
	double _gain;
	unsigned int _cnt_max;  
	Matrix _P_old;
	Matrix _theta_hat_old;
	Matrix _idx;
	double _epsilon;
    float *_x;
	float *_y;
	unsigned int _N; 
	Matrix P;
	Matrix phiT;
	Matrix theta_hat;
	Matrix phiT_trans;
	Matrix dummy;
	Matrix dummy2;
	Matrix dummy3;
	Matrix P_trans;
	Matrix K;       
};

void controller_tf_init(float *num, float *den, int nbr_num_coeffs, int nbr_den_coeffs, int select);
void controller_ss_init(float *F, float *G, float *C, float *D, float *K, float *L, float *xest0, float k1, int no_states);
float obtain_controller_output(float input, float output, int cnt_signal);

#endif /* SYSTEM_CONTROLLER_H */
