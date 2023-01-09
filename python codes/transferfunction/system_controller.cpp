#include "system_controller.h"

controller_struct system_controller;
tf_struct tf_controller;
ss_controller_struct ss_controller;


bool controller_tf_select = 0;
float tsampling = 0;
float *x;
float *y;

float ei = 0;
float ei_old = 0;
float e = 0; 
float u = 0;
float yest = 0;
float yest_old = 0;

float *xest;   
float *xest_old;
float dummy0[1];
float *dummy1;
float *dummy2;

int i, j;

void controller_tf_init(float *num, float *den, int nbr_num_coeffs, int nbr_den_coeffs, int select){
    system_controller.flag = 1;
    system_controller.type = 0;    
    controller_tf_select = select;
    tsampling = 1 / menu.sampling_frequency;
    _create_tf(&tf_controller, num, den, nbr_num_coeffs, nbr_den_coeffs, tsampling);
    x = (float *) malloc(nbr_num_coeffs * sizeof(float));
    y = (float *) malloc(nbr_den_coeffs * sizeof(float));  
}

void controller_ss_init(float *F, float *G, float *C, float *D, float *K, float *L, float *xest0, float k1, int no_states){
    system_controller.flag = 1;
    system_controller.type = 1;  
	ss_controller.F = F;
	ss_controller.G = G;
    ss_controller.C = C;
	ss_controller.D = D;
	ss_controller.no_states = no_states;  
    ss_controller.K = K;
    ss_controller.L = L; 
    ss_controller.k1 = k1;  
    ss_controller.xest0 = xest0;   
	xest = (float *)malloc(no_states * sizeof(float));
	xest_old = (float *)malloc(no_states * sizeof(float)); 
	dummy1 = (float *)malloc(no_states * sizeof(float));  
	dummy2 = (float *)malloc(no_states * sizeof(float));  

    for (i = 0; i < 1; i++) {
        for (j = 0 ; j < no_states; j++) {
            xest_old[i*no_states+j] =xest0[i*no_states+j];
        }
    }
}

float obtain_controller_output(float input, float output, int cnt_signal){
    float controller_output;
    if(system_controller.type == 0){      
        float controller_input;
        if(controller_tf_select == 0){
            controller_input = input;
        }
        else{
            if(cnt_signal == 0)controller_input = input;
            else controller_input = input - output;
        }

        x[cnt_signal % (tf_controller.num_order + 1)] = controller_input;
        controller_output = system_output_per_sample(&tf_controller, x, y, cnt_signal);       
    }
    else if(system_controller.type == 1){      
        if(cnt_signal >= 2){
            e = input - output;
            ei = e + ei_old; 

            mult_matrix(ss_controller.K, xest_old, dummy0, 1, 2, 1);
            u =  ss_controller.k1*ei - dummy0[0];

            mult_matrix(ss_controller.F, xest_old, dummy1, 2, 2, 1);
            mult_matrix_scalar(ss_controller.G, u, dummy2, 2, 1); 
            add_matrices(dummy1, dummy2, dummy2, 2, 1);
            mult_matrix_scalar(ss_controller.L, (output - yest_old), dummy1, 2, 1);
            add_matrices(dummy1, dummy2, xest, 2, 1);

            mult_matrix(ss_controller.C, xest_old, dummy0, 1, 2, 1);
            yest = dummy0[0];

            for (i = 0; i < 1; i++) {
                for (j = 0 ; j < ss_controller.no_states; j++) {
                    xest_old[i*ss_controller.no_states+j] =xest[i*ss_controller.no_states+j];
                }
            }
            system_controller.x1 = xest[0];
            yest_old = yest;
            ei_old = ei;                      
        }
        controller_output = u;      
    }
    return controller_output;
}

RLS_Controller::RLS_Controller(unsigned int na, unsigned int nb, double gain, unsigned int cnt_max) : 
_na(na), _nb(nb), _gain(gain), _cnt_max(cnt_max), _P_old(na+nb+1, na+nb+1), _theta_hat_old(na+nb+1, 1), _idx(na+nb+1, na + nb + 1), P(na+nb+1, na+nb+1), 
phiT(1, na+nb+1), theta_hat(na+nb+1, 1), phiT_trans(na+nb+1, 1), dummy(1, 1), dummy2(na+nb+1, 1), dummy3(na+nb+1, 1), P_trans(na+nb+1, na+nb+1), K(na+nb+1, 1){
    _epsilon = 0;
    _x = new float[_nb + 1];
    _y = new float[_na + 1]; 
    _N = _na + 1;     
    _num = new float[_nb + 1];
    _den = new float[_na + 1]; 
    for (int i = 0; i < _na + _nb + 1; i++){
        for (int j = 0; j < _na + _nb + 1; j++){
            if(i == j){
                _P_old(i, j) = 1;
            }
        }
    }
    _P_old = _P_old * _gain;
    for (int i = 0; i < _na + _nb + 1; i++) {
        for (int j = 0; j < _na + _nb + 1; j++) {
            if (i == j) {
                _idx(i, j) = 1;
            }
        }
    } 
    system_controller.flag = 1;
    system_controller.type = 2;         
}

float RLS_Controller::obtain_controller_output(float xnew, float ynew, int cnt){
    _x[cnt % (_nb+1)] = xnew;
    _y[cnt % (_na+1)] = ynew;
    if(cnt >= _N - 1){
        for(int j = 0; j < _na; j++){
            if((cnt-j)<0){
                phiT(0, j) = 0;
            }
            else {
                int u1 = (cnt % (_na + 1));
                if((u1 - j - 1) < 0){
                    phiT(0, j) = -_y[(cnt % (_na + 1)) - j - 1 + (_na + 1)];
                }
                else {
                    phiT(0, j) = -_y[(cnt % (_na + 1)) - j - 1];
                }
            }
        }
        for(int j = 0; j < (_nb+1); j++){
            if((cnt-j)<0){
                phiT(0, j) = 0;
            }
            else {
                int u2 = (cnt % (_nb + 1));
                if((u2 - j) < 0){
                    phiT(0, (j+_na)) = _x[(cnt % (_nb + 1)) - j + (_nb + 1)];
                }
                else {
                    phiT(0, (j+_na)) = _x[(cnt % (_nb + 1)) - j];
                }
            }
        }
        phiT_trans = phiT.trans();
        dummy = (((phiT * _P_old) * phiT_trans) + 1).inv();
        K = (_P_old * phiT_trans) * dummy;
        dummy = phiT * _theta_hat_old;
        _epsilon = ynew - dummy(0, 0);
        dummy2 = K * _epsilon;
        theta_hat = _theta_hat_old + dummy2;
        dummy3 =  ((K * phiT) * -1);
        P = (_idx + dummy3) * _P_old;
        P_trans = P.trans();
        P = (P + P_trans) * 0.5;
        _P_old = P;
        _theta_hat_old = theta_hat;
    } 
    if(cnt < _cnt_max - 1){
        return -1;
    }
    else{
        _den[0] = 1;
        for (int j = 0; j < _na+_nb+1; j++) {
            if (j < _na) {
                _den[j+1] = theta_hat(j, 0);
            }
            else {
                _num[j - _na] = theta_hat(j, 0);
            }
        }
        return 1;
    }               
}