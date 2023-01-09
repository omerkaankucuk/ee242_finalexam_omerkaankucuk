#ifndef CCONTROL_H
#define CCONTROL_H

#include "mbed.h"

#define PI 3.14159265358979323846

typedef struct {
	float *num; 			// Pointer for numerator array
	float *den;			// Pointer for denerator array
	int num_order;			// Numerator order
	int den_order;		// denerator order
	float sampling_time;	// Sampling time of the system
	float poles_real[20]; 	// Real parts of the system poles
	float poles_imag[20]; 	// Imaginary parts of the system poles
	float zeros_real[20]; 	// Real parts of the system zeros
	float zeros_imag[20]; 	// Imaginary parts of the system zeros
	float poles_mags[20]; 	// Pole magnitudes
	float y_final; 			// Final steady state value of the system output
	float rise_time; 		// Rise time of the system
	float settling_time;	// Settling time of the system
	float settling_min;		// Minimum value of the system output after system response is risen
	float settling_max;		// Maximum value of the system output after system response is risen
	float overshoot;		// Percentage of the overshoot relative to y_final
	float undershoot;		// Percentage of the undershoot relative to y_final
	float peak;				// Absolute peak value of the system output
	float peak_time;		// Time at which peak value occurs
} tf_struct;

typedef struct {
	float *F; 
	float *G;	
	float *C; 
	float *D;	
	int no_states;
    float sampling_time;
} ss_struct;



void create_tf(tf_struct *tf, float *num, float *den, int nbr_num_coeffs, int nbr_den_coeffs, float ts);
extern float unit_pulse_sample(int cnt);
void unit_pulse_signal(int N, float signal[]);
extern float step_sample(int cnt, float Amp);
void step_signal(int N, float Amp, float signal[]);
extern float ramp_sample(int cnt, float Amp);
void ramp_signal(int N, float Amp, float signal[]);
extern float parabolic_sample(int cnt, float Amp);
void parabolic_signal(int N, float Amp, float signal[]);
extern float exponential_sample(int cnt, float a);
void exponential_signal(int N, float a, float signal[]);
extern float sinusoidal_sample(int cnt, float Amp, float Freq, float Phase, float Offset, float Fsample, int select);
void sinusoidal_signal(int N, float Amp, float Freq, float Phase, float Offset, float Fsample, int select, float signal[]);
float damped_sinusoidal_sample(int cnt, float Amp_exp, float Amp_sin, float Freq, float Phase, float Offset, float Fsample, int select);
void damped_sinusoidal_signal(int N, float Amp_exp, float Amp_sin, float Freq, float Phase, float Offset, float Fsample, int select, float signal[]);
float rectangular_sample(int cnt, float Amp, float Period, float Duty, float Offset, float Fsample);
void rectangular_signal(int N, float Amp, float Period, float Duty, float Offset, float Fsample, float signal[]);
float sum_of_sinusoids_sample(int cnt, int No_Sines, float Amps[], float Freqs[], float Phases[], float Offsets[], float Fsample, int select);
void sum_of_sinusoids_signal(int N, int No_Sines, float Amps[], float Freqs[], float Phases[], float Offsets[], float Fsample, int select, float signal[]);
float sweep_sample(int cnt, float Amp, float Freq_start, float Freq_stop, float Offset, float Fsample, float duration);
void sweep_signal(int N, float Amp, float Freq_start, float Freq_stop, float Offset, float Fsample, float duration, float signal[]);
float random_sample(int cnt, float Amp, int Duration, float Offset);
void random_signal(int N, float Amp, int Duration, float Offset, float signal[]);
void _iir_filter(float b[], float a[], float x[], float y[], int size1, int size2, int size3);
void lsim(tf_struct *tf,int N, float x[], float y[]);

void send_data(float signal[],int L);

void _print_stars(int star_number);
void _print_nextline(void);
void _create_tf(tf_struct *tf, float *num, float *den, int nbr_num_coeffs, int nbr_den_coeffs, float ts);
float system_output_per_sample(tf_struct *tf, float x[], float y[], int n);

void add_matrices(float *x1, float *x2, float *result, int size_row, int size_column);
void sub_matrices(float *x1, float *x2, float *result, int size_row, int size_column);
void transpose(float *x, float *result, int size_row, int size_column);
void mult_matrix_scalar(float *x, float A, float *result, int size_row, int size_column);
void mult_matrix(float *x1, float *x2, float *result, int size_row_x1, int size_column_x1, int size_column_x2);
void create_ss(ss_struct *ss, float *F, float *G, float *C, float *D, int no_states, float ts);
void lsim_ss(ss_struct *ss, float *u, float *x0, float *result, int select, int N);
#endif /* CCONTROL_H */