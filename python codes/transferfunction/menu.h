#ifndef MENU_H
#define MENU_H

#include "mbed.h"

#define PWM_freq_initial 20000
#define Sampling_freq_initial 2000

typedef struct {
    float sampling_frequency;
    float PWM_frequency;
    int signal_select;
    float amplitude;
    float frequency;
    float phase;
    float offset;
    int sine_select;
    float amplitude_exp;
    float period;
    float dutycycle;
    int sine_number;
    float *amplitudes;
    float *frequencies;
    float *phases;
    float *offsets;
    float frequency_start;
    float frequency_stop;
    float sweep_duration;
    int random_duration;
    int motor_direction;  
} menu_struct;

extern menu_struct menu;

void get_string_to_buffer(void);
void menu_set_sampling_frequency(void);
void menu_set_PWM_frequency(void);
void menu_input_selection(void);
void constant_amplitude_input(void);
void step_signal_input(void);
void ramp_signal_input(void);
void parabolic_signal_input(void);
void exponential_signal_input(void);
void sinusoidal_signal_input(void);
void damped_exponential_signal_input(void);
void rectangular_signal_input(void);
void sum_of_sinusoids_signal_input(void);
void sweep_signal_input(void);
void random_signal_input(void);
void menu_system_ready(void);
void menu_system_running(void);
void menu_system_stopped(void);
void menu_motor_direction_selection(void);

#endif /* MENU_H */
