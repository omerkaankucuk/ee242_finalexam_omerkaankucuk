#include "system_output.h"

extern InterruptIn signal1;
extern InterruptIn signal2;
extern Timer encoder_timer;
extern AnalogIn ADCin;

output_struct system_output = {
.read_old = 0,
.read_new = 0,
.time_diff = 0
};

int motor_direction_estimate = 0;
int speed_zero_cnt = 0;

void isr_signal1rise(void){
    system_output.read_new = encoder_timer.read_high_resolution_us();
    system_output.time_diff = system_output.read_new - system_output.read_old;
    system_output.read_old = system_output.read_new;
    if(signal2.read() == 1) motor_direction_estimate = 0;
    else motor_direction_estimate = 1;
    speed_zero_cnt = 0;
}

void isr_signal1fall(void){
    if(signal2.read() == 0) motor_direction_estimate = 0;
    else motor_direction_estimate = 1;
}

void isr_signal2rise(void){
    if(signal1.read() == 0) motor_direction_estimate = 0;
    else motor_direction_estimate = 1;
}

void isr_signal2fall(void){
    if(signal1.read() == 1) motor_direction_estimate = 0;
    else motor_direction_estimate = 1;
}

void system_output_init(int type, float gain){
    system_output.flag = 1;
    system_output.gain = gain;
    if(type == 0){
        system_output.type = 0;
    }
    else if(type == 1){
        system_output.type = 1;
        signal1.rise(&isr_signal1rise); 
        signal1.fall(&isr_signal1fall); 
        signal2.rise(&isr_signal2rise); 
        signal2.fall(&isr_signal2fall);
        encoder_timer.start(); 
    }
}

float get_ADC(void){
    float ADC_output;
    ADC_output = ADCin.read() * 3.3 * system_output.gain;
    return ADC_output;  
}

float calculate_frequency(int cycle_difference){
    float freq;
    if(cycle_difference > 0)freq = ((float)1000000 / ((float)(cycle_difference)));
    else freq = 0;
    return freq;
}

float calculate_speed(float freq){
    float motor_speed;
    if(speed_zero_cnt < 10){
        if(motor_direction_estimate == 1) motor_speed = 4*freq*60 / system_output.gain; 
        else motor_speed = -4*freq*60 / system_output.gain; 
        speed_zero_cnt++;
    }
    else {
        motor_speed = 0;   
    }
    return motor_speed;
}

float obtain_system_output(void){
    float output_voltage;
    if(system_output.type == 0){
        output_voltage = get_ADC();
    }
    else if(system_output.type == 1){
        output_voltage = calculate_speed(calculate_frequency(system_output.time_diff));
    }
    return output_voltage;
}





