#include "system_input.h"

extern PwmOut pwm1;
extern PwmOut pwm2;
extern DigitalOut pwm_enable;
extern AnalogOut DACout;

input_struct system_input;

void system_input_init(int type, float gain){
    system_input.gain = gain;
    system_input.flag = 1;
    if(type == 0){
        system_input.type = 0;
    }
    else if(type == 1){  
        float period_us;     
        system_input.type = 1;
        menu_set_PWM_frequency();
        period_us = 1000000 / menu.PWM_frequency;  
        pwm1.period_us(period_us);
        pwm2.period_us(period_us);  
        pwm1 = 0;
        pwm2 = 0;     
        pwm_enable = 1;    
    }
}

float convert_system_input_voltage(float voltage){
        float voltage_mcu;
        voltage_mcu = voltage / system_input.gain;
        if(voltage_mcu > 3.3)voltage_mcu = 3.3;
        if(voltage_mcu < -3.3)voltage_mcu = -3.3;   
        return voltage_mcu;
}

void set_DAC(float voltage_mcu){
    float DAC_input;
    DAC_input = voltage_mcu/3.3;
    if(DAC_input<0)DAC_input = 0;
    DACout.write(DAC_input);  
}  

void set_duty_cycles(float voltage_mcu){
    float duty;
    duty = voltage_mcu/3.3;
    if(menu.motor_direction == 0){
        if(duty<0){
            pwm1 = 0;
            pwm2 = -duty;               
        }            
        else if(duty==0){
            pwm1 = 0;
            pwm2 = 0;                    
        }
        else{
            pwm1 = duty;
            pwm2 = 0; 
            menu.motor_direction = 1;                   
        }
    }
    else if(menu.motor_direction == 1){
        if(duty>0){
            pwm1 = duty;
            pwm2 = 0;               
        }            
        else if(duty==0){
            pwm1 = 0;
            pwm2 = 0;                    
        }
        else{
            pwm1 = 0;
            pwm2 = -duty; 
            menu.motor_direction = 0;                   
        }
    }   
}

void generate_system_input(float voltage){
    float voltage_mcu;
    voltage_mcu = convert_system_input_voltage(voltage);
    if(system_input.type == 0){ 
        set_DAC(voltage_mcu);
    }
    else if(system_input.type == 1){
        set_duty_cycles(voltage_mcu);
    }
}