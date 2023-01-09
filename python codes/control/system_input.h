#ifndef SYSTEM_INPUT_H
#define SYSTEM_INPUT_H

#include "mbed.h"
#include "menu.h"

typedef struct {
    float gain;
    bool flag;
    int type;
} input_struct;

extern input_struct system_input;

void system_input_init(int type, float gain);
float convert_system_input_voltage(float voltage);
void set_DAC(float voltage_mcu);
void set_duty_cycles(float voltage_mcu);
void generate_system_input(float voltage);

#endif /* SYSTEM_INPUT_H */
