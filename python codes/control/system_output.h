#ifndef SYSTEM_OUTPUT_H
#define SYSTEM_OUTPUT_H

#include "mbed.h"

typedef struct {
    float gain;
    bool flag;
    int type;
    uint64_t read_old;
    uint64_t read_new;
    uint64_t time_diff;
} output_struct;

extern output_struct system_output;

void isr_signal1rise(void);
void isr_signal1fall(void);
void isr_signal2rise(void);
void isr_signal2fall(void);
void system_output_init(int type, float gain);
float get_ADC(void);
float calculate_frequency(int cycle_difference);
float calculate_speed(float freq);
float obtain_system_output(void);

#endif /* SYSTEM_OUTPUT_H */
