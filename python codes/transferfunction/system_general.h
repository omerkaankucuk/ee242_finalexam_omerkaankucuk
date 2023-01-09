#ifndef SYSTEM_GENERAL_H
#define SYSTEM_GENERAL_H

#include "mbed.h"
#include "menu.h"
#include "Ccontrol.h"
#include "system_input.h"
#include "system_output.h"
#include "system_controller.h"

typedef union _data {
  float f;
  char  s[4];
} myData;

void isr_button(void);
void user_button_init(void);
void isr_sampling_timer(void);
void sampling_timer_init(void);
void send_data_RT(int select);
void send_variable(float var);
void select_signal(void);
void start_function(void);
void stop_function(void);
void menu_init(void);
void system_loop(void);

#endif /* SYSTEM_GENERAL_H */
