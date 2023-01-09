#include "mbed.h"
#include "system_general.h"
#include "system_input.h"
#include "system_output.h"
#include "system_controller.h"

BufferedSerial pc(USBTX, USBRX, 9600);
UnbufferedSerial device(PD_5, PD_6, 921600);
Ticker sampling_timer;
InterruptIn button(BUTTON1);
Timer debounce;
DigitalOut led(LED1);
PwmOut pwm1(PE_9);
PwmOut pwm2(PE_11);
DigitalOut pwm_enable(PF_15);
InterruptIn signal1(PA_0);
InterruptIn signal2(PB_9);
Timer encoder_timer;
AnalogIn ADCin(PA_3);
AnalogOut DACout(PA_4);

int main()
{
    float num = 0.0923;
    float den = 1;
    user_button_init(); 
    system_input_init(1, 12/3.3);
    system_output_init(1, 3591.84);
    send_data_RT(2);
    sampling_timer_init();
    controller_tf_init(&num, &den, 1, 1, 0);
    menu_init();
    while (true) {
        system_loop();
    }
}