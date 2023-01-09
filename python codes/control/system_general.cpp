#include "system_general.h"

#define RLS_Disable

extern UnbufferedSerial device;
extern Ticker sampling_timer;
extern InterruptIn button;
extern Timer debounce;
extern DigitalOut led;
#ifndef RLS_Disable
    extern RLS_Controller rls;
#endif
bool sample_flag = false;
bool start_flag = false;
bool stop_flag = false;
bool send_data_flag = false;
int cnt_signal = 0;
float system_input_voltage;
float system_output_voltage = 0;
float signal_input;
int send_select;
int signal_select_buf; 
float sampling_period;

void isr_button(void){
    if (debounce.elapsed_time() > 200ms) {
        led = !led;        
        if(led.read() == 1)start_flag = true;
        else if(led.read() == 0)stop_flag = true;
        debounce.reset();
  }
}

void user_button_init(void){
    debounce.start();
    button.rise(&isr_button);
}

void isr_sampling_timer(void){
    sample_flag = true;
}  

void sampling_timer_init(void){
    menu_set_sampling_frequency();
    sampling_period = 1 / menu.sampling_frequency;
    sample_flag = false;
}

void send_data_RT(int select){
    send_data_flag = true;                                   
    send_select = select;
}

void send_variable(float var){
    myData varf;
    varf.f = var;
    device.write(varf.s,sizeof(varf.s));
}

void select_signal(void){
    if(menu.signal_select == 0){
        signal_input = 0;
        system_output_voltage = 0;
    } 
    if(menu.signal_select == 1){ 
        signal_input = menu.amplitude;      
    }   
    else if(menu.signal_select == 2){  
        signal_input = unit_pulse_sample(cnt_signal);  
    }
    else if(menu.signal_select == 3){
        signal_input = step_sample(cnt_signal, menu.amplitude);       
    }
    else if(menu.signal_select == 4){        
        signal_input = ramp_sample(cnt_signal, menu.amplitude); 
    }
    else if(menu.signal_select == 5){    
        signal_input = parabolic_sample(cnt_signal, menu.amplitude);           
    }
    else if(menu.signal_select == 6){
        signal_input = exponential_sample(cnt_signal, menu.amplitude);         
    }
    else if(menu.signal_select == 7){
        signal_input = sinusoidal_sample(cnt_signal, menu.amplitude, menu.frequency, menu.phase, menu.offset, menu.sampling_frequency, menu.sine_select);       
    }
    else if(menu.signal_select == 8){
        signal_input = damped_sinusoidal_sample(cnt_signal, menu.amplitude_exp, menu.amplitude, menu.frequency, menu.phase, menu.offset, menu.sampling_frequency, menu.sine_select);   
    }
    else if(menu.signal_select == 9){
        signal_input = rectangular_sample(cnt_signal, menu.amplitude, menu.period, menu.dutycycle, menu.offset, menu.sampling_frequency);            
    }
    else if(menu.signal_select == 10){
        signal_input = sum_of_sinusoids_sample(cnt_signal, menu.sine_number, menu.amplitudes, menu.frequencies, menu.phases, menu.offsets, menu.sampling_frequency, menu.sine_select);             
    }
    else if(menu.signal_select == 11){
        signal_input = sweep_sample(cnt_signal, menu.amplitude, menu.frequency_start, menu.frequency_stop, menu.offset, menu.sampling_frequency, menu.sweep_duration);            
    }
    else if(menu.signal_select == 12){
        signal_input = random_sample(cnt_signal, menu.amplitude, menu.random_duration, menu.offset);  
    }
}

void start_function(void){
    menu.signal_select = signal_select_buf;
    cnt_signal = 0;
    menu_system_running();         
    start_flag = false;       
    sampling_timer.attach(&isr_sampling_timer, sampling_period); 
    if(system_output.type == 1){
        system_output.time_diff = 0;
        system_output.read_old = 0;
        system_output.read_new = 0;
    } 
}

void stop_function(void){
    signal_select_buf = menu.signal_select;
    menu.signal_select = 0;
    menu_system_stopped();
    stop_flag = false;
    sampling_timer.detach(); 
    sample_flag = 0;
}

void menu_init(void){
    menu_input_selection();
    if(system_input.type == 1)menu_motor_direction_selection();
    signal_select_buf = menu.signal_select;
    menu_system_ready();  
}

void system_loop(void){
    if(start_flag == true){
        start_function();
    }
    if(stop_flag == true){
        generate_system_input(0);
        stop_function();
    } 
    if(sample_flag == true){
        if(system_input.flag == true){
            select_signal();
            if((system_controller.flag == 1)&&(system_controller.type != 2)){
                system_input_voltage = obtain_controller_output(signal_input, system_output_voltage, cnt_signal);
            }
            else {
                system_input_voltage = signal_input;
            }               
            generate_system_input(system_input_voltage);        
        }                         
        if(system_output.flag == 1)system_output_voltage = obtain_system_output(); 
        if(send_data_flag == true){
            if(send_select == 0){
                send_variable(system_input_voltage);
            }
            else if(send_select == 1){
                send_variable(system_output_voltage);
            }
            else if(send_select == 2){                   
                send_variable(system_input_voltage);
                send_variable(system_output_voltage);
            }
            else if(send_select == 3){                   
                send_variable(system_output_voltage);
                send_variable(system_controller.x1);
            }
        }
        #ifndef RLS_Disable
            if(system_controller.type == 2){
                float RLS_check = rls.obtain_controller_output(system_input_voltage, system_output_voltage, cnt_signal);    
                if(RLS_check > 0){
                    generate_system_input(0);
                    stop_function();
                    led = !led;
                    while(device.getc() != 'r');             
                    for(int i = 0; i < rls._nb+1; i++){
                        send_variable(rls._num[i]);
                    }                   
                    for(int i = 0; i < rls._na+1; i++){
                        send_variable(rls._den[i]);
                    }                                                
                } 
            }
        #endif
        cnt_signal++;  
        sample_flag = false; 
    }
}