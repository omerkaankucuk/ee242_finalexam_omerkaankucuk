#include "menu.h"

extern BufferedSerial pc;

char buffer[100];
int buffer_cnt;

menu_struct menu;

void get_string_to_buffer(void){
    buffer_cnt = 0;
    while(1){
        pc.read(buffer + buffer_cnt, 1);
        if((buffer[buffer_cnt] == '\r') || (buffer[buffer_cnt] == '\n')){
            break;
        }
        else if(buffer[buffer_cnt] != 0x7F){
            printf("%c",buffer[buffer_cnt]);
            buffer_cnt++;
        }
        else{
            if(buffer_cnt>0){
                printf("%c",0x08);
                printf("%c",0x20);
                printf("%c",0x08);
                buffer_cnt--;
            }
        }         
    }
}

void menu_set_sampling_frequency(void){
    printf("\nPlease enter the desired sampling frequency in Hz\n"); 
    printf("Note: Samplng frequency must be between 0 and 5000 Hz\n");     
    printf("Note2: Enter 0 for using initial value\n");
    printf("Sampling frequency = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.sampling_frequency); 
    while((menu.sampling_frequency<0)||(menu.sampling_frequency>5000)){
        printf("\nInvalid frequency. Please enter the desired sampling frequency in Hz\n"); 
        printf("Sampling frequency = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.sampling_frequency);
    }
    if(menu.sampling_frequency == 0) menu.sampling_frequency = Sampling_freq_initial;    
    printf("\nSampling frequency is selected as %d Hz\n",(int)menu.sampling_frequency);  
}

void menu_set_PWM_frequency(void){
    printf("\nPlease enter the desired PWM frequency in Hz\n");
    printf("Note: PWM frequency must be greater than 0 Hz\n");     
    printf("Note2: Enter 0 for using initial value\n");
    printf("PWM frequency = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.PWM_frequency);
        while(menu.PWM_frequency<0){
            printf("\nPWM frequency cannot be negative. Please enter the desired PWM frequency in Hz\n"); 
            printf("PWM frequency = ");
            get_string_to_buffer();
            sscanf(buffer, "%g", &menu.PWM_frequency);          
        }
    if(menu.PWM_frequency == 0) menu.PWM_frequency = PWM_freq_initial;
    printf("\nPWM frequency is selected as %d Hz\n",(int)menu.PWM_frequency); 
}

void menu_input_selection(void){
    printf("\nPlease enter the number of desired input signal given below:\n");
    printf("1 = Constant Amplitude\n");
    printf("2 = Unit Pulse Signal\n");
    printf("3 = Step Signal\n");  
    printf("4 = Ramp Signal\n");
    printf("5 = Parabolic Signal\n");
    printf("6 = Exponential Signal\n");  
    printf("7 = Sinusoidal Signal\n");
    printf("8 = Damped Sinusoidal Signal\n");
    printf("9 = Rectangular Signal\n");  
    printf("10 = Sum Of sinusoids Signal\n");
    printf("11 = Sweep Signal\n");
    printf("12 = Random Signal\n");  
    printf("Selection = ");
    get_string_to_buffer();
    sscanf(buffer, "%d", &menu.signal_select);       
    while((menu.signal_select<1)||(menu.signal_select>12)){
            printf("\nInvalid selection. Please enter the number of desired input signal\n"); 
            printf("Selection = ");
            get_string_to_buffer();
            sscanf(buffer, "%d", &menu.signal_select);  
    }
    if(menu.signal_select == 1) printf("\nConstant Amplitude is selected as input\n");
    else if(menu.signal_select == 2) printf("\nUnit Pulse Signal is selected as input\n"); 
    else if(menu.signal_select == 3) printf("\nStep Signal is selected as input\n"); 
    else if(menu.signal_select == 4) printf("\nRamp Signal is selected as input\n"); 
    else if(menu.signal_select == 5) printf("\nParabolic Signal is selected as input\n"); 
    else if(menu.signal_select == 6) printf("\nExponential Signal is selected as input\n"); 
    else if(menu.signal_select == 7) printf("\nSinusoidal Signal is selected as input\n"); 
    else if(menu.signal_select == 8) printf("\nDamped Sinusoidal Signal is selected as input\n"); 
    else if(menu.signal_select == 9) printf("\nRectangular Signal is selected as input\n"); 
    else if(menu.signal_select == 10) printf("\nSum Of sinusoids Signal is selected as input\n"); 
    else if(menu.signal_select == 11) printf("\nSweep Signal is selected as input\n"); 
    else if(menu.signal_select == 12) printf("\nRandom Signal is selected as input\n");  

    if(menu.signal_select == 1)constant_amplitude_input();
    if(menu.signal_select == 3)step_signal_input();
    else if(menu.signal_select == 4)ramp_signal_input();           
    else if(menu.signal_select == 5)parabolic_signal_input();             
    else if(menu.signal_select == 6)exponential_signal_input();          
    else if(menu.signal_select == 7)sinusoidal_signal_input();
    else if(menu.signal_select == 8){
        damped_exponential_signal_input();       
        sinusoidal_signal_input();
    }
    else if(menu.signal_select == 9)rectangular_signal_input();
    else if(menu.signal_select == 10)sum_of_sinusoids_signal_input();
    else if(menu.signal_select == 11)sweep_signal_input();
    else if(menu.signal_select == 12)random_signal_input();        
}

void constant_amplitude_input(void){
    printf("\nPlease enter the desired input amplitude\n");
    printf("Amplitude = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude);
    if(menu.amplitude<0)menu.amplitude = 0;                      
    printf("\nInput amplitude is selected as %d\n",(int)menu.amplitude);  
}

void step_signal_input(void){
    printf("\nPlease enter the desired step amplitude\n");   
    printf("Amplitude = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude); 
    if(menu.amplitude<0)menu.amplitude = 0;                     
    printf("\nStep amplitude is selected as %d\n",(int)menu.amplitude);
}

void ramp_signal_input(void){
    printf("\nPlease enter the desired ramp slope\n");  
    printf("Slope = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude); 
    if(menu.amplitude<0)menu.amplitude = 0;                                 
    printf("\nRamp slope is selected as %d\n",(int)menu.amplitude); 
}

void parabolic_signal_input(void){
    printf("\nPlease enter the desired parabolic slope\n");     
    printf("Slope = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude);
    if(menu.amplitude<0)menu.amplitude = 0;                             
    printf("\nParabolic slope is selected as %d\n",(int)menu.amplitude); 
}

void exponential_signal_input(void){
    printf("\nPlease enter the desired exponential rate\n");
    printf("Rate = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude); 
    if(menu.amplitude<0)menu.amplitude = 0;                              
    printf("\nExponential rate is selected as %d\n",(int)menu.amplitude);   
}

void sinusoidal_signal_input(void){
    printf("\nPlease enter the desired sinusoidal amplitude\n");     
    printf("Amplitude = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude);   
    if(menu.amplitude<0)menu.amplitude = 0;                  
    printf("\nSinusoidal amplitude is selected as %d\n",(int)menu.amplitude); 

    printf("\nPlease enter the desired sinusoidal frequency in Hz\n");
    printf("Note: Sinusoidal frequency must be greater than 0 Hz\n");
    printf("Frequency = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.frequency);  
    while(menu.frequency<=0){
        printf("\nInvalid frequency. Please enter the desired sinusoidal frequency\n"); 
        printf("Frequency = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.frequency); 
    }                     
    printf("\nSinusoidal frequency is selected as %d Hz\n",(int)menu.frequency); 

    printf("\nPlease enter the desired sinusoidal phase in radians\n");   
    printf("Phase = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.phase);                     
    printf("\nSinusoidal phase is selected as %d radians\n",(int)menu.phase);  

    printf("\nPlease enter the desired sinusoidal offset\n");
    printf("Offset = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.offset);
    if(menu.offset<0)menu.offset = 0;                       
    printf("\nSinusoidal offset is selected as %d\n",(int)menu.offset);     

    printf("\nPlease enter the desired sinusoidal type given below:\n");
    printf("0 = Sine wave\n");
    printf("1 = Cosine wave\n");
    printf("Selection = ");
    get_string_to_buffer();
    sscanf(buffer, "%d", &menu.sine_select);          
    while((menu.sine_select<0)||(menu.sine_select>1)){
        printf("\nInvalid selection. Please enter the desired sinusoidal type\n"); 
        printf("Selection = ");
        get_string_to_buffer();
        sscanf(buffer, "%d", &menu.sine_select); 
    }
    if(menu.sine_select == 0)printf("\nSinusoidal type is selected as sine wave\n"); 
    else if(menu.sine_select == 1)printf("\nSinusoidal type is selected as cosine wave\n");   
}

void damped_exponential_signal_input(void){
    printf("\nPlease enter the desired exponential rate\n");
    printf("Note: Exponential rate must be less than 0 in order to obtain damped sinusoidal\n");
    printf("Rate = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude_exp);                               
    while(menu.amplitude_exp>=0){
        printf("\nInvalid exponential rate. Please enter the desired exponential rate\n"); 
        printf("Rate = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.amplitude_exp); 
    }                     
    printf("\nExponential rate is selected as %d\n",(int)menu.amplitude_exp);   
}

void rectangular_signal_input(void){
    printf("\nPlease enter the desired rectangular amplitude\n");    
    printf("Amplitude = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude);   
    if(menu.amplitude<0)menu.amplitude = 0;                     
    printf("\nRectangular amplitude is selected as %d\n",(int)menu.amplitude); 

    printf("\nPlease enter the desired rectangular period in seconds\n");
    printf("Note: Rectangular period must be greater than 0\n");
    printf("Period = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.period);  
    while(menu.period<0){
        printf("\nInvalid period. Please enter the desired rectangular period\n"); 
        printf("Period = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.period); 
    }                     
    printf("\nRectangular period is selected as %d seconds\n",(int)menu.period); 

    printf("\nPlease enter the desired rectangular pulse duration\n");
    printf("Note: Rectangular pulse duration must be selected between 0 and rectangular period\n");
    printf("Pulse duration = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.dutycycle);  
    while((menu.dutycycle<0)||(menu.dutycycle>menu.period)){
        printf("\nInvalid pulse duration. Please enter the desired rectangular pulse duration\n"); 
        printf("Pulse duration = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.dutycycle); 
    }                     
    printf("\nRectangular pulse duration is selected as %d seconds\n",(int)menu.dutycycle);  

    printf("\nPlease enter the desired rectangular offset\n");
    printf("Offset = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.offset);  
    if(menu.offset<0)menu.offset = 0;                  
    printf("\nRectangular offset is selected as %d\n",(int)menu.offset);     
}

void sum_of_sinusoids_signal_input(void){
    int cnt_dummy;
    printf("\nPlease enter the desired number of sinusoids to sum\n");
    printf("Note: Number of sinusoids must be selected between 2 and 10\n");      
    printf("Number = ");
    get_string_to_buffer();
    sscanf(buffer, "%d", &menu.sine_number);  
    while((menu.sine_number<2)||(menu.sine_number>10)){
        printf("\nInvalid number. Please enter the desired number of sinusoids to sum\n"); 
        printf("Number = ");
        get_string_to_buffer();
        sscanf(buffer, "%d", &menu.sine_number); 
    }                     
    printf("\nNumber of sinusoids to sum is selected as %d\n",menu.sine_number); 
    menu.amplitudes = (float*) malloc(menu.sine_number * sizeof(float));
    menu.frequencies = (float*) malloc(menu.sine_number * sizeof(float));
    menu.phases = (float*) malloc(menu.sine_number * sizeof(float));
    menu.offsets = (float*) malloc(menu.sine_number * sizeof(float));
    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        printf("\nPlease enter the desired amplitude for sinusoidal number %d\n", cnt_dummy+1);
        printf("Amplitude %d = ", cnt_dummy+1);
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.amplitudes[cnt_dummy]);  
        if(menu.amplitudes[cnt_dummy]<0)menu.amplitudes[cnt_dummy] = 0;   
    }
    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        if(cnt_dummy == 0){
            printf("\nSinusoidal amplitudes are selected as %d, ",(int)menu.amplitudes[cnt_dummy]);
        }
        else if(cnt_dummy<(menu.sine_number)-1){
            printf("%d, ",(int)menu.amplitudes[cnt_dummy]);
        }
        else{
            printf("and %d Volts\n",(int)menu.amplitudes[cnt_dummy]);
        }
    } 

    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        printf("\nPlease enter the desired frequency for sinusoidal number %d in Hz\n", cnt_dummy+1);
        printf("Note: Sinusoidal frequency must be greater than 0Hz\n");
        printf("Frequency %d = ", cnt_dummy+1);
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.frequencies[cnt_dummy]);  
        while(menu.frequencies[cnt_dummy]<0){
            printf("\nInvalid frequency. Please enter the desired frequency for sinusoidal number %d in Hz\n", cnt_dummy+1); 
            printf("Frequency %d = ", cnt_dummy+1);
            get_string_to_buffer();
            sscanf(buffer, "%g", &menu.frequencies[cnt_dummy]); 
        }       
    }
    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        if(cnt_dummy == 0){
            printf("\nSinusoidal frequencies are selected as %d, ",(int)menu.frequencies[cnt_dummy]);
        }
        else if(cnt_dummy<(menu.sine_number)-1){
            printf("%d, ",(int)menu.frequencies[cnt_dummy]);
        }
        else{
            printf("and %d Hz\n",(int)menu.frequencies[cnt_dummy]);
        }
    } 

    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        printf("\nPlease enter the desired phase for sinusoidal number %d in radians\n", cnt_dummy+1);
        printf("Phase %d = ", cnt_dummy+1);
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.phases[cnt_dummy]);       
    }
    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        if(cnt_dummy == 0){
            printf("\nSinusoidal phases are selected as %d, ",(int)menu.phases[cnt_dummy]);
        }
        else if(cnt_dummy<(menu.sine_number)-1){
            printf("%d, ",(int)menu.phases[cnt_dummy]);
        }
        else{
            printf("and %d radians\n",(int)menu.phases[cnt_dummy]);
        }
    } 

    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        printf("\nPlease enter the desired offset for sinusoidal number %d\n", cnt_dummy+1);
        printf("Offset %d = ", cnt_dummy+1);
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.offsets[cnt_dummy]);  
        if(menu.offsets[cnt_dummy]<0)menu.offsets[cnt_dummy] = 0;   
    }
    for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
        if(cnt_dummy == 0){
            printf("\nSinusoidal offsets are selected as %d, ",(int)menu.offsets[cnt_dummy]);
        }
        else if(cnt_dummy<(menu.sine_number)-1){
            printf("%d, ",(int)menu.offsets[cnt_dummy]);
        }
        else{
            printf("and %d Volts\n",(int)menu.offsets[cnt_dummy]);
        }
    }   
 
    printf("\nPlease enter the desired sinusoidal type given below:\n");
    printf("0 = Sine wave\n");
    printf("1 = Cosine wave\n");
    printf("Selection = ");
    get_string_to_buffer();
    sscanf(buffer, "%d", &menu.sine_select);          
    while((menu.sine_select<0)||(menu.sine_select>1)){
        printf("\nInvalid selection. Please re-enter the desired sinusoidal type\n"); 
        printf("Selection = ");
        get_string_to_buffer();
        sscanf(buffer, "%d", &menu.sine_select); 
    }
    if(menu.sine_select == 0)printf("\nSinusoidal type is selected as sine wave\n"); 
    else if(menu.sine_select == 1)printf("\nSinusoidal type is selected as cosine wave\n");   
}

void sweep_signal_input(void){
    printf("\nPlease enter the desired sweep signal amplitude\n");     
    printf("Amplitude = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude); 
    if(menu.amplitude<0)menu.amplitude = 0;                           
    printf("\nSweep signal amplitude is selected as %d\n",(int)menu.amplitude); 

    printf("\nPlease enter the desired starting frequency for sweep signal in Hz\n");
    printf("Note: Starting frequency for sweep signal must be larger 0 Hz\n");
    printf("Frequency = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.frequency_start);  
    while(menu.frequency_start<0){
        printf("\nInvalid frequency. Please enter the desired starting frequency for sweep signal in Hz\n"); 
        printf("Frequency = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.frequency_start); 
    }                     
    printf("\nStarting frequency for sweep signal is selected as %d Hz\n",(int)menu.frequency_start); 

    printf("\nPlease enter the desired stopping frequency for sweep signal in Hz\n");
    printf("Note: Stopping frequency for sweep signal must be selected between starting frequency and sampling frequency/2 Hz\n");
    printf("Frequency = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.frequency_stop);  
    while((menu.frequency_stop>(menu.sampling_frequency/2))||(menu.frequency_stop<menu.frequency_start)){
        printf("\nInvalid frequency. Please enter the desired stopping frequency for sweep signal in Hz\n"); 
        printf("Frequency = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.frequency_stop); 
    }                     
    printf("\nStopping frequency for sweep signal is selected as %d Hz\n",(int)menu.frequency_stop);  

    printf("\nPlease enter the desired sweep signal offset\n");
    printf("Offset = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.offset);   
    if(menu.offset<0)menu.offset = 0;                  
    printf("\nSweep signal offset is selected as %d\n",(int)menu.offset);    

    printf("\nPlease enter the desired duration in seconds to sweep signal reach its stopping frequency from starting frequency\n");
    printf("Note: Sweep signal duration must be greater than 0\n");
    printf("Duration = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.sweep_duration);  
    while(menu.sweep_duration<0){
        printf("\nInvalid duration. Please enter the desired sweep signal duration in seconds\n"); 
        printf("Duration = ");
        get_string_to_buffer();
        sscanf(buffer, "%g", &menu.sweep_duration); 
    }                     
    printf("\nSweep signal duration is selected as %d seconds\n",(int)menu.sweep_duration);   
}

void random_signal_input(void){
    printf("\nPlease enter the desired maximum random signal amplitude\n");     
    printf("Amplitude = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.amplitude);
    if(menu.amplitude<0)menu.amplitude = 0;                           
    printf("\nMaximum random signal amplitude is selected as %d\n",(int)menu.amplitude); 

    printf("\nPlease enter the desired random signal constant amplitude duration in number of samples\n");
    printf("Note: Random signal constant amplitude duration duration must be greater than 0\n");
    printf("Duration = ");
    get_string_to_buffer();
    sscanf(buffer, "%d", &menu.random_duration);  
    while(menu.random_duration<0){
        printf("\nInvalid duration. Please enter the desired random signal constant amplitude duration in number of samples\n"); 
        printf("Duration = ");
        get_string_to_buffer();
        sscanf(buffer, "%d", &menu.random_duration); 
    }                     
    printf("\nRandom signal constant amplitude duration is selected as %d samples\n",menu.random_duration); 

    printf("\nPlease enter the desired random signal offset\n");
    printf("Offset = ");
    get_string_to_buffer();
    sscanf(buffer, "%g", &menu.offset); 
    if(menu.offset<0)menu.offset = 0;                      
    printf("\nRandom signal offset is selected as %d\n",(int)menu.offset);     
}

void menu_system_ready(void){
    printf("\nSystem is ready. Please push the user button to start system\n"); 
} 

void menu_system_running(void){
    printf("\nSystem is running. Please push the user button to stop the system\n"); 
}

void menu_system_stopped(void){
    printf("\nSystem is stopped. Please push the user button to start the system\n"); 
}

void menu_motor_direction_selection(void){
    int cnt_dummy;
    printf("\nPlease enter the desired motor direction given below:\n");
    printf("0 = Counter-clockwise direction\n");
    printf("1 = Clockwise direction\n");
    printf("Direction = ");
    get_string_to_buffer();
    sscanf(buffer, "%d", &menu.motor_direction);          
    while((menu.motor_direction<0)||(menu.motor_direction>1)){
        printf("\nInvalid direction. Please enter the desired motor direction\n"); 
        printf("Direction = ");
        get_string_to_buffer();
        sscanf(buffer, "%d", &menu.motor_direction); 
    }
    if(menu.motor_direction == 0){                
        menu.amplitude = - menu.amplitude;
        menu.offset = -menu.offset;
        for(cnt_dummy = 0; cnt_dummy<menu.sine_number; cnt_dummy++){
            menu.amplitudes[cnt_dummy] = -menu.amplitudes[cnt_dummy];
            menu.offsets[cnt_dummy] = -menu.offsets[cnt_dummy];
        } 
        printf("\nMotor direction is selected as counter-clockwise\n");  
    }
    else if(menu.motor_direction == 1){
        printf("\nMotor direction is selected as clockwise\n"); 
    } 
}