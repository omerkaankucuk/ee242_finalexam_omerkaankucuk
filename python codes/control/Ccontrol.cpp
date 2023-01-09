#include "Ccontrol.h"

void create_tf(tf_struct *tf, float *num, float *den, int nbr_num_coeffs, int nbr_den_coeffs, float ts){
	tf->num = num;
	tf->den = den;
	tf->num_order = nbr_num_coeffs - 1;
	tf->den_order = nbr_den_coeffs - 1;
	tf->sampling_time = ts;

	int print_cnt;
	int num_cnt = 0;
	int den_cnt = 0;
	int dash_cnt = 0;
	int number_int;
	char num_buf[300];
	char den_buf[300];

	for(print_cnt = 0; print_cnt <= tf->num_order; print_cnt++){
		number_int = (int) (tf->num[print_cnt]*1000000);
		if(print_cnt == tf->num_order) {
			if(number_int != 0){
				if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt," + %g", tf->num[print_cnt]);
				else  num_cnt += sprintf(num_buf+num_cnt," - %g", -tf->num[print_cnt]);
			}
		}
		else if(print_cnt == ((tf->num_order) - 1)){
			if(number_int != 0){
				if((abs(number_int)) == 1000000){
					if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt," + z");
					else  num_cnt += sprintf(num_buf+num_cnt," - z");
				}
				else {
					if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt," + %g z", tf->num[print_cnt]);
					else  num_cnt += sprintf(num_buf+num_cnt," - %g z", -tf->num[print_cnt]);
				}
			}
		}
		else {
			if(number_int != 0){
				if(print_cnt == 0){
					if((abs(number_int)) == 1000000){
						if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt,"z^%d", tf->num_order);
						else num_cnt += sprintf(num_buf+num_cnt,"-z^%d", tf->num_order);
					}
					else {
						if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt,"%g z^%d", tf->num[print_cnt], tf->num_order);
						else num_cnt += sprintf(num_buf+num_cnt,"-%g z^%d", -tf->num[print_cnt], tf->num_order);
					}
				}
				else {
					if((abs(number_int)) == 1000000){
						if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt," + z^%d", (tf->num_order - print_cnt));
						else num_cnt += sprintf(num_buf+num_cnt," - z^%d", (tf->num_order - print_cnt));
					}
					else {
						if((tf->num[print_cnt]) >= 0) num_cnt += sprintf(num_buf+num_cnt," + %g z^%d", tf->num[print_cnt], (tf->num_order - print_cnt));
						else num_cnt += sprintf(num_buf+num_cnt," - %g z^%d", -tf->num[print_cnt], (tf->num_order - print_cnt));
					}
				}
			}
		}
	}

	for(print_cnt = 0; print_cnt <= tf->den_order; print_cnt++){
		number_int = (int) (tf->den[print_cnt]*1000000);
		if(print_cnt == tf->den_order) {
			if(number_int != 0){
				if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt," + %g", tf->den[print_cnt]);
				else  den_cnt += sprintf(den_buf+den_cnt," - %g", -tf->den[print_cnt]);
			}
		}
		else if(print_cnt == ((tf->den_order) - 1)){
			if(number_int != 0){
				if((abs(number_int)) == 1000000){
					if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt," + z");
					else  den_cnt += sprintf(den_buf+den_cnt," - z");
				}
				else {
					if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt," + %g z", tf->den[print_cnt]);
					else  den_cnt += sprintf(den_buf+den_cnt," - %g z", -tf->den[print_cnt]);
				}
			}
		}
		else {
			if(number_int != 0){
				if(print_cnt == 0){
					if((abs(number_int)) == 1000000){
						if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt,"z^%d", tf->den_order);
						else den_cnt += sprintf(den_buf+den_cnt,"-z^%d", tf->den_order);
					}
					else {
						if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt,"%g z^%d", tf->den[print_cnt], tf->den_order);
						else den_cnt += sprintf(den_buf+den_cnt,"-%g z^%d", -tf->den[print_cnt], tf->den_order);
					}
				}
				else {
					if((abs(number_int)) == 1000000){
						if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt," + z^%d", (tf->den_order - print_cnt));
						else den_cnt += sprintf(den_buf+den_cnt," - z^%d", (tf->den_order - print_cnt));
					}
					else {
						if((tf->den[print_cnt]) >= 0) den_cnt += sprintf(den_buf+den_cnt," + %g z^%d", tf->den[print_cnt], (tf->den_order - print_cnt));
						else den_cnt += sprintf(den_buf+den_cnt," - %g z^%d", -tf->den[print_cnt], (tf->den_order - print_cnt));
					}
				}
			}
		}
	}

	if(num_cnt > den_cnt) dash_cnt = num_cnt;
	else dash_cnt = den_cnt;
	if(dash_cnt >= 80) dash_cnt = 80;

	_print_stars(80);
	printf("SYSTEM TRANSFER FUNCTION:");
	_print_nextline();
	_print_nextline();
	printf("%s\n", num_buf);
	for(print_cnt = 0; print_cnt < dash_cnt; print_cnt++){
		if(print_cnt == (dash_cnt - 1)) printf("-\n");
		else printf("-");
	}
	printf("%s\n", den_buf);
	printf("\n");
	printf("Sample time: %g\n", tf->sampling_time);
	_print_stars(80);
}

void _print_stars(int star_number){
	int print_cnt;
	_print_nextline();
	_print_nextline();
	for(print_cnt = 0; print_cnt < star_number; print_cnt++){
		printf("*");
	}
	_print_nextline();
	_print_nextline();
}

void _print_nextline(void){
	printf("\n");
}

float unit_pulse_sample(int cnt){
    float sample;
    if(cnt == 0)sample = 1;       
    else sample = 0;
    return sample;    
}

void unit_pulse_signal(int N, float signal[]){
    int n;
    for(n = 0; n < N; n++){
        signal[n] = unit_pulse_sample(n);    
    }    
}

float step_sample(int cnt, float Amp){
    float sample;
    sample = Amp;
    return sample;   
}

void step_signal(int N, float Amp, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = step_sample(n, Amp);    
    } 
}

float ramp_sample(int cnt, float Amp){
    float sample;
    sample = cnt*Amp;
    return sample;
}

void ramp_signal(int N, float Amp, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = ramp_sample(n, Amp);    
    } 
}

float parabolic_sample(int cnt, float Amp){
    float sample;
    sample = cnt*cnt*Amp/2;
    return sample;
}
    
void parabolic_signal(int N, float Amp, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = parabolic_sample(n, Amp);    
    } 
}

float exponential_sample(int cnt, float a){
    float e = 2.718281828459045;
    float sample;
    sample = pow(e, cnt*a);
    return sample;
}
    
void exponential_signal(int N, float a, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = exponential_sample(n, a);    
    } 
}

float sinusoidal_sample(int cnt, float Amp, float Freq, float Phase, float Offset, float Fsample, int select){
    float sample;   
    if(select == 0) sample = Amp*sin((2*PI*Freq*cnt/Fsample)+Phase) + Offset;
    else if(select == 1) sample = Amp*cos((2*PI*Freq*cnt/Fsample)+Phase) + Offset;
    return sample;
}

void sinusoidal_signal(int N, float Amp, float Freq, float Phase, float Offset, float Fsample, int select, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = sinusoidal_sample(n, Amp, Freq, Phase, Offset, Fsample, select);    
    }     
}

float damped_sinusoidal_sample(int cnt, float Amp_exp, float Amp_sin, float Freq, float Phase, float Offset, float Fsample, int select){
    float sample;
    sample = exponential_sample(cnt, Amp_exp)*sinusoidal_sample(cnt, Amp_sin, Freq, Phase, Offset, Fsample, select);
    return sample;  
}
    
void damped_sinusoidal_signal(int N, float Amp_exp, float Amp_sin, float Freq, float Phase, float Offset, float Fsample, int select, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = damped_sinusoidal_sample(n, Amp_exp, Amp_sin, Freq, Phase, Offset, Fsample, select);    
    }  
}

float rectangular_sample(int cnt, float Amp, float Period, float Duty, float Offset, float Fsample){
    float sample;
    int period_cnt, duty_cnt, cnt_mod, level;
    period_cnt = int(Period * Fsample);
    duty_cnt = int(period_cnt * Duty/100);
    cnt_mod = cnt % period_cnt;
    if(cnt_mod < duty_cnt) level = 1;  
    else level = 0;
    sample = level * Amp + Offset;
    return sample;   
}

void rectangular_signal(int N, float Amp, float Period, float Duty, float Offset, float Fsample, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = rectangular_sample(n, Amp, Period, Duty, Offset, Fsample);    
    } 
}

float sum_of_sinusoids_sample(int cnt, int No_Sines, float Amps[], float Freqs[], float Phases[], float Offsets[], float Fsample, int select){
    float sample;
    int n;
    sample = 0;
    if(select == 0){
        for(n = 0; n < No_Sines; n++){
            sample += Amps[n]*sin((2*PI*Freqs[n]*cnt/Fsample)+Phases[n]) + Offsets[n];
        }
    }
    else if(select == 1){
        for(n = 0; n < No_Sines; n++){
            sample += Amps[n]*cos((2*PI*Freqs[n]*cnt/Fsample)+Phases[n]) + Offsets[n];
        }  
    }
    return sample;
}

void sum_of_sinusoids_signal(int N, int No_Sines, float Amps[], float Freqs[], float Phases[], float Offsets[], float Fsample, int select, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = sum_of_sinusoids_sample(n, No_Sines, Amps, Freqs, Phases, Offsets, Fsample, select);    
    }    
}
   
float sweep_sample(int cntex, float Amp, float Freq_start, float Freq_stop, float Offset, float Fsample, float duration){
    float sample;
    int duration_cnt;
    int cnt;
    duration_cnt = int(duration * Fsample);
    cnt = cntex%duration_cnt;
    sample = Amp*(sin((2*PI*Freq_start*cnt/Fsample)+(2*PI*(Freq_stop-Freq_start)*cnt*cnt/(2*(duration_cnt-1)*Fsample))))+Offset;
    return sample;
}
    
void sweep_signal(int N, float Amp, float Freq_start, float Freq_stop, float Offset, float Fsample, float duration, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = sweep_sample(n, Amp, Freq_start, Freq_stop, Offset, Fsample, duration);    
    }    
}

float rand_old;
float random_sample(int cnt, float Amp, int Duration, float Offset){
    float sample_rand;
    int mod_cnt;
    float rand_num;
    mod_cnt = cnt % Duration;
    rand_num = ((float)(rand()%10000))/10000;
    if(mod_cnt == 0){
        sample_rand = (Amp * rand_num) + Offset;
        rand_old = sample_rand;
    } 
    else sample_rand = rand_old;
    return sample_rand;
}
   
void random_signal(int N, float Amp, int Duration, float Offset, float signal[]){
int n;
    for(n = 0; n < N; n++){
        signal[n] = random_sample(n, Amp, Duration, Offset);    
    }  
}

void _iir_filter(float b[], float a[], float x[], float y[], int size1, int size2, int size3){
    int n, m, k;
    float num[size2], den[size3], sum_iir;
    for(n = 0; n < size2; n++){
        num[n] = b[n] / a[0];    
    } 
    for(n = 0; n < size3; n++){
        den[n] = a[n] / a[0];  
    } 
    den[0] = 1;

    for(n = 0; n < size1; n++){
        sum_iir = 0;
        for(k = 0; k < size3; k++){
            if(k < n + 1){
                if(k < (size3 - size2))sum_iir += 0 * x[n - k];  
                else sum_iir += num[k - (size3 - size2)] * x[n - k + (size3 - size2) - 1];        
            }
        }
        for(m = 1; m < size3; m++){
            if(m < n + 1) sum_iir -= den[m] * y[n - m];       
        }
        y[n] = sum_iir;
    }
}


void lsim(tf_struct *tf,int N, float x[], float y[]){
    _iir_filter(tf->num, tf->den, x, y, N, (tf->num_order + 1), tf->den_order + 1);
}

void _create_tf(tf_struct *tf, float *num, float *den, int nbr_num_coeffs, int nbr_den_coeffs, float ts){
	tf->num = num;
	tf->den = den;
	tf->num_order = nbr_num_coeffs - 1;
	tf->den_order = nbr_den_coeffs - 1;
	tf->sampling_time = ts;
}

float system_output_per_sample(tf_struct *tf, float x[], float y[], int n){
    float sum_iir = 0;
    int k, m;
    float result;
    int n_mod;
    for(k = 0; k<(tf->num_order + 1); k++){
        if(n<(tf->num_order + 1)){
            if(n>=k){
                sum_iir += tf->num[k] * x[n - k];
            }
        }
        else{
            n_mod = n % (tf->num_order + 1);
            if(n_mod - k >= 0){
                sum_iir += tf->num[k] * x[n_mod - k];
            }
            else{
                sum_iir += tf->num[k] * x[n_mod - k + (tf->num_order + 1)];
            }            
        }
    }
    for(k = 1; k<(tf->den_order + 1); k++){
        if(n<(tf->den_order + 1)){
            if(n>=k){
                sum_iir -= tf->den[k] * y[n - k];
            }
        }
        else{
            n_mod = n % (tf->den_order + 1);
            if(n_mod - k >= 0){
                sum_iir -= tf->den[k] * y[n_mod - k];
            }
            else{
                sum_iir -= tf->den[k] * y[n_mod - k + (tf->den_order + 1)];
            }              
        }
    }
    result = sum_iir/tf->den[0];
    y[n % (tf->den_order + 1)] = result;
    return result;
}

void add_matrices(float *x1, float *x2, float *result, int size_row, int size_column){
    int i, j;
      for (i = 0; i < size_row; i++) {
            for (j = 0 ; j < size_column; j++) {
                result[i*size_column+j] = x1[i*size_column+j] + x2[i*size_column+j];
            }
      }
}

void sub_matrices(float *x1, float *x2, float *result, int size_row, int size_column){
    int i, j;
      for (i = 0; i < size_row; i++) {
            for (j = 0 ; j < size_column; j++) {
                result[i*size_column+j] = x1[i*size_column+j] - x2[i*size_column+j];
            }
      }
}

void transpose(float *x, float *result, int size_row, int size_column){
    int i, j;
      for (j = 0 ; j < size_column; j++) { 
            for (i = 0; i < size_row; i++){
                result[j*size_row+i] = x[i*size_column+j];
            }
      }
}

void mult_matrix_scalar(float *x, float A, float *result, int size_row, int size_column){
    int i, j;
      for (i = 0; i < size_row; i++) {
            for (j = 0 ; j < size_column; j++) {
                result[i*size_column+j] = A * x[i*size_column+j];
            }
      }
}

void mult_matrix(float *x1, float *x2, float *result, int size_row_x1, int size_column_x1, int size_column_x2){
    int i, j, k;
    float dummy = 0;
      for (i = 0; i < size_row_x1; i++) {
            for (j = 0 ; j < size_column_x2; j++) {
                for (k = 0 ; k < size_column_x1; k++) {
                    dummy = dummy + x1[i*size_column_x1+k] * x2[k*size_column_x2+j];
                }
                result[i*size_column_x2+j] = dummy;
                dummy = 0;
            }
      }
}

void create_ss(ss_struct *ss, float *F, float *G, float *C, float *D, int no_states, float ts){
	ss->F = F;
	ss->G = G;
    ss->C = C;
	ss->D = D;
	ss->no_states = no_states;
	ss->sampling_time = ts;
}

void lsim_ss(ss_struct *ss, float *u, float *x0, float *result, int select, int N){
    float x_old[20];
    float x_new[20];
    float dummy1[20];
    float dummy2[20];
    int i, j, k;
    float dummy = 0;

    for(i = 0; i < ss->no_states; i++) {
        x_old[i] = x0[i];
    } 
    for(i = 0; i < N; i++) {
        mult_matrix(ss->F, x_old, dummy1, ss->no_states, ss->no_states, 1);
        mult_matrix(ss->G, &u[i], dummy2, ss->no_states, 1, 1);
        add_matrices(dummy1, dummy2, x_new, ss->no_states, 1);
        mult_matrix(ss->C, x_old, dummy1, 1, ss->no_states, 1);
        mult_matrix(ss->D, &u[i], dummy2, 1, 1, 1);
        add_matrices(dummy1, dummy2, &result[i], 1, 1);

        if(select != 0){
            result[i] = x_old[select-1];
        }
        for(j = 0; j < ss->no_states; j++) {
            x_old[j] = x_new[j];
        }   
    }
}

void send_data(float signal[],int L){

typedef union _data {
  float f;
  char  s[4];
} myData;

extern UnbufferedSerial device;
myData signalf;
int n;

    for(n = 0; n < L; n++){
    signalf.f =  signal[n];
    device.write(signalf.s,sizeof(signalf.s));
    }  
}
