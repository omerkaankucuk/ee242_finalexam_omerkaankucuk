import mpcontrol
from math import atan2, factorial, pow, fmod, pi, log10, degrees, sqrt, log, cos
from cmath import exp
from control import tf as tfcontrol
from control import ss as sscontrol
import serial
from ctypes import *


class ss:
    def __init__(self, A, B, C, D, ts=None):
        """
        Create continuous-time or discrete-time state space model.

        Parameters
        ----------
        A : state matrix.
        B : input matrix.
        C : output matrix.
        D : feedthrough matrix.
        ts: sampling time. If it is not defined, continuous-time state space model is created. Otherwise,
            discrete-time state space model is created using this sampling time value.

        Returns
        -------
        Created transfer function.

        Raises
        ------
        ValueError("A matrix must be square!")
            if row and column sizes of A matrix are not eqaul.
        ValueError("A and B matrices must have the same number of rows!")
            if row sizes of A and B matrices are not eqaul.
        ValueError("A and C matrices have the same number of columns!")
            if column sizes of A and C matrices are not eqaul.
         ValueError("C and D matrices must have the same number of rows!")
            if row sizes of C and D matrices are not eqaul.
         ValueError("B and D matrices have the same number of columns!")
            if column sizes of B and D matrices are not eqaul.
        ValueError("Sampling time is not valid!")
            if ts is 0

        Examples
        --------
        >>> ss_continuous = ss([[1.4, 1.9], [1.0, 0.5]], [[1.0], [0.0]], [[0.8, 0.2]], [[0.02]])
        >>> ss_discrete = ss([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)

        """
        if len(A) != len(A[0]):
            raise ValueError("A matrix must be square!")
        if len(A) != len(B):
            raise ValueError("A and B matrices must have the same number of rows!")
        if len(A[0]) != len(C[0]):
            raise ValueError("A and C matrices have the same number of columns!")
        if len(C) != len(D):
            raise ValueError("C and D matrices must have the same number of rows!")
        if len(B[0]) != len(D[0]):
            raise ValueError("B and D matrices have the same number of columns!")

        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.no_states = len(A)
        self.no_inputs = len(D[0])
        self.no_outputs = len(D)

        if ts is None or ts >= 0:
            self.sampling_time = ts
        else:
            raise ValueError("Sampling time is not valid!")

    def is_dtime(self):
        """
        Check if the state space model is discrete-time or not.

        Returns
        -------
        True: if the state space model is discrete-time.
        False: if the state space model is continuous-time.

        Examples
        --------
        >>> ss_discrete = ss([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)
        >>> check = ss_discrete.is_dtime()

        """
        if self.sampling_time is None or self.sampling_time == 0:
            return False
        else:
            return True

    def __str__(self):
        """
        Create a string to print state space model in desired format.

        Returns
        -------
        String representation of state space model.

        Examples
        --------
        >>> ss_discrete = ss([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)
        >>> print(ss_discrete)

        """
        buf_out = ""
        for n in range(len(self.A)):
            if n == 0:
                if (len(self.A) - 1) == 0:
                    buf_out = "A = [" + str(self.A[n]) + "]" + "\n" * 2
                else:
                    buf_out = "A = [" + str(self.A[n]) + "\n"
            elif n == (len(self.A) - 1):
                buf_out += str(self.A[n]) + "]" + "\n" * 2
            else:
                buf_out += str(self.A[n]) + "\n"
        for n in range(len(self.B)):
            if n == 0:
                if (len(self.B) - 1) == 0:
                    buf_out += "B = [" + str(self.B[n]) + "]" + "\n" * 2
                else:
                    buf_out += "B = [" + str(self.B[n]) + "\n"
            elif n == (len(self.B) - 1):
                buf_out += str(self.B[n]) + "]" + "\n" * 2
            else:
                buf_out += str(self.B[n]) + "\n"
        for n in range(len(self.C)):
            if n == 0:
                if (len(self.C) - 1) == 0:
                    buf_out += "C = [" + str(self.C[n]) + "]" + "\n" * 2
                else:
                    buf_out += "C = [" + str(self.C[n]) + "\n"
            elif n == (len(self.C) - 1):
                buf_out += str(self.C[n]) + "]" + "\n" * 2
            else:
                buf_out += str(self.C[n]) + "\n"
        for n in range(len(self.D)):
            if n == 0:
                if (len(self.D) - 1) == 0:
                    buf_out += "D = [" + str(self.D[n]) + "]" + "\n" * 2
                else:
                    buf_out += "D = [" + str(self.D[n]) + "\n"
            elif n == (len(self.D) - 1):
                buf_out += str(self.D[n]) + "]" + "\n" * 2
            else:
                buf_out += str(self.D[n]) + "\n"

        buf_out += str(self.sampling_time) + "\n"
        return buf_out


class tf:
    def __init__(self, num, den, ts=None):
        """
        Create continuous-time or discrete-time transfer function.

        Parameters
        ----------
        num : numerator coefficients of transfer function.
        den : denominator coefficients of transfer function.
        ts: sampling time. If it is not defined, continuous-time transfer function is created. Otherwise,
            discrete-time transfer function is created using this sampling time value.

        Returns
        -------
        Created transfer function.

        Raises
        ------
        ValueError("Denominator must be nonzero vector!")
            if denominator is zero
        ValueError("Sampling time is not valid!")
            if ts is 0

        Examples
        --------
        >>> tf_continuous = tf([10], [1, 10])
        >>> tf_discrete = tf([3], [1, -0.2], 0.01)

        """
        cnt = 0
        for n in range(len(den)):
            if den[n] != 0:
                cnt += 1
            if cnt == 0:
                raise ValueError("Denominator must be nonzero vector!")

        cnt = 0
        for n in range(len(num)):
            if num[n] != 0:
                cnt += 1
            if cnt == 0:
                self.num = [0]
                self.den = [1]
                self.num_order = 0
                self.den_order = 0
            else:
                self.num = num
                self.den = den
                self.num_order = len(num) - 1
                self.den_order = len(den) - 1

        if ts is None or ts >= 0:
            self.sampling_time = ts
        else:
            raise ValueError("Sampling time is not valid!")

    def is_dtime(self):
        """
        Check if the transfer function is discrete-time or not.

        Returns
        -------
        True: if the transfer function is discrete-time.
        False: if the transfer function is continuous-time.

        Examples
        --------
        >>> tf_discrete = tf([3], [1, -0.2], 0.01)
        >>> check = tf_discrete.is_dtime()

        """
        if self.sampling_time is None or self.sampling_time == 0:
            return False
        else:
            return True

    def pole(self):
        """
        Find the poles of transfer function and print them.

        Returns
        -------
        Poles of transfer function.

        Examples
        --------
        >>> tf_discrete = tf([3], [1, -0.2], 0.01)
        >>> poles = tf_discrete.pole()

        """
        pole_arr = get_roots(self.den)
        print('\n' + "SYSTEM POLES:" + '\n')
        if len(pole_arr) == 0:
            print("System has no pole" + '\n')
        else:
            for n in range(len(pole_arr)):
                if pole_arr[n].real == 0:
                    if pole_arr[n].imag == 0:
                        print("0\n")
                    elif pole_arr[n].imag > 0:
                        buf = "%gi\n" % pole_arr[n].imag
                        print(buf)
                    else:
                        buf = "-%gi\n" % (-pole_arr[n].imag)
                        print(buf)
                elif pole_arr[n].real > 0:
                    if pole_arr[n].imag == 0:
                        buf = "%g\n" % pole_arr[n].real
                        print(buf)
                    elif pole_arr[n].imag > 0:
                        buf = "%g + %gi\n" % (pole_arr[n].real, pole_arr[n].imag)
                        print(buf)
                    else:
                        buf = "%g - %gi\n" % (pole_arr[n].real, -pole_arr[n].imag)
                        print(buf)
                else:
                    if pole_arr[n].imag == 0:
                        buf = "-%g\n" % (-pole_arr[n].real)
                        print(buf)
                    elif pole_arr[n].imag > 0:
                        buf = "-%g + %gi\n" % (-pole_arr[n].real, pole_arr[n].imag)
                        print(buf)
                    else:
                        buf = "-%g - %gi\n" % (-pole_arr[n].real, -pole_arr[n].imag)
                        print(buf)
        return pole_arr

    def zero(self):
        """
        Find the zeros of transfer function and print them.

        Returns
        -------
        Zeros of transfer function.

        Examples
        --------
        >>> tf_discrete = tf([3], [1, -0.2], 0.01)
        >>> zeros = tf_discrete.zero()

        """
        zero_arr = get_roots(self.num)
        print('\n' + "SYSTEM ZEROS:" + '\n')
        if len(zero_arr) == 0:
            print("System has no zero" + '\n')
        else:
            for n in range(len(zero_arr)):
                if zero_arr[n].real == 0:
                    if zero_arr[n].imag == 0:
                        print("0\n")
                    elif zero_arr[n].imag > 0:
                        buf = "%gi\n" % zero_arr[n].imag
                        print(buf)
                    else:
                        buf = "-%gi\n" % (-zero_arr[n].imag)
                        print(buf)
                elif zero_arr[n].real > 0:
                    if zero_arr[n].imag == 0:
                        buf = "%g\n" % zero_arr[n].real
                        print(buf)
                    elif zero_arr[n].imag > 0:
                        buf = "%g + %gi\n" % (zero_arr[n].real, zero_arr[n].imag)
                        print(buf)
                    else:
                        buf = "%g - %gi\n" % (zero_arr[n].real, -zero_arr[n].imag)
                        print(buf)
                else:
                    if zero_arr[n].imag == 0:
                        buf = "-%g\n" % (-zero_arr[n].real)
                        print(buf)
                    elif zero_arr[n].imag > 0:
                        buf = "-%g + %gi\n" % (-zero_arr[n].real, zero_arr[n].imag)
                        print(buf)
                    else:
                        buf = "-%g - %gi\n" % (-zero_arr[n].real, -zero_arr[n].imag)
                        print(buf)
        return zero_arr

    def is_stable(self):
        """
        Check if the discrete-time transfer function is stable or not and print the result.

        Returns
        -------
        True: if the discrete-time transfer function is stable.
        False: if the discrete-time transfer function is not stable.

        Raises
        ------
        NotImplementedError
            if transfer function is not discrete-time.

        Examples
        --------
        >>> tf_discrete = tf([3], [1, 0.2], 0.01)
        >>> check = tf_discrete.is_stable()

        """
        if not self.is_dtime():
            raise NotImplementedError("System must be discrete-time!")
        else:
            poles = get_roots(self.den)
            stability = is_system_stable(poles)
            print('\n' + "SYSTEM STABILITY:" + '\n')
            if not stability:
                print("System is not stable!\n")
            else:
                print("System is stable.\n")
            return stability

    def damp(self):
        """
        Find poles, pole magnitudes, pole damping ratios, pole natural frequencies and pole time constants
            of discrete-time transfer function and print the results.

        Returns
        -------
        poles: poles of the discrete-time transfer function.
        pole_mags: pole magnitudes of the discrete-time transfer function.
        pole_damps: pole damping ratios of the discrete-time transfer function.
        pole_naturalfreqs: pole natural frequencies of the discrete-time transfer function in radians/second.
        pole_timeconstants: pole time constants of the discrete-time transfer function in seconds.

        Raises
        ------
        NotImplementedError
            if transfer function is not discrete-time.

        Examples
        --------
        >>> tf_discrete = tf([3], [1, -0.2], 0.01)
        >>> poles, pole_mags, pole_damps, pole_naturalfreqs, pole_timeconstants = tf_discrete.damp()

        """
        if not self.is_dtime():
            raise NotImplementedError("System must be discrete-time!")
        else:
            pole_arr = get_roots(self.den)
            pole_mags = find_pole_magnitudes(pole_arr)
            pole_damps = find_pole_dampings(pole_arr)
            pole_naturalfreqs = find_pole_naturalfreqs(pole_arr, self.sampling_time)
            pole_timeconstants = find_pole_timeconstants(pole_arr, self.sampling_time)

            print('\n' + "SYSTEM POLES:" + '\n')
            if len(pole_arr) == 0:
                print("System has no pole" + '\n')
            else:
                for n in range(len(pole_arr)):
                    if pole_arr[n].real == 0:
                        if pole_arr[n].imag == 0:
                            print("0\n")
                        elif pole_arr[n].imag > 0:
                            buf = "%gi\n" % pole_arr[n].imag
                            print(buf)
                        else:
                            buf = "-%gi\n" % (-pole_arr[n].imag)
                            print(buf)
                    elif pole_arr[n].real > 0:
                        if pole_arr[n].imag == 0:
                            buf = "%g\n" % pole_arr[n].real
                            print(buf)
                        elif pole_arr[n].imag > 0:
                            buf = "%g + %gi\n" % (pole_arr[n].real, pole_arr[n].imag)
                            print(buf)
                        else:
                            buf = "%g - %gi\n" % (pole_arr[n].real, -pole_arr[n].imag)
                            print(buf)
                    else:
                        if pole_arr[n].imag == 0:
                            buf = "-%g\n" % (-pole_arr[n].real)
                            print(buf)
                        elif pole_arr[n].imag > 0:
                            buf = "-%g + %gi\n" % (-pole_arr[n].real, pole_arr[n].imag)
                            print(buf)
                        else:
                            buf = "-%g - %gi\n" % (-pole_arr[n].real, -pole_arr[n].imag)
                            print(buf)
            print('\n' + "SYSTEM POLE MAGNITUDES:" + '\n')
            for n in range(len(pole_mags)):
                buf = "%g\n" % (pole_mags[n])
                print(buf)
            print('\n' + "SYSTEM POLE DAMPINGS:" + '\n')
            for n in range(len(pole_damps)):
                buf = "%g\n" % (pole_damps[n])
                print(buf)
            print('\n' + "SYSTEM POLE NATURAL FREQUENCIES:" + '\n')
            for n in range(len(pole_naturalfreqs)):
                buf = "%g\n" % (pole_naturalfreqs[n])
                print(buf)
            print('\n' + "SYSTEM POLE TIME CONSTANTS:" + '\n')
            for n in range(len(pole_timeconstants)):
                buf = "%g\n" % (pole_timeconstants[n])
                print(buf)
            return pole_arr, pole_mags, pole_damps, pole_naturalfreqs, pole_timeconstants

    def stepinfo(self, stepgain=1, ST=0.02, RT=None):
        """
        Find the step response parameters and print the results.

        Parameters
        ----------
        stepgain: amplitude of the step input. Default is 1 for unit-step input.
        ST : threshold value for defining settling time. Default is 0.02.
        RT : threshold values for defining rise time. Default is [0.1 0.9].

        Returns
        -------
        rise_time: time value for step response to rise from (100*RT[0])% to (100*RT[1])% of the steady-state value.
        settling_time: time value for the error between step response and steady-state value to fall within
            (100*ST)% of steady-state value.
        settling_min: minimum value of step response after rise-time is reached.
        settling_max: maximum value of step response after rise-time is reached.
        overshoot: percentage overshoot.
        undershoot: percentage undershoot.
        peak: peak absolute value of step response.
        peak_time: time value at which the peak occurs.

        Raises
        ------
        NotImplementedError
            if transfer function is not discrete-time.

        Examples
        --------
        >>> tf_discrete = tf([0.4], [1, -0.9], 0.001)
        >>> rt, st, settling_min, settling_max, overshoot, undershoot, peak, peak_time = tf_discrete.stepinfo()

        """
        if not self.is_dtime():
            raise NotImplementedError("System must be discrete-time!")
        else:
            if RT is None:
                RT = [0.1, 0.9]
            settling_time = 0
            settling_min = 0

            poles = get_roots(self.den)
            stability = is_system_stable(poles)

            if stability == 0:
                peak = float('INF')
                peak_time = float('INF')
                overshoot = float('NAN')
                undershoot = float('NAN')
                rise_time = float('NAN')
                settling_time = float('NAN')
                settling_min = float('NAN')
                settling_max = float('NAN')
            else:
                ss_time = find_ss_time(self)
                x = [stepgain for _i in range(ss_time)]
                y = mpcontrol.iir_filter(self.num, self.den, x, ss_time)
                y_max = max(y)
                y_max_index = y.index(y_max)
                peak = abs(y_max)
                peak_time = (float(y_max_index)) * self.sampling_time
                settling_max = y_max
                y_final = stepgain*poly_value(self.num, 1) / poly_value(self.den, 1)
                if y_max > y_final:
                    overshoot = 100 * (y_max - y_final) / y_final
                else:
                    overshoot = 0.0
                y_min = min(y)
                if y_min < 0:
                    undershoot = 100 * abs(y_min) / y_final
                else:
                    undershoot = 0.0

                t_high = 0.0
                t_low = 0.0
                for n in range(y_max_index):
                    if (y[n] - y[0]) >= ((y_final - y[0]) * RT[0]):
                        t_low = float(n) * self.sampling_time
                        break
                for n in range(y_max_index):
                    if (y[n] - y[0]) >= ((y_final - y[0]) * RT[1]):
                        t_high = float(n) * self.sampling_time
                        settling_min = y[n]
                        break
                rise_time = t_high - t_low
                if rise_time == 0:
                    settling_min = y_max
                for n in range(y_max_index, ss_time, 1):
                    if y[n] < settling_min:
                        settling_min = y[n]

                st_flag = 0
                for n in range(ss_time):
                    if bool(st_flag) == 0:
                        if (y[n] < (y_final * (1 + ST))) and (y[n] > (y_final * (1 - ST))):
                            st_flag = 1
                            settling_time = float(n) * self.sampling_time
                    else:
                        if (y[n] >= (y_final * (1 + ST))) or (y[n] <= (y_final * (1 - ST))):
                            st_flag = 0

            print('\n' + "STEP INFO:" + '\n')
            buf = "Rise Time (T_rt): %g sec" % rise_time
            print(buf)
            buf = "Settling Time (T_st): %g sec" % settling_time
            print(buf)
            buf = "Peak Time (T_pt): %g sec" % peak_time
            print(buf)
            buf = "Overshoot (M_p): %g%%" % overshoot
            print(buf)
            buf = "Undershoot: %g%%" % undershoot
            print(buf)
            buf = "Settling Min.: %g" % settling_min
            print(buf)
            buf = "Settling Max.: %g" % settling_max
            print(buf)
            buf = "Peak: %g\n" % peak
            print(buf)

            return rise_time, settling_time, settling_min, settling_max, overshoot, undershoot, peak, peak_time

    def ss_error(self, input_type='Step'):
        """
        Find the steady state error of discrete-time transfer function for selected input type and print the result.

        Parameters
        ----------
        input_type: type of the input signal. It can be 'Step', 'Ramp' or 'Parabolic'. Default is 'Step'.

        Returns
        -------
        Calculated steady-state error value.

        Raises
        ------
        NotImplementedError("System must be discrete-time!")
            if transfer function is not discrete-time
        NotImplementedError("Invalid Input!")
            if input type is not selected 'Step', 'Ramp' or 'Parabolic'

        Examples
        --------
        >>> tf_discrete = tf([3], [1, -0.2], 0.01)
        >>> error = tf_discrete.ss_error('Step')

        """
        if not self.is_dtime():
            raise NotImplementedError("System must be discrete-time!")
        else:
            poles = get_roots(self.den)
            stability = is_system_stable(poles)
            if stability == 0:
                ss_err = float('INF')
            else:
                if self.num_order > self.den_order:
                    tf_num = []
                    for n in range(self.num_order + 1):
                        tf_num.append(0)
                    for n in range(self.num_order + 1):
                        if (self.num_order - self.den_order) > n:
                            tf_num[n] = - self.num[n]
                        else:
                            tf_num[n] = self.den[n - (self.num_order - self.den_order)] - self.num[n]
                elif self.num_order == self.den_order:
                    tf_num = []
                    for n in range(self.num_order + 1):
                        tf_num.append(0)
                    for n in range(self.num_order + 1):
                        tf_num[n] = self.den[n] - self.num[n]
                else:
                    tf_num = []
                    for n in range(self.den_order + 1):
                        tf_num.append(0)
                    for n in range(self.den_order + 1):
                        if (self.den_order - self.num_order) > n:
                            tf_num[n] = self.den[n]
                        else:
                            tf_num[n] = self.den[n] - self.num[n - (self.den_order - self.num_order)]
                if input_type == 'Step':
                    ss_err = poly_value(tf_num, 1) / poly_value(self.den, 1)
                elif input_type == 'Ramp':
                    ss_flag = 0
                    zeros = get_roots(self.num)
                    for n in range(len(zeros)):
                        if zeros[n].imag == 0.0 and zeros[n].real == 1.0:
                            ss_flag = 1
                            break
                    if ss_flag == 0:
                        ss_err = float('INF')
                    else:
                        tf_num_dummy = []
                        for n in range(len(tf_num) - 1):
                            tf_num_dummy.append(0)
                        for n in range(len(tf_num) - 1):
                            tf_num_dummy[n] = self.sampling_time * array_sum(tf_num, n)
                        ss_err = poly_value(tf_num_dummy, 1) / poly_value(self.den, 1)
                elif input_type == 'Parabolic':
                    zeros = get_roots(self.num)
                    ss_flag = 0
                    ss_flag_cnt = 0
                    for n in range(len(zeros)):
                        if zeros[n].imag == 0.0 and zeros[n].real == 1.0:
                            ss_flag_cnt += 1
                            if ss_flag_cnt == 2:
                                ss_flag = 1
                                break
                    if ss_flag == 0:
                        ss_err = float('INF')
                    else:
                        tf_num_dummy = []
                        for n in range(len(tf_num) - 1):
                            tf_num_dummy.append(0)
                        tf_num_dummy2 = []
                        for n in range(len(tf_num) - 2):
                            tf_num_dummy2.append(0)
                        for n in range(len(tf_num) - 1):
                            tf_num_dummy[n] = self.sampling_time * array_sum(tf_num, n)
                        for n in range(len(tf_num) - 2):
                            tf_num_dummy2[n] = self.sampling_time * self.sampling_time * array_sum(tf_num_dummy, n)
                        ss_err = poly_value(tf_num_dummy2, 1) / poly_value(self.den, 1)
                else:
                    raise NotImplementedError("Invalid Method!")
            print('\n' + "STEADY STATE ERROR:" + '\n')
            if input_type == 'Step':
                buf = "Steady State Error for Step Input: %g\n" % ss_err
                print(buf)
            elif input_type == 'Ramp':
                buf = "Steady State Error for Ramp Input: %g\n" % ss_err
                print(buf)
            elif input_type == 'Parabolic':
                buf = "Steady State Error for Parabolic Input: %g\n" % ss_err
                print(buf)
            return ss_err

    def __str__(self):
        """
        Create a string to print transfer function in desired format.

        Returns
        -------
        String representation of transfer function.

        Examples
        --------
        >>> tf_discrete = tf([3], [1, 0.2], 0.01)
        >>> print(tf_discrete)

        """
        buf_num = ""
        buf_den = ""

        if not self.is_dtime():
            buf_out = "CONTINUOUS-TIME SYSTEM TRANSFER FUNCTION:" + "\n" * 2
            x = 's'
        else:
            buf_out = "DISCRETE-TIME SYSTEM TRANSFER FUNCTION:" + "\n" * 2
            x = 'z'
        print('\n')

        if self.num_order == 0 and self.num[0] == 0:
            tf_access = False
            buf_num += '0'
        else:
            tf_access = True

            if self.num_order == 0:
                if self.num[0] > 0:
                    buf_num += "%g" % (self.num[0])
                else:
                    buf_num += "-%g" % (self.num[0])
            elif self.num_order == 1:
                if self.num[0] > 0:
                    if self.num[0] == 1:
                        buf_num += "%s " % x
                    else:
                        buf_num += "%g %s " % (self.num[0], x)
                else:
                    if self.num[0] == -1:
                        buf_num += "-%s " % x
                    else:
                        buf_num += "-%g %s " % (-self.num[0], x)
                if self.num[1] > 0:
                    buf_num += "+ %g" % (self.num[1])
                elif self.num[1] < 0:
                    buf_num += "- %g" % (-self.num[1])
            else:
                for n in range(self.num_order + 1):
                    if self.num[n] != 0:
                        if n == 0:
                            if self.num[n] > 0:
                                if self.num[n] == 1:
                                    buf_num += "%s^%d " % (x, self.num_order - n)
                                else:
                                    buf_num += "%g %s^%d " % (self.num[n], x, self.num_order - n)
                            else:
                                if self.num[n] == -1:
                                    buf_num += "-%s^%d " % (x, self.num_order - n)
                                else:
                                    buf_num += "-%g %s^%d " % (-self.num[n], x, self.num_order - n)
                        elif n < self.num_order - 1:
                            if self.num[n] > 0:
                                if self.num[n] == 1:
                                    buf_num += "+ %s^%d " % (x, self.num_order - n)
                                else:
                                    buf_num += "+ %g %s^%d " % (self.num[n], x, self.num_order - n)
                            else:
                                if self.num[n] == -1:
                                    buf_num += "- %s^%d " % (x, self.num_order - n)
                                else:
                                    buf_num += "- %g %s^%d " % (-self.num[n], x, self.num_order - n)
                        elif n == self.num_order - 1:
                            if self.num[n] > 0:
                                if self.num[n] == 1:
                                    buf_num += "+ %s " % x
                                else:
                                    buf_num += "+ %g %s " % (self.num[n], x)
                            else:
                                if self.num[n] == -1:
                                    buf_num += "- %s " % x
                                else:
                                    buf_num += "- %g %s " % (-self.num[n], x)
                        else:
                            if self.num[n] > 0:
                                buf_num += "+ %g" % (self.num[n])
                            else:
                                buf_num += "- %g" % (-self.num[n])

            if self.den_order == 0:
                if self.den[0] > 0:
                    buf_den += "%g" % (self.den[0])
                else:
                    buf_den += "-%g" % (self.den[0])
            elif self.den_order == 1:
                if self.den[0] > 0:
                    if self.den[0] == 1:
                        buf_den += "%s " % x
                    else:
                        buf_den += "%g %s " % (self.den[0], x)
                else:
                    if self.den[0] == -1:
                        buf_den += "-%s " % x
                    else:
                        buf_den += "-%g %s " % (-self.den[0], x)
                if self.den[1] > 0:
                    buf_den += "+ %g" % (self.den[1])
                elif self.den[1] < 0:
                    buf_den += "- %g" % (-self.den[1])
            else:
                for n in range(self.den_order + 1):
                    if self.den[n] != 0:
                        if n == 0:
                            if self.den[n] > 0:
                                if self.den[n] == 1:
                                    buf_den += "%s^%d " % (x, self.den_order - n)
                                else:
                                    buf_den += "%g %s^%d " % (self.den[n], x, self.den_order - n)
                            else:
                                if self.den[n] == -1:
                                    buf_den += "-%s^%d " % (x, self.den_order - n)
                                else:
                                    buf_den += "-%g %s^%d " % (-self.den[n], x, self.den_order - n)
                        elif n < self.den_order - 1:
                            if self.den[n] > 0:
                                if self.den[n] == 1:
                                    buf_den += "+ %s^%d " % (x, self.den_order - n)
                                else:
                                    buf_den += "+ %g %s^%d " % (self.den[n], x, self.den_order - n)
                            else:
                                if self.den[n] == -1:
                                    buf_den += "- %s^%d " % (x, self.den_order - n)
                                else:
                                    buf_den += "- %g %s^%d " % (-self.den[n], x, self.den_order - n)
                        elif n == self.den_order - 1:
                            if self.den[n] > 0:
                                if self.den[n] == 1:
                                    buf_den += "+ %s " % x
                                else:
                                    buf_den += "+ %g %s " % (self.den[n], x)
                            else:
                                if self.den[n] == -1:
                                    buf_den += "- %s " % x
                                else:
                                    buf_den += "- %g %s " % (-self.den[n], x)
                        else:
                            if self.den[n] > 0:
                                buf_den += "+ %g" % (self.den[n])
                            else:
                                buf_den += "- %g" % (-self.den[n])

        dash_cnt = max(len(buf_num), len(buf_den))
        buf_dashes = '-' * dash_cnt

        if len(buf_num) < dash_cnt:
            buf_num = (' ' * int(round((dash_cnt - len(buf_num)) / 2)) + buf_num)

        if len(buf_den) < dash_cnt:
            buf_den = (' ' * int(round((dash_cnt - len(buf_den)) / 2)) + buf_den)

        buf_out += buf_num + "\n"
        if tf_access:
            buf_out += buf_dashes + "\n"
            buf_out += buf_den + "\n" * 2
        if self.is_dtime():
            buf_out += "Sampling time: " + str(self.sampling_time) + "\n" * 2
        return buf_out


def c2d(tfc, ts, method='bilinear'):
    """
    Convert continuos-time transfer function to discrete-time transfer function using given sampling time
        and conversion method.

    Parameters
    ----------
    tfc : continuous-time transfer function to be converted.
    ts : sampling time for the conversion.
    method: conversion method. 'bilinear' method is used when it is not specified.
        Methods other than 'bilinear' transformation haven't been implemented yet.

    Returns
    -------
    Converted discrete-time transfer function

    Raises
    ------
    ValueError
        if tfc is not continuous-time
    NotImplementedError
        if method is not selected as 'bilinear'

    Examples
    --------
    >>> tf_continuous = tf([10], [1, 10])
    >>> tf_discrete = c2d(tf_continuous, 0.01)

    """
    if tfc.is_dtime():
        raise ValueError("System must be continuous-time!")
    else:
        if method == 'bilinear':
            if tfc.num_order > tfc.den_order:
                tfd_num = [0.0 for _i in range(tfc.num_order + 1)]
                tfd_den = [0.0 for _i in range(tfc.num_order + 1)]

                for n in range(tfc.num_order + 1):
                    for i in range(tfc.num_order + 1):
                        acc = 0
                        for k in range(max(n - tfc.num_order + i, 0), min(i, n) + 1, 1):
                            if fmod(k, 2) == 0:
                                acc += (factorial(i) * factorial(tfc.num_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.num_order - i - n + k))
                            else:
                                acc -= (factorial(i) * factorial(tfc.num_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.num_order - i - n + k))
                        tfd_num[n] += (tfc.num[tfc.num_order - i] * pow((2.0 / ts), i) * acc)

                for n in range(tfc.num_order + 1):
                    for i in range(tfc.den_order + 1):
                        acc = 0
                        for k in range(max(n - tfc.num_order + i, 0), min(i, n) + 1, 1):
                            if fmod(k, 2) == 0:
                                acc += (factorial(i) * factorial(tfc.num_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.num_order - i - n + k))
                            else:
                                acc -= (factorial(i) * factorial(tfc.num_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.num_order - i - n + k))
                        tfd_den[n] += (tfc.den[tfc.den_order - i] * pow((2.0 / ts), i) * acc)

            elif tfc.num_order == tfc.den_order:
                tfd_num = [0.0 for _i in range(tfc.num_order + 1)]
                tfd_den = [0.0 for _i in range(tfc.num_order + 1)]

                for n in range(tfc.num_order + 1):
                    for i in range(tfc.num_order + 1):
                        acc = 0
                        for k in range(max(n - tfc.num_order + i, 0), min(i, n) + 1, 1):
                            if fmod(k, 2) == 0:
                                acc += (factorial(i) * factorial(tfc.num_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.num_order - i - n + k))
                            else:
                                acc -= (factorial(i) * factorial(tfc.num_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.num_order - i - n + k))
                        tfd_num[n] += (tfc.num[tfc.num_order - i] * pow((2.0 / ts), i) * acc)

                for n in range(tfc.den_order + 1):
                    for i in range(tfc.den_order + 1):
                        acc = 0
                        for k in range(max(n - tfc.den_order + i, 0), min(i, n) + 1, 1):
                            if fmod(k, 2) == 0:
                                acc += (factorial(i) * factorial(tfc.den_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.den_order - i - n + k))
                            else:
                                acc -= (factorial(i) * factorial(tfc.den_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.den_order - i - n + k))
                        tfd_den[n] += (tfc.den[tfc.den_order - i] * pow((2.0 / ts), i) * acc)

            else:
                tfd_num = [0.0 for _i in range(tfc.den_order + 1)]
                tfd_den = [0.0 for _i in range(tfc.den_order + 1)]

                for n in range(tfc.den_order + 1):
                    for i in range(tfc.num_order + 1):
                        acc = 0
                        for k in range(max(n - tfc.den_order + i, 0), min(i, n) + 1, 1):
                            if fmod(k, 2) == 0:
                                acc += (factorial(i) * factorial(tfc.den_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.den_order - i - n + k))
                            else:
                                acc -= (factorial(i) * factorial(tfc.den_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.den_order - i - n + k))
                        tfd_num[n] += (tfc.num[tfc.num_order - i] * pow((2.0 / ts), i) * acc)

                for n in range(tfc.den_order + 1):
                    for i in range(tfc.den_order + 1):
                        acc = 0
                        for k in range(max(n - tfc.den_order + i, 0), min(i, n) + 1, 1):
                            if fmod(k, 2) == 0:
                                acc += (factorial(i) * factorial(tfc.den_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.den_order - i - n + k))
                            else:
                                acc -= (factorial(i) * factorial(tfc.den_order - i)) / (
                                            factorial(k) * factorial(i - k) * factorial(n - k) * factorial(
                                                tfc.den_order - i - n + k))
                        tfd_den[n] += (tfc.den[tfc.den_order - i] * pow((2.0 / ts), i) * acc)

            for n in range(1, len(tfd_den), 1):
                tfd_den[n] = tfd_den[n] / tfd_den[0]
            for n in range(len(tfd_num)):
                tfd_num[n] = tfd_num[n] / tfd_den[0]
            tfd_den[0] = 1

            tf_new = tf(tfd_num, tfd_den, ts)
            return tf_new
        else:
            raise NotImplementedError("Invalid Method!")


def step(tfd, stepgain=1, nsample=None):
    """
    Simulate the output of a discrete-time transfer function for step input.

    Parameters
    ----------
    tfd: discrete-time transfer function.
    stepgain: amplitude gain for step input. 1 is used when it is not specified.
    nsample: number of samples. It is calculated automatically when it is not specified.

    Returns
    -------
    Obtained output when step input is applied to given discrete-time transfer function.

    Raises
    ------
    NotImplementedError
        if tf is not discrete-time transfer function
    ValueError
        if nsample is 0

    Examples
    --------
    >>> tf_discrete = tf([3], [1, -0.9], 0.01)
    >>> y1 = step(tf_discrete)
    >>> y2 = step(tf_discrete, stepgain=10, nsample=500)

    """
    if not tfd.is_dtime():
        raise NotImplementedError("System must be discrete-time!")
    else:
        if nsample == 0:
            raise ValueError("Sample size must be nonzero value!")
        else:
            if nsample is None:
                poles = get_roots(tfd.den)
                stability = is_system_stable(poles)
                if stability:
                    length = find_ss_time(tfd)
                else:
                    length = 1000
            else:
                length = nsample
            x = [stepgain for _i in range(length)]
            y = mpcontrol.iir_filter(tfd.num, tfd.den, x, length)
            return y


def series(tf1, tf2):
    """
    Return the series connection of two discrete-time transfer functions, tf1 and tf2.

    Parameters
    ----------
    tf1 : discrete-time transfer function 1.
    tf2 : discrete-time transfer function 2.

    Returns
    -------
    Created discrete-time transfer function as tf1*tf2.

    Raises
    ------
    NotImplementedError
        if both tf1 and tf2 are not discrete-time transfer functions
    TypeError
        if tf1 or tf2 is not tf class
    ValueError
        if sampling times of tf1 and tf2 does not match

    Examples
    --------
    >>> tf_discrete1 = tf([3], [1, 0.2], 0.01)
    >>> tf_discrete2 = tf([1, 2], [1, -2, 0.4], 0.01)
    >>> tf_discrete3 = series(tf_discrete1, tf_discrete2)

    """
    if tf1.is_dtime() is False or tf2.is_dtime() is False:
        raise NotImplementedError("Systems must be discrete-time!")
    else:
        if (type(tf1) != tf) or (type(tf2) != tf):
            raise TypeError("Inputs must be tf classes!")
        else:
            if tf1.sampling_time != tf2.sampling_time:
                raise ValueError("Sampling times must agree!")
            else:
                tf_num = [0.0 for _i in range(tf1.num_order + tf2.num_order + 1)]
                tf_den = [0.0 for _i in range(tf1.den_order + tf2.den_order + 1)]

                for n in range(tf1.num_order + 1):
                    for k in range(tf2.num_order + 1):
                        tf_num[n + k] += tf1.num[n] * tf2.num[k]

                for n in range(tf1.den_order + 1):
                    for k in range(tf2.den_order + 1):
                        tf_den[n + k] += tf1.den[n] * tf2.den[k]

                for n in range(1, len(tf_den), 1):
                    tf_den[n] = tf_den[n] / tf_den[0]
                for n in range(len(tf_num)):
                    tf_num[n] = tf_num[n] / tf_den[0]
                tf_den[0] = 1

                tf_new = tf(tf_num, tf_den, tf1.sampling_time)
                return tf_new


def parallel(tf1, tf2):
    """
    Return the parallel connection of two discrete-time transfer functions, tf1 and tf2.

    Parameters
    ----------
    tf1 : discrete-time transfer function 1.
    tf2 : discrete-time transfer function 2.

    Returns
    -------
    Created discrete-time transfer function as tf1+tf2.

    Raises
    ------
    NotImplementedError
        if both tf1 and tf2 are not discrete-time transfer functions
    TypeError
        if tf1 or tf2 is not tf class
    ValueError
        if sampling times of tf1 and tf2 does not match

    Examples
    --------
    >>> tf_discrete1 = tf([3], [1, 0.2], 0.01)
    >>> tf_discrete2 = tf([1, 2], [1, -2, 0.4], 0.01)
    >>> tf_discrete3 = parallel(tf_discrete1, tf_discrete2)

    """
    if tf1.is_dtime() is False or tf2.is_dtime() is False:
        raise NotImplementedError("Systems must be discrete-time!")
    else:
        if (type(tf1) != tf) or (type(tf2) != tf):
            raise TypeError("Inputs must be tf classes!")
        else:
            if tf1.sampling_time != tf2.sampling_time:
                raise ValueError("Sampling times must agree!")
            else:
                tf_num1 = [0.0 for _i in range(tf1.num_order + tf2.den_order + 1)]
                tf_num2 = [0.0 for _i in range(tf2.num_order + tf1.den_order + 1)]
                if (tf1.num_order + tf2.den_order) > (tf2.num_order + tf1.den_order):
                    tf_num = [0.0 for _i in range(tf1.num_order + tf2.den_order + 1)]
                    for n in range(tf1.num_order + 1):
                        for k in range(tf2.den_order + 1):
                            tf_num1[n + k] += tf1.num[n] * tf2.den[k]
                    for n in range(tf2.num_order + 1):
                        for k in range(tf1.den_order + 1):
                            tf_num2[n + k] += tf2.num[n] * tf1.den[k]

                    for n in range(tf1.num_order + tf2.den_order + 1):
                        if n < ((tf1.num_order + tf2.den_order) - (tf2.num_order + tf1.den_order)):
                            tf_num[n] = tf_num1[n]
                        else:
                            tf_num[n] = tf_num1[n] + tf_num2[
                                n - ((tf1.num_order + tf2.den_order) - (tf2.num_order + tf1.den_order))]

                else:
                    tf_num = [0.0 for _i in range(tf2.num_order + tf1.den_order + 1)]
                    for n in range(tf1.num_order + 1):
                        for k in range(tf2.den_order + 1):
                            tf_num1[n + k] += tf1.num[n] * tf2.den[k]
                    for n in range(tf2.num_order + 1):
                        for k in range(tf1.den_order + 1):
                            tf_num2[n + k] += tf2.num[n] * tf1.den[k]

                    for n in range(tf2.num_order + tf1.den_order + 1):
                        if n < ((tf2.num_order + tf1.den_order) - (tf1.num_order + tf2.den_order)):
                            tf_num[n] = tf_num2[n]
                        else:
                            tf_num[n] = tf_num2[n] + tf_num1[
                                n - ((tf2.num_order + tf1.den_order) - (tf1.num_order + tf2.den_order))]

                tf_den = [0.0 for _i in range(tf1.den_order + tf2.den_order + 1)]
                for n in range(tf1.den_order + 1):
                    for k in range(tf2.den_order + 1):
                        tf_den[n + k] += tf1.den[n] * tf2.den[k]

                for n in range(1, len(tf_den), 1):
                    tf_den[n] = tf_den[n] / tf_den[0]
                for n in range(len(tf_num)):
                    tf_num[n] = tf_num[n] / tf_den[0]
                tf_den[0] = 1

                tf_new = tf(tf_num, tf_den, tf1.sampling_time)
                return tf_new


def feedback(tf1, tf2, sign=-1):
    """
    Return the feedback connection of two discrete-time transfer functions, tf1 and tf2.

    Parameters
    ----------
    tf1 : discrete-time transfer function 1.
    tf2 : discrete-time transfer function 2.
    sign : sign of feedback connection. If it is -1, negative feedback is used. If it is 1, positive feedback
        is used. Negative feedback is used when it is not specified.

    Returns
    -------
    Created discrete-time transfer function as tf1/(1\u00B1tf1*tf2)

    Raises
    ------
    NotImplementedError
        if both tf1 and tf2 are not discrete-time transfer functions
    TypeError
        if tf1 or tf2 is not tf class
    ValueError
        if sampling times of tf1 and tf2 does not match

    Examples
    --------
    >>> tf_discrete1 = tf([3], [1, 0.2], 0.01)
    >>> tf_discrete2 = tf([1, 2], [1, -2, 0.4], 0.01)
    >>> tf_discrete3 = feedback(tf_discrete1, tf_discrete2)     # negative feedback
    >>> tf_discrete4 = feedback(tf_discrete1, tf_discrete2, 1)  # positive feedback

    """
    if tf1.is_dtime() is False or tf2.is_dtime() is False:
        raise NotImplementedError("Systems must be discrete-time!")
    else:
        if (type(tf1) != tf) or (type(tf2) != tf):
            raise TypeError("Inputs must be tf classes!")
        else:
            if tf1.sampling_time != tf2.sampling_time:
                raise ValueError("Sampling times must agree!")
            else:
                tf_num = [0.0 for _i in range(tf1.num_order + tf2.den_order + 1)]
                for n in range(tf1.num_order + 1):
                    for k in range(tf2.den_order + 1):
                        tf_num[n + k] += tf1.num[n] * tf2.den[k]
                tf_den1 = [0.0 for _i in range(tf1.num_order + tf2.num_order + 1)]
                tf_den2 = [0.0 for _i in range(tf1.den_order + tf2.den_order + 1)]
                if (tf1.num_order + tf2.num_order) > (tf1.den_order + tf2.den_order):
                    tf_den = [0.0 for _i in range(tf1.num_order + tf2.num_order + 1)]
                    for n in range(tf1.num_order + 1):
                        for k in range(tf2.num_order + 1):
                            tf_den1[n + k] += tf1.num[n] * tf2.num[k]
                    for n in range(tf1.den_order + 1):
                        for k in range(tf2.den_order + 1):
                            tf_den2[n + k] += tf1.den[n] * tf2.den[k]
                    for n in range(tf1.num_order + tf2.num_order + 1):
                        if n < ((tf1.num_order + tf2.num_order) - (tf1.den_order + tf2.den_order)):
                            tf_den[n] = tf_den1[n]
                        else:
                            if sign == -1:
                                tf_den[n] = tf_den1[n] + tf_den2[
                                    n - ((tf1.num_order + tf2.num_order) - (tf1.den_order + tf2.den_order))]
                            elif sign == 1:
                                tf_den[n] = tf_den1[n] - tf_den2[
                                    n - ((tf1.num_order + tf2.num_order) - (tf1.den_order + tf2.den_order))]
                else:
                    tf_den = [0.0 for _i in range(tf1.den_order + tf2.den_order + 1)]
                    for n in range(tf1.num_order + 1):
                        for k in range(tf2.num_order + 1):
                            tf_den1[n + k] += tf1.num[n] * tf2.num[k]
                    for n in range(tf1.den_order + 1):
                        for k in range(tf2.den_order + 1):
                            tf_den2[n + k] += tf1.den[n] * tf2.den[k]
                    for n in range(tf1.den_order + tf2.den_order + 1):
                        if n < ((tf1.den_order + tf2.den_order) - (tf1.num_order + tf2.num_order)):
                            tf_den[n] = tf_den2[n]
                        else:
                            if sign == -1:
                                tf_den[n] = tf_den2[n] + tf_den1[
                                    n - ((tf1.den_order + tf2.den_order) - (tf1.num_order + tf2.num_order))]
                            elif sign == 1:
                                tf_den[n] = tf_den2[n] - tf_den1[
                                    n - ((tf1.den_order + tf2.den_order) - (tf1.num_order + tf2.num_order))]

                for n in range(1, len(tf_den), 1):
                    tf_den[n] = tf_den[n] / tf_den[0]
                for n in range(len(tf_num)):
                    tf_num[n] = tf_num[n] / tf_den[0]
                tf_den[0] = 1

                tf_new = tf(tf_num, tf_den, tf1.sampling_time)
                return tf_new


def bode(tfd, nsample=None, dB=True, Hz=False, deg=True):
    """
    Find the frequency response of discrete-time transfer function for desired number of frequency points.

    Parameters
    ----------
    tfd: discrete-time transfer function.
    nsample: desired number of frequency points between 0 and 2*pi*fsample.
    dB: calculate magnitude values in dB if it is true. Default is true.
    Hz: calculate frequency values in Hz if it is true, in rad/s if it is false. Default is false.
    deg: calculate phase values in degrees if it is true, in rad if it is false. Default is true.

    Returns
    -------
    mag : magnitude values of frequency response at omega values.
    phase : phase values of frequency response at omega values.
    omega : frequency points of frequency response.

    Raises
    ------
    NotImplementedError
        if tfd is not discrete-time transfer function
    ValueError
        if nsample selected smaller than 2

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> mag1, phase1, omega1 = bode(tf_discrete, nsample=300)
    >>> mag2, phase2, omega2 = bode(tf_discrete, db=False, deg=False)

    """
    if not tfd.is_dtime():
        raise NotImplementedError("System must be discrete-time!")
    else:
        eps = 2.220446049250313e-13
        upper = log10(pi / tfd.sampling_time)
        lower = log10(1)
        if nsample is None:
            length = 200
        elif nsample < 2:
            raise ValueError("Number of frequency points must be larger than 1!")
        else:
            length = nsample
        omega = [(lower + ((i * (upper - lower))/(length-1))) for i in range(length)]
        ejw = [0.0 for _i in range(length)]
        out = [0.0 for _i in range(length)]
        mag = [0.0 for _i in range(length)]
        phase = [0.0 for _i in range(length)]
        for n in range(length):
            omega[n] = pow(10, omega[n])
            ejw[n] = exp(1.j * omega[n] * tfd.sampling_time)
            out[n] = poly_value(tfd.num, ejw[n]) / poly_value(tfd.den, ejw[n])
            mag[n] = abs(out[n])
            if mag[n] == 0:
                if dB:
                    mag[n] = 20 * log10(eps)
            else:
                if dB:
                    mag[n] = 20 * log10(mag[n])
            if(abs(out[n].imag) < eps) and (abs(out[n].real) < eps):
                phase[n] = -pi
            else:
                phase[n] = atan2(out[n].imag, out[n].real)
            if deg:
                phase[n] = degrees(phase[n])
                if phase[n] > 0:
                    phase[n] -= 360
            else:
                if phase[n] > 0:
                    phase[n] -= 2*pi
            if Hz:
                omega[n] = omega[n] / (2*pi)
        return mag, phase, omega


def margin(tfd):
    """
    Find gain margin, phase margin and crossover frequencies of discrete-time transfer function and print them.

    Parameters
    ----------
    tfd: discrete-time transfer function.

    Returns
    -------
    GM : gain margin of discrete-time transfer function.
    PM : phase margin of discrete-time transfer function.
    Wcg : gain crossover frequency of discrete-time transfer function.
    Wcp : phase crossover frequency of discrete-time transfer function.

    Raises
    ------
    NotImplementedError
        if tfd is not discrete-time transfer function

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> GM, PM, Wcg, Wcp = margin(tf_discrete)

    """
    if not tfd.is_dtime():
        raise NotImplementedError("System must be discrete-time!")
    else:
        [mag, angle, freq] = bode(tfd)
        length = len(freq)
        Wcp = 0
        PM = 0
        Wcg = 0
        GM = 0
        eps = 2.220446049250313e-16
        for n in range(1, length, 1):
            if mag[0] <= 0:
                if mag[n - 1] > 0 >= mag[n]:
                    Wcp = (mag[n-1]*log10(freq[n]) - mag[n]*log10(freq[n-1]))/(mag[n-1] - mag[n])
                    Wcp = pow(10, Wcp)
                    ejw = exp(1.j * Wcp * tfd.sampling_time)
                    out = poly_value(tfd.num, ejw) / poly_value(tfd.den, ejw)
                    magx = abs(out)
                    if magx > 0:
                        anglex = atan2(out.imag, out.real)
                    else:
                        anglex = atan2(out.imag, out.real - eps)
                    PM = degrees(anglex) + 180
                    break
                else:
                    Wcp = float('NAN')
                    PM = float('INF')
            else:
                if mag[n] <= 0:
                    Wcp = (mag[n-1]*log10(freq[n]) - mag[n]*log10(freq[n-1]))/(mag[n-1] - mag[n])
                    Wcp = pow(10, Wcp)
                    ejw = exp(1.j * Wcp * tfd.sampling_time)
                    out = poly_value(tfd.num, ejw) / poly_value(tfd.den, ejw)
                    magx = abs(out)
                    if magx > 0:
                        anglex = atan2(out.imag, out.real)
                    else:
                        anglex = atan2(out.imag, out.real - eps)
                    PM = degrees(anglex) + 180
                    break
                else:
                    Wcp = float('NAN')
                    PM = float('INF')
        for n in range(1, length, 1):
            if angle[n] < -180:
                Wcg = ((log10(freq[n-1])*(-180-angle[n]))+(log10(freq[n])*(angle[n-1]+180)))/(angle[n-1]-angle[n])
                Wcg = pow(10, Wcg)
                ejw = exp(1.j * Wcg * tfd.sampling_time)
                out = poly_value(tfd.num, ejw) / poly_value(tfd.den, ejw)
                magx = abs(out)
                if magx == 0:
                    magx = 20 * log10(eps)
                else:
                    magx = 20 * log10(magx)
                GM = -magx
                break
            elif angle[n] == -180:
                Wcg = freq[n]
                ejw = exp(1.j * Wcg * tfd.sampling_time)
                out = poly_value(tfd.num, ejw) / poly_value(tfd.den, ejw)
                magx = abs(out)
                if magx == 0:
                    magx = 20 * log10(eps)
                else:
                    magx = 20 * log10(magx)
                GM = -magx
                break
            else:
                Wcg = float('NAN')
                GM = float('INF')

        print('\n' + "BODE RESPONSE:" + '\n')
        buf = "Gain Margin: %g dB" % GM
        print(buf)
        buf = "Phase Margin: %g deg" % PM
        print(buf)
        buf = "Gain Margin Frequency: %g rad/s" % Wcg
        print(buf)
        buf = "Phase Margin Frequency: %g rad/s" % Wcp
        print(buf)
        return GM, PM, Wcg, Wcp


def pid_tune(tfc, method='ZieglerNichols', pid='PID'):
    """
    Find PID coefficients for given continuous-time transfer function using selected method.

    Parameters
    ----------
    tfc : continuous-time transfer function.
    method: selection for PID tuning method. It can be selected as 'ZieglerNichols',
        'CohenCoon', or 'CHR'. Here, 'CHR' represents Chien-Hrones-Reswich method for 20%
        overshoot with distrubance rejection. Default is 'ZieglerNichols'.
    pid: selection for controller type. It can be selected as 'P', 'PI', 'PID'.
        Default is 'PID'.

    Returns
    -------
    PID coefficients according to pid selection.

    Raises
    ------
    ValueError
        if tfc is not continuous-time
    NotImplementedError
        if pid is not selected as 'P', 'PI' or 'PID'
    NotImplementedError
        if tfc is unstable system

    Examples
    --------
    >>> tf_c = tf([10], [1, 10])
    >>> Kp1, Ki1, Kd1 = pid_tune(tf_c) # calculates PID coefficients using Ziegler-Nichols
    >>> Kp2, Ki2, Kd2 = pid_tune(tf_c, method='CohenCoon', pid='PI')  # calculates PI coefficients using Cohen-Coon

    """
    if tfc.is_dtime():
        raise ValueError("System must be continuous-time!")
    else:
        poles = get_roots(tfc.den)
        pole_mags = find_pole_magnitudes(poles)
        dummy_mags = [0.0 for _i in range(len(poles))]
        for n in range(len(poles)):
            dummy_mags[n] = pole_mags[n]
        for n in range(len(poles)):
            for k in range(len(poles)):
                if dummy_mags[k] < dummy_mags[n]:
                    temp = dummy_mags[n]
                    dummy_mags[n] = dummy_mags[k]
                    dummy_mags[k] = temp
        ts = 1/(100 * (dummy_mags[0]))
        tf_d = c2d(tfc, ts)
        poles = get_roots(tf_d.den)
        stability = is_system_stable(poles)
        if stability is False:
            raise NotImplementedError("System must be stable for PID tuning!")
        else:
            y = step(tf_d)
            t = [0.0 for _i in range(len(y))]
            for n in range(len(t)):
                t[n] = n * tf_d.sampling_time
            y_final = poly_value(tf_d.num, 1) / poly_value(tf_d.den, 1)
            t632 = get_time_near(t, y, 0.632 * y_final)
            t283 = get_time_near(t, y, 0.283 * y_final)
            tau = 1.5 * (t632 - t283)
            L = 1.5 * (t283 - (t632 / 3))
            a = y_final * L / tau
            R = L / tau
            Kp = 0
            Ki = 0
            Kd = 0
            if method == 'ZieglerNichols':
                if pid == 'PID':
                    Kp = 1.2 / a
                    Ti = 2 * L
                    Td = L / 2
                    Ki = Kp / Ti
                    Kd = Kp * Td
                elif pid == 'PI':
                    Kp = 0.9 / a
                    Ti = 3 * L
                    Ki = Kp / Ti
                    Kd = 0
                elif pid == 'P':
                    Kp = 1 / a
                    Ki = 0
                    Kd = 0
                else:
                    raise NotImplementedError("PID type must me P, PI, or PID")
            elif method == 'CohenCoon':
                if pid == 'PID':
                    Kp = (1 / a) * ((4 / 3) + (R / 4))
                    Ti = L * ((32 + 6 * R) / (13 + 8 * R))
                    Td = L * 4 / (11 + 2 * R)
                    Ki = Kp / Ti
                    Kd = Kp * Td
                elif pid == 'PI':
                    Kp = (1 / a) * ((9 / 10) + (R / 12))
                    Ti = L * ((30 + 3 * R) / (9 + 20 * R))
                    Ki = Kp / Ti
                    Kd = 0
                elif pid == 'P':
                    Kp = (1 / a) * (1 + (R / 3))
                    Ki = 0
                    Kd = 0
                else:
                    raise NotImplementedError("PID type must me P, PI, or PID")
            elif method == 'CHR':
                if pid == 'PID':
                    Kp = 1.2 / a
                    Ti = 2 * L
                    Td = 0.42 * L
                    Ki = Kp / Ti
                    Kd = Kp * Td
                elif pid == 'PI':
                    Kp = 0.7 / a
                    Ti = 2.3 * L
                    Ki = Kp / Ti
                    Kd = 0
                elif pid == 'P':
                    Kp = 0.7 / a
                    Ki = 0
                    Kd = 0
                else:
                    raise NotImplementedError("PID type must me P, PI, or PID")
            return Kp, Ki, Kd



def mptf_to_tf(tf_original):
    """
    Convert micropython control library transfer function to original control library transfer function.

    Parameters
    ----------
    tf_original: micropython control library transfer function.

    Returns
    -------
    Original control library tf transfer function.

    Examples
    --------
    >>> tf_mpcontrol = tf([3], [1, 0.2], 0.01)
    >>> tf_control = mptf_to_tf(tf_mpcontrol)

    """
    if not tf_original.is_dtime():
        tf_copy = tfcontrol(tf_original.num, tf_original.den)
    else:
        tf_copy = tfcontrol(tf_original.num, tf_original.den, tf_original.sampling_time)
    return tf_copy


def tf_to_mptf(tf_original):
    """
    Convert original control library transfer function to micropython control library transfer function.

    Parameters
    ----------
    tf_original: original control library transfer function.

    Returns
    -------
    micropython control library transfer function.

    Examples
    --------
    >>> tf_control = control.tf([3], [1, 0.2], 0.01)
    >>> tf_mpcontrol = mptf_to_tf(tf_control)

    """
    if not tf_original.isdtime():
        tf_copy = tf(tf_original.num[0][0], tf_original.den[0][0])
    else:
        tf_copy = tf(tf_original.num[0][0], tf_original.den[0][0], tf_original.dt)
    return tf_copy


def mpss_to_ss(ss_original):
    """
    Convert micropython control library state space model to original control library state space model.

    Parameters
    ----------
    ss_original: micropython control library state space model.

    Returns
    -------
    Original control library state space model.

    Examples
    --------
    >>> ss_mpcontrol = ss([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)
    >>> ss_control = mpss_to_ss(ss_mpcontrol)

    """
    if not ss_original.is_dtime():
        ss_copy = sscontrol(ss_original.A, ss_original.B, ss_original.C, ss_original.D)
    else:
        ss_copy = tfcontrol(ss_original.A, ss_original.B, ss_original.C, ss_original.D, ss_original.sampling_time)
    return ss_copy


def ss_to_mpss(ss_original):
    """
    Convert original control library state space model to micropython control library state space model.

    Parameters
    ----------
    ss_original: original control library state space model.

    Returns
    -------
    micropython control library state space model.

    Examples
    --------
    >>> ss_control = control.ss([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)
    >>> ss_mpcontrol = mpss_to_ss(ss_control)

    """
    if not ss_original.isdtime():
        ss_copy = ss(ss_original.A.tolist(), ss_original.B.tolist(), ss_original.C.tolist(), ss_original.D.tolist())
    else:
        ss_copy = ss(ss_original.A.tolist(), ss_original.B.tolist(), ss_original.C.tolist(), ss_original.D.tolist(), ss_original.dt)
    return ss_copy


class RLS:
    def __init__(self, na, nb, Ts, gain, cnt_max):
        """
        Create RLS model for sample by sample simulation.

        Parameters
        ----------
        na: order of the numerator polynomial.
        nb: order of the denumerator polynomial.
        Ts: sampling time.
        gain: multiplier for the P matrix just for the first step in iteration.
        cnt_max: Maximum number of simulation step.

        Returns
        -------
        created RLS model.

        Examples
        --------
        >>> RLS_controller = RLS(2, 2, 0.0005, 1000000, 2000)

        """
        self.P_old = mpcontrol.mult_matrix_scalar(create_identity_matrix(na + nb + 1), gain)
        self.P = mpcontrol.create_empty_matrix(na + nb + 1, na + nb + 1)
        for i in range(na + nb + 1):
            for j in range(na + nb + 1):
                self.P[i][j] = 50000
        self.theta_hat = mpcontrol.create_empty_matrix(na + nb + 1, 1)
        self.theta_hat[0][0] = -0.6
        self.theta_hat[1][0] = 0.8
        self.theta_hat[2][0] = 0.8
        self.theta_hat_old = mpcontrol.create_empty_matrix(na + nb + 1, 1)
        self.P_trans = mpcontrol.create_empty_matrix(na + nb + 1, na + nb + 1)
        self.K = mpcontrol.create_empty_matrix(na + nb + 1, 1)
        self.phiT = mpcontrol.create_empty_matrix(1, na + nb + 1)
        self.epsilon = 0
        self.x = [0 for _i in range(nb+1)]
        self.y = [0 for _i in range(na+1)]
        self.Ts = Ts
        self.N = max(na+1, nb+1)
        self.na = na
        self.nb = nb
        self.cnt_max = cnt_max
        self.num = [0 for _i in range(nb+1)]
        self.den = [1 for _i in range(na+1)]

    def update(self, xnew, ynew, cnt):
        """
        Update the RLS parameters for each step.

        Parameters
        ----------
        xnew: system input.
        ynew: system output.
        cnt: simulastion step number.

        Returns
        -------
        RLS transfer function.

        Examples
        --------
        >>> RLS_controller = RLS(2, 2, 0.0005, 1000000, 2000)
        >>> RLS_tf = RLS_controller.update(1, 2.5, 300)

        """
        self.x[cnt % (self.nb+1)] = xnew
        self.y[cnt % (self.na+1)] = ynew
        if cnt >= self.N - 1:
            for j in range(self.na):
                if cnt - j < 0:
                    self.phiT[0][j] = 0
                else:
                    if cnt % (self.na + 1) - j - 1 < 0:
                        self.phiT[0][j] = -self.y[cnt % (self.na + 1) - j - 1 + (self.na + 1)]
                    else:
                        self.phiT[0][j] = -self.y[cnt % (self.na + 1) - j - 1]
            for j in range(self.nb + 1):
                if cnt - j < 0:
                    self.phiT[0][j + self.na] = 0
                else:
                    if cnt % (self.nb+1) - j < 0:
                        self.phiT[0][j + self.na] = self.x[cnt % (self.nb+1) - j + (self.nb+1)]
                    else:
                        self.phiT[0][j + self.na] = self.x[cnt % (self.nb+1) - j]
            self.K = mpcontrol.mult_matrices(mpcontrol.mult_matrices(self.P_old, transpose(self.phiT)), inv(
                mpcontrol.add_matrices(mpcontrol.mult_matrices(mpcontrol.mult_matrices(self.phiT, self.P_old), transpose(self.phiT)), [[1]])))
            self.epsilon = ynew - mpcontrol.mult_matrices(self.phiT, self.theta_hat_old)[0][0]
            self.theta_hat = mpcontrol.add_matrices(self.theta_hat_old, mpcontrol.mult_matrix_scalar(self.K, self.epsilon))
            self.P = mpcontrol.mult_matrices(mpcontrol.add_matrices(create_identity_matrix(self.na + self.nb + 1),
                                                mpcontrol.mult_matrix_scalar(mpcontrol.mult_matrices(self.K, self.phiT), -1)), self.P_old)
            self.P_trans = transpose(self.P)
            self.P = mpcontrol.mult_matrix_scalar(mpcontrol.add_matrices(self.P, self.P_trans), 0.5)
            self.P_old = self.P
            self.theta_hat_old = self.theta_hat
        if cnt < self.cnt_max - 1:
            return -1
        else:
            self.num = transpose(self.theta_hat)[0][self.na:(self.na + self.nb + 1)]
            self.den[1:(self.na + 1)] = transpose(self.theta_hat)[0][0:self.na]
            Gz = tf(self.num, self.den, self.Ts)
            return Gz


class RLS_Sim:
    def __init__(self, x, y, na, nb, Ts, gain):
        """
        Create RLS model for array simulation.

        Parameters
        ----------
        x: input array of the system.
        y: output array of the system.
        na: order of the numerator polynomial.
        nb: order of the denumerator polynomial.
        Ts: sampling time.
        gain: multiplier for the P matrix just for the first step in iteration.

        Returns
        -------
        created RLS model.

        Examples
        --------
        >>> RLS_controller = RLS_Sim(x, y, 2, 2, 0.0005, 1000000)

        """
        self.x = x
        self.y = y
        self.na = na
        self.nb = nb
        self.Ts = Ts
        self.gain = gain
        self.size = len(x)
        self.N = max(na+1, nb+1)
        self.epslon = [0 for _i in range(self.size)]
        self.num = mpcontrol.create_empty_matrix(nb + 1, self.size)
        self.den = mpcontrol.create_empty_matrix(na + 1, self.size)

    def Sim(self):
        """
        Simulate the RLS system.

        Returns
        -------
        RLS transfer function.

        Examples
        --------
        >>> RLS_controller = RLS_Sim(x, y, 2, 2, 0.0005, 1000000)
        >>> RLS_tf = RLS_controller.Sim()

        """
        P_old = mpcontrol.mult_matrix_scalar(create_identity_matrix(self.na + self.nb + 1), self.gain)
        phiT = mpcontrol.create_empty_matrix(1, self.na + self.nb + 1)
        theta_hat_old = mpcontrol.create_empty_matrix(self.na + self.nb + 1, 1)
        num_dummy = mpcontrol.create_empty_matrix(self.size, self.nb + 1)
        den_dummy = mpcontrol.create_empty_matrix(self.size, self.na + 1)
        for i in range(self.N - 1, self.size):
            for j in range(self.na):
                if i-j < 0:
                    phiT[0][j] = 0
                else:
                    phiT[0][j] = -self.y[i-j-1]
            for j in range(self.nb + 1):
                if i-j < 0:
                    phiT[0][j+self.na] = 0
                else:
                    phiT[0][j+self.na] = self.x[i-j]
            K = mpcontrol.mult_matrices(mpcontrol.mult_matrices(P_old, transpose(phiT)), inv(mpcontrol.add_matrices(mpcontrol.mult_matrices(mpcontrol.mult_matrices(phiT, P_old), transpose(phiT)), [[1]])))
            self.epslon[i] = self.y[i] - mpcontrol.mult_matrices(phiT, theta_hat_old)[0][0]
            theta_hat = mpcontrol.add_matrices(theta_hat_old, mpcontrol.mult_matrix_scalar(K, self.epslon[i]))
            P = mpcontrol.mult_matrices(mpcontrol.add_matrices(create_identity_matrix(self.na+self.nb+1), mpcontrol.mult_matrix_scalar(mpcontrol.mult_matrices(K, phiT), -1)), P_old)
            P_trans = transpose(P)
            P = mpcontrol.mult_matrix_scalar(mpcontrol.add_matrices(P, P_trans), 0.5)
            P_old = P
            theta_hat_old = theta_hat
            num_dummy[i][:] = transpose(theta_hat)[0][self.na:(self.na+self.nb+1)]
            den_dummy[i][0] = 1
            den_dummy[i][1:(self.na+1)] = transpose(theta_hat)[0][0:self.na]
        self.num = transpose(num_dummy)
        self.den = transpose(den_dummy)
        Gz = tf(num_dummy[self.size-1][:], den_dummy[self.size-1][:], self.Ts)
        return Gz


class MRAC_Sim:
    def __init__(self, x, y_m, tf_system, Ts, gamma):
        self.x = x
        self.y_m = y_m
        self.tf_system = tf_system
        self.Ts = Ts
        self.gamma = gamma
        self.size = len(x)
        self.epsilon = [0 for _i in range(self.size)]
        self.e_m = [0 for _i in range(self.size)]
        self.theta = [0 for _i in range(self.size)]
        self.u = [0 for _i in range(self.size)]
        self.y = [0 for _i in range(self.size)]
        self.tf_adjust = tf([-gamma*Ts, 0], [1, -1], Ts)

    def Sim(self):
        for i in range(self.size):
            if i == 0:
                self.epsilon[i] = -self.y_m[i]
            else:
                self.epsilon[i] = self.y[i - 1] - self.y_m[i]
            self.e_m[i] = self.epsilon[i]*self.y_m[i]
            mpcontrol.system_output_per_sample(self.tf_adjust, self.e_m, self.theta, i)
            self.u[i] = self.x[i] * self.theta[i]
            mpcontrol.system_output_per_sample(self.tf_system, self.u, self.y, i)


def Diophantine(g, d, Am_poles, A0_poles):
    """
    Calculate numerator and denumerator of ISTR controller.

    Parameters
    ----------
    g: transfer function of the system to be controlled.
    d: the system delay parameter
    Am_poles: desired system poles.
    A0_poles: desired observer poles.

    Returns
    -------
    numerator and denumerator of ISTR controller.

    Examples
    --------
    >>> g = tf([0.0194227, 0.0182], [1, -1.81919, 0.822752], 0.0005)
    >>> R, S = Diophantine(g, 1, [0.8+0.2j, 0.8-0.2j], [0.5])

    """
    A = g.den
    B = g.num
    na = len(A) - 1
    nb = len(B) - 1
    nr = nb + d - 1
    nalpha = na + nb + d - 1
    for i in range(nb + 1):
        B[i] = B[i] / A[0]
    for i in range(1, na + 1):
        A[i] = A[i] / A[0]
    A[0] = 1
    D = mpcontrol.create_empty_matrix(nalpha, nalpha)
    if nalpha == (len(Am_poles) + len(A0_poles)):
        Am = [1]
        for i in range(len(Am_poles)):
            Am = conv(Am, [1, -Am_poles[i]])
        A0 = [1]
        for i in range(len(A0_poles)):
            if i <= nalpha - len(Am_poles):
                A0 = conv(A0, [1, -A0_poles[i]])
        Am0 = conv(Am, A0)
        for i in range(nr):
            for j in range(i, len(A) + i):
                D[j][i] = A[j - i]
        for i in range(nalpha - nr):
            for j in range(i, nb + d + i):
                D[j + d - 1][i + nr] = B[j - i]
    alpha = Am0[1:]
    a = [0 for _i in range(nalpha)]
    a[0:na] = A[1:]
    I = mpcontrol.create_empty_matrix(nalpha, 1)
    for i in range(len(a)):
        I[i][0] = alpha[i] - a[i]
    RS = mpcontrol.mult_matrices(inv(D), I)
    R = [0 for _i in range(nr + 1)]
    S = [0 for _i in range(na)]
    R[0] = 1
    for i in range(1, nr + 1):
        R[i] = RS[i - 1][0].real
    for i in range(nr, nalpha):
        S[i - nr] = RS[i][0].real
    return R, S


def get_data_from_MC(NoElements, com_port, selection):
    """
    Obtain desired number of samples for one or two float variables from selected com port.

    Parameters
    ----------
    NoElements: desired number of samples for each variable.
    com_port: selected COM port.
    selection: selection for number of variable. If it is 0, one float variable is obtained. If it is 1,
        two float variables are obtained.

    Returns
    -------
    Obtained variables. Returns one variable if selection is 0, and two variables if selection is 1.

    Examples
    --------
    >>> x1 = get_data_from_MC(2000, 'com10', 0)
    >>> x, y = get_data_from_MC(2000, 'com10', 1)

    """
    port = serial.Serial(com_port, 921600, timeout=10)

    port.write(b'r\r')

    data1 = [0 for _i in range(NoElements)]
    data1f = [0 for _i in range(NoElements)]

    if selection == 1:
        data2 = [0 for _i in range(NoElements)]
        data2f = [0 for _i in range(NoElements)]

    for x in range(NoElements):
        if selection == 0:
            data = port.read(4)
            data1[x] = data[0] | (data[1] << 8) | (data[2] << 16) | (data[3] << 24)
        elif selection == 1:
            data = port.read(4 * 2)
            data1[x] = data[0] | (data[1] << 8) | (data[2] << 16) | (data[3] << 24)
            data2[x] = data[4] | (data[5] << 8) | (data[6] << 16) | (data[7] << 24)

    port.close()

    for x in range(NoElements):
        dummy = hex(data1[x])
        dummy2 = int(dummy, 16)
        cp = pointer(c_int(dummy2))
        fp = cast(cp, POINTER(c_float))
        data1f[x] = fp.contents.value
        if selection == 1:
            dummy = hex(data2[x])
            dummy2 = int(dummy, 16)
            cp = pointer(c_int(dummy2))
            fp = cast(cp, POINTER(c_float))
            data2f[x] = fp.contents.value

    if selection == 0:
        return data1f
    elif selection == 1:
        return data1f, data2f


def get_coefficients_from_MC(NoNums, NoDens, com_port):
    """
    Obtain desired number of RLS transfer function coefficients from selected com port.

    Parameters
    ----------
    NoNums: desired number of numerator coefficient.
    NoDens: desired number of denumerator coefficient.
    com_port: selected COM port.

    Returns
    -------
    Obtained numerator and denumerator coefficients.

    Examples
    --------
    >>> x1 = get_data_from_MC(2000, 'com10', 0)
    >>> num, den = get_data_from_MC(2, 2, 'com10')

    """
    NoElements = NoNums + NoDens
    port = serial.Serial(com_port, 921600, timeout=10)
    port.write(b'r\r')
    data1 = [0 for _i in range(NoElements)]
    data1f = [0 for _i in range(NoElements)]
    for x in range(NoElements):
        data = port.read(4)
        data1[x] = data[0] | (data[1] << 8) | (data[2] << 16) | (data[3] << 24)
    port.close()
    for x in range(NoElements):
        dummy = hex(data1[x])
        dummy2 = int(dummy, 16)
        cp = pointer(c_int(dummy2))
        fp = cast(cp, POINTER(c_float))
        data1f[x] = fp.contents.value

    num = data1f[0:NoNums]
    den = data1f[NoNums:NoElements]
    return num, den


def ma_filter(signal, size):
    """
    Apply moving average filter to desired signal.

    Parameters
    ----------
    signal: signal which moving average filter is applied.
    size: window size for moving average filter.

    Returns
    ----------
    Filtered signal.

    Examples
    --------
    >>> signal = get_data_from_MC(2000, 'com10', 0)
    >>> signal_filtered = ma_filter(signal, 15)

    """
    output = [0 for _i in range(len(signal))]
    k = int((size-1)/2)
    for j in range(len(signal)):
        if j < k:
            for i in range(size):
                if (j - k + i) < 0:
                    output[j] = output[j] + signal[0]
                else:
                    output[j] = output[j] + signal[j - k + i]
        elif k <= j < (len(signal) - k):
            for i in range(size):
                output[j] = output[j] + signal[j - k + i]
        else:
            for i in range(size):
                if (j - k + i) > (len(signal)-1):
                    output[j] = output[j] + signal[len(signal)-1]
                else:
                    output[j] = output[j] + signal[j - k + i]
    for j in range(len(signal)):
        output[j] = output[j] / size
    return output


def save_to_file(data, file_name='data.dat'):
    """
    Save data to file.

    Parameters
    ----------
    data: data to be saved.
    file_name: name of the file. Default is 'data.dat'.

    Examples
    --------
    >>> data = [0.5*n for n in range(500)]
    >>> save_to_file(data, 'dummy.dat')

    """
    file = open(file_name, 'w')
    for x in range(len(data)):
        file.write(str(data[x]))
        file.write("\n")
    file.close()


def read_from_file(file_name='data.dat'):
    """
    Read data from file.

    Parameters
    ----------
    file_name: name of the file. Default is 'data.dat'.

    Examples
    --------
    >>> data = read_from_file('dummy.dat')

    """
    f = open(file_name, 'r')
    x = f.readlines()
    y = [0.0 for _i in range(len(x))]
    for n in range(len(x)):
        y[n] = float(x[n])

    f.close()
    return y


def stepinfo_from_signal(y, ts, ST=0.02, RT=None):
    """
    Find step response parameters from step output signal and print the results.

    Parameters
    ----------
    y : output signal obtained from step input.
    ts : sampling time.
    ST : threshold value for defining settling time. Default is 0.02.
    RT : threshold values for defining rise time. Default is [0.1 0.9].

    Returns
    -------
    rise_time: time value for step response to rise from (100*RT[0])% to (100*RT[1])% of the steady-state value.
    settling_time: time value for the error between step response and steady-state value to fall within
        (100*ST)% of steady-state value.
    settling_min: minimum value of step response after rise-time is reached.
    settling_max: maximum value of step response after rise-time is reached.
    overshoot: percentage overshoot.
    undershoot: percentage undershoot.
    peak: peak absolute value of step response.
    peak_time: time value at which the peak occurs.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> y = step(tf_discrete, nsample=500)
    >>> stepinfo_from_signal(y, 0.01)

    """
    if RT is None:
        RT = [0.1, 0.9]
    settling_time = 0
    settling_min = 0

    y_final = 0
    size = len(y)
    size_avg = int(size / 20)
    for n in range(size - size_avg, size, 1):
        y_final = y_final + y[n]
    y_final = y_final / size_avg

    st_flag = 0
    for n in range(size):
        if st_flag == 0:
            if (y[n] < (y_final * (1 + ST))) and (y[n] > (y_final * (1 - ST))):
                st_flag = 1
                settling_time = float(n) * ts
        else:
            if (y[n] >= (y_final * (1 + ST))) or (y[n] <= (y_final * (1 - ST))):
                st_flag = 0

    if settling_time == 0:
        peak = float('INF')
        peak_time = float('INF')
        overshoot = float('NAN')
        undershoot = float('NAN')
        rise_time = float('NAN')
        settling_time = float('NAN')
        settling_min = float('NAN')
        settling_max = float('NAN')

    else:
        y_max = max(y)
        y_max_index = y.index(y_max)
        peak = abs(y_max)
        peak_time = (float(y_max_index)) * ts
        settling_max = y_max
        if y_max > y_final:
            overshoot = 100 * (y_max - y_final) / y_final
        else:
            overshoot = 0.0
        y_min = min(y)
        if y_min < 0:
            undershoot = 100 * abs(y_min) / y_final
        else:
            undershoot = 0.0

        t_high = 0.0
        t_low = 0.0
        for n in range(y_max_index):
            if (y[n] - y[0]) >= ((y_final - y[0]) * RT[0]):
                t_low = float(n) * ts
                break
        for n in range(y_max_index):
            if (y[n] - y[0]) >= ((y_final - y[0]) * RT[1]):
                t_high = float(n) * ts
                settling_min = y[n]
                break
        rise_time = t_high - t_low
        if rise_time == 0:
            settling_min = y_max
        for n in range(y_max_index, size, 1):
            if y[n] < settling_min:
                settling_min = y[n]

    print('\n' + "STEP INFO:" + '\n')
    buf = "Rise Time (T_rt): %g sec" % rise_time
    print(buf)
    buf = "Settling Time (T_st): %g sec" % settling_time
    print(buf)
    buf = "Peak Time (T_pt): %g sec" % peak_time
    print(buf)
    buf = "Overshoot (M_p): %g%%" % overshoot
    print(buf)
    buf = "Undershoot: %g%%" % undershoot
    print(buf)
    buf = "Settling Min.: %g" % settling_min
    print(buf)
    buf = "Settling Max.: %g" % settling_max
    print(buf)
    buf = "Peak: %g\n" % peak
    print(buf)
    return rise_time, settling_time, settling_min, settling_max, overshoot, undershoot, peak, peak_time


def compare_stepinfo(data, tfd, gain):
    """
    Compare step responses obtained from real system output and discrete-time transfer function.

    Parameters
    ----------
    data: output of the real system when step input is applied.
    tfd: discrete-time transfer function of the real system.
    gain: amplitude of the step input.

    Examples
    --------
    >>> data1 = [0.5*i for i in range(500)]
    >>> data2 = [0.25*i for i in range(500)]

    """
    print("Step info from transfer function simulation")
    tfd.stepinfo(stepgain=gain)
    print("Step info from real time signal output")
    stepinfo_from_signal(data, tfd.sampling_time)
    y = step(tfd, stepgain=gain, nsample=len(data))
    return y


def transpose(x):
    """
    Transpose the input matrix.

    Parameters
    ----------
    x: input matrix.

    Returns
    -------
    Transposed matrix.

    Examples
    --------
    >>> x = mpcontrol.create_empty_matrix(3, 2)
    >>> y = transpose(x)

    """
    y = mpcontrol.create_empty_matrix(len(x[0]), len(x))
    for n in range(len(x[0])):
        for k in range(len(x)):
            y[n][k] = x[k][n]
    return y


def create_identity_matrix(dim):
    """
    Create identity matrix.

    Parameters
    ----------
    dim: dimension of the dim*dim identity matrix.

    Returns
    -------
    Identity matrix.

    Examples
    --------
    >>> x = create_identity_matrix(3)

    """
    x = mpcontrol.create_empty_matrix(dim, dim)
    for n in range(dim):
        x[n][n] = 1
    return x


def minor(x, i, j):
    """
    Calculate minor of matrix.

    Parameters
    ----------
    x: original matrix.

    Returns
    -------
    minor of original matrix.

    Examples
    --------
    >>> F = [[1.8188, -0.8224], [1.0, 0.0]]
    >>> result = det(F, 0, 0)

    """
    y = mpcontrol.create_empty_matrix(len(x[0]) - 1, len(x[0]) - 1)
    for n in range(len(x[0])-1):
        if n < i:
            for k in range(len(x[0])-1):
                if k < j:
                    y[n][k] = x[n][k]
                else:
                    y[n][k] = x[n][k + 1]
        else:
            for k in range(len(x[0]) - 1):
                if k < j:
                    y[n][k] = x[n + 1][k]
                else:
                    y[n][k] = x[n + 1][k + 1]
    return y


def det(x):
    """
    Calculate determinant of matrix.

    Parameters
    ----------
    x: original matrix.

    Returns
    -------
    determinant of original matrix.

    Examples
    --------
    >>> F = [[1.8188, -0.8224], [1.0, 0.0]]
    >>> result = det(F)

    """
    if len(x[0]) == 1:
        return x[0][0]
    elif len(x[0]) == 2:
        return x[0][0]*x[1][1]-x[0][1]*x[1][0]
    else:
        determinant = 0
        for n in range(len(x)):
            determinant += ((-1)**n)*x[n][0]*det(minor(x, n, 0))
        return determinant


def inv(x):
    """
    Calculate inverse of matrix.

    Parameters
    ----------
    x: original matrix.

    Returns
    -------
    inverse of original matrix.

    Examples
    --------
    >>> F = [[1.8188, -0.8224], [1.0, 0.0]]
    >>> result = inv(F)

    """
    if len(x) == len(x[0]):
        determinant = det(x)
        if determinant != 0:
            if len(x[0]) == 1:
                return [[1/determinant]]
            elif len(x[0]) == 2:
                return [[x[1][1] / determinant, -1 * x[0][1] / determinant],
                        [-1 * x[1][0] / determinant, x[0][0] / determinant]]
            cofactors = []
            for n in range(len(x[0])):
                cofactorRow = []
                for k in range(len(x[0])):
                    minor_arr = minor(x, n, k)
                    cofactorRow.append(((-1) ** (n + k)) * det(minor_arr))
                cofactors.append(cofactorRow)
            cofactors = transpose(cofactors)
            for n in range(len(cofactors)):
                for k in range(len(cofactors)):
                    cofactors[n][k] = cofactors[n][k] / determinant
        else:
            raise ValueError("Matrix inverse does not exist, matrix determinant equals to 0.")
    else:
        raise ValueError("Matrix inverse does not exist, matrix is not square.")
    return cofactors


def get_roots(coeffs):
    """
    Find the roots of a polynomial using its coefficients.

    Parameters
    ----------
    coeffs: coefficients of the polynomial.

    Returns
    -------
    Roots of polynomial.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> poles = get_roots(tf_discrete.den)
    >>> zeros = get_roots(tf_discrete.num)

    """
    order = len(coeffs) - 1
    a = [0.0 for _i in range(order + 3)]
    b = [0.0 for _i in range(order + 3)]
    c = [0.0 for _i in range(order + 3)]
    roots = [0.0 for _i in range(order)]
    eps = 1.0e-13
    cnt = 0

    for n in range(order + 1):
        a[n] = coeffs[order - n]
    if order == 0:
        roots = []
    elif order == 1:
        roots[0] = root_first_order(a[1], a[0])
    elif order == 2:
        roots[1], roots[0] = roots_second_order(a[2], a[1], a[0])
    else:
        while order >= 3:
            r = 0.1
            s = 0.1
            eps_r = 1
            eps_s = 1
            while (eps_r > eps) and (eps_s > eps) and (cnt < 50000):
                cnt = cnt + 1
                b[order] = a[order]
                b[order - 1] = a[order - 1] + r * b[order]
                for n in range(order - 2, -1, -1):
                    b[n] = a[n] + r * b[n + 1] + s * b[n + 2]
                c[order] = b[order]
                c[order - 1] = b[order - 1] + r * c[order]
                for n in range(order - 2, 0, -1):
                    c[n] = b[n] + r * c[n + 1] + s * c[n + 2]
                dr = ((b[1] * c[2]) - (b[0] * c[3])) / ((c[1] * c[3]) - (c[2] * c[2]))
                ds = ((b[0] * c[2]) - (b[1] * c[1])) / ((c[1] * c[3]) - (c[2] * c[2]))
                r = r + dr
                s = s + ds
                if r == 0:
                    r = eps
                if s == 0:
                    s = eps
                eps_r = 100 * abs(dr / r)
                eps_s = 100 * abs(ds / s)
            roots[order - 1], roots[order - 2] = roots_second_order(1.0, -r, -s)
            for n in range(order + 1):
                a[n] = b[n]
            for n in range(order + 1):
                a[n] = a[n + 2]
            order = order - 2
        if order == 2:
            roots[1], roots[0] = roots_second_order(a[2], a[1], a[0])
        else:
            roots[0] = root_first_order(a[1], a[0])
    return roots


def root_first_order(a1, a0):
    """
    Find the root of first degree polynomial.

    Parameters
    ----------
    a1: coefficient of x^1.
    a0: coefficient of x^0.

    Returns
    -------
    Root of first degree polynomial.

    Examples
    --------
    >>> den = [1, 0.6]
    >>> root = root_first_order(den[0], den[1])

    """
    real = -a0 / a1
    imag = 0
    root = complex(real, imag)
    return root


def roots_second_order(a2, a1, a0):
    """
    Find the roots of second degree polynomial.

    Parameters
    ----------
    a2: coefficient of x^2.
    a1: coefficient of x^1.
    a0: coefficient of x^0.

    Returns
    -------
    Roots of second degree polynomial.

    Examples
    --------
    >>> den = [1, 1.4, 0.9]
    >>> root1, root2 = roots_second_order(den[0], den[1], den[2])

    """
    disc = a1 * a1 - 4 * a2 * a0
    if disc > 0:
        real1 = (-a1 + sqrt(disc)) / (2 * a2)
        real2 = (-a1 - sqrt(disc)) / (2 * a2)
        imag1 = 0
        imag2 = 0
        root1 = complex(real1, imag1)
        root0 = complex(real2, imag2)
    elif disc == 0:
        real1 = -a1 / (2 * a2)
        real2 = -a1 / (2 * a2)
        imag1 = 0
        imag2 = 0
        root1 = complex(real1, imag1)
        root0 = complex(real2, imag2)
    else:
        real1 = -a1 / (2 * a2)
        real2 = -a1 / (2 * a2)
        imag1 = sqrt(-disc) / (2 * a2)
        imag2 = -sqrt(-disc) / (2 * a2)
        root1 = complex(real1, imag1)
        root0 = complex(real2, imag2)
    return root1, root0


def find_pole_magnitudes(poles):
    """
    Find the pole magnitudes for given set of poles.

    Parameters
    ----------
    poles: system poles.

    Returns
    -------
    Magnitudes of given poles.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> poles = get_roots(tf_discrete.den)
    >>> pole_mags = find_pole_magnitudes(poles)

    """
    pole_mags = [0.0 for _i in range(len(poles))]
    for n in range(len(poles)):
        pole_mags[n] = abs(poles[n])
    return pole_mags


def is_system_stable(poles):
    """
    Check if the discrete-time system is stable or not for given set of poles.

    Returns
    -------
    True: if the discrete-time system is stable for all poles.
    False: if the discrete-time system is not stable for all poles.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> poles = tf_discrete.pole()
    >>> check = is_system_stable(poles)

    """
    pole_mags = find_pole_magnitudes(poles)
    for n in range(len(poles)):
        if pole_mags[n] > 1:
            return False
    return True


def find_pole_dampings(poles):
    """
    Find the pole damping ratios for given set of discrete-time poles.

    Parameters
    ----------
    poles: discrete-time system poles.

    Returns
    -------
    Damping ratios of given discrete-time poles.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> poles = get_roots(tf_discrete.den)
    >>> pole_damps = find_pole_dampings(poles)

    """
    pole_damps = [0.0 for _i in range(len(poles))]
    for n in range(len(poles)):
        ln_pole_real = log(abs(poles[n]))
        ln_pole_imag = atan2(poles[n].imag, poles[n].real)
        pole_damps[n] = -cos(atan2(ln_pole_imag, ln_pole_real))
        if pole_damps[n] < 0:
            pole_damps[n] = 0
    return pole_damps


def find_pole_naturalfreqs(poles, ts):
    """
    Find the pole natural frequencies for given set of discrete-time poles.

    Parameters
    ----------
    poles: discrete-time system poles.
    ts: sampling time for discrete-time poles.

    Returns
    -------
    Natural frequencies of given discrete-time poles.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> poles = get_roots(tf_discrete.den)
    >>> pole_naturalfreqs = find_pole_naturalfreqs(poles, 0.01)

    """
    pole_naturalfreqs = [0.0 for _i in range(len(poles))]
    for n in range(len(poles)):
        ln_pole_real = log(abs(poles[n])) / ts
        ln_pole_imag = atan2(poles[n].imag, poles[n].real) / ts
        pole_naturalfreqs[n] = sqrt((ln_pole_real * ln_pole_real) + (ln_pole_imag * ln_pole_imag))
    return pole_naturalfreqs


def find_pole_timeconstants(poles, ts):
    """
    Find the pole time constants for given set of discrete-time poles.

    Parameters
    ----------
    poles: discrete-time system poles.
    ts: sampling time for discrete-time poles.

    Returns
    -------
    Time constants of given discrete-time poles.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> poles = get_roots(tf_discrete.den)
    >>> pole_timeconstants = find_pole_timeconstants(poles, 0.01)

    """
    pole_damps = find_pole_dampings(poles)
    pole_naturalfreqs = find_pole_naturalfreqs(poles, ts)
    pole_timeconstants = [0.0 for _i in range(len(pole_damps))]
    for n in range(len(pole_damps)):
        if pole_naturalfreqs[n] == 0 or pole_damps[n] == 0:
            pole_timeconstants[n] = float('INF')
        else:
            pole_timeconstants[n] = 1 / (pole_damps[n] * pole_naturalfreqs[n])
    return pole_timeconstants


def find_ss_time(tfd):
    """
    Find the estimated number of samples needed for a discrete-time system reaches its steady state.

    Parameters
    ----------
    tfd: discrete-time transfer function.

    Returns
    -------
    The estimated number of samples needed for a discrete-time system reaches its steady state.

    Examples
    --------
    >>> tf_discrete = tf([3], [1, 0.2], 0.01)
    >>> ss_time = find_ss_time(tf_discrete)

    """
    poles = get_roots(tfd.den)
    pole_timeconstants = find_pole_timeconstants(poles, tfd.sampling_time)
    for n in range(len(poles)):
        for k in range(len(poles)):
            if pole_timeconstants[k] < pole_timeconstants[n]:
                temp = pole_timeconstants[n]
                pole_timeconstants[n] = pole_timeconstants[k]
                pole_timeconstants[k] = temp
    if pole_timeconstants[0] == float('INF') or pole_timeconstants[0] > 1.e3:
        result = int(pole_timeconstants[1] * 20 / tfd.sampling_time)
    else:
        result = int(pole_timeconstants[0] * 20 / tfd.sampling_time)
    return result


def poly_value(arr, x):
    """
    Find the result of polynomial at desired x value.

    Parameters
    ----------
    arr: coefficients of the polynomial.
    x: desired calculation point.

    Returns
    -------
    Result of the polynomial at x.

    Examples
    --------
    >>> den = [1, 0.2]
    >>> y = poly_value(den, 1)

    """
    result = 0
    for n in range(len(arr)):
        result = result * x + arr[n]
    return result


def array_sum(arr, k):
    """
    Find the sum of array elements starting from first element to kth element.

    Parameters
    ----------
    arr: array of interest.
    k: desired index value.

    Returns
    -------
    Result of the summation.

    Examples
    --------
    >>> den = [1, 0.2, 1.6, 2.54, 3.6]
    >>> result = array_sum(den, 2)

    """
    if k < 0:
        return 0
    else:
        return arr[k] + array_sum(arr, k - 1)


def get_time_near(t, y, y_desired):
    """
    Find the time value when desired y value is obtained.

    Parameters
    ----------
    t : time array contains time steps.
    y: amplitude array contains y values for given t array.
    y_desired: desired y value.

    Returns
    -------
    Time value when y_desired is obtained.

    Examples
    --------
    >>> time_arr = [0.01*n for n in range(1000)]
    >>> y_arr = [0.5*n for n in range(1000)]
    >>> get_time_near(time_arr, y_arr, 3.86)

    """
    my_time = 0
    for i in range(len(y)):
        if y[i] >= y_desired >= y[i-1]:
            my_time = t[i-1]+(((y_desired-y[i-1])*(t[i]-t[i-1]))/(y[i]-y[i-1]))
            break
    return my_time


def conv(u, v):
    """
    Return the convolution of two arrays.

    Parameters
    ----------
    u : first array.
    v: second array.

    Returns
    -------
    Result of the convolution operation.

    Examples
    --------
    >>> x = [1 for n in range(10)]
    >>> y = [2 for n in range(10)]
    >>> z = conv(x, y)

    """
    m = len(u)
    n = len(v)
    l = m + n -1
    w = [0 for _i in range(l)]
    for k in range(l):
        for j in range(k+1):
            if j >= m:
                u_dum = 0
            else:
                u_dum = u[j]
            if k-j >= n:
                v_dum = 0
            else:
                v_dum = v[k-j]
            w[k] += u_dum * v_dum
    return w

