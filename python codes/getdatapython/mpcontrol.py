from math import cos, pi, sin
from random import random


class ss_simple:
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
        >>> ss_continuous = ss_simple([[1.4, 1.9], [1.0, 0.5]], [[1.0], [0.0]], [[0.8, 0.2]], [[0.02]])
        >>> ss_discrete = ss_simple([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)

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
        >>> ss_discrete = ss_simple([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)
        >>> check = ss_discrete.is_dtime()

        """
        if self.sampling_time is None or self.sampling_time == 0:
            return False
        else:
            return True


class tf_simple:
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
        >>> tf_continuous = tf_simple([10], [1, 10])
        >>> tf_discrete = tf_simple([3], [1, -0.2], 0.01)

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
        >>> tf_discrete = tf_simple([3], [1, -0.2], 0.01)
        >>> check = tf_discrete.is_dtime()

        """
        if self.sampling_time is None or self.sampling_time == 0:
            return False
        else:
            return True


def iir_filter(b, a, x, size):
    """
    Find the system output using the difference equation.

    Parameters
    ----------
    b: numerator coefficients of discrete-time system.
    a: denominator coefficients of discrete-time system.
    x: input signal.
    size: desired number of samples for output signal.

    Returns
    -------
    System output.

    Examples
    --------
    >>> num = 3
    >>> den = [1, 0.2]
    >>> x = [1 for i in range(1000)]
    >>> y = iir_filter(num, den, x, 1000)

    """
    size2 = len(b)
    size3 = len(a)
    y = [0.0 for _i in range(size)]
    num = [0.0 for _i in range(size2)]
    den = [0.0 for _i in range(size3)]
    for n in range(1, size3, 1):
        den[n] = a[n] / a[0]
    for n in range(size2):
        num[n] = b[n] / a[0]
    den[0] = 1

    for n in range(size):
        sum_iir = 0
        for k in range(size3):
            if (size3 - size2) <= k < (n + 1):
                sum_iir += num[k - (size3 - size2)] * x[n - k]
        for m in range(1, size3, 1):
            if m < n + 1:
                sum_iir -= den[m] * y[n - m]
        y[n] = sum_iir
    return y


def lsim(tfd, x, nsample=None):
    """
    Simulate the output of a discrete-time transfer function for given input.

    Parameters
    ----------
    tfd: discrete-time transfer function.
    x: input signal.
    nsample: desired number of output samples. Default is None, which means size
        of the input signal is used to generate the output signal.

    Returns
    -------
    Calculated output signal.

    Raises
    ------
    NotImplementedError
        if transfer function is not discrete-time.
    ValueError("Sample size must be nonzero value!")
        if nsample is selected as 0
    ValueError("Size of the output cannot be larger than size of the input!")
        if nsample is selected larger than the size of input signal

    Examples
    --------
    >>> tf_discrete = tf_simple([3], [1, 0.2], 0.01)
    >>> x = [1 for i in range(1000)]
    >>> y = lsim(tf_discrete, x, 400)

    """
    if not tfd.is_dtime():
        raise NotImplementedError("System must be discrete-time!")
    else:
        if nsample is None:
            length = len(x)
        elif nsample == 0:
            raise ValueError("Sample size must be nonzero value!")
        elif nsample > len(x):
            raise ValueError("Size of the output cannot be larger than size of the input!")
        else:
            length = nsample
        y = iir_filter(tfd.num, tfd.den, x, length)
        return y


def unit_pulse_sample(k):
    """
    Calculate the kth sample of unit pulse signal.

    Parameters
    ----------
    k: index of the desired sample.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = unit_pulse_sample(50)

    """
    if k == 0:
        sample = 1
    else:
        sample = 0
    return sample


def unit_pulse_signal(N):
    """
    Calculate the unit pulse signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.

    Returns
    -------
    Unit pulse signal.

    Examples
    --------
    >>> y = unit_pulse_signal(200)

    """
    signal = [unit_pulse_sample(n) for n in range(N)]
    return signal


def step_sample(Amp):
    """
    Calculate a sample of step signal.

    Parameters
    ----------
    Amp: amplitude of step signal.

    Returns
    -------
    Calculated sample value.

    Examples
    --------
    >>> y = step_sample(2)

    """
    sample = Amp
    return sample


def step_signal(N, Amp):
    """
    Calculate the step signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: amplitude of step signal.

    Returns
    -------
    Step signal.

    Examples
    --------
    >>> y = step_signal(200, 2)

    """
    signal = [step_sample(Amp) for _i in range(N)]
    return signal


def ramp_sample(k, Amp):
    """
    Calculate the kth sample of ramp signal.

    Parameters
    ----------
    k: index of the desired sample.
    Amp: slope of the ramp signal.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = ramp_sample(50, 2)

    """
    sample = k*Amp
    return sample


def ramp_signal(N, Amp):
    """
    Calculate the ramp signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: slope of the ramp signal.

    Returns
    -------
    Ramp signal.

    Examples
    --------
    >>> y = ramp_signal(200, 2)

    """
    signal = [ramp_sample(n, Amp) for n in range(N)]
    return signal


def parabolic_sample(k, Amp):
    """
    Calculate the kth sample of parabolic signal.

    Parameters
    ----------
    k: index of the desired sample.
    Amp: slope of the parabolic signal.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = parabolic_sample(50, 2)

    """
    sample = k*k*Amp/2
    return sample


def parabolic_signal(N, Amp):
    """
    Calculate the parabolic signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: slope of the parabolic signal.

    Returns
    -------
    Parabolic signal.

    Examples
    --------
    >>> y = parabolic_signal(200, 2)

    """
    signal = [parabolic_sample(n, Amp) for n in range(N)]
    return signal


def exponential_sample(k, a):
    """
    Calculate the kth sample of exponential signal.

    Parameters
    ----------
    k: index of the desired sample.
    a: rate of the exponential signal.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = exponential_sample(50, 0.02)

    """
    e = 2.718281828459045
    try:
        sample = e**(k*a)
    except OverflowError:
        sample = float('inf')
    return sample


def exponential_signal(N, a):
    """
    Calculate the exponential signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    a: rate of the exponential signal.

    Returns
    -------
    Exponential signal.

    Examples
    --------
    >>> y = exponential_signal(200, 0.02)

    """
    signal = [exponential_sample(n, a) for n in range(N)]
    return signal


def sinusoidal_sample(k, Amp, Freq, Phase, Offset, Fsample, select):
    """
    Calculate the kth sample of sinusoidal signal.

    Parameters
    ----------
    k: index of the desired sample.
    Amp: amplitude of the sinusoidal signal.
    Freq: frequency of the sinusoidal signal in Hz.
    Phase: phase of the sinusoidal signal in radians.
    Offset: offset of the sinusoidal signal.
    Fsample: sampling frequency of the sinusoidal signal in Hz.
    select: selection of sinusoidal type. It must be selected 0 for sine wave and 1 for cosine wave.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = sinusoidal_sample(50, 1, 50, 0.2, 0, 1000, 0)

    """
    sample = 0
    if select == 0:
        sample = Amp*sin((2*pi*Freq*k/Fsample)+Phase) + Offset
    elif select == 1:
        sample = Amp*cos((2*pi*Freq*k/Fsample)+Phase) + Offset
    return sample


def sinusoidal_signal(N, Amp, Freq, Phase, Offset, Fsample, select):
    """
    Calculate the sinusoidal signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: amplitude of the sinusoidal signal.
    Freq: frequency of the sinusoidal signal in Hz.
    Phase: phase of the sinusoidal signal in radians.
    Offset: offset of the sinusoidal signal.
    Fsample: sampling frequency of the sinusoidal signal in Hz.
    select: selection of sinusoidal type. It must be selected 0 for sine wave and 1 for cosine wave.

    Returns
    -------
    Sinusoidal signal.

    Examples
    --------
    >>> y = sinusoidal_signal(200, 1, 50, 0.2, 0, 1000, 0)

    """
    signal = [sinusoidal_sample(n, Amp, Freq, Phase, Offset, Fsample, select) for n in range(N)]
    return signal


def damped_sinusoidal_sample(k, Amp_exp, Amp_sin, Freq, Phase, Offset, Fsample, select):
    """
    Calculate the kth sample of damped sinusoidal signal.

    Parameters
    ----------
    k: index of the desired sample.
    Amp_exp: rate of the exponential decay.
    Amp_sin: amplitude of the damped sinusoidal signal.
    Freq: frequency of the damped sinusoidal signal in Hz.
    Phase: phase of the damped sinusoidal signal in radians.
    Offset: offset of the damped sinusoidal signal.
    Fsample: sampling frequency of the damped sinusoidal signal in Hz.
    select: selection of sinusoidal type. It must be selected 0 for sine wave and 1 for cosine wave.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = damped_sinusoidal_sample(50, -0.02, 1, 50, 0.2, 0, 1000, 0)

    """
    sample = exponential_sample(k, Amp_exp)*sinusoidal_sample(k, Amp_sin, Freq, Phase, Offset, Fsample, select)
    return sample


def damped_sinusoidal_signal(N, Amp_exp, Amp_sin, Freq, Phase, Offset, Fsample, select):
    """
    Calculate the damped sinusoidal signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp_exp: rate of the exponential decay.
    Amp_sin: amplitude of the damped sinusoidal signal.
    Freq: frequency of the damped sinusoidal signal in Hz.
    Phase: phase of the damped sinusoidal signal in radians.
    Offset: offset of the damped sinusoidal signal.
    Fsample: sampling frequency of the damped sinusoidal signal in Hz.
    select: selection of sinusoidal type. It must be selected 0 for sine wave and 1 for cosine wave.

    Returns
    -------
    Damped sinusoidal signal.

    Examples
    --------
    >>> y = damped_sinusoidal_signal(200, -0.02, 1, 50, 0.2, 0, 1000, 0)

    """
    signal = [damped_sinusoidal_sample(n, Amp_exp, Amp_sin, Freq, Phase, Offset, Fsample, select) for n in range(N)]
    return signal


def rectangular_sample(k, Amp, Period, Duty, Offset):
    """
    Calculate the kth sample of rectangular signal.

    Parameters
    ----------
    k: index of the desired sample.
    Amp: amplitude of the rectangular signal.
    Period: period of the rectangular signal in samples.
    Duty: Percentage duty cycle of the rectangular signal
    Offset: offset of the rectangular signal.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = rectangular_sample(50, 1, 20, 50, 0)

    """
    duty_cnt = int(Period*Duty/100)
    cnt_mod = k % Period
    if cnt_mod < duty_cnt:
        level = 1
    else:
        level = 0
    sample = level * Amp + Offset
    return sample


def rectangular_signal(N, Amp, Period, Duty, Offset):
    """
    Calculate the rectangular signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: amplitude of the rectangular signal.
    Period: period of the rectangular signal in samples.
    Duty: Percentage duty cycle of the rectangular signal
    Offset: offset of the rectangular signal.

    Returns
    -------
    Rectangular signal.

    Examples
    --------
    >>> y = rectangular_signal(200, 1, 20, 50, 0)

    """
    signal = [rectangular_sample(n, Amp, Period, Duty, Offset) for n in range(N)]
    return signal


def sum_of_sinusoids_sample(k, No_Sines, Amps, Freqs, Phases, Offsets, Fsample, select):
    """
    Calculate the kth sample of sum of sinusoids.

    Parameters
    ----------
    k: index of the desired sample.
    No_Sines: number of sinusoids to sum.
    Amps: amplitudes of the sinusoids.
    Freqs: frequencies of the sinusoids in Hz.
    Phases: phases of the sinusoids in radians.
    Offsets: offsets of the sinusoids.
    Fsample: sampling frequency of the sinusoids in Hz.
    select: selection of sinusoidal type. It must be selected 0 for sine wave and 1 for cosine wave.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> amps = [0.2, 0.4, 0.5]
    >>> freqs = [50, 80, 120]
    >>> phases = [0.1, 0, 0.4]
    >>> offsets = [0, 0.5, 1]
    >>> y = sum_of_sinusoids_sample(50, 3, amps, freqs, phases, offsets, 1000, 0)

    """
    sample = 0
    if select == 0:
        for n in range(No_Sines):
            sample = sample + Amps[n]*sin((2*pi*Freqs[n]*k/Fsample)+Phases[n]) + Offsets[n]
    elif select == 1:
        sample = 0
        for n in range(No_Sines):
            sample = sample + Amps[n]*cos((2*pi*Freqs[n]*k/Fsample)+Phases[n]) + Offsets[n]
    return sample


def sum_of_sinusoids_signal(N, No_Sines, Amps, Freqs, Phases, Offsets, Fsample, select):
    """
    Calculate the sum of sinusoids for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    No_Sines: number of sinusoids to sum.
    Amps: amplitudes of the sinusoids.
    Freqs: frequencies of the sinusoids in Hz.
    Phases: phases of the sinusoids in radians.
    Offsets: offsets of the sinusoids.
    Fsample: sampling frequency of the sinusoids in Hz.
    select: selection of sinusoidal type. It must be selected 0 for sine wave and 1 for cosine wave.

    Returns
    -------
    Sum of sinusoids.

    Examples
    --------
    >>> amps = [0.2, 0.4, 0.5]
    >>> freqs = [5, 25, 100]
    >>> phases = [0.1, 0, 0.4]
    >>> offsets = [0, 0.5, 1]
    >>> y = sum_of_sinusoids_signal(200, 3, amps, freqs, phases, offsets, 1000, 0)

    """
    signal = [sum_of_sinusoids_sample(n, No_Sines, Amps, Freqs, Phases, Offsets, Fsample, select) for n in range(N)]
    return signal


def sweep_sample(cnt, Amp, Freq_start, Freq_stop, Offset, Fsample, Duration):
    """
    Calculate the kth sample of sweep signal.

    Parameters
    ----------
    cnt: index of the desired sample.
    Amp: amplitude of the sweep signal.
    Freq_start: starting frequency of the sweep signal in Hz.
    Freq_stop: stopping frequency of the sweep signal in Hz.
    Offset: offset of the sweep signal.
    Fsample: sampling frequency of the sweep signal in Hz.
    Duration: duration of the sweep signal from starting to stopping frequency in seconds.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = sweep_signal(50, 1, 1, 10, 0, 1000, 2)

    """
    duration_cnt = int(Duration * Fsample)
    k = cnt % duration_cnt
    sample = Amp*(sin((2*pi*Freq_start*k/Fsample)+(2*pi*(Freq_stop-Freq_start)*k*k/(2*(duration_cnt-1)*Fsample))))+Offset
    return sample


def sweep_signal(N, Amp, Freq_start, Freq_stop, Offset, Fsample, Duration):
    """
    Calculate the sweep signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: amplitude of the sweep signal.
    Freq_start: starting frequency of the sweep signal in Hz.
    Freq_stop: stopping frequency of the sweep signal in Hz.
    Offset: offset of the sweep signal.
    Fsample: sampling frequency of the sweep signal in Hz.
    Duration: duration of the sweep signal from starting to stopping frequency in seconds.

    Returns
    -------
    Sweep signal.

    Examples
    --------
    >>> y = sweep_signal(2000, 1, 1, 10, 0, 1000, 2)

    """
    signal = [sweep_sample(n, Amp, Freq_start, Freq_stop, Offset, Fsample, Duration) for n in range(N)]
    return signal


sample_rand = 0
def random_sample(k, Amp, Duration, Offset):
    """
    Calculate the kth sample of random signal.

    Parameters
    ----------
    k: index of the desired sample.
    Amp: maximum amplitude of the random signal.
    Duration: number of samples for random signal to remain unchanged.
    Offset: offset of the random signal.

    Returns
    -------
    Calculated kth sample value.

    Examples
    --------
    >>> y = random_sample(50, 1, 5, 0)

    """
    global sample_rand
    mod_cnt = k % Duration
    if mod_cnt == 0:
        sample_rand = (Amp * random()) + Offset
    return sample_rand


def random_signal(N, Amp, Duration, Offset):
    """
    Calculate the random signal for desired number of samples.

    Parameters
    ----------
    N: number of samples.
    Amp: maximum amplitude of the random signal.
    Duration: number of samples for random signal to remain unchanged.
    Offset: offset of the random signal.

    Returns
    -------
    Random signal.

    Examples
    --------
    >>> y = random_sample(200, 1, 5, 0)

    """
    signal = [random_sample(n, Amp, Duration, Offset) for n in range(N)]
    return signal


def system_output_per_sample(tfd, x, y, k):
    """
    Calculate the kth output sample of the discrete-time transfer function.

    Parameters
    ----------
    tfd: discrete-time transfer function.
    x: input signal array.
    y: output signal array.
    k: index of the desired output sample.

    Returns
    -------
    Calculated kth output sample value.

    Examples
    --------
    >>> tf_discrete = tf_simple([0.6], [1, 0.8], 0.0001)
    >>> x = [1, 1, 1, 1, 1]
    >>> y = [0, 0.6, 0.12, 0.504, 0]
    >>> system_output_per_sample(tf_discrete, x, y, 4)

    """
    sum_iir = 0
    for n in range(tfd.num_order+1):
        if k < (tfd.num_order+1):
            if k >= n:
                sum_iir += tfd.num[n] * x[k - n]
        else:
            k_mod = k % (tfd.num_order + 1)
            if (k_mod-n) >= 0:
                sum_iir += tfd.num[n] * x[k_mod - n]
            else:
                sum_iir += tfd.num[n] * x[k_mod - n + (tfd.num_order + 1)]
    for n in range(1, tfd.den_order+1):
        if k < (tfd.den_order+1):
            if k >= n:
                sum_iir -= tfd.den[n] * y[k - n]
        else:
            k_mod = k % (tfd.den_order + 1)
            if (k_mod-n) >= 0:
                sum_iir -= tfd.den[n] * y[k_mod - n]
            else:
                sum_iir -= tfd.den[n] * y[k_mod - n + (tfd.den_order + 1)]
    result = sum_iir/tfd.den[0]
    y[k % (tfd.den_order + 1)] = result
    return result


def lsim_ss(ssd, u, nsample=None, select=0, x0=None):
    """
    Simulate the output or states of a discrete-time state-space model for given input and initial states.

    Parameters
    ----------
    ssd: discrete-time state space model.
    u: input signal.
    nsample: desired number of output or state samples. Default is None, which means size
        of the input signal is used to generate the output or state signal.
    select: selection for simulation result. If it is selected as 0, the output is given as a result.
        If it is selected as other than 0, it is used to select a state. For example if it is selected
        as 1, first state is given as a result. If it is selected as 2, second state is given as a result.
        Default is 0.
    x0: initial states of the state space model. Default is None, which means initial states are zero.

    Returns
    -------
    Selected simulation result.

    Raises
    ------
    NotImplementedError
        if transfer function is not discrete-time.
    ValueError("Sample size must be nonzero value!")
        if nsample is selected as 0
    ValueError("Size of the output cannot be larger than size of the input!")
        if nsample is selected larger than the size of input signal

    Examples
    --------
    >>> ss_discrete = ss_simple([[1.9, -0.9], [1.0, 0.0]], [[1.0], [0.0]], [[0.012, 0.0003]], [[0.003]], 0.0005)
    >>> x = [1 for i in range(1000)]
    >>> y1 = lsim_ss(ss_discrete, x)
    >>> y2 = lsim_ss(ss_discrete, x, 400, 1, [100, 100])

    """
    if not ssd.is_dtime():
        raise NotImplementedError("System must be discrete-time!")
    else:
        if ssd.no_inputs != 1 or ssd.no_outputs != 1:
            raise ValueError("System must be SISO!")
        if nsample is None:
            length = len(u)
        elif nsample == 0:
            raise ValueError("Sample size must be nonzero value!")
        elif nsample > len(u):
            raise ValueError("Size of the output cannot be larger than size of the input!")
        else:
            length = nsample
        result = [0.0 for _i in range(length)]
        if x0 is None:
            xold = create_empty_matrix(ssd.no_states, 1)
        else:
            xold = create_empty_matrix(len(x0), 1)
            for i in range(len(x0)):
                xold[i][0] = x0[i]
        for n in range(length):
            xnew = add_matrices(mult_matrices(ssd.A, xold), mult_matrix_scalar(ssd.B, u[n]))
            y = add_matrices(mult_matrices(ssd.C, xold), mult_matrix_scalar(ssd.D, u[n]))
            if select == 0:
                result[n] = y[0][0]
            else:
                result[n] = xold[select-1][0]
            xold = xnew
        return result


def mult_matrices(x1, x2):
    """
    Calculate the multiplication of two matrices.

    Parameters
    ----------
    x1: first matrix as a list.
    x2: second matrix as a list.

    Returns
    -------
    Multiplication of two matrices as a list.

    Raises
    ------
    ValueError
        if the column size of first matrix and row size of second matrix are not equal.

    Examples
    --------
    >>> x1 = [[1.0], [2.0]]
    >>> x2 = [[2.0, 4.0]]
    >>> y = mult_matrices(x1, x2)

    """
    result = 0
    if len(x1[0]) == len(x2):
        y = create_empty_matrix(len(x1), len(x2[0]))
        for n in range(len(x1)):
            for k in range(len(x2[0])):
                for m in range(len(x1[0])):
                    result = result + x1[n][m] * x2[m][k]
                y[n][k] = result
                result = 0
    else:
        raise ValueError("Matrices are not suitable for multiplication.")
    return y


def mult_matrix_scalar(x, A):
    """
    Calculate the multiplication of a matrix with scalar.

    Parameters
    ----------
    x: Matrix as a list.
    A: scalar value.

    Returns
    -------
    Calculated matrix as a list.

    Examples
    --------
    >>> x = [[1.0], [2.0]]
    >>> y = mult_matrix_scalar(x, 2)

    """
    y = create_empty_matrix(len(x), len(x[0]))
    for n in range(len(x)):
        for k in range(len(x[0])):
            y[n][k] = A*x[n][k]
    return y


def add_matrices(x1, x2):
    """
    Calculate the addition of two matrices.

    Parameters
    ----------
    x1: first matrix as a list.
    x2: second matrix as a list.

    Returns
    -------
    Sum of two matrices as a list.

    Raises
    ------
    ValueError
        if the row sizes or column sizes of the matrices are not equal.

    Examples
    --------
    >>> x1 = [[1.0], [2.0]]
    >>> x2 = [[2.0], [4.0]]
    >>> y = add_matrices(x1, x2)

    """
    if (len(x1) == len(x2)) and (len(x1[0]) == len(x2[0])):
        y = create_empty_matrix(len(x1), len(x1[0]))
        for n in range(len(x1)):
            for k in range(len(x1[0])):
                y[n][k] = x1[n][k] + x2[n][k]
    else:
        raise ValueError("Matrices are not the same size.")
    return y


def create_empty_matrix(row, column):
    """
    Create empty matrix as a list with given dimensions.

    Parameters
    ----------
    row: row size of the empty matrix.
    column: column size of the empty matrix.

    Returns
    -------
    Created empty matrix as a list.

    Examples
    --------
    >>> x = create_empty_matrix(3, 2)

    """
    x = [[0.0 for _n in range(column)] for _k in range(row)]
    return x
