import numpy as np 

# define class ESControl
# using class will eliminate the need of global variables

class ESControl:
    '''
        Ts      : sampling time
        u_min   : min of control signal
        u_max   : max of control signal
        int0    : initial of control signal (should be the control signal value when ESControl starts)
        Asine   : amplitude of the pertubation sine wave (usually ~5% of max-min)
        Tsine   : period of sine wave (in seconds)
        fLow    : cut-off frequency of  low-pass filter (Hz), usually << 1/Tsine
        fHigh   : cut-off frequency of high-pass filter (Hz), usually << 1/Tsine
        kg      : convergence gain (see Bryan's thesis for tuning), negative: maximum seeking, positive: minimum seeking
    ''' 
    def __init__(self, **kwargs):
        if not 'Ts' in kwargs:
            print(' ERROR: ESControl needs Ts')
            return
        if not 'u_min' in kwargs:
            print(' ERROR: ESControl needs u_min')
            return
        if not 'u_max' in kwargs:
            print(' ERROR: ESControl needs u_max')
            return
        if not 'int0' in kwargs:
            print(' ERROR: ESControl needs int0 (initial value of integrator or output)')
            return
        if not 'Asine' in kwargs:
            print(' ERROR: ESControl needs Asine (sinewave amplitue, suggestion: 5% of control range)')
            return
        if not 'Tsine' in kwargs:
            print(' ERROR: ESControl needs Tsine (sinewave period in seconds)')
            return
        if not 'fLow' in kwargs:
            print(' ERROR: ESControl needs fLow (cutoff Freq for low-pass filter (Hz))')
            return
        if not 'fHigh' in kwargs:
            print(' ERROR: ESControl needs fLow (cutoff Freq for high-pass filter (Hz))')
            return
        if not 'kg' in kwargs:
            print (' ERROR: ESControl needs kg (integrator gain, becareful with this!!!, should be positive)')
            return 
        if not 'delay_start' in kwargs:
            self.delay_start = 2 # delay start iteration to collect data
        else:
            self.delay_start = kwargs['delay_start']
        self.Asine   = kwargs['Asine']
        self.intk    = kwargs['int0']
        self.umin    = kwargs['u_min']
        self.umax    = kwargs['u_max']
        self.Ts      = kwargs['Ts']
        self.Tsine   = kwargs['Tsine']
        self.fLow    = kwargs['fLow']
        self.fHigh   = kwargs['fHigh']
        self.kg      = kwargs['kg']
        self.hpfk    = 0
        self.lpfk    = 0
        self.k       = 0
        # self.ahpf    = np.exp( -self.Ts * self.fHigh  )
        # self.alpf    = np.exp( -self.Ts * self.fLow  )
        self.ahpf    = np.exp( -2 * np.pi * self.Ts * self.fHigh  )
        self.alpf    = np.exp( -2 * np.pi * self.Ts * self.fLow  )
        
        self.omega   = 2*np.pi/self.Tsine

    def getControlSingal(self, J):
        if self.k > self.delay_start:
            self.yk = J
            self.hpfk = self.ahpf * self.hpfk + (1-self.ahpf)*J
            
            productk  = 2.0/self.Asine * np.sin(self.k * self.Ts *self.omega ) * (J-self.hpfk)
            self.lpfk = self.alpf * self.lpfk + (1-self.alpf) * productk;
            
            self.intk = self.intk - self.kg * self.lpfk * self.Ts + self.Ts * self.Asine * (self.omega) * np.cos( self.k * self.Ts *self.omega)
            self.intk = max( self.intk, self.umin )
            self.intk = min( self.intk, self.umax )
        else:
            self.hpfk = J
        self.uk   = self.intk
        self.k   += 1

        return self.uk

