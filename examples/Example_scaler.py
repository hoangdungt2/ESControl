'''
    Simple example for a scalar static function
    input u in range of (2500,6000)
    by jose.vu (jose.vu@kaer.com)
'''
import numpy as np 
from ESControl import ESControl as esc 
import matplotlib.pylab as plt 
from scipy.optimize import minimize_scalar

# function to be optimized
p = [  1.58951919e-05,  -1.38389530e-01,   1.18363410e+03]
# p = [2e-5, 1e-1, 1e3]
def f(x):
    return np.polyval(p,x)

# parameters for ESControl
Ts   = 5     # seconds
umin = 2500  # minimum input to f
umax = 6000  # maximum input to f
int0 = 3000  # initial input to f 

# injection sine wave
asine = (umax-umin) * .03  # amplitude: usually 3% of the range
Tsine = Ts * 30            # in seconds, this period needs to be slow enough to capture the plant dynamic

# for low/high pass filters
fLow  = (1/Tsine)/10  # cut-off frequency in Hz, should be << 1/Tsine
fHigh = fLow

# these are for kg tuning, read Bryan thesis for detail
dtime = 15 * Tsine                              # converge to minimum point after 15 cycles of injection sine
kg    = (1 - np.exp( -Ts/dtime ) )/(p[0])

# create an ESC instance
esc0 = esc(  Ts = Ts, u_min = umin, u_max = umax,  \
             int0 = int0, Asine = asine, \
             Tsine = Tsine, \
             fLow  = fLow, \
             fHigh = fHigh, \
             kg = kg )

d0    = {'u':[] ,'y':[]}     # data dictionary to store input, output of f
u0    = int0
n     = round((dtime*1.5)/Ts)  # simulation time (3 times the dtime)
d0['u'].append(u0)

#simulation time
k  = 0
while (k<n):
    J  = f( d0['u'][k] )
    u  = esc0.getControlSingal( J )
    d0['u'].append( u )
    d0['y'].append( J )
    k += 1

# use scipy to search for optimal point for plotting 
res  = minimize_scalar( f, bounds=(umin,umax), method='bounded' )
uopt = res.x
yopt = res.fun


# plotting time
t = np.arange(0,n)*Ts
d0['u'] = d0['u'][:-1]
d0['u'] = np.asarray(d0['u']).flatten()
d0['y'] = np.asarray(d0['y']).flatten()

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot( t,d0['u'],label='esc signal (input to f)' )
axarr[0].axhline(y=uopt,color='r',linestyle=':', label = 'optimal input to f')
axarr[0].grid()
axarr[0].set_title('u')
axarr[0].legend(loc='best')

axarr[1].plot( t,d0['y'],label='J signal (output of f)' )
axarr[1].axhline(y=yopt,color='r',linestyle=':', label = 'optimal output of f')
axarr[1].set_title('y')
axarr[1].grid()
axarr[1].set_xlabel('Time (s)')
axarr[1].legend(loc='best')
plt.show()

