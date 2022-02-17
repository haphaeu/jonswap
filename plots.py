"""
Lazy Python interface to make plots of JONSWAP spectrum and timeserie.
"""
import numpy
import subprocess
import matplotlib.pyplot as plt

Hs = 4.5
Tp = 10
duration = 300

# Timeserie

result = subprocess.run(
    [
        './jonswap',
        str(Hs),
        str(Tp),
        '-d', str(duration),
        '-h'
    ],
    stdout=subprocess.PIPE
)

data = numpy.loadtxt(
    result.stdout.decode('utf-8').split('\n'),
    skiprows=4
)

time, elev = data.transpose()

plt.close()
plt.plot(time, elev)
plt.grid()
plt.title(f'JONSWAP Hs {Hs} Tp {Tp}')
plt.xlabel('Time [s]')
plt.ylabel('Elevation [m]')
plt.savefig('timeserie.png')
plt.close()

# Spectrum

result = subprocess.run(
    [
        './jonswap',
        str(Hs),
        str(Tp)
    ],
    stdout=subprocess.PIPE
)

data = numpy.loadtxt(
    result.stdout.decode('utf-8').split('\n'),
    skiprows=5
)

T, w, PM, JS, amp, phi = data.transpose()

plt.plot(T, JS)
plt.grid()
plt.title(f'JONSWAP Hs {Hs} Tp {Tp}')
plt.xlabel('Period [s]')
plt.ylabel('Energy [m2s/rd]')
plt.savefig('spectrum.png')
plt.close()
