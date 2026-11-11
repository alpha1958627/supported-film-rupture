"""Supported liquid film with no-slip wall and frozen flat-film mobility."""

import numpy as np


class WallFilm2D:
    def __init__(self, n, length, dt, attraction=0., theta=0.):
        self.n, self.length, self.dt = n, length, dt
        self.attraction, self.theta = attraction, theta
        dx = length/n
        k = 2*np.pi*np.fft.fftfreq(n, dx)
        self.q = np.hypot(k[:, None], k)
        qd = (2/dx)*np.sin(k*dx/2)
        qd2 = qd[:, None]**2 + qd**2
        mobility = qd2/3.
        stiffness = qd2
        self.mobility = mobility
        self.rate = mobility*(3*attraction - stiffness)
        half = .5*dt*self.rate
        if np.any(half >= 1.):
            raise ValueError('nonpositive midpoint denominator')
        self.gain = (1+half)/(1-half)
        self.pressure = -dt*mobility/(1-half)
        self.noise = np.sqrt(2*theta*mobility*dt)/(dx*(1-half))

    @staticmethod
    def spectral(field, symbol):
        return np.fft.ifft2(symbol*np.fft.fft2(field)).real

    def advance(self, height, white=None):
        old = height - 1.
        base = self.spectral(old, self.gain)
        if white is not None:
            base += self.spectral(white, self.noise)
        a = self.attraction
        if a == 0.:
            return 1.+base

        updated = base.copy()
        for _ in range(8):
            u = .5*(old+updated)
            h = 1.+u
            pressure = a*u*u*(6+u*(8+3*u))/h**3
            residual = updated-base-self.spectral(pressure, self.pressure)
            if np.max(abs(residual)) <= 2e-13:
                return 1.+updated
            derivative = 3*a*(1-h**-4)
            correction = residual.copy()
            for _ in range(12):
                error = residual-correction+self.spectral(.5*derivative*correction, self.pressure)
                if np.max(abs(error)) <= 2e-14:
                    break
                correction += error
            else:
                raise RuntimeError('JAM midpoint Jacobian solve did not converge')
            updated -= correction
        raise RuntimeError('JAM midpoint residual did not converge')
