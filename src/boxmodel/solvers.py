from numpy import zeros


def evol_ratio(self, ratio, t):
    """ The evolution function that is used for isotopic ratio evolution"""
    rationew = zeros(ratio.size)
    for ii in range(ratio.size):
        outflux = 0
        influx = 0
        for jj in range(ratio.size):
            outflux = outflux + self._Flux[ii][jj] / self._Mass[ii] * \
                self._Partcoeff[ii][jj] * ratio[ii]
            influx = influx + self._Flux[jj][ii] / self._Mass[ii] * \
                self._Partcoeff[jj][ii] * ratio[jj]
        rationew[ii] = influx - outflux
    return rationew