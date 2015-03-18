'''
Created on Dec 10, 2012

@author: william
'''

import numpy as np


class n_component(object):  # Define a N-component SFH vector.

    def __init__(self, ages):
        """
        Calculate a parametric SFH history model.

        Parameters
        ----------
        ages : array_like
            Array with the ages of the simple stellar populations (in yr) .
        """

        self.ages = ages
        # Calculate distances between ages components
        self.DeltaT = (ages[1:] - ages[:-1]) / 2
        #Eval left and right side of each element of the base
        left = np.insert(self.DeltaT, 0, 0)
        right = np.append(self.DeltaT, 0)
        #Eval start and end of each element of the base
        self.ages_start = ages - left
        self.ages_end = ages + right
        self._individual_sfh_curves = []  # Here we will append one by one the ugly components...
        self._total_M_fraction = 0.0  # Will be used to store the total mass fraction of each component to be used afterwards on the normalization.

    def _exp(self, t_eb, tau):
        mf_ = np.zeros_like(self.ages, dtype=np.float)
        for i_age in range(len(self.ages)):
            if t_eb > self.ages_start[i_age]:
                if t_eb >= self.ages_end[i_age]:
                    mf_[i_age] = tau * (
                    np.exp((self.ages_end[i_age] - t_eb) / tau) - np.exp((self.ages_start[i_age] - t_eb) / tau) )
                elif t_eb < self.ages_end[i_age]:
                    mf_[i_age] = tau * ( 1 - np.exp((self.ages_start[i_age] - t_eb) / tau) )

        Me = np.sum(mf_)
        return mf_ / Me

    def _gauss(self, t_b, tau):
        mf_ = np.zeros_like(self.ages, dtype=np.float)
        for i_age in range(len(self.ages)):
            if t_b > self.ages_start[i_age]:
                if t_b >= self.ages_end[i_age]:
                    mf_[i_age] = np.exp(-(t_b - self.ages_end[i_age]) ** 2 / (2 * tau ** 2)) - np.exp(
                        -(t_b - self.ages_start[i_age]) ** 2 / (2 * tau ** 2))
                elif t_b < self.ages_end[i_age]:
                    mf_[i_age] = 1 - np.exp(-(t_b - self.ages_start[i_age]) ** 2 / (2 * tau ** 2))

        Me = np.sum(mf_)
        return mf_ / Me

    def _square(self, t0, l):
        # Set an array with all the elements between t0 and t0+l to 1
        mf_ = np.zeros_like(self.ages, dtype=np.float)
        mf_[np.bitwise_and(self.ages_end >= t0, self.ages_start <= (t0 + l))] = 1
        args = np.argwhere(mf_ == 1)
        #Deal with the boundaries
        #lower
        mf_[args[0]] = (self.ages_end[mf_ == 1][0] - t0) / (self.ages_end[mf_ == 1][0] - self.ages_start[mf_ == 1][0])
        #upper
        if len(args) > 1:
            mf_[args[-1]] = (t0 + l - self.ages_start[mf_ == 1][-1]) / (
            self.ages_end[mf_ == 1][-1] - self.ages_start[mf_ == 1][-1])

        mf_ = mf_ * (self.ages_end - self.ages_start)
        # Total Mass
        Ms = np.sum(mf_)
        return mf_ / Ms

    def add_exp(self, t_eb, tau, frac):
        """
        Add a exponential :math:`A e^{- (t_0 - t) / \\tau}` to the SFH curves.

        Parameters
        ----------
        t_eb : float
            :math:`t_0` of the exponential burst.

        tau : float
            :math:`\\tau` of the exponential burst.

        frac : float
            Light fraction :math:`A` of this burst.
        """
        self._individual_sfh_curves.append(self._exp(t_eb, tau) * frac)
        self._total_M_fraction += frac

    def add_gauss(self, t_b, tau, frac):
        """
        Add a gaussian :math:`t e^{-(t_0 - t)^2/2 \\tau^2}` to the SFH curves.

        Parameters
        ----------
        t_b : float
            Time of the center of the burst.

        tau : float
           Related to the burst width. See eq above.

        frac : float
            Light fraction of this burst.

        """
        self._individual_sfh_curves.append(self._gauss(t_b, tau) * frac)
        self._total_M_fraction += frac

    def add_square(self, t0, l, frac, t_max=None):
        """
        Add a square to the SFH curves.

        Parameters
        ----------
        t0 : float
            :math:`t_0` of the burst.

        l : float
            Burst duration

        frac : float
            Light fraction of this burst.

        t_max : float
            Limit the bust to no older than this parameter.
        """
        if t_max:
            square = self._square(t0, l)
            square[self.ages > t_max] = 0.0
            self._individual_sfh_curves.append(square / np.trapz(square, self.ages) * frac)
        else:
            self._individual_sfh_curves.append(self._square(t0, l) * frac)
        self._total_M_fraction += frac

    @property
    def individual_sfh_curves(self):
        """
        Returns a list with the individual SFH curves.
        """
        return self._individual_sfh_curves

    def get_sfh(self):
        """
        Add a square to the SFH curves.

        Returns
        -------
        sfh : array-like
            Array with the final SFH curve. The curve is normalized to 1.

        See Also
        --------
        individual_sfh_curves

        """
        sfh = np.sum(self.individual_sfh_curves, axis=0)
        sfh /= np.sum(sfh)

        return sfh

    def plot(self, log=True):
        """
        Do a plot of the SFH curves.

        Parameters
        ----------
        log : bool
            If True, will plot the x-axis in log.
        """
        import matplotlib.pyplot as plt

        plt.figure(1)
        plt.clf()
        if log:
            x = np.log10(self.ages)
        else:
            x = self.ages
        for i_sfh in range(len(self.individual_sfh_curves)):
            plt.plot(x, self.individual_sfh_curves[i_sfh], label='Curve %s' % np.int(i_sfh))
        plt.plot(x, self.get_sfh(), lw=2, alpha=0.8, label='Total SFH curve')
        plt.legend(loc=2)
        if log:
            plt.xlabel('$\log t$')
        else:
            plt.xlabel('$t$')
        plt.ylabel('SFH')
