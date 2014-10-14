
from astropy.cosmology import WMAP9 as cosmo
from astropy import units
from matplotlib import cm

__author__ = 'william'

import io
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter

from ..util.matchs import get_zslice
from ..core.exceptions import MAGALException

import ConfigParser
import ast

class plotMagal():

    def __init__(self, Magal, Input, Library=None):
        self.i = Input
        self.m = Magal
        self.l = Library
        self._check_paths()
        self._check_magal_config()

    def _check_paths(self):
        if self.i.path != self.m.path:
            raise 'Error: Library and MagalFit paths must be the same!'

    def plot_inout_params(self, param, type='hexbin', stat='AVG', pNaN_treshold=None, cmap=cm.hot_r,
                          plot_errorbars=False, plot_errorlines=False, xscale='linear', yscale='linear', wait=False,
                          autosave=False, autosave_dir='./', xlim=None, ccds=[1], plot_errorbars_count=False, grid=False):

        if param is None:
            for param in self.m.props:
                self.plot_inout_param(str(param), type, stat, pNaN_treshold, cmap, plot_errorbars, plot_errorlines,
                                      xscale, yscale, wait, autosave, autosave_dir, xlim, ccds, plot_errorbars_count, grid)
        elif isinstance(param, basestring):
            self.plot_inout_param(str(param), type, stat, pNaN_treshold, cmap, plot_errorbars, plot_errorlines,
                                      xscale, yscale, wait, autosave, autosave_dir, xlim, ccds, plot_errorbars_count, grid)
        elif not isinstance(param, (tuple, list)):
            raise 'param must be string, a tuple or a list or None for all props of magal model.'

        for p_ in param:
            self.plot_inout_param(p_, type, stat, pNaN_treshold, cmap, plot_errorbars, plot_errorlines, xscale,
                                  yscale, wait, autosave, autosave_dir, xlim, ccds, plot_errorbars_count, grid)


    def _check_magal_config(self):
        self.model_config = ConfigParser.RawConfigParser(allow_no_value=True)
        self.model_config.readfp(io.BytesIO(self.m.ini_file))
        try:
            self.Nobj_max = int(self.model_config.get('FitGeneral', 'Nobj_max'))
        except ConfigParser.NoOptionError:
            self.Nobj_max = None # None means run over all objects!
        try:
            self.filters_exclude = ast.literal_eval(self.model_config.get('FitGeneral', 'filters_exclude'))
        except ConfigParser.NoOptionError:
            self.filters_exclude = None

    def plot_fit(self, i_obj, ax=None, is_simulation=False, starlight_input=None, percentile=16):

        if ax is None:
            fig = plt.figure(1)
            fig.clf()
            ax = fig.add_axes((0.1,0.1,.8,.8))

        if self.l is None:
            raise MAGALException('You must init the object w/ a Library!')

        bmx_z, bmx_T = self.m.i_BMX[i_obj]
        if bmx_z == 0:  #FIXME:
            i_z = np.argmin((self.l.z - self.i.z[0])**2)
        else:
            i_z = bmx_z

        aux_model = self.l.library[i_z, bmx_T]['m_ab'] + self.m.s[i_obj, bmx_z, bmx_T]

        if self.filters_exclude:
            filterset_mask = ~np.sum([k == self.l.filterset['ID_filter'] for k in self.filters_exclude], axis=0, dtype=np.bool)
        else:
            filterset_mask = np.ones(len(self.l.filterset), dtype=np.bool)

        if is_simulation: #FIXME
            aux_inp = self.i.library[0, i_obj]['m_ab'] - 2.5 * self.i.properties[i_obj]['Mini_fib'] #FIXME
            aux_err = self.i.library[0, i_obj]['e_ab'] #FIXME
        else:
            aux_inp = self.i.data.value[i_obj]['m_ab']
            aux_err = self.i.data.value[i_obj]['e_ab']

        filterset_mask = np.bitwise_and(filterset_mask, np.isfinite(aux_model))
        filterset_mask = np.bitwise_and(filterset_mask, np.isfinite(aux_inp))

        aux_lib = np.ma.masked_invalid(self.l.library['m_ab'])

        aux_model_tot = [
            ((aux_lib[i_z,:,i_mag] + self.m.s[i_obj, i_z]) * self.m.likelihood[i_obj, i_z]).sum() / self.m.likelihood[i_obj, i_z][np.argwhere(~aux_lib[i_z,:,i_mag].mask)].sum()
            for i_mag in np.argwhere(filterset_mask)]

        aux_l = self.l.filterset['wl_central'][filterset_mask]

        z_lib = self.l.z[i_z]
        if z_lib == 0:
            d_L_lib = 3.08567758e19  # 10 parsec in cm
        else:
            d_L_lib = cosmo.luminosity_distance(z_lib).to('cm')
        # k_cosmo = 1/(4 * np.pi * np.power(d_L, 2)) #FIXME:

        # aux_l2 = aux_l**-2
        # ax.plot(aux_l, aux_l2 * 10**(-.4 * (aux_inp + 2.41)) * k_cosmo.value, '.', color='blue')
        # ax.plot(aux_l, aux_l2 * 10**(-.4 * (aux_model[filterset_mask] + 2.41)) * k_cosmo.value, '.', color='red')

        # ax.plot(aux_l, aux_inp+aux_err, color='black')
        # ax.plot(aux_l, aux_inp-aux_err, color='black')

        plt.plot(aux_l, aux_model_tot)
        n_percentile = self.m.get('/Alhambra_24/1/statistics/model/n_percentiles')[i_obj][str(percentile)]
        # print np.argsort(self.m.likelihood[0])[0][::-n_percentile]
        print n_percentile
        # for i_T in np.argsort(self.m.likelihood[0])[0][::-1][1:n_percentile]:
        #     aux_model = self.l.library[i_z, i_T]['m_ab'] + self.m.s[i_obj, bmx_z, i_T]
        #     plt.plot(aux_l, aux_model[filterset_mask], '.', color='black', linewidth=1, alpha=.1)
        #     # print 'i', i


        if starlight_input:
            z_sl = self.i.properties['z'][i_obj]
            if z_sl == 0:
                d_L = 3.08567758e19  # 10 parsec in cm
            else:
                d_L = cosmo.luminosity_distance(z_sl).to('cm')
            k_cosmo = (1+z_sl)**3 / (4 * np.pi * np.power(d_L.value, 2)) #FIXME:
            k_cosmo *= (4 * np.pi * np.power(d_L_lib.value, 2))

            try:
                print 'Trying to load', starlight_input
                aux_slight = np.loadtxt(starlight_input, dtype=np.dtype([('wl', '<f8'), ('flux', '<f8'), ('error', '<f8'), ('flag', '<i8')])) #, )
                wl_ = aux_slight['wl'] * (1+z_lib)
                fl_ = aux_slight['flux'] * np.power(z_sl + 1., 2) / (1+z_lib)
                # plt.plot(wl_, -2.5 * np.log10(aux_slight['flux'] * 1e-17 * wl_**2 / k_cosmo ) - 2.41, alpha=.3)
            except IOError:
                pass

        zp_error = self.m.input_config.getfloat('FitGeneral', 'zp_error') #FIXME
        plt.errorbar(aux_l, aux_inp[filterset_mask], fmt='o', color='green', linewidth=2, yerr=np.sqrt(aux_err[filterset_mask]**2+zp_error**2))
        plt.plot(aux_l, aux_model[filterset_mask], 'o', color='red', linewidth=2)



            # plt.plot(aux_l)
                #plt.plot(aux_l, 10**-.4*(aux_inp + 2.41) * (aux_l**-2))

        # print aux_l, aux_inp, aux_model


    def plot_inout_param(self, param, plt_type, stat, pNaN_treshold, cmap, plot_errorbars, plot_errorlines, xscale, yscale,
                         wait, autosave, autosave_dir, xlim, ccds, plot_errorbars_count, grid):
        """
        bla

        Parameters
        ----------

        Returns
        -------

        Examples
        --------

        See Also
        --------

        Notes
        -----
        """

        plt.clf()

        if isinstance(param, basestring):
            p1 = p2 = save_name = param
            p_color = None
        elif len(param) == 2:
            p1, p2 = param
            p_color = None
            save_name = '_'.join(param)  # Will be used in the end, to savefig()
        elif len(param) == 3:
            p1, p2, p_color = param
            if type(p_color).__name__ == 'ndarray': param = [p1, p2, 'custom']
            save_name = '_'.join(param[0:2])  # Will be used in the end, to savefig()

        print param

        print 'Plotting properties %s and %s' % (p1, p2)
        if ccds is None:
            ccds = self.i.ccds[self.i.fsys]
        for ccd in ccds:
            print '%s/%s' % (self.i.fsys, ccd)
            self.i.get_filtersys(self.i.fsys, ccd)
            self.m.get_filtersys(self.i.fsys, ccd)
            try:
                l_ = len(input_param)
                input_param = np.ma.masked_invalid(np.append(input_param, np.ma.masked_invalid(np.copy(self.i.properties[p1][:self.Nobj_max]))))
                if isinstance(p_color, basestring):
                    z = np.ma.masked_invalid(np.append(z, np.ma.masked_invalid(np.copy(self.i.properties[p_color][:self.Nobj_max]))))
                elif p_color:
                    z = np.ma.masked_invalid(np.append(z, np.ma.masked_invalid(p_color[:self.Nobj_max])))
                if xscale == 'log':
                    if stat == 'AVG':
                        dp = np.append(dp, np.power(10,np.copy(self.m.__getattribute__(p2)[stat])) - np.power(10,input_param[l_:]))
                    else:
                        dp = np.append(dp, np.power(10,np.copy(self.m.__getattribute__(p2)['percentiles'][stat])) - np.power(10,input_param[l_:]))
                else:
                    if stat == 'AVG':
                        dp = np.append(dp, np.copy(self.m.__getattribute__(p2)[stat]) - input_param[l_:])
                    else:
                        dp = np.append(dp, np.copy(self.m.__getattribute__(p2)['percentiles'][stat]) - input_param[l_:])
                    dp = np.ma.masked_array(dp, mask = input_param.mask)
            except UnboundLocalError:
                input_param = np.ma.masked_invalid(np.copy(self.i.properties[p1][:self.Nobj_max]))
                if isinstance(p_color, basestring):
                    z = np.ma.masked_invalid(np.copy(self.i.properties[p_color][:self.Nobj_max]))
                elif np.any(p_color):
                    z = np.ma.masked_invalid(np.copy(p_color[:self.Nobj_max]))
                    p_color = True
                if xscale == 'log':
                    if stat == 'AVG':
                        dp = np.ma.masked_array(np.copy(np.power(10,self.m.__getattribute__(p2)[stat]))) - np.power(10,input_param)
                    else:
                        dp = np.ma.masked_array(np.copy(np.power(10,self.m.__getattribute__(p2)['percentiles'][stat]))) - np.power(10,input_param)
                else:
                    if stat == 'AVG':
                        dp = np.ma.masked_array(np.copy(self.m.__getattribute__(p2)[stat])) - input_param
                    else:
                        dp = np.ma.masked_array(np.copy(self.m.__getattribute__(p2)['percentiles'][stat])) - input_param

        if pNaN_treshold:
            mask = np.ma.greater_equal(self.m.__getattribute__(p2)['pNaN'], pNaN_treshold)
            print '%i less objects due to pNaN_threshold' % mask.sum()
            input_param.mask[mask] = True
            dp.mask[mask] = True
            z.mask[mask] = True

        y = np.ma.copy(dp)
        x = np.ma.copy(input_param)
        del input_param
        del dp

        if xscale == 'log':
            x = 10**x

        if xlim:
            mask = np.bitwise_or(x <= xlim[0], x >= xlim[1])
            x.mask[mask] = True
            y.mask[mask] = True
            z.mask = x.mask

        x = x.compressed()
        y = y.compressed()
        # if p_color:
        z = z.compressed()

        if len(x) <= 20:
            print 'Too few points to plot! Skipping...'
            return

        plt.clf()
        #
        print x.shape,y.shape,z.shape
        if plt_type == 'hexbin':
            if xscale == 'log':
                if p_color:
                    plt.hexbin(x, y, C=z, cmap=cmap, xscale='log')
                    plt.colorbar(label=param[-1])
                else:
                    plt.hexbin(x, y, bins='log', cmap=cmap, xscale='log')
            else:
                if p_color:
                    plt.hexbin(x, y, C=z, cmap=cmap, xscale='log')
                    plt.colorbar(label=param[-1])
                else:
                    plt.hexbin(x, y, bins='log', cmap=cmap)
        else:
            if xscale == 'log':
                ax = plt.gca()
                ax.set_xscale('log')
                ax.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
            if p_color:
                plt.scatter(x,y,c=z, cmap=cmap)
                plt.colorbar(label=param[-1], cmap=cmap)
            else:
                plt.scatter(x, y)
        xlim, ylim = plt.xlim(), plt.ylim()

        if pNaN_treshold:
            plt.title('%s Mean = %3.2f, Stdev = %3.2f, N_obj = %i, stat = %s, pNaN <= %3.2f percent' % (param, np.mean(y), np.std(y), len(x), stat, 100*pNaN_treshold))
        else:
            plt.title('%s Mean = %3.2f, Stdev = %3.2f, N_obj = %i, stat = %s' % (param, np.mean(y), np.std(y), len(x), stat))
        plt.xlabel('Input')
        plt.ylabel('Output - Input')

        ylim = np.sqrt(np.max(np.power(ylim,2)))

        plt.ylim((-ylim, ylim))

        if plot_errorbars:
            n_bins = 10
            extent = [np.min(x), np.max(x), np.min(y), np.max(y)]
            if xscale == 'linear':
                bins_cut = np.linspace(extent[0], extent[1], n_bins)
            else:
                bins_cut = np.logspace(np.log10(extent[0]), np.log10(extent[1]), 10)
            aux_x = []; aux_y = []; aux_error = []; aux_count = [];
            for i_bin in range(len(bins_cut) - 1):
                mask = np.bitwise_and(x >= bins_cut[i_bin], x < bins_cut[i_bin+1])
                if mask.sum() > 0:
                    avg, std = np.average(y[mask]), np.std(y[mask])
                    aux_x.append(bins_cut[i_bin] + (bins_cut[i_bin+1] - bins_cut[i_bin])/2)
                    aux_y.append(avg)
                    aux_error.append(std)
                    aux_count.append(mask.sum())
                    xlim, ylim = plt.xlim(), plt.ylim()
            plt.errorbar(aux_x, aux_y,yerr=aux_error, fmt='o', color='green')
            plt.xlim(extent[0], extent[1])
            plt.ylim(extent[2], extent[3])
            if plot_errorbars_count:
                for i in range(len(aux_x)):
                    plt.text(aux_x[i], aux_y[i]+aux_error[i], aux_count[i], color='red')

            # Dotted guide lines when log xscale plot. TODO: Move this to def as an array
            if plot_errorlines:
                aux_error_lines = [.05, .1]
                aux = np.linspace(xlim[0], xlim[1])
                for el_ in aux_error_lines:
                    plt.plot(aux, aux*el_, 'b--')
                    plt.plot(aux, aux*-el_, 'b--')
                    # plt.text(xlim[1]*.5, xlim[1]*.5*el_, '%f ' % el_*100, color='red')

        ylim = np.sqrt(np.max(np.array(plt.ylim())**2))
        plt.ylim(-ylim, ylim)

        if grid:
            plt.grid()
            if xscale == 'log':
                plt.grid(which='minor')

        if autosave:
            plt.savefig('%s/%s.png' % (autosave_dir, save_name))

        if wait:
            raw_input('Press ENTER for next plot...')

    def plot_pdf(self, i_obj, prop, bins=50):

        aux = np.cumsum(np.sort(self.m.likelihood), axis=-1)
        perc = [2.5, 16, 50, 84, 97.5]
        # aux_n_perc = np.array([np.greater_equal(aux, p_/100.).sum(axis=-1) for p_ in perc])
        # print aux_n_perc[:, i_obj]

        # Plot the Prior
        plt.figure(1)
        plt.clf()
        aux_hst, aux_bins = np.histogram(self.l.properties[prop], bins=bins, normed=True)
        aux_hst /= np.max(aux_hst)

        left = np.array(aux_bins[:-1])
        right = np.array(aux_bins[1:])
        pts = left + (right - left) / 2
        plt.plot(pts, aux_hst, linewidth=2, color='blue')

        # Plot the Posterior
        aux_hst, aux_bins = np.histogram(self.l.properties[prop], weights=self.m.likelihood_T[i_obj], bins=bins)
        aux_hst /= np.max(aux_hst)

        left = np.array(aux_bins[:-1])
        right = np.array(aux_bins[1:])
        pts = left + (right - left) / 2
        plt.plot(pts, aux_hst, linewidth=2, color='green')

        aux_cumhst = np.cumsum(aux_hst)
        aux_cumhst /= aux_cumhst.max()
        plt.plot(pts, aux_cumhst, linewidth=2, color='black', alpha=0.3)

        aux_percentiles = [float(p_) for p_ in self.m.A_V.percentiles.dtype.names]
        aux_x = list(self.m.__getattribute__(prop).percentiles[i_obj])
        aux_y = np.interp(aux_x, pts, aux_cumhst)
        plt.plot(aux_x, aux_y, 'o', color='red')

        # Plot the AVG value
        ylim = (0., 1.01)
        plt.plot([self.m.__getattribute__(prop)['AVG'][i_obj], self.m.__getattribute__(prop)['AVG'][i_obj]], ylim, '--', color='red', linewidth=2, alpha=.5)

        # AND the correct value...
        plt.plot([self.i.properties[i_obj][prop], self.i.properties[i_obj][prop]], ylim, '--', color='green', linewidth=2, alpha=.5)
        plt.ylim(ylim)

        # plt.figure(2)
        # plt.clf()
        # plt.plot(self.l.properties[prop], self.m.likelihood_T[i_obj], '.')

    # def plot_bmx(self, i_gal, i_fig=1):
    #
    #     plt.figure(i_fig)
    #     plt.clf()
    #     plt.plot(self.i.filterset['wl_central'], self.i.library[0,i_gal]['m_ab'], '.', color='blue', label='Data') #TODO: do this for all kinds of input objects
    #     plt.plot(self.i.filterset['wl_central'], self.l.library[self.m.i_BMX[i_gal][0], m.i_BMX[i_gal][1]]['m_ab'], '.', color='red', label='Best match')

