import numpy as np
import scipy.ndimage.filters as sc
from astropy import constants as c
from pystarlight.util.velocity_fix import SpectraVelocityFixer

lick_dt = np.dtype([('index_low', np.float), ('index_upp', np.float), ('blue_low', np.float), ('blue_upp', np.float),
                    ('red_low', np.float), ('red_upp', np.float)])

lick_indexes = {'CN_1': np.array((4142.125, 4177.125, 4080.125, 4117.625, 4244.125, 4284.125), dtype=lick_dt),
                'CN_2': np.array((4142.125, 4177.125, 4083.875, 4096.375, 4244.125, 4284.125), dtype=lick_dt),
                'Ca4227': np.array((4222.250, 4234.750, 4211.000, 4219.750, 4241.000, 4251.000), dtype=lick_dt),
                'G4300': np.array((4281.375, 4316.375, 4266.375, 4282.625, 4318.875, 4335.125), dtype=lick_dt),
                'Fe4383': np.array((4369.125, 4420.375, 4359.125, 4370.375, 4442.875, 4455.375), dtype=lick_dt),
                'Ca4455': np.array((4452.125, 4474.625, 4445.875, 4454.625, 4477.125, 4492.125), dtype=lick_dt),
                'Fe4531': np.array((4514.250, 4559.250, 4504.250, 4514.250, 4560.500, 4579.250), dtype=lick_dt),
                'Fe4668': np.array((4633.999, 4720.250, 4611.500, 4630.250, 4742.750, 4756.500), dtype=lick_dt),
                'H_beta': np.array((4847.875, 4876.625, 4827.875, 4847.875, 4876.625, 4891.625), dtype=lick_dt),
                'Fe5015': np.array((4977.750, 5054.000, 4946.500, 4977.750, 5054.000, 5065.250), dtype=lick_dt),
                'Mg_1': np.array((5069.125, 5134.125, 4895.125, 4957.625, 5301.125, 5366.125), dtype=lick_dt),
                'Mg_2': np.array((5154.125, 5196.625, 4895.125, 4957.625, 5301.125, 5366.125), dtype=lick_dt),
                'Mg_b': np.array((5160.125, 5192.625, 5142.625, 5161.375, 5191.375, 5206.375), dtype=lick_dt),
                'Fe5270': np.array((5245.650, 5285.650, 5233.150, 5248.150, 5285.650, 5318.150), dtype=lick_dt),
                'Fe5335': np.array((5312.125, 5352.125, 5304.625, 5315.875, 5353.375, 5363.375), dtype=lick_dt),
                'Fe5406': np.array((5387.500, 5415.000, 5376.250, 5387.500, 5415.000, 5425.000), dtype=lick_dt),
                'Fe5709': np.array((5696.625, 5720.375, 5672.875, 5696.625, 5722.875, 5736.625), dtype=lick_dt),
                'Fe5782': np.array((5776.625, 5796.625, 5765.375, 5775.375, 5797.875, 5811.625), dtype=lick_dt),
                'Na_D': np.array((5876.875, 5909.375, 5860.625, 5875.625, 5922.125, 5948.125), dtype=lick_dt),
                'TiO_1': np.array((5936.625, 5994.125, 5816.625, 5849.125, 6038.625, 6103.625), dtype=lick_dt),
                'TiO_2': np.array((6189.625, 6272.125, 6066.625, 6141.625, 6372.625, 6415.125), dtype=lick_dt),
                'Hdelta_A': np.array((4083.500, 4122.250, 4041.600, 4079.750, 4128.500, 4161.000), dtype=lick_dt),
                'Hgamma_A': np.array((4319.750, 4363.500, 4283.500, 4319.750, 4367.250, 4419.750), dtype=lick_dt),
                'Hdelta_F': np.array((4090.999, 4112.250, 4057.250, 4088.500, 4114.750, 4137.250), dtype=lick_dt),
                'Hgamma_F': np.array((4331.250, 4352.250, 4283.500, 4319.750, 4354.750, 4384.750), dtype=lick_dt),
                'D4000': np.array((np.nan, np.nan, 3850., 3950., 4000., 4100.), dtype=lick_dt),
                'MgFe': np.array((np.nan, np.nan, np.nan, np.nan, np.nan, np.nan), dtype=lick_dt)}


class LickIndex(object):
    def __init__(self):
        # Define some constants:
        self._light_spped = c.c.to('km/s').value  # Speed of light.
        self.lambda_vazdekis = 5000.  # Ref: Vazdekis.et.al.2010 sec 4.1.3

    @staticmethod
    def _calc_index(indexes, spec_error, spec_flags, spec_flux, spec_lambda):

        aux_return = {}
        for index_name in indexes:

            # Define where we will cut the spectra
            blue_range = np.bitwise_and(spec_lambda >= lick_indexes[index_name]['blue_low'],
                                        spec_lambda <= lick_indexes[index_name]['blue_upp'])
            red_range = np.bitwise_and(spec_lambda >= lick_indexes[index_name]['red_low'],
                                       spec_lambda <= lick_indexes[index_name]['red_upp'])
            window_range = np.bitwise_and(spec_lambda >= lick_indexes[index_name]['index_low'],
                                          spec_lambda <= lick_indexes[index_name]['index_upp'])
            # Define where spectra is good
            if spec_flags is not None:
                blue_good = spec_flags[blue_range] <= 1
                red_good = spec_flags[red_range] <= 1
                window_good = spec_flags[window_range] <= 1
            else:
                blue_good = np.ones(blue_range.sum(), dtype=np.bool)
                red_good = np.ones(red_range.sum(), dtype=np.bool)
                window_good = np.ones(window_range.sum(), dtype=np.bool)
            # ... and count how many lambdas do we have
            blue_ok = float(blue_good.sum())
            red_ok = float(red_good.sum())
            window_ok = float(window_good.sum())

            if blue_ok / len(blue_good) > 0.4 and red_ok / len(red_good) > 0.4 and window_ok / len(window_good) > 0.4:
                mean_blue = np.mean(spec_flux[blue_range])
                mean_red = np.mean(spec_flux[red_range])

                lambda_blue = (lick_indexes[index_name]['blue_low'] + lick_indexes[index_name]['blue_upp']) / 2
                lambda_red = (lick_indexes[index_name]['red_low'] + lick_indexes[index_name]['red_upp']) / 2

                # y = ax + b
                continuum = ((mean_blue - mean_red) / (lambda_blue - lambda_red)) * spec_lambda[window_range] + \
                            (mean_red * lambda_blue - mean_blue * lambda_red) / (lambda_blue - lambda_red)
                ind = np.sum((continuum - spec_flux[window_range]) / continuum)
                if spec_error is not None:
                    err_ind = np.sqrt(np.sum(spec_error[window_range] ** 2 / continuum ** 2))
                else:
                    err_ind = np.nan
                aux_return.update({index_name: (ind, err_ind)})
            else:
                aux_return.update({index_name: (np.nan, np.nan)})

        return aux_return

    def lick_index(self, indexes, spec_lambda, spec_flux, spec_error, spec_flags=None, v0=0, instrument_res=0, vd=0,
                   instrument_res_final=6, vd_final=357):
        """
        Calculate lick indexes on the input spectrum.

        Parameters
        ----------
        indexes : tuple
            Index names that will be evaluated.

        spec_lambda : array-like
            Array with wavelenghts of the spectrum in Angstrom.

        spec_flux : array-like
            Flux of the spectrum. The units will correspond to the output units of the indices.

        spec_error : array-like. Optional.
            Flux errors associated to each lambda. Same units of spec_flux.

        spec_flags : array-like
            Flags to mask the spectrum. Points with spec_flags != 0 will be masked.

        v0 : float. Default: 0.
            TODO: docme

        instrument_res : float. Default: 0.
            TODO: docme

        vd : float. Default: 0.
            TODO: docme

        instrument_res_final : float. Default: 6.
            TODO: docme

        vd_final : float. Default: 357.
            TODO: docme

        Returns
        -------
        index : dict
            Dictionary with the indexes as keywords and value[0] is the index and value[1] is its associated error.
            Missing values are represented as np.nan .
        """

        if not np.allclose(spec_lambda, np.linspace(np.min(spec_lambda), np.max(spec_lambda), len(spec_lambda))):
            # Lambdas must be equally spaced to the vd correction
            spec_lambda_new = np.arange(np.min(spec_lambda), np.max(spec_lambda), 1.0)
            spec_flux = np.interp(spec_lambda_new, spec_lambda, spec_flux)
            if spec_flags is not None:
                spec_flags = np.interp(spec_lambda_new, spec_lambda, spec_flags)
            if spec_error is not None:
                spec_error = np.interp(spec_lambda_new, spec_lambda, spec_error)
            spec_lambda = spec_lambda_new

        aux_return = {}

        # ahora degradamos la resolucion. Primero en lambda por si queremos pasar por ejemplo 3A (FWHM) FIXME:
        # de los modelos a 6A (FWHM) de la resolucion instrumental de califa    #FIXME:
        sigma_instrument = instrument_res / 2.3548200450309493  # 2 * np.sqrt(2 * np.log(2))
        sigma_l_fin = instrument_res_final / 2.3548200450309493
        sigma_dif = np.sqrt(sigma_l_fin ** 2 - sigma_instrument ** 2)
        spec_flux = sc.gaussian_filter1d(spec_flux, sigma_dif)

        sigma_instrument = instrument_res_final / 2.3548200450309493
        sigma_d = (vd / self._light_spped) * self.lambda_vazdekis
        sigma_ini = np.sqrt(sigma_instrument ** 2 + sigma_d ** 2)
        vel_ini = (sigma_ini / self.lambda_vazdekis) * self._light_spped

        # #####bloque para corregir de velocidad de dispersion (funcion de Andre). Dando v_0 a la funcion tambien FIXME:
        # corrige de velocidad de rotacion.  # FIXME:
        if vel_ini < vd_final:
            vfix = SpectraVelocityFixer(spec_lambda, v0, vel_ini, nproc=1)
            spec_flux = vfix.fix(spec_flux, vd_final)

        for index_name in indexes:
            if index_name in ('D4000', 'MgFe'):
                if index_name == 'D4000':
                    blue_range = np.bitwise_and(spec_lambda >= lick_indexes[index_name]['blue_low'],
                                                spec_lambda <= lick_indexes[index_name]['blue_upp'])
                    red_range = np.bitwise_and(spec_lambda >= lick_indexes[index_name]['red_low'],
                                               spec_lambda <= lick_indexes[index_name]['red_upp'])

                    if spec_flags is not None:
                        # Define where spectra is good
                        blue_good = spec_flags[blue_range] <= 1
                        red_good = spec_flags[red_range] <= 1
                        # ... and count how many lambdas do we have
                    else:
                        blue_good = np.ones(blue_range.sum(), dtype=np.bool)
                        red_good = np.ones(red_range.sum(), dtype=np.bool)

                    blue_ok = float(blue_good.sum())
                    red_ok = float(red_good.sum())

                    if blue_ok / len(blue_good) > 0.4 and red_ok / len(red_good) > 0.4:
                        flux_blue = np.sum(spec_flux[blue_range])
                        flux_red = np.sum(spec_flux[red_range])
                        d4000 = flux_red / flux_blue
                        if spec_error is not None:
                            err_d4000 = d4000 * np.sqrt(np.sum(spec_error[blue_range] ** 2) / (flux_blue ** 2)
                                                        + np.sum(spec_error[red_range] ** 2) / (flux_red ** 2))
                        else:
                            err_d4000 = np.nan
                        aux_return.update({'D4000': (d4000, err_d4000)})
                elif index_name == 'MgFe':
                    aux = self._calc_index(('Mg_b', 'Fe5270', 'Fe5335'), spec_error, spec_flags, spec_flux, spec_lambda)
                    aux_fe = 0.72 * aux['Fe5270'][0] + 0.28 * aux['Fe5335'][0]
                    MgFe = np.sqrt(aux['Mg_b'][0] * aux_fe)
                    if spec_error is not None:
                        MgFe_err = (aux['Mg_b'][1] ** 2 * aux_fe ** 2 +
                                    aux['Fe5270'][1] ** 2 * (0.72 * aux['Mg_b'][0]) ** 2 +
                                    aux['Fe5335'][1] ** 2 * (0.28 * aux['Mg_b'][0]) ** 2) * .25 * MgFe ** -2
                        MgFe_err = np.sqrt(MgFe_err)
                    else:
                        MgFe_err = np.nan
                    aux_return.update({'MgFe': (MgFe, MgFe_err)})
            else:
                aux_return.update(self._calc_index((index_name,), spec_error, spec_flags, spec_flux, spec_lambda))

        return aux_return

    def test(self):
        import os
        import magal
        import pprint

        data_example_dir = '%s/../test/data_example/' % os.path.dirname(magal.__file__)
        dt = np.dtype([('lambda', np.float), ('flux', np.float), ('error', np.float), ('flag', np.int)])
        for specfile in ('0266.51630.004.7xt.bz2', '0266.51630.014.7xt.bz2', '0266.51630.022.7xt.bz2'):
            spec = np.loadtxt('%s/starlight/input/%s' % (data_example_dir, specfile), dtype=dt)
            print specfile + ':'
            pprint.pprint(self.lick_index(lick_indexes.keys(), spec['lambda'], spec['flux'], spec['error'],
                                          spec['flag'], v0=50, instrument_res=3, vd=100))
        specfile = 'bc2003_hr_m62_chab_ssp_100.spec'
        dt = np.dtype([('lambda', np.float), ('flux', np.float)])
        spec = np.loadtxt('%s/%s' % (data_example_dir, specfile), dtype=dt)
        print specfile + ':'
        pprint.pprint(self.lick_index(lick_indexes.keys(), spec['lambda'], spec['flux'], None, None, v0=50,
                                      instrument_res=3, vd=100))


if __name__ == '__main__':
    t = LickIndex()
    t.test()
