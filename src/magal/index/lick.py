import numpy as np
import scipy.ndimage.filters as sc
import scipy.interpolate as sci
from pystarlight.util.velocity_fix import SpectraVelocityFixer

Lick = np.array([['index_ind','name', 'Index band', 'blue continuum', 'red continuum'],\
                [1, 'CN_1', 4142.125, 4177.125, 4080.125, 4117.625, 4244.125, 4284.125],\
                [2, 'CN_2', 4142.125, 4177.125, 4083.875, 4096.375, 4244.125, 4284.125],\
                [3, 'Ca4227', 4222.250, 4234.750, 4211.000, 4219.750, 4241.000, 4251.000],\
                [4, 'G4300', 4281.375, 4316.375, 4266.375, 4282.625, 4318.875, 4335.125],\
                [5, 'Fe4383', 4369.125, 4420.375, 4359.125, 4370.375, 4442.875, 4455.375],\
                [6, 'Ca4455', 4452.125, 4474.625, 4445.875, 4454.625, 4477.125, 4492.125],\
                [7, 'Fe4531', 4514.250, 4559.250, 4504.250, 4514.250, 4560.500, 4579.250],\
                [8, 'Fe4668', 4633.999, 4720.250, 4611.500, 4630.250, 4742.750, 4756.500],\
                [9, 'H_beta', 4847.875, 4876.625, 4827.875, 4847.875, 4876.625, 4891.625],\
                [10, 'Fe5015', 4977.750, 5054.000, 4946.500, 4977.750, 5054.000, 5065.250],\
                [11, 'Mg_1', 5069.125, 5134.125, 4895.125, 4957.625, 5301.125, 5366.125],\
                [12, 'Mg_2', 5154.125, 5196.625, 4895.125, 4957.625, 5301.125, 5366.125],\
                [13, 'Mg_b', 5160.125, 5192.625, 5142.625, 5161.375, 5191.375, 5206.375],\
                [14, 'Fe5270', 5245.650, 5285.650, 5233.150, 5248.150, 5285.650, 5318.150],\
                [15, 'Fe5335', 5312.125, 5352.125, 5304.625, 5315.875, 5353.375, 5363.375],\
                [16, 'Fe5406', 5387.500, 5415.000, 5376.250, 5387.500, 5415.000, 5425.000],\
                [17, 'Fe5709', 5696.625, 5720.375, 5672.875, 5696.625, 5722.875, 5736.625],\
                [18, 'Fe5782', 5776.625, 5796.625, 5765.375, 5775.375, 5797.875, 5811.625],\
                [19, 'Na_D', 5876.875, 5909.375, 5860.625, 5875.625, 5922.125, 5948.125],\
                [20, 'TiO_1', 5936.625, 5994.125, 5816.625, 5849.125, 6038.625, 6103.625],\
                [21, 'TiO_2', 6189.625, 6272.125, 6066.625, 6141.625, 6372.625, 6415.125],\
                [22, 'Hdelta_A', 4083.500, 4122.250, 4041.600, 4079.750, 4128.500, 4161.000],\
                [23, 'Hgamma_A', 4319.750, 4363.500, 4283.500, 4319.750, 4367.250, 4419.750],\
                [24, 'Hdelta_F', 4090.999, 4112.250, 4057.250, 4088.500, 4114.750, 4137.250],\
                [25, 'Hgamma_F', 4331.250, 4352.250, 4283.500, 4319.750, 4354.750, 4384.750],\
                ])

def Lick_index(lambda_obs,flux_obs,error_obs, flags_obs, index,MgFe,D4000,rot_vel, res_ini, disp_vel, res_inst_fin=6, v_disp_fin=357):
    
    '''
    l_obs: 'float', array con lambda observado
    f_obs: 'float', array con flujo observado
    index: 'float', array con los indices a calcular
    MgFe: 'yes' or 'no'
    D4000: 'yes' or 'no
    galaxies: 'yes' or 'no'. si degradamos la resolucion a 14A o no
    Por defecto ponemos res_inst_fin=6 que es lo que tenemos en Califa. Si estamos midiendo sobre modelos, supongamos con res_ini = 3A FWHM, primero hacemos convolucion en lambda para poner en 6.
    Si queremos tambien podemos dejar la resolucion en 3FWHM, pero el programa da la opcion de cambiar.
    Si estamos midiendo en CALIFA, el espectro tendra res_ini = 6, res_inst_fin = 6, con lo que no hara convolucion en lambda.
    Por defecto colocaremos v_disp_fin = 357 km/s. Esto equivale a 14A (FWHM) en 5000 A
    '''
    
    if flags_obs != 'no':
        flags = flags_obs
        ind_flag = np.where(flags==0)
        l_obs = lambda_obs[ind_flag]
        f_obs = flux_obs[ind_flag]
        error = error_obs[ind_flag]
    else:
        l_obs = lambda_obs
        f_obs = flux_obs
        flags = np.zeros_like(f_obs)
        error = np.zeros_like(f_obs)
    lamb = np.arange(3800, 6801)
    delta_lambda_obs = 1
    flux = np.interp(lamb, l_obs, f_obs)
    flags = np.interp(lamb, lambda_obs, flags)
    error = np.interp(lamb, l_obs, error)
    indices = []
    errores = []
    name_indices = []
    
    #ahora degradamos la resolucion. Primero en lambda por si queremos pasar por ejemplo 3A (FWHM) de los modelos a 6A (FWHM) de la resolucion instrumental de califa
    sigma_instrumental = res_ini/(2*np.sqrt(2*np.log(2)))
    sigma_l_fin = res_inst_fin/(2*np.sqrt(2*np.log(2)))
    sigma_dif = np.sqrt(sigma_l_fin**2 - sigma_instrumental**2)
    flux = sc.gaussian_filter1d(flux,sigma_dif)
    
    sigma_instrumental = res_inst_fin/(2*np.sqrt(2*np.log(2)))
    sigma_d = (disp_vel/300000.)*5000.
    sigma_ini = np.sqrt(sigma_instrumental**2 + sigma_d**2)
    vel_ini = (sigma_ini/5000.)*300000
    
    ######bloque para corregir de velocidad de dispersion (funcion de Andre). Dando v_0 a la funcion tambien corrige de velocidad de rotacion.
    if vel_ini < v_disp_fin:
        vfix = SpectraVelocityFixer(lamb, rot_vel,vel_ini, nproc=1)
        flux = vfix.fix(flux, v_disp_fin)
        #flags = vfix.fix(flags, v_disp_fin)[0]
        #error = vfix.fix(flags, v_disp_fin)[0]
    else:
        flux = flux
        flags = flags
        error = error
    
    if MgFe == 'yes':
        aux = np.array([13,14,15])
        list = np.append(index,aux)
    else:
        list = index
    
    M_b_fac ='no'
    Fe5270_fac = 'no'
    Fe5335_fac = 'no'
    for ii in list:
        X_pos_b = np.where(flags[np.bitwise_and(lamb > Lick[ii][4], lamb < Lick[ii][5])] != 0)
        X_pos_r = np.where(flags[np.bitwise_and(lamb > Lick[ii][6], lamb < Lick[ii][7])] != 0)
        X_pos_l = np.where(flags[np.bitwise_and(lamb > Lick[ii][2], lamb < Lick[ii][3])] != 0)
        
        len_b = Lick[ii][5]-Lick[ii][4]
        len_r = Lick[ii][7]-Lick[ii][6]
        len_l = Lick[ii][3]-Lick[ii][2]
        
        if (len(X_pos_b) < 0.4*len_b/delta_lambda_obs) & (len(X_pos_r) < 0.4*len_r/delta_lambda_obs) & (len(X_pos_l) < 0.4*len_l/delta_lambda_obs):
            mean_blue = np.mean(flux[np.bitwise_and(lamb > Lick[ii][4], lamb < Lick[ii][5])])
            mean_red = np.mean(flux[np.bitwise_and(lamb > Lick[ii][6], lamb < Lick[ii][7])])
    
            p_b = (Lick[ii][4] + Lick[ii][5])/2
            p_r = (Lick[ii][6] + Lick[ii][7])/2
    
            # recta definida por los dos continuos:
            # ((mean_blue - mean_red)/(p_b - p_r))*x + (mean_red*p_b - mean_blue*p_r)/(p_b-p_r)
            l_low = 0
            if int(Lick[ii][2]) == Lick[ii][2]:
                l_low = Lick[ii][2]
            else:
                l_low = int(Lick[ii][2]) + 1
    
            range = np.arange(l_low, int(Lick[ii][3])+1)
            line = ((mean_blue - mean_red)/(p_b - p_r))*range + (mean_red*p_b - mean_blue*p_r)/(p_b-p_r)
            ind = np.sum((line-flux[np.bitwise_and(lamb >= Lick[ii][2], lamb <= Lick[ii][3])])/line)
            if ind < 0:
                if (Lick[ii][1] == 'Mg_b')|(Lick[ii][1] == 'Fe5270')|(Lick[ii][1] == 'Fe5335'):
                    ind = 0
            err_ind = np.sum((error[np.bitwise_and(lamb >= Lick[ii][2], lamb <= Lick[ii][3])]**2)/(line**2))
        
            if Lick[ii][1] == 'Mg_b':
                index_Mg_b = ind
                M_b_fac = 'yes'
                err_Mg_b = np.sqrt(err_ind)
            if Lick[ii][1] == 'Fe5270':
                index_Fe5270 = ind
                Fe5270_fac = 'yes'
                err_Fe5270 = np.sqrt(err_ind)
            if Lick[ii][1] == 'Fe5335':
                index_Fe5335 = ind
                Fe5335_fac = 'yes'
                err_Fe5335= np.sqrt(err_ind)

            
            if ii in index:
                indices = np.append(indices,ind)
                errores = np.append(errores,np.sqrt(err_ind))
                name_indices = np.append(name_indices,Lick[ii][1])
        else:
            print 'No se puede calcular el indice ' + Lick[ii][1]
                    
    
    if (MgFe == 'yes'):
        if (M_b_fac == 'yes')&(Fe5270_fac == 'yes')&(Fe5335_fac == 'yes'):
            x = index_Mg_b
            y = index_Fe5270
            z = index_Fe5335
            Ex = err_Mg_b
            Ey = err_Fe5270
            Ez = err_Fe5335
        
            MgFe =  np.sqrt(index_Mg_b*(0.72 * index_Fe5270 + 0.28 * index_Fe5335))
            if MgFe == 0:
                MgFe_e = 0.0001
                error_MgFe = 0.01
            else:
                MgFe_e = MgFe
                error_MgFe = (1/(4*(MgFe_e**2)))*((Ex**2)*(0.72 * y + 0.28 * z)**2 + (Ey**2)*(0.72 * x)**2 + (Ez**2)*(0.28 * x)**2)
            indices = np.append(indices,MgFe_e)
            errores = np.append(errores,np.sqrt(error_MgFe))
            name_indices = np.append(name_indices,'MgFe')
        else:
            print 'No se puede calcular MgFe'
    
    if D4000 == 'yes':
        Y_pos_b = np.where(flags[np.bitwise_and(lamb > 3850, lamb < 3950)] != 0)
        Y_pos_r = np.where(flags[np.bitwise_and(lamb > 4000, lamb < 4100)] != 0)
        if (len(Y_pos_b) < 0.4*100/delta_lambda_obs) & (len(Y_pos_r) < 0.4*100/delta_lambda_obs):
            d4_b = np.sum(flux[np.bitwise_and(lamb >= 3850, lamb <= 3950)])
            d4_r = np.sum(flux[np.bitwise_and(lamb >= 4000, lamb <= 4100)])
            d4000 = d4_r/d4_b
            sigma_2_b = np.sum(error[np.bitwise_and(lamb >= 3850, lamb <= 3950)]**2)/(d4_b**2)
            sigma_2_r = np.sum(error[np.bitwise_and(lamb >= 4000, lamb <= 4100)]**2)/(d4_r**2)
            err_d4000 = d4000*np.sqrt(sigma_2_b + sigma_2_r)
            indices = np.append(indices,d4000)
            errores = np.append(errores, err_d4000)
            name_indices = np.append(name_indices,'D4000')
        else: 
            print 'No se puede calcular D4000'
    
    if flags_obs != 'no':
        return name_indices, indices, errores
    else:
        return name_indices, indices


def photometry(l_obs, f_obs, l_filter, f_filter, redshift, magnitudeAB, starlight):
    
    '''
    l_obs: 'float', array con lambda observado
    f_obs: 'float', array con flujo obervado
    dir_filter: 'str', direccion del filtro, incluido nombre del fichero
    redshift: 'float', redshift del espectro observado
    magnitudeAB: 'str', 'yes' o 'no', para saber si queremos en magnitud o flujo, se espera que el flujo observado este en erg/s/cm^2/A
    starlight: 'str', 'yes' o 'no', si es un input de starlight, solo va multiplicado por (1+z), no por (1+z)**3
    '''
    
    z_pho = redshift
    lambda_min = np.min(l_filter)
    lambda_max = np.max(l_filter)
    lambda_filter = np.arange(np.int(lambda_min),np.int(lambda_max)+1,1)
    filter_curve = np.interp(lambda_filter, l_filter,f_filter)
    l_obs = l_obs*(1+z_pho)
    f_int = sci.interp1d(l_obs,f_obs)
    flux_filter = f_int(lambda_filter)
    if starlight == 'yes':
        flux_filter = flux_filter/(1+z_pho)
    elif starlight == 'no':
        flux_filter = flux_filter/((1+z_pho)**3)
    
    if magnitudeAB == 'no':
        PHO = np.sum(filter_curve*flux_filter)/np.sum(filter_curve)
    elif magnitudeAB == 'yes':
        A = np.trapz(lambda_filter*flux_filter*filter_curve,lambda_filter)
        B = np.trapz(filter_curve/lambda_filter,lambda_filter)
    
    PHO = -2.5*np.log10(A/B) - 2.41
    
    return PHO