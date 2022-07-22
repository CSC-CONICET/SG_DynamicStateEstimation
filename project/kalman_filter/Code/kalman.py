from project.kalman_filter.Functions.ukf_UI import ukf_UI_sinc
import numpy as np


def define_kalman(dse_param, gen_params, measurements):
    ## Algorithm

    kalman_vars, x, P = define_kalman_filter_ini(measurements)

    for i in np.arange(0,measurements['l']).reshape(-1):
        u = define_kalman_input(measurements, gen_params,i)
        z = measurements['vector'][:,i]

        x,P,residuals,residual_cov, d, P_d = ukf_UI_sinc(x, P, z, u, kalman_vars['Qk'], kalman_vars['Rk'], gen_params, dse_param)

        kalman_vars['state'][:,i]=x
        kalman_vars['state_UI'][:,i]= d.T
        kalman_vars['residuals_cov'][:,:,i]=residual_cov
        kalman_vars['residuals'][:,i]=residuals

    return kalman_vars


def define_kalman_filter_ini(measurements):
    kalman_vars = {}

    kalman_vars.update({'x_ini': np.array([np.angle(measurements['volt'][0]), 1, 0.5, 0.9, 1, 0.26])})

    kalman_vars.update({'state': np.zeros((np.size(kalman_vars['x_ini']), measurements['l']))})
    kalman_vars.update({'residuals': np.zeros((measurements['m'], measurements['l']))})

    kalman_vars.update({'residuals_cov': np.zeros((measurements['m'], measurements['m'], measurements['l']))})
    kalman_vars.update({'xk': kalman_vars['x_ini']})
    kalman_vars.update({'Pk': np.dot(1e-08, np.eye(np.size(kalman_vars['xk'])))})
    kalman_vars.update({'Qk': np.dot(1e-12, np.eye(np.size(kalman_vars['xk'])))})
    kalman_vars['Qk'][1][1] = 1e-14
    kalman_vars['Qk'][3][3] = 1e-12

    # Measurement covariance
    kalman_vars.update({'Rk': np.dot(measurements['Meas_NoiseVar'], np.eye(measurements['m']))})
    kalman_vars['Rk'][2, 2] = np.dot(measurements['f_NoiseStd'] ** 2, 100)

    x = kalman_vars['x_ini']
    P = kalman_vars['Pk']
    d = np.array([0.8, 2])

    kalman_vars['state'][:, 0] = kalman_vars['xk']
    kalman_vars.update({'state_UI': np.zeros((np.size(d), measurements['l']))})
    kalman_vars['state_UI'][:, 0] = d

    return kalman_vars,x, P


def define_kalman_input(measurements, gen_params,i):
    kalman_input = {}

    kalman_input.update({'vk': measurements['volt'][i]})
    kalman_input.update({'ik': np.dot(measurements['cur'][i], gen_params['sbase']) / gen_params['mbase']})
    kalman_input.update({'Pe': np.real(np.dot(kalman_input['vk'], np.conj(kalman_input['ik'])))})
    kalman_input.update({'Qe': np.imag(np.dot(kalman_input['vk'], np.conj(kalman_input['ik'])))})
    kalman_input.update({'ir': np.real(kalman_input['ik'])})
    kalman_input.update({'ii': np.imag(kalman_input['ik'])})
    kalman_input.update({'V': abs(kalman_input['vk'])})
    kalman_input.update({'theta': np.angle(kalman_input['vk'])})
    kalman_input.update({'I': abs(kalman_input['ik'])})
    kalman_input.update({'phi': np.angle(kalman_input['ik'])})

    return kalman_input