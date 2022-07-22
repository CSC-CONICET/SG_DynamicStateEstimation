import numpy as np


def define_measurements(input_data,input_gen_type,dse_param):

    """

    Parameters
    ----------
    input_data {dict: 23}
    input_gen_type {str}
    dse_param {dict: 12}

    Returns
    -------
    measurements {dict: 30}
    """

    measurements = {}
    measurements.update({'l' : len(input_data['angV'][np.arange(dse_param['offset1'], dse_param['offset2'])])})
    measurements.update({'TVE_max': dse_param['TVE_max'] / 100})
    measurements.update({'FE_max': dse_param['FE_max'] / dse_param['f0']})
    measurements.update({'RFE_max': dse_param['RFE_max'] / dse_param['f0']})
    measurements.update({'Meas_NoiseVar': (measurements['TVE_max'] / (np.sqrt(2) * 3)) ** 2})
    measurements.update({'f_NoiseStd': measurements['FE_max'] / 3})
    measurements.update({'Meas_NoiseStd': np.sqrt(measurements['Meas_NoiseVar'])})

    # Current and voltage phasors.
    measurements.update({'volt': np.multiply(input_data['V'][np.arange(dse_param['offset1'], dse_param['offset2'])],np.exp(np.dot(1j,input_data['angV'][np.arange(dse_param['offset1'],dse_param['offset2'])])))})
    measurements.update({'volt_ref': measurements['volt']})
    measurements.update({'cur': np.multiply(input_data['I'][np.arange(dse_param['offset1'], dse_param['offset2'])],np.exp(np.dot(1j,input_data['angI'][np.arange(dse_param['offset1'],dse_param['offset2'])])))})
    measurements.update({'cur_ref': measurements['cur']})
    # Derivative of the voltage/current amplitude/phase
    measurements.update({'abs_voltp': np.diff(abs(measurements['volt'])) / dse_param['diff_t']})
    measurements.update({'abs_voltp': np.hstack((measurements['abs_voltp'],measurements['abs_voltp'][-1]))})
    measurements.update({'phase_voltp': np.diff(np.angle(measurements['volt'])) / dse_param['diff_t']})
    measurements.update({'phase_voltp': np.hstack((measurements['phase_voltp'],measurements['phase_voltp'][-1]))})
    measurements.update({'abs_curp': np.diff(abs(measurements['cur'])) / dse_param['diff_t']})
    measurements.update({'abs_curp': np.hstack((measurements['abs_curp'],measurements['abs_curp'][-1]))})
    measurements.update({'phase_curp': np.diff(np.angle(measurements['cur'])) / dse_param['diff_t']})
    measurements.update({'phase_curp': np.hstack((measurements['phase_curp'], measurements['phase_curp'][-1]))})

    # Phasor derivative
    measurements.update({'volt_prime':  np.multiply(measurements['abs_voltp'], np.cos(np.angle(measurements['volt']))) - np.multiply(np.multiply(abs(measurements['volt']),np.sin(np.angle(measurements['volt']))),measurements['phase_voltp']) + np.dot(1j,(np.multiply(measurements['abs_voltp'],np.sin(np.angle(measurements['volt']))) + np.multiply(np.multiply(abs(measurements['volt']),np.cos(np.angle(measurements['volt']))),measurements['phase_voltp'])))})
    measurements.update({'volt_prime_ref': measurements['volt_prime']})
    measurements.update({'cur_prime': np.multiply(measurements['abs_curp'], np.cos(np.angle(measurements['cur']))) - np.multiply(np.multiply(abs(measurements['cur']),np.sin(np.angle(measurements['cur']))),measurements['phase_curp']) + np.dot(1j,(np.multiply(measurements['abs_curp'],np.sin(np.angle(measurements['cur']))) + np.multiply(np.multiply(abs(measurements['cur']),np.cos(np.angle(measurements['cur']))),measurements['phase_curp'])))})
    measurements.update({'cur_prime_ref': measurements['cur_prime']})
    measurements.update({'f_est': 1 + np.dot(1 / (np.dot(np.dot(2,np.pi), dse_param['f0'])),np.imag(measurements['volt_prime'] / measurements['volt_ref']))})
    measurements.update({'f_est_ref': measurements['f_est']})
    measurements.update({'rocof': np.diff(measurements['f_est'] / dse_param['diff_t'])})
    measurements.update({'rocof': np.hstack((measurements['rocof'], measurements['rocof'][-1]))})

    # Additive Noise
    np.random.seed(1)
    measurements.update({'volt': measurements['volt'] + np.dot(np.dot(np.mean(abs(measurements['volt'])),measurements['Meas_NoiseStd']),(np.random.randn(measurements['l'],) + np.dot(1j,np.random.randn(measurements['l'],))))})
    measurements.update({'cur': measurements['cur'] + np.dot(np.dot(np.mean(abs(measurements['cur'])),measurements['Meas_NoiseStd']),(np.random.randn(measurements['l'],) + np.dot(1j,np.random.randn(measurements['l'],))))})
    measurements.update({'f_est': measurements['f_est'] + np.dot(measurements['f_NoiseStd'],np.random.randn(measurements['l'],))})

    # Data interpolation
    measurements.update({'volt': np.interp(dse_param['t_interp'],input_data['t'][np.arange(dse_param['offset1'], dse_param['offset2'])],np.real(measurements['volt'])) + np.dot(1j,np.interp(dse_param['t_interp'],input_data['t'][np.arange(dse_param['offset1'],dse_param['offset2'])],np.imag(measurements['volt'])))})
    measurements.update({'cur':  np.interp(dse_param['t_interp'],input_data['t'][np.arange(dse_param['offset1'], dse_param['offset2'])],np.real(measurements['cur'] )) + np.dot(1j,np.interp(dse_param['t_interp'],input_data['t'][np.arange(dse_param['offset1'],dse_param['offset2'])],np.imag(measurements['cur'] )))})
    measurements.update({'Pe': np.real(np.multiply(measurements['volt'],np.conj(measurements['cur'])))})
    measurements.update({'Qe': np.imag(np.multiply(measurements['volt'],np.conj(measurements['cur'])))})
    measurements.update({'frequency': np.interp(dse_param['t_interp'],input_data['t'][np.arange(dse_param['offset1'], dse_param['offset2'])], measurements['f_est'])})
    measurements.update({'rocof': np.interp(dse_param['t_interp'],input_data['t'][np.arange(dse_param['offset1'], dse_param['offset2'])], measurements['rocof'])})

    # Measurement vector
    measurements.update({'vector': np.vstack((np.real(measurements['volt']), np.imag(measurements['volt']), measurements['frequency']))})

    measurements.update({'m': np.size(measurements['vector'], 0)})
    measurements.update({'l': np.size(measurements['vector'], 1)})

    return measurements
