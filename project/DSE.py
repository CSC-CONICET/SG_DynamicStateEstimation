
import numpy as np
from project.kalman_filter.Code.gen_params import define_sync_gen_params,define_renewable_gen_params
from project.kalman_filter.Code.kalman import define_kalman



def define_DSE_params(input_data):

    """

    Parameters
    ----------
    input_data {dict: 23}

    Returns
    -------
    dse_param {dict: 12}
    """

    dse_param={}
    dse_param.update({'f0': 60})

    dse_param.update({'diff_t': input_data['t'][1] - input_data['t'][0]})
    dse_param.update({'ReportRate' : round(1 / dse_param['diff_t'])})
    dse_param.update({'offset1' : 0})
    dse_param.update({'offset2' : len(input_data['V'])-1})
    dse_param.update({'interp_fact' : 1})
    dse_param.update({'delta_t' : dse_param['diff_t'] / dse_param['interp_fact']})
    dse_param.update({'TVE_max' : 1e-1})
    dse_param.update({'FE_max' : 5e-4})
    dse_param.update({'RFE_max' : 0.1})
    dse_param.update({'Graphs_boolean' : 1})
    dse_param.update({'t_interp' : np.arange(input_data['t'][dse_param['offset1']], input_data['t'][dse_param['offset2']],dse_param['delta_t'])})

    return dse_param


def perform_DSE(input_gen_type, measurements, dse_param):

    """

    Parameters
    ----------
    input_gen_type
    input_gen_case
    measurements
    dse_param

    Returns
    -------

    """

    # Generator model parameters are defined
    gen_params = define_sync_gen_params(input_gen_type)

    # UKF
    kalman_vars = define_kalman(dse_param, gen_params, measurements)

    return kalman_vars


