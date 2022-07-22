# Generated with SMOP  0.41-beta
# Gen Parameters
import numpy as np
from sympy import *


def define_sync_gen_params(input_gen_type):

    gen_param = {}
    gen_param.update({'Pm': 0.8089})
    gen_param.update({'omega0': 1})
    gen_param.update({'xl': 0.1})
    gen_param.update({'ra': 0})
    gen_param.update({'H': 2.15})
    gen_param.update({'xd': 1.4})
    gen_param.update({'xd0': 1.4})
    gen_param.update({'xdp': 0.3})
    gen_param.update({'xdpp': 0.2})
    gen_param.update({'xq': 1.35})
    gen_param.update({'xqp': 0.6})
    gen_param.update({'xqpp': 0.2})
    gen_param.update({'D': 0})
    gen_param.update({'Tdp': 6})
    if (gen_param['Tdp'] == 0):
        gen_param['Tdp'] = float('inf')

    gen_param.update({'Tdpp': 0.5})
    if (gen_param['Tdpp'] == 0):
        gen_param['Tdpp'] = float('inf')

    gen_param.update({'Tqp': 1})
    if (gen_param['Tqp'] == 0):
        gen_param['Tqp'] = float('inf')

    gen_param.update({'Tqpp': 0.05})
    if (gen_param['Tqpp'] == 0):
        gen_param['Tqpp'] = float('inf')

    if input_gen_type == 'H' :
        gen_param.update({'mbase': 2541})
    elif input_gen_type == 'B':
        gen_param.update({'mbase': 122})
    elif input_gen_type == 'G':
        gen_param.update({'mbase': 1272})
    else:
        raise ValueError('The MBASE value is missing')

    gen_param.update({'sbase': 100})

    gen_param.update({'S10': 0.03})
    gen_param.update({'S12': 0.4})
    A,B = symbols('A B')
    # Cuadratic Saturation Function
    #eqns = [B*(1-A)^2 == gen_param.S10, B*(1.2-A)^2 == gen_param.S12];
    # Scaled Quadratic Saturation Function
    eqn1 = Eq(B * (1 - A) ** 2, gen_param['S10'])
    eqn2 = Eq(B * (1.2 - A) ** 2 / 1.2, gen_param['S12'])
    # DefineGenParams.m:53
    # Exponential Saturation Function
    # eqns = [B*1^A == gen_param.S10, B*1.2^A == gen_param.S12];

    S=nonlinsolve([eqn1,eqn2],[A,B])
    gen_param.update({'A': float(S.args[0][0])})
    gen_param.update({'B': float(S.args[0][1])})

    gen_param.update({'Efd_max': 3})
    gen_param.update({'Efd_min': 0})
    gen_param.update({'Pmech_max': 10})
    gen_param.update({'Pmech_min': 0})
    gen_param.update({'K1': gen_param['xqp'] - gen_param['xl']})
    gen_param.update({'K2': gen_param['xdp'] - gen_param['xl']})
    gen_param.update({'K3': gen_param['xd'] - gen_param['xdp']})
    gen_param.update({'K4': (gen_param['xdp'] - gen_param['xdpp']) / (gen_param['xdp'] - gen_param['xl']) ** 2})
    gen_param.update({'K5': gen_param['xq'] - gen_param['xqp']})
    gen_param.update({'K6': (gen_param['xqp'] - gen_param['xdpp']) / (gen_param['xqp'] - gen_param['xl']) ** 2})
    gen_param.update({'K7': (gen_param['xdpp'] - gen_param['xl']) / (gen_param['xdp'] - gen_param['xl'])})
    gen_param.update({'K8': (gen_param['xdp'] - gen_param['xdpp']) / (gen_param['xdp'] - gen_param['xl'])})
    gen_param.update({'K9': (gen_param['xdpp'] - gen_param['xl']) / (gen_param['xqp'] - gen_param['xl'])})
    gen_param.update({'K10': (gen_param['xqp'] - gen_param['xdpp']) / (gen_param['xqp'] - gen_param['xl'])})
    gen_param.update({'K11': (gen_param['xq'] - gen_param['xl']) / (gen_param['xd'] - gen_param['xl'])})
    gen_param.update({'Kprop': 0})

    return gen_param

def define_renewable_gen_params(measurements):

    gen_param = {}

    gen_param['base2'] = 100
    gen_param.update({'base1': 1666})
    gen_param.update({'REPCA1_P_ref': np.dot(np.mean(measurements.Pe(np.arange(100,2000))),gen_param['base2']) / gen_param['base1']})
    gen_param.update({'REPCA1.Q_ref': np.dot(np.mean(measurements.Qe(np.arange(100,2000))),gen_param['base2']) / gen_param['base1']})
    gen_param.update({'REGCA1_Tg': 0.02})
    gen_param.update({'REECB1_Tiq': 0.02})
    gen_param.update({'REECB1_Trv': 0.02})
    gen_param.update({'REECB1_Tpord': 0.02})
    gen_param.update({'REECB1_Iqmax': 1.11})
    gen_param.update({'REECB1_Iqmin': - 1.11})
    gen_param.update({'REECB1_Ipmax': gen_param.REECB1.Iqmax})
    gen_param.update({'REECB1_Pmax': 1})
    gen_param.update({'REECB1_Pmin': 0})
    gen_param.update({'REPCA1_Tg': 0.1})
    gen_param.update({'REPCA1_Tfltr': 0.02})
    gen_param.update({'REPCA1_Tfv': 0.05})
    gen_param.update({'REPCA1_Tp': 0.25})
    gen_param.update({'REPCA1_Kp': 18})
    gen_param.update({'REPCA1_Ki': 5})
    gen_param.update({'REPCA1_Kpg': 1})
    gen_param.update({'REPCA1_Kig': 0.05})
    gen_param.update({'REPCA1_dbd1': - 1})
    gen_param.update({'REPCA1_dbd2': 1})
    gen_param.update({'REPCA1_emax': 0.1})
    gen_param.update({'REPCA1_emin': - 0.1})
    gen_param.update({'REPCA1_Pmax': 1})
    gen_param.update({'REPCA1_Pmin': 0})
    gen_param.update({'REPCA1_Qmax': 0.43})
    gen_param.update({'REPCA1_Qmin': - 0.43})

    measurements['vector'] = np.dot(measurements['vector'],gen_param['base2']) / gen_param['base1']

    return gen_param