import numpy as np
import scipy.io


def flatten(d):
    for a, b in d.items():
        try:
            d[a] = d[a].flatten()
        except AttributeError:
            pass
    return d


def load_data(input_file, input_file_playback, input_gen_type):

    # To make sure that there is no “nested lists" we use the flatten function
    input_data= flatten(scipy.io.loadmat('DataIn/'+input_file))

    IT = (input_data['PT'] - np.dot(1j, input_data['QT'])) / np.conj(np.multiply(input_data['V'], np.exp(np.dot(1j, (input_data['angV'])))))

    input_data2 = flatten(scipy.io.loadmat('DataIn/'+input_file_playback))

    IB = (input_data2['Pe_B'] - np.dot(1j, input_data2['Qe_B'])) / np.conj(np.multiply(input_data['V'],np.exp(np.dot(1j, (input_data['angV'])))))
    IG = (input_data2['Pe_G'] - np.dot(1j, input_data2['Qe_G'])) / np.conj(np.multiply(input_data['V'],np.exp(np.dot(1j, (input_data['angV'])))))
    IH = (input_data2['Pe_H'] - np.dot(1j, input_data2['Qe_H'])) / np.conj(np.multiply(input_data['V'], np.exp(np.dot(1j, (input_data['angV'])))))
    IS = (input_data2['Pe_S'] - np.dot(1j, input_data2['Qe_S'])) / np.conj(np.multiply(input_data['V'],np.exp(np.dot(1j, (input_data['angV'])))))
    IW = (input_data2['Pe_W'] - np.dot(1j, input_data2['Qe_W'])) / np.conj(np.multiply(input_data['V'],np.exp(np.dot(1j, (input_data['angV'])))))

    if input_gen_type == 'H':
        I = IT - IB - IG - IS - IW
    elif input_gen_type == 'B':
        I = IT - IH - IG - IS - IW
    elif input_gen_type == 'G':
        I = IT - IB - IH - IS - IW
    else:
        raise ValueError('Select a correct generator type: options {H,B,G}')

    input_data.update({'angI' : np.angle(I)})
    input_data.update({'I' : abs(I)})

    return input_data