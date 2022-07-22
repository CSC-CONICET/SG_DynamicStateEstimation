from project.import_data import load_data
from project.DSE import define_DSE_params, perform_DSE
from project.measurements import define_measurements
from project.kalman_filter.Code.graphics import plot_graphics_sync

case = 1
input_file_name = ['Efd_FO.mat', 'Pmech_FO']
input_file_playback = ['Efd_FO_playback.mat', 'Pmech_FO_playback']

input_gen_type = 'H' # Select a generator type - options = {H: Hydro, B: Biomass ,G: Gas}

input_data = load_data(input_file_name[case], input_file_playback[case], input_gen_type)
dse_param = define_DSE_params(input_data)
measurements = define_measurements(input_data, input_gen_type, dse_param)
kalman_vars = perform_DSE(input_gen_type, measurements, dse_param)

if dse_param['Graphs_boolean']:
    plot_graphics_sync(dse_param, kalman_vars, input_data)
