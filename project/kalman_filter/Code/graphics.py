# graphics.py
import matplotlib.pyplot as plt

def plot_graphics_sync(dse_param,kalman_vars,input_data):

    fig_path='/home/pablo/Dropbox/PosDoc/Codigo/DSE_FO_playback_v2 (copia)/'

    fig1, axs1 = plt.subplots(2, 4)
    fig1.suptitle('Generator state variables')
    axs1[0,0].plot(input_data['t'],input_data['delta'], 'r--',label='True')
    axs1[0,0].plot(dse_param['t_interp'],kalman_vars['state'][0,:],label='FMDD')
    axs1[0,0].legend()
    axs1[0,0].set_ylabel('$\delta$')
    axs1[1,0].plot(input_data['t'], input_data['speed'], 'r--', label='True')
    axs1[1,0].plot(dse_param['t_interp'],kalman_vars['state'][1,:])
    axs1[1,0].set_ylabel('$\omega$')
    axs1[1,0].set_xlabel('Time[s]')
    axs1[0,1].plot(input_data['t'], input_data['psikq'], 'r--', label='True')
    axs1[0,1].plot(dse_param['t_interp'],kalman_vars['state'][2,:])
    axs1[0,1].set_ylabel('$\Psi_q$')
    axs1[1,1].plot(input_data['t'], input_data['psikd'], 'r--', label='True')
    axs1[1,1].plot(dse_param['t_interp'],kalman_vars['state'][3,:])
    axs1[1,1].set_xlabel('Time[s]')
    axs1[1,1].set_ylabel('$\Psi_d$')
    axs1[0,2].plot(input_data['t'], input_data['epq'], 'r--', label='True')
    axs1[0,2].plot(dse_param['t_interp'],kalman_vars['state'][4,:])
    axs1[0,2].set_ylabel('$E\'_q$')
    axs1[1,2].plot(input_data['t'], input_data['epd'], 'r--', label='True')
    axs1[1,2].plot(dse_param['t_interp'],kalman_vars['state'][5,:])
    axs1[1,2].set_xlabel('Time[s]')
    axs1[1,2].set_ylabel('$E\'_d$')
    axs1[0,3].plot(dse_param['t_interp'],kalman_vars['state_UI'][0,:])
    axs1[0,3].plot(input_data['t'], input_data['pmech'], 'r--', label='True')
    axs1[0,3].set_ylabel('$P_{mech}$')
    axs1[1,3].plot(dse_param['t_interp'],kalman_vars['state_UI'][1,:])
    axs1[1,3].plot(input_data['t'], input_data['efd'], 'r--', label='True')
    axs1[1,3].set_xlabel('Time[s]')
    axs1[1,3].set_ylabel('$E_{fd}$')

    fig2, axs2 = plt.subplots(3, 1)
    fig2.suptitle('Kalman filter residuals')
    axs2[0].plot(dse_param['t_interp'],kalman_vars['residuals'][0,:])
    axs2[0].set_ylabel('$V_{re}$')
    axs2[1].plot(dse_param['t_interp'],kalman_vars['residuals'][1,:])
    axs2[1].set_ylabel('$V_{im}$')
    axs2[2].plot(dse_param['t_interp'],kalman_vars['residuals'][2,:])
    axs2[2].set_ylabel('$f$')
    axs2[2].set_xlabel('Time[s]')


def plot_graphics_renewables(dse_param, kalman_vars,input_data):
    fig_path = '/home/pablo/Dropbox/PosDoc/Codigo/DSE_FO_playback_v2 (copia)/'

    fig1, axs1 = plt.subplots(2, 4)
    axs1[0, 0].plot(dse_param['t_interp'], kalman_vars['state'][0, :])
    axs1[1, 0].plot(dse_param['t_interp'], kalman_vars['state'][1, :])
    axs1[0, 1].plot(dse_param['t_interp'], kalman_vars['state'][2, :])
    axs1[1, 1].plot(dse_param['t_interp'], kalman_vars['state'][3, :])
    axs1[0, 2].plot(dse_param['t_interp'], kalman_vars['state'][4, :])
    axs1[1, 2].plot(dse_param['t_interp'], kalman_vars['state'][5, :])
    axs1[0, 3].plot(dse_param['t_interp'], kalman_vars['state_UI'][0, :])
    axs1[1, 3].plot(dse_param['t_interp'], kalman_vars['state_UI'][1, :])

    fig2, axs2 = plt.subplots(3, 1)
    axs2[0].plot(dse_param['t_interp'], kalman_vars['residuals'][0, :])
    axs2[1].plot(dse_param['t_interp'], kalman_vars['residuals'][1, :])
    axs2[2].plot(dse_param['t_interp'], kalman_vars['residuals'][2, :])