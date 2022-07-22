# ukf_UI
from numpy import *
import copy


def ukf_UI_sinc(x=None, P=None, z=None, u=None, Q=None, R=None, gen_param=None, dse_param=None, *args,
                **kwargs):
    # UKF   Unscented Kalman Filter for nonlinear dynamic systems
    # [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P
    # for nonlinear dynamic system (for simplicity, noises are assumed as additive):
    #           x_k+1 = f(x_k) + G d_k + w_k
    #           z_k   = h(x_k) + v_k
    # where
    #           w ~ N(0,Q) meaning w is gaussian noise with covariance Q
    #           v ~ N(0,R) meaning v is gaussian noise with covariance R
    #           d is the unknown input
    #           G is the gradient of g respect to x (where g(x)=d)
    # Inputs:   f: function handle for f(x)
    #           x: "a priori" state estimate
    #           P: "a priori" estimated state covariance
    #           h: function handle for h(x)
    #           z: current measurement
    #           Q: process noise covariance
    #           R: measurement noise covariance
    # Output:   x: "a posteriori" state estimate
    #           P: "a posteriori" state covariance
    #           y_tilde: residuals vector
    #           P2: Covariance of residuals vector
    #           d: "a posteriori" unknown input
    #           P_d: "a posteriori" covariance of the unknown input

    L = size(x)
    m = size(z)
    alpha = 0.001
    ki = 0
    beta = 2
    lambda_ = dot(alpha ** 2, (L + ki)) - L

    c = L + lambda_

    Wm = hstack((reshape(array([lambda_ / c]), (1, 1)), 0.5 / c + zeros((1, 2 * L))))
    Wm = Wm[0]
    Wc = copy.deepcopy(Wm)
    Wc[0] = Wc[0] + (1 - alpha ** 2 + beta)

    c = sqrt(c)
    G = zeros((size(x), 2))
    G[1, 0] = dse_param['delta_t'] / (2 * gen_param['H'] * x[1])
    G[4, 1] = dse_param['delta_t'] / gen_param['Tdp']

    # Estimation of the UI
    X = sigmas(x, P, c)
    d = transpose(zeros(G.ndim))
    x1, X1, P1, X2 = ut_fstate(X, Wm, Wc, L, Q, gen_param, dse_param, u, G, d)
    z1, __, P2, Z2 = ut_hmeas(X1, Wm, Wc, m, R, gen_param, u)

    P12 = dot(dot(X2, diagflat(Wc)), Z2.T)
    H = dot(P12.T, linalg.inv(P1))
    y_tilde = z - z1
    S = dot(H, G)
    R_tilde = dot(dot(H, P1), H.T) + R
    W = linalg.inv(R_tilde)
    P_d = linalg.inv((dot(dot(S.T, W), S)))
    M = dot(dot(P_d, S.T), W)
    d = dot(M, y_tilde.T)
    d[0] = min(gen_param['Pmech_max'], max(gen_param['Pmech_min'], d[0]))
    d[1] = min(gen_param['Efd_max'], max(gen_param['Efd_min'], d[1]))

    # Updating state vector
    x1, X1, P1, X2 = ut_fstate(X, Wm, Wc, L, Q, gen_param, dse_param, u, G, d)
    z1, __, P2, Z2 = ut_hmeas(X1, Wm, Wc, m, R, gen_param, u)
    P12 = dot(dot(X2, diagflat(Wc)), Z2.T)
    K = dot(P12, linalg.inv(P2))
    y_tilde = z - z1

    x = x1 + transpose(dot(K, y_tilde.T))
    P = P1 - dot(dot(K, P2), K.T)

    return x.flatten(), P, y_tilde, P2, d, P_d


def ut_fstate(X=None, Wm=None, Wc=None, n=None, Q=None, gen_param=None, dse_param=None, u=None, G=None,
               d=None, *args, **kwargs):
    # Unscented Transformation
    # Input:
    #        f: nonlinear map
    #        X: sigma points
    #       Wm: weights for mean
    #       Wc: weights for covraiance
    #        n: numer of outputs of f
    #        Q: additive covariance
    # Output:
    #        y: transformed mean
    #        Y: transformed smapling points
    #        P: transformed covariance
    #       Y1: transformed deviations

    L = X.shape[1]

    Y = fstate_sat(X, gen_param, dse_param, u) + reshape(repeat(dot(G, d),L,axis=0),(X.shape[0],X.shape[1]))

    y = reshape(dot(Y,Wm),(1,n))
    Y1 = Y - dot(y.T, ones((1, L)))
    P = dot(dot(Y1, diagflat(Wc)), Y1.T) + Q

    return y, Y, P, Y1


def ut_hmeas(X=None, Wm=None, Wc=None, n=None, R=None, gen_param=None, u=None, *args, **kwargs):
    # Unscented Transformation
    # Input:
    #        f: nonlinear map
    #        X: sigma points
    #       Wm: weights for mean
    #       Wc: weights for covraiance
    #        n: numer of outputs of f
    #        R: additive covariance
    # Output:
    #        y: transformed mean
    #        Y: transformed smapling points
    #        P: transformed covariance
    #       Y1: transformed deviations

    L = X.shape[1]

    Y = hmeas(X, gen_param, u)

    y = reshape(dot(Y, Wm), (1, n))
    Y1 = Y - dot(y.T, ones((1, L)))
    P = dot(dot(Y1, diagflat(Wc)), Y1.T) + R

    return y, Y, P, Y1


def sigmas(x=None, P=None, c=None, *args, **kwargs):
    # Sigma points around reference point
    # Inputs:
    #       x: reference point
    #       P: covariance
    #       c: coefficient
    # Output:
    #       X: Sigma points

    A = reshape(c * linalg.cholesky(P), (x.size, x.size))
    Y = dot(x.reshape(size(x), 1), ones((1, A.shape[1])))
    X = hstack((hstack((x.reshape(x.size, 1), Y + A)), Y - A))
    return X


def fstate_sat(x=None, gen_param=None, dse_param=None, u=None, *args, **kwargs):
    # f: nonlinear map - x_k+1 = f(x_k,u_k)
    # state variables: x = [delta, omega, psi_q, psi_d, ep_q, ep_d]

    delta = x[0,:]
    omega = x[1,:]
    psi_q = x[2,:]
    psi_d = x[3,:]
    ep_q  = x[4,:]
    ep_d  = x[5,:]

    delta_t = dse_param['delta_t']
    sin_delta = sin(delta)
    cos_delta = cos(delta)
    H     = gen_param['H']
    D     = gen_param['D']
    ir    = u['ir']
    ii    = u['ii']
    K1    = gen_param['K1'] # Generator parameters constants
    K2 = gen_param['K2']
    K3 = gen_param['K3']
    K4 = gen_param['K4']
    K5 = gen_param['K5']
    K6 = gen_param['K6']
    K7    = gen_param['K7']
    K8 = gen_param['K8']
    K9 = gen_param['K9']
    K10= gen_param['K10']
    K11 = gen_param['K11']

    x = array(
        [delta + 2 * pi * dse_param['f0'] * (omega - 1) * delta_t,
         omega + 1 / (H * 2) * ((0 - D * (omega - 1)) / omega - ((ep_q * K7 + psi_d * K8) * (ii *  sin_delta + ir * cos_delta) + (ep_d * K9 + psi_q * K10) * (ir *  sin_delta - ii * cos_delta))) * delta_t,
         psi_q + (-psi_q + K1 * (ii *  sin_delta + ir * cos_delta) + ep_d) / gen_param['Tqpp'] * delta_t,
         psi_d + (-psi_d - K2 * (ir *  sin_delta - ii * cos_delta) + ep_q) / gen_param['Tdpp'] * delta_t,
         ep_q + (0 - (ep_q + K7 * K3 * (ir *  sin_delta - ii * cos_delta) + K4 * K3 * (ep_q - psi_d) + gen_param['B'] * (abs(K8 * psi_d + K7 * ep_q + 1j * (- K10 * psi_q - K9 * ep_d)) - gen_param['A']) ** 2 / abs(K8 * psi_d + K7 * ep_q + 1j * (- K10 * psi_q - K9 * ep_d)) ** 1 * (  K8  * psi_d + K7 * ep_q))) / gen_param['Tdp'] * delta_t,
         ep_d + (  -  ep_d + K9 * K5 * (ii *  sin_delta + ir * cos_delta) - K6 * K5 * (ep_d - psi_q) + gen_param['B'] * (abs(K8 * psi_d + K7 * ep_q + 1j * (- K10 * psi_q - K9 * ep_d)) - gen_param['A']) ** 2 / abs(K8 * psi_d + K7 * ep_q + 1j * (- K10 * psi_q - K9 * ep_d)) ** 1 * (- K10 * psi_q - K9 * ep_d) * K11) / gen_param['Tqp'] * delta_t])

    return x


def hmeas(x=None, gen_param=None, kalman_input=None, *args, **kwargs):
    # h: nonlinear map - z_k = h(x_k,u_k)
    # state variables: x = [delta, omega, psi_q, psi_d, ep_q, ep_d]
    # measurement vector: z = [V_re V_im f]

    sin1 = sin(x[0,:])
    cos1 = cos(x[0,:])

    z = array([((x[2,:] * gen_param['K10'] + x[5,:] * gen_param['K9']) * sin1 +
                (x[3,:] * gen_param['K8'] + x[4,:] * gen_param['K7']) * cos1) * x[1,:] -
               kalman_input['ir'] * gen_param['ra'] + kalman_input['ii'] * gen_param['xdpp'],
               ((-x[2,:] * gen_param['K10'] - x[5,:] * gen_param['K9']) * cos1 +
                (x[3,:] * gen_param['K8'] + x[4,:] * gen_param['K7']) * sin1) * x[1,:] -
               kalman_input['ii'] * gen_param['ra'] - kalman_input['ir'] * gen_param['xdpp'],
               x[1,:]])
    return z
