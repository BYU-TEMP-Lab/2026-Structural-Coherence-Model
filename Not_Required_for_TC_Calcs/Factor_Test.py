import numpy as np

def P_factor(NI, KF, PH, DF):
    
    P_ph = NI * KF * PH
    P_dif = NI * KF * (1 - PH) * DF

    return P_ph, P_dif

# Example usage for LiF
x_range = np.linspace(0,14,100)
transfer_len = 2.5
transfer_points = np.arange(transfer_len, np.max(x_range),transfer_len)
cc_peak = 3.8
cc_transfer = False

NI = 0.95
KF = 0.86
ph = [0.8, 0.3, 0.3, 0.3, 0.3]
DF = ( (transfer_points[0]/cc_peak)**3 - (1/2)**3 ) / (1 - (1/2)**3)
if DF < 0:
    DF = 0  # Ensure DF is non-negative

transfer_points_with_cc = []
if cc_peak > transfer_points[0] and cc_peak < transfer_points[1]:
    cc_transfer = True
    cc_transfer_pnts = np.arange(1, len(transfer_points)) * cc_peak
    for i in range(len(transfer_points)):
        transfer_points_with_cc.append(transfer_points[i])
        if i < len(transfer_points) - 1:
            transfer_points_with_cc.append(cc_transfer_pnts[i])
    print(transfer_points_with_cc)
    transfer_points = transfer_points_with_cc

b_disruption = np.ones(len(transfer_points))  # Initialize b_disruption array
S_y = np.zeros(len(transfer_points))
beta_i = 0
beta_i_integral = 0

i_ca = 0
for i in range(len(transfer_points)-1):
    PH = ph[i_ca]
    P_ph, P_dif = P_factor(NI, KF, PH, DF)
    
    if cc_transfer and i % 2 == 1:

        b_disruption[i] += - P_dif

        beta_i = b_disruption[i] / (1 - b_disruption[i])
        beta_i_integral += beta_i * (transfer_points[i+1] - transfer_points[i])

        S_y[i] = (np.exp(-beta_i_integral))

    else:

        b_disruption[i] += - P_ph - P_dif

        beta_i = b_disruption[i] / (1 - b_disruption[i])
        beta_i_integral += beta_i * (transfer_points[i+1] - transfer_points[i])

        S_y[i] = (np.exp(-beta_i_integral))

        i_ca += 1



np.plot(transfer_points, S_y, marker='o')
    #print(f"Total normalized failure probability: {failure_probability:.4f}")
