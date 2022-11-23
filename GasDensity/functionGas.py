from sklearn.feature_extraction import DictVectorizer
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import cmath
import imageio as iio
import os
from matplotlib.ticker import LinearLocator
from scipy.signal import savgol_filter


def dict2mat(Dict, dimensional):
    """
    :param Dict: Dictionary data
    :param dimensional: integer, 1 or 2
    :return: A matrix with dimension in input
    """
    if dimensional == 1:
        dictvectorizer = DictVectorizer(sparse=False)
        Matrix1 = dictvectorizer.fit_transform(Dict)
        return Matrix1
    if dimensional == 2:
        temp_x = max([cord[0] for cord in Dict.keys()])
        temp_y = max([cord[1] for cord in Dict.keys()])
        Dict2 = [[0] * (temp_y + 1) for ele in range(temp_x + 1)]
        for (i, j), val in Dict.items():
            Dict2[i][j] = val
        Matrix2 = np.array(Dict2).reshape(temp_x + 1, -1)

        return Matrix2


def get_structure(fid, title):
    """
    :param fid: file ID
    :return: Img: Read in image.
        width: The fringe line width in the image.
        k: Fringe number in image.
        X: Fringe coordinate value in X axis.
        Y: Fringe coordinate value in X axis.
        Z: Fringe value.
        a_single: New image with fringe in 1 pixel.
        fig_singlePix: Plot the image with single pixel.
    """

    Img = iio.v3.imread(fid)
    count = []
    imax = Img.max()  # RGB image max     background color
    imin = Img.min()  # RGB image min     fringe color
    a_single = np.double(Img)
    a_single[:, :] = imax
    X = {}
    Y = {}
    Z = {}
    for i in np.arange(Img.shape[0]).reshape(-1):
        countn = 0
        k = 1
        loc = np.asarray(np.where(Img[i, :] == imin))
        for n in np.arange(1, len(loc[0, :])):
            a_single[i, loc[0, 0]] = 0
            X[i, 0] = loc[0, 0]
            Y[i, 0] = i
            Z[i, 0] = (2 * np.pi)
            if loc[0, n] - loc[0, n - 1] > 2:
                a_single[i, loc[0, n]] = 0
                X[i, k] = loc[0, n]
                Y[i, k] = i
                Z[i, k] = (2 * np.pi * (k + 1))
                k = k + 1
        for j in np.arange(Img.shape[1]).reshape(-1):
            temp = Img[i, j]
            if temp == 0:
                countn = countn + 1
        count.append(countn)
    count = np.array(count).reshape(-1)
    sum = np.sum(count) / Img.shape[0]
    width = int(sum / k) + 1
    X = dict2mat(X, 2)
    Y = dict2mat(Y, 2)
    Z = dict2mat(Z, 2)

    fig_singlePix = plt.imshow(a_single, cmap='gray')
    plt.title(title, fontweight='bold', fontsize=14)
    plt.xlabel(r'x [Pixels]', fontsize=14)
    plt.ylabel(r'y [Pixels]', fontsize=14)
    plt.show()

    return Img, width, k, X, Y, Z, a_single, fig_singlePix


def cutoff(Img_single_bk):
    """
    :param Img_single_bk: Image of Background fringe in single pixel.
    :return: Left cutoff and right cutoff.
    """
    loc = np.zeros([Img_single_bk.shape[0], Img_single_bk.shape[1] - np.count_nonzero(Img_single_bk[0, :])])
    for i in np.arange(Img_single_bk.shape[0]):
        loc_temp = np.asarray(np.where(Img_single_bk[i, :] == 0))
        loc[i, :] = loc_temp
    cutoff_L = int(loc[:, 0].max())
    cutoff_R = int(loc[:, -1].min())

    return cutoff_L, cutoff_R


def PhaseShift(Fringe, Background, Img_Single_bk):
    """
    :param Fringe: The fringe X coordinate value with Gas introduce fringe shift.
    :param Background: The fringe X coordinate value with background fringe.
    :param Img_Single_bk: Single pixel background fringe data array.
    :return: phase_shift: Fringe -  Background, or Background - Fringe, it depends on X coordinate value.
        phase_shift_interp: phase_shift interpolation data.
        phase_shift_plot: Plot the phase shift interpolation result.
    """
    phase_shift = Background - Fringe
    delta_ph = {}  # phase shift difference between adjacent fringes
    delta_fr = {}  # background fringe difference
    phase_shift_interp = {}
    temp_x = phase_shift.shape[0]
    row, col = Img_Single_bk.shape
    k = len(Fringe[0, :])
    for i in np.arange(temp_x):
        for kk in np.arange(k - 1):
            delta_ph[i, kk] = phase_shift[i, kk + 1] - phase_shift[i, kk]
            delta_fr[i, kk] = Background[i, kk + 1] - Background[i, kk]
    delta_ph = dict2mat(delta_ph, 2)
    delta_fr = dict2mat(delta_fr, 2)
    phase_step = np.nan_to_num(delta_ph / delta_fr)
    delta_fr_avg = np.mean(np.mean(delta_fr))
    phase_shift_norm = phase_shift / delta_fr_avg * (2 * np.pi)
    phase_step_norm = phase_step / delta_fr_avg * (2 * np.pi)

    for i in np.arange(row).reshape(-1):
        loc = np.asarray(np.where(Img_Single_bk[i, :] == 0))
        k = len(loc[0, :])
        m = 1
        for j in np.arange(col):
            if Img_Single_bk[i, j] == 0:
                phase_shift_interp[i, j] = phase_shift_norm[i, m - 1]
                if m < k:
                    for n in np.arange(1, delta_fr[i, m - 1]):
                        phase_shift_interp[i, j + n] = phase_shift_interp[i, j] + n * phase_step_norm[i, m - 1]
                m = m + 1

    phase_shift_interp = dict2mat(phase_shift_interp, 2)
    cutoff_L, cutoff_R = cutoff(Img_Single_bk)
    phase_shift_plot = plt.pcolor(phase_shift_interp[:, cutoff_L: cutoff_R])
    plt.title('Phase Shift Map', fontweight='bold', fontsize=14)
    plt.xlabel(r'x [Pixels]', fontsize=14)
    plt.ylabel(r'y [Pixels]', fontsize=14)
    plt.show()

    return phase_shift, phase_shift_interp, phase_shift_plot


def find_center(phase_shift, ori_image, fringe_num):
    """
    :param phase_shift: Fringe -  Background, or Background - Fringe, it depends on X coordinate value.
    :param ori_image: Original background fringe image.
    :param fringe_num: Fringes number in original image.
    :return: center: Locate fringe peaks.
    """
    up_signal = {}
    down_signal = {}
    temp_x = phase_shift.shape[0]
    row, col = ori_image.shape
    for j in np.arange(fringe_num):
        for i in np.arange(1, int(row / 2) + 1):
            if np.abs(phase_shift[i, j]) <= 5:
                up_signal[j] = i

        for i in np.arange(int(row / 2), temp_x):
            if np.abs(phase_shift[i, j]) >= 5:
                down_signal[j] = i + 1

    channel_top = dict2mat(up_signal, 1)
    channel_bottom = dict2mat(down_signal, 1)
    center_avg = np.around((channel_top[:] + channel_bottom[:]) / 2)

    ## Average 5 points to locate phaseshift peak
    center_avg2 = {}
    for j in np.arange(phase_shift.shape[1]):

        sum_3row = {}
        ini_num = 0
        for i in np.arange(phase_shift.shape[0]):
            sum_3row[i] = np.sum(phase_shift[i - 2:i + 3, j])
        sum_3row2 = dict2mat(sum_3row, 1)
        if np.amax(sum_3row2) > ini_num:
            ini_num = np.where(sum_3row2 == sum_3row2.max())
        center_avg2[j] = int(np.array(ini_num[-1][-1]))

    center_avg2 = dict2mat(center_avg2, 1)

    for i in np.arange(phase_shift.shape[1]):
        if center_avg2[:, i] > center_avg[:, i]:
            center_avg[:, i] = center_avg2[:, i]

    center = np.amin(center_avg)

    return center


def abel(phi, y, coeff):
    """
    :param phi: Phase
    :param y: Radial distance
    :param coeff: Micrometer/pixels * Contracting factor from image processing
    :return: fgas: Gas density.
            r: Radian.
    """
    n = len(y[0, :])
    r_t = y[0, 0:n]
    ones_v = np.ones(len(y[0, :]) - 1)
    ones_v = np.array([ones_v]).T
    dphi = np.diff(phi)
    y_smooth = {}
    for jj in np.arange(0, n - 1):
        y_smooth[jj] = ((y[0, jj] + y[0, jj + 1]) / 2)
    y_sm_arr = dict2mat(y_smooth, 1)
    dphi_mat = ones_v * dphi
    y_sm_arr_mat = ones_v * y_sm_arr
    r_mat = (ones_v * r_t).T

    coeff_phase = 2 * np.pi
    Cn = 1.16e-23                    # Nitrogen gas refractive index coefficient
    wavelength = 0.532               # Wavelength for Green (um)
    coeff_n = wavelength / Cn / coeff_phase / coeff
    temp = {}
    for i in np.arange(len(y_sm_arr_mat[:, 0])):
        for j in np.arange(len(y_sm_arr_mat[0, :])):
            temp[i, j] = (1 / cmath.sqrt(y_sm_arr_mat[i, j] ** 2 - r_mat[i, j] ** 2))
    temp = dict2mat(temp, 2)[0:j, 1::]
    fgas = 1 / np.pi * (np.matmul(np.triu(dphi_mat[0:j, 1::] * temp), ones_v[0:j])) * coeff_n
    r = r_t * coeff

    return fgas, r


def density(phase_shift_interp, center, step, coeff):
    """
    :param phase_shift_interp: phase_shift map.
    :param center: Fringe peaks location
    :param step: Pick one frame every step (e.g. pick 1 frame out of 10 frames, step=10)
    :param coeff: Micrometer/pixels * Contracting factor from image processing.
    :return: F_r_arr: Assign Smooth density to each frame.
    """
    temp_x = phase_shift_interp.shape[0]
    edge = step                 # cut edge
    peakrange = temp_x / 2 / 3
    end = phase_shift_interp.shape[1]
    st = 0
    frame1 = {}
    frame2 = {}
    frame = {}
    F_smooth = {}
    x = {}
    y = {}
    xz = {}
    F_r = {}
    r0 = {}
    p = 50          # displacement to separate plots

    for z in np.arange(st, end, step).reshape(-1):
        for i in np.arange(temp_x).reshape(-1):
            frame1[i] = np.nan_to_num(np.mean(phase_shift_interp[i, z:z + edge]))
        fr1_arr = dict2mat(frame1, 1)

        for i in np.arange(temp_x - edge + 1).reshape(-1):
            frame2[i] = p - np.mean(fr1_arr[0, i:i + edge])
        fr2_arr = dict2mat(frame2, 1)

        ref = 500
        x1 = int(center - peakrange / 2)
        x2 = int(center + peakrange / 2)

        for i in np.arange(x1, x2):
            s = 0
            for j in np.arange(1, int(peakrange / 2)):
                s = s + np.abs(fr2_arr[0, i + j] - fr2_arr[0, i - j])
            if s < ref:
                ref = s
                r0 = i
        print(r0)
        for i in np.arange(int(center)):
            r1 = int(center)
            if r0 < center:
                r1 = r0
            frame[i] = fr2_arr[0, r1 - i - 1]
        fr_arr = dict2mat(frame, 1)

        for i in np.arange(len(fr_arr[0, :])):
            x[i] = fr_arr[0, i] - np.amin(fr_arr)
            y[i] = i
        x_arr = dict2mat(x, 1)
        y_arr = dict2mat(y, 1)
        xz[(z - st) / step] = x_arr
        xz_arr = np.array(list(xz.values()))[:].reshape(-1, len(x_arr[0, :]))
        F, r = abel(x_arr, y_arr, coeff)
        ini = 0
        for i in np.arange(ini, len(F) - step):
            F_smooth[i] = np.mean(F[np.max([1, int(i - edge)]):int(i + edge)])
        F_sm_arr = dict2mat(F_smooth, 1)
        F_r[(z - st) / step] = F_sm_arr
        F_r_arr = np.array(list(F_r.values()))[:].reshape(-1, len(F_sm_arr[0, :]))

    return F_r_arr


def Axis_plot(phase_shift_interp, step, F_r_arr, coeff):
    """
    :param phase_shift_interp: phase_shift map.
    :param step: Grab one frame every step (e.g. grab 1 frame out 10 frames, step=10)
    :param F_r_arr: Smooth density data.
    :param coeff: Micrometer/pixels * Contracting factor from image processing.
    :return: On-axis plot of phase, density, and real density
    """
    end = phase_shift_interp.shape[1]
    st = 0
    sl = (end - st) / step  # slice
    r_avg = 10                             # smooth window length = 2*r_avg.
    F_r_abs = np.abs(F_r_arr)
    ne = {}
    z_real = {}
    for i in np.arange(1, int(sl) + 1):
        ne[i] = np.mean(F_r_abs[i - 1, 0:r_avg - 1])
        z_real[i] = coeff * (st + step * (i - 1))
    ne_arr = dict2mat(ne, 1)
    z_real = dict2mat(z_real, 1)
    ne_z = {}
    for i in np.arange(int(sl)):
        ne_z[i] = np.mean(ne_arr[0, np.amax([0, i - r_avg]): np.amin([i + r_avg, int(sl) - 1])])
    ne_z = dict2mat(ne_z, 1)
    fig = plt.plot(z_real.T, (ne_z / 1e19).T)
    plt.title('On-axis Density (Smooth)', fontweight='bold', fontsize=14)
    plt.xlabel(r'z [$\mu$m]', fontsize=14)
    plt.ylabel(r'Density [$1e19/cm^3$]', fontsize=14)
    plt.show()

    return fig


def density2(F_r_arr, center, coeff, Img_Single_bk):
    """
    :param F_r_arr: Smooth density data.
    :param center: Fringe peaks location
    :param coeff: Micrometer/pixels * Contracting factor from image processing.
    :param Img_Single_bk: Single pixel background fringe data array
    :return: fig: 2D density distribution
            fig2: Enhanced 2D density distribution
    """
    F_r_arr = np.abs(F_r_arr)
    F_rT_arr = F_r_arr[:, :].T
    f_temp = np.zeros([len(F_rT_arr[:, 0]) * 2, len(F_rT_arr[0, :])])
    for i in np.arange(len(F_r_arr[0, :])):
        f_temp[i] = F_rT_arr[len(F_r_arr[0, :]) - i - 1, :]
        f_temp[i + len(F_r_arr[0, :]) - 1] = F_rT_arr[i, :]

    f_reshape = f_temp
    a, b = f_reshape.shape
    b1 = np.arange(0, b, 0.1)
    a1 = np.arange(0, a)
    fun = interpolate.interp2d(np.arange(b), np.arange(a), f_reshape, kind='linear')
    finalf = fun(b1, a1)
    finalf = savgol_filter(finalf, 10, 2)

    cutoff_L, cutoff_R = cutoff(Img_Single_bk)
    xt, yt = finalf.shape
    fig = plt.pcolor(finalf[:, cutoff_L:cutoff_R])
    yticks11 = (np.arange(xt) - center) * coeff
    yticksMax = yticks11.max() - yticks11.max() % 100
    if yticksMax * 2 / 100 > 7:
        pointsY = int(yticksMax * 2 / 200 + 1)
    if yticksMax * 2 / 100 <= 7:
        pointsY = int(yticksMax * 2 / 100 + 1)
    yticlabel = np.linspace(yticksMax, -yticksMax, pointsY, endpoint=True).astype(int)

    xticks11 = np.arange(cutoff_L, cutoff_R) * coeff / 1000
    xticksMax = xticks11.max()
    xticksMin = xticks11.min()
    if (xticksMax - xticksMin) / 0.2 > 7:
        pointsX = int(xticksMax / 0.5 + 1)
        xticlabel = np.round(np.arange(xticksMin, xticksMax, 0.5), 2)
    if (xticksMax - xticksMin) / 0.2 <= 7:
        pointsX = int(xticksMax / 0.2 + 1)
        xticlabel = np.round(np.arange(xticksMin, xticksMax, 0.2), 2)

    ax = plt.gca()
    ax.set_xticklabels(xticlabel)
    ax.set_yticklabels(yticlabel)
    ax.xaxis.set_major_locator(LinearLocator(pointsX))
    ax.yaxis.set_major_locator(LinearLocator(pointsY))
    plt.title('Density Distribution', fontweight='bold', fontsize=14)
    plt.xlabel(r'z [mm]', fontsize=14)
    plt.ylabel(r'y [$\mu$m]', fontsize=14)
    plt.colorbar(fig)
    plt.show()

    fig2 = plt.pcolor(finalf[:, cutoff_L:cutoff_R].T)
    ax2 = plt.gca()
    ax2.set_yticklabels(xticlabel)
    ax2.set_xticklabels(yticlabel)
    ax2.yaxis.set_major_locator(LinearLocator(pointsX))
    ax2.xaxis.set_major_locator(LinearLocator(pointsY))
    plt.gca().invert_yaxis()
    plt.title('Nozzle Entrance', fontweight='bold', fontsize=14)
    plt.xlabel(r'y [$\mu$m]', fontsize=14)
    plt.ylabel(r'z [mm]', fontsize=14)
    plt.colorbar(fig2)
    plt.clim(0e19, 6e19)
    plt.show()

    return fig, fig2
