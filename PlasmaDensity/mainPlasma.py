from functionPlasma import *

path = "C:/Users/adminRen/PycharmProjects/PlasmaDensity"
fid_bk = os.path.join(path, 'plasma2bk.bmp')
fid = os.path.join(path, 'plasma2.bmp')
cali_fac = 1.3                                      # Micrometer/pixels
pic_conv_fac = 2.2                                  # Contracting factor from image processing
coeff = cali_fac * pic_conv_fac
step = 10                                           # Pick 1 out of 10 pixels to accelerate final density mapping

#########################################################################################
# Extract Background Fringe from fid_bk to get the fringe number k and fringe value Zbk #
#########################################################################################
title1 = 'Vacuum Phase Lines'
Img_bk, width_bk, k, Xbk, Ybk, Zbk, a_single_bk, plot_bk = get_structure(fid_bk, title1)

#################################
# Plasma Introduce Fringe Shift #
#################################
title2 = 'Phase Lines Shift from Plasma'
Img, width, k, X, Y, Z, a_single, plot_p = get_structure(fid, title2)

######################################################
# 2D Phase Mapping from phase differences on fringes #
######################################################

phase_shift, phase_shift_interp, Ph_shift_plot = PhaseShift(X, Xbk, a_single_bk)

################################################
# Locate Fringe Peaks to set up axial symmetry #
################################################

center = find_center(phase_shift, Img, k)

###############################
# Inverse Abel Transformation #
###############################

F_r_arr = density(phase_shift_interp, center, step, coeff)

########################
# On-axis Density Plot #
########################

plot_real_dens = Axis_plot(phase_shift_interp, step, F_r_arr, coeff)


######################
# 2D Density Mapping #
######################
dens_plot = density2(F_r_arr, center, coeff, a_single_bk)