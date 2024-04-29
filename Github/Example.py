# coding: utf-8

#<editor-fold desc="库引用声明">
import matplotlib as mpl
mpl.use('TkAgg')
from scipy.integrate import trapz
import os
from lidar_params import *
# </editor-fold>
'''数据保存路径'''
current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, 'data')

Is_sea_fog = 1
c = 3e8
band = ['S', 'C', 'X', 'Ka', 'W']
freq = [2.885e9, 5.43e9, 9.38e9, 34.88e9, 93.75e9] #Hz
wave_length = [1e3*c/freq_0 for freq_0 in freq]    #mm

def threshold(n_wave, R_max, ele, wind, visibility, rainrate, S_min_dBm):

    # 求解bsca和ext的数值范围
    # <editor-fold desc="晴空，与仰角有关，只能现场计算,已与Pr_in_aerosol对比验证bsca和ext计算结果">
    R_set = np.arange(10, R_max, 20)#m
    Wind_set = np.arange(wind[0], wind[1], 2)#m/s
    concerned_height = R_set * np.tan(ele)
    rd = np.arange(0.01, 50, 0.1)  # rd 0.04, 12.59, 0.01
    bsca_aerosol = np.zeros((len(Wind_set), len(R_set))) #不同风速，不同高度上的消光截面
    ext_aerosol = np.zeros((len(Wind_set), len(R_set)))
    scattering_aerosol = np.zeros((len(Wind_set), len(R_set)))

    data = np.load(data_path+'/'+'cross_section_aerosol in ' + band[n_wave] + '.npz')
    bsca = data['bsca']
    ext = data['ext']
    U = data['U']
    Z = data['Z']
    for i in range(0, len(Wind_set)):
        for j in range(0, len(concerned_height)):
            index_U = np.where(U == min(U, key=lambda x: abs(x - Wind_set[i])))[0][0]
            index_Z = np.where(Z == min(Z, key=lambda x: abs(x - concerned_height[j])))[0][0]
            bsca_aerosol[i, j] = bsca[index_U, index_Z]
            ext_aerosol[i, j] = ext[index_U, index_Z]
            scattering_aerosol[i, j] = bsca_aerosol[i, j] * 10 ** ( -0.2 * trapz(ext_aerosol[i, 0:j + 1], R_set[0:j + 1]))
    # < / editor - fold >

    #<editor-fold desc="雾">
    data = np.load(data_path + '/' + 'IsSeaFog=' + str(Is_sea_fog) + ' in ' + band[n_wave] + ' band.npz')
    sigma_bsca = data['sigma_bsca']
    sigma_ext = data['sigma_ext']
    data = np.load(data_path+'/' + 'Visibility.npz')
    Visibility = data['Visibility']
    # 创建一个布尔数组，表示 B 中的每个元素是否在阈值范围内
    mask = (Visibility >= visibility[0]) & (Visibility <= visibility[1])
    # 使用布尔数组来索引 A 和 B 中的元素
    bsca_fog = sigma_bsca[mask]
    ext_fog = sigma_ext[mask]
    # </editor-fold>

    #<editor-fold desc="雨">
    data = np.load(data_path + '/' + band[n_wave] + ' band-rain-RainRatePSD.npz')
    sigma_bsca = data['sigma_bsca']
    sigma_ext = data['sigma_ext']
    data = np.load(data_path + '/' + 'RainRate.npz')
    RainRate = data['RainRate']
    # 创建一个布尔数组，表示 B 中的每个元素是否在阈值范围内
    mask = (RainRate >= rainrate[0]) & (RainRate <= rainrate[1])
    # 使用布尔数组来索引 A 和 B 中的元素
    bsca_rain = sigma_bsca[mask]
    ext_rain = sigma_ext[mask]
    # </editor-fold>

    bsca_precipitation = np.concatenate((bsca_fog, bsca_rain), axis=0)
    ext_precipitation = np.concatenate((ext_fog, ext_rain), axis=0)
    scattering_precipitation = bsca_precipitation * 10**(-0.2*R_max*ext_precipitation)

    scattering_all = np.concatenate((scattering_aerosol.flatten(), scattering_precipitation), axis=0)
    scattering_min = np.min(scattering_all)
    min_level = 10**(S_min_dBm/10)*1e-3 * 2**10 * np.log(2) * np.pi**2 * R_max**2 / ((wave_length[n_wave]*1e-3)**2 * c * scattering_min)

    return min_level

if __name__ == '__main__':
    R_max = 500 # m, ranging from 1 to 5e3
    wind = [20, 25] # m/s, ranging from 0.5 to 30
    visibility = [1, 5] # km, ranging from 0.05 to 10
    rainrate = [1, 15] # mm/h, ranging from 0.1 to 50
    ele = np.deg2rad(3.5) #rad
    n_wave = 4 # ranging from 0 to 4 for S, C, X, Ka, W
    S_min_dBm = -104 #dBm
    min_level = threshold(n_wave, R_max, ele, wind, visibility, rainrate, S_min_dBm)
    print('The lower limit is ', str(min_level))
    S=1