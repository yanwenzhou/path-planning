"""
Path planning with Bezier curve.Atsushi Sakai
yanwenzhou
"""
import sys
print(sys.executable)
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import matplotlib.image as mpimg
import matplotlib.patches as patches
import math
from scipy import interpolate
from scipy.interpolate import CubicSpline
from scipy import integrate
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.interpolate import interp1d
import cv2
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontManager
import matplotlib.cm as cm
import matplotlib.colors as mcolors
# font = FontProperties(fname=r'D:\1111\times.ttf')
# 添加字体文件到matplotlib的字体库
fm = FontManager()
fm.addfont('D:\\1111\\times.ttf')
# 设置全局字体为"Times New Roman"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
# 设置全局的作图文字为宋体
# plt.rcParams['font.sans-serif'] = ['SimHei']
# 设置全局的数字和字母为 Times New Roman
plt.rcParams['font.serif'] = ['times.ttf']
# 如果需要显示负号，需要设置
plt.rcParams['axes.unicode_minus'] = False
show_animation = True                           #是否显示动画
plt.close("all")                                #关闭所有图形窗口


def calc_bezier_path(control_points, n_points=100):                                 #计算贝塞尔曲线路径
    """
    Compute bezier path (trajectory) given control points.                          #给定控制点，计算贝塞尔路径（轨迹）
    :param control_points: (numpy array)                                        #控制点
    :param n_points: (int) number of points in the trajectory                   #轨迹中的点数
    :return: (numpy array)                                                      #返回值是numpy数组
    """
    traj = []                                                           #轨迹
    for t in np.linspace(0, 1, n_points):         #在0到1之间均匀取n_points个点,这里的t是一个参数，用来控制贝塞尔曲线的形状
        traj.append(bezier(t, control_points))                          #将每个点添加到轨迹中
    return np.array(traj)                                               #返回轨迹


def bernstein_poly(n, i, t):                 #伯恩斯坦多项式，是用来计算贝塞尔曲线的
    """         
    Bernstein polynomial.          #伯恩斯坦多项式
    :param n: (int) polynom degree          #多项式的次数
    :param i: (int)                         #整数
    :param t: (float)                       #浮点数
    :return: (float)                    #返回值是浮点数
    """
    return scipy.special.comb(n, i) * t ** i * (1 - t) ** (n - i)   #返回值是组合数乘以t的i次方乘以（1-t）的n-i次方


def bezier(t, control_points):                  #贝塞尔曲线
    """
    Return one point on the bezier curve.       #返回贝塞尔曲线上的一个点
    :param t: (float) number in [0, 1]          #t是一个参数，用来控制贝塞尔曲线的形状
    :param control_points: (numpy array)        #控制点
    :return: (numpy array) Coordinates of the point                     #返回值是一个numpy数组，表示点的坐标
    """
    n = len(control_points) - 1                     #控制点的个数，n是控制点的个数减1
    return np.sum([bernstein_poly(n, i, t) * control_points[i] for i in range(n + 1)], axis=0)      #返回值是一个numpy数组，表示点的坐标


def bezier_derivatives_control_points(control_points, n_derivatives):       #计算贝塞尔曲线的导数
    """
    Compute control points of the successive derivatives of a given bezier curve.   #计算给定贝塞尔曲线的导数的控制点
    A derivative of a bezier curve is a bezier curve.
    See https://pomax.github.io/bezierinfo/#derivatives
    for detailed explanations
    :param control_points: (numpy array)
    :param n_derivatives: (int)
    e.g., n_derivatives=2 -> compute control points for first and second derivatives
    :return: ([numpy array])
    """
    w = {0: control_points}             #w是一个字典，key是0，value是控制点
    for i in range(n_derivatives):      #计算n_derivatives次导数
        n = len(w[i])                   #控制点的个数, n是w[i]的长度,w[i]是一个numpy数组
        w[i + 1] = np.array([(n - 1) * (w[i][j + 1] - w[i][j])          
                            for j in range(n - 1)])
    return w


def curvature(dx, dy, ddx, ddy):
    """
    Compute curvature at one point given first and second derivatives.      #计算给定点处的曲率，给定一阶和二阶导数

    :param dx: (float) First derivative along x axis
    :param dy: (float)
    :param ddx: (float) Second derivative along x axis
    :param ddy: (float)
    :return: (float)
    """
    return (dx * ddy - dy * ddx) / (dx ** 2 + dy ** 2) ** (3 / 2)


def plot_arrow(x, y, yaw, length=1.0, width=0.5, fc="r", ec="k",ax =None):  # pragma: no cover   #绘制箭头，length=1.0指的是箭头的长度，width=0.5指的是箭头的宽度，fc="r"指的是箭头的颜色，ec="k"指的是箭头的边缘颜色
    if not isinstance(x, float):                                                        #如果x不是浮点数
        for (ix, iy, iyaw) in zip(x, y, yaw):                                           #遍历x,y,yaw
            plot_arrow(ix, iy, iyaw)                                                    #绘制箭头
    else:
        plt.arrow(x, y, length * np.cos(yaw), length * np.sin(yaw),                     #绘制箭头
                fc=fc, ec=ec, head_width=width*4, head_length=width*4)                #箭头的颜色，边缘颜色，头的宽度，头的长度
        plt.plot(x, y)
    if ax is None:
        ax = plt.gca()

# def create_berth(ax):

#     # 定义泊位的顶点
#     berth_vertices = [(10, 10), (30, 20), (10, 5), (30, 5)]

#     # 创建一个多边形
#     berth = patches.Polygon(berth_vertices, closed=True, fill=True, color='black')

#     # 将多边形添加到图形中
#     ax.add_patch(berth)

#     # # 显示图形
#     # plt.show()

wpt_t = []
wpt_t1 = []
wpt_t2 = []
wpt_t3 = []
scenarios = [
        {"start_x": 0.0, "start_y": 0.0, "start_yaw": np.radians(30.0), "wpt_t": wpt_t,"n1":0.3,"n2":0.3,"n3":1.6,"n4":1.6},  # 场景1
        {"start_x": 200.0, "start_y": 0.0, "start_yaw": np.radians(110.0), "wpt_t": wpt_t1,"n1":0.9,"n2":0.9,"n3":0.9,"n4":0.9},  # 场景2
        {"start_x": 200, "start_y": 500.0, "start_yaw": np.radians(-10.0), "wpt_t": wpt_t2,"n1":1.2,"n2":13,"n3":1,"n4":1},  # 场景3
        {"start_x": 200.0, "start_y": 500.0, "start_yaw": np.radians(-65.0), "wpt_t": wpt_t3,"n1":2.2,"n2":2.6,"n3":1,"n4":1},  # 场景4
    ]



# 使用场景1的参数
scenario = scenarios[3]  # 索引为0表示第一个场景




#贝塞尔曲线是一种参数化的曲线，它是由一系列控制点决定的，这些控制点决定了曲线的形状。calc_4points_bezier_path函数计算了四个控制点的贝塞尔曲线路径。
def calc_4points_bezier_path(sx, sy, syaw, ex, ey, eyaw, offset):               #这里的意思是   通过起始点和终点的位置和方向，计算出控制点和路径
    """
    Compute control points and path given start and end position.               #计算控制点和路径
    :param sx: (float) x-coordinate of the starting point                           #起始点的x坐标
    :param sy: (float) y-coordinate of the starting point                       #起始点的y坐标
    :param syaw: (float) yaw angle at start                                     #起始点的偏航角
    :param ex: (float) x-coordinate of the ending point                         #终点的x坐标    
    :param ey: (float) y-coordinate of the ending point                         #终点的y坐标
    :param eyaw: (float) yaw angle at the end                                   #终点的偏航角
    :param offset: (float)                                                      #偏移量
    :return: (numpy array, numpy array)                                         #返回值是两个numpy数组
    """
    dist = np.hypot(sx - ex, sy - ey) / offset                                  #计算两点之间的距离
    control_points = np.array(                                                      #控制点
        [[sx, sy],                                                                  #起始点
           [sx + scenario["n1"]*dist * np.cos(syaw), sy + scenario["n2"]*dist * np.sin(syaw)],                          #情况一
           [ex - scenario["n3"]*dist * np.cos(eyaw), ey - scenario["n4"]*dist * np.sin(eyaw)], 
        #  [sx + 0.9*dist * np.cos(syaw), sy + 0.9*dist * np.sin(syaw)],                              #情况二
        #  [ex - 0.9*dist * np.cos(eyaw), ey - 0.9*dist * np.sin(eyaw)],   
        #    [sx + 0.6*dist * np.cos(syaw), sy + 0.6*dist * np.sin(syaw)],                          #情况三                    
        #    [sx + 1.2*dist * np.cos(syaw), sy + 13*dist * np.sin(syaw)],                          
        # #    [ex - 1.0*dist * np.cos(eyaw), ey - 1.0*dist * np.sin(eyaw)], 
            # [sx + 1.2*dist * np.cos(syaw), sy + 1.6*dist * np.sin(syaw)],                          #情况四                  
            # [sx + 2.2*dist * np.cos(syaw), sy + 2.6*dist * np.sin(syaw)],                          
            # [ex - 1.0*dist * np.cos(eyaw), ey - 1.0*dist * np.sin(eyaw)], 
        [ex, ey]])                                                                 #终点
    path = calc_bezier_path(control_points, n_points=100)                           #计算贝塞尔曲线路径
    return path, control_points                                                     #返回路径和控制点

def main(): 
    # 定义一个列表，保存四个场景的参数
    # 定义一个列表，保存四个场景的参数

    start_x = scenario["start_x"]
    start_y = scenario["start_y"]
    start_yaw = scenario["start_yaw"]

    m = 2                                                            #主函数
    """Plot an example bezier curve."""
    # start_x = 10.0  # [m]                                               #起始点的x坐标            情况四
    # start_y = 90.0 # [m]                                                #起始点的y坐标
    # start_yaw = np.radians(-65.0)  # [rad]                               #起始点的偏航角，radians()函数将角度转换为弧度

    # start_x = -10.0  # [m]                                               #起始点的x坐标            情况三
    # start_y = 80.0 # [m]                                                #起始点的y坐标
    # start_yaw = np.radians(-10.0)  # [rad]                               #起始点的偏航角，radians()函数将角度转换为弧度

    # start_x = 00.0  # [m]                                               #起始点的x坐标            情况二
    # start_y = 00.0 # [m]                                                #起始点的y坐标
    # start_yaw = np.radians(110.0)  # [rad] 

    # start_x = 0.0  # [m]                                               #起始点的x坐标            情况一
    # start_y = 0.0 # [m]                                                #起始点的y坐标
    # start_yaw = np.radians(30.0)  # [rad] 

    end_x =  300.0  # [m]                                             #终点的x坐标    
    end_y = 400.0  # [m]                                                 #终点的y坐标
    end_yaw = np.radians(110.0)  # [rad]                                #终点的偏航角
    # end_yaw = np.radians(30.0)
    offset = 2                                                           #偏移量
    # 切线的长度
    length = 10

    # 计算切线的结束点
    tangent_end_x = end_x + length * np.cos(end_yaw)
    tangent_end_y = end_y + length * np.sin(end_yaw)
    path, control_points = calc_4points_bezier_path(
        start_x, start_y, start_yaw, end_x, end_y, end_yaw, offset)
    min_radius = float('inf')
    # Display the tangent, normal and radius of curvature at the point with the minimum radius          #在曲率最小的点处显示切线、法线和曲率半径
    derivatives_cp = bezier_derivatives_control_points(control_points, 2)       #计算贝塞尔曲线的导数 这里的2指的是计算二阶导数   

    # 遍历曲线上的所有点
    t1 = 0.5
    vs_array = []  # 创建空数组
    ud = np.linspace(0.8, 0, num=100)
    ud_values = []  # 创建空列表来存储ud的值
    wpt = []
    a=1
    v_start = 1.6
    for t in np.linspace(0, a, 100):
        # if t < v_start/4:
        #     ud = v_start - t * (3/4 * v_start) / (a/5)
        # else:
        #     ud = (1/4*v_start)/(4*a/5) - t * (1/4*v_start)/(4*a/5)
        
        # if t < 5*v_start/8:
        #     ud = v_start - t * (3/4 * v_start) / (a/2)
        # else:
        #     ud = (v_start/4)/(1*a/2) - t * (v_start/4)/(1*a/2)

        ud=v_start-v_start*t

        # ud = 0.008/(t+0.01)
        # ud = 0.8/np.log(t+2)
        # ud = 0.8 * np.exp(-0.01 * t)
        # ud = 0.008/np.exp(t/100)

        ud_values.append(ud)  # 将ud的值添加到列表中
        point = bezier(t, control_points)
        wpt.append([point[1], point[0]])       #将点添加到wpt中
        dt = bezier(t, derivatives_cp[1])
        dt_length = np.linalg.norm(dt)          #计算矢量的长度
        # Vs = ud / np.sqrt(dt_length) ** (1/4)             #计算速度
        Vs = ud / dt_length            #计算速度
        # print(dt_length )
        vs_array.append(Vs)  # Append Vs value to the array
        ddt = bezier(t, derivatives_cp[2])
        radius = 1 / curvature(dt[0], dt[1], ddt[0], ddt[1])
        # Update the minimum radius and corresponding parameter value
        if not np.isinf(radius) and radius < min_radius:
            min_radius = radius
            t1 = t
    vs_scaled = (vs_array - np.min(vs_array)) / (np.max(vs_array) - np.min(vs_array)) * (np.max(ud_values) - np.min(ud_values)) + np.min(ud_values)
    for i in range(100):
        speed = vs_scaled[i]
        wpt[i] = list(wpt[i])
        wpt[i].append(speed)
    # print(vs_array)
    vs_scaled_smooth = savgol_filter(vs_scaled, 49, 3)          #平滑处理
    # plt.plot(np.linspace(0, 1, 100), vs_scaled_smooth)               #绘制速度，linspace(0,1,100)的意思是在0到1之间生成100个点
    # # plt.plot(np.linspace(0, 1, 100), vs_scaled)            
    # plt.xlabel('x')
    # plt.ylabel('Vs')
    # plt.title('trajectory speed')
    # plt.show()
    # # print(wpt)
    # np.savetxt('wpt.txt', wpt)
    # 假设 x 是你的 x 轴数据
    x = np.linspace(0, 10, len(vs_scaled_smooth))
    # 计算面积
    area = integrate.trapz(vs_scaled_smooth, x)
    # 打印结果
    print("area of vs*x(即走过的距离，没有缩放为物理上的时间，数值较小)=",area)


    # plt.title('speed assignment')

    print("min radius:", min_radius)  #输出最小半径
    point = bezier(t1, control_points)                                           #计算贝塞尔曲线上的一个点
    dt = bezier(t1, derivatives_cp[1])                                          #计算一阶导数
    print("tangent:", dt)                                                  
    ddt = bezier(t1, derivatives_cp[2])                                          #计算二阶导数
    radius = 1 / curvature(dt[0], dt[1], ddt[0], ddt[1])                #计算曲率半径
    dt /= np.linalg.norm(dt, 2)                             #计算单位切线
    # tangent = np.array([point, point + dt])             #计算切线
    # normal = np.array([point, point + [- dt[1], dt[0]]])        #计算法线
    # curvature_center = point + np.array([- dt[1], dt[0]]) * radius      #计算曲率中心 
    # circle = plt.Circle(tuple(curvature_center), radius,        #绘制圆
                        # color=(0, 0.8, 0.8), fill=False, linewidth=1)       #圆的颜色，填充，线宽

    assert path.T[0][0] == start_x, "path is invalid"                   #断言
    assert path.T[1][0] == start_y, "path is invalid"                   
    assert path.T[0][-1] == end_x, "path is invalid"
    assert path.T[1][-1] == end_y, "path is invalid"


    # 计算路径点,计算总的路程
    path_points = np.array([path.T[0], path.T[1]]).T
    # 计算相邻路径点之间的距离
    distances = np.sqrt(np.sum(np.diff(path_points, axis=0)**2, axis=1))
    # 检查和处理零值
    vs_scaled_no_zero = vs_scaled[1:][vs_scaled[1:] != 0]
    # 获取对应的距离值
    distances_no_zero = distances[vs_scaled[1:] != 0]
    # 计算每个路径点处的物理时间
    times = distances_no_zero / vs_scaled_no_zero
    # 计算总的物理时间  
    total_time = np.sum(times)
    # 计算总路程
    total_distance = np.sum(distances)
    # 计算累积时间
    cumulative_times = np.cumsum(times)
    # 获取 cumulative_times 的最后一个值
    last_value = cumulative_times[-1]
    # 将最后一个值添加到 cumulative_times 两次
    cumulative_times = np.append(cumulative_times, [last_value, last_value])
    print("Total distance:", total_distance)
    # print("path_points:", path_points)
    # print(wpt)
    # print(vs_scaled)
    print("Total time:", total_time)            #这个不准
    # print("Cumulative times:", cumulative_times)
    # print("Length of wpt:", len(wpt))
    # print("Length of Cumulative time:", len(cumulative_times))
    # print("Length of vs_scaled:",len(vs_scaled))

    #  最新速度分配方法  （改进后）
    # 计算当前面积
    current_area = integrate.trapz(vs_scaled_smooth, np.linspace(0, 1, len(vs_scaled_smooth)))
    # 计算缩放比例
    scale_factor = total_distance / current_area
    # 缩放横坐标
    x_scaled = np.linspace(0, 1, len(vs_scaled_smooth)) * scale_factor
    # 绘制缩放后的曲线
    plt.plot(x_scaled, vs_scaled_smooth)
    plt.xlabel('t')
    plt.ylabel('Vs')
    plt.title('trajectory speed')
    plt.show()

    #把数据添加到wpt_t中
    # 创建插值函数
    f = interp1d(x_scaled, vs_scaled_smooth)
    # 创建一个空列表来保存结果
    vs_t = []
    # 遍历每一秒
    for t in range(int(x_scaled[-1]) + 1):
        # 计算对应的速度
        v = f(t)
        # 将结果保存到 vs_t 中
        vs_t.append(v)
    # 打印结果
    # print(vs_t)
    # 计算每秒走过的路程
    distances = [v * 1 for v in vs_t]  # 假设每个步长对应 1 秒
    # 计算总路程
    total_distances = np.cumsum(distances)
    # 打印结果
    print(total_distances)
    # 计算路径的长度
    path_distances = [0]
    for i in range(1, len(path)):
        path_distances.append(path_distances[-1] + np.sqrt((path[i][0] - path[i - 1][0]) ** 2 + (path[i][1] - path[i - 1][1]) ** 2))
    # 创建插值函数
    fx = interp1d(path_distances, [p[0] for p in path], bounds_error=False, fill_value='extrapolate')
    fy = interp1d(path_distances, [p[1] for p in path], bounds_error=False, fill_value='extrapolate')
    # 创建一个空列表来保存结果
    points_t = []
    # 遍历每一秒
    for total_distance in total_distances:
        # 计算对应的点
        x = fx(total_distance)
        y = fy(total_distance)
        # 将结果保存到 points_t 中
        points_t.append([x, y])
    # 打印结果
    # print('points_t',points_t)
    # print("Length of points_t:", len(points_t))
    # print("Length of total_distances:", len(total_distances))
    # print("Length of vs_t:",len(vs_t))
    # 创建一个空列表来保存结果
    # 遍历 points_t 和 vs_t
    for t, (point, v) in enumerate(zip(points_t, vs_t)):
        # 将 y 坐标、x 坐标、t 和 v 添加到 wpt_t1
        # wpt_t11.append([point[1], point[0], t, v])
        scenario["wpt_t"].append([point[1], point[0], t, v])
    # 创建一个空列表来保存结果
    tangent_directions = []
    # 遍历每个点
    for i in range(1, len(points_t)):
        # 计算 dx 和 dy
        dx = points_t[i][0] - points_t[i - 1][0]
        dy = points_t[i][1] - points_t[i - 1][1]
        # 计算弧度
        radian = np.arctan2(dy, dx)
        # 将结果保存到 tangent_directions 中
        tangent_directions.append(radian)
    # 添加初始点的切线方向
    tangent_directions.insert(0, tangent_directions[0])

    for i in range(len(scenario["wpt_t"])):
        scenario["wpt_t"][i].append(tangent_directions[i])
    # 打印结果
    print('wpt_t',scenario["wpt_t"])
    # 根据 scenario 的值来决定保存到哪个文件
    if scenario == scenarios[0]:
        filename = 'changjing1.txt'
    elif scenario == scenarios[1]:
        filename = 'changjing2.txt'
    elif scenario == scenarios[2]:
        filename = 'changjing3.txt'
    elif scenario == scenarios[3]:
        filename = 'changjing4.txt'
    # 保存文件
    np.savetxt(filename,scenario["wpt_t"] )

    #     # 定义起点和终点
    # start = np.array([61, 60])
    # end = np.array([78, 66])
    # # 计算距离
    # distance = np.linalg.norm(end - start)
    # # 定义速度
    # speed = 0.04  # m/s
    # # 计算所需时间
    # time = distance / speed
    # # 创建时间点
    # t = np.linspace(0, time, int(time) + 1)
    # # 创建插值函数
    # f_x = interpolate.interp1d([0, time], [start[0], end[0]])
    # f_y = interpolate.interp1d([0, time], [start[1], end[1]])
    # # 计算每一秒走到哪个轨迹点
    # points_second = np.array([f_x(t), f_y(t)]).T
    # # 定义恒定的速度和psi值
    # v = 0.04
    # psi = 110 * np.pi / 180
    # # 创建一个新的列表
    # path_second = []
    # # 遍历每个轨迹点和对应的时间点
    # for i, point in enumerate(points_second):
    #     # 添加轨迹点、时间点和恒定的速度、psi值到新的列表中
    #     path_second.append([point[1], point[0], i, v, psi])
    # # print(path_second)
    # path_second = np.array(path_second)
    # wpt_t = np.array(wpt_t)                                        #将wpt_t转换为NumPy数组
    # # 将path_second中的第三列统一加上wpt_t最后一行第三列的值
    # path_second[:, 2] += wpt_t[-1, 2]
    # # 创建一个新的变量来存储合并后的结果
    # wpt_t_new = np.concatenate((wpt_t, path_second), axis=0)   
    # np.savetxt('wpt_t_new.txt',wpt_t_new) 
    # print('加上第二阶段的轨迹',wpt_t_new)

    if show_animation:  # pragma: no cover
        img = mpimg.imread('D://1111//inference//environment_port.png')

        # 创建一个新的图形，包含两个子图
        fig, ax = plt.subplots(1, 2, figsize=(12,4.5))

        # 设置第一个子图的背景图像
        ax[0].imshow(img, extent=[0, 650, 0, 500])  # 显示图片
        ax[0].plot(path.T[0], path.T[1], color='red', linewidth=1.4, label="Bezier Path")  # 绘制贝塞尔路径，设置颜色为红色，线宽为2
        # 在第一个子图中绘制箭头
        plot_arrow([start_x], [start_y], [start_yaw], ax=ax[0], length=6)  # 绘制箭头
        plot_arrow([end_x], [end_y], [end_yaw], ax=ax[0], length=6)

        # 设置第二个子图的背景图像
        ax[1].imshow(img, extent=[0, 600, 0, 400])  # 显示图片

        # 将速度分配可视化
        cmap = mcolors.LinearSegmentedColormap.from_list("mycmap", ["red", "green"])
        colors = [cmap(i) for i in np.linspace(0, 1, 10)]
        wpt_t = np.loadtxt(filename)
        segments_wpt_t = np.array_split(wpt_t, 10)
        segments_path_T_0 = np.array_split(path.T[0], 10)
        segments_path_T_1 = np.array_split(path.T[1], 10)
        for i, (segment_wpt_t, segment_path_T_0, segment_path_T_1) in enumerate(zip(segments_wpt_t, segments_path_T_0, segments_path_T_1)):
            color = colors[i % len(colors)]
            ax[1].plot(segment_path_T_0, segment_path_T_1, color=color, label="{:.2f} m/s".format(abs(segment_wpt_t[-1][3])))
            ax[1].plot(np.nan, np.nan)
            ax[1].legend(fontsize=12, borderaxespad=0.1)

        # ax.plot(control_points.T[0], control_points.T[1],       #绘制控制点
        #         '--o', label="Control Points")      
        fig, ax = plt.subplots(2)
        ax[0].plot(control_points.T[0], control_points.T[1],'--o', label="Control Points")            
        # print(path.T[0], path.T[1])
        # ax.plot([end_x + 2, tangent_end_x +2], [end_y, tangent_end_y],linewidth=2, color='black')
        # ax.plot(tangent[:, 0], tangent[:, 1], label="Tangent")  #绘制切线
        # ax.plot(normal[:, 0], normal[:, 1], label="Normal")     #绘制法线
        # ax.add_artist(circle)                                   #添加圆
        # # 在起始点的位置添加一个代表小船的矩形
        # boat_width = 1
        # boat_height = 0.4
        # start_boat = Rectangle((start_x, start_y), boat_width, boat_height, angle=start_yaw)
        # ax.add_patch(start_boat)
        plot_arrow(start_x, start_y, start_yaw,length=6)                 #绘制箭头
        plot_arrow(end_x, end_y, end_yaw,length=6)                       #绘制箭头

        # 在第一个子图上添加图例
        ax[0].legend()
        # 在第二个子图上添加图例
        ax[1].legend()

        # 设置第一个子图的坐标轴比例
        ax[0].axis("equal")
        # 设置第二个子图的坐标轴比例
        ax[1].axis("equal")

        # 在第一个子图上显示网格
        ax[0].grid(True)
        # 在第二个子图上显示网格
        ax[1].grid(True)

        # 给第一个子图的X轴添加标签
        ax[0].set_xlabel('Y-East / m')
        # 给第二个子图的X轴添加标签
        ax[1].set_xlabel('Y-East / m')

        # 给第一个子图的Y轴添加标签
        ax[0].set_ylabel('X-North / m')
        # 给第二个子图的Y轴添加标签
        ax[1].set_ylabel('X-North / m')

        # 设置第一个子图的比例
        ax[0].set_aspect('equal')
        # 设置第二个子图的比例
        ax[1].set_aspect('equal')
        plt.show()                                 #显示图形


def main2():                                #主函数2
    """Show the effect of the offset."""    #显示偏移的效果
    start_x = 10.0  # [m]                   #起始点的x坐标
    start_y = 1.0  # [m]                    #起始点的y坐标                      
    start_yaw = np.radians(180.0)  # [rad]  #起始点的偏航角

    end_x = -0.0  # [m]                             #终点的x坐标
    end_y = -3.0  # [m]                             #终点的y坐标
    end_yaw = np.radians(-45.0)  # [rad]            #终点的偏航角

    for offset in np.arange(1.0, 5.0, 1.0):         #遍历1.0到5.0之间的数，步长为1.0
        path, control_points = calc_4points_bezier_path(        #计算四个控制点的贝塞尔曲线路径
            start_x, start_y, start_yaw, end_x, end_y, end_yaw, offset)  #通过起始点和终点的位置和方向，计算出控制点和路径
        assert path.T[0][0] == start_x, "path is invalid"       #断言,这里的意思是如果path.T[0][0]不等于start_x,则抛出异常，path.T[0][0]是路径的x坐标，T[0][0】表示第一个点的x坐标
        assert path.T[1][0] == start_y, "path is invalid"       #assert的意思是断言，如果后面的条件为真，则继续执行，否则抛出异常
        assert path.T[0][-1] == end_x, "path is invalid"
        assert path.T[1][-1] == end_y, "path is invalid"

        if show_animation:  # pragma: no cover              
            plt.plot(path.T[0], path.T[1], label="Offset=" + str(offset))   #绘制从起始点到终点的路径，label="Offset=" + str(offset)表示偏移量，str(offset)表示将offset转换为字符串

    if show_animation:  # pragma: no cover            #如果show_animation为真
        plot_arrow(start_x, start_y, start_yaw)
        plot_arrow(end_x, end_y, end_yaw)
        plt.legend()
        plt.axis("equal")
        plt.grid(True)
        plt.show()


if __name__ == '__main__':
    main()
    # main2()

