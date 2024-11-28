import sympy as sy
import numpy as np

#MDH
alpha=[0,-np.pi/2,0,-np.pi,-np.pi/2,np.pi/2]
a=[0,0,368,316,0,0]
offset=[-3*np.pi/2,-np.pi/3,-np.pi/2,np.pi/2,np.pi/2,0]
d=[122,140.5,0,0,-102.5,94]


"""推导转换矩阵"""
#z轴旋转角度
theta1=sy.Symbol("theta1")
#z轴平移
d1=sy.Symbol("d1")

#x轴旋转角度
al1=sy.Symbol("al1")
#x轴平移

a1=sy.Symbol("a1")

def MTi(alpha,a,d,theta,DHtype="MDH"):

    #x轴旋转变化阵，右手系，逆时针为正
    Rx=sy.Matrix([
        [1, 0,        0,        0],
        [0, sy.cos(alpha), -sy.sin(alpha), 0],
        [0, sy.sin(alpha), sy.cos(alpha),  0],
        [0, 0,        0,         1]
    ])
    #沿着x轴平移，有方向
    Tx=sy.Matrix([
        [1, 0, 0,  a],
        [0, 1, 0,  0],
        [0, 0,1,   0],
        [0, 0,0,   1]
    ])
    #z轴旋转阵，右手系，逆时针为正
    Rz=sy.Matrix([
        [sy.cos(theta),  -sy.sin(theta),   0,  0],
        [sy.sin(theta), sy.cos(theta),   0,  0],
        [0,             0,            1,  0],
        [0,             0,            0,  1]
    ])
    #沿着z轴平移，有方向
    Tz=sy.Matrix([
        [1, 0, 0,  0],
        [0, 1, 0,  0],
        [0, 0,1,   d],
        [0, 0,0,   1]
    ])

    print("MDH——TR",Rx@Tx@Rz@Tz)
    #右乘原则，先旋转X轴alpha，再沿X轴平移a，再旋转Z轴theta（这里为offset值），最后沿着Z轴平移D
    return Rx@Tx@Rz@Tz
#MDH的转换矩阵如下：
MTi(al1,a1,d1,theta1)