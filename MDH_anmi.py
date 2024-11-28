#MDH方法
import sympy as sy
import numpy as np
import math
from math import atan2

#MDH参数
alpha=[0,-np.pi/2,0,-np.pi,-np.pi/2,np.pi/2]
a=[0,0,368,316,0,0]
offset=[-3*np.pi/2,-np.pi/3,-np.pi/2,np.pi/2,np.pi/2,0]
d=[122,140.5,0,0,-102.5,94]

def MTi(alpha,a,d,theta,DHtype="MDH"):

    #x轴变化阵，右手系，顺针旋负方向
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
    #z轴旋转阵，右手系，顺针旋负方向
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

    #print("MDH——TR",Rx@Tx@Rz@Tz)
    return Rx@Tx@Rz@Tz

def Aiv(alpha,a,d,offset,thetas=[0,0,0,0,0,0],Tnum=6):
    """
    依次右乘，推导出机械臂位姿矩阵
    """
    print("theta:",thetas)

    for i in range(Tnum):
        print("tnumi",offset[i]+thetas[i])
        Ti=MTi(alpha[i],a[i],d[i],offset[i]+thetas[i])
        if i==0:
            Ai=Ti
        else:
            Ai=Ai@Ti
    return Ai
#机械臂位姿矩阵如下
A06= Aiv(alpha,a,d,offset)
#print("位姿矩阵：",A06)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm,animation
def DoARM(alpha,a,d,offset,thetas=[0,0,0,0,0,0]):
    X=[0]
    Y=[0]
    Z=[0]
    for i in range(6):
        A0i=Aiv(alpha,a,d,offset,thetas,Tnum=i+1)
        #print("A",i,":",A0i)
        X.append(A0i[:3, 3][0])
        Y.append(A0i[:3, 3][1])
        Z.append(A0i[:3, 3][2])
    return X,Y,Z
# Plot the arm segments
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X,Y,Z=DoARM(alpha,a,d,offset)
an=ax.scatter(X,Y,Z, c=['r', 'g', 'b', 'c', 'm', 'y','g'], s=50)
links,=ax.plot(X, Y,Z,'-',linewidth=3.0)
ax.set_xlim3d(-200, 500)
ax.set_ylim3d(-200, 500)
ax.set_zlim3d(-200, 500)
# Set the axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Show the plot
plt.show()
ii=0
jj=0
def update_point(aa):
    global ii ,an,jj
    ii=ii+1
    if ii>=180:
        ii=0
        jj=0
    if ii>=90:
        jj=180-ii
    else:
        jj=jj+1

    X,Y,Z=DoARM(alpha,a,d,offset,thetas=[jj/180*np.pi,0,jj/180*np.pi,jj/180*np.pi,jj/180*np.pi,0])
    X=np.array(X)
    Y=np.array(Y)
    Z=np.array(Z)
    links.set_data(X, Y)
    links.set_3d_properties(Z, 'z')
    an.remove()
    an=ax.scatter(X, Y,Z, c=['r', 'g', 'b', 'c', 'm', 'y','g'], s=50)



ani=animation.FuncAnimation(fig, update_point,interval=10, fargs=())



