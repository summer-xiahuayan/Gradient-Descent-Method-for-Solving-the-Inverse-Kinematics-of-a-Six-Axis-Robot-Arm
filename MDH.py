#MDH方法
import sympy as sy
import math
import matplotlib.pyplot as plt
import numpy as np


#MDH参数
alpha=[0,-np.pi/2,0,-np.pi,-np.pi/2,np.pi/2]
a=[0,0,368,316,0,0]
offset=[-3*np.pi/2,-np.pi/3,-np.pi/2,np.pi/2,np.pi/2,0]
d=[122,140.5,0,0,-102.5,94]
point_list=[]
target=[142.400000000000,-1.88901768768479e-15,491.810000000000]



def plot_3D(point_list):
    # 创建一个新的figure
    fig = plt.figure()
    M = ['red', 'blue', 'yellow', 'orange', 'green', 'moccasin']
    # 添加一个3D坐标轴
    ax = fig.add_subplot(111, projection='3d')
    x=[point_i[0,0] for point_i in point_list]
    y=[point_i[1,0] for point_i in point_list]
    z=[point_i[2,0] for point_i in point_list]
    x.insert(0,0)
    y.insert(0,0)
    z.insert(0,0)

    an=ax.scatter(x,y,z, c=['r', 'g', 'b', 'c', 'm', 'y','g'], s=50)
    links,=ax.plot(x, y,z,'-',linewidth=3.0)
    ax.set_xlim3d(-200, 500)
    ax.set_ylim3d(-200, 500)
    ax.set_zlim3d(-200, 500)
    # Set the axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # 显示图形
    plt.show()



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

def Aiv(alpha,a,d,offset,Tnum=6):
    """
    依次右乘，推导出机械臂位姿矩阵
    """
    base_point=sy.Matrix([
        [0],
        [0],
        [0],
        [1]
    ])
    for i in range(Tnum):
        theta=sy.Symbol("q"+str(i+1))
        Ti=MTi(alpha[i],a[i],d[i],offset[i]+theta)
        ti=MTi(alpha[i],a[i],d[i],offset[i])
        if i==0:
            Ai=Ti
            ai=ti
        else:
            Ai=Ai@Ti
            ai=ai@ti

        point_list.append(ai@base_point)

    return Ai,ai



def loss_function(t1, t2, t3, t4, t5, t6, target):
    # 假设target是一个包含六个目标值的张量
    theate=[t1+offset[0], t2+offset[1], t3+offset[2], t4+offset[3], t5+offset[4], t6+offset[5]]
    A06= Aiv(alpha,a,d,theate)

    loss = math.sqrt((A06[:3, 3][0]-target[0])**2+
                     (A06[:3, 3][1]-target[1])**2+
                     (A06[:3, 3][2]-target[2])**2)
    return loss


if __name__=="__main__":
    #机械臂位姿矩阵如下
    A06,a06= Aiv(alpha,a,d,offset)
    print("位姿矩阵：",A06)
    print(f"目标点位置： x:{a06[:3, 3][0]},y:{a06[:3, 3][1]},z:{a06[:3, 3][2]}")
    plot_3D(point_list)







