'''2Dアニメーション'''
import numpy as np
import matplotlib.pyplot as plt

im=500
jm=500
km=2
dx_m=10.e3
dy_m=10.e3

#- ファイル読み込み
f=open('../OUT/test2/test2.ss0500','rb')
nstep = np.fromfile(f,'i4',1)
t_sec = np.fromfile(f,'f4',1)
dat = np.fromfile(f,'f8',(im+2)*(jm+2)*km*3).reshape((im+2,jm+2,km,3),order='F')
f.close()

#- 海面高度
k = 0
eta_cm = dat[:,:,k,2] * 1.e2
x_km = (np.arange(im+2)-0.5)*dx_m * 1.e-3
y_km = (np.arange(jm+2)-0.5)*dy_m * 1.e-3

#- 流速 (ベクトルの長さは同じにして、色で流速を示す)
dec = 20  # decimate (間引く間隔)
u = np.zeros( (im+2,jm+2) )
v = np.zeros( (im+2,jm+2) )
u[1:,:] = 0.5*( dat[1:,:,k,0] + dat[:-1,:,k,0] )
v[:,1:] = 0.5*( dat[:,1:,k,1] + dat[:,:-1,k,1] )
u = u[::dec,::dec]
v = v[::dec,::dec]
u_abs=np.sqrt( pow(u,2) + pow(v,2) )
u = u / u_abs
v = v / u_abs
xu_km = x_km[::dec]
yu_km = y_km[::dec]

#- キャンバス
fig, ax = plt.subplots()

#for i in range(1):  # 3つ目の引数が描画の間隔
#    ax.cla()
cs = ax.contour(x_km,y_km,eta_cm.T)
    #- index order of eta is reversed, so eta.T (transpose) are used.
Q = ax.quiver(xu_km, yu_km, u[:,:].T, v[:,:].T, u_abs[:,:].T, cmap='jet', pivot='mid' )
#    plt.colorbar(Q, label='Velocity [m/s]', shrink=0.6, ax=ax) うまくいかない
ax.clabel(cs)
ax.set_title( f'SSH [cm] t = {t_sec}',fontsize=10)
ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
#    plt.pause(0.1)

plt.show()


