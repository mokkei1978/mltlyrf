#!/usr/bin/env python
"""水位・流速の水平分布を描く

Usage: eta_vel.py IREC KL

Arguments:
  IREC  record number (4 digits)
  KL    layer number (0: surface)
"""


'''2Dアニメーション'''
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt

im=60
jm=50
km=2
dx_m=91.2e3
dy_m=111.e3

args = docopt(__doc__)
irec = args.get('IREC')
ckl = args.get('KL')

#- ファイル読み込み
f=open('../OUT/rect_mricom04/rect_mricom.ss'+irec,'rb')
nstep = np.fromfile(f,'i4',1)
t_sec = np.fromfile(f,'f4',1)
dat = np.fromfile(f,'f8',(im+2)*(jm+2)*km*3).reshape((im+2,jm+2,km,3),order='F')
f.close()

#- 海面高度
kl = int(ckl)
eta_cm = dat[:,:,kl,2] * 1.e2
x_km = (np.arange(im+2)-0.5)*dx_m * 1.e-3
y_km = (np.arange(jm+2)-0.5)*dy_m * 1.e-3

#- 流速 (ベクトルの長さは同じにして、色で流速を示す)
dec = 2  # decimate (間引く間隔)
u = np.zeros( (im+2,jm+2) )
v = np.zeros( (im+2,jm+2) )
u[1:,:] = 0.5*( dat[1:,:,kl,0] + dat[:-1,:,kl,0] )
v[:,1:] = 0.5*( dat[:,1:,kl,1] + dat[:,:-1,kl,1] )
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
plt.colorbar(Q, label='Velocity [m/s]', shrink=0.6, ax=ax)
ax.clabel(cs)
t_day = int( t_sec[0] / 3600. / 24.)
ax.set_title( f'eta [cm] k={kl} t = {t_day} [day]',fontsize=10)
ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
#    plt.pause(0.1)

#plt.show()
plt.savefig('temp.png', bbox_inches='tight')
plt.savefig('eta'+ckl+'_'+irec+'.png', bbox_inches='tight')

