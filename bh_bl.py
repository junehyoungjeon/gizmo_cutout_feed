from process import *
import random
from scipy.optimize import curve_fit
import matplotlib
import mpl_toolkits.mplot3d
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from miniquasar import tGW, coaltime, sigma_inf, rho_inf, Macc_edd, binary_gen, ms_mbh_default
import sys
from metalicity import *
from gwform import *

fac0 = 0.3783187*0.3497985*0.399293243
fac = 0.3783187*0.3497985*0.399293243
#fac = 1
bs = 4
V = (bs/0.6774)**3 * fac
Zsun = 0.02

Mup = lambda z: 2.5e7*((1+z)/10)**-1.5
Mdown = lambda z: 1e6*((1+z)/10)**-2

def sumy(l):
	if len(l)>1:
		return np.array([np.min(l), np.max(l), np.average(l), np.median(l), np.std(l)])
	else:
		return np.array([np.min(l), np.max(l), np.average(l), np.median(l)])

def psfr(l, lb, ls = '-', mode = 0, V = V):
	if mode==0:
		plt.plot(1./l[0]-1., l[2]/V, ls=ls, label=lb)
	else:
		plt.plot(1./l[0]-1., l[3], ls=ls, label=lb)

def psfr_tot(l1, l2, lb, ls = '-', mode = 0, V = V):
	if mode==0:
		plt.plot(1./l1[0]-1, (l1[2]+l2[2])/V, ls=ls, label=lb)
	else:
		plt.plot(1./l1[0]-1, (l1[3]+l2[3]), ls=ls, label=lb)

def sfrmd(t, t0, a=0.182, b=1.26, c=1.865, d=0.071):
	return a*(t**b*np.exp(-t/c)+d*np.exp(d*(t-t0)/c))

def sfrbf(z):
	return 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)

def sfrtf(z, a, b, c, d):
	#t = 1/(z+1)
	#return a*(t**b*np.exp(-t/c)+d*np.exp(d*(t-1)/c)) 
	return a*(1+z)**b/(1+((1+z)/c)**d)

def LWbg(z, sfrd3, sfrd2, eta = 2e4):
	return (eta/1e4)*(sfrd3 + sfrd2/5.) / 1e-2 * ((1+z)/10)**3

def countBH(ds, Zth = 100.):
	z = ds['Redshift']
	keys = ds.field_list
	tag = np.sum([x[0] == 'PartType3' for x in keys])
	if tag>0:
		ad = ds.all_data()
		poss = np.array(ad[('PartType3','Coordinates')].to('kpccm/h'))
		sid = ad[('PartType3','ParticleIDs')]
		#star = np.array(ad['PartType4', 'Star Type'])
		bht = np.array(ad['PartType3', 'BH Type'])
		Z = np.array(ad['PartType3', 'Metallicity_00'])/Zsun
		jlw = np.array(ad[('PartType3','J21_SF')])
		h2sf = np.array(ad[('PartType3','H2_SF')])
		#print('JLW_min = {}, JLW_max = {}, JLW_bar = {}, JLW_med = {}'.format(*sumy(jlw)))
		select = (bht==2)#*(Z<Zth)
		if np.sum(select)>0:
			print('JLW_SF:', sumy(jlw[select]))
			print('H2_SF:', sumy(h2sf[select]))
			print('Z_SF:', sumy(Z[select]))
			print(poss[select])
			#print('Star Type:', star[select])
		n1 = np.sum(bht==1) #len(bht[bht==1])
		n2 = np.sum(bht==2) #len(bht[bht==2])
		#n3 = len(bht[select])
		return [z, n1, n2, poss[select], sid[select]]#, n3]
	else:
		return [z, 0, 0]#, 0]

def smass(ds, Zcr = 0.1):
	ad = ds.all_data()
	keys = ds.field_list	
	bflag = np.sum([x[0] == 'PartType3' for x in keys])
	lmass = np.array(ad[('PartType4','Masses')].to('Msun'))
	stype = ad[('PartType4', 'Star Type')]
	lm2 = lmass[stype == 2]
	if bflag>0:
		#lmbh = np.array(ad[('PartType3','Masses')].to('Msun'))
		bht = ad[('PartType3','BH Type')]
		sel = bht==1
		Mbh = bht[sel].shape[0] * lmass[0] #np.sum(lmbh[bht == 1])
	M2 = np.sum(lm2)
	M1 = np.sum(lmass) - M2 + Mbh
	print('Mass in PopIII BH particles: {:.0f} Msun Mpc^-3'.format(Mbh/V))
	return [M1, M2]

new_3Dax = lambda x: plt.subplot(x,projection='3d')

def plot3d(fn, ID, R, ds, nump=1e5, alp=0.5, age = 1e7):
	llb = [r'$x\ [h^{-1}\mathrm{kpc}]$', r'$y\ [h^{-1}\mathrm{kpc}]$', r'$z\ [h^{-1}\mathrm{kpc}]$']
	ad0 = ds.all_data()
	keys = ds.field_list
	
	bflag = np.sum([x[0] == 'PartType3' for x in keys])
	if bflag > 0:
		bid = ad0[('PartType3','ParticleIDs')]
		posb = np.array(ad0[('PartType3','Coordinates')].to('kpccm/h'))
		sel0 = bid==ID
	else:
		sel0 = [0]
	
	sflag = np.sum([x[0] == 'PartType4' for x in keys])
	if sflag > 0:
		sid = ad0[('PartType4','ParticleIDs')]
		poss = np.array(ad0[('PartType4','Coordinates')].to('kpccm/h'))
		sel1 = sid==ID
	else:
		sel1 = [0]
		
	gid = ad0[('PartType0','ParticleIDs')]
	posg = np.array(ad0[('PartType0','Coordinates')].to('kpccm/h'))
	sel2 = gid==ID
	if np.sum(sel0)>0:
		pos = posb[sel0][0]
	if np.sum(sel1)>0:
		pos = poss[sel1][0]
	if np.sum(sel2)>0:
		pos = posg[sel2][0]
	ad = ds.sphere(pos, (R, 'kpccm/h'))
	if sflag > 0:
		poss = np.array(ad[('PartType4','Coordinates')].to('kpccm/h'))
		lts = np.array(ad[('PartType4', 'Stellar Age')])
	else:
		poss = []
	if bflag > 0:
		posb = np.array(ad[('PartType3','Coordinates')].to('kpccm/h'))
	else:
		posb = []
	posg = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))
	posd = np.array(ad[('PartType1','Coordinates')].to('kpccm/h'))
	T = temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')])
	print('Ngas_ngb = {} for ID = {} with hBH = {} kpc'.format(len(T), ID, R))
	lpd = np.random.choice(posd.shape[0], min(int(nump),posd.shape[0]),replace=False)
	lpg = np.random.choice(posg.shape[0], min(int(nump),posg.shape[0]),replace=False)
	posd = posd[lpd].T
	posg = posg[lpg]
	lT = T[lpg]
	hot = lT>1e4
	cold = lT<=1e4
	posgh = posg[hot].T
	posgc = posg[cold].T
	ms = (2*nump/(posg.shape[1] + posd.shape[1]))**0.5 * 0.05
	plt.figure(figsize=(5,5))
	ax = new_3Dax(111)
	#ax.set_aspect('equal','box')
	ax.plot(*posd, '.', markersize=ms, label='DM', alpha = alp, zorder=0)
	ax.plot(*posgc, '.', markersize=ms, label='Cold gas', alpha = alp, zorder=2)
	if len(poss)>0:
		lps = np.random.choice(poss.shape[0], min(int(nump),poss.shape[0]),replace=False)
		poss = poss[lps]
		lts = lts[lps]
		old = lts > age
		young = lts <= age
		po = poss[old].T
		py = poss[young].T
		mss = min((2*nump/poss.shape[0])**0.5 * 0.1, 7)
		ax.plot(*po, '*', label='Old star', markersize=1.5*mss, alpha = alp, zorder=3)
		ax.plot(*py, 'h', label='Young star', markersize=1.5*mss, alpha = alp, zorder=4)
	if len(posb)>0:
		posb = posb.T
		ax.plot(*posb, '^', label='BH', markersize=ms, alpha = alp, zorder=4)
	ax.plot(*posgh, '.', markersize=ms, label='Hot gas', color='r', alpha = alp, zorder=1)
	#ax.plot(*[[x] for x in pos], 'o', color='g', label='Target particle', alpha=0.5, zorder=5)
	ax.set_xlabel(llb[0])
	ax.set_ylabel(llb[1])
	ax.set_zlabel(llb[2])
	ax.set_xlim(pos[0]-R, pos[0]+R)
	ax.set_ylim(pos[1]-R, pos[1]+R)
	ax.set_zlim(pos[2]-R, pos[2]+R)
	plt.legend()
	plt.tight_layout()
	plt.savefig(fn)
	#plt.show()
	plt.close()

def LW_bg(sn, ds, le = [1500, 1500, 1950], re = [2500, 2500, 2050], ax = [0, 1], rep = './', mode=0, nump=1e6):
	llb = [r'$x\ [h^{-1}\mathrm{kpc}]$', r'$y\ [h^{-1}\mathrm{kpc}]$', r'$z\ [h^{-1}\mathrm{kpc}]$']
	ad = ds.box(le, re)
	posg = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))
	sel = np.random.choice(posg.shape[0], min(int(nump),posg.shape[0]),replace=False)
	posg = posg[sel].T
	jlw = np.array(ad[('PartType0','J21')][sel])
	bht = np.array(ad['PartType3', 'BH Type'])
	poss = np.array(ad[('PartType3','Coordinates')].to('kpccm/h'))
	posbh1 = poss[bht==1].T
	posbh2 = poss[bht==2].T
	plt.figure()
	axis = plt.subplot(111)
	axis.set_aspect(aspect=1)
	plt.tripcolor(posg[ax[0]], posg[ax[1]], np.log10(jlw),cmap=plt.cm.binary)
	cb = plt.colorbar()
	plt.clim(-.5,2.5)
	cb.set_label(r'$\log(J_{21})$')
	if mode!=0:
		plt.plot(*posbh1[ax], '.', color='g', label='Stellar BH', alpha=0.5)
		plt.plot(*posbh2[ax], '.', color='purple', label='DCBH', alpha=0.5)
		plt.legend()
	plt.xlim(le[ax[0]], re[ax[0]])
	plt.ylim(le[ax[1]], re[ax[1]])
	plt.xlabel(llb[ax[0]])
	plt.ylabel(llb[ax[1]])
	plt.tight_layout()
	plt.savefig(rep+'LW_bg_'+str(sn)+'_'+str(np.sum(ax))+'.png',dpi=300)
	plt.close()
	print(sumy(jlw))

def temp_bg(sn, ds, le = [1500, 1500, 1950], re = [2500, 2500, 2050], ax = [0, 1], rep = './', mode=0, nump=1e6, lflag=0, wh=6.5, lss=[1, 8, 64]):
	llb = [r'$x\ [h^{-1}\mathrm{kpc}]$', r'$y\ [h^{-1}\mathrm{kpc}]$', r'$z\ [h^{-1}\mathrm{kpc}]$']
	ad = ds.box(le, re)
	posg = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))
	sel = np.random.choice(posg.shape[0], min(int(nump),posg.shape[0]),replace=False)
	posg = posg[sel].T
	lT = temp(ad[('PartType0','InternalEnergy')][sel],ad[('PartType0','Primordial HII')][sel])
	bht = np.array(ad['PartType3', 'BH Type'])
	poss = np.array(ad[('PartType3','Coordinates')].to('kpccm/h'))
	posst = np.array(ad[('PartType4','Coordinates')].to('kpccm/h')).T
	posbh1 = poss[bht==1].T
	posbh2 = poss[bht==2].T
	if lflag==0:
		asp = (re[0]-le[0])/(re[1]-le[1])
		fsz = (wh, wh/asp)
		plt.figure(figsize=fsz)
	else:
		plt.figure()
	axis = plt.subplot(111)
	axis.set_aspect(aspect=1)
	plt.tripcolor(posg[ax[0]], posg[ax[1]], np.log10(lT),cmap=plt.cm.seismic)
	plt.clim(2.5,5.5)
	if lflag>0:
		cb = plt.colorbar()
		cb.set_label(r'$\log(T\ [\mathrm{K}])$')
	if mode!=0:
		if lflag==0:
			nc = 1
			plt.scatter(*posst[ax], marker='.', c='orange', s=lss[0], alpha=1.0, label='Star')
		else:
			nc = 0
		plt.scatter(*posbh1[ax], marker='.', color='g', label='Stellar BH', s=lss[1], alpha=0.5)
		plt.scatter(*posbh2[ax], marker='.', color='purple', label='DCBH', alpha=0.5, zorder=2, s=lss[2])
		if lflag>0:
			lgnd = plt.legend()#markerscale=3)
			#lgnd.legendHandles[0]._sizes = [32]
			lgnd.legendHandles[nc]._sizes = [64]
			lgnd.legendHandles[nc+1]._sizes = [256]
		else:
			plt.tick_params(top=False, bottom=False, left=False, right=False, 
			labelleft=False, labelbottom=False)
		#for legend_handle in lgd.legendHandles:
		#	legend_handle._legmarker.set_markersize(9)
	plt.xlim(le[ax[0]], re[ax[0]])
	plt.ylim(le[ax[1]], re[ax[1]])
	if lflag>0:
		plt.xlabel(llb[ax[0]])
		plt.ylabel(llb[ax[1]])
	plt.tight_layout()
	if lflag==0:
		plt.savefig(rep+'temp_'+str(sn)+'.png',dpi=666)
	else:
		plt.savefig(rep+'temp_'+str(sn)+'_'+str(np.sum(ax))+'.png',dpi=300)
	plt.close()
	print(sumy(lT))

def denspro(sn, ds, le = [1500, 1500, 1950], re = [2500, 2500, 2050], ax = [0, 1], nb=100, rep = './', norm = 1e6, mode=0):
	h = ds['HubbleParam']
	llb = [r'$x\ [h^{-1}\mathrm{kpc}]$', r'$y\ [h^{-1}\mathrm{kpc}]$', r'$z\ [h^{-1}\mathrm{kpc}]$']
	ad = ds.box(le, re)
	posg = np.array(ad[('PartType0','Coordinates')].to('kpccm/h')).T
	bht = np.array(ad['PartType3', 'BH Type'])
	poss = np.array(ad[('PartType3','Coordinates')].to('kpccm/h'))
	posbh1 = poss[bht==1].T
	posbh2 = poss[bht==2].T
	obj = recaesar(sn, rep)
	posh = np.array([x.pos.to('kpccm/h') for x in obj.halos]).T
	lselh = [(posh[i]>le[i]) * (posh[i]<re[i]) for i in range(3)]
	selh = lselh[0] * lselh[1] * lselh[2]
	posh = posh.T
	posh = posh[selh]
	posh = posh.T
	lm = np.array([(1e10*x.virial_quantities['circular_velocity']**2*x.radii['virial'].to('kpc/h')/G) for x in obj.halos])[selh]
	
	plt.figure()
	axis = plt.subplot(111)
	axis.set_aspect(aspect=1)
	plt.hist2d(posg[ax[0]],posg[ax[1]],bins=nb,norm=LogNorm(),cmap=plt.cm.coolwarm)
	cb = plt.colorbar()
	plt.clim(1, 5e3)
	cb.set_label(r'$\log(N)$')
	plt.scatter(posh[ax[0]],posh[ax[1]],s=lm/norm,edgecolors='g',facecolors='none',lw=1,label='Halo')
	#if mode!=0:
		#plt.plot(*posbh1[ax], '.', color='g', label='Stellar BH', alpha=0.5)
		#plt.plot(*posbh2[ax], '.', color='purple', label='DCBH', alpha=0.5)
	plt.xlim(le[ax[0]], re[ax[0]])
	plt.ylim(le[ax[1]], re[ax[1]])
	plt.xlabel(llb[ax[0]])
	plt.ylabel(llb[ax[1]])
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'denspro_'+str(sn)+'_'+str(np.sum(ax))+'.png')
	plt.close()

def pbd(x, Bxs):
	# make sure that the displacement between two particles
	# is calculated with their nearest images 
	return x - Bxs*(x>Bxs/2.) + Bxs*(x<-Bxs/2.)

def pbd2(x, Bxs):
	# convert all particles outside the box to their images inside the box
	return x + Bxs*(x<0) - Bxs*(x>Bxs)

def profile0(pos, cen, R, nbin = 20, fac=1e-2, mode = 1, Bxs = 4000):
	lpart = np.array(pos)
	rcen = np.array(cen)
	lr0 = np.array([np.linalg.norm(pbd(x-rcen, Bxs)) for x in lpart])
	within = lr0<R
	lr = lr0[within]
	lr = sorted(lr, key=lambda x: x)
	if mode==0:
		lrb1 = np.linspace(fac*R, R, nbin)
	else:
		lrb1 = np.geomspace(lr[0]*2, R, nbin)
	lrb0 = np.hstack([[0.0], lrb1])
	lrb = 0.5*(lrb0[1:] + lrb0[:-1])
	lvs = (np.pi*4/3)*(lrb0[1:]**3-lrb0[:-1]**3)
	lp = np.zeros(nbin)
	count = 0
	for r in lr:
		while r>=lrb0[count+1]:
			count += 1
		if count>=nbin:
			break
		lp[count] += 1
	sel = lp>0
	d = {}
	d['pro'] = [lrb[sel], (lp/lvs)[sel]]+ \
		[(lp**0.5/lvs)[sel], np.sum(lp), lvs[sel]]
		# radius bins, density, 
		# Poisson error, total particle number, volume of shells
	d['part'] = within # particles used in the profle calculation 
		# denoted by a boolean array
	#d['cen'] = rcen # CoM
	return d
	
def profile1(pos, cen, yin, R, nbin = 20, fac=1e-2, mode=1, Bxs = 4000):
	lpart = np.array(pos)
	rcen = np.array(cen)
	lr0 = np.array([np.linalg.norm(pbd(x-rcen, Bxs)) for x in lpart])
	within = lr0<R
	lr = lr0[within]
	ly0 = yin[within]
	data = sorted(np.array([lr, ly0]).T, key=lambda x: x[0])
	data = np.array(data).T
	lr, ly = data
	if mode==0:
		lrb1 = np.linspace(fac*R, R, nbin)
	else:
		lrb1 = np.geomspace(lr[0]*2, R, nbin)
	lrb0 = np.hstack([[0.0], lrb1])
	lrb = 0.5*(lrb0[1:] + lrb0[:-1])
	lvs = (np.pi*4/3)*(lrb0[1:]**3-lrb0[:-1]**3)
	lp = np.zeros(nbin)
	lcol = [[] for i in range(nbin)]
	count = 0
	for r, y in zip(lr, ly):
		while r>=lrb0[count+1]:
			count += 1
		if count>=nbin:
			break
		lp[count] += 1
		lcol[count].append(y)
	sel = lp>0
	lden = np.zeros(nbin)
	lstd = np.zeros(nbin)
	for i in range(nbin):
		col = lcol[i]
		if lp[i]==0:
			lden[i], lstd[i] = 0, 0
			continue
		lden[i] = np.median(col)
		#lden[i] = np.average(col)
		if lp[i]==1:
			lstd[i] = lden[i]/2**0.5
		else:
			lstd[i] = np.std(col)
	d = {}
	d['pro'] = [lrb[sel], lden[sel], lstd[sel], np.sum(lp), lvs[sel]]
	d['data'] = data
	return d
	
def profile2(pos, cen, yin, R, nbin = 20, fac=1e-2, mode=1, Bxs = 4000):
	lpart = np.array(pos)
	rcen = np.array(cen)
	dis0 = np.array([pbd(x-cen, Bxs) for x in lpart])
	lr0 = np.array([np.linalg.norm(x) for x in dis0])
	within = lr0<R
	lr = lr0[within]
	dis = dis0[within]
	ly_raw = yin[within]
	ymean = np.average(ly_raw, axis=0)
	ly0 = -np.array([np.dot(x, y-ymean) for x, y in zip(dis, ly_raw)])/lr
	data = sorted(np.array([lr, ly0]).T, key=lambda x: x[0])
	data = np.array(data).T
	lr, ly = data
	if mode==0:
		lrb1 = np.linspace(fac*R, R, nbin)
	else:
		lrb1 = np.geomspace(lr[0]*2, R, nbin)
	lrb0 = np.hstack([[0.0], lrb1])
	lrb = 0.5*(lrb0[1:] + lrb0[:-1])
	lvs = (np.pi*4/3)*(lrb0[1:]**3-lrb0[:-1]**3)
	lp = np.zeros(nbin)
	lcol = [[] for i in range(nbin)]
	count = 0
	for r, y in zip(lr, ly):
		while r>=lrb0[count+1]:
			count += 1
		if count>=nbin:
			break
		lp[count] += 1
		lcol[count].append(y)
	lden = np.zeros(nbin)
	lstd = np.zeros(nbin)
	for i in range(nbin):
		col = lcol[i]
		if lp[i]==0:
			lden[i], lstd[i] = 0, 0
			continue
		lden[i] = np.average(col)
		if lp[i]==1:
			lstd[i] = lden[i]/2**0.5
		else:
			lstd[i] = np.std(col)
	sel = lp>0
	d = {}
	d['pro'] = [lrb[sel], lden[sel], lstd[sel], np.sum(lp), lvs[sel]]
	d['data'] = data
	return d

oi = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#000000']
mst = ['o', '^', 's', 'D', 'h', 'p']

aplog = lambda x: np.sign(x)*(np.log10(np.abs(x)+1))
rga = lambda l: [np.min(l), np.max(l)]
nth_bh, nth_sf, nth_hii = 2e3, 1e2, 1e3
T1, T2 = 7e3, 1e4
Zth = 2e-4
h2th = 1e-6

def main(sn, nbh = 3, R = 2, dm = 1, rep = './', x0 = 1e-2, nbin = 30, fac=1e-2, mode=1, Bxs = 4000, alp = 0.1, ms = 0.5, cs = oi, mst = mst):
	ds = yt.load(rep+'snapshot_'+str(sn).zfill(3)+'.hdf5')
	z = ds['Redshift']
	h = ds['HubbleParam']
	a = 1./(1.+z)
	ad0 = ds.all_data()
	keys = ds.field_list
	tag = np.sum([x[0] == 'PartType3' for x in keys])
	if tag<=0:
		print('There are no BHs!')
		return -1
	poss = np.array(ad0[('PartType3','Coordinates')].to('kpccm/h'))
	bht = np.array(ad0['PartType3', 'BH Type'])
	sel = bht==2
	n0 = np.sum(sel)
	if n0==0:
		print('There is no DCBH!')
		return -2
	n = int(min(n0, nbh))
	idbh = np.array(ad0['PartType3', 'ParticleIDs'])[sel]#
	posbh = poss[sel]#[0:n]
	Zbh = np.array(ad0['PartType3', 'Metallicity_00'])[sel]/Zsun #[0:n]/Zsun
	h2sf = np.array(ad0[('PartType3','H2_SF')])[sel]#[0:n]
	data_bh = [[i, p, z, h] for i, p, z, h in zip(idbh, posbh, Zbh, h2sf)]
	data_bh = np.array(sorted(data_bh, key = lambda x: x[0]))
	idbh = data_bh[:,0][0:n]
	posbh = np.array(data_bh[:,1])[0:n]
	Zbh = np.array(data_bh[:,2])[0:n]
	h2sf = np.array(data_bh[:,3])[0:n]
	if dm==0:
		norm = np.array(ad0[('PartType0', 'Masses')][0].to('g')/(KPC*a/h)**3)*Hfrac/PROTON
	else:
		norm = Hfrac/PROTON
	#"""	
	plt.figure(figsize=(8, 6))
	ax0 = plt.subplot(221)
	ax1 = plt.subplot(222)
	ax2 = plt.subplot(223)
	ax3 = plt.subplot(224)
	for i in range(n):
		cen = posbh[i]
		ad = ds.sphere(cen, (R, 'kpccm/h'))
		pos = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))
		T = temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')])
		Z = np.array(ad['PartType0', 'Metallicity_00'])/Zsun
		xh2 = np.array(ad['PartType0', 'Primordial H2'])
		if dm==0:
			prho = profile0(pos, cen, R, nbin, fac, mode, Bxs)
		else:
			rho = np.array(ad[('PartType0','Density')].to('g/cm**3'))
			prho = profile1(pos, cen, rho, R, nbin, fac, mode, Bxs)
		ptemp = profile1(pos, cen, T, R, nbin, fac, mode, Bxs)
		pZ = profile1(pos, cen, Z, R, nbin, fac, mode, Bxs)
		ph2 = profile1(pos, cen, xh2, R, nbin, fac, mode, Bxs)
		p0, p1 = prho['pro'], ptemp['pro']
		p2, p3 = pZ['pro'], ph2['pro']
		if dm!=0:
			d0, d1 = prho['data'], ptemp['data']
		d2, d3 = pZ['data'], ph2['data']
		ax0.errorbar(p0[0], p0[1]*norm, p0[2]*norm, ls=lls[i], zorder = 2, color = cs[i])
		ax1.errorbar(p1[0], p1[1], p1[2], ls=lls[i], label=r'DCBH_'+str(i), zorder = 2, color = cs[i])
		if dm!=0:
			ax0.plot(d0[0], d0[1]*norm, '.', alpha = alp, zorder = 1, color = cs[i])
			ax1.plot(d1[0], d1[1], '.', alpha = alp, zorder = 1, color = cs[i])
		
		ax2.plot([x0*2], [Zbh[i]], mst[i], color = cs[i], alpha = 0.5)
		ax3.plot([x0*2], [h2sf[i]], mst[i], color = cs[i], alpha = 0.5)
		ax2.plot(p2[0], p2[1], ls=lls[i], zorder = 2, color = cs[i])
		ax3.plot(p3[0], p3[1], ls=lls[i], zorder = 2, color = cs[i])
		ax2.plot(d2[0], d2[1], '.', alpha = alp, zorder = 1, color = cs[i])
		ax3.plot(d3[0], d3[1], '.', alpha = alp, zorder = 1, color = cs[i])

	ax0.plot([x0, R], [nth_bh, nth_bh], 'k--', label=r'$n_{\mathrm{th}}$(BH)')
	ax0.plot([x0, R], [nth_sf, nth_sf], 'k-.', label=r'$n_{\mathrm{th}}$(SF)')
	ax0.plot([x0, R], [nth_hii, nth_hii], 'k:', label=r'$n_{\mathrm{th}}$(HII)')
	ax1.fill_between([x0, R], [T1, T1], [T2, T2], facecolor='k', zorder = 0, label=r'$T$(BH)', alpha = alp*3)
	ax2.plot([x0, R], [Zth, Zth], 'k--', label=r'$Z_{\mathrm{th}}$')
	ax3.plot([x0, R], [h2th, h2th], 'k--', label=r'$x_{\mathrm{H_{2},th}}$')

	ax0.set_ylabel(r'$n_{\mathrm{H}}\ [\mathrm{cm^{-3}}]$')
	ax1.set_ylabel(r'$T\ [\mathrm{K}]$')
	ax2.set_ylabel(r'$Z\ [\mathrm{Z_{\odot}}]$')
	ax3.set_ylabel(r'$\mathrm{[H_{2}/H]}$')
	ax0.set_yscale('log')
	ax1.set_yscale('log')
	ax2.set_yscale('log')
	ax3.set_yscale('log')
	#ax0.set_ylim(5e-1, 5e3)
	#ax1.set_ylim(1e2, 3e4)
	#ax2.set_ylim(1e-5, 1)
	#ax3.set_ylim(9e-11, 1e-4)
	
	ax0.set_xlim(x0, R)
	ax1.set_xlim(x0, R)
	ax2.set_xlim(x0, R)
	ax3.set_xlim(x0, R)
	ax0.set_xlabel(r'$r\ [h^{-1}\mathrm{kpc}]$')
	ax1.set_xlabel(r'$r\ [h^{-1}\mathrm{kpc}]$')
	ax2.set_xlabel(r'$r\ [h^{-1}\mathrm{kpc}]$')
	ax3.set_xlabel(r'$r\ [h^{-1}\mathrm{kpc}]$')
	ax0.set_xscale('log')
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')
	ax0.legend(loc=3)
	ax1.legend(loc=4)
	ax2.legend(loc=2)
	ax3.legend(loc=4)	

	plt.tight_layout()
	plt.savefig(rep+'DCBHpro_'+str(sn)+'.pdf')
	#plt.show()
	plt.close()	
	#"""	
	#umass = np.array(ad0[('PartType0', 'Masses')][0].to('Msun'))
	uvel = 1e5
	ul = a*KPC/h
	plt.figure()
	for i in range(n):
		cen = posbh[i]
		ad = ds.sphere(cen, (R, 'kpccm/h'))
		pos = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))
		vel = np.array(ad[('PartType0','Velocities')].to('km/s'))
		pro = profile2(pos, cen, vel, R, nbin, fac, mode, Bxs)
		if dm==0:
			prho = profile0(pos, cen, R, nbin, fac, mode, Bxs)
		else:
			rho = np.array(ad[('PartType0','Density')].to('g/cm**3'))
			prho = profile1(pos, cen, rho, R, nbin, fac, mode, Bxs)
		p = pro['pro']
		pr = prho['pro']
		lr = p[0]
		lm = np.cumsum(p[4]*pr[1]) * ul**3 *1e10/UM
		lacc = 4*np.pi*lr**2*p[1]*pr[1] * uvel * ul**2 *1e10/UM * YR
		lacc_er = 4*np.pi*lr**2*p[1]*pr[2] * uvel * ul**2 *1e10/UM * YR
		plt.errorbar(lm, lacc, lacc_er, ls=lls[i], zorder = 2, color = cs[i], label=r'DCBH_'+str(i))
		#plt.plot(lm, lacc, ls=lls[i], zorder = 2, color = cs[i], label=r'DCBH_'+str(i))
	xr = [1e4, 1e8]
	plt.plot(xr, [0, 0], 'k--', lw=0.5)
	plt.plot(xr, [0.04, 0.04], 'k-.', label=r'$\mathrm{\dot{M}_{crit}}$')
	plt.xlim(xr)
	plt.ylim(-0.5, 1.)
	plt.xscale('log')
	#plt.yscale('symlog')
	plt.xlabel(r'Enclosed mass [$\mathrm{M_{\odot}}$]')
	plt.ylabel(r'Accretion rate [$\mathrm{M_{\odot}\ yr^{-1}}$]')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig(rep+'Mdot_M_'+str(sn)+'.pdf')
	plt.close()
	print('Redshift: {:.3f}'.format(z))
	print(Zbh)
	print(h2sf)
	print(posbh)
	print(idbh)

def findhost(pos, obj, fac = 2., Bxs = 4000, h = 0.6774, mode=0):
	lh = obj.halos
	out = lh[0]
	rmin = fac
	lr = []
	for halo in lh:
		r200 = np.array(halo.radii['r200'].to('kpccm/h'))
		hpos = np.array(halo.pos.to('kpccm/h'))
		dis = pbd(hpos - pos, Bxs)
		r = np.linalg.norm(dis)
		if mode==0:
			rat = r/r200
		else:
			rat = r
		lr.append(r)
		if rat <= rmin:
			rmin = rat
			out = halo
	return [out, rmin]
	
def track(halo, sns, snf, fac = 5., rep = './', base = 'snapshot', ext = '.hdf5', post = 'caesar', Bxs  =4000):
	ds0 = yt.load(rep+base+'_'+str(sns).zfill(3)+ext)
	hub = ds0['HubbleParam']
	z0 = ds0['Redshift']
	ad0 = ds0.all_data()
	IDlist = ad0[('PartType1', 'ParticleIDs')][halo.dmlist]
	vel = np.linalg.norm(np.array(halo.vel))
	pos = np.array(halo.pos.to('kpccm/h'))
	r200 = np.array(halo.radii['r200'].to('kpccm/h'))
	sn = sns - 1
	lfrac = [1.0]
	lz = [z0]
	lhalo = [halo]
	while sn >= snf:
		obj = caesar.load(rep+post+'_'+str(sn).zfill(3)+ext)
		ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
		ad = ds.all_data()
		IDs = ad[('PartType1', 'ParticleIDs')]
		z = ds['Redshift']
		dt = np.abs(TZ(z) - TZ(z0))
		R = fac * max(dt * vel * 1e5 * (1+z)/KPC * hub, r200)
		frac = 0.0
		flag = 0
		for h in obj.halos:
			pos_ = np.array(h.pos.to('kpccm/h'))
			r = np.linalg.norm(pbd(pos_ - pos, Bxs))
			if r > R:
				continue
			flag = 1
			IDlist_ = IDs[h.dmlist]
			overlap = np.intersect1d(IDlist, IDlist_)
			frac_ = len(overlap)/len(IDlist)
			if frac_ > frac:
				frac = frac_
				sel = h
		if flag==0:
			break;
		sn -= 1
		print(frac)
		lfrac.append(frac)
		lhalo.append(sel)
		lz.append(z)
		z0 = z
		IDlist = IDs[sel.dmlist]
		vel = np.linalg.norm(np.array(sel.vel))
		pos = np.array(sel.pos.to('kpccm/h'))
		r200 = np.array(sel.radii['r200'].to('kpccm/h'))
	return np.array([lhalo, lz, lfrac])
	
def haloget(his, eps = 0.1):
	d = {}
	lh, lz, lfrac = his
	d['z'] = lz
	d['pos'] = np.array([h.pos.to('kpccm/h') for h in lh])
	d['mtot'] = np.array([h.masses['total'] for h in lh])
	d['mdm'] = np.array([h.masses['dm'] for h in lh])
	d['mgas'] = np.array([h.masses['gas'] for h in lh])
	d['mstr'] = np.array([h.masses['stellar'] for h in lh])*0.1
	d['temp'] = np.array([h.temperatures['mass_weighted'] for h in lh])
	d['vtemp'] = np.array([h.temperatures['virial'] for h in lh])
	d['r200'] = np.array([h.radii['r200'].to('kpccm/h') for h in lh])
	d['frac'] = lfrac
	return d

def SFgas(sn, rep = './', base = 'snapshot', ext = '.hdf5', X = 0.76):
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	ad = ds.all_data()
	gas = np.array(ad[('PartType0', 'Masses')].to('Msun'))
	mgas = np.max(gas)
	star = np.array(ad[('PartType4', 'Masses')].to('Msun'))
	mstar = star[0]
	mstar_tot = np.sum(star)
	lm = gas[gas<mgas]
	lms = mgas - lm
	ms_tot = np.sum(lms)
	Ng = int(mgas/mstar+0.1)
	lfrac = np.zeros(Ng)
	lfrac_ = np.zeros(Ng)
	for i in range(1, Ng):
		lfrac[i-1] = np.sum(lms[lms/mstar+0.1>i])
		lfrac_[i-1] = np.sum(lms[lms/mstar+0.1>i] - i*mstar)
	lfrac += mstar_tot - ms_tot
	return lfrac/mstar_tot, lfrac_/mstar_tot
	
def BHmass_spec(sn, rep = './', tp = 3, tbh = 1, base = 'snapshot', ext = '.hdf5', X = 0.76, mode = 0, nbm=100):
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	z = ds['Redshift']
	zf = '{:.2f}'.format(z)
	h = ds['HubbleParam']
	a = 1./(1.+z)
	keys = ds.field_list
	#print(keys)
	plab = 'PartType'+str(tp)
	tag = np.sum([x[0] == plab for x in keys])
	if tag==0:
		print('No BH particles!')
		return 0
	ad = ds.all_data()
	bht = ad[(plab, 'BH Type')]
	sel = bht == tbh
	mass = np.array(ad[(plab, 'BH Mass')][sel])
	massm = np.array(ad[(plab, 'BH Max mass')][sel])
	mdot = np.array(ad[plab, 'BH_Mdot'])
	mdotmax = 0
	for i in range(len(mdot)):
		if mdot[i]>mdotmax:
			bhind = i
			mdotmax = mdot[i]
	print('Mdot: ', sumy(mdot))
	print('Maximum Mdot: {} [Msun yr^-1]'.format(mdotmax))
	bhm = np.array(ad[(plab, 'BH_Mass')])
	bhmacc = np.array(ad[(plab, 'BH_Macc')])
	medd = Macc_edd(bhm, 0.125)
	print('Macc/MBH: ', sumy(bhmacc/bhm))
	rat = mdot/medd
	print('Mdot/Medd: ', sumy(rat))
	bhmacc += (bhm - bhmacc) * (bhm < bhmacc)
	#print('Mdot: ', mdot)
	bhdens = np.array(ad[(plab, 'Density')].to('g/cm**3'))[sel]
	#print('Density: ', bhdens * X / PROTON)
	nb = min(int(len(mass)/30),nbm)
	#print(mass.shape, nb)
	plt.figure()
	plt.hist(np.log10(rat[rat>0]), nb, alpha = 0.5, label=r'$z='+zf+'$', density = True)
	plt.legend()
	plt.xlabel(r'$\log(\dot{M}_{\mathrm{acc}}/\dot{M}_{\mathrm{Edd}})$')
	plt.ylabel(r'Probability density')
	#plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'BHrat_spec_'+str(sn)+'.pdf')
	plt.close()
	
	plt.figure()
	plt.hist(np.log10(bhdens[bhdens>0]), nb, alpha = 0.5, label=r'$z='+zf+'$', density = True)
	plt.legend()
	plt.xlabel(r'$\log(\rho_{\mathrm{gas}}\ [\mathrm{g\ cm^{-3}}])$')
	plt.ylabel(r'Probability density')
	#plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'BHdens_spec_'+str(sn)+'.pdf')
	plt.close()
	
	plt.figure()
	plt.hist(np.log10(mdot[mdot>0]), nb, alpha = 0.5, label=r'$z='+zf+'$', density = True)
	plt.legend()
	plt.xlabel(r'$\log(\dot{M}_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ yr^{-1}}])$')
	plt.ylabel(r'Probability density')
	#plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'BHmdot_spec_'+str(sn)+'.pdf')
	plt.close()

	plt.figure()
	plt.hist(np.log10(bhmacc[bhmacc>0]), nb, alpha = 0.5, label=r'$z='+zf+'$', density = True)
	plt.legend()
	plt.xlabel(r'$\log(M_{\mathrm{acc}}\ [\mathrm{M}_{\odot}])$')
	plt.ylabel(r'Probability density')
	#plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'BHmass_spec_acc_'+str(sn)+'.pdf')
	plt.close()
	
	plt.figure()
	if mode==0:
		plt.hist(bhm, nb, alpha = 0.5, label=r'$z='+zf+'$', density = True)
		plt.xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M}_{\odot}]$')
	else:
		plt.hist(np.log10(bhm), nb, alpha = 0.5, label=r'$z='+zf+'$', density = True)
		plt.xlabel(r'$\log(M_{\mathrm{BH}}\ [\mathrm{M}_{\odot}])$')
	plt.legend()
	plt.ylabel(r'Probability density')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'BHmass_spec_tot_'+str(sn)+'.pdf')
	plt.close()
	
	m1, m2 = 40, 600
	bins = np.geomspace(m1, m2, nb)
	
	plt.figure()
	plt.hist(np.log10(mass), np.log10(bins), alpha = 0.5, label=r'Hseed', density = True)
	plt.hist(np.log10(massm), np.log10(bins), alpha = 0.5, label=r'Lseed', density = True)
	plt.legend()
	#plt.xlim(m1, m2)
	plt.xlabel(r'$\log(M_{\mathrm{BH, seed}}\ [\mathrm{M}_{\odot}])$')
	plt.ylabel(r'Probability density')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'BHmass_spec_seed_'+str(sn)+'.pdf')
	plt.close()
	
	return bhind

def BH_mass(ds, X = 0.76):
	z = ds['Redshift']
	h = ds['HubbleParam']
	ad = ds.all_data()
	gas = np.sum(ad[('PartType0', 'Masses')].to('Msun'))
	star = np.sum(ad[('PartType4', 'Masses')].to('Msun'))
	keys = ds.field_list
	tag = np.sum([x[0] == 'PartType3' for x in keys])
	if tag>0:
		bh = np.sum(ad[('PartType3', 'Masses')].to('Msun'))
		print('Mtot: ', gas+star+bh)
		sid = ad[('PartType3','ParticleIDs')]
		mass = np.array(ad[('PartType3', 'Masses')].to('Msun'))
		bhmseed = np.array(ad[('PartType3', 'BH Mass')])
		bhm = np.array(ad[('PartType3', 'BH_Mass')])
		bhnseed = np.array(ad[('PartType3', 'BH_NProgs')])
		#bhsound = np.array(ad[('PartType3', 'BH_sound')])
		#print('cs: ', sumy(bhsound**0.5))
		#bhgasv = np.array(ad[('PartType3', 'BH_gasvel')])
		#print('vgas: ', sumy(bhgasv))
		#bhvtheta = np.array(ad[('PartType3', 'BH_vtheta')])
		#print('vtheta: ', sumy(bhvtheta))
		print('Nseed: ', sumy(bhnseed))
		print('MBH/MBH_seed: ', sumy(bhm/bhmseed))
		print('MBH/Mdyn: ', sumy(bhm/mass))
		print('MBH_seed/Mdyn', sumy(bhmseed/mass))
		bhacrb = np.array(ad[('PartType3', 'BH_AccretionLength')])
		print('Hsml', sumy(bhacrb))
		#bhmgas = np.array(ad[('PartType3', 'BH_Ngas')])#/h * UM/Msun
		#print('Ngas', sumy(bhmgas))
		#bhmstar = np.array(ad[('PartType3', 'BH_Nstar')])#/h * UM/Msun
		#print('Ngstar', sumy(bhmstar))
		bhsinf = np.array(ad[('PartType3', 'BH_sigma_inf')])/1e5
		print('sigma_inf', sumy(bhsinf))
		bhrinf = np.array(ad[('PartType3', 'BH_rho_inf')]) * PC**3/Msun #X / PROTON
		print('rho_inf', sumy(bhrinf))
		print('MBH/MBH_seed (total): {}'.format(np.sum(bhm) / np.sum(bhmseed)))
		print('NBH_seed/NBH: {}'.format(np.sum(bhnseed) / bhm.shape[0]))
		return [z, np.sum(bhmseed), np.sum(bhm), sid, bhacrb]
	else:
		return [z, 0, 0]

def details(sn, z, rep = './', ncol = 13, npro = 168, base = 'bhmergers', rep0 = 'blackhole_details/', ext = '.txt'):
	out = []
	for i in range(npro):
		fname = rep+rep0+base+'_'+str(i)+ext
		l = retxt_nor(fname, ncol, 0, 0)
		out.append(l)
	out = np.array(np.hstack(out))
	sel = out[0] <= 1/(1+z)
	#print(np.sum(sel))
	out = out.T[sel]
	totxt(rep+base+'_'+str(sn)+ext, out.T,0,0,0)
	return out.T

nb = 20

def splothis(d, rep, s = 128, alp = 0.7, mode=0, selab=0, ne = nb, lsens=[], lins=[], fage=0, LIGO = []):
	gamma = d['gamma']
	if mode==0:
		rg = d['rg']
		le = d['le0']
		plt.figure()
		plt.hist(le, ne, alpha=alp, density=True)
		plt.xlabel(r'$e$')
		plt.ylabel(r'Probability density')
		plt.yscale('log')
		plt.tight_layout()
		plt.savefig(rep+'ecc_dis.pdf')
		plt.close()
	else:
		fsk = d['sk']
	dz = d['dz']
	if fage>0:
		ldage = d['Dage']
		sel = ldage>0
		print('fraction of age difference: {:.3f}'.format(np.sum(sel)/len(ldage)), sumy(ldage[sel]))
		print('fraction of common parent: {:.3f}'.format(d['com_parent']/len(ldage)))
		plt.figure()
		plt.hist(np.log10(ldage[sel])-6, ne, alpha=alp, density=True)
		plt.xlabel(r'$\log(\Delta t_{\mathrm{BH}}\ [\mathrm{Myr}])$')
		plt.ylabel(r'Probability density')
		plt.yscale('log')
		plt.tight_layout()
		plt.savefig(rep+'Dage_dis.pdf')
		plt.close()
	
	lx, ly, yer = d['Mrate']
	plt.figure()
	plt.errorbar(lx[ly>0], ly[ly>0]*1e9, yerr=yer[ly>0]*1e9, xerr=dz/2, fmt='*')
	plt.xlabel(r'$z_{\mathrm{BBH}}$')
	plt.ylabel(r'$\dot{n}_{\mathrm{BBH}}\ [\mathrm{Gpc^{-3}\ yr^{-1}}]$')
	plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'merger_rate.pdf')
	else:
		plt.savefig(rep+'merger_rate_col.pdf')
	plt.close()

	lzGW, lz, lm, lm1, lm2, le = d['events']
	
	#if mode==0:
	tickloc = np.linspace(np.min(lzGW), np.max(lzGW), 5)[:-1]
	ticklab = ['{:.2f}'.format(np.log10(TZ(z)/1e9/YR)) for z in tickloc]
	plt.figure(figsize=(12,4))
	ax = plt.subplot(131)
	plt.scatter(lzGW, lm, c=lz, cmap = plt.cm.coolwarm, s=s, alpha=alp)
	cb = plt.colorbar()
	cb.ax.set_title('$z_{\mathrm{BBH}}$')#, y=-0.1)
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.ylabel(r'$M_{1}+M_{2}\ [\mathrm{M_{\odot}}]$')
	plt.yscale('log')
	#plt.xlim(xlim)
	#if mode==0:
	ax2 = ax.twiny()
	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(tickloc)
	ax2.set_xticklabels(ticklab)
	ax2.set_xlabel(r'$\log(t\ [\mathrm{Gyr}])$')
	
	ax = plt.subplot(132)
	plt.scatter(lzGW, lm, c=lm2/lm1, cmap = plt.cm.cool, s=s, alpha=alp, label=r'$\gamma={}$'.format(gamma))
	cb = plt.colorbar()
	cb.ax.set_title('$M_{2}/M_{1}$')#, y=-0.11)
	plt.clim(0, 1)
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.yscale('log')
	#plt.xlim(xlim)
	xmax = np.max(lzGW)
	ymax = np.max(lm)
	if (mode!=0) and (fsk>0):
		plt.text(xmax*0.2, ymax*0.9, r'$\gamma={}$'.format(gamma)+r', $f_{\mathrm{fk}}='+str(fsk)+'$')
	else:
		plt.text(xmax*0.7, ymax*0.9, r'$\gamma={}$'.format(gamma))
	#if mode==0:
	ax2 = ax.twiny()
	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(tickloc)
	ax2.set_xticklabels(ticklab)
	ax2.set_xlabel(r'$\log(t\ [\mathrm{Gyr}])$')
	
	ax = plt.subplot(133)
	plt.scatter(lzGW, lm, c=le, s=s, alpha=alp)
	cb = plt.colorbar()
	if mode==0:
		cb.ax.set_title('$e$')#, y=-0.1)
		plt.clim(0, 1)
	else:
		cb.ax.set_title('$a_{0}\ [\mathrm{AU}]$')
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.yscale('log')
	#plt.xlim(xlim)
	#if mode==0:
	ax2 = ax.twiny()
	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(tickloc)
	ax2.set_xticklabels(ticklab)
	ax2.set_xlabel(r'$\log(t\ [\mathrm{Gyr}])$')
	
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'bhmerger_'+str(selab)+'_'+str(rg)+'.pdf')
	else:
		#plt.tight_layout(rect=[0, 0.03, 1, 0.95])
		plt.savefig(rep+'bhmerger_col.pdf')
	#plt.show()
	plt.close()
	
	if len(lsens)>0:
		lfp = d['fpeak']
		lstrain = 2*d['strain']*d['fpeak']
		f1, f2 = 1e-4, 1e3
		y1, y2 = 1e-24, 1e-16
		plt.figure()
		plt.scatter(lfp, lstrain, c=lzGW, s=s, alpha=alp, cmap=plt.cm.cool, zorder=2)
		cb = plt.colorbar()
		cb.ax.set_title(r'$z_{\rm{GW}}$')#, y=-0.11)
		for lf, lh22, m1, m2 in zip(d['lf'], d['lh22'], lm1, lm2):
			if (m1<9e3) and (m2<9e3):
				plt.plot(lf, 2*lf*lh22, alpha=0.2, color='gray', zorder=0)
			elif m2>9e3:
				plt.plot(lf, 2*lf*lh22, alpha=0.4, color='r', zorder=1)
			else:
				plt.plot(lf, 2*lf*lh22, alpha=0.4, color='orange', zorder=1)
		i = 0
		if len(LIGO)>0:
			plt.loglog(LIGO[0], LIGO[1]*LIGO[0]**0.5, 'k-', label='LIGO-Livingston', lw=0.5)
		for sens, ins in zip(lsens,lins):
			plt.plot(sens[0], sens[1]*sens[0]**0.5, label=ins, ls=lls[i], color='k')
			i += 1
		plt.xlabel(r'$f\ [\rm{Hz}]$')
		plt.ylabel(r'Characteristic strain')
		plt.xscale('log')
		plt.yscale('log')
		plt.xlim(f1, f2)
		plt.ylim(y1, y2)
		#plt.legend()
		plt.legend(bbox_to_anchor=(0., 0.79, 1., 1.), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		plt.tight_layout()
		plt.savefig(rep+'GWevents_'+str(selab)+'.pdf')
		plt.close()

lfmt = ['^', 'o', 's', 'D'] * 2

def mplothis(ld, rep, z1 = 4, z2 = 20, nz = nb, alp = 0.7, cs = oi, m1 = 6e2, m2 = 3e3, nm = nb, ne = nb, mode = 0, lref = [], llab=[], Mtot=1e5, lMs=[],dind=0,lins=[]):
	plt.figure()
	i = 0
	y1 = 100
	y2 = 0
	x2 = 0
	for d in ld:
		if mode==2:
			fsk = d['sk']
		else:
			gamma = d['gamma']
		if False:#len(d['snr'])>0:
			lx, ly, yer = d['snr'][dind]['GWrate']
		else:
			lx, ly, yer = d['GWrate']
		dz = d['dz']
		if y1>np.min(ly[ly>0]):
			y1 = np.min(ly[ly>0])
		if y2<np.max(ly):
			y2 = np.max(ly)
		if x2<np.max(lx):
			x2 = np.max(lx)
		if mode==0:
			plt.errorbar(lx[ly>0], ly[ly>0]*1e9, yerr=yer[ly>0]*1e9, xerr=dz/2, fmt=lfmt[i], label=r'$\gamma={}$'.format(gamma), alpha = alp)
		elif mode==1:
			plt.plot(lx[ly>0], ly[ly>0]*1e9, ls=lls[i], label=r'$\gamma={}$'.format(gamma))
			#plt.errorbar(lx[ly>0], ly[ly>0]*1e9, yerr=yer[ly>0]*1e9, ls=lls[i], label=r'$\gamma={}$'.format(gamma))
		else:
			plt.errorbar(lx[ly>0], ly[ly>0]*1e9, yerr=yer[ly>0]*1e9, xerr=dz/2, fmt=lfmt[i], label=r'$f_{\rm{sk}}='+str(fsk)+'$', alpha= alp)
		i += 1
	#xmin = np.min(lx[ly>0])
	#ymax = np.max(ly)
	#plt.text(xmin, ymax*1, r'$\gamma={}$'.format(gamma))
	y1, y2 = min(0.01, y1*1e9*0.5), max(y2*1e9 * 2, 1e3)
	x1, x2 = 0, min(x2*1.2, 10)
	#plt.fill_between([x1, x2], [1, 1], [1e2, 1e2], facecolor='gray', alpha=0.5, label=r'Kinugawa14, $z_{\mathrm{GW}}\simeq 0$')
	if len(lref)>0:
		i = 0
		for ref in lref:
			print('Pop III stellar mass density ({}): {:.0f} Msun Mpc^-3, ratio: {:.3f}'.format(llab[i], lMs[i], lMs[i]/Mtot))
			plt.plot(ref[0], ref[1], color='k', ls=lls[i], label=llab[i])
			plt.plot(ref[0], ref[1]*Mtot/lMs[i], color='k', ls=lls[i], lw=0.5)
			i += 1
	plt.xlim(x1, x2)
	plt.ylim(y1, y2)
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.ylabel(r'$\dot{n}_{\mathrm{GW}}\ [\mathrm{yr^{-1}\ Gpc^{-3}}]$')
	plt.yscale('log')
	plt.legend(loc=2)
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'GW_rate_dens_'+str(selab)+'.pdf')
	else:
		plt.savefig(rep+'GW_rate_dens_col_'+str(selab)+'.pdf')
	plt.close()
	
	zlim = ld[0]['Drate'][1]
	y2 = 0
	plt.figure()
	i = 0
	for d in ld:
		if mode==2:
			fsk = d['sk']
		else:
			gamma = d['gamma']
		dz = d['dz']
		if len(d['snr'])>0:
			lx = d['snr'][dind]['Dzbase']
			ly = d['snr'][dind]['Dlr']
		else:
			lx = d['Dzbase']
			ly = d['Dlr']
		if y2<np.max(ly):
			y2 = np.max(ly)
		if mode==2:
			plt.plot(lx[ly>0], ly[ly>0], ls = lls[i], label=r'$f_{\rm{sk}}='+str(fsk)+'$')
		else:
			plt.plot(lx[ly>0], ly[ly>0], ls = lls[i], label=r'$\gamma={}$'.format(gamma))
		i += 1
	y1, y2 = 0.1, y2 * 2
	x1, x2 = 0, zlim
	#plt.fill_between([x1, x2], [1, 1], [1e2, 1e2], facecolor='gray', alpha=0.5, label=r'Kinugawa14, $z_{\mathrm{GW}}\simeq 0$')
	plt.xlim(x1, x2)
	plt.ylim(y1, y2)
	plt.xlabel(r'$z_{\mathrm{lim}}$')
	if len(d['snr'])>0:
		inst = lins[dind]
		plt.ylabel(r'$\dot{N}_{\mathrm{'+inst+r'}}(z_{\rm{GW}}<z_{\rm{lim}})\ [\mathrm{yr^{-1}}]$')
	else:
		plt.ylabel(r'$\dot{N}_{\mathrm{GW}}(z_{\rm{GW}}<z_{\rm{lim}})\ [\mathrm{yr^{-1}}]$')
	plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'GW_rate_'+str(selab)+'_'+inst+'.pdf')
	else:
		plt.savefig(rep+'GW_rate_col_'+str(selab)+'.pdf')
	plt.close()
		
	bins = np.linspace(z1, z2, nz+1)
	plt.figure()
	i = 0
	for d in ld:
		gamma = d['gamma']
		if len(d['snr'])>0:
			lz = d['snr'][dind]['Dlz']
		else:
			lz = d['Dlz']
		print('z_merge(D): ', sumy(lz))
		plt.hist(lz, bins, alpha = alp*2/3, color=cs[i], label=r'$\gamma={}$'.format(gamma), density=True)
		i += 1
	plt.xlabel(r'$z_{\mathrm{BBH}}$')
	plt.ylabel(r'Probability density')
	#plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'dis_Dzmerge_'+str(selab)+'.pdf')
	else:
		plt.savefig(rep+'dis_Dzmerge_col_'+str(selab)+'.pdf')
	plt.close()
	
	bins = np.linspace(0, zlim, nz+1)
	plt.figure()
	i = 0
	for d in ld:
		gamma = d['gamma']
		if len(d['snr'])>0:
			lz = d['snr'][dind]['DlzGW']
		else:
			lz = d['DlzGW']
		print('z_GW(D): ', sumy(lz))
		plt.hist(lz, bins, alpha = alp*2/3, color=cs[i], label=r'$\gamma={}$'.format(gamma), density=True)
		i += 1
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.ylabel(r'Probability density')
	#plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'dis_DGW_'+str(selab)+'.pdf')
	else:
		plt.savefig(rep+'dis_DGW_col_'+str(selab)+'.pdf')
	plt.close()
	
	m1 = 1e5
	m2 = 0
	for d in ld:
		if len(d['snr'])>0:
			lm = d['snr'][dind]['Dlm']
		else:
			lm = d['Dlm']
		ml, mh = np.min(lm), np.max(lm)
		if ml<m1:
			m1 = ml
		if mh>m2:
			m2 = mh
	bins = np.log10(np.geomspace(m1, m2, nm+1))
	plt.figure()
	i = 0
	for d in ld:
		gamma = d['gamma']
		if len(d['snr'])>0:
			lm = d['snr'][dind]['Dlm']
		else:
			lm = d['Dlm']
		print('M(D): ', sumy(lm))
		plt.hist(np.log10(lm), bins, alpha = alp*2/3, color=cs[i], label=r'$\gamma={}$'.format(gamma), density=True)
		i += 1
	plt.xlabel(r'$\log(M\ [\mathrm{M}_{\odot}])$')
	plt.ylabel(r'Probability density')
	plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep+'dis_DM_'+str(selab)+'.pdf')
	else:
		plt.savefig(rep+'dis_DM_col_'+str(selab)+'.pdf')
	plt.close()
	
	if mode==0:
		bins = np.linspace(0, 1, ne+1)
		plt.figure()
		i = 0
		for d in ld:
			gamma = d['gamma']
			if len(d['snr'])>0:
				le = d['snr'][dind]['Dle']
			else:
				le = d['Dle']
			print('ecc(D): ', sumy(le))
			plt.hist(le, bins, alpha = alp*2/3, color=cs[i], label=r'$\gamma={}$'.format(gamma), density=True)
			i += 1
		plt.xlabel(r'$e$')
		plt.ylabel(r'Probability density')
		plt.yscale('log')
		plt.legend()
		plt.tight_layout()
		plt.savefig(rep+'dis_Decc_'+str(selab)+'.pdf')
		plt.close()
		
		bins = np.linspace(0, 1, ne+1)
		plt.figure()
		i = 0
		for d in ld:
			gamma = d['gamma']
			if len(d['snr'])>0:
				lq = d['snr'][dind]['Dlq']
			else:
				lq = d['Dlq']
			print('q(D): ', sumy(lq))
			plt.hist(lq, bins, alpha = alp*2/3, color=cs[i], label=r'$\gamma={}$'.format(gamma), density=True)
			i += 1
		plt.xlabel(r'$q\equiv M_{2}/M_{1}$')
		plt.ylabel(r'Probability density')
		#plt.yscale('log')
		plt.legend()
		plt.tight_layout()
		plt.savefig(rep+'dis_Dq_'+str(selab)+'.pdf')
		plt.close()
		
	if len(ld[0]['snr'])>0:
		for i in range(len(ld[0]['snr'])):
			inst = lins[i]
			plt.figure()
			bins = np.linspace(0, 5, ne+1)
			j = 0
			for d in ld:
				gamma = d['gamma']
				lsnr = d['snr'][i]['snr']
				lsnr = np.log10(lsnr[lsnr>0])
				plt.hist(lsnr, bins, alpha = alp*2/3, color=cs[j], label=r'$\gamma={}$'.format(gamma), density=True)
				j += 1
			plt.xlabel(r'$\log(\rm{SNR})$')
			plt.ylabel(r'Probability density')
			plt.yscale('log')
			plt.legend()
			plt.tight_layout()
			plt.savefig(rep+'dis_DSNR_'+inst+'_'+str(selab)+'.pdf')
			plt.close()
				
			
def mergerhis(l0, ds, mode = 0, gamma = 1.5, H = 17.5, nth = 1e1, beta = 2.0, X = 0.76, bfrac=1, rg = 0, zlim = 5, dz = 0.5, mrel = ms_mbh_default, zsn=4, lsens=[], snrlim=7, ffac=0.9, fage=0, V=V):
	h = ds['HubbleParam']
	lz = 1/l0[0] - 1
	ncol = len(l0)
	sel = lz >= zsn
	l = l0.T[sel]
	Nbbh = len(l)
	l = l.T
	lm1 = l[2]*UM/h/Msun
	lm2 = l[7]*UM/h/Msun
	le = l[12]
	out = {}
	out['Nbbh'] = Nbbh
	out['le0'] = le
	out['gamma'] = gamma
	out['rg'] = rg
	out['dz'] = dz
	out['snr'] = []
	
	if fage>0:
		lz1 = 1/l[ncol-2]-1
		lz2 = 1/l[ncol-1]-1
		lid1 = l[ncol-4]
		lid2 = l[ncol-3]
		out['com_parent'] = sum(lid1==lid2)
		Dage = np.array([abs(TZ(z1)-TZ(z2))/YR for z1, z2 in zip(lz1, lz2)])
		out['Dage'] = Dage
		
	lz = 1/l[0] - 1
	tlim = TZ(0)/1e9/YR
	lm = lm1 + lm2
	
	if len(lz)==0:
		out['Mrate'] = [[], [], []]
	else:	
		#z1, z2 = np.min(lz), np.max(lz)
		z1 = max(np.min(lz)-dz/2, 0)
		z2 = max(np.max(lz)+dz/2, z1 + dz)
		nzb = min(int((z2-z1)/dz), 20)
		if nzb<=0:
			nzb = 1
		zbase = np.linspace(z1, z2, nzb+1)
		dt = np.array([TZ(zbase[i]) - TZ(zbase[i+1]) for i in range(nzb)])/YR
		ly, lx0 = np.histogram(lz, zbase)
		yer = ly**0.5/dt/V
		ly = ly/dt/V
		lx = (lx0[1:] + lx0[:-1])/2
		out['Mrate'] = [lx, ly, yer]
	
	if mode==0:
		lt = (np.array([TZ(z) for z in lz])/YR + np.array([tGW(m1+m2, m1, m2, e, H, gamma) for m1, m2, e in zip(lm1, lm2, le)]))/1e9
	else:
		lsig0 = sigma_inf(lm)
		lrho0 = rho_inf(lm, gamma, mrel)
		"""
		lrho1 = rho_inf(lm, gamma, mrel, rmax)
		lrho2 = rho_inf(lm, gamma, mrel, rmin)
		lrhou, lrhod = lrho1, lrho2
		if np.sum(lrhou)<np.sum(lrhod):
			lrhou, lrhod = lrho2, lrho1
		#print(sumy(lrhou/lrhod))
		rhoth = nth * PROTON / X
		lb = (l[14]/rhoth)**beta * l[14] * (l[14]>rhoth) + l[14] * (l[14]<=rhoth)
		lb = lb * (lb<=lrhou) + lrhou * (lb>lrhou)
		lb += (lrhod-lb) * (lb<lrhod)
		lrho = lb + lrhod * (lb<=0)
		"""
		lrho = rho_inf(lm, gamma, mrel, bfrac)
		lsig = l[13] + lsig0 * (l[13]<=0)
		lt = (np.array([TZ(z) for z in lz])/YR + np.array([coaltime(m1, m2, sig, rho, e, H) for m1, m2, sig, rho, e in zip(lm1, lm2, lsig, lrho, le)]))/1e9
	sel = lt<tlim
	#print(lt, tlim, le)

	"""
	if mode!=0:
		rat0 = lsig/lsig0
		sel0 = rat0>=0
		rat1 = lrho/lrho0
		sel1 = rat1>=0
		#print('sigma_inf: ', l[13], sumy(l[13]))
		#print('rho_inf: ', lb, sumy(l[14]))
		#print('sig_ratio: ', sumy(rat0[sel0]))
		#print('rho_ratio: ', sumy(rat1[sel1]))
		plt.figure()
		plt.scatter(lm[sel0], rat0[sel0], c=lm2[sel0]/lm1[sel0], label=r'$\sigma_{\mathrm{sim}}/\sigma_{\mathrm{est}}$', s=128, alpha=0.7, marker='o')
		plt.scatter(lm[sel1], rat1[sel1], c=np.log10(lm2/lm1)[sel1], label=r'$\rho_{\mathrm{sim}}/\rho_{\mathrm{est}}$', s=128, alpha=0.7, marker='^')
		cb = plt.colorbar()
		cb.ax.set_title('$M_{2}/M_{1}$')#, y=-0.11)
		plt.yscale('log')
		plt.xscale('log')
		plt.ylim(1e-3, 1e3)
		plt.xlabel(r'$M_{1}+M_{2}\ [\mathrm{M_{\odot}}]$')
		plt.ylabel(r'Ratio')
		plt.legend()
		plt.tight_layout()
		plt.savefig(rep+'ratio_gamma'+str(gamma)+'.pdf')
		plt.close()
	"""

	d = np.array([l[1], lt, lz, lm1, lm2, le, np.ones(len(lz))]).T
	d = d[sel]
	if rg!=0:
		d = np.array(sorted(d, key=lambda x:x[0]))
		nme = len(d)
		for i in range(0, nme):
			ID = d[i][0]
			t = d[i][1]
			z = d[i][2]
			if d[i][6]<=0:
				continue
			for j in range(i+1, nme):
				if ID != d[j][0]:
					break
				if ((d[j][1]>t) or (d[j][6]<=0)) and (d[j][2]>z):
					temp = d[i][3]
					d[i][3] -= d[j][4]
					d[i][1] = (TZ(z)/YR + tGW(temp+d[i][4], d[i][3], d[i][4], d[i][5], H, gamma))/1e9
					t = d[i][1]
					if rg > 1:
						if d[i][4] > d[j][4] * (rg-1):
							d[j][6] = 0
					d[j][3] = temp
					d[j][1] = (TZ(d[j][2])/YR + tGW(temp+d[j][4], d[j][3], d[j][4], d[j][5], H, gamma))/1e9
		ltag = d[:,6]
		d = d[ltag>0]
	lID, lt, lz, lm1, lm2, le, ltag = d.T
	sel = lt<tlim
	d = d[sel]
	lID, lt, lz, lm1, lm2, le, ltag = d.T
	lm = lm1 + lm2
	lzGW = np.array([ZT(t) for t in np.log10(lt)])
	Nmg = len(lm)
	out['Nmg'] = Nmg
			
	#Nbbh = d['Nbbh']
	if Nmg>0:
		out['Nhh'] = np.sum((lm1>9e3)*(lm2>9e3))
		out['Nhl'] = np.sum((lm1>9e3)*(lm2<9e3))
		out['Nll'] = Nmg-out['Nhh']-out['Nhl']
	else:
		out['Nhh'] = out['Nhl'] = out['Nll'] = 0
			
	if len(lzGW)==0:
		out['Drate'] = [0, zlim]
	else:
		out['strain'] = np.array([PMwf(np.array([fRD_a(afin(0,0,m1,m2),m1+m2)*ffac]), m1, m2, DZ(z)*(1+z))[0][0] for m1, m2, z in zip(lm1, lm2, lzGW)])
		out['fpeak'] = np.array([fRD_a(afin(0,0,m1,m2), m1+m2)*ffac/(1+z) for m1, m2, z in zip(lm1, lm2, lzGW)])
		lf = np.geomspace(1e-4, 1e3, 100)
		out['lh22'] = np.array([PMwf(lf[lf*(1+z)<2*fRD_a(afin(0,0,m1,m2),m1+m2)]*(1+z), m1, m2, DZ(z)*(1+z))[0] for m1, m2, z in zip(lm1, lm2, lzGW)])
		out['lf'] = np.array([lf[lf*(1+z)<2*fRD_a(afin(0,0,m1,m2),m1+m2)] for m1, m2, z in zip(lm1, lm2, lzGW)])
		#print(len(lzGW))
		z1 = max(np.min(lzGW)-dz/2, 0)
		z2 = max(np.max(lzGW)+dz/2, z1 + dz)
		nzb = int((z2-z1)/dz)#, 30)
		if nzb<=0:
			nzb = 1
		zbase = np.linspace(z1, z2, nzb+1)
		#t2, t1 = np.max(lt), np.min(lt)
		#nzb = int((t2-t1)*1e9/ts)
		#if nzb<=0:
		#	nzb = 1
		#ltedge = np.linspace(min(t2+ts/1e15,tlim), t1-ts/1e15, nzb+1)
		#zbase = np.array([ZT(t) for t in np.log10(ltedge)])
		dt = np.array([TZ(zbase[i]) - TZ(zbase[i+1]) for i in range(nzb)])/YR
		ly, lx0 = np.histogram(lzGW, zbase)
		yer = ly**0.5/dt/V
		ly = ly/dt/V
		lx = (lx0[1:] + lx0[:-1])/2
		out['GWrate'] = [lx, ly, yer]
	
		lV = np.array([DZ(x2)**3 - DZ(x1)**3 for x2, x1 in zip(lx0[1:], lx0[:-1])])*4*np.pi/3/MPC**3
		intg = ly/(1+lx) * lV
		lrate = np.cumsum(intg[lx<=zlim])
		if len(lrate)>0:
			rate = lrate[-1]#np.sum(ly[lx<zlim] / (1+lx[lx<zlim]) * lV[lx<zlim])
			out['Dle'] = le[lzGW<zlim]
			out['Drate'] = [rate, zlim]
			out['Dlz'] = lz[lzGW<zlim]
			out['DlzGW'] = lzGW[lzGW<zlim]
			out['Dzbase'] = lx[lx<=zlim]
			out['Dlr'] = lrate
			out['Dlm'] = lm[lzGW<zlim]
			lq = lm2[lzGW<zlim]/lm1[lzGW<zlim]
			lq = lq * (lq<=1) + 1/lq * (lq>1)
			out['Dlq'] = lq
		else:
			print(lrate, lx, zlim, intg[lx<=zlim], intg)
			out['Drate'] = [0, zlim]
	
	out['events'] = [lzGW, lz, lm, lm1, lm2, le]
	
	lsnr = []
	if len(lsens)>0:
		dd0 = np.array([lzGW, lm1, lm2, lz, le]).T
		for sens in lsens:
			lzGW, lm1, lm2, lz, le = dd0.T
			snr = np.array([SNR(m1, m2, sens, z) for m1, m2, z in zip(lm1, lm2, lzGW)])
			sel = snr>=snrlim
			dd = dd0[sel].T
			dsnr = {}
			dsnr['snr'] = snr
			lzGW, lm1, lm2, lz, le = dd
			if len(lzGW)<=0:
				dsnr['Drate'] = [0, zlim]
				dsnr['frac'] = 0
				lsnr.append(dsnr)
				continue
			lm = lm1 + lm2
			z1 = max(np.min(lzGW)-dz/2, 0)
			z2 = max(np.max(lzGW)+dz/2, z1 + dz)
			nzb = int((z2-z1)/dz)#, 30)
			if nzb<=0:
				nzb = 1
			zbase = np.linspace(z1, z2, nzb+1)
			dt = np.array([TZ(zbase[i]) - TZ(zbase[i+1]) for i in range(nzb)])/YR
			ly, lx0 = np.histogram(lzGW, zbase)
			yer = ly**0.5/dt/V
			ly = ly/dt/V
			lx = (lx0[1:] + lx0[:-1])/2
			dsnr['GWrate'] = [lx, ly, yer]
			lV = np.array([DZ(x2)**3 - DZ(x1)**3 for x2, x1 in zip(lx0[1:], lx0[:-1])])*4*np.pi/3/MPC**3
			intg = ly/(1+lx) * lV
			lrate = np.cumsum(intg)
			rate = lrate[-1]
			Nmg_ = len(lzGW)
			dsnr['Nmg'] = Nmg_
			dsnr['frac'] = Nmg_/Nmg
			dsnr['Dle'] = le
			dsnr['Drate'] = [rate, z2]
			dsnr['Dlz'] = lz
			dsnr['DlzGW'] = lzGW
			dsnr['Dzbase'] = lx
			dsnr['Dlr'] = lrate
			dsnr['Dlm'] = lm
			lq = lm2/lm1
			lq = lq * (lq<=1) + 1/lq * (lq>1)
			dsnr['Dlq'] = lq
			lsnr.append(dsnr)
		out['snr'] = lsnr
	return out

def stellarb(sn, rep = './', base = 'snapshot', ext = '.hdf5', gamma = 1.5, sk = 0.1, zlim = 5, dz = 0.1, frm = 0.57, seed = 1234, mode=0):#, fbh = 0.8):
	np.random.seed(seed)
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	ad = ds.all_data()
	ms = np.array(ad[('PartType3','Masses')].to('Msun'))
	lz = 1/np.array(ad[('PartType3','StellarFormationTime')]) - 1
	sel = np.array(ad[('PartType3','BH Type')])==1
	ms = ms[sel]
	lz = lz[sel]
	#mmaxbh = np.array(ad[('PartType3', 'BH Max mass')][sel])
	#minput = ms*fbh
	out = []
	for m, z in zip(ms, lz):
		if sk>0:
			l = binary_gen(m, gamma=gamma, sk=sk, frm=frm, adis=0, mode=mode)
		else:
			l = binary_gen(m, gamma=gamma, frm=frm, adis=1, mode=mode)
		out.append([z*np.ones(len(l[0])), *l])
	d = np.hstack(out)
	out = {}
	lz, lt, lm1, lm2, la0, le = d
	lm = lm1 + lm2
	ltGW = (lt + np.array([TZ(z) for z in lz])/YR)/1e9
	tlim = TZ(0)/1e9/YR
	sel = ltGW < tlim
	d = np.array([ltGW, lz, lm, lm1, lm2, la0]).T[sel]
	ltGW, lz, lm, lm1, lm2, la0 = d.T
	lzGW = np.array([ZT(t) for t in np.log10(ltGW)])
	out['events'] = [lzGW, lz, lm, lm1, lm2, la0]
	Nmg = len(lzGW)
	out['Nmg'] = Nmg
	out['gamma'] = gamma
	out['sk'] = sk
	out['dz'] = dz
	
	if len(lzGW)==0:
		out['Drate'] = [0, zlim]
	else:
		z1 = max(np.min(lzGW)-dz/2, 0)
		z2 = max(np.max(lzGW), z1 + dz)
		nzb = min(int((z2-z1)/dz), 300)
		if nzb<=0:
			nzb = 1
		zbase = np.linspace(z1, z2, nzb+1)
		dt = np.array([TZ(zbase[i]) - TZ(zbase[i+1]) for i in range(nzb)])/YR
		ly, lx0 = np.histogram(lzGW, zbase)
		yer = ly**0.5/dt/V
		ly = ly/dt/V
		lx = (lx0[1:] + lx0[:-1])/2
		out['GWrate'] = [lx, ly, yer]
		
		lV = np.array([DZ(x2)**3 - DZ(x1)**3 for x2, x1 in zip(lx0[1:], lx0[:-1])])*4*np.pi/3/MPC**3
		intg = ly/(1+lx) * lV
		lrate = np.cumsum(intg[lx<=zlim])
		rate = lrate[-1]#np.sum(ly[lx<zlim] / (1+lx[lx<zlim]) * lV[lx<zlim])
		out['Dla0'] = la0[lzGW<zlim]
		out['Drate'] = [rate, zlim]
		out['Dlz'] = lz[lzGW<zlim]
		out['DlzGW'] = lzGW[lzGW<zlim]
		out['Dzbase'] = lx[lx<=zlim]
		out['Dlr'] = lrate
		out['Dlm'] = lm[lzGW<zlim]

	if len(lz)==0:
		out['Mrate'] = [[], [], []]
	else:		
		z1 = max(np.min(lzGW)-dz/2, 0)
		z2 = max(np.max(lzGW), z1 + dz)
		nzb = min(int((z2-z1)/dz), 300)
		if nzb<=0:
			nzb = 1
		zbase = np.linspace(z1, z2, nzb+1)
		dt = np.array([TZ(zbase[i]) - TZ(zbase[i+1]) for i in range(nzb)])/YR
		ly, lx0 = np.histogram(lz, zbase)
		yer = ly**0.5/dt/V
		ly = ly/dt/V
		lx = (lx0[1:] + lx0[:-1])/2
		out['Mrate'] = [lx, ly, yer]
	return out

def bhtrack(sni, snf, rep = './', base = 'snapshot', ext = '.hdf5', sp = 1):
	ds0 = yt.load(rep+base+'_'+str(sni).zfill(3)+ext)
	ad0 = ds0.all_data()
	bhm = np.array(ad0[('PartType3', 'BH_Mass')])
	bhmacc = np.array(ad0[('PartType3', 'BH_Macc')])
	bhid = ad0[('PartType3','ParticleIDs')]
	maxm = np.max(bhm)
	sel = bhm==maxm
	ID = bhid[sel][0]
	lz = []
	lmbh = []
	lmacc = []
	lz.append(ds0['Redshift'])
	lmbh.append(maxm)
	lmacc.append(bhmacc[sel][0])
	sn = sni - sp
	while sn>=snf:
		ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
		keys = ds.field_list
		tag = np.sum([x[0] == 'PartType3' for x in keys])
		if tag>0:
			ad = ds.all_data()
			bhid = ad[('PartType3','ParticleIDs')]
			sel = bhid==ID
			if np.sum(sel)<=0:
				print('Lose track at z = {:.3f}, sn = {}'.format(ds['Redshift'],sn))
				break
			bhm = np.array(ad[('PartType3', 'BH_Mass')])[sel][0]
			bhmacc = np.array(ad[('PartType3', 'BH_Macc')])[sel][0]
			lz.append(ds['Redshift'])
			lmbh.append(bhm)
			lmacc.append(bhmacc)
			sn = sn - sp
		else:
			print('Lose track at z = {:.3f}, sn = {}'.format(ds['Redshift'],sn))
			break
	return [lz, lmbh, lmacc]

def bhhis(l, sni, rep = './', base = 'snapshot', ext = '.hdf5'):
	ds0 = yt.load(rep+base+'_'+str(sni).zfill(3)+ext)
	h = ds0['HubbleParam']
	ad0 = ds0.all_data()
	bhm = np.array(ad0[('PartType3', 'BH_Mass')])
	bhm0 = np.array(ad0[('PartType3', 'BH Mass')])
	bhmacc = np.array(ad0[('PartType3', 'BH_Macc')])
	bhid = ad0[('PartType3','ParticleIDs')]
	mdyn = np.array(ad0['PartType3', 'Masses'].to('Msun'))
	maxm = np.max(bhm)
	sel = bhm==maxm
	mdynmax = mdyn[sel][0]
	print('MBH_max = {} Msun (Mdyn = {} Msun)'.format(maxm, mdynmax))
	ID = bhid[sel][0]
	mseed = bhm0[sel][0]
	sel = l[1]==ID
	d = l.T[sel]
	d = np.array(sorted(d, key=lambda x:x[0]))
	d = d.T
	lz = 1/d[0]-1
	lm1 = d[2]*UM/h/Msun
	lm2 = d[7]*UM/h/Msun
	lm = lm1 + lm2
	lmm = np.cumsum(lm2)
	lm2_ = np.hstack([[0], lm2[:-1]])
	lm += (lm - lmm < mseed) * lm2_
	lmacc = lm - lmm - mseed #lm1[1:] - (lm1[:-1]+lm2[:-1])
	#print(d[0][1:]-d[0][:-1])
	
	plt.figure()
	plt.plot(lz, lm, label='Total')
	plt.plot(lz, lmm, '--', label='Merge', marker='.')
	plt.plot(lz, lmacc, '-.', label='Accreted')
	plt.yscale('log')
	plt.ylim(1e2, np.max(lm)*1.2)
	plt.legend()
	plt.xlabel(r'$z$')
	plt.ylabel(r'$M\ [\mathrm{M_{\odot}}]$')
	plt.tight_layout()
	plt.savefig(rep+'Maxhis.pdf')
	plt.close()
	return ID
	
#Mup = lambda z: 2.5e7*((1+z)/10)**-1.5
#Mdown = lambda z: 1.54e5*((1+z)/31)**-2.074 #1e6*((1+z)/10)**-2
	
def group_info(ds, obj, rep = './', mode = 0, mfac = 1.0, rfac = 0.2, mlab = 0):
	z = ds['Redshift']
	if mlab==0:
		Mmin = Mdown(z)
	else:
		Mmin = Mup(z)
	lh = obj.halos
	datah = []
	for halo in lh:
		M = halo.masses['total'].to('Msun')*mfac
		if M<Mmin:
			continue
		Ms_ = halo.masses['stellar'].to('Msun')
		Mg = halo.masses['gas'].to('Msun')
		if len(halo.galaxy_index_list)>0:
			gal = halo.central_galaxy
			pos = gal.pos.to('kpccm/h')
			Ms = gal.masses['stellar'].to('Msun')
			R = gal.radii['baryon_half_mass'].to('kpccm/h')
		else:
			pos = halo.pos.to('kpccm/h')
			R = halo.radii['baryon_half_mass'].to('kpccm/h')
			if R<=0:
				R = halo.radii['total_half_mass'].to('kpccm/h')
			Ms = Ms_
		ad = ds.sphere(pos, (R*rfac, 'kpccm/h'))
		lbh = ad[('PartType3', 'BH_Mass')]
		if len(lbh)>0:
			if mode==0:
				Mbh = np.max(lbh)
			else:
				Mbh = np.sum(lbh)
		else:
			Mbh = 0
		datah.append([Mbh, Ms, R.to('kpc'), M, Ms_, Mg])
	datah = np.array(datah).T
	d = {}
	d['z'] = z
	d['data'] = datah
	return d
	
def gal_info(ds, obj, rep = './', mode = 0, mfac = 1.0, rfac = 0.2, Mmin = 32*586):
	z = ds['Redshift']
	lg = obj.galaxies
	datag = []
	for gal in lg:
		Ms = gal.masses['stellar'].to('Msun')*mfac
		if Ms<Mmin:
			continue
		Mg = gal.masses['gas'].to('Msun')
		pos = gal.pos.to('kpccm/h')
		R = gal.radii['baryon_half_mass'].to('kpccm/h')
		ad = ds.sphere(pos, (R*rfac, 'kpccm/h'))
		lbh = ad[('PartType3', 'BH_Mass')]
		if len(lbh)>0:
			if mode==0:
				Mbh = np.max(lbh)
			else:
				Mbh = np.sum(lbh)
		else:
			Mbh = 0
		datag.append([Mbh, Ms, R.to('kpc'), Mg])
	datag = np.array(datag).T
	d = {}
	d['z'] = z
	d['data'] = datag
	return d

#Mbh_Ms0 = lambda x: 10**(1.25*np.log10(x)-4.8)
Mbh_Ms1 = lambda x: 10**(1.4*np.log10(x)-6.45)
Mbh_Ms2 = lambda x: 10**(1.05*np.log10(x)-4.1)
	
def linear(x, b0, b1):
	return b0 + x*b1
	
def get_range(l, fac = 0.7):
	return np.min(l)*fac, np.max(l)/fac
	
def halo_star(d, rep, sn, s=128, Ob = 0.048, Om=0.315, Mmin = 32*586):
	z = d['z']
	mbh = d['data'][0]
	ms = d['data'][4]
	mh = d['data'][3]
	mg = d['data'][5]
	sel = ms>Mmin
	ngal = np.sum(sel)
	nh = len(mh)
	occ = ngal/nh
	lab = r'$z={:.2f}$'.format(z)+r', $f_{\mathrm{occ}}'+r'={:.3f}'.format(occ)+r'$'

	ms = ms[sel]
	mh = mh[sel]
	mbh = mbh[sel]
	mg = mg[sel]
	mb = ms + mg
	#print('Ms/Mh: ', sumy(ms/mh))
	print('Mb/Mh: ', sumy(mb/mh))
	print('Ms/Mb: ', sumy(ms/mb))
	x1, x2 = get_range(mh)
	y1, y2 = get_range(ms)
	xr = np.array([x1, x2])
	plt.figure()
	plt.scatter(mh, ms, s=s, c = np.log10(mb/mh*Om/Ob), cmap=plt.cm.cool, label=lab)
	cb = plt.colorbar()
	cb.ax.set_title(r'$\log\left(\frac{M_{\rm{baryon}}}{f_{\rm{baryon}}M_{\rm{halo}}}\right)$')#, y=-0.11)
	plt.plot(xr, xr*Ob/Om, 'k-', label=r'$M_{\mathrm{\star}}/M_{\mathrm{halo}}=f_{\mathrm{baryon}}$')
	plt.plot(xr, xr*Ob/Om*0.1, 'k--', label=r'$M_{\mathrm{\star}}/M_{\mathrm{halo}}=0.1f_{\mathrm{baryon}}$')
	plt.plot(xr, xr*Ob/Om*0.01, 'k-.', label=r'$M_{\mathrm{\star}}/M_{\mathrm{halo}}=0.01f_{\mathrm{baryon}}$')
	plt.plot(xr, xr*Ob/Om*0.001, 'k:', label=r'$M_{\mathrm{\star}}/M_{\mathrm{halo}}=10^{-3}f_{\mathrm{baryon}}$')
	plt.xlim(xr)
	plt.ylim(y1, y2)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M_{\mathrm{halo}}\ [\mathrm{M}_{\odot}]$')
	plt.ylabel(r'$M_{\star}\ [\mathrm{M}_{\odot}]$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'Ms_Mh_'+str(sn)+'.pdf')
	plt.close()	
	
def plotrela(d, rep, mode, base, sn, s = 128, Mmin = 32*586, y1=40, y2=1e5, f1=1e5, f2=1e8, alp=0.7):
	z = d['z']
	mbh = d['data'][0]
	ms = d['data'][1]
	lr = d['data'][2]
	sel = (mbh>0) * (ms>Mmin)
	nbh = np.sum(sel)
	ngal = np.sum(ms>Mmin)
	occ = nbh/ngal
	lab = r'$z={:.2f}$'.format(z)+r', $f_{\mathrm{occ}}'+r'={:.3f}'.format(occ)+r'$'
	
	ms = ms[sel]
	mbh = mbh[sel]
	lr = lr[sel]
	x1, x2 = get_range(ms)
	#y1, y2 = get_range(mbh)
	xr = np.array([x1, x2])
	#sel = ms>0
	#sel = (ms>f1) * (ms<f2)
	fit = curve_fit(linear, np.log10(ms), np.log10(mbh))
	b0, b1 = fit[0]
	if b0>=0:
		flab = r'$\log M_{\mathrm{BH}}='+'{:.2f}'.format(b1)+r'\log M_{\star}'+r'+{:.2f}$'.format(b0)
	else:
		flab = r'$\log M_{\mathrm{BH}}='+'{:.2f}'.format(b1)+r'\log M_{\star}'+r'-{:.2f}$'.format(-b0)
	perr = np.sqrt(np.diag(fit[1]))
	#print(perr)
	Mbh_Ms0 = lambda x: 10**(np.log10(x)*b1 + b0)
	Mbh_Msu = lambda x: 10**(np.log10(x)*(b1+perr[1]) + b0+perr[0])
	Mbh_Msd = lambda x: 10**(np.log10(x)*(b1-perr[1]) + b0-perr[0])
	plt.figure()
	plt.scatter(ms, mbh, s=s, c = np.log10(lr), cmap=plt.cm.cool, label=lab, alpha=alp)
	cb = plt.colorbar()
	cb.ax.set_title('$\log(R_{1/2}\ [\mathrm{kpc}])$')#, y=-0.11)
	plt.plot(xr, Mbh_Ms0(xr), 'k-', label=flab)
	plt.fill_between(xr, Mbh_Msd(xr), Mbh_Msu(xr), facecolor='gray', alpha=0.5, zorder=0)
	plt.plot(xr, Mbh_Ms1(xr), 'k--', label=r'$\log M_{\mathrm{BH}}=1.4\log M_{\star}-6.45$')
	plt.plot(xr, Mbh_Ms2(xr), 'k:', label=r'$\log M_{\mathrm{BH}}=1.05\log M_{\star}-4.1$')
	plt.xlim(xr)
	plt.ylim(y1, y2)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M_{\star}\ [\mathrm{M}_{\odot}]$')
	plt.ylabel(r'$M_{\mathrm{BH}}\ [\mathrm{M}_{\odot}]$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'Mbh_Ms_'+str(sn)+'_'+base+str(mode)+'.pdf')
	plt.close()
	
if __name__ == "__main__":
	sn = 23
	npro = 168
	#npro = 1152
	#npro = 672
	#npro = 408
	#rlab = 'NSFDBKzoom'
	#rlab = 'FDzoom_HR'
	rlab = 'FDzoom'
	#rlab = 'FDbox'
	#rep = 'NSFDBKzoom_Hseed/'
	#rep = 'NSFDBKzoom_Lseed/'
	#rep = 'FDzoom_NBH/'
	rep = 'FDzoom_Lseed/'
	#rep = 'FDzoom_Lseed_comp/'
	#rep = 'FDzoom_Lseed_ADC/'#NDF/'#local/'
	#rep = 'FDzoom_Lseed_HR/'
	#rep = 'FDzoom_Hseed/'
	#rep = 'SFtest/'
	#rep = 'LRbox_FDWIND_Hseed/'
	#rep = 'FDbox_Lseed/'
	cr = 'comp/'
	gamma = 1.5
	
	zlim = 10
	snrlim = 10
	dz = 0.2
	selab = 1
	glab = 0
	ins = 2
	bfrac = 1.0
	fage = 1
	ncol = 19
	nr = 20+1
	
	fitlab = 1

	Mh_Ms_Hseed = lambda x: 10**((np.log10(x)-1.45)/0.24)
	Mh_Ms_Lseed = lambda x: 10**((np.log10(x)-1)/0.22)
	if selab==0:
		rf = ms_mbh_default
	else:
		rf = Mh_Ms_Lseed

	mode = int(sys.argv[1])
	if mode==-1:
		data = loaddata(sn, sn, rep)
		ds0 = data['ds']
		#plt.figure()
		#ax = new_3Dax(111)
		#dis3D(data, 0, ax)
		#plt.tight_layout()
		#plt.show()
		#plt.close()
	else:
		ds0 = yt.load(rep+'snapshot_'+str(sn).zfill(3)+'.hdf5')
	zsn = ds0['Redshift']
	h = ds0['HubbleParam']
	#print(SFgas(sn, rep))
	#print(ds0.field_list)
	xlm = [25, 4]

	#if mode==2 or mode==0:
	M1, M2 = smass(ds0)
	print('Snapshot: M(PopIII) = {:.0f} Msun Mpc^-3, M(PopII) = {:.0f} Msun Mpc^-3'.format(M1/V, M2/V))
	MpopIII = M1

	sfr0III = np.array(retxt_nor(rep+'popIII_sfr.txt', 4, 0, 0))
	sfr1III = np.array(retxt_nor(cr+'popIII_sfr.txt', 4, 0, 0))
	sfr0II = np.array(retxt_nor(rep+'popII_sfr.txt', 4, 0, 0))
	sfr1II = np.array(retxt_nor(cr+'popII_sfr.txt', 4, 0, 0))

	lm1 = sfr0III[1] * sfr0III[2]
	lm2 = sfr0II[1] * sfr0II[2]
	M1 = np.sum(lm1[sfr0III[0]<=1/(1+zsn)])/V
	M2 = np.sum(lm2[sfr0II[0]<=1/(1+zsn)])/V
	print('OTF: M(PopIII) = {:.0f} Msun Mpc^-3, M(PopII) = {:.0f} Msun Mpc^-3'.format(M1, M2))
	lm1 = sfr1III[1] * sfr1III[2]
	lm2 = sfr1II[1] * sfr1II[2]
	M1 = np.sum(lm1[sfr1III[0]<=1/(1+zsn)])/(4/0.6774)**3/fac0
	M2 = np.sum(lm2[sfr1II[0]<=1/(1+zsn)])/(4/0.6774)**3/fac0
	print('Old: M(PopIII) = {:.0f} Msun Mpc^-3, M(PopII) = {:.0f} Msun Mpc^-3'.format(M1, M2))
	
	ref0 = np.array(retxt(cr+'SF16.txt', 3, 0, 0))
	ref1 = np.array(retxt(cr+'edges_bottom.txt', 2, 0, 0))
	ref2 = np.array(retxt(cr+'edges_top.txt', 2, 0, 0))
	ref3 = np.array(retxt(cr+'sf18_sfrd.txt', 3, 0, 0)) 
	ref4 = np.array(retxt(cr+'sf18_sfrd_evo.txt', 3, 0, 0)) 
	ref5 = np.array(retxt(cr+'popIII_sfr_s2_skx.txt', 4, 0, 0))
	ref6 = np.array(retxt(cr+'popII_sfr_s2_skx.txt', 4, 0, 0))

	lm1 = ref5[1] * ref5[2]
	lm2 = ref6[1] * ref6[2]
	M1 = np.sum(lm1[ref5[0]<=1/(1+zsn)])/(4/0.6774)**3
	M2 = np.sum(lm2[ref6[0]<=1/(1+zsn)])/(4/0.6774)**3
	print('JJ: M(PopIII) = {:.0f} Msun Mpc^-3, M(PopII) = {:.0f} Msun Mpc^-3'.format(M1, M2))
	
	des0 = np.array(retxt(cr+'deS_popIII_sfr.txt',2,0,0))
	des1 = np.array(retxt(cr+'deS_popIII1_sfr.txt',2,0,0))
	i0 = interp1d(des0[0],np.log(des0[1]))
	i1 = interp1d(des1[0],np.log(des1[1]))
	zb0 = np.linspace(3, 12, 30)
	zb1 = np.linspace(12.1, 29, 50)
	sfrds0 = np.exp(i0(zb0))
	sfrds1 = np.exp(i0(zb1)) + np.exp(i1(zb1))
	zb = np.hstack([zb0, zb1])
	sfrds = np.hstack([sfrds0, sfrds1]) * 0.5

	htwg0 = np.array(retxt(cr+'Htwg_popIII_sfr.txt',2,0,0))
	htwg1 = np.array(retxt(cr+'Htwg_popII_sfr.txt',2,0,0))
	lz0 = np.array([ZT(np.log10(x)) for x in htwg0[0]])
	lz1 = np.array([ZT(np.log10(x)) for x in htwg1[0]])
	#print(sumy(lz0), sumy(lz1))
	hsfr0 = interp1d(lz0, np.log(htwg0[1]))
	hsfr1 = interp1d(lz1, np.log(htwg1[1]))
	zh0 = np.linspace(7.1, 28,  50)
	zh1 = np.linspace(5.7, 7, 20)
	hsfrd0 = np.exp(hsfr0(zh0)) + np.exp(hsfr1(zh0))
	hsfrd1 = np.exp(hsfr1(zh1))
	zh = np.hstack([zh1, zh0])
	hsfrd = np.hstack([hsfrd1, hsfrd0])

	lz = 1/sfr0II[0]-1
	JLW = LWbg(lz, sfr0III[3], sfr0II[3])
	plt.figure()
	plt.plot(lz, JLW)
	plt.xlabel(r'$z$')
	plt.ylabel(r'$J_{\mathrm{LW,bg}}\ [10^{-21}\mathrm{erg\ s^{-1}\ cm^{-2}\ Hz^{-1}\ sr^{-1}}]$')
	plt.yscale('log')
	plt.xlim(xlm)
	plt.tight_layout()
	plt.savefig(rep+'JLW_bg.pdf')
	plt.close()

	if fitlab!=0:
		fit = curve_fit(sfrtf, 1/sfr0III[0]-1, sfr0III[3])#, [2, -4, 12, -6])
		para = fit[0]
		print('Pop III:', para)
		lt = np.linspace(10, TZ(0)/YR/1e6, 1000)/1e3
		lz = np.array([ZT(t) for t in np.log10(lt)])
		lsfrd = sfrtf(lz, *para)
		Mextra = np.trapz(lsfrd, lt)*1e9
		print('Extrapolated total Pop III mass density: {:.0f} Msun Mpc^-3'.format(Mextra))
	else:
		Mextra = 1e5
	
	lMs = []
	lt = np.array([TZ(z) for z in zb])/YR
	Ms = abs(np.trapz(sfrds, lt)*2)
	lMs.append(Ms)
	lMs.append(Ms)
	lt = np.array([TZ(z) for z in lz0])/YR
	lMs.append( abs(np.trapz(htwg0[-1], lt)) )
	
	plt.figure()
	psfr(sfr0III, rlab, mode=1)
	#psfr(sfr1III, 'LB19', '--', V=fac0*(4/0.6774)**3)
	plt.plot(1/ref5[0]-1, ref5[3], 'r--', label='JJ18')
	plt.plot(lz0, htwg0[-1], 'k-.', label='HT16')
	plt.plot(zb, sfrds, 'k:', label='deS11')
	if fitlab!=0:
		plt.plot(zb, sfrtf(zb, *para), lw=3, alpha = 0.5, label=r'Fit')#: $\frac{a(1+z)^{b}}{[1+(1+z)/c]^{d}}$'+'\n'+r'$a={:.2f}$, $b={:.2f}$, $c={:.2f}$, $d={:.2f}$'.format(*para))
	plt.yscale('log')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\mathrm{SFRD}\ [\mathrm{M}_{\odot}\ \mathrm{yr^{-1}\ Mpc^{-3}}]$')
	plt.xlim(xlm)
	plt.ylim(1e-6, 5e-3)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'popIIIsfrd.pdf')
	plt.close()

	lm0 = np.cumsum(sfr0III[1]*sfr0III[3])
	plt.figure()
	plt.plot(1/sfr0III[0]-1, lm0)
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\rho_{\star,\mathrm{PopIII}}\ [\mathrm{M_{\odot}\ Mpc^{-3}}]$')
	plt.yscale('log')
	plt.xlim(xlm)
	plt.tight_layout()
	plt.savefig(rep+'rho_popIII.pdf')
	plt.close()
	
	if fitlab!=0:
		fit = curve_fit(sfrtf, 1/sfr0II[0]-1, sfr0II[3])
		para = fit[0]
		print('Pop II:', para)
	plt.figure()
	psfr(sfr0II, rlab, mode=1)
	#psfr(sfr1II, 'LB19', '--', V=fac0*(4/0.6774)**3)
	plt.plot(1/ref6[0]-1, ref6[3], 'r--', label='JJ18')
	plt.plot(lz1, htwg1[-1], 'k-.', label='HT16')
	if fitlab!=0:
		plt.plot(zb, sfrtf(zb, *para), lw=3, alpha = 0.5, label=r'Fit: $\frac{a(1+z)^{b}}{[1+(1+z)/c]^{d}}$'+'\n'+r'$a={:.2f}$, $b={:.2f}$, $c={:.2f}$, $d={:.2f}$'.format(*para))
	plt.yscale('log')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\mathrm{SFRD}\ [\mathrm{M}_{\odot}\ \mathrm{yr^{-1}\ Mpc^{-3}}]$')
	plt.xlim(xlm)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'popIIsfrd.pdf')
	plt.close()
	
	lm1 = np.cumsum(sfr0II[1]*sfr0II[3])
	plt.figure()
	plt.plot(1/sfr0II[0]-1, lm1)
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\rho_{\star,\mathrm{PopII}}\ [\mathrm{M_{\odot}\ Mpc^{-3}}]$')
	plt.yscale('log')
	plt.xlim(xlm)
	plt.tight_layout()
	plt.savefig(rep+'rho_popII.pdf')
	plt.close()
	
	e1 = interp1d(*ref1)
	e2 = interp1d(*ref2)
	lz = np.linspace(5, 25, 100)
	
	zs = np.linspace(15, xlm[1], 100)
	ts = np.array([TZ(z)/1e9/YR for z in zs])
	t0 = TZ(0)/1e9/YR
	sfrds = sfrmd(ts, t0)
	
	if fitlab!=0:
		fit = curve_fit(sfrtf, 1/sfr0II[0]-1, sfr0II[3]+sfr0III[3])
		para = fit[0]
		print('Total:', para)	
	plt.figure()
	psfr_tot(sfr0II, sfr0III, rlab, mode=1)
	#if fitlab!=0:
	#	plt.plot(zb, sfrtf(zb, *para), lw=3, alpha = 0.5, label=r'Fit')#: $\frac{a(1+z)^{b}}{[1+(1+z)/c]^{d}}$'+'\n'+r'$a={:.2f}$, $b={:.2f}$, $c={:.2f}$, $d={:.2f}$'.format(*para))
	plt.plot(1/ref6[0]-1, ref5[3]+ref6[3], 'r--', label='JJ18')
	plt.plot(zh, hsfrd, 'k-.', label='HT16')
	plt.plot(zs, sfrds, 'k:', label='SL04')
	plt.fill_between(zs, sfrbf(zs)*10**0.2, sfrbf(zs)/10**0.2, facecolor=oi[1], label='MP14', alpha=0.3)
	plt.fill_between(lz, e1(lz), e2(lz), facecolor=oi[0], alpha=0.3, label='MJ19')
	plt.plot(ref0[0], 10**ref0[1], 'o', label='SF16: observed', color='k')
	plt.plot(ref0[0], 10**ref0[2], '^', label='SF16: estimated', color='k')
	#plt.fill_between(ref3[0], 10**ref3[1], 10**ref3[2], facecolor=oi[1], alpha=0.5, label=r'SF18: $\alpha=-2.35$')
	#plt.fill_between(ref4[0], 10**ref4[1], 10**ref4[2], facecolor=oi[2], alpha=0.5, label=r'SF18: $\alpha>-2.35$')
	plt.yscale('log')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\mathrm{SFRD}\ [\mathrm{M}_{\odot}\ \mathrm{yr^{-1}\ Mpc^{-3}}]$')
	plt.xlim(xlm)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'sfrd.pdf')
	plt.close()

	if mode==-2:
		obj = recaesar(sn, rep)
		bhmode = 1
		rfac = 0.5
		rela = gal_info(ds0, obj, rep, bhmode, rfac=rfac)
		plotrela(rela, rep, bhmode, 'gal', sn)
		mlab = 1
		rela = group_info(ds0, obj, rep, bhmode, rfac=rfac, mlab=mlab)
		plotrela(rela, rep, bhmode, 'halo', sn)
		halo_star(rela, rep, sn)
		
		#rela = group_info(ds0, obj, rep, bhmode, rfac=rfac, mlab=0)
		#halo_star(rela, './', sn, Mmin = 0)
	
	demo = 0
	wh = 6.5
	ax = [0, 1]
	#ax = [1, 2]
	#ax = [2,1]
	thk = 400
	cen = np.array([1940, 1960, 1850])
	if demo>0:
		lflag = 0
		ext = 800
		asp = 4/3
		lss = [1, 8, 64]
	else:
		lflag = 1
		ext = 1000
		asp = 1
		lss = [1, 32, 256]
	lext = np.array([ext, ext/asp, thk])
	#le, re = [0, 0, 675-thk/2], [4000, 4000, 675+thk/2]
	le, re = cen - lext/2, cen + lext/2
	#le, re = [1500, 1500, 1850], [2500, 2500, 2150]
	#le, re = [1950, 1850, 1940], [2050, 1950, 1970]
	
	if mode==0:		
		"""
		pm = 1
		met = metal(sn, rep, rlab, Zsun=Zsun, mode=1, nump=1e7)
		print(met)
		ds0 = phase(sn, rep, indm=1, nump=1e7, mode=pm)
		"""
		bh = 1
		#temp_bg(sn, ds0, le, re, ax, rep = rep, mode=bh, lss=lss, lflag=lflag, wh=wh)
		if demo==0:
			LW_bg(sn, ds0, le, re, ax, rep = rep, mode=0)
		
		#denspro(sn, ds0, le, re, [0,1], nb = 100, rep = rep, mode=1)
		#denspro(sn, ds0, le, re, [1,2], nb = 100, rep = rep, mode=1)
		#denspro(sn, ds0, le, re, [2,0], nb = 100, rep = rep, mode=1)
		#denspro(sn, ds0, le, re, ax, nb = 100, rep = rep)
		#main(sn, nbh = 3, R = 3, alp = 0.1, ms = 0.3, rep=rep)
		
	#d = loaddata(sn, sn, rep)
	#ds0 = d['ds']

	#"""
	
	if mode==1:
		llab = ['BK17: FS1', 'BK17: FS2', 'HT16: IMF 1->100']	
		BK1 = np.array(retxt(cr+'BK_GWR1.txt',2,0,0))
		BK2 = np.array(retxt(cr+'BK_GWR2.txt',2,0,0))
		HT0 = np.array(retxt(cr+'HT_GWR.txt',2,0,0))
		lt0 = HT0[0]
		tmax = TZ(0)/1e9/YR
		sel = lt0<=tmax
		lt = lt0[sel]
		lz = [ZT(np.log10(t)) for t in lt]
		ly = HT0[1][sel]
		HT = np.array([lz, ly])
	
		LIGO_ = np.array(retxt(cr+'AdLIGOcurrent.txt',2,0,0))
		LIGO = np.array(retxt(cr+'AdLIGOdesign.txt',2,0,0))
		ET2 = np.array(retxt(cr+'ETxylophone.txt',2,0,0))
		deciHz = np.array(retxt(cr+'deciHz_op.txt',2,0,0))
		dhz = np.array([deciHz[0], deciHz[1]/deciHz[0]**0.5])
		LISA = np.array(retxt(cr+'LISA.txt',2,0,0))
		lisa = np.array([LISA[0], LISA[1]/LISA[0]**0.5])
		lsens = [LIGO, ET2, dhz, lisa]
		lsens0 = [LIGO, ET2, dhz, lisa]
		lins = ['AdLIGO', 'ETxylophone', 'DOoptimal','LISA']
	
		lg = [1.5, 1.0, 0.5]
		l = details(sn, zsn, rep, ncol, npro)# base='bhmergers_'+str(sn))
		#bhhis(l, sn, rep)
		ld = [mergerhis(l, ds0, rg = mode, gamma = g, dz=dz, mode=selab, mrel=rf, zsn=zsn, lsens=lsens,zlim=zlim,snrlim=snrlim, bfrac=bfrac, fage=fage) for g in lg]
		print('bfrac = {}'.format(bfrac))
		for d in ld:
			g = d['gamma']
			Nmg = d['Nmg'] #d['snr'][0]['Nmg']
			Nbbh = d['Nbbh']
			Nhh = d['Nhh']
			Nhl = d['Nhl']
			Nll = d['Nll']
			#le = d['snr'][0]['Dle']
			#lsnr = [dd['snr'] for dd in d['snr']]
			#[print('SNR ({}): '.format(lins[i]), sumy(lsnr[i])) for i in range(len(lins))]
			#print('Eccentricity: ', sumy(le))
			print('gamma = {}: Nmerger: {} ({} Mpc^-3), fraction: {}, efficiency: {} ({}) Msun^-1'.format(g, Nmg, Nmg/V, Nmg/Nbbh, Nll/MpopIII, Nmg/MpopIII))
			print('Nhh = {} ({:.3f}%), Nhl = {} ({:.3f}%), Nll = {} ({:.3f}%)'.format(Nhh, 1e2*Nhh/Nmg, Nhl, Nhl/Nmg*1e2, Nll, Nll/Nmg*1e2))
			i = 0
			for dsnr in d['snr']:
				rate, zmax = dsnr['Drate']
				frac = dsnr['frac']
				print('Detection rate (for z < {}): {:.3f} yr^-1, fraction: {:.3f}, for {}'.format(zmax, rate, frac, lins[i]))
				i += 1
		splothis(ld[0], rep, 128, 0.7, selab=selab, lsens=lsens0, lins=lins, fage=fage, LIGO=LIGO_)
		#for ins in range(len(lins)):
		mplothis(ld, rep, z1=zsn, lref=[BK1, BK2, HT], llab=llab, Mtot=Mextra, 
					lMs=lMs, lins=lins, dind = ins)
	
		if glab>0:
			#"""
			#nr = 20
			lg = np.linspace(0.5, 1.5, nr)
			ld = [mergerhis(l, ds0, rg = mode, gamma = g, dz=dz, mode=selab, mrel=rf, zsn=zsn, lsens=lsens, zlim=zlim,snrlim=snrlim) for g in lg]
			lrate0 = [d['snr'][0]['Drate'][0] for d in ld]
			lrate1 = [d['snr'][1]['Drate'][0] for d in ld]
			lrate2 = [d['snr'][2]['Drate'][0] for d in ld]
			lrate3 = [d['snr'][3]['Drate'][0] for d in ld]
			plt.figure()
			plt.plot(lg, lrate0, label='AdLIGO')
			plt.plot(lg, lrate1, '--', label='ETxylophone')
			plt.plot(lg, lrate2, '-.', label='DOoptimal')
			plt.plot(lg, lrate3, ':', label='LISA')
			plt.legend(loc=2)
			plt.xlabel(r'$\gamma$')
			plt.ylabel(r'$\dot{N}_{\mathrm{detection}}\ [\mathrm{yr^{-1}}]$')#)(z_{\mathrm{GW}}<'+str(zlim)+')\ [\mathrm{yr^{-1}}]$')
			plt.yscale('log')
			plt.tight_layout()
			plt.savefig(rep+'rate_gamma_'+str(selab)+'.pdf')
			plt.close()
			#"""
			
			lrd = np.array([d['GWrate'][1][0]*(d['GWrate'][0][0]<dz) for d in ld])
			plt.figure()
			plt.plot(lg, lrd*1e9)
			plt.xlabel(r'$\gamma$')
			plt.ylabel(r'$\dot{n}_{\mathrm{GW}}(z_{\rm{GW}}\simeq 0)\ [\mathrm{yr^{-1}\ Gpc^{-3}}]$')
			plt.yscale('log')
			plt.tight_layout()
			plt.savefig(rep+'rate_dens_gamma'+str(selab)+'.pdf')
			plt.close()
			
			if selab==1:
				lg = [1.5, 1.0, 0.5]
				#nr = 20
				lbf = np.geomspace(0.2, 1.0, nr)
				md = [[mergerhis(l, ds0, rg = mode, gamma = g, dz=dz, mode=selab, mrel=rf, bfrac=bf, zsn=zsn, lsens=lsens,zlim=zlim,snrlim=snrlim) for bf in lbf] for g in lg]
				for ins in range(len(lins)):
					lrate = np.array([[d['snr'][ins]['Drate'][0] for d in ld] for ld in md])
			
					plt.figure()
					i = 0
					for r, g in zip(lrate, lg):
						plt.plot(lbf, r, ls=lls[i], label=r'$\gamma='+str(g)+'$')
						i += 1
					plt.xlabel(r'$f_{\rm{bulge}}$')
					plt.ylabel(r'$\dot{N}_{\mathrm{'+lins[ins]+r'}}\ [\mathrm{yr^{-1}}]$')
					#plt.ylabel(r'$\dot{N}_{\mathrm{GW}}(z_{\mathrm{GW}}<'+str(zlim)+')\ [\mathrm{yr^{-1}}]$')
					plt.yscale('log')
					#plt.xscale('log')
					plt.legend()
					plt.tight_layout()
					plt.savefig(rep+'rate_bf_'+lins[ins]+'.pdf')
					plt.close()
			
				lrd = np.array([[d['GWrate'][1][0]*(d['GWrate'][0][0]<dz) for d in ld] for ld in md])
				plt.figure()
				i = 0
				for r, g in zip(lrd*1e9, lg):
					plt.plot(lbf, r, ls=lls[i], label=r'$\gamma='+str(g)+'$')
					i += 1
				plt.xlabel(r'$f_{\rm{bulge}}$')
				plt.ylabel(r'$\dot{n}_{\mathrm{GW}}(z_{\rm{GW}}\simeq 0)\ [\mathrm{yr^{-1}\ Gpc^{-3}}]$')
				#plt.xscale('log')
				plt.yscale('log')
				plt.legend()
				plt.tight_layout()
				plt.savefig(rep+'rate_dens_bf.pdf')
				plt.close()
	#"""
	
	if mode==2:
		nbh = countBH(ds0)
		print('z = {:.2f}, N(SBH) = {}, N(DCBH) = {}'.format(*nbh))
		bhind = BHmass_spec(sn, rep, mode = 1)
		mass = BH_mass(ds0)
		#print('z = {:.2f}, M_BH(seed) = {} Msun, M_BH = {} Msun'.format(*mass))
		#ds1 = phase(sn, cr, indm=0)
		#bhind = 0
		ID = mass[3][bhind]
		hsml = mass[4][bhind]
	#if mode==2:
		#plot3d(rep+'dis_'+str(sn)+'.png', ID, hsml, ds0)
	
		bh1 = np.array(retxt_nor(rep+'BH1_fr.txt', 4, 0, 0))
		sel = bh1[0]<=1/(1+zsn)
		Mseed = np.sum(bh1[1][sel] * bh1[2][sel])
		bh2 = np.array(retxt_nor(rep+'BH2_fr.txt', 4, 0, 0))
		sel = bh2[0]<=1/(1+zsn)
		Mseed += np.sum(bh2[1][sel] * bh2[2][sel])	
		print('Total seed mass: {} Msun Mpc^-3'.format(Mseed/V))
	
		bhinfo = np.array(retxt_nor(rep+'blackholes.txt',7,0,0))
		lz = 1/bhinfo[0]-1
		sel = lz>=zsn
		m = bhinfo[2] * UM/Msun/h
		mdot = bhinfo[4]
		print('Total BH mass: {} Msun Mpc^-3'.format(m[sel][-1]/V))
		lt = np.array([TZ(z)/YR for z in lz])
		macc = np.sum(mdot[sel][:-1] * (lt[sel][1:]-lt[sel][:-1])) #np.trapz(mdot[sel], lt[sel])
		macc_ = m[sel][-1] - Mseed
		print('Total accreted mass: {} ({}) Msun Mpc^-3'.format(macc/V, macc_/V))
		print('Accreted ratio: {} ({})'.format(macc/m[sel][-1], macc_/m[sel][-1]))
		#print('NBH_seed/NBH: {}'.format(bhinfo[1][-1]/(nbh[1]+nbh[2])))
	
		plt.figure()
		plt.plot(lz[m>0], m[m>0]/V)
		plt.xlabel(r'$z$')
		plt.ylabel(r'$\rho_{\mathrm{BH}}\ [\mathrm{M_{\odot}\ Mpc^{-3}}]$')
		plt.yscale('log')
		plt.xlim(xlm)
		plt.tight_layout()
		plt.savefig(rep+'rhoBH_z.pdf')
		plt.close()
		plt.figure()
		plt.plot(lz[mdot>0], mdot[mdot>0]/V)
		plt.xlabel(r'$z$')
		plt.ylabel(r'$\dot{\rho}_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ yr^{-1}\ Mpc^{-3}}]$')
		plt.yscale('log')
		plt.xlim(xlm)
		plt.tight_layout()
		plt.savefig(rep+'rhodotacc_z.pdf')
		plt.close()
		#mdot_ = (mdot[:-1] + mdot[1:])/2
		lmacc = np.cumsum(mdot[:-1] * (lt[1:]-lt[:-1]))
		lz = lz[1:]
		plt.figure()
		y1, y2 = 1e-5, 1e5
		plt.plot(lz[lmacc>0], lmacc[lmacc>0]/V)
		plt.fill_between(xlm, [0.66e4, 0.66e4], [1.4e4, 1.4e4], facecolor='gray', label='Salvaterra12: $z=5$')
		plt.plot([5, 5], [y1, y2], 'k:')
		plt.xlabel(r'$z$')
		plt.ylabel(r'$\rho_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ Mpc^{-3}}]$')
		plt.yscale('log')
		plt.xlim(xlm)
		plt.ylim(y1, y2)
		plt.legend(loc=3)
		plt.tight_layout()
		plt.savefig(rep+'rhoacc_z.pdf')
		plt.close()
		rat = lmacc/m[1:]
		plt.figure()
		plt.plot(lz[rat>0], rat[rat>0])
		plt.plot(xlm, [1,1], 'k--')
		plt.xlabel(r'$z$')
		plt.ylabel(r'$\rho_{\mathrm{acc}}/\rho_{\mathrm{BH}}$')
		plt.yscale('log')
		plt.xlim(xlm)
		plt.tight_layout()
		plt.savefig(rep+'ratio_acc_BH_z.pdf')
		plt.close()
	
	"""
		dz = 0.1
		#selab = 1
		ld = [stellarb(sn, rep, gamma=g, dz=dz, sk=sk, frm=frm, mode=selab) for g in [1.5, 1.0, 0.5]]
		#ld = [stellarb(sn, rep, gamma=gamma, dz=dz, sk=k, frm=frm) for k in [1e-3, 1e-2, 1e-1]]
		for d in ld:
			k = d['sk']
			g = d['gamma']
			rate, zlim = d['Drate']
			Nmg = d['Nmg']
			la0 = d['Dla0']
			print('a0: ', sumy(la0))
			print('Nmerger: {} ({} Mpc^-3)'.format(Nmg, Nmg/V))
			print('Detection rate (with gamma = {}, fsk = {}, for z < {}): {} yr^-1'.format(g, k, zlim, rate))
		splothis(ld[0], rep, 32, 0.5, mode=1)
		mplothis(ld, rep, z1=zsn, mode=1, m1 = 70, m2 = 300)
	
		
		nr = 20
		lg = np.linspace(0.5, 1.5, nr)
		lrate = [stellarb(sn, rep, gamma=g, dz=dz, sk=sk, frm=frm)['Drate'][0] for g in lg]
		plt.figure()
		if sk>0:
			plt.plot(lg, lrate, label=r'$f_{\mathrm{sk}}='+str(sk)+'$')
		else:
			plt.plot(lg, lrate)#, label=r'$f_{\mathrm{sk}}='+str(sk)+'$')
		plt.xlabel(r'$\gamma$')
		plt.ylabel(r'$\dot{N}_{\mathrm{GW}}(z_{\mathrm{GW}}<'+str(zlim)+')\ [\mathrm{yr^{-1}}]$')
		plt.legend()
		plt.tight_layout()
		plt.savefig(rep+'rate_gamma_col.pdf')
		plt.close()
		

	lg = [1.5, 1.0, 0.5]
	if mode==3:
		nr = 20
		x1, x2 = 0.001, 0.01
		lsk = np.geomspace(x1, x2, nr)
		out = [lsk]
		for g in lg:
			lrate = [stellarb(sn, rep, gamma=g, dz=dz, sk=k, frm=frm)['Drate'][0] for k in lsk]
			out.append(lrate)
		totxt(rep+'rate_fsk.txt', out, ['fsk']+[str(g) for g in lg], 1, 0)
	if mode==4:
		out = retxt(rep+'rate_fsk.txt', len(lg)+1, 1, 0)
		lsk = out[0]
		x1, x2 = np.min(lsk), np.max(lsk)
		#y1, y2 = 0.2, 1e3
		plt.figure()	
		for i in range(len(lg)):
			lrate = out[i+1]
			plt.plot(lsk, lrate, label=r'$\gamma='+str(lg[i])+'$', ls=lls[i])
			i += 1
		plt.xlabel(r'$f_{\mathrm{sk}}$')
		plt.ylabel(r'$\dot{N}_{\mathrm{GW}}(z_{\mathrm{GW}}<'+str(zlim)+')\ [\mathrm{yr^{-1}}]$')
		plt.legend()
		plt.yscale('log')
		plt.xscale('log')
		plt.xlim(x1, x2)
		#plt.ylim(y1, y2)
		plt.tight_layout()
		plt.savefig(rep+'rate_sk_col.pdf')
		plt.close()
	"""

	#print(fac, V*fac)
	
	
	"""
	sni, snf = 41, 21
	bht = bhtrack(sni, snf, rep, sp = 5)
	plt.figure()
	plt.plot(bht[0], bht[1], label=r'$M_{\mathrm{BH}}$', marker='.')
	#plt.plot(bht[0], bht[2], '--', label=r'$M_{\mathrm{acc}}$')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$M\ [\mathrm{M_{\odot}}]$')
	plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'Maxhis.pdf')
	plt.close()
	"""
	
	"""
	sn1, sn2 = 73, 74
	R = 10
	ds1 = yt.load(rep+'snapshot_'+str(sn1).zfill(3)+'.hdf5')
	ds2 = yt.load(rep+'snapshot_'+str(sn2).zfill(3)+'.hdf5')
	pos = np.array(countBH(ds2)[4][bhind])
	plot3d(rep+'pre.pdf', pos, R, ds1)
	plot3d(rep+'post.pdf', pos, R, ds2)
	"""
	
	"""
	lnbh = []
	for n in range(40, 50):
		ds0 = yt.load(rep+'snapshot_'+str(n).zfill(3)+'.hdf5')
		nbh = countBH(ds0)
		lnbh.append(nbh)
	for x in lnbh:
		print('z = {:.2f}, N(SBH) = {}, N(DCBH) = {}'.format(*x))#N(DCBH,Z<Zth) = {}'.format(*x))
	"""
	
	"""
	ID = np.array(countBH(ds0)[4][bhind])
	for sn in range(72, 75):
		R = 3
		ds = yt.load(rep+'snapshot_'+str(sn).zfill(3)+'.hdf5')
		plot3d(rep+'dis_'+str(sn)+'.pdf', ID, R, ds)
		main(sn, nbh = 3, R = 3, alp = 0.1, ms = 0.3, rep=rep)
	"""
	
	"""
	bhind = 2
	mode = 0
	fac1 = 10.
	fac2 = 5.
	sns = 49
	snf = 47
	ds0 = yt.load(rep+'snapshot_'+str(sns).zfill(3)+'.hdf5')
	pos = np.array(countBH(ds0)[3][bhind])
	obj = recaesar(sns, rep)
	halo, rat = findhost(pos, obj, fac1, mode=mode)
	his = track(halo, sns, snf, fac2, rep=rep)
	prop = haloget(his)
	print('DCBH position: ',  pos)
	print('Host halo center: ', prop['pos'][0])
	print('r/r200: ', rat)
	print('Match fraction: ', his[2])
	print('Stellar masses: ', prop['mstr'])
	totxt(rep+'HostInfo.txt', [prop[k] for k in prop if k is not 'pos'], [k for k in prop if k is not 'pos'], 1, 0)
	
	plt.figure()
	plt.plot(his[1], prop['mtot'])
	plt.plot(his[1], Mup(his[1]), 'k--', label='Atomic')
	plt.legend()
	plt.yscale('log')
	plt.xlim(np.max(his[1]), np.min(his[1]))
	plt.xlabel(r'Redshift')
	plt.ylabel(r'$M_{h}\ [\mathrm{M_{\odot}}]$')
	plt.tight_layout()
	#plt.show()
	plt.savefig(rep+'AccHis_'+str(bhind)+'.pdf')
	plt.close()
	"""
	
	#"""
	if mode==3:
		lze = np.array([8, 7, 6, 5.5, 5, 4.5, 4])
		z0 = 24
		ze = 4
		dt = 10e6
		ntot = 24
		ns = 10
		#nr = int((ntot-ns)/2)
		nr = ntot - ns - len(lze)
		z1 = 9
		a1 = 1/(1+z1)
		t1 = TZ(z1)/YR
		t2 = t1 + dt
		lt = np.linspace(t1, t2, ns+1)
		lam = np.array([1/(ZT(np.log10(t/1e9))+1) for t in lt])
		la0 = np.linspace(1/(1+z0), lam[0], nr+1)
		#lae = np.linspace(lam[-1], 1/(1+ze), nr+1)
		lae = 1/(1+lze)
		la = np.hstack([la0, lam[1:-1], lae])
		totxt('output_{}_{}_{}.txt'.format(z0, z1, ze), [la], 0, 0, 0)
	#"""
