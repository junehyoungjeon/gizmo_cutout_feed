from bh_bl import *

def findstar(ds, stype = 2):
	ad = ds.all_data()
	st = ad[('PartType4', 'Star Type')]
	if stype==2:
		sel = st==2
	if stype==3:
		sel = st>=30000
	age = np.array(ad[('PartType4', 'Stellar Age')])[sel]
	lid = ad[('PartType4', 'ParticleIDs')][sel]
	mdt = 1e7
	for t, sid in zip(age, lid):
		if abs(t)<mdt:
			mdt = t
			ID = sid
	return [ID, mdt, lid[age==0]]
	
oi = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#000000']
mst = ['o', '^', 's', 'D', 'h', 'p']
lls = ['-', '--', '-.', ':'] * 2
	
def ngbstar(ad, cen, t0, R = 0.2, Bxs = 4000):
	count = 0
	Rmin = 2
	pos = np.array(ad[('PartType4','Coordinates')].to('kpccm/h'))
	age = np.array(ad[('PartType4', 'Stellar Age')])
	st = ad[('PartType4', 'Star Type')]
	lty = []
	for x, t, ty in zip(pos, age, st):
		dis = pbd(x-cen, Bxs)
		r = np.linalg.norm(dis)
		if (r<R) and t <= t0:
			count += 1
			lty.append(ty)
			if (r<Rmin) and (r>0):
				Rmin = r
	if len(lty)>0:
		print(sumy(lty))
	return [count, Rmin]
	
TSF = 1e3
def main(ID, lsn = [0], rep = './', R = 2, x00 = 1e-2, nbin = 20, Bxs = 4000, alp = 0.1, ms = 0.5, cs = oi, mst = mst, lls = lls, dm = 1, lab = 'popII', Rs = 0.2, sind = 0):
	norm = Hfrac/PROTON
	plt.figure(figsize=(8,6))
	ax0 = plt.subplot(221)
	ax1 = plt.subplot(222)
	ax2 = plt.subplot(223)
	ax3 = plt.subplot(224)
	i = 0
	for sn in lsn:
		ds = yt.load(rep+'snapshot_'+str(sn).zfill(3)+'.hdf5')
		z = ds['Redshift']
		if i==0:
			z0 = z
			h = ds['HubbleParam']
			x0 = x00/(1+z)/h
		ad0 = ds.all_data()
		lid = ad0[('PartType4', 'ParticleIDs')]
		sel = lid==ID
		if np.sum(sel)>0:
			age0 = np.array(ad0[('PartType4', 'Stellar Age')])[sel][0]
			cen = np.array(ad0[('PartType4','Coordinates')].to('kpccm/h'))[sel][0]
			h2sf = np.array(ad0[('PartType4','H2_SF')])[sel][0]
			Zsf = np.array(ad0['PartType4', 'Metallicity_00'])[sel][0]/Zsun
		else:
			lid = ad0[('PartType3', 'ParticleIDs')]
			sel = lid==ID
			age0 = (TZ(z) - TZ(z0))/YR
			cen = np.array(ad0[('PartType3','Coordinates')].to('kpccm/h'))[sel][0]
			h2sf = np.array(ad0[('PartType3','H2_SF')])[sel][0]
			Zsf = np.array(ad0['PartType3', 'Metallicity_00'])[sel][0]/Zsun
		age = np.round(age0/1e5)/1e1
		if age<0.1:
			age = 0.1
		ad = ds.sphere(cen, (R, 'kpccm/h'))
		print('Nstar = {}, rmin = {} kpc (for r < {} kpc)'.format(*ngbstar(ad, cen, age0, Rs, Bxs), Rs))
		pos = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))
		vel = np.array(ad[('PartType0','Velocities')].to('km/s'))
		T = temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')])
		rho = np.array(ad[('PartType0','Density')].to('g/cm**3'))
		Z = np.array(ad['PartType0', 'Metallicity_00'])/Zsun
		xh2 = np.array(ad['PartType0', 'Primordial H2'])
		prho = profile1(pos, cen, rho, R, nbin, fac, 1, Bxs)
		ptemp = profile1(pos, cen, T, R, nbin, fac, 1, Bxs)
		#pZ = profile1(pos, cen, Z, R, nbin, fac, 1, Bxs)
		pv = profile2(pos, cen, vel, R, nbin, fac, 1, Bxs)
		ph2 = profile1(pos, cen, xh2, R, nbin, fac, 1, Bxs)
		p0, p1 = prho['pro'], ptemp['pro']
		p2, p3 = pv['pro'], ph2['pro']
		#if dm!=0:
		d0, d1 = prho['data'], ptemp['data']
		d2, d3 = pv['data'], ph2['data']
		#ax0.errorbar(p0[0], p0[1]*norm, p0[2]*norm, ls=lls[i], color=cs[i], zorder=2)
		#ax1.errorbar(p1[0], p1[1], p1[2], ls=lls[i], color = cs[i], zorder=2)
		ax0.plot(p0[0]/(1+z)/h, p0[1]*norm, ls=lls[i], color=cs[i], zorder=2)
		ax1.plot(p1[0]/(1+z)/h, p1[1], ls=lls[i], color = cs[i], zorder=2)
		ax2.plot(p2[0]/(1+z)/h, p2[1], ls=lls[i], zorder = 2, color = cs[i])
		ax3.plot(p3[0]/(1+z)/h, p3[1], ls=lls[i], zorder = 2, color = cs[i], label=r'$t='+str(age)+'$ Myr')
		#if dm!=0:
		ax0.plot(d0[0]/(1+z)/h, d0[1]*norm, '.', alpha = alp, zorder = 1, color = cs[i])
		ax1.plot(d1[0]/(1+z)/h, d1[1], '.', alpha = alp, zorder = 1, color = cs[i])
		ax2.plot(d2[0]/(1+z)/h, d2[1], '.', alpha = alp, zorder = 1, color = cs[i])
		ax3.plot(d3[0]/(1+z)/h, d3[1], '.', alpha = alp, zorder = 1, color = cs[i])
		if i==0:
			print('H2_SF: {}, Z_SF: {}'.format(h2sf, Zsf))
		#	ax3.plot([x0*2], [h2sf], mst[i], color = cs[i], alpha = 0.5)
		i += 1
			
	ax0.plot([x0, R/(1+z)/h], [nth_sf, nth_sf], 'k-.', label=r'$n_{\mathrm{th}}$')
	#ax0.plot([x0, R/(1+z)/h], [nth_hii, nth_hii], 'k:', label=r'$n_{\mathrm{th}}$(HII)')
	ax1.plot([x0, R/(1+z)/h], [TSF, TSF], 'k--', label=r'$T_{\mathrm{th}}$')
	ax2.plot([x0, R/(1+z)/h], [0, 0], 'k:')
	#ax1.fill_between([x0, R], [T1, T1], [T2, T2], facecolor='k', zorder = 0, label=r'$T$(BH)', alpha = alp*3)
	
	ax0.set_ylabel(r'$n_{\mathrm{H}}\ [\mathrm{cm^{-3}}]$')
	ax1.set_ylabel(r'$T\ [\mathrm{K}]$')
	#ax2.set_ylabel(r'$Z\ [\mathrm{Z_{\odot}}]$')
	ax2.set_ylabel(r'$-v_{r}\ [\mathrm{km\ s^{-1}}]$')
	ax3.set_ylabel(r'$\mathrm{[H_{2}/H]}$')
	ax0.set_yscale('log')
	ax1.set_yscale('log')
	#ax2.set_yscale('log')
	ax3.set_yscale('log')
	ax0.set_xscale('log')
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')
	ax0.set_xlim(x0, R/(1+z)/h)
	ax1.set_xlim(x0, R/(1+z)/h)
	ax2.set_xlim(x0, R/(1+z)/h)
	ax3.set_xlim(x0, R/(1+z)/h)
	ax0.set_xlabel(r'$r\ [\mathrm{kpc}]$')
	ax1.set_xlabel(r'$r\ [\mathrm{kpc}]$')
	ax2.set_xlabel(r'$r\ [\mathrm{kpc}]$')
	ax3.set_xlabel(r'$r\ [\mathrm{kpc}]$')
	ax0.legend(loc=3)
	ax1.legend()
	ax3.legend()
	plt.tight_layout()
	plt.savefig(rep+'Spro_'+lab+'_'+str(sind)+'.pdf')
	#plt.show()
	plt.close()	
	
	plot3d(rep+'dis_'+lab+'_'+str(sind)+'.pdf', ID, Rs, ds, age = age0)
	
	return 0
	
if __name__ == "__main__":
	rep = 'SFtest/'
	sn = 1
	lsn1 = np.array([0, 3, 12, 24]) + sn
	lsn2 = np.array([0, 3, 6, 9]) + sn
	ds0 = yt.load(rep+'snapshot_'+str(sn).zfill(3)+'.hdf5')
	print(ds0.field_list)
	f1 = findstar(ds0, 3)
	print(f1)
	ID0 = f1[0]
	sind = sn
	R = 10
	main(ID0, lsn1, rep=rep, lab='popIII', alp = 0.2, sind = sind, x00 = 5e-2, R = R)
	"""
	f2 = findstar(ds0, 2)
	print(f2)
	Ns = len(f2[2])
	print(Ns)
	sind = 0
	ID1 = f2[2][sind]
	main(ID1, lsn2, rep=rep, lab='popII', dm = 0, sind = sind)
	"""
