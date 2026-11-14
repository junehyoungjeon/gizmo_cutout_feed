from bh_bl import *
import os
fac = 0.3783187*0.3497985*0.399293243
bs = 4
V = (bs/0.6774)**3 * fac

if __name__=="__main__":
	fcp = int(sys.argv[1])
	lab1 = 'FDzoom'
	if fcp==0:
		sn1, sn2 = 78, 78
		rep = 'NSFDBKcomp/'
		rep1 = 'FDzoom_Hseed_comp/'
		rep2 = 'NSFDBKzoom_Hseed/'
		lab1 = 'FDzoom'
		lab2 = 'NSFDBKzoom'
		clab = 'NSFDBK'
	elif fcp==1:
		sn1, sn2 = 80, 80
		#rep = 'HRcomp_new/'
		rep = 'HRcomp/'
		rep1 = 'FDzoom_Lseed_comp/'
		#rep1 = 'FDzoom_Lseed_ADC/'
		rep2 = 'FDzoom_Lseed_HR/'
		lab2 = 'FDzoomHR'
		clab = 'HR'
	else:
		sn1, sn2 = 100, 23
		rep = 'ADCcomp/'
		rep1 = 'FDzoom_Lseed_comp/'
		rep2 = 'FDzoom_Lseed/'#_ADC/'
		lab2 = 'FDzoomADC'
		clab = 'ADC'
		
	cr = 'comp/'
	snap = 1
	acc = 1
	dind = 1

	ncol = 15
	dz = 0.2
	gamma = 1.0
	bfrac = 0.6
	fage = 0
	zlim = 5
	snrlim = 10
	alp = 0.7
	
	if not os.path.exists(rep):
		os.makedirs(rep)
	
	ds1 = yt.load(rep1+'snapshot_'+str(sn1).zfill(3)+'.hdf5')
	ds2 = yt.load(rep2+'snapshot_'+str(sn2).zfill(3)+'.hdf5')
	zsn = ds1['Redshift']
	h = ds1['HubbleParam']
	print('zsn = {}'.format(zsn))
	
	l1 = np.array(retxt(rep1+'bhmergers_'+str(sn1)+'.txt',ncol,0.0))
	l2 = np.array(retxt(rep2+'bhmergers_'+str(sn2)+'.txt',ncol,0.0))
	
	ref0 = np.array(retxt(cr+'SF16.txt', 3, 0, 0))
	ref1 = np.array(retxt(cr+'edges_bottom.txt', 2, 0, 0))
	ref2 = np.array(retxt(cr+'edges_top.txt', 2, 0, 0))
	
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
	
	sfr13 = np.array(retxt_nor(rep1+'popIII_sfr.txt', 4, 0, 0))
	sfr23 = np.array(retxt_nor(rep2+'popIII_sfr.txt', 4, 0, 0))
	sfr12 = np.array(retxt_nor(rep1+'popII_sfr.txt', 4, 0, 0))
	sfr22 = np.array(retxt_nor(rep2+'popII_sfr.txt', 4, 0, 0))
	lm13 = sfr13[1]*sfr13[2]
	lm23 = sfr23[1]*sfr23[2]
	lm12 = sfr12[1]*sfr12[2]
	lm22 = sfr22[1]*sfr22[2]
	M13 = np.sum(lm13[sfr13[0]<=1/(1+zsn)])
	M23 = np.sum(lm23[sfr23[0]<=1/(1+zsn)])
	M12 = np.sum(lm12[sfr12[0]<=1/(1+zsn)])
	M22 = np.sum(lm22[sfr22[0]<=1/(1+zsn)])
	print('OTF ({}/{}): M(PopIII) = {:.0f}/{:.0f}={:.4f} [Msun Mpc^-3], M(PopII) = {:.0f}/{:.0f}={:.4f} [Msun Mpc^-3]'.format(lab1, lab2, M13/V, M23/V, M13/M23, M12/V, M22/V, M12/M22))
	
	if snap>0:
		M13, M12 = smass(ds1)
		M23, M22 = smass(ds2)
		print('snapshot ({}/{}): M(PopIII) = {:.0f}/{:.0f}={:.4f} [Msun Mpc^-3], M(PopII) = {:.0f}/{:.0f}={:.4f} [Msun Mpc^-3]'.format(lab1, lab2, M13/V, M23/V, M13/M23, M12/V, M22/V, M12/M22))
	
	xlm = [25, 4]
	
	plt.figure()
	psfr(sfr13, lab1, mode=1)
	psfr(sfr23, lab2, '--', mode=1)
	plt.plot(lz0, htwg0[-1], 'k-.', label='HT16')
	plt.plot(zb, sfrds, 'k:', label='deS11')
	plt.yscale('log')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\mathrm{SFRD}\ [\mathrm{M}_{\odot}\ \mathrm{yr^{-1}\ Mpc^{-3}}]$')
	plt.xlim(xlm)
	plt.ylim(1e-6, 5e-3)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'popIIIsfrd_'+clab+'.pdf')
	plt.close()
	
	e1 = interp1d(*ref1)
	e2 = interp1d(*ref2)
	lz = np.linspace(5, 25, 100)
	
	zs = np.linspace(15, xlm[1], 100)
	ts = np.array([TZ(z)/1e9/YR for z in zs])
	t0 = TZ(0)/1e9/YR
	sfrds = sfrmd(ts, t0)
	
	plt.figure()
	psfr_tot(sfr13, sfr12, lab1, mode=1)
	psfr_tot(sfr23, sfr22, lab2, '--', mode=1)
	plt.plot(zh, hsfrd, 'k-.', label='HT16')
	plt.plot(zs, sfrds, 'k:', label='SL04')
	plt.fill_between(zs, sfrbf(zs)*10**0.2, sfrbf(zs)/10**0.2, facecolor=oi[1], label='MP14', alpha=0.3)
	plt.fill_between(lz, e1(lz), e2(lz), facecolor=oi[0], alpha=0.3, label='MJ19')
	plt.plot(ref0[0], 10**ref0[1], 'o', label='SF16: observed', color='k')
	plt.plot(ref0[0], 10**ref0[2], '^', label='SF16: estimated', color='k')
	plt.yscale('log')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\mathrm{SFRD}\ [\mathrm{M}_{\odot}\ \mathrm{yr^{-1}\ Mpc^{-3}}]$')
	plt.xlim(xlm)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'sfrd_'+clab+'.pdf')
	plt.close()
	
	if acc>0:
		bhinfo1 = np.array(retxt_nor(rep1+'blackholes.txt',7,0,0))
		bhinfo2 = np.array(retxt_nor(rep2+'blackholes.txt',7,0,0))
		lz1 = 1/bhinfo1[0]-1
		lt1 = np.array([TZ(z)/YR for z in lz1])
		mdot1 = bhinfo1[4]
		lz2 = 1/bhinfo2[0]-1
		lt2 = np.array([TZ(z)/YR for z in lz2])
		mdot2 = bhinfo2[4]
	
		lmacc1 = np.cumsum(mdot1[:-1] * (lt1[1:]-lt1[:-1]))
		lz1 = lz1[1:]
		lmacc2 = np.cumsum(mdot2[:-1] * (lt2[1:]-lt2[:-1]))
		lz2 = lz2[1:]
		plt.figure()
		y1, y2 = 1e-5, 1e5
		plt.plot(lz1[lmacc1>0], lmacc1[lmacc1>0]/V, label=lab1)
		plt.plot(lz2[lmacc2>0], lmacc2[lmacc2>0]/V, '--', label=lab2)
		plt.fill_between(xlm, [0.66e4, 0.66e4], [1.4e4, 1.4e4], facecolor='gray', label='Salvaterra12: $z=5$')
		plt.plot([5, 5], [y1, y2], 'k:')
		plt.xlabel(r'$z$')
		plt.ylabel(r'$\rho_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ Mpc^{-3}}]$')
		plt.yscale('log')
		plt.xlim(xlm)
		plt.ylim(y1, y2)
		plt.legend(loc=3)
		plt.tight_layout()
		plt.savefig(rep+'rhoacc_z_'+clab+'.pdf')
		plt.close()
	
	Mh_Ms_Hseed = lambda x: 10**((np.log10(x)-0.93)/0.34)
	Mh_Ms_Hseed_ = lambda x: 10**((np.log10(x)-0.47)/0.42)
	Mh_Ms_Lseed = lambda x: 10**((np.log10(x)-1)/0.22)

	if fcp==0:
		rf1 = Mh_Ms_Hseed
		rf2 = Mh_Ms_Hseed_
	else:
		rf1 = Mh_Ms_Lseed
		rf2 = Mh_Ms_Lseed
	
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
	lref = [BK1, BK2, HT]
	
	LIGO = np.array(retxt(cr+'AdLIGOdesign.txt',2,0,0))
	ET2 = np.array(retxt(cr+'ETxylophone.txt',2,0,0))
	deciHz = np.array(retxt(cr+'deciHz_op.txt',2,0,0))
	dhz = np.array([deciHz[0], deciHz[1]/deciHz[0]**0.5])
	LISA = np.array(retxt(cr+'LISA.txt',2,0,0))
	lisa = np.array([LISA[0], LISA[1]/LISA[0]**0.5])
	lsens = [LIGO, ET2, dhz, lisa]
	lsens0 = [LIGO, ET2, dhz, lisa]
	lins = ['AdLIGO', 'ETxylophone', 'DOoptimal','LISA']
	
	rf = ms_mbh_default
	d10 = mergerhis(l1, ds1, rg =1, gamma=gamma, dz=dz, mode=0, mrel=rf, zsn=zsn, lsens=lsens, zlim=zlim, snrlim=snrlim, bfrac=bfrac, fage=fage, V=V)
	d20 = mergerhis(l2, ds2, rg =1, gamma=gamma, dz=dz, mode=0, mrel=rf, zsn=zsn, lsens=lsens, zlim=zlim, snrlim=snrlim, bfrac=bfrac, fage=fage, V=V)

	d11 = mergerhis(l1, ds1, rg =1, gamma=gamma, dz=dz, mode=1, mrel=rf1, zsn=zsn, lsens=lsens, zlim=zlim, snrlim=snrlim, bfrac=bfrac, fage=fage, V=V)
	d21 = mergerhis(l2, ds2, rg =1, gamma=gamma, dz=dz, mode=1, mrel=rf2, zsn=zsn, lsens=lsens, zlim=zlim, snrlim=snrlim, bfrac=bfrac, fage=fage, V=V)
	
	lm3 = [M13, M23, M13, M23]
	ld = [d10, d20, d11, d21]
	
	gwlab = [lab1+'_obs-based', lab2+'_obs-based', lab1+'_sim-based', lab2+'_sim-based']
	lfs = ['full', 'full', 'none', 'none']
	plt.figure()
	i = 0
	y1 = 100
	y2 = 0
	x2 = 0
	for d, m3 in zip(ld, lm3):
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
		print('gamma = {}: Nmerger: {} ({} Mpc^-3), fraction: {}, efficiency: {} ({}) Msun^-1'.format(g, Nmg, Nmg/V, Nmg/Nbbh, Nll/m3, Nmg/m3))
		print('Nhh = {} ({:.3f}%), Nhl = {} ({:.3f}%), Nll = {} ({:.3f}%)'.format(Nhh, 1e2*Nhh/Nmg, Nhl, Nhl/Nmg*1e2, Nll, Nll/Nmg*1e2))
		j = 0
		for dsnr in d['snr']:
			rate, zmax = dsnr['Drate']
			frac = dsnr['frac']
			print('Detection rate (for z < {}): {:.3f} yr^-1, fraction: {:.3f}, for {}'.format(zmax, rate, frac, lins[j]))
			j += 1
			
		lx, ly, yer = d['GWrate']
		if y1>np.min(ly[ly>0]):
			y1 = np.min(ly[ly>0])
		if y2<np.max(ly):
			y2 = np.max(ly)
		if x2<np.max(lx):
			x2 = np.max(lx)
		plt.errorbar(lx[ly>0], ly[ly>0]*1e9, yerr=yer[ly>0]*1e9, xerr=dz/2, fmt=lfmt[i], label=gwlab[i], alpha =alp, fillstyle=lfs[i])
		i += 1
	y1, y2 = min(0.01, y1*1e9*0.5), max(y2*1e9 * 2, 1e3)
	x1, x2 = 0, min(x2*1.05, 5)
	if len(lref)>0:
		i = 0
		for ref in lref:
			plt.plot(ref[0], ref[1], color='k', ls=lls[i], label=llab[i])
			i += 1
	plt.xlim(x1, x2)
	plt.ylim(y1, y2)
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.ylabel(r'$\dot{n}_{\mathrm{GW}}\ [\mathrm{yr^{-1}\ Gpc^{-3}}]$')
	plt.yscale('log')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig(rep+'GW_rate_dens_'+clab+'.pdf')
	plt.close()
	
	lslab = [lab1, lab2]
	
	i = 0
	y2 = 0
	plt.figure()
	for d in ld:
		lx = d['snr'][dind]['Dzbase']
		ly = d['snr'][dind]['Dlr']
		if y2<np.max(ly):
			y2 = np.max(ly)
		plt.plot(lx[ly>0], ly[ly>0], ls = lls[i], label=gwlab[i])
		i += 1
	y1, y2 = 0.1, y2 * 2
	x1, x2 = 0, zlim
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
	plt.savefig(rep+'GW_rate_'+inst+'_'+clab+'.pdf')
	plt.close()
	
	i = 0
	plt.figure()
	for d in [d10, d20]:
		lx, ly, yer = d['Mrate']
		plt.errorbar(lx[ly>0], ly[ly>0]*1e9, yerr=yer[ly>0]*1e9, xerr=dz/2, fmt=lfmt[i], label=lslab[i], alpha = alp, fillstyle=lfs[i])
		i += 1
	plt.legend()
	plt.xlabel(r'$z_{\mathrm{BBH}}$')
	plt.ylabel(r'$\dot{n}_{\mathrm{BBH}}\ [\mathrm{Gpc^{-3}\ yr^{-1}}]$')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'merger_rate_'+clab+'.pdf')
	plt.close()
	
