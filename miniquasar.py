from cosmology import *
import os
import multiprocessing as mp
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from txt import *
#plt.style.use('test2')

def Macc_edd(MBH, eps = 0.1, mu = 1.22):
	return 4*np.pi*GRA*PROTON*mu/SPEEDOFLIGHT/SIGMATH/eps * YR * MBH
	#return 2.2e-6 * (MBH/100)

def LBH(Mdot, eps = 0.1):
	return Mdot * Msun/YR * SPEEDOFLIGHT**2 * eps
	
def Tinner(MBH, Mdot, fin = 3):
	rs = 2*GRA*MBH*Msun/SPEEDOFLIGHT**2
	rin = fin*rs
	Tin = (3*GRA*MBH*Mdot * Msun**2/YR/ (8*np.pi*STEFAN*rin**3))**(1./4)
	return Tin
	
e1, e2, e3 = 11.2, 13.6, 10e3
nu1, nu2, nu3 = e1*eV/PLANCK, e2*eV/PLANCK, e3*eV/PLANCK
nu_NT = 0.2e3*eV/PLANCK
print('Transition temperature: {:.0f} K'.format(nu_NT*PLANCK/BOL))
print('Eddington rate: {} (MBH/(100 Msun)) [Msun yr^-1]'.format(Macc_edd(100)))

def MCD_spec(nu, Tin, rin, rout, p, xmax = 200, xmin = 1e-10):
	xin = PLANCK*nu/(BOL*Tin)
	xout = xin * (rout/rin)**p
	out = 2*PLANCK/SPEEDOFLIGHT**2 * (BOL*Tin/(PLANCK*nu))**(2./p) * nu**3
	def integrand(x):
		return (x**(2./p)-1)/(np.e**x-1)
	x1 = max(xin, xmin)
	x2 = min(xout, xmax)
	flag = 0
	if x1>xmax:
		res = integrate.quad(integrand, xmax*(xin/xout), xmax)[0] * np.exp(xmax-x1)
		flag = 1
	if x2<xmin:
		res = integrate.quad(integrand, xmin, xmin*(xout/xin))[0] * (x2/xmin)**(2./p)	
		flag = 1
	if flag==0:
		res = integrate.quad(integrand, x1, x2)[0]
	#if xout<=lim:
	#	res = integrate.quad(integrand, xin, xout)[0]
	#else:
	#	res = integrate.quad(integrand, xin, lim)[0] #* np.exp(-(xout-lim))
	return out * res * 16*np.pi**2 * rin**2

def BH_spec(MBH, Mdot, fac = 0.5, beta = 1, 
			fin = 3, fout = 1e4, p = 3./4, eps = 0.1, Nnu = 10000, minfac = 1e-4, Tmin = 1e4):
	Ltot = LBH(Mdot, eps)
	rs = 2*GRA*MBH*Msun/SPEEDOFLIGHT**2
	rin = fin*rs
	rout = fout*rs
	Tin = max(Tmin, (3*GRA*MBH*Mdot * Msun**2/YR/ (8*np.pi*STEFAN*rin**3))**(1./4))
	nu_min = min(nu1, minfac * Tin*BOL/PLANCK)
	if nu_min>=nu_NT:
		nu_min = nu_NT/1e3
	lnu0 = np.linspace(nu_min, nu_NT, Nnu)
	spec0 = np.array([MCD_spec(nu, Tin, rin, rout, p) for nu in lnu0])
	norm0 = Ltot * fac / np.trapz(spec0, lnu0)
	if beta==1:
		norm1 = np.log(nu3/nu_NT)
	else:
		norm1 = 1./(1.-beta)*( nu3**(1-beta) - nu_NT**(1-beta) )
	norm1 = Ltot * (1.-fac) * norm1
	lnu1 = np.geomspace(nu_NT, nu3, Nnu)
	spec1 = lnu1**(-beta)
	
	lnu_LW = np.linspace(nu1, nu2, Nnu)
	spec_LW = np.array([MCD_spec(nu, Tin, rin, rout, p) for nu in lnu_LW])
	L_LW = np.trapz(spec_LW, lnu_LW) * norm0
	N_LW = np.trapz(spec_LW/(lnu_LW*PLANCK), lnu_LW) * norm0
	
	lnu_ion = np.linspace(nu2, nu_NT, Nnu)
	spec_ion = np.array([MCD_spec(nu, Tin, rin, rout, p) for nu in lnu_ion])
	L_ion = np.trapz(spec_ion, lnu_ion) * norm0 + Ltot * (1-fac)
	N_ion = np.trapz(spec_ion/(lnu_ion*PLANCK), lnu_ion) * norm0
	N_ion += Ltot * (1-fac)/PLANCK /beta * (nu_NT**-beta - nu3**-beta)
	
	d = {}
	specMCD = spec0*norm0
	specNT = spec1*norm1
	lnu = np.hstack([lnu0[1:], lnu1])
	lspec = np.hstack([specMCD[1:], specNT])
	d['MCD_orig'] = [lnu0[1:]*PLANCK/eV, spec0[1:]]
	d['MCD'] = [lnu0[1:]*PLANCK/eV, specMCD[1:]]
	d['NT'] = [lnu1*PLANCK/eV, specNT]
	d['tot'] = [lnu*PLANCK/eV, lspec]
	d['LW'] = [L_LW, N_LW]
	d['ion'] = [L_ion, N_ion]
	return d
	
def main(m1 = 10, m2 = 1e9, a1 = 1e-10, a2 = 1e4, nbin = 32, ncore = 8, 
		Nnu = 1000, fac = 0.5, up = 100, Tmin = 1e4):
	lm = np.geomspace(m1, m2, nbin)
	la = np.geomspace(a1, a2, nbin)
	#print(lm)
	#print(la)
	X, Y = np.meshgrid(lm, la, indexing='ij')
	Nion = np.zeros(X.shape)
	Lion = np.zeros(X.shape)
	LLW = np.zeros(X.shape)
	NLW = np.zeros(X.shape)
	np_core = int(nbin/ncore)
	lpr = [[i*np_core, (i+1)*np_core] for i in range(ncore-1)] + [[(ncore-1)*np_core, nbin]]
	manager = mp.Manager()
	def sess(pr0, pr1, j):
		out = []
		for i in range(pr0, pr1):
			if la[j] > Macc_edd(lm[i]) * up:
				out.append([0, 0, 0, 0])
			else:
				d = BH_spec(MBH=lm[i], Mdot=la[j], fac = fac, Nnu = Nnu, Tmin = Tmin)
				out.append(d['ion']+d['LW'])
		output.put((pr0, np.array(out).T))
	for i in range(nbin):
		output = manager.Queue()
		pros = [mp.Process(target=sess, args=(lpr[k][0], lpr[k][1], i)) for k in range(ncore)]
		for p in pros:
			p.start()
		for p in pros:
			p.join()
		out = [output.get() for p in pros]
		out.sort()
		Lion[:,i] = np.hstack([x[1][0] for x in out])
		Nion[:,i] = np.hstack([x[1][1] for x in out])
		LLW[:,i] = np.hstack([x[1][2] for x in out])
		NLW[:,i] = np.hstack([x[1][3] for x in out])
	return X, Y, Lion, Nion, LLW, NLW
	
def linear(x, a, b):
	return a*x + b
	
c = 8.5e47
print('Nion = {} s^-1 [Mdot/(Msun yr^-1)]'.format(c/Macc_edd(1)))
def model_Nion(mbh, eta, a=1, b=1, c=c):
	return mbh**a * eta**b * c

def time_isobhb_circ(a0, M1, M2, e0 = 0):
	m = (M1 + M2)*Msun
	m1m2 = M1 * M2 * Msun**2
	e = e0
	Fe = (1-e**2)**(-3.5)*(1+(73./24)*e**2+(37./96)*e**4)
	B = 64*GRA**3*m1m2*m*Fe/(5*SPEEDOFLIGHT**5)
	return a0**4/(4*B)/YR

def eint(e):
	return e**(29/19)*(1+(121/304)*e**2)**(1181/2299)/(1-e**2)**(3/2)

def time_isobhb(a0, M1, M2, e0 = 0):
	m = (M1 + M2)*Msun
	m1m2 = M1 * M2 * Msun**2
	beta = 64*GRA**3*m1m2*m/(5*SPEEDOFLIGHT**5)
	if e0==0:
		return a0**4/(4*beta)/YR
	else:
		c0 = a0*(1-e0**2)/e0**(12/19)*(1+121/304 * e0**2)**(-870/2299)
		t = 12*c0**4/(19*beta)
		t *= integrate.quad(eint, 0, e0)[0]
		return t / YR

def coaltime(M1, M2, sig, rho, e0, H = 17.5, emax = 0.99, a0 = 0.04*KPC):
	m = (M1 + M2)*Msun
	m1m2 = M1 * M2 * Msun**2
	e = e0 * (e0<emax) + emax*(e0>=emax)
	Fe = (1-e**2)**(-3.5)*(1+(73./24)*e**2+(37./96)*e**4)
	A = GRA*H*rho/sig
	B = 64*GRA**3*m1m2*m*Fe/(5*SPEEDOFLIGHT**5)
	agw = (B/A)**(1/5) 
	agw = agw * (agw<a0) + a0 * (agw>=a0)
	return (1/(A*agw) - 1/(A*a0)) / YR + time_isobhb(agw, M1, M2, e0) #agw**4/(4*B) / YR
	
def Reff_m(m):
	r1 = 2.95*(m/1e6)**0.596 * PC
	r2 = 34.8*(m/1e6)**0.399 * PC
	return r1 * (r1>r2) + r2 * (r1<=r2)
	
def rho_inf(lm, gamma = 1.5):
	lms = (lm/1e9/0.49)**(1/1.16) * 1e11
	Reff = Reff_m(lms)
	r0 = Reff/0.75 * (2**(1/(3-gamma))-1)
	rinf = r0/((lms/2/lm)**(1/(3-gamma))-1)
	lrho = (3-gamma)*lms/(4*np.pi) * r0/(rinf**gamma*(rinf+r0)**(4-gamma)) * Msun
	return lrho

def sigma_inf(lm):
	lsig = (lm/1e9/0.309)**(1/4.38) * 2e2*1e5
	return lsig
	
def tGW(m, m1, m2, e, H = 17.5, gamma = 1.5, a0 = 0.04*KPC, mode=0):
	sig = sigma_inf(m)
	rho = rho_inf(m, gamma)
	if mode==0:
		return coaltime(m1, m2, sig, rho, e, H, a0 = a0)
	else:
		return time_isobhb(a0, m1, m2, e)
	
def gen_Ms(alp = 0.17, M1 = 40, M2 = 140):
	D = M2**(1-alp) - M1**(1-alp)
	N = M1**(1-alp)
	p = np.random.random()
	return (p*D + N)**(1/(1-alp))
	
def gen_q(q1, q2 = 1.0, alp = 0.75):
	p = np.random.random()
	D = q2**alp - q1**alp
	N = q1**alp
	return (p*D + N)**(1/alp)
	
def gen_P(alp = 2.54, p1 = 1.0, p2 = 4.0):
	p = np.random.random()
	D = p2**alp - p1**alp
	N = p1**alp
	P = (p*D + N)**(1/alp)
	return 10**P
	
def binary_gen(Ms = 500, fb = 0.36, Mmin = 40, gamma = 1.5, sk = 0.1, fbh = 0.8, frm = 0.57, mode=0):
	mb = Ms*2*fb/(1+fb) * fbh * frm
	mtot=  0.0
	lm1 = []
	lm2 = []
	la0 = []
	lt = []
	while mtot<mb:
		M1 = gen_Ms()
		q = gen_q(Mmin/M1)
		M2 = q*M1
		M = M1 + M2		
		if M+mtot>Ms:
			M2 = Ms-mtot-M1
			if M2<Mmin:
				break
		P = gen_P()
		a = (P*YR*(GRA*M*Msun)**0.5/(2*np.pi))**(2/3) * sk
		t = tGW(M, M1, M2, 0, gamma=gamma, a0 = a, mode=mode)
		lm1.append(M1)
		lm2.append(M2)
		la0.append(a/AU)
		lt.append(t)
		mtot += M
	return np.array([lt, lm1, lm2, la0])
	
def testb(N = 1000, Ms = 586, gamma = 1.5, sk = 0.01, mode=0):
	out = []
	for i in range(N):
		out.append(binary_gen(Ms, gamma = gamma, sk = sk, mode=mode))
	return np.hstack(out)
	
if __name__ == "__main__":
	gamma, sk = 1.5, 0.01
	Np = 1000
	mode = 0
	test = testb(Np, gamma=gamma, sk=sk, mode=mode)
	lt, lm1, lm2, la0 = test
	Nb = len(lt)
	tmax = TZ(0)/1e9/YR
	y1, y2 = np.min(lm1+lm2)*0.9, np.max(lm1+lm2)*1.2
	plt.figure()
	plt.scatter(lt/1e9, lm1+lm2, c=lm2/lm1, cmap = plt.cm.cool, s=128, alpha=0.7, label=r'{} binaries from {} stellar populations'.format(Nb, Np))#,vmin=0, vmax=1)
	plt.plot([tmax, tmax], [y1, y2], 'k--', label=r'Age of the Universe')
	cb = plt.colorbar()
	cb.ax.set_title('$M_{2}/M_{1}$', y=-0.1)
	plt.ylim(y1, y2)
	plt.legend(loc=2)
	#plt.xscale('log')
	plt.xlabel(r'$t_{\mathrm{GW}}\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$M_{1}+M_{2}\ [\mathrm{M_{\odot}}]$')
	plt.title(r'$\gamma='+str(gamma)+'$, $f_{\mathrm{sk}}='+str(sk)+'$, $m_{\star}=586\ \mathrm{M}_{\odot}$, $f_{\mathrm{B}}=0.36$, $f_{\mathrm{BH}}=0.8$, $f_{\mathrm{rm}}=0.57$')
	plt.tight_layout()
	plt.savefig('ColBHB_test.pdf')
	plt.close()
	
	plt.figure()
	plt.hist(np.log10(la0), 30, density = True, alpha = 0.7)
	plt.xlabel(r'$\log(a_{0}\ [\mathrm{AU}])$')
	plt.ylabel(r'Probability density')
	plt.tight_layout()
	plt.savefig('dis_a0.pdf')
	plt.close()
	
	plt.figure()
	plt.hist(np.log10(lt/1e9), 30, density=True, alpha=0.7)
	plt.xlabel(r'$\log(t_{\mathrm{GW}}\ [\mathrm{Gyr}])$')
	plt.ylabel(r'Probability density')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig('dis_tGW.pdf')
	plt.close()
	
	tlim = TZ(0)/1e9/YR
	ltd = lt[lt/1e9<tlim]/1e9
	lzGW = np.array([ZT(t) for t in np.log10(ltd)])
	plt.figure()
	plt.hist(lzGW, 30, density=True, alpha=0.7)
	plt.xlabel(r'$z_{\mathrm{GW}}$')
	plt.ylabel(r'Probability density')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig('dis_zGW.pdf')
	plt.close()
	
	gamma = 1.5
	M1 = M2 = 100
	la0 = np.geomspace(0.1, 1e2, 10)
	for a0 in la0:
		t = tGW(M1+M2, M1, M2, 0, a0 = a0*AU, gamma = gamma) / 1e9
		print('tGW(a0 = {} AU) = {} Gyr'.format(a0, t))

	"""
	lls = ['-', '--', '-.', ':']
	le = [0, 0.25, 0.5, 0.9]
	z = 15.
	boost = 1e5
	delta = 200.
	sig = Vcir(1e6, z)
	rho = delta*rhom(1/(z+1.)) * boost
	print('Msun pc^-3 mp^-1: {} (cgs)'.format(Msun/PC**3/PROTON))
	print('sigma = {:.3f} km s^-1, rho = {} Msun pc^-3'.format(sig/1e5, rho / Msun * PC**3))
	M1, M2 = 1e2, 1e2
	print('eccentricity: ', le)
	print('tGW:', [coaltime(M1, M2, sig, rho, e)/1e9 for e in le])
	
	gamma = 1.5
	lm = np.geomspace(1e2, 1e7, 100)
	lsig = sigma_inf(lm)
	lrho = rho_inf(lm, gamma)
	lt = [tGW(lm, lm/2, lm/2, e, gamma=gamma)/1e9 for e in le]#[coaltime(lm/2, lm/2, lsig, lrho, e)/1e9 for e in le]
	
	plt.figure()
	plt.loglog(lm, lsig/1e5)
	plt.xlabel(r'$M_{1+2}\ [\mathrm{M}_{\odot}]$')
	plt.ylabel(r'$\sigma_{\mathrm{inf}}\ [\mathrm{km\ s^{-1}}]$')
	plt.tight_layout()
	plt.savefig('siginf_MB.pdf')
	plt.close()
	
	plt.figure()
	plt.loglog(lm, lrho/Msun * PC**3)
	plt.xlabel(r'$M_{1+2}\ [\mathrm{M}_{\odot}]$')
	plt.ylabel(r'$\rho_{\mathrm{inf}}\ [\mathrm{M_{\odot}\ pc^{-3}}]$')
	plt.tight_layout()
	plt.savefig('rhoinf_MB.pdf')
	plt.close()
	
	plt.figure()
	[plt.plot(lm, t, label=r'$e='+str(e)+'$', ls=s) for t, e, s in zip(lt, le, lls)]
	plt.ylabel(r'$t_{\mathrm{GW}}\ [\mathrm{Gyr}]$')
	plt.xlabel(r'$M_{1+2}\ [\mathrm{M}_{\odot}]$')
	plt.xscale('log')
	plt.yscale('log')
	#plt.ylim(1e-2, 1)
	plt.legend()
	plt.tight_layout()
	plt.savefig('tGW_MB.pdf')
	plt.close()
	"""
	
	"""
	ncore = 8
	nbin = 64
	rep = 'miniquasar_fdbk/'
	if not os.path.exists(rep):
		os.makedirs(rep)

	m1, m2 = 40, 1e9
	a1, a2 = 1e-10, 1e4
		
	lmbh = np.array([3, 5, 7, 9])
	ler = [[1e-4, 5], [1e-4, 1e2], [1e-4, 1e2], [5e-3, 1e2]]
	leta = np.geomspace(1e-4, 1e2, 200)
	lmdot = Macc_edd(10**lmbh)
	lspec = [[BH_spec(m, eta*md, Nnu=100) for eta in leta] for m, md in zip(10**lmbh, lmdot)]
	lion = [np.array([x['ion'] for x in y]).T for y in lspec]
	lLW = [np.array([x['LW'] for x in y]).T for y in lspec]
	#print(lion)
	plt.figure(figsize=(10,7))
	ax1 = plt.subplot(221)
	ax2 = plt.subplot(222)
	ax3 = plt.subplot(223)
	ax4 = plt.subplot(224)
	for i in range(len(lmbh)):
		Lion = lion[i][0]
		Nion = lion[i][1]
		Nion_mod = model_Nion(10**lmbh[i], leta)
		LLW = lLW[i][0]
		NLW = lLW[i][1]
		fLion = curve_fit(linear, np.log(leta), np.log(Lion))
		fNion = curve_fit(linear, np.log(leta), np.log(Nion))
		sel = (leta > ler[i][0]) * (leta < ler[i][1])
		fLLW = curve_fit(linear, np.log(leta[sel]), np.log(LLW[sel]))
		fNLW = curve_fit(linear, np.log(leta[sel]), np.log(NLW[sel]))
		print('MBH = 10^{} Msun'.format(lmbh[i]))
		print('Lion fit parameters: slope: {:.3f}, norm: {}'
				.format(fLion[0][0], np.exp(fLion[0][1])))
		print('Nion fit parameters: slope: {:.3f}, norm: {}'
				.format(fNion[0][0], np.exp(fNion[0][1])))
		print('LLW fit parameters: slope: {:.3f}, norm: {}'
				.format(fLLW[0][0], np.exp(fLLW[0][1])))
		print('NLW fit parameters: slope: {:.3f}, norm: {}'
				.format(fNLW[0][0], np.exp(fNLW[0][1])))
		ax1.loglog(leta, Lion, ls=lls[i], label=r'$M_{\mathrm{BH}}=10^'+str(lmbh[i])+'\ \mathrm{M_{\odot}}$')
		ax2.loglog(leta, Nion, ls=lls[i], label=r'$M_{\mathrm{BH}}=10^'+str(lmbh[i])+'\ \mathrm{M_{\odot}}$')
		ax2.loglog(leta, Nion_mod, ls=lls[i], color='gray')#, label=r'$M_{\mathrm{BH}}=10^'+str(lmbh[i])+'\ \mathrm{M_{\odot}}$')
		ax3.loglog(leta, LLW, ls=lls[i], label=r'$M_{\mathrm{BH}}=10^'+str(lmbh[i])+'\ \mathrm{M_{\odot}}$')
		ax4.loglog(leta, NLW, ls=lls[i], label=r'$M_{\mathrm{BH}}=10^'+str(lmbh[i])+'\ \mathrm{M_{\odot}}$')
	ax1.legend()
	#ax2.legend()
	#ax3.legend()
	#ax4.legend()
	ax1.set_xlabel(r'$\eta$')
	ax2.set_xlabel(r'$\eta$')
	ax3.set_xlabel(r'$\eta$')
	ax4.set_xlabel(r'$\eta$')
	ax1.set_ylabel(r'$\log(L_{\mathrm{ion}}\ [\mathrm{erg\ s^{-1}}])$')
	ax2.set_ylabel(r'$\log(N_{\mathrm{ion}}\ [\mathrm{s^{-1}}])$')
	ax3.set_ylabel(r'$\log(L_{\mathrm{LW}}\ [\mathrm{erg\ s^{-1}}])$')
	ax4.set_ylabel(r'$\log(N_{\mathrm{LW}}\ [\mathrm{s^{-1}}])$')
	plt.tight_layout()
	plt.savefig(rep+'eta_MBH.pdf')
	
	
	leta = [1e-4, 1e-2, 1, 100]
	lmr = [[40, 1e7], [40, 1e9], [2e2, 1e9], [1e5, 1e9]]
	lm = np.geomspace(m1, m2, 200)
	lmdot = Macc_edd(lm)
	lspec = [[BH_spec(m, eta*md, Nnu=100) for m, md in zip(lm, lmdot)] for eta in leta]
	lion = [np.array([x['ion'] for x in y]).T for y in lspec]
	lLW = [np.array([x['LW'] for x in y]).T for y in lspec]
	plt.figure(figsize=(10,7))
	ax1 = plt.subplot(221)
	ax2 = plt.subplot(222)
	ax3 = plt.subplot(223)
	ax4 = plt.subplot(224)
	for i in range(len(leta)):
		Lion = lion[i][0]
		Nion = lion[i][1]
		Nion_mod = model_Nion(lm, leta[i])
		LLW = lLW[i][0]
		NLW = lLW[i][1]
		fLion = curve_fit(linear, np.log(lm), np.log(Lion))
		fNion = curve_fit(linear, np.log(lm), np.log(Nion))
		sel = (lm > lmr[i][0]) * (lm < lmr[i][1])
		fLLW = curve_fit(linear, np.log(lm[sel]), np.log(LLW[sel]))
		fNLW = curve_fit(linear, np.log(lm[sel]), np.log(NLW[sel]))
		print('eta = {}'.format(leta[i]))
		print('Lion fit parameters: slope: {:.3f}, norm: {}'
				.format(fLion[0][0], np.exp(fLion[0][1])))
		print('Nion fit parameters: slope: {:.3f}, norm: {}'
				.format(fNion[0][0], np.exp(fNion[0][1])))
		print('LLW fit parameters: slope: {:.3f}, norm: {}'
				.format(fLLW[0][0], np.exp(fLLW[0][1])))
		print('NLW fit parameters: slope: {:.3f}, norm: {}'
				.format(fNLW[0][0], np.exp(fNLW[0][1])))
		ax1.loglog(lm, Lion, ls=lls[i], label=r'$\eta='+str(leta[i])+'$')
		ax2.loglog(lm, Nion, ls=lls[i], label=r'$\eta='+str(leta[i])+'$')
		ax2.loglog(lm, Nion_mod, ls=lls[i], color='gray')#, label=r'$\eta='+str(leta[i])+'$')
		ax3.loglog(lm, LLW, ls=lls[i], label=r'$\eta='+str(leta[i])+'$')
		ax4.loglog(lm, NLW, ls=lls[i], label=r'$\eta='+str(leta[i])+'$')
	#ax1.legend()
	#ax2.legend()
	#ax3.legend()
	ax4.legend()
	ax1.set_xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	ax2.set_xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	ax3.set_xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	ax4.set_xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	ax1.set_ylabel(r'$\log(L_{\mathrm{ion}}\ [\mathrm{erg\ s^{-1}}])$')
	ax2.set_ylabel(r'$\log(N_{\mathrm{ion}}\ [\mathrm{s^{-1}}])$')
	ax3.set_ylabel(r'$\log(L_{\mathrm{LW}}\ [\mathrm{erg\ s^{-1}}])$')
	ax4.set_ylabel(r'$\log(N_{\mathrm{LW}}\ [\mathrm{s^{-1}}])$')
	plt.tight_layout()
	plt.savefig(rep+'MBH_eta.pdf')
	"""
	
	"""
	
	tag = 1
	if tag==0:
		X, Y, Lion, Nion, LLW, NLW = main(m1, m2, a1, a2, ncore=ncore, nbin=nbin, Tmin = 1e4)
		totxt(rep+'X'+'.txt',X,0,0,0)
		totxt(rep+'Y'+'.txt',Y,0,0,0)
		totxt(rep+'Lion'+'.txt',Lion,0,0,0)
		totxt(rep+'Nion'+'.txt',Nion,0,0,0)
		totxt(rep+'LLW'+'.txt',LLW,0,0,0)
		totxt(rep+'NLW'+'.txt',NLW,0,0,0)
	else:
		X = np.array(retxt(rep+'X.txt',nbin,0,0))
		Y = np.array(retxt(rep+'Y.txt',nbin,0,0))
		Lion = np.array(retxt(rep+'Lion.txt',nbin,0,0))
		Nion = np.array(retxt(rep+'Nion.txt',nbin,0,0))
		LLW = np.array(retxt(rep+'LLW.txt',nbin,0,0))
		NLW = np.array(retxt(rep+'NLW.txt',nbin,0,0))
	
	plt.figure()
	plt.plot(X[:,0], Macc_edd(X[:,0]), 'k--', label=r'$\dot{M}_{\mathrm{Edd}}$')
	ctf = plt.contourf(X, Y, np.log10(Lion), 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(L_{\mathrm{ion}}\ [\mathrm{erg\ s^{-1}}])$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	plt.ylabel(r'$\dot{M}_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ yr^{-1}}]$')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig(rep+'LionMap'+'.pdf')
	plt.close()
	
	plt.figure()
	plt.plot(X[:,0], Macc_edd(X[:,0]), 'k--', label=r'$\dot{M}_{\mathrm{Edd}}$')
	ctf = plt.contourf(X, Y, np.log10(Nion), 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(N_{\mathrm{ion}}\ [\mathrm{s^{-1}}])$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	plt.ylabel(r'$\dot{M}_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ yr^{-1}}]$')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig(rep+'NionMap'+'.pdf')
	plt.close()
	
	plt.figure()
	plt.plot(X[:,0], Macc_edd(X[:,0]), 'k--', label=r'$\dot{M}_{\mathrm{Edd}}$')
	ctf = plt.contourf(X, Y, np.log10(LLW), 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(L_{\mathrm{LW}}\ [\mathrm{erg\ s^{-1}}])$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	plt.ylabel(r'$\dot{M}_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ yr^{-1}}]$')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig(rep+'LLWMap'+'.pdf')
	plt.close()
	
	plt.figure()
	plt.plot(X[:,0], Macc_edd(X[:,0]), 'k--', label=r'$\dot{M}_{\mathrm{Edd}}$')
	ctf = plt.contourf(X, Y, np.log10(NLW), 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(N_{\mathrm{LW}}\ [\mathrm{s^{-1}}])$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M_{\mathrm{BH}}\ [\mathrm{M_{\odot}}]$')
	plt.ylabel(r'$\dot{M}_{\mathrm{acc}}\ [\mathrm{M_{\odot}\ yr^{-1}}]$')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig(rep+'NLWMap'+'.pdf')
	plt.close()
	"""
	
	"""
	x1, x2 = 10, 10e3
	y1, y2 = 1e21, 4e23
	MBH = 1e2
	Mdot = 1e-6 #Macc_edd(MBH)*100
	test = BH_spec(MBH, Mdot)
	print('LM luminosities: L = {} erg s^-1, Ndot = {} s^-1'.format(*test['LW']))
	print('Ionization luminosities: L = {} erg s^-1, Ndot = {} s^-1'.format(*test['ion']))

	plt.figure()
	plt.loglog(*test['tot'], label='Total')
	plt.loglog(*test['MCD'], '--', label='MCD')
	plt.loglog(*test['NT'], '-.', label='NT')
	plt.loglog(*test['MCD_orig'], ':', label='MCD original')
	plt.xlim(x1, x2)
	plt.ylim(y1, y2)
	plt.xlabel(r'$h\nu\ [\mathrm{eV}]$')
	plt.ylabel(r'$L_{\nu}\ [\mathrm{erg\ s^{-1}\ Hz^{-1}}]$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep+'spec_test.pdf')
	#plt.show()
	plt.close()
	"""
