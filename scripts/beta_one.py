
import numpy
from srxraylib.plot.gol import plot

def get_fwhm(bins,h):
    tt = numpy.where(h >= max(h) * 0.5)
    if h[tt].size > 1:
        binSize = bins[1] - bins[0]
        fwhm = binSize * (tt[0][-1] - tt[0][0])
        # print('fwhm = %g'%fwhm)
        return fwhm
    else:
        return 0.0

# def phi(x,c,mode=0):
#
#     if mode==0:
#         return (2*c/numpy.pi)**(1.0/4) * numpy.exp(-c*x**2)
#
#     if mode==1:
#         return (2*c/numpy.pi)**(1.0/4) * 2 * numpy.sqrt(c) * x * numpy.exp(-c*x**2)
#
#     if mode==2:
#         return (2*c/numpy.pi)**(1.0/4) / numpy.sqrt(2) * (4 * x - 1) * numpy.exp(-c*x**2)
#
# x = numpy.linspace(-600e-6,600e-6,200)
# sigma = 100e-6
# wavelength = 2e-9
# c = numpy.sqrt(5./16) / sigma**2
#
# # plot(x,numpy.exp(-c*x**2))

# lambda0 = numpy.sqrt(4*numpy.pi/(3.0+numpy.sqrt(5)))
# lambda1 = lambda0 * (2.0/(3.0+numpy.sqrt(5)))
# lambda2 = lambda0 * (2.0/(3.0+numpy.sqrt(5)))**2
# print(lambda0,lambda1,lambda2)
# print(lambda0/lambda0,lambda1/lambda0,lambda2/lambda0)
#
# # plot(x,lambda2*phi(x,c,mode=2)**2)
#
#
# spectral_density  = lambda0 * phi(x,c,mode=0)**2
# spectral_density += lambda1 * phi(x,c,mode=1)**2
# spectral_density += lambda2 * phi(x,c,mode=2)**2
#
# # plot(x*1e6,spectral_density,title="source spectral density",xtitle="x [um]")
#
# source_fwhm_intensity = get_fwhm(x,spectral_density)
#
#
#
# # plot(1e6*x/SIGMA*SIGMAP*118.0,spectral_density)
#
#
#
# # plot(xp*1e6*118.0,spectral_density,title="propagated (118m) spectral density",xtitle="x [um]")
# print("c: ",c,gs.c())

from wofry.propagator.util.gaussian_schell_model import GaussianSchellModel1D

sigma = 100e-6
wavelength = 2e-9

gs = GaussianSchellModel1D(1.0,sigma,sigma)
x = numpy.linspace(-600e-6,600e-6,600)


print("eigenvalues",gs.beta(0),gs.beta(1),gs.beta(2))
print("eigenvalues/eigenvalue0",gs.beta(0)/gs.beta(0),gs.beta(1)/gs.beta(0),gs.beta(2)/gs.beta(0))

sdensity = gs.beta(0) * gs.phi(0, x) ** 2
source_fwhm_intensity = get_fwhm(x,sdensity)
sdensity += gs.beta(1) * gs.phi(1, x) ** 2
sdensity += gs.beta(2) * gs.phi(2, x) ** 2


SIGMA = 67e-6 # source_fwhm_intensity / 2.355
SIGMAP = wavelength / 4 / numpy.pi / SIGMA
print("0th MODE SIGMA, SIGMAP, FWHM [um]= ",1e6*SIGMA,1e6*SIGMAP,1e6*2.35*SIGMA)
print("propagated 0th mode FWHM at 0 m= %f um"%(1e6*SIGMA*2.355))
print("propagated 0th mode FWHM at 118 m= %f um"%(SIGMAP*2.355*118.0*1e6))
print("propagated 0th mode FWHM at 218 m= %f um"%(SIGMAP*2.355*218.0*1e6))


print("propagated FWHM at 0 m= %f um"%(get_fwhm(1e6*x,sdensity)))
print("propagated FWHM at 118 m= %f um"%(get_fwhm(1e6*x/SIGMA*SIGMAP*118.0,sdensity)))
print("propagated FWHM at 218 m= %f um"%(get_fwhm(1e6*x/SIGMA*SIGMAP*218.0,sdensity)))

plot(1e6 * x,sdensity,
     1e6 * x / SIGMA * SIGMAP * 118.0,sdensity,
     1e6 * x / SIGMA * SIGMAP * 218.0,sdensity,
     legend=["source FWHM=%d um"%(get_fwhm(x*1e6,sdensity)),
             "118m FWHM=%d um"%(get_fwhm(1e6*x/SIGMA*SIGMAP*118.0,sdensity)),
             "218m FWHM=%d um"%(get_fwhm(1e6*x/SIGMA*SIGMAP*218.0,sdensity))],
             xtitle="x/um")