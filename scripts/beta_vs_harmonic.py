import numpy
import scipy.constants as codata

m2ev = codata.h * codata.c / codata.e

def get_sigmas_radiation(photon_energy,undulator_length):
    lambdan = m2ev / photon_energy
    return 2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),0.69*numpy.sqrt(lambdan/undulator_length),

def get_sigmas_ALSU(photon_energy,undulator_length):
    sr, srp = get_sigmas_radiation(photon_energy, undulator_length)

    sx, sz, sxp, szp = 12.1e-6, 14.7e-6, 5.7e-6, 4.7e-6

    Sx = numpy.sqrt(sx ** 2 + sr ** 2)
    Sz = numpy.sqrt(sz ** 2 + sr ** 2)
    Sxp = numpy.sqrt(sxp ** 2 + srp ** 2)
    Szp = numpy.sqrt(szp ** 2 + srp ** 2)

    return Sx,Sz,Sxp,Szp

if __name__ == "__main__":
    from scipy import interpolate


    beta = numpy.linspace(0,10,100)
    # https://www.wolframalpha.com/input/?i=0.99+%3D+1+%2F+%281+%2B+1%2F+%280.5+x%5E2+%2B+x+sqrt%28%28x%2F2%29%5E2+%2B+1%29%29%29
    CF = 1 / (1 + 1/ (0.5 * beta**2 + beta * numpy.sqrt((beta/2)**2 + 1)))
    from srxraylib.plot.gol import plot
    # plot(beta,CF,xtitle="beta",ytitle="CF")

    f = interpolate.interp1d(CF,beta)


    #
    # flexon paper
    #
    energy0 = 230.888 # 806.0
    L = 137 * 0.0288

    for factor in [1,3,5]:
        energy = energy0 * factor
        sr, srp = get_sigmas_radiation(energy, L)
        Sx, Sz, Sxp, Szp = get_sigmas_ALSU(energy, L)
        # print("\n Radiation Sizes L=%f m at E=%f eV: \nSr=%5.3f um \nSrp=%5.3f urad: " % (L, energy, 1e6 * sr, 1e6 * srp))
        CFh = sr * srp / (Sx * Sxp)
        print("CF (horizontal) at E=%f eV is %f, beta: %f, Sx=%f um, Smu=%f um" % \
              (energy, CFh, f(CFh), 1e6 * Sx, 1e6 * Sx * f(CFh) ))

    # h1m0 = numpy.loadtxt("h1m0.txt")
    # h1m1 = numpy.loadtxt("h1m1.txt")
    # plot(h1m0[:, 0], h1m0[:, 1],
    #      h1m1[:, 0], h1m1[:, 1],
    #      h1m0[:, 0], h1m0[:, 1] + h1m1[:, 1],
    #      xtitle="x [um]", ytitle="intensity [a.u.]",
    #      title="CF = 0.91, beta=3.07, Sx = 34 um",
    #      legend=["mode 0", "mode 1", "modes 0+1"],
    #      show = 0,
    #      )
    #
    # h1m0_image = numpy.loadtxt("h1m0_image.txt")
    # h1m1_image = numpy.loadtxt("h1m1_image.txt")
    # plot(h1m0_image[:, 0], h1m0_image[:, 1],
    #      h1m1_image[:, 0], h1m1_image[:, 1],
    #      h1m0_image[:, 0], h1m0_image[:, 1] + h1m1_image[:, 1],
    #      xtitle="x [um]", ytitle="intensity [a.u.]",
    #      legend = ["mode 0", "mode 1", "modes 0+1"],
    #      )


    h5m0 = numpy.loadtxt("h5m0.txt")
    h5m1 = numpy.loadtxt("h5m1.txt")
    plot(h5m0[:, 0], h5m0[:, 1],
         h5m1[:, 0], h5m1[:, 1],
         h5m0[:, 0], h5m0[:, 1] + h5m1[:, 1],
         xtitle="x [um]", ytitle="intensity [a.u.]",
         title="CF = 0.68, beta=1.2, Sx = 18.6 um",
         legend=["mode 0", "mode 1", "modes 0+1"],
         show = 0,
         )

    h5m0_image = numpy.loadtxt("h5m0_image.txt")
    h5m1_image = numpy.loadtxt("h5m1_image.txt")
    plot(h5m0_image[:, 0], h5m0_image[:, 1],
         h5m1_image[:, 0], h5m1_image[:, 1],
         h5m0_image[:, 0], h5m0_image[:, 1] + h5m1_image[:, 1],
         xtitle="x [um]", ytitle="intensity [a.u.]",
         legend = ["mode 0", "mode 1", "modes 0+1"],
         )


