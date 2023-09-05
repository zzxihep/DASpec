#!/usr/bin/python

"""
Author: Pu Du @ IHEP
Date: 2016.07
"""

import os
import sys
import numpy as np
from scipy import optimize
from PyAstronomy.pyasl import unred
from . import carray
from . import swigDASpec
from . import dasreadspec

template_PATH = os.getenv("DASpec_TEMPLATE_PATH")

line_centers_air = {
        "Hbeta": 4861.325,
        "Hgamma": 4340.464,
        "[OIII]4959": 4958.911,
        "[OIII]5007": 5006.843,
        "[OIII]4363": 4363.210,
        "HeII": 4685.71,
        "MgII": 2797.963
        }

line_centers_vacuum = {
        "Hbeta": 4862.683,
        "Hgamma": 4341.684,
        "[OIII]4959": 4960.295,
        "[OIII]5007": 5008.240,
        "[OIII]4363": 4364.436,
        "HeII": 4687.02,
        "MgII": 2798.788
        }

fitwindows = [
    [2470, 2625],
    [2675, 2755],
    [2855, 3010],
    [3625, 3645],
    [4170, 4260],
    # [4430, 4770],
    # [5080, 5550],
    [4430, 5550],
    [6050, 6200],
    [6890, 7010]]

fitwindows_continuum = [
    [2470, 2625],
    [2675, 2755],
    [2855, 3010],
    [3625, 3645],
    [4170, 4260],
    # [4430, 4770],
    # [5080, 5550],
    [4430, 5550],
    [6050, 6200],
    [6890, 7010]]

model_name = [
    "powerlaw",
    "powerlaw_d",
    "balmer_continuum",
    "template_spec",
    "template_spec_reddened",
    "line_gaussian",
    "line_lorentzian",
    "line_dgaussian",
    "line_gh4",
    "ccm_reddening"
]

model_name_p = [
    "powerlaw",
    "powerlaw_d",
    "balmer_continuum",
    "template_spec_gaussian",
    "template_spec_dgaussian",
    "template_spec_lorentzian",
    "template_spec_gh4",
    "template_spec_reddened_gaussian",
    "template_spec_reddened_dgaussian",
    "template_spec_reddened_lorentzian",
    "template_spec_reddened_gh4",
    "line_gaussian",
    "line_lorentzian",
    "line_dgaussian",
    "line_gh4",
    "ccm_reddening"
]

tie_name = [
    "tie",
    "profile",
    "flux_profile"
]

def str_to_model(arrays):

    """ string array to model """

    model = compcontainer()

    for i in arrays:
        text = 'model.add(' + i + ')'
        #print text
        eval(text)
    return model

class array(object):

    """ array object used for transmission with C++ """

    def __init__(self, p):

        """
        initialize object

        p is python array or numpy array
        """

        self.len = len(p)
        self.c_array = carray.newarray(self.len)
        self.py_array = p

        for i in range(self.len):
            # print(type(p[i]))
            carray.setelement(self.c_array, i, p[i])

    def carray(self):
        return self.c_array

    def pyarray(self):
        p = []
        for i in range(self.len):
            p.append(carray.getelement(self.c_array, i))
        return np.array(p)

    def __del__(self):
        carray.delarray(self.c_array)


class compcontainer(object):

    """
    component container

    int npar: number of parameters (all of the models)
    int nmodel: number of models
    component *model: array containing pointers to model

    info: print information of the component

    calc: calculate the model including all of the components
    calc_totpar: calculate the model including all of the components and long par array
    int num: length of x and y
    double *x: wavelength array (Angstrom)
    double *y: flux array (erg/s/cm2/Angstrom)
    double *p: parameter array

    addfix: add fix
    c: component
    p: parameter
    val: value

    addtie: add tie
    c: component
    p: parameter
    ct: target component
    pt: target parameter
    type: tie type ("ratio":0, "offset":1)
    val: value (ratio or offset)
    """



    def __init__(self):
        self.container = swigDASpec.compcontainer()
        self.comps = []

    def add(self, func):
        self.comps.append(func)  # put it into comps to avoid memory cleanning
        self.container.add(self.comps[-1].func)  # put it into c class

    def info(self):
        self.container.info()

    def calc_totpar(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.container.calc_totpar(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.container.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

    def pars2l(self, p):
        ptot = [0.0] * self.container.npar
        cptot = array(ptot)
        cp = array(p)
        self.container.pars2l(cp.carray(), cptot.carray())
        return np.array(cptot.pyarray())

    def parl2s(self, ptot):
        p = [0.0] * self.container.netnpar()
        cp = array(p)
        cptot = array(ptot)
        self.container.parl2s(cptot.carray(), cp.carray())
        return np.array(cp.pyarray())

    def clean(self):
        self.container.clean()
        comps = []

    def addfix(self, c, p, val):
        self.container.addfix(c, p, val)

    def addtie(self, c, p, ct, pt, t, val):
        self.container.addtie(c, p, ct, pt, t, val)

    def addtie_profile(self, c, ct):
        self.container.addtie_profile(c, ct)

    def addtie_flux_profile(self, c, ct, val):
        self.container.addtie_flux_profile(c, ct, val)

    def netnpar(self):
        return self.container.netnpar()

    def calc_comp(self, n, x, par):

        """ calculate component n using parameter """

        n = n - 1

        if n < 0 or n >= len(self.comps):
            print("Error! n is invalid!")
            sys.exit()

        pfit = par

        # print self.model.comps[n].calc()
        nnp = 0
        if n == 0:
            # print 'npar', self.model.comps[n].func.npar
            p = pfit[0: self.comps[n].func.npar]
        else:
            for i in range(n):
                nnp += self.comps[i].func.npar
            # print np
            p = pfit[nnp: nnp + self.comps[n].func.npar]
        # print p

        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.comps[n].func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())


class compcontainer2(compcontainer):
    def __init__(self):
        super(compcontainer2, self).__init__()
        self.filename = ''
        self.funclst = []
        self.parlst = []
        self.parlimitlst = []
        self.locklst = []
        self.window = []
        self.tielst = []  # format ['tiefuncname', ID1, ID2, val(option)]
        self.fixlst = []  # format [compID, p, val]
        self.unit = 1.0
        self.errlst = []
        # self.redshift = 0.0
        # self.redshift = 0.05500
        # self.redshift = 0.158339
        # self.redshift = 0.129
        # self.redshift = 0.0
        self.redshift = 0.0

    def index(self, func):
        """
        return the component func subscript+1 in compcontainer2.
        if func not in this compcontainer2, return None.
        if func type is int, return func.
        """
        if type(func) is int:
            return func
        for i in range(len(self.funclst)):
            if self.funclst[i] == func:
                return i+1
        return None

    def add(self, func, par=None, parlimit=None):
        super(compcontainer2, self).add(func)
        npar = func.func.npar
        self.funclst.append(func)
        ipar = [1.0 for i in range(npar)]
        self.parlst.append(ipar)
        if par is not None:
            self.set_par(func, par)
        iparlimit = [float('-inf') for i in range(2*npar)]
        for i in range(len(iparlimit)//2):
            iparlimit[2*i+1] = float('inf')
        if parlimit is not None:
            for i, val in enumerate(parlimit):
                iparlimit[i] = val
        self.parlimitlst.append(iparlimit)
        ilocklst = [False for i in range(npar)]
        self.locklst.append(ilocklst)
        ierrlst = [0 for i in range(npar)]
        self.errlst.append(ierrlst)
        return func

    def addfix(self, fc, p, val):
        """
        fc should be func or int(comp subscript)
        other same to addfix in compcontainer.
        """
        ind = self.index(fc)
        super(compcontainer2, self).addfix(ind, p, val)
        self.locklst[ind-1][p-1] = True
        fixinfo = [ind, p, val]
        self.fixlst.append(fixinfo)

    def addtie(self, fc, p, fct, pt, t, val):
        """
        fc, fct should be func or int(comp subscript)
        others same to addfix in compcontainer.
        t: type of tie('0' or 'ratio': ratio, '1' or 'offset': offset)
        type: string or int
        """
        c = self.index(fc)
        ct = self.index(fct)
        t = str(t)
        if t == '0':
            t = 'ratio'
        elif t == '1':
            t = 'offset'
        super(compcontainer2, self).addtie(c, p, ct, pt, t, val)
        self.locklst[c-1][p-1] = True
        self.tielst.append(['tie', c, p, ct, pt, t, val])

    def addtie_profile(self, fc, fct):
        c = self.index(fc)
        ct = self.index(fct)
        super(compcontainer2, self).addtie_profile(c, ct)
        for i in range(1, len(self.locklst[c-1])):
            self.locklst[c-1][i] = True
        self.tielst.append(['profile', c, ct])

    def addtie_flux_profile(self, fc, fct, val):
        c = self.index(fc)
        ct = self.index(fct)
        super(compcontainer2, self).addtie_flux_profile(c, ct, val)
        for i in range(len(self.locklst[c-1])):
            self.locklst[c-1][i] = True
        self.tielst.append(['flux_profile', c, ct, val])

    def set_filename(self, fn):
        """
        set the spectrum data file name.
        """
        self.filename = fn

    def set_unit(self, val):
        """
        set the unit of flux
        """
        self.unit = val

    def set_par(self, fc, par):
        """
        set parameter of component fc
        fc should be func or int
        """
        c = self.index(fc)
        siz = min(len(self.parlst[c-1]), len(par))
        for i in range(siz):
            self.parlst[c-1][i] = par[i]

    def set_par_all(self, par):
        """
        set all free parameter.
        par: free parameter list or a all parameter list
        """
        u = 0
        if len(par) < len(self.parlst):
            for i in range(len(self.parlst)):
                for j in range(len(self.parlst[i])):
                    if self.locklst[i][j] is False:
                        self.parlst[i][j] = par[u]
                        u += 1
        else:
            for i in range(len(self.parlst)):
                for j in range(len(self.parlst[i])):
                    self.parlst[i][j] = par[u]
                    u += 1

    def set_parlimit(self, fc, parlimit):
        """
        set parameter limit of component fc
        fc should be func or int
        parlimit: list
        """
        c = self.index(fc)
        self.parlimitlst[c-1] = parlimit

    def set_parlimit_all(self, parlimit):
        """
        set parameter limit.
        parlimit: all free parameter limit list
            [par1l, par1u, par2l, par2u, par3l, par3u...]
        """
        u = 0
        # print 'run set_parlimit_all'
        # print len(parlimit)
        # print 'len par func lst = ', len(self.funclst)
        for i in range(len(self.parlimitlst)):
            for j in range(len(self.parlimitlst[i])):
                if self.locklst[i][j//2] is False:
                    self.parlimitlst[i][j] = parlimit[u]
                    u += 1

    def freepar_all(self):
        """
        get all free parameter in compcontainer2
        return: list
        """
        freelst = []
        for i, parcomp in enumerate(self.parlst):
            for j, v in enumerate(parcomp):
                if self.locklst[i][j] is not True:
                    if j != 0:
                        freelst.append(v)
                    else:
                        # freelst.append(v*self.unit)
                        freelst.append(v)
        return freelst

    def parlimit_all(self):
        """
        get all free parameter limit in compcontainer2
        return: list
        """
        # print 'run parlimit_all'
        limitlst = []
        for i, parlcomp in enumerate(self.parlimitlst):
            for j, v in enumerate(parlcomp):
                if self.locklst[i][j//2] is not True:
                    limitlst.append(v)
        # print len(limitlst)
        # print limitlst
        return limitlst

    def add_window(self, w1, w2):
        """
        add the fit window.
        """
        self.window.append([w1, w2])

    def select(self, wave, flux, err=None):
        """
        select the wave, flux, err in fit windows
        wave, flux, err should be numpy.array
        return: selected wave, flux, err
        type: array, array, array
        """
        arg = np.zeros(wave.shape, dtype=bool)
        for win in self.window:
            w1, w2 = win
            arg = arg | ((wave > w1) & (wave < w2))
        nwave = wave[arg]
        nflux = flux[arg]
        if err is not None:
            nerr = err[arg]
            return nwave, nflux, nerr
        return nwave, nflux

    def readspectrum(self, fn=None):
        """
        read spectrum and return wave, flux, err
        fn: spectrum file name
            default=None, return default fn spectrum saved in self.filename
        return wave(numpy.array), flux(numpy.array), err(numpy.array)
        """
        if fn is None:
            fn = self.filename
        # print fn
        wave, flux, err = dasreadspec.readspec(fn)
        if hasattr(self, 'ebv') and self.ebv is not None:
            # print('ebv = ', self.ebv)
            flux = unred(wave, flux, self.ebv)
            err = unred(wave, err, self.ebv)
        flux = np.zeros(flux.shape) + flux
        err = np.zeros(err.shape) + err
        arg = np.where(err > 0)
        wave = wave[arg]
        flux = flux[arg]
        err = err[arg]
        arg = np.where(flux > 0)
        wave = wave[arg]
        flux = flux[arg]
        err = err[arg]
        return wave, flux, err

    def readspectrum_unit(self, fn=None):
        """
        read spectrum and return wave, flux, err
            the flux and err are divided by self.unit
        fn: spectrum file name
            default=None, return default fn spectrum saved in self.filename
        return wave(numpy.array),
               flux(numpy.array)/self.unit,
               err(numpy.array)/self.unit
        """
        # add temp code
        # redshift = 0.158339
        wave, flux, err = self.readspectrum(fn)
        wave /= (1+self.redshift)
        flux = flux / self.unit
        err = err / self.unit
        flux *= (1+self.redshift)
        err *= (1+self.redshift)
        # print 'wave type = ', wave.dtype
        # print 'flux type = ', flux.dtype
        # print 'error type = ', err.dtype
        return wave, flux, err

    def func2txt(self, ind):
        """
        convert a func in self.funclst to string
            exmaple: 'line_gaussian(4861.0)'
        ind: func subscript in self.funclst, base from 1
        return: string, line_gaussian(4861.0), powerlaw(5100.0) etc.
        """
        def funcpar(func):
            """
            get func par
            powerlaw: keyword = func.ref
            gaussian: keyword = func.center
            template_spec: keyword =
                templatename, func.kernel, func.flux_llim, func.flux_rlim
            return: string or string tuple
            type: string or (string, string...)
            """
            # if type(func) == type(line_gaussian()):
            if isinstance(func, line_gaussian):
                return str(func.func.center)
            # if type(func) == type(powerlaw()):
            if isinstance(func, line_dgaussian):
                return str(func.func.center)
            if isinstance(func, powerlaw):
                return str(func.func.ref)
            if isinstance(func, powerlaw_d):
                return str(func.func.ref)
            # if type(func) == type(template_spec('fetemplate_no3')):
            if isinstance(func, template_spec):
                tempname = func.templatename
                kernelname = func.func.kernel
                llim = str(func.func.flux_llim)
                rlim = str(func.func.flux_rlim)
                return tempname, kernelname, llim, rlim
            if isinstance(func, line_gh4):
                return str(func.func.center)
            return None

        func = self.funclst[ind-1]
        par = funcpar(func)
        funcname = func.func.name
        # if type(func) == type(template_spec('fetemplate_no3')):
        if isinstance(func, template_spec):
            text = "%s('%s','%s',%s,%s)" % tuple([funcname]+list(par))
            return text
        else:
            if par is None:
                text = "%s()" % funcname
                return text
            else:
                text = "%s(%s)" % (funcname, par)
                return text

    def set_err(self, fc, par):
        """
        set err of par
        """
        ind = self.index(fc)
        self.errlst[ind-1] = par

    def get_err(self, fc):
        """
        get component fc fitting parameter err list.
        fc: func or int
        """
        ind = self.index(fc)
        return self.errlst[ind-1]

    def get_par(self, fc):
        """
        get component fc fitting parameter.
        fc: func or int
        """
        ind = self.index(fc)
        return self.parlst[ind-1]


class curvefit(object):

    """
    fit model (compcontainer) to the data
    """

    def __init__(self):
        self.curvefit = swigDASpec.curvefit()
        self.model = 0
        self.cx = 0
        self.cy = 0
        self.cerr = 0
        self.cp = 0
        self.fitdone = 0
        self.limits = []

    def set_data(self, x, y, err):
        self.cx = array(x)
        self.cy = array(y)
        self.cerr = array(err)
        self.curvefit.setdata(self.cx.len, self.cx.carray(), self.cy.carray(), self.cerr.carray())

    def set_model(self, m):
        self.model = m
        self.curvefit.setmodel(self.model.container)

    def set_model2(self, m):
        self.set_model(m)
        self.set_init(m.freepar_all())
        self.set_limit_tot(m.parlimit_all())
        wave, flux, err = m.readspectrum_unit()
        nwave, nflux, nerr = m.select(wave, flux, err)
        self.set_data(nwave, nflux, nerr)

    def set_init(self, p):
        self.cp = array(p)
        self.curvefit.setinitp(self.cp.carray())

    def set_limit(self, n, p, limit):
        """
        n: index of parameter
        p: limit
        limit: 0-lower, 1-upper
        """
        self.limits.append([n, p, limit])
        self.curvefit.setlimit(n, p, limit)

    def set_limit_tot(self, p):
        if self.cp.len != len(p) // 2:
            print('error!! number of limits is not equal to the number of parameters!')
            sys.exit()
        for i in range(len(p) // 2):
            #print i + 1, i * 2, i * 2 + 1
            self.set_limit(i + 1, p[i * 2], 0)
            self.set_limit(i + 1, p[i * 2 + 1], 1)

    def lmfit(self, nitermax = 200):
        if self.model == 0:
            print("Error! please set model for fitting!")
            sys.exit()
        if self.cp == 0 or self.cp.len != self.model.netnpar():
            print("Error! please set correct initial parameters for fitting!")
            sys.exit()
        if self.cx == 0 or self.cy == 0 or self.cerr == 0:
            print("Error! please set data for fitting!")
            sys.exit()
        self.curvefit.lmfit(nitermax)
        self.fitdone = 1

    def siman(self, ntmax = 50000, ninner = 5000, jump = 0.2,
            ninit = 2000, Tratio = 0.001, delta = 1.0e-6, nstable = 20):
        if self.model == 0:
            print("Error! please set model for fitting!")
            sys.exit()
        if self.cp == 0 or self.cp.len != self.model.netnpar():
            print("Error! please set correct initial parameters for fitting!")
            sys.exit()
        if self.cx == 0 or self.cy == 0 or self.cerr == 0:
            print("Error! please set data for fitting!")
            sys.exit()
        self.curvefit.siman(ntmax, ninner, jump, ninit, Tratio, delta, nstable)
        self.fitdone = 1

    def mix_fit(self, ntmax = 50000, ninner = 200, jump = 0.2,
            ninit = 500, Tratio = 0.01, delta = 1.0e-5, nstable = 20):
        if self.model == 0:
            print("Error! please set model for fitting!")
            sys.exit()
        if self.cp == 0 or self.cp.len != self.model.netnpar():
            print("Error! please set correct initial parameters for fitting!")
            sys.exit()
        if self.cx == 0 or self.cy == 0 or self.cerr == 0:
            print("Error! please set data for fitting!")
            sys.exit()
        self.curvefit.mix_fit(ntmax, ninner, jump, ninit, Tratio, delta, nstable)
        self.fitdone = 1

    def info(self):
        self.curvefit.info()

    def calc(self, x):

        """ calculate model using the fitted parameterr """

        if self.fitdone == 0:
            print("Error! please do the fitting first!")
            sys.exit()

        p = self.par()
        return self.model.calc(x, p)

    def calc_comp(self, n, x):

        """ calculate component n using the fitted parameterr """

        n = n - 1

        if self.fitdone == 0:
            print("Error! please do the fitting first!")
            sys.exit()

        if n < 0 or n >= len(self.model.comps):
            print("Error! n is invalid!")
            sys.exit()

        pfit = self.par_tot()

        #print self.model.comps[n].calc()
        nnp = 0
        if n == 0:
            #print 'npar', self.model.comps[n].func.npar
            p = pfit[0: self.model.comps[n].func.npar]
        else:
            for i in range(n):
                nnp += self.model.comps[i].func.npar
            #print np
            p = pfit[nnp: nnp + self.model.comps[n].func.npar]
        #print p

        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.model.comps[n].func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

    def par(self):
        """ return final parameters """
        p = []
        for i in range(self.curvefit.npout):
            p.append(carray.getelement(self.curvefit.pout, i))
        return np.array(p)

    def parerr(self):
        """ return error of parameters """
        p = []
        for i in range(self.curvefit.npout):
            p.append(carray.getelement(self.curvefit.perrout, i))
        return np.array(p)

    def par_tot(self):
        """ return final parameters (tot) """
        p = []
        for i in range(self.curvefit.npout_tot):
            p.append(carray.getelement(self.curvefit.pout_tot, i))
        return np.array(p)

    def parerr_tot(self):
        """ return error of parameters (tot) """
        p = []
        for i in range(self.curvefit.npout_tot):
            p.append(carray.getelement(self.curvefit.perrout_tot, i))
        return np.array(p)

    def par_comp(self, ind):
        """return final parameters of nth element, ind begin from 1
            if ind < 0 or ind greater than element number, return None
        """
        nn = ind - 1
        if nn < 0 or nn >= len(self.model.comps):
            print("par_comp ERROR!!! ind is invalid!")
            return None
        fromm = 0
        for i in range(nn):
            fromm += self.model.comps[i].func.npar
        parnum = self.model.comps[nn].func.npar
        ret = []
        for i in range(fromm, fromm+parnum):
            ret.append(carray.getelement(self.curvefit.pout_tot, i))
        return np.array(ret)

    def parerr_comp(self, ind):
        """return final parameters error of nth element, ind begin from 1
            if ind < 0 or ind greater than element number, return None
        """
        nn = ind - 1
        if nn < 0 or nn >= len(self.model.comps):
            print("par_comp ERROR!!! ind is invalid!")
            return None
        fromm = 0
        for i in range(nn):
            fromm += self.model.comps[i].func.npar
        parnum = self.model.comps[nn].func.npar
        ret = []
        for i in range(fromm, fromm+parnum):
            ret.append(carray.getelement(self.curvefit.perrout_tot, i))
        return np.array(ret)

    # def limit_comp(self, ind):
    #     """return the limit of nth component, ind begin from 1
    #         if ind < 0 or ind greater than element number, return None
    #     """
    #     nn = ind - 1
    #     if nn < 0 or nn >= len(self.model.comps):
    #         print "par_comp ERROR!!! ind is invalid!"
    #         return None
    #     fromm = 0
    #     for i in range(nn):
    #         fromm += self.model.comps[i].func.npar
    #     parnum = self.model.comps[nn].func.npar
    #     limit = []
    #     for i in range(fromm, fromm+parnum):
    #         for flag, v in enumerate(self.limits):
    #             if v[0] == i+1 and v[2] == 0:
    #                 limit.append(v[1])
    #                 break
    #             if flag == len(self.limits)-1:
    #                 limit.append(float('-inf'))
    #         for flag, v in enumerate(self.limits):
    #             if v[0] == i+1 and v[2] == 1:
    #                 limit.append(v[1])
    #                 break
    #             if flag == len(self.limits)-1:
    #                 limit.append(float('inf'))
    #     # print 'put out limit list'
    #     # for v in self.limits:
    #         # print v[0], v[1], v[2]
    #     return limit

    def par_limit_all(self):
        """return all limit
        type: list
        """
        sortlimits = sorted(self.limits)
        sortlimits = [i[1] for i in sortlimits]
        sortlimits = [[sortlimits[2*i], sortlimits[2*i+1]]
                      for i in range(len(sortlimits)//2)]
        # print sortlimits
        sortlimits = [sorted(i) for i in sortlimits]
        retlst = []
        for comp in sortlimits:
            for v in comp:
                retlst.append(v)
        return retlst

    def chisq(self):
        """ return chisq """
        return self.curvefit.chisq

    def reduced_chisq(self):
        """ return reduced chisq """
        return self.curvefit.reduced_chisq

    def DOF(self):
        """ return DOF """
        return self.curvefit.DOF

    def status(self):
        """ return fit status """
        return self.curvefit.status

    def iternum(self):
        """ return iteration number """
        return self.curvefit.iternum


class line_gaussian(object):

    """
    Gaussian line profile

    center: center of the line

    info: print information of the component

    calc: calculate the component
    p[0]: flux (erg/s/cm2)
    p[1]: FWHM (km/s)
    p[2]: shift (km/s)
    """

    def __init__(self, center = 4861.0):
        self.func = swigDASpec.line_gaussian(center)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class line_dgaussian(object):

    """
    double Gaussian line profile

    center: center of the line

    info: print information of the component

    calc: calculate the component
    p[0]: flux (erg/s/cm2)
    p[1]: FWHM (km/s)
    p[2]: shift (km/s)
    p[3]: FWHM (km/s)
    p[4]: shift (km/s)
    p[5]: flux ratio (1/(1+2))
    """

    def __init__(self, center = 4861.0):
        self.func = swigDASpec.line_dgaussian(center)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class line_lorentzian(object):

    """
    Lorentzian line profile

    center: center of the line

    info: print information of the component

    calc: calculate the component
    p[0]: flux (erg/s/cm2)
    p[1]: FWHM (km/s)
    p[2]: shift (km/s)
    """

    def __init__(self, center = 4861.0):
        self.func = swigDASpec.line_lorentzian(center)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class powerlaw(object):

    """
    power law
    a * (x / c)^b

    ref: c

    info: print information of the component

    calc: calculate the component
    p[0]: flux (erg/s/cm2/A)
    p[1]: power index
    """

    def __init__(self, ref = 5100.0):
        self.func = swigDASpec.powerlaw(ref)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())


class powerlaw_d(object):

    """
    double power law
    y = A * (x / ref)^b1 for x < x_edge
    y = B * (x / ref)^b2 for x >= x_edge
    so B = A * (x_edge/ref)^(b1-b2)
    But in fact, we use flux_ref, b1, b2, x_edge as free parameters

    ref: c

    info: print information of the component

    calc: calculate the component
    p[0]: flux (erg/s/cm2/A)
    p[1]: power index1
    p[2]: power index2
    p[3]: wave edge of index1 and index2
    """

    def __init__(self, ref=5100.0):
        self.func = swigDASpec.powerlaw_d(ref)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class line_gh4(object):

    """
    Gauss-Hermite polynomials line profile

    center: center of the line

    info: print information of the component

    calc: calculate the component
    p[0]: flux (erg/s/cm2)
    p[1]: sigma (km/s)
    p[2]: shift (km/s)
    p[3]: skewness
    p[4]: kurtosis
    """

    def __init__(self, center = 4861.0):
        self.func = swigDASpec.line_gh4(center)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class template_spec(object):

    """
    template convolved by a certain kernel

    p: parameters of kernel
    """

    cx = 0
    cy = 0

    def __init__(self, filename, kernel = "gaussian", f_llim = 4434.0, f_rlim = 4684.0):
        # read file
        l = open(template_PATH + filename).readlines()
        x = [float(i.split()[0]) for i in l if i[0] != '#']
        y = [float(i.split()[1]) for i in l if i[0] != '#']
        self.cx = array(x)
        self.cy = array(y)
        self.func = swigDASpec.template_spec(self.cx.len, self.cx.carray(), self.cy.carray(),
            kernel, f_llim, f_rlim)
        self.templatename = filename  # appended variable

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class template_spec_reddened(object):

    """
    reddened template convolved by a certain kernel

    p: parameters of kernel
    """

    cx = 0
    cy = 0

    def __init__(self, filename, kernel = "gaussian", f_llim = 4434.0, f_rlim = 4684.0, r = 3.1):
        # read file
        l = open(template_PATH + filename).readlines()
        x = [float(i.split()[0]) for i in l if i[0] != '#']
        y = [float(i.split()[1]) for i in l if i[0] != '#']
        self.cx = array(x)
        self.cy = array(y)
        self.func = swigDASpec.template_spec_reddened(self.cx.len, self.cx.carray(), self.cy.carray(),
            kernel, f_llim, f_rlim, r)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class balmer_continuum(object):

    """
    balmer continuum

    p[0]: flux at Balmer edge
    p[1]: optical depth of Balmer edge
    """

    def __init__(self):
        self.func = swigDASpec.balmer_continuum()

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [0.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())

class ccm_reddening(object):

    """
    ccm reddening

    p[0]: ebv
    """

    def __init__(self, r = 3.1):
        self.func = swigDASpec.ccm_reddening(r)

    def info(self):
        self.func.info()

    def calc(self, x, p):
        y = [1.0 for i in x]
        cx = array(x)
        cy = array(y)
        cp = array(p)
        self.func.calc(len(x), cx.carray(), cy.carray(), cp.carray())
        return np.array(cy.pyarray())


def test():

    from astropy.io import fits as pyfits
    from pyastrolib import astro
    import matplotlib.pyplot as plt

    fit = pyfits.open("spSpec-52368-0881-064.fit")
    #print fit[0].header
    wave0 = fit[0].header["CRVAL1"]
    nwave = fit[0].header["NAXIS1"]
    dwave = fit[0].header["CD1_1"]
    wave = wave0 + dwave * np.arange(nwave)
    wave = 10.0**wave
    data = fit[0].data
    flux = data[0]
    err = data[2]

    index = np.where(err > 0)
    wave = wave[index[0]]
    flux = flux[index[0]]
    err = err[index[0]]
    (wave, flux, err) = np.array((wave, flux, err), dtype = np.double)

    #astro.ccm_unred(wave, flux, 0.112 - 0.085)
    z = 0.153375
    wave = wave / (1.0 + z)
    flux = flux * (1.0 + z)
    err = err * (1.0 + z)

    index = np.where(#((wave >= 2470) & (wave <= 2625)) |
            #((wave >= 2675) & (wave <= 2755)) |
            #((wave >= 2855) & (wave <= 3010)) |
            #((wave >= 3625) & (wave <= 3645)) |
            #((wave >= 4170) & (wave <= 4260)) |
            #((wave >= 4430) & (wave <= 4770)) |
            #((wave >= 5080) & (wave <= 5550)) |
            ((wave >= 4430) & (wave <= 5550)))# |
            #((wave >= 6050) & (wave <= 6200)) |
            #((wave >= 6890) & (wave <= 7010)) )

    wave1 = wave[index[0]]
    flux1 = flux[index[0]]
    err1 = err[index[0]]

    model = compcontainer()
    model.add(powerlaw())
    #model.add(balmer_continuum())
    model.add(template_spec("fetemplate_no3"))
    model.add(line_dgaussian(4861))
    model.add(line_gaussian(4959))
    model.add(line_gaussian(5007))
    model.addtie_flux_profile(5, 4, 3.0)
    model.info()
    #model.comps[0].info()

    fit = curvefit()
    fit.set_model(model)
    fit.set_data(wave1, flux1, err1)
    #fit.set_init([1.0, -2.0, 1.0, 0.2, 100.0, 100.0, 0.0, 1.0, 1000.0, 0.0, 1.0, 200.0, 0, 1, 200, 0, 1, 200, 0])
    fit.set_init([0.79978568690512, -0.067457407191308,# 3.520216701612777, 0.1,
        5.85184328147245, 3000.48748920527652, -800.88713428317556,
        0, 1000.0, 0, 3000, 1000, 0.2, 1, 500, 0])
    #fit.set_limit(4, 0.01, 0)
    #fit.set_limit(4, 2.00, 1)
    #fit.set_limit(6, 10.0, 0)
    #fit.set_limit(6, 10000.0, 1)
    #fit.set_limit(7, -10000.0, 0)
    #fit.set_limit(7, 10000.0, 1)

    #fit.lmfit()

    #print fit.par()
    #print fit.parerr()
    #print fit.par_tot()
    #print fit.parerr_tot()
    #y = fit.calc_comp(1, [1, 2, 3])
    #print fit.parerr()

    fit.set_limit(1, 0, 0)
    fit.set_limit(1, 1000, 1)
    fit.set_limit(2, -3, 0)
    fit.set_limit(2, 0, 1)
    #fit.set_limit(3, 0, 0)
    #fit.set_limit(3, 100, 1)
    #fit.set_limit(4, 0.01, 0)
    #fit.set_limit(4, 0.5, 1)
    fit.set_limit(3, 0.0, 0)
    fit.set_limit(3, 1000, 1)
    fit.set_limit(4, 100.0, 0)
    fit.set_limit(4, 2000, 1)
    fit.set_limit(5, -2000.0, 0)
    fit.set_limit(5, 2000, 1)
    fit.set_limit(6, 0.0, 0)
    fit.set_limit(6, 10000, 1)
    fit.set_limit(7, 1.0, 0)
    fit.set_limit(7, 5000, 1)
    fit.set_limit(8, -2000.0, 0)
    fit.set_limit(8, 2000, 1)
    fit.set_limit(9, 1.0, 0)
    fit.set_limit(9, 5000, 1)
    fit.set_limit(10, -2000.0, 0)
    fit.set_limit(10, 2000, 1)
    fit.set_limit(11, 0.0, 0)
    fit.set_limit(11, 1.0, 1)
    fit.set_limit(12, 0.0, 0)
    fit.set_limit(12, 1e4, 1)
    fit.set_limit(13, 1.0, 0)
    fit.set_limit(13, 5000, 1)
    fit.set_limit(14, -2000.0, 0)
    fit.set_limit(14, 2000, 1)

    #fit.siman()
    fit.mix_fit()
    #fit.lmfit()
    print(fit.par())
    print(fit.parerr())
    #print fit.chisq()

    #fit.lmfit()
    #print fit.par()
    #print fit.parerr()

    plt.plot(wave, flux)
    plt.plot(wave, fit.calc(wave))
    plt.plot(wave, fit.calc_comp(1, wave))
    plt.plot(wave, fit.calc_comp(2, wave))
    plt.plot(wave, fit.calc_comp(3, wave))
    plt.plot(wave, fit.calc_comp(4, wave))
    plt.plot(wave, fit.calc_comp(5, wave))
    #plt.plot(wave, fit.calc_comp(6, wave))
    #plt.show()



    #a = compcontainer()
    #a.add(line_gaussian())
    #a.add(template_spec("fetemplate_no3"))
    #a.add(template_spec_reddened("fetemplate_no3"))
    #a.info()

    #x = np.linspace(4500.0, 5200.0, 1000)
    #da = line_gaussian().calc(x, [1.0, 1000.0, 0.0])

    #b = curvefit()
    #b.set_model(a)
    #b.set_init([10.0, 100.0, -100.0])
    #b.set_data(x, da, np.ones(len(x)))
    #b.info()
    #b.lmfit()
    #b.info()

if __name__ == "__main__":
    test()
