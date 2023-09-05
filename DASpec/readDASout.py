#!/usr/bin/env python
# coding=utf-8

# @Author: zhixiang zhang <zzx>
# @Date:   2017-07-27T21:54:18+08:00
# @Email:  zhangzx@ihep.ac.cn
# @Filename: readDASout.py
# @Last modified by:   zzx
# @Last modified time: 2017-07-27T21:55:00+08:00

import os
import sys
import numpy as np
from PyAstronomy import pyasl
from matplotlib.backends.backend_pdf import PdfPages
# from termcolor import colored
import DASpec
from DASpec import *
from tqdm import tqdm


def log2model2(loglst, redshift=0.0):
    """
    convert DASpec fitting log to DASpec.compcontainer2 object
    loglst: string list
    return: DASpec.compcontainer2 object
    """
    model = DASpec.compcontainer2()
    if redshift is not None:
        model.redshift = redshift
    log = loglst[:]
    if '##########' in log[0]:
        log = log[1:]
    fname = log[0]
    model.set_filename(fname)
    log = log[2:]
    unit = float(log[0])
    model.set_unit(unit)
    log = log[2:]
    window = log[0]
    window = window.replace('[', ' ').replace(']', ' ')
    window = [float(i) for i in window.split()]
    window = np.array(window)
    window = window.reshape((len(window) // 2, 2))
    for win in window:
        model.add_window(win[0], win[1])
    log = log[2:]
    while log[0][0] != '#':
        ID, func = log[0].split(' ', 1)
        model.add(eval(func))
        log = log[1:]
    log = log[1:]  # # fix
    while log[0][0] != '#':
        parlst = log[0].split()[1:]
        parlst = [i.split(':')[1] for i in parlst]
        comp = int(parlst[0])
        par = int(parlst[1])
        val = float(parlst[2])
        model.addfix(comp, par, val)
        log = log[1:]
    log = log[1:]  # # tie
    while log[0][0] != '#':
        # print log[0]
        if 'profile' in log[0]:
            ID, other = log[0].split(' ', 1)
            func, other = other.split(' ', 1)
            if 'value:' in other:
                other, val = other.split('value:', 1)
                val = float(val)
            comp1, comp2 = other.split('->')
            compID1 = int(comp1.split(':')[1])
            compID2 = int(comp2.split(':')[1])
            if func.strip() == 'flux_profile':
                model.addtie_flux_profile(compID1, compID2, val)
            elif func.strip() == 'profile':
                model.addtie_profile(compID1, compID2)
        else:
            pre, value = log[0].split('value:')
            value = float(value)
            pre, tietype = pre.split('type:')
            tietype = tietype.strip()
            pre = pre.split(' ', 1)[1]
            part1, part2 = pre.split('->')
            comp1, par1 = part1.split()
            comp2, par2 = part2.split()
            comp1 = int(comp1.split(':')[1])
            par1 = int(par1.split(':')[1])
            comp2 = int(comp2.split(':')[1])
            par2 = int(par2.split(':')[1])
            model.addtie(comp1, par1, comp2, par2, tietype, value)
        log = log[1:]
    while '# result' not in log[0]:
        log = log[1:]
    log = log[1:]
    while 'par of comp' in log[0]:
        ID, par = log[0].split(':')
        ID = int(ID.split()[-1])
        par = [float(i) for i in par.split()]
        model.set_par(ID, par)
        log = log[1:]
    while len(log) > 0 and 'err of comp' not in log[0]:
        log = log[1:]
    while len(log) > 0 and 'err of comp' in log[0]:
        ID, par = log[0].split(':')
        ID = int(ID.split()[-1])
        par = [float(i) for i in par.split()]
        model.set_err(ID, par)
        log = log[1:]
    while len(log) > 0 and 'limit of comp' not in log[0]:  # skip err of comp
        log = log[1:]
#    while len(log) > 0 and 'limit of comp' in log[0]:
#        ID, par = log[0].split(':', 1)
#        ID = int(ID.split()[-1])
#        par = [float(i) for i in par.split()]
#        model.set_parlimit(ID, par)
#        log = log[1:]
    while len(log) > 0 and 'limit of comp' in log[0]:
        par = log[0].split(':')[1]
        par = [float(i) for i in par.split()]
        model.set_parlimit_all(par)
        log = log[1:]
    return model


def DASout2compmodel2(stlst):
    """
    convert a DASpec_GUI out log to a compcontainer2 object.
    return fname, compcontainer2 object.
    """
    model = log2model2(stlst)
    fname = model.filename
    return fname, model


def getlog(fit):
    """
    fit: DASpec.curvefit object
    get the model component and parameter log list, the log format like this
    rucawftboYFyk240591.txt
    # unit:
    1.000000e-13
    # window:
    [4170.000 4260.000] [4430.000 5550.000]
    # model:
    1 powerlaw(5100.0)
    2 line_gaussian(4861.0)
    3 line_gaussian(4861.0)
    4 line_gaussian(4861.0)
    5 line_gaussian(5007.0)
    6 line_gaussian(4959.0)
    7 line_gaussian(4861.0)
    8 line_gaussian(4686.0)
    9 template_spec('fetemplate_no3','gaussian',4434.0,4684.0)
    10 template_spec('fetemplate_no3','gaussian',4434.0,4684.0)
    11 template_spec('ssp_11Gyr_z02','gaussian',4434.0,4684.0)
    # fix:
    # tie:
    1 flux_profile comp:6 -> comp:5 value:3.333000e-01
    2 profile comp:7 -> comp:5
    3 profile comp:8 -> comp:5
    # reduced chisq: 3.023444
    # result:
    par of comp 1: 8.588537e-02 -1.434816e+00
    par of comp 2: 5.230844e+00 1.645516e+04 2.377541e+02
    par of comp 3: 5.549721e+00 4.119392e+03 -3.608210e+02
    par of comp 4: 1.864929e+00 2.632614e+03 2.561207e+03
    par of comp 5: 1.223337e+00 7.750321e+02 2.431640e-06
    par of comp 6: 4.077381e-01 7.750321e+02 2.431640e-06
    par of comp 7: 1.013139e-01 7.750321e+02 2.431640e-06
    par of comp 8: 4.998204e-02 7.750321e+02 2.431640e-06
    par of comp 9: 5.355106e+00 1.927529e+03 1.892953e+03
    par of comp 10: 4.938511e+00 2.403345e+03 -9.218627e+02
    par of comp 11: 7.645054e-01 1.000000e+02 -1.160580e+02
    err of comp 1: 1.038241e-03 3.222571e-02
    err of comp 2: 1.553819e-01 4.969910e+02 1.626463e+02
    err of comp 3: 2.820423e-01 1.284230e+02 9.044870e+01
    err of comp 4: 2.438830e-01 1.296748e+02 8.221022e+01
    err of comp 5: 1.973197e-02 1.141233e+01 0.000000e+00
    err of comp 6: 0.000000e+00 0.000000e+00 0.000000e+00
    err of comp 7: 2.430528e-02 0.000000e+00 0.000000e+00
    err of comp 8: 1.425685e-02 0.000000e+00 0.000000e+00
    err of comp 9: 2.501325e-01 9.660500e+01 4.556856e+01
    err of comp 10: 2.527012e-01 1.608643e+02 7.131541e+01
    err of comp 11: 2.347847e-01 0.000000e+00 5.146854e+01
    limit of comp: 0.00000e+00 1.00000e+02 -4.00000e+00 0.00000e+00
                   0.00000e+00 1.00000e+04 1.00000e+00 3.00000e+04
                   -4.00000e+03 4.00000e+03 0.00000e+00 1.00000e+04
                   1.00000e+00 1.00000e+04 -4.00000e+03 2.00000e+03
                   0.00000e+00 1.00000e+04 1.00000e+00 1.00000e+04
                   -2.00000e+03 4.00000e+03 0.00000e+00 1.00000e+04
                   1.00000e+00 1.00000e+04 -2.00000e+03 2.00000e+03
                   0.00000e+00 1.00000e+04 0.00000e+00 1.00000e+04
                   0.00000e+00 1.00000e+04 1.00000e+02 1.00000e+04
                   -2.00000e+03 4.00000e+03 0.00000e+00 1.00000e+04
                   1.00000e+02 1.00000e+04 -2.00000e+03 4.00000e+03
                   0.00000e+00 1.00000e+04 1.00000e+02 1.00000e+04
                   -2.00000e+03 2.00000e+03
    """
    loglst = []
    header = '#' * 20
    loglst.append(header)
    loglst.append(fit.model.filename)  # filename
    loglst.append('# unit:')
    unittext = '%e' % fit.model.unit
    loglst.append(unittext)  # unit
    loglst.append('# window:')
    windowtext = ''
    for win in fit.model.window:
        windowtext += '[%f %f] ' % (win[0], win[1])
    loglst.append(windowtext)  # window
    loglst.append('# model:')
    for i in range(len(fit.model.funclst)):
        ind = i + 1
        functext = fit.model.func2txt(ind)
        functext = '%d %s' % (ind, functext)
        loglst.append(functext)  # func
    loglst.append('# fix:')  # skip fix part
    for indfix, fixinfo in enumerate(fit.model.fixlst):
        textfix = '%d comp:%d par:%d value:%.3e' % (indfix+1, fixinfo[0], fixinfo[1], fixinfo[2])
        loglst.append(textfix)
    loglst.append('# tie:')
    for i in range(len(fit.model.tielst)):
        ID = i + 1
        tie = fit.model.tielst[i]
        if len(tie) == 4:
            text = ('%d %s comp:%d -> comp:%d value:%e' %
                    (ID, tie[0], tie[1], tie[2], tie[3]))
        elif len(tie) == 3:
            text = '%d %s comp:%d -> comp:%d' % (ID, tie[0], tie[1], tie[2])
        elif len(tie) == 7:
            text = ('%d comp:%d par:%d -> comp:%d par:%d type:%s value:%e' %
                    (ID, tie[1], tie[2], tie[3], tie[4], tie[5], tie[6]))
        loglst.append(text)
    chisqtext = '#reduced chisq: %f' % fit.reduced_chisq()
    loglst.append(chisqtext)
    loglst.append('# result:')
    for i in range(len(fit.model.comps)):
        ID = i+1
        parlst = fit.par_comp(ID)
        text = 'par of comp %d:' % ID
        for val in parlst:
            text += ' %e' % val
        loglst.append(text)
    for i in range(len(fit.model.comps)):
        ID = i+1
        parlst = fit.parerr_comp(ID)
        text = 'err of comp %d:' % ID
        for val in parlst:
            text += ' %e' % val
        loglst.append(text)
    limitlst = fit.par_limit_all()
    text = 'limit of comp:'
    for val in limitlst:
        text += ' %e' % val
    loglst.append(text)
    return loglst


def plot_model(model, fig):
    fig.clf()
    ax1 = fig.add_axes([0.1, 0.35, 0.8, 0.55])
    ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
    ax1.set_title(model.filename)
    wave, flux, err = model.readspectrum_unit()
    lwidth = 1
    ax1.step(wave, flux, linewidth=lwidth, color='black', alpha=0.7)
    fitflux = np.zeros(wave.shape)
    # ncom = len(model.comps)
    ax1.axhline(0, color='black', linewidth=lwidth)

    for func in model.comps:
        par = model.get_par(func)
        tempflux = func.calc(wave, par)
        fitflux = fitflux + tempflux
        arg = np.where(np.abs(tempflux) > 0.001)
        nwave = wave[arg]
        nflux = tempflux[arg]
        ax1.plot(nwave, nflux, linewidth=lwidth)

    ax1.plot(wave, fitflux, color='red', linewidth=lwidth)

    yl, yr = ax1.get_ylim()
    if yl < -0.1:
        yl = -0.05
    for win in model.window:
        ax1.fill_between(win, [yl, yl], [yr, yr], alpha=0.2, color='black')
    ax1.set_ylim(yl, yr)

    freepar = model.freepar_all()
    dim = len(wave) - len(freepar)
    chisq = np.sum(((flux - fitflux) / err)**2) / dim

    text = 'chisq = %f' % chisq
    xp = wave[0]
    yp = np.max(flux)*0.95
    ax1.text(xp, yp, text, color='red', fontsize=14)
    ax1.set_ylim(yl, yr)

    resid = flux - fitflux
    ax2.step(wave, resid, linewidth=lwidth, color='black')
    ax2.axhline(0, color='blue', linestyle='--')
    sewave, seresid = model.select(wave, resid)
    miny = np.min(seresid)
    maxy = np.max(seresid)
    hight = maxy - miny
    marg = 0.309 * hight
    yr = maxy + marg
    yl = miny - marg
    # tmpyl, tmpyr = ax2.get_ylim()
    for win in model.window:
        ax2.fill_between(win, [yl, yl], [yr, yr],
                         alpha=0.2, color='black')
    ax2.set_ylim(yl, yr)
    return ax1, ax2


def plot(fit, wave=None, flux=None, ax=None):
    if wave is None or flux is None:
        wave, flux, err = fit.model.readspectrum_unit()
    import matplotlib
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_axes([0.1, 0.35, 0.8, 0.55])
        ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
    elif isinstance(ax, matplotlib.figure.Figure):
        fig = ax
        fig.clf()
        ax1 = fig.add_axes([0.1, 0.35, 0.8, 0.55])
        ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
    else:
        ax1 = ax
        ax1.cla()
    ax1.set_title(fit.model.filename)
    lwidth = 1
    ax1.step(wave, flux, linewidth=lwidth, color='black', alpha=0.7)
    
    ax1.plot(wave, fit.calc(wave), color='red', linewidth=lwidth)
    ncom = len(fit.model.comps)
    powerflux = fit.calc_comp(1, wave)
    ax1.plot(wave, powerflux, linewidth=lwidth)
    #################   temp add

    fil = open('powerlaw_data.txt', 'w')
    fpow = powerflux / (1+fit.model.redshift)
    fpow = fpow * fit.model.unit
    powwave = wave * (1+fit.model.redshift)
    for ind in range(len(powwave)):
        text = '%.2f  %.6e' % (powwave[ind], fpow[ind])
        fil.write(text+'\n')
    fil.close()

    ######################
    ax1.axhline(0, color='black', linewidth=lwidth)
    for i in range(2, ncom + 1):
        tempflux = fit.calc_comp(i, wave)
        arg = np.where(np.abs(tempflux) > 0.001)
        nwave = wave[arg]
        nflux = tempflux[arg]
        npowflux = powerflux[arg]
        ax1.plot(nwave, nflux, linewidth=lwidth)
    yl, yr = ax1.get_ylim()
    for win in fit.model.window:
        ax1.fill_between(win, [yl, yl], [yr, yr], alpha=0.2, color='black')
    chisq = fit.reduced_chisq()
    text = 'chisq = %f' % chisq
    xp = wave[0]
    yp = np.max(flux)*0.95
    ax1.text(xp, yp, text, color='red', fontsize=14)
    ax1.set_ylim(yl, yr)
    if isinstance(ax, matplotlib.axes._axes.Axes):
        return
    resid = flux - fit.calc(wave)
    ax2.step(wave, resid, linewidth=lwidth, color='black')
    ax2.axhline(0, color='blue', linestyle='--')
    sewave, seresid = fit.model.select(wave, resid)
    miny = np.min(seresid)
    maxy = np.max(seresid)
    hight = maxy - miny
    marg = 0.309 * hight
    yr = maxy + marg
    yl = miny - marg
    tmpyl, tmpyr = ax2.get_ylim()
    for win in fit.model.window:
        ax2.fill_between(win, [tmpyl, tmpyl], [tmpyr, tmpyr],
                         alpha=0.2, color='black')
    ax2.set_ylim(tmpyl, tmpyr)
    if ax is None:
        plt.show()


def read_outfile(fn):
    """
    read DASpec_GUI out file, split the log to individual models.
    for a object, just save its last fit result log.
    return: list [logobj1, logobj2, logobj3,...]
    """
    objnamelst = set()
    logdic = dict()
    alllst = [i.strip() for i in open(fn)]
    alllst = alllst[1:]

    objtmp = []
    logdic = {}
    if len(alllst) > 10000:
        alllst = tqdm(alllst)
    for line in alllst:
        if '##########' in line:
            if len(objtmp) > 0:
                objname = objtmp[0]
                logdic[objname] = objtmp
                objnamelst.add(objname)
                objtmp = []
        else:
            objtmp.append(line)
    if len(objtmp) > 0:
        objname = objtmp[0]
        logdic[objname] = objtmp
        objnamelst.add(objname)
        objtmp = []
    retlst = [logdic[name] for name in objnamelst]
    return retlst

    # while len(alllst) > 0:
    #     obj = []
    #     while len(alllst) > 0 and '#' * 10 not in alllst[0]:
    #         obj.append(alllst[0])
    #         alllst = alllst[1:]
    #     alllst = alllst[1:]
    #     objname = obj[0]
    #     if objname not in objnamelst:
    #         objnamelst.append(objname)
    #     logdic[objname] = obj
    # retlst = [logdic[name] for name in objnamelst]
    # return retlst


def read_outfile_noreject_repeat(fn):
    """
    read DASpec_GUI out file, split the log to individual models.
    for a object, just save its last fit result log.
    return: list [logobj1, logobj2, logobj3,...]
    """
    objnamelst = []
    loglst = []
    alllst = [i.strip() for i in open(fn)]
    alllst = alllst[1:]
    while len(alllst) > 0:
        obj = []
        while len(alllst) > 0 and '#' * 10 not in alllst[0]:
            obj.append(alllst[0])
            alllst = alllst[1:]
        alllst = alllst[1:]
        objname = obj[0]
        if objname not in objnamelst:
            objnamelst.append(objname)
        loglst.append(obj)
    return loglst


def read_spec(fn):
    """
    read text spectrum
    return: wave, flux, err
    """
    filetype = os.path.splitext(fn)[1]
    if filetype == '.fits':
        wave, flux = pyasl.read1dFitsSpec(fn)
        rerr = np.sqrt(wave*flux)
        arg = np.where(wave > 5000)[0]
        con = 0.01*flux[arg[0]]/rerr[arg[0]]
        err = con*rerr
    else:
        data = np.loadtxt(fn)
        wave = data[:, 0]
        flux = data[:, 1]
        err = data[:, 2]
    return wave, flux, err


def combine_flux(lst):
    """
    combine component flux in lst, just sum(lst)
    """
    return sum(lst)


def combine_err(lst):
    """
    combine error
    """
    tmplst = np.array(lst)
    err = np.sum(tmplst*tmplst)**0.5 / len(tmplst)
    return err


def get_chisq(log):
    for line in log:
        if 'chisq' in line:
            val = float(line.strip().split(':')[-1])
            return val


def get_result(log, ID, index=None):
    """
    get fitting result from log.
    log: fitting log
    type: string list
    ID: the ID of component in fitting model
    type: int
    index: the component parameter[index]
        if index is None, return all parameter of component in a list
    type: int or None, default = None
    return: fitting result
    type: float or float list
    """
    nlog = log[:]
    while 'par of comp' not in nlog[0]:
        nlog = nlog[1:]
    result = {}
    while 'par of comp' in nlog[0]:
        iID, par = nlog[0].split(':')
        iID = int(iID.split()[-1])
        par = [float(i) for i in par.split()]
        result[iID] = par
        nlog = nlog[1:]
    if index is None:
        return result[ID]
    else:
        return result[ID][index]


def get_result_err(log, ID, index=None):
    """
    get fitting result error from log.
    log: fitting log
    type: string list
    ID: the ID of component in fitting model
    type: int
    index: the component parameter[index] error
        if index is None, return all parameter error of component in a list
    type: int or None, default = None
    return: fitting result error
    type: float or float list
    """
    nlog = log[:]
    while 'err of comp' not in nlog[0]:
        nlog = nlog[1:]
    result = {}
    while len(nlog) > 0 and 'err of comp' in nlog[0]:
        iID, par = nlog[0].split(':')
        iID = int(iID.split()[-1])
        par = [float(i) for i in par.split()]
        result[iID] = par
        nlog = nlog[1:]
    if index is None:
        return result[ID]
    else:
        return result[ID][index]


def main():
    """
    run as a script， used to plot logfile.
    parameter 1: log file name
    """
    import matplotlib.pyplot as plt
    logfn = sys.argv[1]
    if len(sys.argv) > 2:
        redshift = float(sys.argv[2])
    else:
        redshift = 0.0
    loglst = read_outfile(logfn)
    plt.ion()
    pdf = PdfPages('fitresult.pdf')
    fig = plt.figure()
    for log in loglst:
        print(log)
        model = log2model2(log, redshift)
        fit = DASpec.curvefit()
        fit.set_model2(model)
        fit.lmfit()
        plot(fit, ax=fig)
        fig.canvas.draw()
        fig.show()
        pdf.savefig(fig)
        raw_input('press any key to continue')
    pdf.close()


if __name__ == '__main__':
    main()
