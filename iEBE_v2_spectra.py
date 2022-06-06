#!/usr/bin/env python
"""
    This module setups a environment for later convenient calculations using the
    EbeDBReader class.

    One way to use this module from command line to evaluate a single
    expression, for example:
    uhg.py database_filename "e_2(ed)"

    Another way to use this module by importing all its contents to an
    interactive shell, for example:
    python -ic "from uhg import *"

    In the interactive mode use the "use" function to connect to a database, use
    "h" function to print out a short help, and use the "e" function to evaluate
    an expression.

"""
#import py_compile
#py_compile.compile('EbeCollector.py')
from numpy import *
import EbeCollector
import DBR
import matplotlib.pyplot as plt
#import ROOT
from threading import Thread
Ntreds = 20

first_meth = True

second_meth = True
plot_hist = False
ifrec_yields = False



pt = [0.8, 1.0, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.7, 5.0, 6.0, 8.0]
pt_m = [0.9, 1.1, 1.45, 1.95, 2.45, 2.95, 3.45, 3.95, 4.45, 4.85, 5.5, 7.0]
delta_pt = [0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 1.0, 2.0]
pt_r = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
pt_r_m = [0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25]
Npart = [8.2, 85.5, 314, 524]

#namepp = ["collected_pp.db", "collected_pp1.db", "collected_pp2.db"]
#namepp = ["/home/yoren/bnl/iEBE/iEBE-stable/RESULTS/collected.db"]
#names = ["collected.db", "collected1.db"]#20-60
#names = ["/home/yoren/bnl/iEBE/collected_020.db"]
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS020/collected.db"]#0-20
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS2040/collected.db"]#20-40
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS4060/collected.db"]#4060
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS5060/collected.db"]#5060
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS2060/collected.db"]#2060
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS01/collected.db"]#01
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTSUU/collected.db"]#UU
#names = ["/home/yoren/bnl/iEBE/collected4060.db","/home/yoren/bnl/iEBE/collected4060_1.db"]
#names = ["/home/yoren/bnl/iEBE/collectedUU2.db"]
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/RESULTS/collected.db"]
#names = ["/home/yoren/bnl/iEBE/collected_test.db"]
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/RESULTS/collected.db"]
#names = ["/home/yoren/bnl/iEBE/collected_020.db","/home/yoren/bnl/iEBE/collected4060.db","/home/yoren/bnl/iEBE/collected4060_1.db"]
#names = ["/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS020/collected.db",
#         "/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS4060/collected.db",
#         "/home/yoren/bnl/iEBE/iEBE-stable/NEW/RESULTS2060/collected.db"]

names = ["/home/yoren/bnl/iEBE/iEBE-stable/FINAL/CuAu093/collected6.db"]
for i in range(7, 11):
    names.append("/home/yoren/bnl/iEBE/iEBE-stable/FINAL/CuAu093/collected{}.db".format(i))

#names = ["/home/yoren/bnl/iEBE/iEBE-stable/FINAL/UU093/collected0.db"]
#for i in range(1, 5):
#    names.append("/home/yoren/bnl/iEBE/iEBE-stable/FINAL/UU093/collected{}.db".format(i))



print(names)

def first(name, particle, result, itred, start, stop, normall):
    res = [0, 0, 0, 0, 0, 0, 0]
    res_cen_mult = [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]]
    counts = [0, 0, 0, 0, 0, 0, 0]
    counts_mult = [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]]
    for iname in range(len(name)):
        reader = EbeCollector.EbeDBReader(name[iname])
        Nev = reader.getNumberOfEvents()
        mults0 = reader.getMultiplicities(particleName="charged_eta_hydro")
        mults = reader.getMultiplicities(particleName="charged_eta_hydro", orderBy="N")
        boardes = [0., 0, 0, 0]
        for i in range(1, 5):
            boardes[i-1] = int(Nev*(5-i/1)/5)
        mult_borders = [mults[boardes[0]], max(mults), mults[boardes[1]], mults[boardes[0]],
                        mults[boardes[2]], mults[boardes[1]], mults[boardes[2]], mults[boardes[0]]]
        if itred == 0: print(iname, Nev, mult_borders)
        counts_bad = 1
        ecc = abs(reader.get_V_n(particleName="phi_thermal"))
        ecc_av = sum(ecc) / len(ecc)
        ecc1 = sorted(ecc)
        rms = 0.0
        for iecc in ecc:
            rms = rms + (iecc-ecc_av)**2
        rms = sqrt(rms) / len(ecc)
        if itred == 0: print("eccen = ", ecc_av, rms/ecc_av)
        mults_phi = reader.getMultiplicities(particleName=particle)

        for ievent in range(1+int(Nev*start), int(Nev*stop) + counts_bad):

            try:

                #mults_phi_inner = reader.getInterpretedSpectraForOneEvent(event_id=ievent, particleName=particle, pTs=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.6])
                mults_phi_inner = mults_phi[ievent-1]

                v2_iner = abs(reader.getInterpretedComplexDifferentialFlowForOneEvent(particleName=particle,
                                                                                   pTs=[0.5, 1.0, 1.5, 2.0, 2.5,
                                                                                        3.0, 4.6], order=2,
                                                                                   event_id=ievent))
                #print(reader.getDifferentialFlowDataForOneEvent(event_id=ievent, particleName=particle, order=2, pT_range=None, where="", orderBy="pT"))
                v2_iner = v2_iner * mults_phi_inner
                res += v2_iner
                counts += mults_phi_inner

                mult = mults0[ievent-1]
                for iecc in range(4):
                    if mult > mult_borders[2*iecc] and mult <= mult_borders[2*iecc + 1]:
                        res_cen_mult[iecc] += v2_iner
                        counts_mult[iecc] += mults_phi_inner

            except:
                #print(ievent)
                counts_bad += 1

    #print(res / counts)
    #for iecc in range(4):
    #    res_cen_mult[iecc] = res_cen_mult[iecc] / counts_mult[iecc]
    result[itred] = res_cen_mult
    normall[itred] = counts_mult

def second(name, particle):
    try:
        v2 = 0
        yields = 0
        Ncolls = 0
        for i in range(len(name)):
            reader = EbeCollector.EbeDBReader(name[i])
            ncolls_inner = reader.getNumberOfEvents()
            Ncolls += ncolls_inner
            inname = "<v_2([0.5,1.0,1.5,2.0,2.5,3.0,4.6])({})>".format(particle)
            res1 = reader.evaluateExpression(inname)
            v2 += ncolls_inner * res1[0]
            inname = "<dN/dydpT([0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 7.25])({})>".format(particle)
            res1 = reader.evaluateExpression(inname)
            yields += ncolls_inner * res1[0]

        yields = [yields[i]/Ncolls for i in range(len(yields))]
        print(v2 / Ncolls, yields)
        return yields
    except:
        print(name)
        return -9999



try:

    if second_meth:
        yield_AB = second(names, "phi_thermal")
        # yield_pp = second(namepp, "phi_thermal")
        # yield_AB = [yield_AB[i]/yield_pp[i]/Npart[1] for i in range(len(yield_AB))]
        # print ("RAB", pt_r_m)
        # print(yield_AB)

    if first_meth:

        threads = [None] * Ntreds
        results = [None] * Ntreds
        normall = [None] * Ntreds

        for itred in range(len(threads)):
            threads[itred] = Thread(target=first, args=(names, "phi_thermal", results, itred,
                                                        (0.+itred)/Ntreds, (itred+1.)/Ntreds, normall))
            threads[itred].start()

        print(Ntreds, " trends have started")

        for i in range(len(threads)):
            threads[i].join()


        print("yolo ", results[0][0])
        print("yolo ", normall[0][0])
        res = results[0]
        norm = normall[0]
        for i in range(1, Ntreds):
            for j in range(4):
                res[j] = res[j] + results[i][j]
                norm[j] = norm[j] + normall[i][j]
        for j in range(4):
            res[j] = res[j]/norm[j]
        print("shiiiiiiiiiiiiiiiiiiiiiit", res)
        #yield_AB = second(names, "phi_thermal")
        # yield_pp = second(namepp, "phi_thermal")
        # yield_AB = [yield_AB[i]/yield_pp[i]/Npart[1] for i in range(len(yield_AB))]
        # print ("RAB", pt_r_m)




except:
    print("lox")
