import matplotlib.pyplot as plt

from ocelot.gui import *
from ocelot import *
from pylab import *
import numpy as np
import cmath
from scipy.constants import c as cLight

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def OptimizeLattice(Dists,Quads,MonitorDictionary,dbzTotal, varsTotal, beam, MSquareMax, SigmaLevel, k, InitialLattice='Parallel',Initialize=True):
    
    # d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11 = [Drift(l=length, eid=f'd{ii+1:.0f}') for ii, length in enumerate(d_lengths)]
    
    QuadLength = Quads[0].l
    num_quads = len(Dists)-1
    quad_eids = [f'q{i+1}' for i in range(num_quads)]
    # quads = [None] * num_quads

    if Initialize == True:
        for i in range(num_quads):
            # if i == 1:
            #     f = 1/(1/(d_lengths[i]) + 1/(d_lengths[i+1])) /2
            if InitialLattice == 'FOFO':
                f = 1/(1/(Dists[i].l) + 1/(Dists[i+1].l)) /2
            elif InitialLattice == 'Parallel':
                f = 1/(1/(Dists[i].l) )
            elif InitialLattice == 'Relay':
                if i == 0:
                    f = 1/(1/(Dists[0].l) + 2/(Dists[1].l)) 
                else:
                    if i%2 == 1:
                        f = 0.5 * Dists[i].l # collimating
                    else:
                        f = 0.5 * Dists[i+1].l # focusing
                
            # quads[i] = Quadrupole(l=QuadLength, k1=1/(f*QuadLength), eid=quad_eids[i])
            Quads[i].k1 = 1/(f*QuadLength)
    
    # m1, mSmall, mBig = Monitors


    # print(quads)

    # q1, q2, q3, q4, q5, q6, q7, q8, q9, q10 = quads

    # dbzTotal = (m1,d1,mSmall,q1,d2,mSmall,q2,d3,mBig,q3,d4,mBig,q4,d5,mSmall,q5,d6,mSmall,q6,d7,mSmall,q7,d8,mSmall,q8,d9,mBig,q9,d10,mBig,q10,d11,mBig)
    # dbzTotal = (m1,d1,q1,d2,mSmall,q2,d3,mBig,q3,d4,mBig,q4,d5,mBig,q5,d6,mBig,q6,d7,mBig,q7,d8,mBig,q8,d9,mBig,q9,d10,mBig,q10,d11,mBig)
    # dbzTotal = (m1,d1,q1,d2,mSmall,q2,d3,mBig,q3,d4,mBig,q4,d5,mBig,q5,d6,mBig,q6,d7,mBig,q7,d8,mBig,q8,d9,mBig,q9,d10,mBig,q10,d11,mBig)
    # varsTotal = [q1,q2,q3,q4,q5,q6,q7,q8,q9,q10]

    
    emittance = beam.emit_x

    latTotal = MagneticLattice(dbzTotal)
    tw0 = Twiss(beam)

    # max_rBig = (1/SigmaLevel) * mirror_radiusBig
    # max_rSmall = (1/SigmaLevel) * mirror_radiusSmall
    # max_rWindow = (1/SigmaLevel) * window_radius
    
    print('emittance', (MSquareMax*(1/(2*k)))*1e6, 'mm mrad')
    
    emittance_max = (MSquareMax*(1/(2*k)))
    # max_betaBig = max_rBig**2 / emittance_max
    # max_betaSmall = max_rSmall**2 / emittance_max
    # max_betaWindow = max_rWindow**2 / emittance_max
    # max_beta_300GHz = max_rBig**2 / emittance
    Monitors = np.fromiter(MonitorDictionary.keys(),dtype = type(MonitorDictionary.keys()))
    # print(Monitors)
    MonitorSizes = np.fromiter(MonitorDictionary.values(), dtype=float)
    # print(MonitorSizes)
    MonitorBetas = (MonitorSizes/SigmaLevel)**2 / emittance
    MonitorBetasMax = (MonitorSizes/SigmaLevel)**2 / emittance_max
    # print(Monitors)
    # print(MonitorBetas)


    # constrTotal = {m1:{'Dx':0.0, 'Dxp':0.0, 'beta_x':beam.beta_x, 'alpha_x':beam.alpha_x, 'beta_y':beam.beta_y, 'alpha_y':beam.alpha_y},mSmall:{'beta_x':['<',max_betaSmall]},mBig:{'beta_x':['<',max_betaBig]}}
    # constrTotal = {m1:{'Dx':0.0, 'Dxp':0.0, 'beta_x':beam.beta_x, 'alpha_x':beam.alpha_x, 'beta_y':beam.beta_y, 'alpha_y':beam.alpha_y},
    #                 mSmall:{'beta_x':['<',max_betaSmall]},mBig:{'beta_x':['<',max_betaBig]},mWindow:{'beta_x':['<',max_betaWindow]}}
    constrTotal = {Monitors[0]:{'Dx':0.0, 'Dxp':0.0, 'beta_x':beam.beta_x, 'alpha_x':0.0, 'beta_y':beam.beta_y, 'alpha_y':0.0}}
    constrTotal.update({Monitor:{'beta_x':['<',beta]} for Monitor, beta in zip(Monitors[1:],MonitorBetasMax[1:])})
    # constrTotal = {Monitors[ii]:{'Dx':0.0, 'Dxp':0.0, 'beta_x':MonitorBetas[ii], 'alpha_x':0.0, 'beta_y':MonitorBetas[ii], 'alpha_y':0.0} for ii, Monitor in enumerate(Monitors)}
    if varsTotal:
        match(latTotal, constrTotal, varsTotal, tw0, max_iter=20000)

    return latTotal, tw0, Quads


def PlotLattice(lat,beam, MonitorDictionary,k,MSquareMax,SigmaLevel,extension,quads_tot,FixedQuadsList=[],mWindow=None):
    f = plt.figure(figsize=(8, 15))
    ax = f.add_subplot(411)

    tws = twiss(lat, Twiss(beam), nPoints = 1000)
    sigma_0 = (beam.beta_x*beam.emit_x)**0.5

    pos = 0

    # pos = d1.l 
    ii = 0

    QuadList = []
    PosList = []
    fList = []
    sizeList = []
    betaList = []

    lambda_thz = 2*np.pi/k

    emittance = beam.emit_x
    MSquare = emittance*2*k

    emittance_max = (MSquareMax*(1/(2*k)))
    # emittance = (MSquare*(1/(2*k)))
    # beta_MaxBig = (MirrorSizeBig/SigmaLevel)**2 / emittance_max
    # beta_MaxSmall = (MirrorSizeSmall/SigmaLevel)**2 / emittance_max
    # beta_MaxWindow = (window_radius/SigmaLevel)**2 / emittance_max
    # # print(beta_MaxBig, beta_MaxSmall)
    # beta_Big = (MirrorSizeBig/SigmaLevel)**2 / emittance
    # beta_Small = (MirrorSizeSmall/SigmaLevel)**2 / emittance
    # beta_Window = (window_radius/SigmaLevel)**2 / emittance
    Monitors = np.fromiter(MonitorDictionary.keys(),dtype = type(MonitorDictionary.keys()))
    MonitorSizes = np.fromiter(MonitorDictionary.values(), dtype=float)
    # print(MSquareMax)
    MonitorBetas = (MonitorSizes/SigmaLevel)**2 / emittance
    # print(lat.)
    # print(lat.sequence[2].k1)

    for element in lat.sequence:
        pos += element.l
        # print(element)
        if element.__class__ == Quadrupole:
            print("Position magnet",ii, pos, 'm')  
            print( "Focal length magnet ", ii, ' = ' , 1/(element.k1*element.l), 'm')
            print( "Strength magnet ", ii, ' = ' , element.k1, '1/m^2')
            plt.axvline(x=pos, color='r', alpha=0.5, linestyle='--')
            # PosList.append(pos+d1.l)
            PosList.append(pos)
            fList.append(1/(element.k1*element.l))
            # QuadList.append(vars[ii//2])
            ii += 1
            QuadList.append(element)
            # PosList.append(pos)
            # fList.append(element.k1)
            # if ii == 1:
            #     sizeList.append(MirrorSizeSmall*1e3 )
            #     betaList.append(beta_Small)
        elif element.__class__ == Monitor:
            if element != Monitors[0] and element != mWindow:
                MonitorSize = MonitorDictionary[element]
                sizeList.append(MonitorSize*1e3)
                betaList.append( (MonitorSize/SigmaLevel)**2 / emittance )
            # if element == mSmall:
            #     sizeList.append(MirrorSizeSmall*1e3 )
            #     betaList.append(beta_Small)
            # elif element == mBig:
            #     sizeList.append(MirrorSizeBig*1e3)
            #     betaList.append(beta_Big)
            elif element == mWindow:
                MonitorSize = MonitorDictionary[element]
                sizeList.append(MonitorSize*1e3)
                betaList.append( (MonitorSize/SigmaLevel)**2 / emittance )
                PosList.append(pos)
                fList.append(0)
        # pos += element.l

    # sizeListPrint = sizeList[0:len(sizeList)-1]
    # print(lat)
    # print(PosList)
    # print(sizeList)
    # print(betaList)

    np.savetxt(f'flist_{extension}.txt', np.column_stack((PosList, fList, sizeList)), delimiter=',', header='PosList (m), fList (m), mirror size (mm)', comments='')

    sizeList = [MonitorSizes[0]*1e3] + sizeList
    betaList = [MonitorBetas[0]] + betaList
    PosList = [0] + PosList
    PosList = np.array(PosList)
    # print( 'PosList', PosList)
    fList = np.array(fList)
    # sizeList = sizeList[0:len(fList)]
    sizeList = np.array(sizeList)
    # betaList = betaList[0:len(fList)]

    betaList = np.array(betaList)
    

    s = [p.s for p in tws]
    # PosList = np.append(PosList, s[-1])
    # fList = np.append(fList, 0)
    # print(len(s))
    ax.set_xlim(0, lat.totalLen)

    # plt.title(''.join(['Optimised transport line for ',f'{cLight/lambda_thz/1e12:.1f}',r' THz and $\sigma_{r,0}$ =',f'{sigma_0*1e3:.1f}','mm']))
    plt.title(''.join(['Optimised transport line for ',f'{cLight/lambda_thz/1e9:.1f}',r' GHz, $\sigma_{r,0}$ = ',f'{sigma_0*1e6:.1f}',r' µm and $M^2$ = ',f'{MSquare:.2f}']))
    p1, = plt.plot(s, [p.beta_x for p in tws], lw=2.0,label=r"$\beta_x$")
    #p2, = plt.plot(s, [p.beta_y for p in tws], lw=2.0)
    plt.grid(True)
    # plt.axhline(y=max_beta, color='b', alpha=0.7, linestyle='--',label='Max beta')
    plt.plot(PosList,betaList*MSquare/MSquareMax, color='k', alpha=0.5, linestyle='--',label=fr'Max beta for $M^2$={MSquareMax:.1f}')
    plt.plot(PosList,betaList, color='b', alpha=0.7, linestyle='--',label=fr'Max beta $M^2$={MSquare:.1f}')
    # plt.plot([PosList[1],s[-1]],[max_beta,max_beta], color='k', alpha=0.5, linestyle='--',label='__nolegend__')
    # plt.plot([0,PosList[1],PosList[1]],[max_beta_300GHz*(2/3)**2,max_beta_300GHz*(2/3)**2,max_beta_300GHz], color='b', alpha=0.7, linestyle='--',label='Max beta 300GHz')
    # plt.plot([PosList[1],s[-1]],[max_beta_300GHz,max_beta_300GHz], color='b', alpha=0.7, linestyle='--',label='__nolegend__')
    # plt.axhline(y=min_beta, color='b', alpha=0.7, linestyle='--',label='Min beta')
    plt.ylabel(r'$\beta_x$ [m]', rotation=0, labelpad=20,color='b')
    # plt.legend(loc='upper left')
    plt.legend(loc=(1.1,0.6))

    ax.twinx()
    p3,=plt.plot(s, [-2* p.alpha_x for p in tws], 'r', alpha=0.7,lw=2.0,label=r"$\beta_x'$")

    #plt.legend([p1,p2,p3], [r'$\beta_x$',r'$\beta_y$', r'$D_x$'])
    # plt.legend(loc='upper right')
    plt.ylabel(r"$\beta_x'$", rotation=0,labelpad=10,color='r')
    plt.legend(loc=(1.1,0.45))

    ax2 = f.add_subplot(412)
    # plot_lattice(lat, ax2, alpha=0.5)

    # add beam size (arbitrary scale)

    scale = 1

    sig_x = scale * np.array([np.sqrt(p.beta_x*emittance) for p in tws]) # 0.03 is for plotting same scale
    # sig_y = scale * np.array([np.sqrt(p.beta_y*beam.emit_y) for p in tws])

    x = scale * np.array([p.x for p in tws])
    # y = scale * np.array([p.y for p in tws])

    plt.plot()
    plt.plot(s, (x + sig_x)*1e3, color='r', lw=2.0,label='Beam size')
    plt.plot(s, (x-sig_x)*1e3, color='r', lw=2.0,label='__nolegend__')

    TempList = (x + SigmaLevel*sig_x)*1e3

    plt.plot(s, (x + SigmaLevel*sig_x)*1e3, color='r', alpha=0.6, lw=2.0,label=fr'Beam size ${SigmaLevel:.1f} \sigma$')
    plt.plot(s, (x-SigmaLevel*sig_x)*1e3, color='r', alpha=0.6, lw=2.0,label='__nolegend__')

    # plt.plot(s, sig_y, color='g', lw=2.0)
    # plt.plot(s, -sig_y, color='g', lw=2.0)

    #f=plt.figure()
    plt.plot(s, x*1e3, 'r--', lw=2.0)

    # plt.axhline(y=mirror_radius*1e3, color='k', alpha=0.7, linestyle='--',label='Mirror edge')
    # plt.axhline(y=-mirror_radius*1e3, color='k', alpha=0.7, linestyle='--',label='__nolegend__')
    # plt.plot([0,PosList[1],PosList[1]],[mirror_radius*2/3 *1e3,mirror_radius*2/3 *1e3,mirror_radius *1e3], color='k', alpha=0.7, linestyle='--',label='Beam pipe')
    # plt.plot([PosList[1],s[-1]],[mirror_radius*1e3,mirror_radius*1e3], color='k', alpha=0.7, linestyle='--',label='__nolegend__')
    # plt.plot([0,PosList[1],PosList[1]],[-mirror_radius*2/3 *1e3,-mirror_radius*2/3 *1e3,-mirror_radius *1e3], color='k', alpha=0.7, linestyle='--',label='__nolegend__')
    # plt.plot([PosList[1],s[-1]],[-mirror_radius*1e3,-mirror_radius*1e3], color='k', alpha=0.7, linestyle='--',label='__nolegend__')

    plt.plot(PosList,sizeList, color='k', alpha=0.7, linestyle='--',label='Mirror size')
    plt.plot(PosList,-sizeList, color='k', alpha=0.7, linestyle='--',label='__nolegend__')

    for ii in range(len(quads_tot)):
        # plt.axvline(x=PosList[ii], color='r', alpha=0.7, linestyle='--',label='Magnet')
        # plt.axhline(y=fList[ii]*1e3, color='g', alpha=0.7, linestyle='--',label='Focal length')
        # plt.arrow(PosList[ii], mirror_radius*1e3 /2, fList[ii],0, head_width=2, head_length=2, fc='g', ec='g')
        # plt.arrow(PosList[ii], TempList[find_nearest(s, PosList[ii])], fList[ii],0, head_width=4, head_length=0.2, fc='g', ec='g')
        if quads_tot[ii] in FixedQuadsList:
            plt.plot( [PosList[ii]-2*fList[ii], PosList[ii],PosList[ii]+2*fList[ii]], [0,TempList[find_nearest(s, PosList[ii])],0], 'm',linestyle='--',alpha=0.9, lw=1.0)
            plt.plot( [PosList[ii]-2*fList[ii], PosList[ii],PosList[ii]+2*fList[ii]], [0,-TempList[find_nearest(s, PosList[ii])],0], 'm',linestyle='--',alpha=0.9, lw=1.0)
        else:
            if fList[ii] < 5:
                plt.plot( [PosList[ii]-2*fList[ii], PosList[ii],PosList[ii]+2*fList[ii]], [0,TempList[find_nearest(s, PosList[ii])],0], 'g',linestyle='--',alpha=0.9, lw=1.0)
                plt.plot( [PosList[ii]-2*fList[ii], PosList[ii],PosList[ii]+2*fList[ii]], [0,-TempList[find_nearest(s, PosList[ii])],0], 'g',linestyle='--',alpha=0.9, lw=1.0)
        # plt.arrow(PosList[ii], mirror_radius*1e3 /2, -fList[ii],0, head_width=2, head_length=2, fc='g', ec='g')
    #plt.plot(s, y, 'r--', lw=2.0)

    ax2.set_xlim(0, lat.totalLen)

    plt.grid(True)
    plt.plot([],[],'g--',alpha=0.9, lw=1.0,label='Optimised mirror')
    if FixedQuadsList:
        plt.plot([],[],'m--',alpha=0.9, lw=1.0,label='Fixed mirror')
    plt.ylabel('x [mm]',rotation=0,labelpad=20)
    plt.xlabel('s [m]')
    plt.legend(loc=(1.05,0.2))
    plt.savefig(f'twiss_{extension}.png', dpi=300, bbox_inches='tight')
    plt.close('all')

    lat.save_as_py_file('lat_{extension}.py', tws0=tws, remove_rep_drifts=False, power_supply=False)
    # betaSecondMirror = betaList[np.argmin(np.abs(sList-d1.l-d2.l))]
    # alphaSecondMirror = alphaList[np.argmin(np.abs(sList-d1.l-d2.l))]
    # with open(f'lat_{extension}.py', 'a') as file:
    #     file.write(f"betaSecondMirror = {betaSecondMirror}\n")
    #     file.write(f"alphaSecondMirror = {alphaSecondMirror}\n")

    # betaThirdMirror = betaList[np.argmin(np.abs(sList-d1.l-d2.l-d3.l))]
    # alphaThirdMirror = alphaList[np.argmin(np.abs(sList-d1.l-d2.l-d3.l))]
    # with open(f'lat_{extension}.py', 'a') as file:
    #     file.write(f"betaThirdMirror = {betaThirdMirror}\n")
    #     file.write(f"alphaThirdMirror = {alphaThirdMirror}\n")
    return PosList, fList, sizeList, betaList

def PlotAcceptance(lat,beam, PosList, sizeList, SigmaLevel, k,SigmaList,MsquareList,extension):
    AcceptedList = np.zeros((len(SigmaList),len(MsquareList)))

    ii = 0
    for sigma_0 in SigmaList:
        jj = 0
        for Msquare in MsquareList:
            beam = Beam()
            beam.E = 16.0
            emittance = Msquare*(1/(2*k))
            beam.emit_x = emittance
            beam.emit_y = emittance

            beam.beta_x = sigma_0**2 /emittance
            beam.beta_y = sigma_0**2 /emittance
            beam.alpha_x = 0.0
            beam.alpha_y = 0.0
            beam.disp_x = 0.0
            beam.disp_y = 0.0
            beam.disp_dx = 0.0
            beam.disp_dy = 0.0

            tws = twiss(lat, Twiss(beam), nPoints = 1000)

            sig_x = np.array([np.sqrt(p.beta_x*emittance) for p in tws]) 
            s = [p.s for p in tws]
            Accepted = 1
            for kk in range(len(PosList)):
                index = find_nearest(s, PosList[kk])
                if sig_x[index]*SigmaLevel > sizeList[kk]*1e-3:
                    Accepted = 0
                    break
            AcceptedList[ii,jj] = Accepted

            jj += 1
        ii += 1
    
    fig, ax = plt.subplots()
    
    # print(AcceptedList)
    contour = ax.contourf(MsquareList,SigmaList*1e3, AcceptedList,levels=200, cmap='viridis')
    fig.colorbar(contour)
    # ax.set_xticklabels(['']+MsquareList)
    # ax.set_yticklabels(['']+SigmaList)
    ax.set_xlabel('M^2')
    ax.set_ylabel(r'$\sigma_{r,0}$ [mm]')
    plt.title('Acceptance')
    plt.savefig(f'Acceptance_{extension}.png', dpi=300, bbox_inches='tight')
    plt.close('all')
    return


