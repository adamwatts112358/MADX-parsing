import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

def parsetwiss(filename):
    names = []
    keywords = []
    s = []
    betx = []
    alphax = []
    mux = []
    bety = []
    alphay = []
    muy = []
    dx = []
    dy = []
    L = []

    with open(filename,'r') as f:
        lines = f.readlines()
        del lines[:48]
        for line in lines:
            names.append(line.split()[0].strip('\"'))
            keywords.append(line.split()[1].strip('\"'))
            s.append(float(line.split()[2]))
            betx.append(float(line.split()[3]))
            alphax.append(float(line.split()[4]))
            mux.append(float(line.split()[5]))
            bety.append(float(line.split()[6]))
            alphay.append(float(line.split()[7]))
            muy.append(float(line.split()[8]))
            dx.append(float(line.split()[15]))
            dy.append(float(line.split()[17]))
            L.append(float(line.split()[34]))

    return names,keywords,s,betx,alphax,mux,bety,alphay,muy,dx,dy,L

def draw_quad(s_start,length,name,color):
    plt.gca().add_patch(mpatches.Rectangle((s_start, -0.025), length, .05, facecolor=color,zorder=10))
    plt.text(s_start+length/2,-0.12,name,horizontalalignment='center',fontsize=7,rotation='vertical',color=color)
def draw_dipole(s_start,length,name,color):
    plt.gca().add_patch(mpatches.Rectangle((s_start, -0.025), length, .05, facecolor=color,zorder=10))
    plt.text(s_start+length/2,-0.12,name,horizontalalignment='center',fontsize=7,rotation='vertical',color=color)
def draw_corrector(s_start,length,name,color):
    plt.gca().add_patch(mpatches.Rectangle((s_start, -0.025), length, .05, facecolor=color,zorder=10))
    plt.text(s_start+length/2,0.35,name,horizontalalignment='center',fontsize=7,rotation='vertical',color=color)
def draw_instrumentation(s_start,length,name,color):
    plt.axvline(x=s_start,ymin=0.5,ymax=0.75,color=color,linestyle='--')
    plt.text(s_start+length/2,0.75,name,horizontalalignment='center',fontsize=7,rotation='vertical',color=color)

def plotMarker(name,height):
    s_val=s[names.index(name)]
    plt.axvline(x=s_val,color='red',linestyle='--')
    plt.text(s_val+0.1,height,name,color='red',rotation=90)

def plot(s,betx,bety,dx,dy,names,L):
    #s_Hbpms = []
    #for i in range(len(names)):
    #    if names[i].startswith('HP'):
    #        s_Hbpms.append(float(s[i]))
    #        plotMarker(names[i],1.0)
    #s_Vbpms = []
    #for i in range(len(names)):
    #    if names[i].startswith('VP'):
    #        s_Vbpms.append(float(s[i]))
    #        plotMarker(names[i],1.0)
    ## Measured dispersion
    #eta_H = [-0.287195409405,-0.575898406261,-1.49376789931,-1.73196671656,1.2329552703,1.29401256206,-3.05638229222,0.914100524424,1.02716958324,2.79984116098,1.39100068807,1.29652520781,1.08948319788]
    #eta_V = [0.273627122347,0.69273643371,-1.04551189723,-0.501524092015,-1.02415440834,-2.18449421642,-5.80270409869,3.05261332359,0.720124272402,0.158296682349,0.692987698285]

    plt.figure(figsize=(16,8))
    plt.suptitle('NuMI Beamline',fontsize=22)

    plt.subplot(311)
    plt.plot(s,[0]*len(s),'k-')
    plt.xlim(0,s[-1])
    plt.ylim(-1,1)
    plt.gca().set_xlim(left=0)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.gca().set_axis_bgcolor('#f1f1f1')
    # Plot magnets
    for i in range(len(keywords)):
        if 'QUADRUPOLE' in keywords[i]:
            draw_quad(s[i]-L[i],L[i],names[i],'#CB6015')
        elif 'RBEND' in keywords[i]:
            draw_dipole(s[i]-L[i],L[i],names[i],'#004C97')
        elif 'KICKER' in keywords[i]:
            draw_corrector(s[i]-L[i],L[i],names[i],'#4C8C2B')
        elif names[i].startswith("PM"):
            draw_instrumentation(s[i]-L[i],L[i],names[i],'#EAAA00')

    plt.subplot(312)
    plt.plot(s,betx,label=r'$\beta_x$')
    plt.plot(s,bety,label=r'$\beta_y$')
    plt.xlabel('s [m]')
    plt.xlim(0,s[-1])
    plt.ylabel(r'$\beta$ [m]',fontsize=14)
    plt.legend(loc='upper left')
    plt.gca().set_axis_bgcolor("#F8F8F8")
    plt.grid()

    plt.subplot(313)
    plt.plot(s,dx,label=r'$D_x [m]$')
    plt.plot(s,dy,label=r'$D_y [m]$')
    #plt.plot(s_Hbpms,eta_H,'bo')
    #plt.plot(s_Vbpms,eta_V,'go')
    plt.xlabel('s [m]')
    plt.xlim(0,s[-1])
    plt.ylabel(r'$D [m]$',fontsize=14)
    plt.legend(loc='upper left')
    plt.gca().set_axis_bgcolor("#F8F8F8")
    plt.grid()

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('plot_twiss.png')

def plotSize(s,betx,bety,dx,dy,names,L):

    # Compute Lorentz factors
    E = 120.0E9 # Beam total energy in eV
    gamma = E/938E6
    beta = np.sqrt(1-(1/gamma)**2)

    # Normalized emittances, 95%
    ex_N = 18E-6 # pi*m*r
    ey_N = 18E-6 # pi*m*r

    # Geometric emittances, 95%
    ex = ex_N/(beta*gamma)
    ey = ey_N/(beta*gamma)

    # Momentum spread
    dpp = 2*1.0649888E-3 # Given in RMS, so multiply by 2 to get 95%

    xsig = (np.asarray(betx)*ex + np.asarray(dx)**2*dpp**2)**0.5*1000.0
    ysig = (np.asarray(bety)*ey + np.asarray(dy)**2*dpp**2)**0.5*1000.0

    s_profiles = []
    for i in range(len(names)):
        if names[i].startswith('PM'):
                if (names[i] != 'PM115') and (names[i] != 'PM112'):
                    s_profiles.append(float(s[i]))
    Hsigmas = []
    Vsigmas = []
    with open('profiles.csv','r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            Hsigmas.append(float(line.split(',')[1]))
            Vsigmas.append(float(line.split(',')[2]))

    # Convert beam widths to two-sigma, 95%
    Hsigmas = np.asarray(Hsigmas)*2.0
    Vsigmas = np.asarray(Vsigmas)*2.0

    plt.figure(figsize=(16,8))
    plt.suptitle('NuMI Beamline',fontsize=22)

    plt.subplot(311)
    plt.plot(s,[0]*len(s),'k-')
    plt.xlim(0,s[-1])
    plt.ylim(-1,1)
    plt.gca().set_xlim(left=0)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.gca().set_axis_bgcolor('#f1f1f1')
    # Plot magnets
    for i in range(len(keywords)):
        if 'QUADRUPOLE' in keywords[i]:
            draw_quad(s[i]-L[i],L[i],names[i],'#CB6015')
        elif 'RBEND' in keywords[i]:
            draw_dipole(s[i]-L[i],L[i],names[i],'#004C97')
        elif 'KICKER' in keywords[i]:
            draw_corrector(s[i]-L[i],L[i],names[i],'#4C8C2B')
        elif names[i].startswith("PM"):
            draw_instrumentation(s[i]-L[i],L[i],names[i],'#EAAA00')

    plt.subplot(312)
    plt.plot(s,xsig,label=r'$\sigma_x$',color='blue')
    plt.plot(s_profiles,Hsigmas,label=r'PM $\sigma_x$',marker='o',linestyle='None',color='red')
    plt.xlabel('s [m]')
    plt.xlim(0,s[-1])
    plt.ylabel(r'2$\sigma_x$ [mm]',fontsize=14)
    plt.legend(loc='upper center', bbox_to_anchor=(0.25, 1.10),
          ncol=3, fancybox=True, shadow=True)
    plt.gca().set_axis_bgcolor("#F8F8F8")
    plt.grid()

    plt.subplot(313)
    plt.plot(s,ysig,label=r'$\sigma_y$',color='green')
    plt.plot(s_profiles,Vsigmas,label=r'PM $\sigma_y$',marker='o',linestyle='None',color='red')
    plt.xlabel('s [m]')
    plt.xlim(0,s[-1])
    plt.ylabel(r'2$\sigma_y$ [mm]',fontsize=14)
    plt.legend(loc='upper center', bbox_to_anchor=(0.25, 1.10),
          ncol=3, fancybox=True, shadow=True)
    plt.gca().set_axis_bgcolor("#F8F8F8")
    plt.grid()

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig('plot_size.png')

names,keywords,s,betx,alphax,mux,bety,alphay,muy,dx,dy,L = parsetwiss('twiss.out')
plot(s,betx,bety,dx,dy,names,L)
plotSize(s,betx,bety,dx,dy,names,L)

with open('stations.txt','w') as f:
    for i in range(len(s)):
        f.write('%s \t %f \t %f \n'%(names[i], s[i], L[i]))

NAME = []
X = []
Y = []
Z = []
with open('survey.txt', 'r') as f:
    lines = f.readlines()
for line in lines[8:]:
    NAME.append(line.split()[0])
    X.append(float(line.split()[4]))
    Y.append(float(line.split()[5]))
    Z.append(float(line.split()[6]))

with open('parsed_survey.txt','w') as f:
    for i in range(len(NAME)):
        f.write('%s \t %f \t %f \t %f \n'%(NAME[i], X[i], Y[i], Z[i]))
