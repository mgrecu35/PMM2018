using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = [12]
using PyPlot
unshift!(PyVector(pyimport("sys")["path"]), "")
@pyimport readPickle as r
@pyimport numpy as np
(zsL,t1sL,citL,fname)=r.readPickle()
include("SDSU_Wrap.jl")
iprof=1000
zw=zsL[iprof,3,:];
piaW=0.;
zwC=Array{Float32}(80)*0.
zwC2=Array{Float32}(80)*0.
zkuP=Array{Float32}(80)*0.
wrainL=Array{Float64,1}(0)
kextL=Array{Float64,1}(0)
dbzL=Array{Float64,1}(0)
for k=1:100
    imelt=0;
    wrain=0.01*k;
    wsnow=0.;
    wgrpl=0.;
    icldw=0.;
    temp=273.;
    fmelt=0.;
    whail=0.;
    fmelt=0.;
    dBZ,kexttot,dBZm,kexttotm=getScatt(imelt,wrain,wsnow,wgrpl,whail,icldw,temp,fmelt);
    push!(wrainL,wrain)
    push!(dbzL,dBZ[4])
    push!(kextL,kexttot[4])
end

kextL2=Array{Float64,1}(0)
dbzL2=Array{Float64,1}(0)
for k=1:100
    imelt=0;
    wsnow=0.01*k;
    wrain=0.;
    icldw=0.;
    temp=273.;
    fmelt=0.;
    whail=0.;
    wgrpl=0.;
    fmelt=0.;
    dBZ,kexttot,dBZm,kexttotm=getScatt(imelt,wrain,wsnow,wgrpl,whail,icldw,temp,fmelt)
    #push!(wrainL,wrain)
    push!(dbzL2,dBZ[4])
    push!(kextL2,kexttot[4])
end

res=np.polyfit(dbzL,log10.(1.0*4.343*kextL),1)
sRes=np.polyfit(dbzL2,log10.(1.0*4.343*kextL2),1)

z_clw=getcloud(94.,0.15)
print(z_clw)

dr=0.125;
kextW=getKExtWint()
zwC2D=zeros(2001,80)
piaWL=Array{Float64}(0)
qKu=Array{Float64}(0)
qW=Array{Float64}(0)
for iprof=500:2500
    piaW=0
    zetaS=0.
    beta=sRes[1]*10.
    betav=zeros(80)
    q=0.2*log(10.)*beta
    alpha=10.^sRes[2]
    zeta=zeros(80)
    zw=zsL[iprof,3,:]+0.
    
    for k=80:-1:1
        piaW+=kextW[k]*4.343*0.125
        zw[k]+=piaW
        piaW+=kextW[k]*4.343*0.125
    end
    #print(piaW)
    piaWC=0
    for k=15:-1:1
        if k>=8+3
            piaWC+=z_clw*4.343*0.15
            zw[k]+=piaWC
            piaWC+=z_clw*4.343*0.125
        else
            zw[k]+=piaWC
        end
    end
    for k=80:-1:12
        if (zw[k]==zw[k]) & (zw[k]>0) 
            zeta[k]=zetaS+q*alpha*10.^(0.1*zw[k]*beta)*dr
            zetaS=zeta[k]
            betav[k]=beta
        end
    end
    
    for k=11:-1:8
        if (zw[k]==zw[k]) & (zw[k]>0)
            f=(k-8)/4.
            beta=f*sRes[1]*10.+(1-f)*res[1]*10.
            q=0.2*log(10.)*beta
            alpha=f*10.^sRes[2]+(1-f)*10.^res[2]
            zeta[k]=zetaS+q*alpha*10.^(0.1*zw[k]*beta)*dr
            zetaS=zeta[k]
            betav[k]=beta
        end
    end

    for k=7:-1:4
        if (zw[k]==zw[k]) & (zw[k]>0)
            beta=res[1]*10.
            q=0.2*log(10.)*beta
            alpha=10.^res[2]
            zeta[k]=zetaS+q*alpha*10.^(0.1*zw[k]*beta)*dr
            zetaS=zeta[k]
            betav[k]=beta
        end
    end
    zwc=zeros(80)
    if(zeta[4]>0.98) 
        eps=0.98/zeta[4]
    else
        eps=1.
    end
    for k=80:-1:4
        if (zw[k]==zw[k]) & (zw[k]>0) 
            zwc[k]=zw[k]-log10(1-eps*zeta[k])*10/betav[k]
        end
    end
    push!(piaWL,-log10(1-eps*zeta[4])*10/betav[4])
    print(iprof)
    zwC2D[iprof-500+1,:]=zwc
    q1,e1=gauss_newton(zwc[6],4)
    q2,e2=gauss_newton(zsL[iprof-35,1,6],2)
    push!(qKu,q2)
    push!(qW,q1)
    #print(e1)
    #print(e2)
    #println()
end
#z2=np.ma[:array](zwC2D,mask=zwC2D.<1)
r.plotret(t1sL[500:2500],range(1,80)*0.125,zwC2D,30.,
          "Corrected W-band Reflectivity", "zCWband.png")
r.plotret(t1sL[500:2500],range(1,80)*0.125,zsL[500-35:2500-35,1,:],40.,
          "Observed Ku-band Reflectivity", "zKuband.png")
r.plotret(t1sL[500:2500],range(1,80)*0.125,zsL[500:2500,3,:],30.,
          "Observed W-band Reflectivity", "zWband.png")
#PyPlot.pcolormesh(t1sL[1000:2000],range(1,80)*0.25,z2',vmin=0,vmax=40)
PyPlot.figure()
PyPlot.plot(t1sL[500:2500],39.5-zsL[500:2500,3,1])
PyPlot.plot(t1sL[500:2500],piaWL)
PyPlot.legend(["SRT PIA","HB PIA"])
PyPlot.xlabel("Time")
PyPlot.ylabel("dB")
PyPlot.savefig("piaW.png")

r.plotWC(qKu,qW)


print(np.corrcoef(37-zsL[500:2500,3,1],piaWL))
print(np.corrcoef(zwC2D[:,4],zsL[500-35:2500-35,1,4]))
print(np.corrcoef(zsL[500:2500,3,4],zsL[500-35:2500-35,1,4]))
