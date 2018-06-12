using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = [12]
using PyPlot
@pyimport matplotlib.colors as col

unshift!(PyVector(pyimport("sys")["path"]), "")

#@pyimport testSDSU as sdsu
#@pyimport testSDSU as sdsu
#sdsu.init()
imd=2
ccall((:setparams_, "SDSU"), Void, () )    
ccall((:simradarlut_, "SDSU"), Void, (Ptr{Int32},), &imd)

function getKExtWint()
    kext=[0.08216629, 0.07987797, 0.07557384, 0.07255906, 0.06967789,
          0.06681726, 0.06410733, 0.0579389 , 0.0535338 , 0.05497366,
          0.04758603, 0.04699508, 0.04475723, 0.04301621, 0.0412752 ,
          0.03960225, 0.03795199, 0.03521832, 0.03287926, 0.03057081,
          0.02841735, 0.02642137, 0.02526997, 0.0246093 , 0.02303219,
          0.02217608, 0.0213438 , 0.02051153, 0.01967925, 0.01884698,
          0.01811531, 0.01741973, 0.01671269, 0.01577828, 0.01484388,
          0.01403029, 0.01323571, 0.0125283 , 0.011831  , 0.01113371,
          0.01050953, 0.01017607, 0.00990444, 0.00963282, 0.00926262,
          0.00873582, 0.00813506, 0.00680915, 0.00657739, 0.00638145,
          0.00615319, 0.00587261, 0.00559203, 0.00531145, 0.00503086,
          0.00458927, 0.00440378, 0.00420007, 0.00394702, 0.00380493,
          0.00365312, 0.00348672, 0.00336368, 0.00330825, 0.00325281,
          0.00319737, 0.00306326, 0.00292377, 0.00278428, 0.00264478,
          0.0025104 , 0.00241081, 0.00231121, 0.00223822, 0.00216396,
          0.00209095, 0.00201794, 0.00193494, 0.00183698, 0.00174985]
    return kext
end
function getcloud(freq,cldw)
    t=Float32(273.)
    z_clw=Array{Float32}(1)
    cldwi=Float32(cldw)
    freqi=Float32(freq)
    #print(z_clw)
    ccall((:gcloud_,"SDSU"),Void,(Ptr{Float32},Ptr{Float32},
                                  Ptr{Float32},Ptr{Float32}),
          &freqi,&t,&cldwi,z_clw)
    z_clw1=z_clw[1]
    #print(z_clw1)
    #print(x)
    #SUBROUTINE GCLOUD(FREQY,T,CLW,Z_CLW)
    return z_clw1
end

function getScatt(imelt,wrain,wsnow,wgrpl,whail,icldw,temp,fmelt)
    wcldw=0.
    nrdfreq=6
    kexttot=Array{Float32}(nrdfreq)
    asymtot=Array{Float32}(nrdfreq)
    salbtot=Array{Float32}(nrdfreq)
    kexttotm=Array{Float32}(nrdfreq)
    asymtotm=Array{Float32}(nrdfreq)
    salbtotm=Array{Float32}(nrdfreq)
    asymtotm=Array{Float32}(nrdfreq)
    salbtotm=Array{Float32}(nrdfreq)
    dBZ=Array{Float32}(nrdfreq)
    dBZm=Array{Float32}(nrdfreq)
  
    ccall((:fromlktables2_,"SDSU"),Void, (Ptr{Int32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                       Ptr{Float32},Ptr{Float32},Ptr{Int32}),
           &imelt,&wrain,&wsnow,&wgrpl,&whail,&wcldw,&icldw,
           &temp, kexttot,salbtot,asymtot,kexttotm,
          salbtotm, asymtotm, dBZ, dBZm,&fmelt,&nrdfreq)
    return dBZ,kexttot,dBZm,kexttotm
          
end


function gauss_newton(zobs,ifreq)
    imelt=0
    wrain=0.1
    wsnow=0.
    wgrpl=0.
    whail=0.
    icldw=0.
    temp=274.
    fmelt=0.
    dn=0.4
    eps=1e3
    for it=1:5
        dBZ,kexttot,dBZm,kexttotm=getScatt(imelt,wrain/dn,wsnow,wgrpl,whail,
                                           icldw,temp,fmelt)
        dBZ=dBZ+log10(dn)*10.
        kexttot=kexttot*dn
        wrain1=1.1*wrain
        dBZ1,kexttot1,dBZm1,kexttotm1=getScatt(imelt,wrain1/dn,wsnow,wgrpl,
                                               whail,icldw,temp,fmelt)
        dBZ1=dBZ1+log10(dn)*10.
        kexttot1=kexttot1*dn
        
        eps=dBZ[ifreq]-zobs
        gradZ=(dBZ1[ifreq]-dBZ[ifreq])/(0.1*wrain)
        dwrain=eps*gradZ/(gradZ*gradZ+0.001)
        wrain=wrain-dwrain
        if wrain<0.01
            wrain=0.01
        end
    end
  
    return wrain, eps
end
