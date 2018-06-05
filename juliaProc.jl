using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = [12]
using PyPlot
@pyimport matplotlib.colors as col

unshift!(PyVector(pyimport("sys")["path"]), "")

@pyimport readFileList as r
@pyimport readData as rd

include("cfad.jl")
cfad=zeros(40,40)
fList=r.readFileNames()
for if1 in range(1,16)
    println(fList[if1])
    fname=fList[if1]
    i=1
    (zFactKu,node,flag,precipType,sfcType,zeroDeg,binClut,
     lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,
     binBBPeak,piaHB,sfcRain,binSf,reliabF,locZAngle,s0,eps,dm,n)=
    rd.readKuR(fname,(i-1)*300,i*300,12,37)
    
    
    
    
    nchunks= trunc(Int,n/300)
    for i in range(1,nchunks)
        (zFactKu,node,flag,precipType,sfcType,zeroDeg,binClut,
         lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,
         binBBPeak,piaHB,sfcRain,binSf,reliabF,locZAngle,s0,eps,dm,n)=
        rd.readKuR(fname,(i-1)*300,i*300,12,37)
        iret=cfad_Fill0(zFactKu,hZero,zeroDeg,precipType,binClut,binBBPeak,cfad)
    end

end
x=linspace(12,51,40)
y=linspace(1,40,40)*0.25-5
y=y[end:-1:1]+0.125
PyPlot.pcolormesh(x,y,cfad',norm=col.LogNorm(),cmap="jet")
PyPlot.ylim(-1,4)
PyPlot.title("Ku-band reflectivity CFAD")
PyPlot.xlabel("Reflectivity")
PyPlot.ylabel("Height relative to zero isotherm (km)")
c=PyPlot.colorbar()
c[:set_label]("Counts")
PyPlot.savefig("stratiform_noBB_CFAD.png")
