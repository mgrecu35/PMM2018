function cfad_Fill(zku,hzero,zeroDeg,precipType,binClut,binBBPeak,cfad)
    (nscan,nray,nh)=size(zku)
    for j in range(1,nray)
        for i in range(1,nscan)
            if trunc(Int,precipType[i,j]/1e7)==1
                #print((i,j))
                nc=binClut[i,j]
                #println(nc)
                for k in range(1,nc)
                    if (zku[i,j,k]>12) & (binBBPeak[i,j]>0) & (hzero[i,j]<2000)
                        #println(zku[i,j,k])
                        ih=trunc(Int,(k-binBBPeak[i,j])*0.25)+21
                        iz=trunc(Int,(zku[i,j,k]-12)/1.)+1
                        #println((ih,iz))
                        if ((iz>0) & (iz<=40) & (ih>=1) & (ih<=40))
                            cfad[iz,ih]=cfad[iz,ih]+1
                        end
                    end
                end
            end
        end
    end
    return 1
end

function cfad_Fill0(zku,hzero,zeroDeg,precipType,binClut,binBBPeak,cfad)
    (nscan,nray,nh)=size(zku)
    for j in range(1,nray)
        for i in range(1,nscan)
            if trunc(Int,precipType[i,j]/1e7)==1
                #print((i,j))
                nc=binClut[i,j]
                #println(nc)
                for k in range(1,nc)
                    if (zku[i,j,k]>12) & (zeroDeg[i,j]>0) & (hzero[i,j]<2000)
                        #println(zku[i,j,k])
                        ih=trunc(Int,(k-zeroDeg[i,j])*0.25)+21
                        iz=trunc(Int,(zku[i,j,k]-12)/1.)+1
                        #println((ih,iz))
                        if ((iz>0) & (iz<=40) & (ih>=1) & (ih<=40))
                            cfad[iz,ih]=cfad[iz,ih]+1
                        end
                    end
                end
            end
        end
    end
    return 1
end
