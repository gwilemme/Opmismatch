
using NBInclude
## IVER = 5
@nbinclude("startup.ipynb") ## choose IVER before, or IVER = 3 by default

println("---*** Begin STEP1 on version $IVER ***---")
let 
    dist = sum(abs2,(PAR.emp.mom2-Diagonal(PAR.emp.mom2)) .* PAR.wgt[2])  /(22*22) + 
        sum(abs2,PAR.emp.mom3 .*PAR.wgt[3])/(22*21)
    println("The distance cannot go below $dist with segmented equilibrium.")
end

function compute_step1(ES)
    
    comp1 = 60*30 #60*30 30min
    comp2 = 60*3 #60*3 3min
    
    ## calibrate a segmented equilibrium
    (Asegm, Msegm) = calibrate_segm_equilibrium(PAR, RTax, Para_es(ES,M_init), N,
                    toestim=[true,true,true],verbose=false, maxti=comp1 )

    ## fix a benchmark equilibrium with arbitrary alpha and xi
    (A0,M0) = let
            alpha_cons = 0.9
            xi = 1.1
            vec= vcat(diag(Msegm.y),fill(alpha_cons,N.J),Msegm.es, xi, Asegm.theta, nu2nuc(Msegm.nu))
            compute_equilibrium_theta(vec,RTax,Msegm,N,Xtol=1e-2)
    end
    M0.k[:] = compute_k(A0,RTax,M0,N)            
    println("Distance benchmark alpha=0.95 and xi=1: ", distance_mom(PAR, A0,RTax,M0,N))

    
    ## find good priors for alpha and xi
    (A1, M1) = calibrate_equilibrium_step1(PAR,RTax,M0,N, 
            A0=A0,verbose = false,  Xtol=1e-2, maxti= comp2, algo=:LN_SBPLX) 
    M1.k[:] = compute_k(A1,RTax,M1,N)   
    
    dist = distance_mom(PAR, A1,RTax,M1,N)
    return( dist, A1, M1 )
end

rangeES = range(0.4,0.6,length=10)
distvec = fill(1000.,length(rangeES))
index = 0
for ES in rangeES
    global index += 1
    
    esr = round(ES,digits=2)
    println("- Loop $index, for ES=$esr")
    dist = compute_step1(ES)[1]
    
    distvec[index] = dist
end    
;

ES = let 
    dist_noNaN = distvec[.!isnan.(distvec)] ## exclude NaN
    rangeES_noNaN = rangeES[.!isnan.(distvec)]
    rangeES_noNaN[findmin(dist_noNaN)[2]]
end

dist, A1, M1 = compute_step1(ES)

## save
save_struc(A1, M1, joinpath(dir_calib,"step1_v$IVER.csv"), "v$IVER", N) 

esr = round(M1.es,digits = 2)
xir = round(M1.xi, digits=2)
alphar = round((M1.y[1,1]-M1.y[2,1])/(M1.y[1,1]-M1.h),digits=2)
println("---*** End STEP1 on version $IVER with ES=$esr, XI=$xir, alph=$alphar ***---")
