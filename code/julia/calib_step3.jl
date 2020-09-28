
using NBInclude
## IVER = 5
@nbinclude("startup.ipynb") ## choose IVER before, or IVER = 3 by default

if "COMP" in keys(ENV)
    COMP = parse(Int,ENV["COMP"])
else
    COMP = 1/60 
end

println("---*** Begin STEP3 on version $IVER, for maxtime=$(COMP)h***---")

A0, M0 = open_struc(joinpath(dir_calib,"step2_v$IVER.csv"),N)

println("Distance before: ",distance_mom(PAR, A0,RTax,M0,N))

@time (A1, M1) = calibrate_equilibrium(PAR,RTax,M0,N, 
        toestim=[false,false,true,true,true,false],
        A0=A0,verbose = false,  Xtol=1e-2, maxti= 60*60*COMP, algo=:LN_SBPLX) 
   
M1.k[:] = compute_k(A1,RTax,M1,N)

println("Distance after: ",distance_mom(PAR, A1,RTax,M1,N))

## save
save_struc(A1, M1, joinpath(dir_calib,"step3_v$IVER.csv"), "v$IVER", N)

esr = round(M1.es, digits = 2)
xir = round(M1.xi, digits=2)
println("---*** End STEP3 on version $IVER with ES=$esr, XI=$xir***---") 
