using NBInclude
@nbinclude("startup.ipynb")
timing = 60*30 ## time per computation of equilibrium/optimum 
println("Code takes ", timing*2*3/60/60, " hours." )

df_stats = DataFrame(CSV.read(joinpath(dir_mom,"stats.csv")))
println("The unemployment rate is at ", 
    round(sum((df_stats.emp1 .== 0) .* df_stats.N) / sum(df_stats.N),sigdigits=2)) 
println("The share of unemployed workers with known occupation or armed forces is at ", 
    round(df_stats.N[(df_stats.emp1 .== 0) .& (df_stats.occ1 .==-1)][1] / sum((df_stats.emp1 .== 0) .* df_stats.N),digits=2)) 
println("The av. wage in pop is ", 
    round(weightedmean(df_stats.w[df_stats.emp1 .== 1],df_stats.Nw[df_stats.emp1 .== 1]),digits=1))


df_rates = DataFrame(CSV.read(joinpath(dir_mom,"rates.csv")))
(shUU,shUE,shEU,shEEnew,shEEsame) =share_transit(df_rates,N)
(nUU,nUE,nEU,nEEnew,nEEsame) = n_transit(df_rates,N)

ntot = sum(nUU)+sum(nUE)+sum(nEU)+sum(nEEnew)+sum(nEEsame)
nU = sum(nUU)+sum(nUE)
println("The number of transitions observed (for 23 months) is ", Int(ntot))
println("The share of transitions from U is ", round(nU/ntot,sigdigits=2))
println("The share of unemployed workers getting employed is ", round((sum(nUE))/nU,sigdigits=2))
println("The share of unemployed workers now employed in a new occ. is ", 
            round((sum(nUE)-sum(diag(nUE)))/sum(nUE),sigdigits=2))
println("The share of employed workers becoming unemployed is ", round(sum(nEU)/(ntot-nU),sigdigits=2))
println("The share of employed workers changing job is ", round(sum(nEEnew)/(ntot-nU),sigdigits=2))
println("The share of employed workers changing occupation is ", round((sum(nEEnew)-sum(diag(nEEnew)))/(ntot-nU),sigdigits=2))
println("The share of E2E with occupation change employed is ", round((sum(nEEnew)-sum(diag(nEEnew)))/sum(nEEnew),sigdigits=2))

## stats on moments
mom1 = emp.mom1
mom2=emp.mom2
mom3=emp.mom3
mom4 = emp.mom4
let
    println("--- Description of the moments.")
    mi = round(minimum(mom4),sigdigits=3)
    omi = collect(1:length(mom4))[mom4.==minimum(mom4)][1] 
    ma = round(maximum(mom4),sigdigits=3)
    oma = collect(1:length(mom4))[mom4.==maximum(mom4)][1] 
    println("The most frequent occupation is ", oma, " with ", ma, ". ", 
        "The least frequent occupation is ", omi, " with ", mi, ".")
    me = round(mean(mom1),digits=1)
    mi = round(minimum(mom1),digits=1)
    omi = collect(1:length(mom1))[mom1.==minimum(mom1)][1] 
    ma = round(maximum(mom1),digits=1)
    oma = collect(1:length(mom1))[mom1.==maximum(mom1)][1] 
    println("The occupational wage has mean ", me , 
        ". The min is at ", mi , " for occ. ", omi ,
        ". The max is at ", ma , " for occ. ", oma , ".")
    me = round(mean(jdr),sigdigits=3)
    mi = round(minimum(jdr),sigdigits=3)
    omi = collect(1:length(jdr))[jdr.==minimum(jdr)][1] 
    ma = round(maximum(jdr),sigdigits=3)
    oma = collect(1:length(jdr))[jdr.==maximum(jdr)][1] 
    println("The job separation rate has mean ", me , 
        ". The min is at ", mi , " for occ. ", omi ,
        ". The max is at ", ma , " for occ. ", oma , ".")
    jfr = [sum(mom2[j0,:]) for j0 in 1:N.J]
    me = round(mean(jfr),sigdigits=3)
    mi = round(minimum(jfr),sigdigits=3)
    omi = collect(1:length(jfr))[jfr.==minimum(jfr)][1] 
    ma = round(maximum(jfr),sigdigits=3)
    oma = collect(1:length(jfr))[jfr.==maximum(jfr)][1] 
    println("The job-finding rate has mean ", me , 
        ". The min is at ", mi , " for occ. ", omi ,
        ". The max is at ", ma , " for occ. ", oma , ".")
    osr = [sum(mom3[j0,:]) for j0 in 1:N.J]
    me = round(mean(osr),sigdigits=3)
    mi = round(minimum(osr),sigdigits=3)
    omi = collect(1:length(osr))[osr.==minimum(osr)][1] 
    ma = round(maximum(osr),sigdigits=3)
    oma = collect(1:length(osr))[osr.==maximum(osr)][1] 
    println("The occupation-switching rate has mean ", me , 
        ". The min is at ", mi , " for occ. ", omi ,
        ". The max is at ", ma , " for occ. ", oma , ".")
end
;

table = zeros(27,3)

for i in 1:3

    (Acal,Mcal) = open_struc(joinpath(dir_calib,"step4_v$i.csv"),N)
    
    Mcal.eta
    ## table
    table[1,i] = Mcal.eta
    table[2,i] = Mcal.beta
    index = 3
    
    
    ## calibrated parameters
    table[index,i] = round(mean(diag(Mcal.y)),digits=1)
    table[index+1,i] = round(100* mean( vcat( 1-Mcal.y[2,1]/Mcal.y[1,1], [1-Mcal.y[1,i]/Mcal.y[i,i] for i in 2:N.J]  ) ))
    table[index+2,i] = round(Mcal.es,digits=3)
    table[index+3,i] = round(Mcal.xi,digits=3)
    index += 4
    
    
    ## fit
    theo = theo_moments(Acal,RTax,Mcal,N)
    table[index,i] = round(sqrt(0.25*distance_mom(PAR.emp,theo,PAR.wgt,N)),digits=2)
    index += 1
    
    ## moments 4 
    dist = (theo.mom4 - PAR.emp.mom4) * PAR.wgt[4]
    table[index,i] = round(sqrt(sum(abs2,dist ) / N.J), digits = 2) ## contribution
    table[index+1,i] = round(mean(dist),sigdigits = 3)
    table[index+2,i] = round(std(dist), sigdigits = 2)
    table[index+3,i] = round(minimum(dist), digits = 2) ## min
    table[index+4,i] = round(maximum(dist),digits = 2) ## max
    index += 5
    
    ## moments 1 
    dist = (theo.mom1 - PAR.emp.mom1) * PAR.wgt[1]
    table[index,i] = round(sqrt(sum(abs2,dist ) / N.J), digits = 2) ## contribution
    table[index+1,i] = round(mean(dist),sigdigits = 2)
    table[index+2,i] = round(std(dist), sigdigits = 2)
    table[index+3,i] = round(minimum(dist), digits = 2) ## min
    table[index+4,i] = round(maximum(dist),digits = 2) ## max
    index += 5
    
    ## moments 2 
    dist = (theo.mom2 - PAR.emp.mom2) * PAR.wgt[2]
    table[index,i] = round(sqrt(sum(abs2,dist ) / (N.J*N.J)), digits = 2) ## contribution
    table[index+1,i] = round(mean(dist),sigdigits = 2)
    table[index+2,i] = round(std(dist), sigdigits = 2)
    table[index+3,i] = round(minimum(dist), digits = 2) ## min
    table[index+4,i] = round(maximum(dist),digits = 2) ## max
    index += 5
    
    ## moments 3 
    dist = (theo.mom3 - PAR.emp.mom3) * PAR.wgt[3]
    table[index,i] = round(sqrt(sum(abs2,dist ) / (N.J*(N.J-1))), digits = 2) ## contribution
    table[index+1,i] = round(mean(dist),sigdigits = 2)
    table[index+2,i] = round(std(dist), sigdigits = 2)
    table[index+3,i] = round(minimum(dist), digits = 2) ## min
    table[index+4,i] = round(maximum(dist),digits = 2) ## max

end

function printlineLatex(firstcol, index)
    print("$firstcol")
    for i in 1:3
        print(" & ",table[index,i])
    end
    println("\\\\")
end

println("Fixed parameters &&&\\\\")
printlineLatex("\\qquad\\qquad  matching elasticity \$\\eta\$", 1)
printlineLatex("\\qquad\\qquad  bargaining power \$\\beta\$", 2)
index = 3

## calibrated par
println("Calibrated parameters &&&\\\\")
printlineLatex("\\qquad\\qquad  mean perfect-match productivity", index)
printlineLatex("\\qquad\\qquad  mean mismatching penalty (\\%)", index+1)
printlineLatex("\\qquad\\qquad  search elasticity \$\\epsilon\$", index+2)
printlineLatex("\\qquad\\qquad  search parameter \$\\xi\$", index+3)
index += 4

println("\\hline\\hline")
printlineLatex("Distance to moments", index)
index += 1

## mom1
println("Occupational employment \$\\mathcal{M}_1\$ &&&\\\\")
printlineLatex("\\qquad\\qquad distance", index)
printlineLatex("\\qquad\\qquad mean", index+1)
printlineLatex("\\qquad\\qquad standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)
index += 5

## mom2
println("Occupational wages \$\\mathcal{M}_2\$ &&&\\\\")
printlineLatex("\\qquad\\qquad distance", index)
printlineLatex("\\qquad\\qquad mean", index+1)
printlineLatex("\\qquad\\qquad standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)
index += 5

## mom3
println("Occupational job-finding rates \$\\mathcal{M}_3\$ &&&\\\\")
printlineLatex("\\qquad\\qquad distance", index)
printlineLatex("\\qquad\\qquad mean", index+1)
printlineLatex("\\qquad\\qquad standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)
index += 5

## mom4
println("Occupation-switching rates \$\\mathcal{M}_4\$ &&&\\\\")
printlineLatex("\\qquad\\qquad distance", index)
printlineLatex("\\qquad\\qquad mean", index+1)
printlineLatex("\\qquad\\qquad standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)

table = zeros(18,3)

for i in 1:3
    println(i)
    
    (Acal,Mcal) = open_struc(joinpath(dir_calib,"step4_v$i.csv"),N)
    
    ## table
    table[1,i] = Mcal.eta
    table[2,i] = Mcal.beta
    index = 3
    
    ## op vs eq
    A0 = compute_equilibrium(Notax,Mcal,N,verbose=false,theta0=Acal.theta,maxti=timing,
        valf=0.001,Xtol=1e-3)
    opA= compute_optimum(Mcal,N,verbose=false,theta0=Acal.theta,maxti=timing,
        valf=0.001,Xtol=1e-3)
    oT1,oT2 = optax(opA,Mcal,N)
    println("Compare welfare with different formula:")
    println("decentralised with flat tax:", welfare_dec(Acal,RTax,Mcal,N), " vs ",
                    welfare2(Acal,Mcal,N))
    println("decentralised laissez-faire:", welfare_dec(A0,Notax,Mcal,N), " vs ",
                    welfare2(A0,Mcal,N))
    println("oT1 vs oT2: ", welfare_dec(opA,oT1,Mcal,N)," vs " , welfare_dec(opA,oT2,Mcal,N) )
    println("welfare opt: ", welfare_opt(opA,Mcal,N) )
    println("welfare2 opt: ", welfare2(opA,Mcal,N) )
    oT=oT2 ## choose oT2
    #opA = compute_equilibrium(oT,Mcal,N,verbose=false,theta0=Acal.theta,maxti=timing,
    #    valf=0.001,Xtol=1e-3)
    #println("welfare opt: ", welfare_dec(opA,oT,Mcal,N) )
    
    ## efficiency loss
    table[index,i] = round((1 - welfare_dec(Acal,RTax,Mcal,N) / welfare_opt(opA,Mcal,N))*100,digits=1)
    table[index+1,i] = round((1 - welfare_dec(A0,Notax,Mcal,N) / welfare_opt(opA,Mcal,N))*100,digits=1)
    index += 2
    
    ## vacancies
    vAcal = vacancies(Acal,Mcal,N)
    vA0 = vacancies(A0,Mcal,N)
    vopA = vacancies(opA,Mcal,N)
    table[index,i] = round((sum(vAcal)/sum(vopA)-1)*100, sigdigits = 3)
    table[index+1,i] = round((sum(vA0)/sum(vopA)-1)*100, sigdigits = 3)
    println((1 - mean(Acal.theta ./ opA.theta)) *100)
    println((1 - mean(A0.theta ./ opA.theta)) *100)
    index += 2
    
    ## unemployment
    table[index,i] = round(sum(Acal.u)*100,digits=1)
    table[index+1,i] = round(sum(A0.u)*100,digits=1)
    table[index+2,i] = round(sum(opA.u)*100,digits=1)
    index += 3
    
    ## mismatch
    table[index,i] = round((1-sum([Acal.n[j,j] for j in 1:N.J])/sum(Acal.n))*100,digits=1)
    table[index+1,i] = round((1-sum([A0.n[j,j] for j in 1:N.J])/sum(A0.n))*100,digits=1)
    table[index+2,i] = round((1-sum([opA.n[j,j] for j in 1:N.J])/sum(opA.n))*100,digits=1)
    
    # Mm ratios
    MmJ, MmI = meanmin(Acal,RTax,Mcal,N)
    table[index+3,i] = round(mean(MmI),digits=1)
    MmJ, MmI = meanmin(Acal,Notax,Mcal,N)
    table[index+4,i] = round(mean(MmI),digits=1)
    MmJ, MmI = meanmin(opA,oT,Mcal,N)
    table[index+5,i] = round(mean(MmI),digits=1)
    
    ## redistrib
    table[index+6,i] = round(ty_stat2(opA,oT,Mcal,N),digits=2)
    table[index+7,i] = round(ty_stat(opA,oT,Mcal,N),digits=2)
    thetaRD = [oT.f1[j0] * (Mcal.r+Mcal.delta[j1]) for j0 in 1:N.J, j1 in 1:N.J] 
    ## Theta^{j_1}_{j0} *(r+\delta_{j_1})
    table[index+8,i] = round(mean(thetaRD),sigdigits=3)
end

function printlineLatex(firstcol, index)
    print("$firstcol")
    for i in 1:3
        print(" & ",table[index,i])
    end
    println("\\\\")
end


println("Fixed parameters &&&\\\\")
printlineLatex("\\qquad\\qquad  matching elasticity \$\\eta\$", 1)
printlineLatex("\\qquad\\qquad  bargaining power \$\\beta\$", 2)
println("\\hline\\hline")
index = 3

## efficiency
println("Efficiency loss (\\%)&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
index += 2

## vacancies
println("Vacancies compared to optimum (\\%)&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
index += 2

## unemployment
println("Unemployment rate (\\%) &&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
printlineLatex("\\qquad\\qquad  optimum", index+2)
index += 3

## mismatch
println("Share of mismatches (\\%) &&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
printlineLatex("\\qquad\\qquad  optimum", index+2)
index += 3

## Mm
println("Mean-min ratio &&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
printlineLatex("\\qquad\\qquad  optimum", index+2)

println("Taxes-production coeff. &&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index+3)
printlineLatex("\\qquad\\qquad  laissez-faire", index+4)
printlineLatex("Tax on job-to-job transition", index+5)


