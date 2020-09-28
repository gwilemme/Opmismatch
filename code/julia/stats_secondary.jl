using NBInclude
@nbinclude("startup.ipynb")
timing = 60*30 ## time per computation of equilibrium/optimum 
println("Code takes ", (timing*2*6)/60/60, " hours." )

table = zeros(27,6)

for j in 4:9

    (Acal,Mcal) = open_struc(joinpath(dir_calib,"step4_v$j.csv"),N)
    i = j-3
    
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
    for i in 1:6
        print(" & ",table[index,i])
    end
    println("\\\\")
end

spa = "\\qquad" #"\\qquad\\qquad"

println("Fixed parameters &&&&&&\\\\")
printlineLatex("$spa  matching elasticity \$\\eta\$", 1)
printlineLatex("$spa  bargaining power \$\\beta\$", 2)
index = 3

## calibrated par
println("Calibrated parameters &&&&&&\\\\")
printlineLatex("$spa  mean perfect-match productivity", index)
printlineLatex("$spa  mean mismatching penalty (\\%)", index+1)
printlineLatex("$spa  search elasticity \$\\epsilon\$", index+2)
printlineLatex("$spa  search parameter \$\\xi\$", index+3)
index += 4

println("\\hline\\hline")
printlineLatex("Distance to moments", index)
index += 1

## mom1
println("Occupational employment \$\\mathcal{M}_1\$ &&&&&&\\\\")
printlineLatex("$spa distance", index)
printlineLatex("$spa mean", index+1)
printlineLatex("$spa standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)
index += 5

## mom2
println("Occupational wages \$\\mathcal{M}_2\$ &&&&&&\\\\")
printlineLatex("$spa distance", index)
printlineLatex("$spa mean", index+1)
printlineLatex("$spa standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)
index += 5

## mom3
println("Occupational job-finding rates \$\\mathcal{M}_3\$ &&&&&&\\\\")
printlineLatex("$spa distance", index)
printlineLatex("$spa mean", index+1)
printlineLatex("$spa standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)
index += 5

## mom4
println("Occupation-switching rates \$\\mathcal{M}_4\$ &&&&&&\\\\")
printlineLatex("$spa distance", index)
printlineLatex("$spa mean", index+1)
printlineLatex("$spa standard deviation", index+2)
#printlineLatex("\\qquad\\qquad min", index+3)
#printlineLatex("\\qquad\\qquad max", index+4)

table = zeros(18,6)

for j in 4:9
    
    (Acal,Mcal) = open_struc(joinpath(dir_calib,"step4_v$j.csv"),N)
    i = j-3
    
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
    println("Compare welfare depending on computation of optimal taxes:")
    println("oT1 vs oT2: ", welfare_dec(opA,oT1,Mcal,N)," vs " , welfare_dec(opA,oT2,Mcal,N) )
    println("welfare opt: ", welfare_opt(opA,Mcal,N) )
    oT=oT2 ## choose oT2
    
    ## efficiency loss
    table[index,i] = round((1 - welfare_dec(Acal,RTax,Mcal,N) / welfare_opt(opA,Mcal,N))*100,digits=1)
    table[index+1,i] = round((1 - welfare_dec(A0,Notax,Mcal,N) / welfare_opt(opA,Mcal,N))*100,digits=1)
    index += 2
    
    ## vacancies
    vAcal = vacancies(Acal,Mcal,N)
    vA0 = vacancies(A0,Mcal,N)
    vopA = vacancies(opA,Mcal,N)
    table[index,i] = round((sum(vAcal)/sum(vopA)-1) *100, digits = 1)
    table[index+1,i] = round((sum(vA0)/sum(vopA)-1) *100, digits = 1)
    index += 2
    
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
    for i in 1:6
        print(" & ",table[index,i])
    end
    println("\\\\")
end

println("Fixed parameters &&&&&&\\\\")
printlineLatex("\\qquad\\qquad  matching elasticity \$\\eta\$", 1)
printlineLatex("\\qquad\\qquad  bargaining power \$\\beta\$", 2)
println("\\hline\\hline")
index = 3

## efficiency
println("Efficiency loss (\\%)&&&&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
index += 2

## vacancies
println("Loss of vacancies (\\%)&&&&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
index += 2

## unemployment
println("Unemployment rate (\\%) &&&&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
printlineLatex("\\qquad\\qquad  optimum", index+2)
index += 3

## mismatch
println("Share of mismatches (\\%) &&&&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
printlineLatex("\\qquad\\qquad  optimum", index+2)
index += 3

## Mm
println("Mean-min ratio &&&&&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index)
printlineLatex("\\qquad\\qquad  laissez-faire", index+1)
printlineLatex("\\qquad\\qquad  optimum", index+2)

println("Taxes-production coeff. &&&\\\\")
printlineLatex("\\qquad\\qquad  benchmark", index+3)
printlineLatex("\\qquad\\qquad  laissez-faire", index+4)
printlineLatex("Tax on job-to-job transition", index+5)
