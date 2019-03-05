###
###  Run Julia on the entire file
###########################
include("LS-method-7s-init.jl");


## first take the bound from the AL method
w_bounds=(-1.6292720455478487, 7.18257249498523);
begin
    Lg=2000;
    grda=[w_bounds[1]+i*(w_bounds[end]-w_bounds[1])/(Lg-1) for i=0:(Lg-1)];
    β=exp(-ρ);
end


demands=[]

## trial interst rates
intr_trials=[0.03701851068729933,0.03,0.02,0.025,0.0275,0.026250000000000002,0.025625000000000002,0.025937500000000002,0.025781250000000002,0.025703125]

intr=intr_trials[end]
#intr=0.03701851068729933

## load previously calculated values
@load "LS_7s_last_saved.jld2" LS_7s_last_saved;
Λ1,tbl2=LS_7s_last_saved[1],LS_7s_last_saved[2];
intr

## those are important checkpoints
sum(Λ1)
[sum(Λ1[i,:]) for i=1:nos]-BSP

begin
    demand=[grda[Int64(tbl2[kk,ii])] for kk=1:nos, ii=1:Lg];
    demands=vcat(demands,sum(Λ1.*demand));
end

demands[end]
## asset holdings computed with the trial interst rates
results=[2.5962400597605844,2.1641391789354287,-1.1063492587202628,-0.9458937477350954,1.3390600887510085,0.9675183139405708,-0.9219084160642973,0.8951184111096828,0.8437074739560597,-0.9169339221019778]


scatter(results,intr_trials,label="",color=:black,border=true)
plot!((0.02:0.0001:0.035)*0,0.02:0.0001:0.035,color=:black,label="")
xlabel!("average assets")
ylabel!("interest rate")

A=hcat(intr_trials,results)
Aplus=A[A[:,2].>0,:]
Aminus=A[A[:,2].<=0,:]

iloc1=findmin(Aplus[:,2])[2]
iloc2=findmax(Aminus[:,2])[2]
Aplus[iloc1,1],Aminus[iloc2,1]
Aplus[iloc1,2],Aminus[iloc2,2]
#the gap in the wealth is
Aplus[iloc1,2]-Aminus[iloc2,2]

#the gap in the interest rate is
Aplus[iloc1,1]-Aminus[iloc2,1]


begin
    bond_price=ι/(1+intr)
    cons=[(1+intr)*grda[i]+wage*hours[k]-grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
    shares=[grda[Int64(tbl2[k,i])]/bond_price for k=1:nos, i=1:Lg];
end

begin
    SL=[sum(Λ1[state,:]) for state=1:nos];
    LF=[[Λ1[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+Λ1[state,i]/SL[state])
        end
    end
end


begin
    plot(aassets[1,:],LF[1],label="(LS)",color=:black)
    plot!(aassets[2,:],LF[2],label="",color=:black)
    plot!(aassets[3,:],LF[3],label="",color=:black)
    plot!(aassets[4,:],LF[4],label="",color=:black)
    plot!(aassets[5,:],LF[5],label="",color=:black)
    plot!(aassets[6,:],LF[6],label="",color=:black)
    plot!(aassets[7,:],LF[7],label="",color=:black)
end


begin
    plot(aassets[1,1:1200],LF[1][1:1200],label="",color=:black)
    plot!(aassets[2,1:1200],LF[2][1:1200],label="",color=:black)
    plot!(aassets[3,1:1200],LF[3][1:1200],label="",color=:black)
    plot!(aassets[4,1:1200],LF[4][1:1200],label="",color=:black)
    plot!(aassets[5,1:1200],LF[5][1:1200],label="",color=:black)
    plot!(aassets[6,1:1200],LF[6][1:1200],label="",color=:black)
    plot!(aassets[7,1:1200],LF[7][1:1200],label="",color=:black)
end
