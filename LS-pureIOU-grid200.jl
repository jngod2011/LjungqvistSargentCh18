####
### Run the blocks manually
###########################
include("LS-method-7s-init.jl");


demands=[]

## REPEAT FROM HERE
##  manually insert on the next line the ansatz values for r one by one
## the very first rate is the one obtained by the AL nethod
intr_trials=[0.037018453344873636,0.03,0.025,0.0275,0.028749999999999998,0.029375,0.029062499999999998,0.028906249999999998,0.028984375,0.0289453125,0.028925781249999998,0.028935546875,0.0289306640625]
### use the next line to compute the average of two previous candidate rates
### (one too big the other one too small)
###
#(intr_trials[end]+intr_trials[end-2])/2
intr=intr_trials[end]

begin
    b=3.0;
    minkap=-b;
    if intr<=0.0
        minkap = -b;
    else
        minkap=-min(b, wage*hours[1]/intr);
    end;
    maxkap=16.0;
end
minkap
-2.0821797333204723

begin
    Lg=200;
    grda=[minkap+i*(maxkap-minkap)/(Lg-1) for i=0:(Lg-1)];
    β=exp(-ρ);
end
## 75 seconds
@time Λ1,tbl2=LS_method(intr,5,5,nos,hours,wage,β,grda,Lg,CPROB,R);

intr

## those are important checkpoints
sum(Λ1)
[sum(Λ1[i,:]) for i=1:nos]-BSP

begin
    demand=[grda[tbl2[kk,ii]] for kk=1:nos, ii=1:Lg];
    demands=vcat(demands,sum(Λ1.*demand));
end

## produce the aggregate demand corresponding to the last interest rate
demands[end]
## Copy the output from the above line, add it to the list below, and go back to the REPEAT point
## The latest demand value may already be added
results=[7.369655374922196,0.5100668459460974,-1.0654201141611441,-0.505665361435287,-0.1113888746946741,0.16311610386135325,0.02379514486718422,-0.04304169618186418,0.01870022538380409,0.015794907576101826,-0.04179338316455158,0.015397628445206781,-0.04146948169368391]
intr_trials[1]-intr_trials[end-1]
0.008082906469873637
#demands is the same as results


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

scatter(results,intr_trials,label="",color=:black,border=true)
plot!((0.025:0.0001:0.035)*0,0.025:0.0001:0.035,color=:black,label="")
xlabel!("average assets")
ylabel!("interest rate")



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
    plot(aassets[1,2:92],LF[1][2:92],label="(LS)",color=:black)
    plot!(aassets[2,2:92],LF[2][2:92],label="",color=:black)
    plot!(aassets[3,2:92],LF[3][2:92],label="",color=:black)
    plot!(aassets[4,2:92],LF[4][2:92],label="",color=:black)
    plot!(aassets[5,2:92],LF[5][2:92],label="",color=:black)
    plot!(aassets[6,2:92],LF[6][2:92],label="",color=:black)
    plot!(aassets[7,2:92],LF[7][2:92],label="",color=:black)
end
