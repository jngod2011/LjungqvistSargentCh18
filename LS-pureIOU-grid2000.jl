####
### Run the blocks manually
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

## REPEAT FROM HERE
##  manually insert on the next line the ansatz values for r one by one
## the very first rate is the one obtained by the AL nethod
intr_trials=[0.03701851068729933,0.03,0.02,0.025,0.0275,0.026250000000000002,0.025625000000000002,0.025937500000000002,0.025781250000000002,0.025703125]
### use the next line to compute the average of two previous candidate rates
### (one too big the other one too small)
###
#(intr_trials[end]+intr_trials[end-2])/2
intr=intr_trials[end]
#intr=0.03701851068729933

#the next line could take up to an hour
@time Λ1,tbl2=LS_method(intr,5,5,nos,hours,wage,β,grda,Lg,CPROB,R);
## uncomment the next two lines to save the output
#LS_7s_last_saved=[Λ1,tbl2];
#@save "LS_7s_last_saved.jld2" LS_7s_last_saved;
intr

## those are important checkpoints
sum(Λ1)
[sum(Λ1[i,:]) for i=1:nos]-BSP

begin
    demand=[grda[Int64(tbl2[kk,ii])] for kk=1:nos, ii=1:Lg];
    demands=vcat(demands,sum(Λ1.*demand));
end

## produce the aggregate demand corresponding to the last interest rate
demands[end]
## Copy the output from the above line, add it to the list below, and go back to the REPEAT point
## The latest demand value may already be added
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
