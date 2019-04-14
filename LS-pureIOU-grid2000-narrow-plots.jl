include("LS-method-7s-init.jl");

begin
    @load "LS_pureIOU_2000_narrow_saved.jld2" LS_pureIOU_2000_narrow_saved;
    r_list=LS_pureIOU_2000_narrow_saved[1]
    a_list=LS_pureIOU_2000_narrow_saved[2]
    Λ1=LS_pureIOU_2000_narrow_saved[3]
    tbl2=LS_pureIOU_2000_narrow_saved[4]
end;

## Two important checkpoints!
sum(Λ1) # should print 1
[sum(Λ1[i,:]) for i=1:nos]-BSP # should print a list of zeros

begin
    A=hcat(r_list,a_list)
    Aplus=A[A[:,2].>0,:]
    Aminus=A[A[:,2].<=0,:]
    iloc1=findmin(Aplus[:,2])[2]
    iloc2=findmax(Aminus[:,2])[2]
end;

#the nearest interest rates
Aplus[iloc1,1],Aminus[iloc2,1]

#the respective aggregate demand
Aplus[iloc1,2],Aminus[iloc2,2]

#the gap in the demand a(r)
Aplus[iloc1,2]-Aminus[iloc2,2]

#the gap in the interest rate is
Aplus[iloc1,1]-Aminus[iloc2,1]

begin
    scatter(a_list,r_list,label="",color=:black,border=true)
    plot!((0.02:0.0001:0.032)*0,0.02:0.0001:0.032,color=:black,label="")
    xlabel!("average assets")
    ylabel!("interest rate")
end


begin
    r=r_list[end]
    w_bounds=(-1.6292720455478487, 7.18257249498523);
    Lg=2000 #the only place to change the grid size
    grda=[w_bounds[1]+i*(w_bounds[end]-w_bounds[1])/(Lg-1) for i=0:(Lg-1)]
    β=exp(-ρ)
end;


begin
    bond_price=ι/(1+r)
    cons=[(1+r)*grda[i]+wage*hours[k]-grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
    shares=[grda[Int64(tbl2[k,i])]/bond_price for k=1:nos, i=1:Lg];
end;

begin
    SL=[sum(Λ1[state,:]) for state=1:nos];
    LF=[[Λ1[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+Λ1[state,i]/SL[state])
        end
    end
end;


begin
    plot(aassets[1,2:920],LF[1][2:920],label="(LS)",color=:black)
    plot!(aassets[2,2:920],LF[2][2:920],label="",color=:black)
    plot!(aassets[3,2:920],LF[3][2:920],label="",color=:black)
    plot!(aassets[4,2:920],LF[4][2:920],label="",color=:black)
    plot!(aassets[5,2:920],LF[5][2:920],label="",color=:black)
    plot!(aassets[6,2:920],LF[6][2:920],label="",color=:black)
    plot!(aassets[7,2:920],LF[7][2:920],label="",color=:black)
end
