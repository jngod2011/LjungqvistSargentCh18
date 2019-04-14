## solving by the L&S method
begin
    using JLD2, FileIO
    using LinearAlgebra
    using Interpolations
    using Roots
    using Plots
    using Optim
    pyplot()
    #CPROB=[1.0/3.0 2.0/3.0;1.0/4.0 3.0/4.0];
    #CPROB=[0.35 0.65;0.22 0.78];
    CPROB=[1.0/3.0 2.0/3.0;1.0/4.0 3.0/4.0];

    nos=size(CPROB)[1];
    [sum(CPROB[i,:]) for i=1:nos];

    mM=CPROB'-I;
    mM[end,:]=[1 for i=1:nos];
    rhs=[floor(i/nos) for i=1:nos];
    BSP=mM\rhs; #the list of steady-state probabilities

    hours = [0.0,1.0];
    wage=1.0;
    inc=wage*hours;
    ι=inc'*BSP; #payoff from the bond (= the aggregate output in the economy)

    #ρ=0.05
    ρ=-log(0.96)
    R=3.0;
end

function U(ξ::Float64,R::Float64)
    if ξ <= 0.
        -9e10
        else ((1/ξ^(R-1.0))/(1.0 - R))
    end
end

function II(aa::Int64,k::Int64,a::Int64,tbl2::Array{Int64,2})
    if aa==tbl2[k,a]
        return 1
    else
        return 0
    end
end

## the step on the grid created above

function LS_method(intr::Float64,accu_pow1::Int64,accu_pow2::Int64,nos::Int64,hours::Array{Float64,1},wage::Float64,β::Float64,grda::Array{Float64,1},Lg::Int64,CPROB::Array{Float64,2},R::Float64)
    local TBL,TBL1,TBL2,tbl,tbl1,tbl2,Λ0,Λ1
    currentU=[U((1+intr)*grda[i]+wage*hours[k]-grda[j],R) for k=1:nos, i=1:Lg, j=1:Lg];
    nextV=[U((1+intr)*grda[j]+wage*hours[k],R) for k=1:nos, j=1:Lg];

    TBL=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*nextV)') for k=1:nos, i=1:Lg];
    TBL1=first.(TBL);TBL2=last.(TBL);

    tbl=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*TBL1)') for k=1:nos, i=1:Lg];
    tbl1=first.(tbl);tbl2=last.(tbl);

    while maximum(abs.(TBL2-tbl2))>0||maximum(abs.(TBL1-tbl1))>1/(10^accu_pow1)
        TBL1=tbl1;
        TBL2=tbl2;
        tbl=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*TBL1)') for k=1:nos, i=1:Lg];
        tbl1=first.(tbl);tbl2=last.(tbl);
    end

    Λ0=[Float64(1/(nos*Lg)) for k=1:nos, i=1:Lg];

    Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];

    Λ0=Λ1;
    Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    while maximum(abs.(Λ1-Λ0))>1/10^accu_pow2
        Λ0=Λ1;
        Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    end
    return Λ1,tbl2
end
