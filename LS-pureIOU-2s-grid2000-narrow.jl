#
# A single execution of this file generates a single pair (a(r),r).
# The pair is then stored in "LS_pureIOU_2s_2000_narrow_saved.jld2", together with
# all previously computed pairs.
#      If running for the first time (no "LS_pureIOU_2s_2000_narrow_saved.jld2" exists)
#      uncomment lines 24-29. Put the #'s back for all subsequent runs.
#
include("LS-method-2s-init.jl");


# Assume the bounds on the wealth generated by the dual method and create the grid.
#
w_bounds=(-3.1534299198619835, 3.828835045437629);

begin
    Lg=2000 #the only place to change the grid size
    grda=[w_bounds[1]+i*(w_bounds[end]-w_bounds[1])/(Lg-1) for i=0:(Lg-1)]
    β=exp(-ρ)
end;

## Uncomment and execute the next block the very *first* time,
## *only* if "LS_pureIOU_2s_2000_narrow_saved.jld2" does not exist in the current path.
## Put the #'s back for all subsequent runs.

#begin
#    a_list=Array{Float64}(undef, 0)
#    r_list=Array{Float64}(undef, 0)
#    LS_pureIOU_2s_2000_narrow_saved=[r_list,a_list]
#    @save "LS_pureIOU_2s_2000_narrow_saved.jld2" LS_pureIOU_2s_2000_narrow_saved
#end;

##
## The next block guesses the next candidate for r.
## It can be bypassed by setting r=<value> explicitly. See line 68.
##


begin      
    @load "LS_pureIOU_2s_2000_narrow_saved.jld2" LS_pureIOU_2s_2000_narrow_saved;
    r_list=LS_pureIOU_2s_2000_narrow_saved[1]
    a_list=LS_pureIOU_2s_2000_narrow_saved[2]
    if length(r_list)==0
        r=0.03
        elseif length(r_list)==1
        r=0.02
        else 
        r=0.5*(r_list[findall(x->x<0,a_list)[end]]+r_list[findall(x->x>0,a_list)[end]])
    end
    push!(r_list,r)
end;



# the next line could take up to an hour
# arguments 2 and 3are the requested digits of convergence
#
@time Λ1,tbl2=LS_method(r,5,5,nos,hours,wage,β,grda,Lg,CPROB,R);

## Two important checkpoints!
sum(Λ1) # should print 1
[sum(Λ1[i,:]) for i=1:nos]-BSP # should print a list of zeros

begin
    policy=[grda[Int64(tbl2[kk,ii])] for kk=1:nos, ii=1:Lg]
    a=sum(Λ1.*policy) # the aggregate demand a(r)
end;

##
## If running with an explicit r (see line 32)
## be sure that a_list and r_list exist before running the next block.
##

begin
    a_list=push!(a_list,a)
    LS_pureIOU_2s_2000_narrow_saved=[r_list,a_list,Λ1,tbl2]
    @save "LS_pureIOU_2s_2000_narrow_saved.jld2" LS_pureIOU_2s_2000_narrow_saved
end;