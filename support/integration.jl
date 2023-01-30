using Symbolics
using SymbolicNumericIntegration

@variables x_m
@variables x_l
@variables t_0
@variables N
@variables l
@variables m

function P(x)
    return  (x-t_0)/(1-t_0)
end
function p(x)
    return 1/(1-t_0)
end

integrate(p(x_m))