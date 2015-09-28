#include <TMB.hpp>

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
    pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
    return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
DATA_SCALAR(eps);
 DATA_VECTOR(x);

 PARAMETER(p);
 PARAMETER(Dummy);
 Type pen;
 pen=0;

 Type nll;
 nll=0;

 Type var;
 var = p*(1-p);
 var = posfun(var, eps, pen);
 nll += pen;

 nll += -sum( dnorm(x, p, sqrt(var),true) );
 nll += -dnorm(Dummy, Type(0.0), Type(1.0), true);
 return nll;

}
