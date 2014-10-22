/** 	\file
	\brief Macros to do vectorization.
	*/

// Function body type declarations
// V=vector, T=scalar, I=integer, N=none
#define declareV(arg) const vector<Type> &arg
#define declareT(arg) Type arg
#define declareI(arg) int arg
#define declareN(arg)
// How to extract elementwise subset of the four types
#define elementV(arg,i) arg[i]
#define elementT(arg,i) arg
#define elementI(arg,i) arg
#define elementN(arg,i)
// How to place comma in front of the types
#define commaV ,
#define commaT ,
#define commaI ,
#define commaN
// Update output vector size
#define outputsizeV(n,arg) n = (arg.size()>n ? arg.size() : n)
#define outputsizeT(n,arg)
#define outputsizeI(n,arg)
#define outputsizeN(n,arg)
/** \brief General vectorize macro up to 6 arguments

    Applied type abbreviations: V=vector, T=scalar, I=integer, N=none.
    The longest vector input determines the length of the output.
    Arguments are not re-cycled; unequal vector lengths should result
    in a crash.
*/
#define GVECTORIZE(FUN,Type1,Type2,Type3,Type4,Type5,Type6)		\
template <class Type>							\
vector<Type> FUN( declare##Type1(arg1) comma##Type2			\
		  declare##Type2(arg2) comma##Type3			\
		  declare##Type3(arg3) comma##Type4			\
		  declare##Type4(arg4) comma##Type5			\
		  declare##Type5(arg5) comma##Type6			\
		  declare##Type6(arg6) )				\
{									\
  int n = 1;								\
  outputsize##Type1(n,arg1);						\
  outputsize##Type2(n,arg2);						\
  outputsize##Type3(n,arg3);						\
  outputsize##Type4(n,arg4);						\
  outputsize##Type5(n,arg5);						\
  outputsize##Type6(n,arg6);						\
  vector<Type> res(n);							\
  for(int i=0;i<n;i++) res[i] = FUN( element##Type1(arg1,i) comma##Type2 \
				     element##Type2(arg2,i) comma##Type3 \
				     element##Type3(arg3,i) comma##Type4 \
				     element##Type4(arg4,i) comma##Type5 \
				     element##Type5(arg5,i) comma##Type6 \
				     element##Type6(arg6,i) );		\
  return res;								\
}

/** \brief Vectorize 1-argument functions. */
#define VECTORIZE1_t(FUN)			\
  GVECTORIZE(FUN,V,N,N,N,N,N)

/** \brief Vectorize 2-argument functions.

    For two-arguments functions (Type, Type),
    vectorize both arguments.
*/
#define VECTORIZE2_tt(FUN)			\
  GVECTORIZE(FUN,V,T,N,N,N,N);			\
  GVECTORIZE(FUN,T,V,N,N,N,N);			\
  GVECTORIZE(FUN,V,V,N,N,N,N);

/** \brief Vectorize 3-argument functions.

    For three-arguments functions (Type, Type, int),
    vectorize first two arguments.
*/
#define VECTORIZE3_tti(FUN)			\
  GVECTORIZE(FUN,V,T,I,N,N,N);			\
  GVECTORIZE(FUN,T,V,I,N,N,N);			\
  GVECTORIZE(FUN,V,V,I,N,N,N);

/** \brief Vectorize 3-argument functions.

    For three-arguments functions (Type, Type, Type),
    vectorize all three arguments.
*/
#define VECTORIZE3_ttt(FUN)			\
  GVECTORIZE(FUN,V,T,T,N,N,N);			\
  GVECTORIZE(FUN,T,V,T,N,N,N);			\
  GVECTORIZE(FUN,T,T,V,N,N,N);			\
  GVECTORIZE(FUN,V,V,T,N,N,N);			\
  GVECTORIZE(FUN,T,V,V,N,N,N);			\
  GVECTORIZE(FUN,V,T,V,N,N,N);			\
  GVECTORIZE(FUN,V,V,V,N,N,N);

/** \brief Vectorize 4-argument functions.

    For Four-arguments functions (Type, Type, Type, int),
    vectorize first three arguments.
*/
#define VECTORIZE4_ttti(FUN)			\
  GVECTORIZE(FUN,V,T,T,I,N,N);			\
  GVECTORIZE(FUN,T,V,T,I,N,N);			\
  GVECTORIZE(FUN,T,T,V,I,N,N);			\
  GVECTORIZE(FUN,V,V,T,I,N,N);			\
  GVECTORIZE(FUN,T,V,V,I,N,N);			\
  GVECTORIZE(FUN,V,T,V,I,N,N);			\
  GVECTORIZE(FUN,V,V,V,I,N,N);

/** \brief Vectorize 6-argument functions.

    For Six-arguments functions (Type, Type, Type, Type, Type, int),
    vectorize first three arguments.
*/
#define VECTORIZE6_ttttti(FUN)			\
  GVECTORIZE(FUN,V,T,T,T,T,I);			\
  GVECTORIZE(FUN,T,V,T,T,T,I);			\
  GVECTORIZE(FUN,T,T,V,T,T,I);			\
  GVECTORIZE(FUN,V,V,T,T,T,I);			\
  GVECTORIZE(FUN,T,V,V,T,T,I);			\
  GVECTORIZE(FUN,V,T,V,T,T,I);			\
  GVECTORIZE(FUN,V,V,V,T,T,I);


using CppAD::abs;
VECTORIZE1_t(abs);
VECTORIZE1_t(acos);
VECTORIZE1_t(asin);
VECTORIZE1_t(atan);
VECTORIZE1_t(cos);
VECTORIZE1_t(erf);
VECTORIZE1_t(exp);
VECTORIZE1_t(log);
VECTORIZE1_t(log10);
VECTORIZE1_t(sin);
VECTORIZE1_t(sqrt);   //mangler atan2, pow
VECTORIZE2_tt(pow);

/* max and min of vector */
double max(const vector<double> &x)
{
  double res=x[0];
  for(int i=0;i<x.size();i++){
    //res=CondExpGt(res,x[i],res,x[i]);
    if(res<x[i])res=x[i];
  }
  return res;
}
template <class Type>
Type max(const vector<Type> &x)
{
  Type res=x[0];
  for(int i=0;i<x.size();i++){
    res=CondExpGt(res,x[i],res,x[i]);
  }
  return res;
}
double min(const vector<double> &x)
{
  double res=x[0];
  for(int i=0;i<x.size();i++){
    //res=CondExpGt(res,x[i],res,x[i]);
    if(res>x[i])res=x[i];
  }
  return res;
}
template <class Type>
Type min(const vector<Type> &x)
{
  Type res=x[0];
  for(int i=0;i<x.size();i++){
    res=CondExpLt(res,x[i],res,x[i]);
  }
  return res;
}


/* Normalize rows and colums in a matrix */
// template<class Type>
// void normalizeRows(matrix<Type> &x){
//   for(int i=0;i<x.rows();i++){
//     matrix_row<matrix<Type> > xr(x,i);
//     xr/=sum(xr);
//   }
// }
// template<class Type>
// void normalizeCols(matrix<Type> &x){
//   for(int i=0;i<x.cols();i++){
//     matrix_column<matrix<Type> > xc(x,i);
//     xc/=sum(xc);
//   }
// }


/* Tape the choice of a matrix-row  */
template<class Type>
vector<Type> select_row(const matrix<Type> &x, Type index)
{
   int n1=x.rows();
   int n2=x.cols();
   vector<Type> res(n2);
   //res.clear();
   for(int j=0;j<n2;j++)res(j)=Type(0);
   for(int i=0;i<n1;i++)
     for(int j=0;j<n2;j++)
        res(j)+=CppAD::CondExpEq(index,Type(i),x(i,j),Type(0));
   return res;
}

/* Tape the choice of a vector-element  */
template<class Type>
Type select_elt(const vector<Type> &x, Type i)
{
   int n=x.size();
   Type res(0.0);
   for(int j=0;j<n;j++)
     {
        res+=CppAD::CondExpEq(i,Type(j),x[j],Type(0));
     }
   return res;
}

/* Generate piecewice constant function-object based on n discontinuity-points
   and n+1 function values */
template<class Type>
class piecewice{
private:
  Type y0;
  vector<Type> x,dy;
  int n;
public:
  bool leftcontinuous;
  piecewice(const vector<Type> &x_, const vector<Type> &y,
	    bool leftcontinuous_=true):
    x(x_),leftcontinuous(leftcontinuous_){
    n=x.size();
    if(n!=y.size()-1)error("piecewice: y must be one longer than x");
    dy.resize(n);
    y0=y(0);
    for(int i=0;i<n;i++)dy(i)=y(i+1)-y(i);
  }
  inline Type operator()(const Type &t){
    Type res=y0;
    if(leftcontinuous){
      for(int i=0;i<n;i++)res+=CppAD::CondExpGt(t,x(i),dy(i),Type(0));
    } else {
      for(int i=0;i<n;i++)res+=CppAD::CondExpGe(t,x(i),dy(i),Type(0));
    }
    return res;
  }
}; 
