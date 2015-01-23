/**
   \file R_inla.cpp
   Makes SPDE methods from INLA R-package available in TMB.
*/

/** \brief SPDE methods from INLA R-package .

*/
using namespace Eigen;
using namespace tmbutils;

namespace R_inla {
	
template<class Type>
struct spde_t{  
  SparseMatrix<Type> G0;
  SparseMatrix<Type> G1;
  SparseMatrix<Type> G2;
  spde_t(SEXP x){  /* x = List passed from R */
  G0 = asSparseMatrix<Type>(getListElement(x,"G0"));
  G1 = asSparseMatrix<Type>(getListElement(x,"G1"));
  G2 = asSparseMatrix<Type>(getListElement(x,"G2"));
}
};

template<class Type>
struct spde_aniso_t{
  int n_s;
  int n_tri;
  vector<Type> Tri_Area;
  matrix<Type> E0;
  matrix<Type> E1;
  matrix<Type> E2;
  vector<int> TV;
  SparseMatrix<Type> G0;
  SparseMatrix<Type> G0_inv;
  spde_aniso_t(SEXP x){  /* x = List passed from R */
  n_s = 	CppAD::Integer(asVector<Type>(getListElement(x,"n_s"))[0]);  
  n_tri = 	CppAD::Integer(asVector<Type>(getListElement(x,"n_tri"))[0]);  
  Tri_Area = asVector<Type>(getListElement(x,"Tri_Area"));
  E0 = asMatrix<Type>(getListElement(x,"E0"));
  E1 = asMatrix<Type>(getListElement(x,"E1"));
  E2 = asMatrix<Type>(getListElement(x,"E2"));  
  TV = asVector<int>(getListElement(x,"TV")); 
  G0 = asSparseMatrix<Type>(getListElement(x,"G0"));
  G0_inv = asSparseMatrix<Type>(getListElement(x,"G0_inv"));
}
};

/** Precission matrix eqn (10) in Lindgren et al. (2011) */    
template<class Type>
  SparseMatrix<Type> Q_spde(spde_t<Type> spde, Type kappa){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  	
  return kappa_pow4*spde.G0 + Type(2.0)*kappa_pow2*spde.G1 + spde.G2;    
}


/** Precission matrix for the anisotropic case, eqn (20) in Lindgren et al. (2011) */    
template<class Type>
  SparseMatrix<Type> Q_spde(spde_aniso_t<Type> spde, Type kappa, vector<Type> ln_H_input){

  int i;
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  int n_s = spde.n_s;
  int n_tri = spde.n_tri;
  vector<Type> Tri_Area = spde.Tri_Area;
  matrix<Type> E0 = spde.E0;
  matrix<Type> E1 = spde.E1;
  matrix<Type> E2 = spde.E2;
  vector<int> TV = spde.TV;
  SparseMatrix<Type> G0 = spde.G0;
  SparseMatrix<Type> G0_inv = spde.G0_inv;
	  	  
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  Type H_trace = H(0,0)+H(1,1);
  Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  SparseMatrix<Type> G1_aniso(n_s,n_s); 
  SparseMatrix<Type> G2_aniso(n_s,n_s); 
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  // Calculate new SPDE matrices

  // Calculate G1 - pt. 1
  array<Type> Gtmp(n_tri,3,3);
  for(i=0; i<n_tri; i++){    
    // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
    Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
  }
  // Calculate G1 - pt. 2
  for(i=0; i<n_tri; i++){
    int i0 = i;
    int i1 = i + n_tri; 
    int i2 = i + 2*n_tri; 
    G1_aniso.coeffRef(TV(i1),TV(i0)) = G1_aniso.coeffRef(TV(i1),TV(i0)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i0),TV(i1)) = G1_aniso.coeffRef(TV(i0),TV(i1)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i2),TV(i1)) = G1_aniso.coeffRef(TV(i2),TV(i1)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i1),TV(i2)) = G1_aniso.coeffRef(TV(i1),TV(i2)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i2),TV(i0)) = G1_aniso.coeffRef(TV(i2),TV(i0)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i0),TV(i2)) = G1_aniso.coeffRef(TV(i0),TV(i2)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i0),TV(i0)) = G1_aniso.coeffRef(TV(i0),TV(i0)) + (Gtmp(i,0,0));  
    G1_aniso.coeffRef(TV(i1),TV(i1)) = G1_aniso.coeffRef(TV(i1),TV(i1)) + (Gtmp(i,1,1));  
    G1_aniso.coeffRef(TV(i2),TV(i2)) = G1_aniso.coeffRef(TV(i2),TV(i2)) + (Gtmp(i,2,2));  
  }
  G2_aniso = G1_aniso * G0_inv * G1_aniso; 

  return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
}


} // end namespace R_inla
