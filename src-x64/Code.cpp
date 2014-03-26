// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

#include "Code.h"

using namespace std; 
using namespace arma;
using namespace Rcpp; 


uvec myRank( colvec x){
  
   int n = x.n_elem; 
   uvec ind(n); 
   ind(0) = 1; 
   for( uvec::iterator i = ind.begin()+1; i != ind.end(); i++){
     *i = *(i-1)+1;  
   }  
   
   uvec sind = stable_sort_index(stable_sort_index(x));
    
   return ind.elem(sind); 
}



mat CSumI( colvec yy, int FUN, colvec Yi, mat Vi, bool v){
  //FUN == 0 for < 
  //    == 1 for <=
  //    == 2 for > 
  //    == 4 for >=
  
  int yy_n = yy.n_elem; 
  int Yi_n = Yi.n_elem; 


  if(FUN==0 || FUN == 4){
     yy = -yy; 
     Yi = -Yi; 
  }

  colvec all_y(yy_n+Yi_n); 
  for( int i=0; i<yy_n; i++){  
    all_y(i) = yy(i); 
  }

  for( int j=yy_n; j < Yi_n + yy_n; j++){
    all_y(j) = Yi(j-yy_n); 
  }

  uvec pos = myRank(all_y).subvec(0, yy_n-1) - myRank(yy);   

  //if there is an equality
  if( FUN==1 || FUN ==4){
    pos = Yi_n - pos; 
  }

  if(v){
     uvec tmpind = stable_sort_index(Yi); 

     if(FUN ==1 || FUN == 4){ 
        tmpind = stable_sort_index(-Yi); 
     }
        
     
      
      mat newVi = cumsum(Vi.rows(tmpind));
 
      newVi.insert_rows(0,1);
      

      return newVi.rows(pos); 

   }else return conv_to<colvec>::from(pos); 

}

mat Vec2Mat(colvec yy, int Nrows){
   int yy_n = yy.n_elem;
 
   mat out(Nrows, yy_n); 
   out.each_row() = yy.t(); 

   return out; 

}

   
colvec myPmin( colvec myvec, double x){
    
   myvec(find(myvec > x)).fill(x); 
   return myvec; 
}



SEXP getWEXP(SEXP datar_SEXP, SEXP Yr_SEXP, SEXP N_SEXP, SEXP RT_outr_SEXP, SEXP predictTime_SEXP, SEXP resid_sco_SEXP, SEXP fitvar_SEXP){
   //convert SEXP to arma
 
   arma::mat data = Rcpp::as<arma::mat>(datar_SEXP);
   arma::mat Y = Rcpp::as<arma::mat>(Yr_SEXP);
   arma::mat RT_out = Rcpp::as<arma::mat>(RT_outr_SEXP);
   arma::vec resid_sco = Rcpp::as<arma::vec >(resid_sco_SEXP);
   
   int N = Rcpp::as<int>(N_SEXP);
   double predictTime = Rcpp::as<double>(predictTime_SEXP); 
   double fitvar = Rcpp::as<double>(fitvar_SEXP); 
   
   //create dataD
   arma::mat dataDuns = data.rows(find(data.col(1)==1));
   arma::mat dataD = dataDuns.rows(sort_index(dataDuns.col(0)));


   int n = data.n_rows, nD = dataD.n_rows, np = Y.n_cols; 
   //guide
    colvec times = data.col(0); 
   // status= data.col(1); 
   // Dtimes = dataD.col(0); 
   // weights = data.col(4); 
   //  Dweights = dataD.col(4); 

   uvec tmpind = find((data.col(0) <= predictTime)%data.col(1)); // indices of (data$times<=predict.time)&(data$status==1)
   mat myempty(1,1); 
   uvec tmpindT = conv_to<uvec>::from(CSumI(times.elem(tmpind), 4, dataD.col(0), myempty,  FALSE)); 


   colvec rrk =  exp(data.col(6)); //data.col(6)  is data$linearY
   
   // build riskmat
   mat riskmat(n, nD); 
   mat::iterator riskmat_it = riskmat.begin(); 
 
for(colvec::iterator j = dataD.begin_col(0); j != dataD.end_col(0); j++){   
   for( colvec::iterator i = data.begin_col(0); i != data.end_col(0); i++){
       *riskmat_it = *i >= *j; riskmat_it++;
      }  
   }   
   //s0 and s1
   colvec s0 = riskmat.t()*(rrk%data.col(4)); 
   colvec s1 = riskmat.t()*trans(Vec2Mat(rrk%data.col(4), np)%Y.t()); 

   //haz0 and cumhaz0
   
   colvec haz0 = dataD.col(4)/sum(riskmat%trans(Vec2Mat(rrk%data.col(4), nD))).t(); 
   colvec cumhaz0 = cumsum(haz0); 
   colvec ptvec(1); ptvec(0) = predictTime; 
   colvec cumhaz_t0 = CSumI(ptvec, 4, dataD.col(0), haz0, TRUE); 

   //Wexp
   colvec Wexp_beta = resid_sco*fitvar*N; 

   colvec WexpLam1(n); WexpLam1.zeros(n); 
   WexpLam1(tmpind) =  N/s0(tmpindT - 1);     
   WexpLam1 = WexpLam1 - CSumI( myPmin(data.col(0), predictTime), 4, dataD.col(0), haz0/s0, TRUE)%rrk*N; 
 
   colvec WexpLam2 = Wexp_beta*CSumI(ptvec, 4, dataD.col(0), haz0%s1/trans(Vec2Mat(s0, np)), TRUE); 
   colvec WexpLam = WexpLam1 - WexpLam2; 

   //Fyk = Pr(Sy < c)
   colvec Fyk = CSumI( data.col(6), 4, data.col(6), data.col(4), TRUE)/sum(data.col(4)); 

   colvec dFyk(Fyk.n_elem); dFyk(0) = 0;
   dFyk(span(1, Fyk.n_elem - 1)) = Fyk(span(0, Fyk.n_elem-2));
   dFyk = Fyk - dFyk; 
   
   colvec Sy = data.col(5); 
   colvec St0_Fyk = cumsum(Sy%dFyk); 
   double St0 = max(St0_Fyk); 
   colvec St0_Syk = St0 - St0_Fyk; 
 //
  mat Wexp_Cond_Stc =  -Vec2Mat(Sy%rrk, n)%(trans(Vec2Mat(WexpLam, n)) +as_scalar(cumhaz_t0)*Wexp_beta*Y.t()); 

  mat tmpmat = conv_to<mat>::from(trans(Vec2Mat(data.col(6), n)) > Vec2Mat(data.col(6), n)); 
  mat Wexp_Stc = trans(CSumI(data.col(6), 0, data.col(6), Wexp_Cond_Stc.t()%Vec2Mat(dFyk, n).t(), TRUE)) + trans(Vec2Mat(Sy,n))%tmpmat - Vec2Mat(St0_Syk, n); 

   colvec Wexp_St = sum(trans(Wexp_Cond_Stc)%trans(Vec2Mat(dFyk, dFyk.n_elem))).t() + Sy - St0; 

   mat Wexp_Fc = 1-tmpmat - Vec2Mat(Fyk, n); 

   //assemble for classic performance measures, given linear predictor

   List out(8); 

   out[0] = -Wexp_Cond_Stc; 
   out[1] =  Wexp_Fc;
   mat Wexp_St_mat = Vec2Mat(Wexp_St, n).t(); 

   out[2] = (-Wexp_St_mat%Vec2Mat(RT_out.col(3), n) + Wexp_Stc)/St0; 
   out[3] = (Wexp_St_mat%Vec2Mat(RT_out.col(4), n) -Wexp_Fc - Wexp_Stc)/(1-St0); 

   out[4] = -Wexp_St; 
   out[5] = (Wexp_St_mat - Wexp_Stc - Vec2Mat(RT_out.col(6), n)%Wexp_Fc)/Vec2Mat(Fyk, n); 
   out[6] = (Vec2Mat(RT_out.col(5)-1, n)%Wexp_Fc - Wexp_Stc)/Vec2Mat(1-Fyk, n);
   out[7] = Wexp_beta;

   return out; 
}








