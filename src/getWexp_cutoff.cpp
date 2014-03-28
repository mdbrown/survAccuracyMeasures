// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "Code.h"

using namespace std; 
using namespace Rcpp;
using namespace arma; 


// [[Rcpp::export]]
Rcpp::List getWEXPcutoff(arma::mat data, arma::mat subdata, arma::mat Y, arma::mat subY, int N, arma::mat RT_out, double predictTime, arma::vec resid_sco, double fitvar, arma::colvec cutoffs){

   
   //create dataD
   arma::mat dataDuns = data.rows(find(data.col(1)==1));
   arma::mat dataD = dataDuns.rows(sort_index(dataDuns.col(0)));


   int n = data.n_rows, nD = dataD.n_rows, np = Y.n_cols, ncuts = subdata.n_rows; 
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
   colvec Fyk = CSumI( cutoffs, 4, data.col(6), data.col(4), TRUE)/sum(data.col(4)); 

   colvec dFyk(Fyk.n_elem); dFyk(0) = 0;
   dFyk(span(1, Fyk.n_elem - 1)) = Fyk(span(0, Fyk.n_elem-2));
   dFyk = Fyk - dFyk; 
   
   colvec Sy = subdata.col(5); 
   colvec Syall = data.col(5); 
   colvec St0_Fyk = cumsum(Sy%dFyk); 
   double St0 = max(St0_Fyk); 
   colvec St0_Syk = St0 - St0_Fyk; 
 //
  mat Wexp_Cond_Stc =  -Vec2Mat(Sy%exp(subdata.col(6)), n)%(trans(Vec2Mat(WexpLam, ncuts)) +as_scalar(cumhaz_t0)*Wexp_beta*subY.t()); 

  mat tmpmat = conv_to<mat>::from(trans(Vec2Mat(data.col(6), ncuts)) > Vec2Mat(cutoffs, n)); 
  mat Wexp_Stc = trans(CSumI(cutoffs, 0, subdata.col(6), Wexp_Cond_Stc.t()%Vec2Mat(dFyk, n).t(), TRUE)) + trans(Vec2Mat(Syall,ncuts))%tmpmat - Vec2Mat(St0_Syk, n); 

  colvec Wexp_St = sum(trans(Wexp_Cond_Stc)%trans(Vec2Mat(dFyk, n))).t() + Syall - St0; 

   mat Wexp_Fc = 1-tmpmat - Vec2Mat(Fyk, n); 

   //assemble for classic performance measures, given linear predictor

   List out(8); 

   out[0] = -Wexp_Cond_Stc; 
   out[1] =  Wexp_Fc;
   mat Wexp_St_mat = Vec2Mat(Wexp_St, ncuts).t(); 

   out[2] = (-Wexp_St_mat%Vec2Mat(RT_out.col(3), n) + Wexp_Stc)/St0; 
   out[3] = (Wexp_St_mat%Vec2Mat(RT_out.col(4), n) -Wexp_Fc - Wexp_Stc)/(1-St0); 

   out[4] = -Wexp_St; 
   out[5] = (Wexp_St_mat - Wexp_Stc - Vec2Mat(RT_out.col(6), n)%Wexp_Fc)/Vec2Mat(Fyk, n); 
   out[6] = (Vec2Mat(RT_out.col(5)-1, n)%Wexp_Fc - Wexp_Stc)/Vec2Mat(1-Fyk, n);
   out[7] = Wexp_beta;

   return out; 
}

