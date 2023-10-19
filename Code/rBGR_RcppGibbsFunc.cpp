# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export()]]
double Gibbs_trunNormC(arma::vec LB_Vec, arma::vec f_LB_a1, arma::vec f_LB_a2,
                       arma::vec UB_Vec, arma::vec g_UB_b1, arma::vec g_UB_b2, 
                       double prior_a1, double prior_a2) {
  int LB_len = LB_Vec.n_elem, UB_len = UB_Vec.n_elem;
  arma::vec B_pt = join_cols(LB_Vec, UB_Vec);
  B_pt = sort(B_pt);
  
  int segNum = LB_len + UB_len + 1;
  arma::vec D_vec(segNum), E_vec(segNum), sgmts_cnst(segNum);
  
  for(int iIdx=0; iIdx < segNum; iIdx++){
    double LB=0, UB=0;
    if(iIdx == 0){
      LB=R_NegInf;
      UB=B_pt(iIdx);
    }else if(iIdx == (segNum - 1)){
      LB=B_pt(iIdx-1);
      UB=R_PosInf;
    }else{
      LB=B_pt(iIdx-1);
      UB=B_pt(iIdx);
    }
    
    double Di = prior_a1, Ei = prior_a2;
    
    arma::uvec LB_idx = (LB_Vec <= LB), UB_idx = (UB_Vec >= UB);
    arma::uvec LB_nIdx = find(LB_idx==1), UB_nIdx = find(UB_idx==1);
    Di = Di + sum(f_LB_a1(LB_nIdx))+sum(g_UB_b1(UB_nIdx));
    Ei = Ei + sum(f_LB_a2(LB_nIdx))+sum(g_UB_b2(UB_nIdx));
    D_vec(iIdx) = Di;
    E_vec(iIdx) = Ei;
    
    double in_mean = -Ei/(2*Di), in_sd= sqrt(-1/(2*Di));
    
    double logp_UB = R::pnorm(UB, in_mean, in_sd, 1, 1);
    double logp_LB = R::pnorm(LB, in_mean, in_sd, 1, 1);
    
    double chkVal = exp(logp_UB - logp_LB);
    double cnst;
    if(iIdx==1 | chkVal == std::numeric_limits<double>::infinity()){
      // Rcout << "iIdx " << iIdx << std::endl;
      cnst = logp_UB;
    }else if(chkVal - 1 ==0){
      cnst = R_NegInf;
    }else{
      cnst = logp_LB + log(exp(logp_UB - logp_LB)-1);
    }
    sgmts_cnst(iIdx) = cnst;
  }
  double maxVal = max(sgmts_cnst);
  arma::vec maxVal_vec(segNum);
  maxVal_vec.fill(maxVal);
  double tmpSum = sum(exp(sgmts_cnst - maxVal_vec));
  double logSumExp = log(tmpSum) + maxVal;
  // Rcout << "LSE: " << logSumExp <<std::endl;

  arma::vec prob_(segNum);
  for(int tmpIdx = 0 ; tmpIdx < segNum ; tmpIdx++){
    prob_(tmpIdx) = exp(sgmts_cnst(tmpIdx) - logSumExp);
  }

  // Function set_seed("set.seed");
  // set_seed(Named("seed")=seed_);
  IntegerVector intSeq = seq_len(segNum);
  NumericVector prob_NV = wrap(prob_);
  IntegerVector sel_idx = sample(intSeq, 1, false, prob_NV);
  

  int sIdx = sel_idx[0] - 1;
  double S_Di = D_vec(sIdx), S_Ei = E_vec(sIdx), S_LB, S_UB;
  double S_in_mean = -S_Ei/(2*S_Di), S_in_sd= sqrt(-1/(2*S_Di));
  if(sIdx == 0){
    S_LB=R_NegInf;
    S_UB=B_pt(sIdx);
  }else if(sIdx == (segNum - 1)){
    S_LB=B_pt(sIdx-1);
    S_UB=R_PosInf;
  }else{
    S_LB=B_pt(sIdx-1);
    S_UB=B_pt(sIdx);
  }

  Function rtrunN("rtruncnorm");
  // Rcout << "sIdx: " << sIdx << "; S_LB: " << S_LB << "; S_UB: " << S_UB << "; mean: " << S_in_mean << "; sd: " << S_in_sd << std::endl;
  // set_seed(Named("seed")=seed_);
  NumericVector tmp_res = rtrunN(1, Named("a") = S_LB, Named("b") = S_UB, Named("mean") = S_in_mean, Named("sd") = S_in_sd);
  double res = tmp_res[0];
  return res;
}

// [[Rcpp::export()]]
arma::mat beta_df(arma::mat u_df, arma::mat coef_df, double thre){
  arma::mat theta_df = u_df * coef_df.t();
  arma::mat res = theta_df % (abs(theta_df)>thre);
  return res;
}

// [[Rcpp::export()]]
arma::mat Gibbs_l_etaC(arma::vec y_, arma::mat x_, arma::mat u_, arma::vec d0_, double t_, 
                       arma::mat l_eta_, arma::mat l_xi_, arma::mat l_Z_,  arma::mat l_nu_){
  int n = x_.n_rows;
  int p = x_.n_cols;
  int q = u_.n_cols;
  // Rcout << "n: " << n << "; p: " << p << "; q:" <<q <<std::endl;
  
  arma::mat l_s_ = l_Z_ % l_nu_;
  arma::mat l_eta_old = l_eta_;
  arma::mat l_alpha_old = l_xi_ % l_eta_;
  arma::mat beta_old = beta_df(u_, l_alpha_old, t_);
  
  // int k=0, j=9;
  for(int k=0; k<q; k++){
    arma::uvec pos_u_idx = u_.col(k)>=0, neg_u_idx = u_.col(k)<0;
    arma::uvec pos_u_nIdx = find(pos_u_idx==1), neg_u_nIdx = find(neg_u_idx==1);
    for(int j=0; j<p; j++){
      // if(k==0 & (j<=3)){
      //   Rcout << "k: " << k << "; j: " << j <<std::endl;
      //   Rcout << l_alpha_old << std::endl;
      // }
      IntegerVector min_kTmp =  seq_len(q), min_jTmp =  seq_len(p);
      min_kTmp.erase(k);
      arma::uvec min_k = as<arma::uvec>(min_kTmp) - 1;
      min_jTmp.erase(j);
      arma::uvec min_j = as<arma::uvec>(min_jTmp) - 1;
      
      arma::uvec jVec(1);
      jVec.fill(j);
      
      arma::vec tmp_term = u_.cols(min_k) * l_alpha_old.submat(jVec,min_k).t();
      arma::vec a1_vec(n), a2_vec(n), T1_vec(n), T2_vec(n);
      
      a1_vec.fill(-0.5 * l_xi_(j,k) * l_xi_(j,k)); 
      a1_vec = a1_vec % u_.col(k) % u_.col(k) % x_.col(j) % x_.col(j) / d0_;
      
      a2_vec.fill(-1 * l_xi_(j,k));
      a2_vec = a2_vec % u_.col(k) / d0_;
      a2_vec = a2_vec % (x_.col(j) % x_.col(j) % tmp_term  + x_.col(j) % (sum(beta_old.cols(min_j) % x_.cols(min_j),1) - y_) );
      
      T1_vec.fill(t_);
      T1_vec = T1_vec - tmp_term; 
      T1_vec = T1_vec /(l_xi_(j,k) * u_.col(k));
      
      T2_vec.fill(-t_);
      T2_vec = T2_vec - tmp_term; 
      T2_vec = T2_vec /(l_xi_(j,k) * u_.col(k));
      
      arma::vec L_vec = join_cols(T1_vec.elem(pos_u_nIdx), T2_vec.elem(neg_u_nIdx)),
        U_vec = join_cols(T2_vec.elem(pos_u_nIdx), T1_vec.elem(neg_u_nIdx)) ;
      L_vec = sort(L_vec);
      U_vec = sort(U_vec);
      
      double a1_prior = -0.5/l_s_(j,k), a2_prior=0.0;
      double l_eta_prop=0.0;
      l_eta_prop = Gibbs_trunNormC(L_vec, a1_vec, a2_vec, U_vec, a1_vec, a2_vec, a1_prior, a2_prior);
      // List l_eta_prop = List::create(Named("L_vec") = L_vec, Named("a1_vec") = a1_vec, Named("a2_vec") = a2_vec,
      //                                Named("U_vec") = U_vec, Named("a1_prior") = a1_prior, Named("a2_prior") = a2_prior);
      // IntegerVector l_eta_prop = Gibbs_trunNormC(L_vec, a1_vec, a2_vec, U_vec, a1_vec, a2_vec, a1_prior, a2_prior, s_);
      // return l_eta_prop;
      
      l_eta_old(j,k) = l_eta_prop;
      l_alpha_old(j,k) =l_xi_(j,k) * l_eta_prop;
      beta_old = beta_df(u_, l_alpha_old, t_);
    }
  }
  return l_eta_old;
}




// [[Rcpp::export()]]
double max_GibbsThreC(arma::vec y_, arma::mat x_, arma::mat u_, arma::vec d0_, double sigProp_,
                      arma::mat l_eta_, arma::mat l_xi_){
  int n = x_.n_rows;
  int p = x_.n_cols;
  arma::mat l_alpha_ = l_xi_ % l_eta_;
  arma::mat theta_ = u_ * l_alpha_.t();
  arma::mat thetaAbs = abs(theta_);
  arma::vec thetaVec = reshape(thetaAbs, n*p, 1);
  arma::vec sigPropVec = { 1.0-sigProp_ };
  arma::vec sigLBVec = quantile(thetaVec, sigPropVec);
  double sigLB = sigLBVec(0);
  // Rcout << "LB: " <<sigLB << std::endl;

  arma::vec tmpCndtTheta  = thetaAbs.elem(find(thetaAbs > sigLB));
  arma::vec cndtTheta  = tmpCndtTheta.elem(find_unique(tmpCndtTheta));
  cndtTheta = sort(cndtTheta);
  arma::vec cndtTheta_ = cndtTheta.elem(find(cndtTheta >= pow(10.0, -13.0)));
  
  int nCndt = cndtTheta_.n_elem;
  arma::vec llh(nCndt), regressor_(n);
  arma::mat beta_(n,p);
  for(int cndtIdx=0; cndtIdx < nCndt; cndtIdx++){
    double t_ = cndtTheta_(cndtIdx), llhTmp=0.0;
    beta_ = theta_ % (thetaAbs > t_);
    regressor_ = sum(beta_ % x_,1);
    for(int i=0; i<n; i++){
      llhTmp = llhTmp + R::dnorm(y_(i), regressor_(i), d0_(i) ,1);
    }
    llh(cndtIdx) = llhTmp;
  }
  int maxIdx = llh.index_max();
  double t_LB = cndtTheta_(maxIdx), res=0.0, t_UB=0.0;
  if(maxIdx == (llh.n_elem-1)){
    res=t_LB;
  }else{
    t_UB = cndtTheta_(maxIdx+1);
    // Function set_seed("set.seed");
    // set_seed(Named("seed")=s_);
    res = R::runif(t_LB, t_UB);
  }
  // Rcout << "t_LB: " << t_LB << "; t_UB: " << t_UB <<std::endl;
  return res;
}

// [[Rcpp::export()]]
arma::mat Gibbs_l_xiC(arma::vec y_, arma::mat x_, arma::mat u_, arma::vec d0_, double t_,
                      arma::mat l_eta_, arma::mat l_xi_, arma::mat l_Z_,  arma::mat l_nu_, arma::mat l_m_){
  int n = x_.n_rows;
  int p = x_.n_cols;
  int q = u_.n_cols;
  // Rcout << "n: " << n << "; p: " << p << "; q:" <<q <<std::endl;

  arma::mat l_s_ = l_Z_ % l_nu_;
  arma::mat l_xi_old = l_xi_;
  arma::mat l_alpha_old = l_xi_ % l_eta_;
  arma::mat beta_old = beta_df(u_, l_alpha_old, t_);

  for(int k=0; k<q; k++){
    // Rcout << "k: " << k << std::endl;
    arma::uvec pos_u_idx = u_.col(k)>=0, neg_u_idx = u_.col(k)<0;
    arma::uvec pos_u_nIdx = find(pos_u_idx==1), neg_u_nIdx = find(neg_u_idx==1);
    for(int j=0; j<p; j++){
      // Rcout << "j: " << j << std::endl;
      IntegerVector min_kTmp =  seq_len(q), min_jTmp =  seq_len(p);
      min_kTmp.erase(k);
      arma::uvec min_k = as<arma::uvec>(min_kTmp) - 1;
      min_jTmp.erase(j);
      arma::uvec min_j = as<arma::uvec>(min_jTmp) - 1;

      arma::uvec jVec(1);
      jVec.fill(j);

      arma::vec tmp_term = u_.cols(min_k) * l_alpha_old.submat(jVec,min_k).t();
      arma::vec a1_vec(n), a2_vec(n), T1_vec(n), T2_vec(n);

      a1_vec.fill(-0.5 * l_eta_(j,k) * l_eta_(j,k));
      a1_vec = a1_vec % u_.col(k) % u_.col(k) % x_.col(j) % x_.col(j) / d0_;

      a2_vec.fill(-1 * l_eta_(j,k));
      a2_vec = a2_vec % u_.col(k) / d0_;
      a2_vec = a2_vec % (x_.col(j) % x_.col(j) % tmp_term  + x_.col(j) % (sum(beta_old.cols(min_j) % x_.cols(min_j),1) - y_) );

      T1_vec.fill(t_);
      T1_vec = T1_vec - tmp_term;
      T1_vec = T1_vec /(l_eta_(j,k) * u_.col(k));

      T2_vec.fill(-t_);
      T2_vec = T2_vec - tmp_term;
      T2_vec = T2_vec /(l_eta_(j,k) * u_.col(k));

      arma::vec L_vec = join_cols(T1_vec.elem(pos_u_nIdx), T2_vec.elem(neg_u_nIdx)),
        U_vec = join_cols(T2_vec.elem(pos_u_nIdx), T1_vec.elem(neg_u_nIdx)) ;
      L_vec = sort(L_vec);
      U_vec = sort(U_vec);

      double a1_prior = -0.5, a2_prior = l_m_(j,k);
      double l_xi_prop = 0.0;
      l_xi_prop = Gibbs_trunNormC(L_vec, a1_vec, a2_vec, U_vec, a1_vec, a2_vec, a1_prior, a2_prior);

      l_xi_old(j,k) = l_xi_prop;
      l_alpha_old(j,k) = l_xi_prop * l_eta_(j,k);
      beta_old = beta_df(u_, l_alpha_old, t_);
    }
  }
  return l_xi_old;
}

