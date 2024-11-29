// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MartG_TestStatC
double MartG_TestStatC(int n, int k, Eigen::MatrixXd X, Eigen::MatrixXd Xj);
RcppExport SEXP _HDTSA_MartG_TestStatC(SEXP nSEXP, SEXP kSEXP, SEXP XSEXP, SEXP XjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Xj(XjSEXP);
    rcpp_result_gen = Rcpp::wrap(MartG_TestStatC(n, k, X, Xj));
    return rcpp_result_gen;
END_RCPP
}
// MartG_ftC
Eigen::MatrixXd MartG_ftC(int n, int k, int p, int d, Eigen::MatrixXd X, Eigen::MatrixXd Xj);
RcppExport SEXP _HDTSA_MartG_ftC(SEXP nSEXP, SEXP kSEXP, SEXP pSEXP, SEXP dSEXP, SEXP XSEXP, SEXP XjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Xj(XjSEXP);
    rcpp_result_gen = Rcpp::wrap(MartG_ftC(n, k, p, d, X, Xj));
    return rcpp_result_gen;
END_RCPP
}
// MartG_bootc
std::vector<double> MartG_bootc(const int n, const int k, const int p, const int d, const int B, double bn, int method, Eigen::MatrixXd ft, Eigen::MatrixXd Xi_temp);
RcppExport SEXP _HDTSA_MartG_bootc(SEXP nSEXP, SEXP kSEXP, SEXP pSEXP, SEXP dSEXP, SEXP BSEXP, SEXP bnSEXP, SEXP methodSEXP, SEXP ftSEXP, SEXP Xi_tempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type bn(bnSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ft(ftSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Xi_temp(Xi_tempSEXP);
    rcpp_result_gen = Rcpp::wrap(MartG_bootc(n, k, p, d, B, bn, method, ft, Xi_temp));
    return rcpp_result_gen;
END_RCPP
}
// sigmak
Eigen::MatrixXd sigmak(Eigen::MatrixXd Y, Eigen::MatrixXd Y_mean, int k, int n);
RcppExport SEXP _HDTSA_sigmak(SEXP YSEXP, SEXP Y_meanSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y_mean(Y_meanSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmak(Y, Y_mean, k, n));
    return rcpp_result_gen;
END_RCPP
}
// thresh_C
Eigen::MatrixXd thresh_C(Eigen::MatrixXd mat, double delta);
RcppExport SEXP _HDTSA_thresh_C(SEXP matSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(thresh_C(mat, delta));
    return rcpp_result_gen;
END_RCPP
}
// MatMult
SEXP MatMult(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _HDTSA_MatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// WN_teststatC
Rcpp::List WN_teststatC(Eigen::MatrixXd X, int n, int p, int k);
RcppExport SEXP _HDTSA_WN_teststatC(SEXP XSEXP, SEXP nSEXP, SEXP pSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(WN_teststatC(X, n, p, k));
    return rcpp_result_gen;
END_RCPP
}
// WN_ftC
Eigen::MatrixXd WN_ftC(int n, int k, int p, Eigen::MatrixXd X, Eigen::MatrixXd X_mean);
RcppExport SEXP _HDTSA_WN_ftC(SEXP nSEXP, SEXP kSEXP, SEXP pSEXP, SEXP XSEXP, SEXP X_meanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X_mean(X_meanSEXP);
    rcpp_result_gen = Rcpp::wrap(WN_ftC(n, k, p, X, X_mean));
    return rcpp_result_gen;
END_RCPP
}
// WN_bootc
std::vector<double> WN_bootc(const int n, const int k, const int p, const int B, double bn, int method, Eigen::MatrixXd ft, Eigen::MatrixXd X, Eigen::VectorXd sigma_zero, Eigen::MatrixXd Xi_temp);
RcppExport SEXP _HDTSA_WN_bootc(SEXP nSEXP, SEXP kSEXP, SEXP pSEXP, SEXP BSEXP, SEXP bnSEXP, SEXP methodSEXP, SEXP ftSEXP, SEXP XSEXP, SEXP sigma_zeroSEXP, SEXP Xi_tempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type bn(bnSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ft(ftSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sigma_zero(sigma_zeroSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Xi_temp(Xi_tempSEXP);
    rcpp_result_gen = Rcpp::wrap(WN_bootc(n, k, p, B, bn, method, ft, X, sigma_zero, Xi_temp));
    return rcpp_result_gen;
END_RCPP
}
// bandwith
double bandwith(Eigen::MatrixXd ft, int k, int p, int d, int ken_type);
RcppExport SEXP _HDTSA_bandwith(SEXP ftSEXP, SEXP kSEXP, SEXP pSEXP, SEXP dSEXP, SEXP ken_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ft(ftSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type ken_type(ken_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bandwith(ft, k, p, d, ken_type));
    return rcpp_result_gen;
END_RCPP
}
// TaperQsC
double TaperQsC(double x);
RcppExport SEXP _HDTSA_TaperQsC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(TaperQsC(x));
    return rcpp_result_gen;
END_RCPP
}
// TaperBartC
double TaperBartC(double x);
RcppExport SEXP _HDTSA_TaperBartC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(TaperBartC(x));
    return rcpp_result_gen;
END_RCPP
}
// TaperFlatC
double TaperFlatC(double x, double c);
RcppExport SEXP _HDTSA_TaperFlatC(SEXP xSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(TaperFlatC(x, c));
    return rcpp_result_gen;
END_RCPP
}
// CmpGammaC
Rcpp::List CmpGammaC(Eigen::MatrixXd Vt);
RcppExport SEXP _HDTSA_CmpGammaC(SEXP VtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Vt(VtSEXP);
    rcpp_result_gen = Rcpp::wrap(CmpGammaC(Vt));
    return rcpp_result_gen;
END_RCPP
}
// EvalGammaJC
Eigen::MatrixXd EvalGammaJC(Rcpp::List Gamma, int j, int len);
RcppExport SEXP _HDTSA_EvalGammaJC(SEXP GammaSEXP, SEXP jSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(EvalGammaJC(Gamma, j, len));
    return rcpp_result_gen;
END_RCPP
}
// CmpRhoC
Rcpp::List CmpRhoC(Rcpp::List Gamma, int len);
RcppExport SEXP _HDTSA_CmpRhoC(SEXP GammaSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(CmpRhoC(Gamma, len));
    return rcpp_result_gen;
END_RCPP
}
// EvalRhoMC
Eigen::MatrixXd EvalRhoMC(Rcpp::List Rho, int m, int len);
RcppExport SEXP _HDTSA_EvalRhoMC(SEXP RhoSEXP, SEXP mSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(EvalRhoMC(Rho, m, len));
    return rcpp_result_gen;
END_RCPP
}
// CmpHatSC
Eigen::MatrixXd CmpHatSC(Rcpp::List Rho, double C0, int KT, double cef, int p, int len);
RcppExport SEXP _HDTSA_CmpHatSC(SEXP RhoSEXP, SEXP C0SEXP, SEXP KTSEXP, SEXP cefSEXP, SEXP pSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< double >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< int >::type KT(KTSEXP);
    Rcpp::traits::input_parameter< double >::type cef(cefSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(CmpHatSC(Rho, C0, KT, cef, p, len));
    return rcpp_result_gen;
END_RCPP
}
// SpecEstC
Rcpp::List SpecEstC(Rcpp::List Gamma, int n, int p, int r, int K, Eigen::MatrixXd cross_indices, Eigen::VectorXd J_set, double l_band, double flag_c);
RcppExport SEXP _HDTSA_SpecEstC(SEXP GammaSEXP, SEXP nSEXP, SEXP pSEXP, SEXP rSEXP, SEXP KSEXP, SEXP cross_indicesSEXP, SEXP J_setSEXP, SEXP l_bandSEXP, SEXP flag_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type cross_indices(cross_indicesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type J_set(J_setSEXP);
    Rcpp::traits::input_parameter< double >::type l_band(l_bandSEXP);
    Rcpp::traits::input_parameter< double >::type flag_c(flag_cSEXP);
    rcpp_result_gen = Rcpp::wrap(SpecEstC(Gamma, n, p, r, K, cross_indices, J_set, l_band, flag_c));
    return rcpp_result_gen;
END_RCPP
}
// TestStatC
double TestStatC(Rcpp::List Gamma, int n, int p, int r, int K, Eigen::MatrixXd cross_indices, Eigen::VectorXd J_set, double l_band, double flag_c);
RcppExport SEXP _HDTSA_TestStatC(SEXP GammaSEXP, SEXP nSEXP, SEXP pSEXP, SEXP rSEXP, SEXP KSEXP, SEXP cross_indicesSEXP, SEXP J_setSEXP, SEXP l_bandSEXP, SEXP flag_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type cross_indices(cross_indicesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type J_set(J_setSEXP);
    Rcpp::traits::input_parameter< double >::type l_band(l_bandSEXP);
    Rcpp::traits::input_parameter< double >::type flag_c(flag_cSEXP);
    rcpp_result_gen = Rcpp::wrap(TestStatC(Gamma, n, p, r, K, cross_indices, J_set, l_band, flag_c));
    return rcpp_result_gen;
END_RCPP
}
// CEst2C
Rcpp::List CEst2C(Eigen::MatrixXd x, Rcpp::List Gamma, int n_tilde, int n, int p, int r, Eigen::MatrixXd cross_indices, int l_band);
RcppExport SEXP _HDTSA_CEst2C(SEXP xSEXP, SEXP GammaSEXP, SEXP n_tildeSEXP, SEXP nSEXP, SEXP pSEXP, SEXP rSEXP, SEXP cross_indicesSEXP, SEXP l_bandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type n_tilde(n_tildeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type cross_indices(cross_indicesSEXP);
    Rcpp::traits::input_parameter< int >::type l_band(l_bandSEXP);
    rcpp_result_gen = Rcpp::wrap(CEst2C(x, Gamma, n_tilde, n, p, r, cross_indices, l_band));
    return rcpp_result_gen;
END_RCPP
}
// CEst3C
Eigen::MatrixXd CEst3C(Eigen::MatrixXd x, Rcpp::List Gamma, int n_tilde, int n, int p, int r, Eigen::MatrixXd cross_indices, int l_band);
RcppExport SEXP _HDTSA_CEst3C(SEXP xSEXP, SEXP GammaSEXP, SEXP n_tildeSEXP, SEXP nSEXP, SEXP pSEXP, SEXP rSEXP, SEXP cross_indicesSEXP, SEXP l_bandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type n_tilde(n_tildeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type cross_indices(cross_indicesSEXP);
    Rcpp::traits::input_parameter< int >::type l_band(l_bandSEXP);
    rcpp_result_gen = Rcpp::wrap(CEst3C(x, Gamma, n_tilde, n, p, r, cross_indices, l_band));
    return rcpp_result_gen;
END_RCPP
}
// BandEstC
double BandEstC(Eigen::MatrixXd Chat, int n_tilde, int r, int l_band, int type);
RcppExport SEXP _HDTSA_BandEstC(SEXP ChatSEXP, SEXP n_tildeSEXP, SEXP rSEXP, SEXP l_bandSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Chat(ChatSEXP);
    Rcpp::traits::input_parameter< int >::type n_tilde(n_tildeSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type l_band(l_bandSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(BandEstC(Chat, n_tilde, r, l_band, type));
    return rcpp_result_gen;
END_RCPP
}
// etaC
Eigen::MatrixXd etaC(int n, int p, int B, int n_tilde, double bn, int type);
RcppExport SEXP _HDTSA_etaC(SEXP nSEXP, SEXP pSEXP, SEXP BSEXP, SEXP n_tildeSEXP, SEXP bnSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n_tilde(n_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type bn(bnSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(etaC(n, p, B, n_tilde, bn, type));
    return rcpp_result_gen;
END_RCPP
}
// LongCovEstC
Eigen::MatrixXd LongCovEstC(int n_tilde, int ln, int r, Eigen::VectorXi Shat_c, Eigen::MatrixXd Chat, int Kern);
RcppExport SEXP _HDTSA_LongCovEstC(SEXP n_tildeSEXP, SEXP lnSEXP, SEXP rSEXP, SEXP Shat_cSEXP, SEXP ChatSEXP, SEXP KernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_tilde(n_tildeSEXP);
    Rcpp::traits::input_parameter< int >::type ln(lnSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type Shat_c(Shat_cSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Chat(ChatSEXP);
    Rcpp::traits::input_parameter< int >::type Kern(KernSEXP);
    rcpp_result_gen = Rcpp::wrap(LongCovEstC(n_tilde, ln, r, Shat_c, Chat, Kern));
    return rcpp_result_gen;
END_RCPP
}
// TestStarC
Eigen::VectorXd TestStarC(Eigen::MatrixXd x, Rcpp::List GhatC, int n_tilde, int n, int p, int r, int K, double flag_c, Eigen::MatrixXd cross_indices, Eigen::VectorXd J_set, int l_band, int B_monte, int type);
RcppExport SEXP _HDTSA_TestStarC(SEXP xSEXP, SEXP GhatCSEXP, SEXP n_tildeSEXP, SEXP nSEXP, SEXP pSEXP, SEXP rSEXP, SEXP KSEXP, SEXP flag_cSEXP, SEXP cross_indicesSEXP, SEXP J_setSEXP, SEXP l_bandSEXP, SEXP B_monteSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type GhatC(GhatCSEXP);
    Rcpp::traits::input_parameter< int >::type n_tilde(n_tildeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type flag_c(flag_cSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type cross_indices(cross_indicesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type J_set(J_setSEXP);
    Rcpp::traits::input_parameter< int >::type l_band(l_bandSEXP);
    Rcpp::traits::input_parameter< int >::type B_monte(B_monteSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(TestStarC(x, GhatC, n_tilde, n, p, r, K, flag_c, cross_indices, J_set, l_band, B_monte, type));
    return rcpp_result_gen;
END_RCPP
}
// minor_P
Eigen::VectorXd minor_P(Eigen::MatrixXd Wr, Eigen::MatrixXd Ws, int d1, int d2);
RcppExport SEXP _HDTSA_minor_P(SEXP WrSEXP, SEXP WsSEXP, SEXP d1SEXP, SEXP d2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Wr(WrSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Ws(WsSEXP);
    Rcpp::traits::input_parameter< int >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< int >::type d2(d2SEXP);
    rcpp_result_gen = Rcpp::wrap(minor_P(Wr, Ws, d1, d2));
    return rcpp_result_gen;
END_RCPP
}
// Vech2Mat_new
Eigen::MatrixXd Vech2Mat_new(Eigen::VectorXd P, int d);
RcppExport SEXP _HDTSA_Vech2Mat_new(SEXP PSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(Vech2Mat_new(P, d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HDTSA_MartG_TestStatC", (DL_FUNC) &_HDTSA_MartG_TestStatC, 4},
    {"_HDTSA_MartG_ftC", (DL_FUNC) &_HDTSA_MartG_ftC, 6},
    {"_HDTSA_MartG_bootc", (DL_FUNC) &_HDTSA_MartG_bootc, 9},
    {"_HDTSA_sigmak", (DL_FUNC) &_HDTSA_sigmak, 4},
    {"_HDTSA_thresh_C", (DL_FUNC) &_HDTSA_thresh_C, 2},
    {"_HDTSA_MatMult", (DL_FUNC) &_HDTSA_MatMult, 2},
    {"_HDTSA_WN_teststatC", (DL_FUNC) &_HDTSA_WN_teststatC, 4},
    {"_HDTSA_WN_ftC", (DL_FUNC) &_HDTSA_WN_ftC, 5},
    {"_HDTSA_WN_bootc", (DL_FUNC) &_HDTSA_WN_bootc, 10},
    {"_HDTSA_bandwith", (DL_FUNC) &_HDTSA_bandwith, 5},
    {"_HDTSA_TaperQsC", (DL_FUNC) &_HDTSA_TaperQsC, 1},
    {"_HDTSA_TaperBartC", (DL_FUNC) &_HDTSA_TaperBartC, 1},
    {"_HDTSA_TaperFlatC", (DL_FUNC) &_HDTSA_TaperFlatC, 2},
    {"_HDTSA_CmpGammaC", (DL_FUNC) &_HDTSA_CmpGammaC, 1},
    {"_HDTSA_EvalGammaJC", (DL_FUNC) &_HDTSA_EvalGammaJC, 3},
    {"_HDTSA_CmpRhoC", (DL_FUNC) &_HDTSA_CmpRhoC, 2},
    {"_HDTSA_EvalRhoMC", (DL_FUNC) &_HDTSA_EvalRhoMC, 3},
    {"_HDTSA_CmpHatSC", (DL_FUNC) &_HDTSA_CmpHatSC, 6},
    {"_HDTSA_SpecEstC", (DL_FUNC) &_HDTSA_SpecEstC, 9},
    {"_HDTSA_TestStatC", (DL_FUNC) &_HDTSA_TestStatC, 9},
    {"_HDTSA_CEst2C", (DL_FUNC) &_HDTSA_CEst2C, 8},
    {"_HDTSA_CEst3C", (DL_FUNC) &_HDTSA_CEst3C, 8},
    {"_HDTSA_BandEstC", (DL_FUNC) &_HDTSA_BandEstC, 5},
    {"_HDTSA_etaC", (DL_FUNC) &_HDTSA_etaC, 6},
    {"_HDTSA_LongCovEstC", (DL_FUNC) &_HDTSA_LongCovEstC, 6},
    {"_HDTSA_TestStarC", (DL_FUNC) &_HDTSA_TestStarC, 13},
    {"_HDTSA_minor_P", (DL_FUNC) &_HDTSA_minor_P, 4},
    {"_HDTSA_Vech2Mat_new", (DL_FUNC) &_HDTSA_Vech2Mat_new, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_HDTSA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
