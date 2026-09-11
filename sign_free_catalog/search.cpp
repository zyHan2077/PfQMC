#include "checks.h"
#include "catalog/search_models.h"
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace catalog;
static MatType eye(int n){return MatType::Identity(n,n);}
static void quaternionControls(Recorder& r) {
 // Check the representation against direct quaternion congruence, without spectra.
 for(int n:{2,3})for(bool anti:{false,true}) {
  auto basis=quaternionBasis(n,anti),herm=quaternionBasis(n,false);
  MatType k=.19*(herm[0]-herm[1])+.13*herm[n]+.11*herm[n+1];
  MatType l=-.09*(herm[0]-herm[1])+.17*herm[n+2];
  MatType g=matrixExp(l)*matrixExp(k),direct(basis.size(),basis.size());
  for(int b=0;b<int(basis.size());++b)for(int a=0;a<int(basis.size());++a) {
   require(std::abs(.5*(basis[a].adjoint()*basis[b]).trace().real()-(a==b))<1e-12,"Quaternion orthonormal basis");
   direct(a,b)=.5*(basis[a].adjoint()*g*basis[b]*g.adjoint()).trace().real();
  }
  MatType product=matrixExp(quaternionAction(l,anti))*matrixExp(quaternionAction(k,anti));
  require((product-direct).norm()<1e-10,"Quaternion congruence oracle");
  require((quaternionAction(k,anti)*quaternionAction(l,anti)-quaternionAction(l,anti)*quaternionAction(k,anti)).norm()>.01,"Quaternion noncommuting directions");
 }
 // Herm_H(2) is the six-orbital Lorentz family, and does not count twice.
 auto basis=quaternionBasis(2,false);MatType u=eye(6);double c=1./std::sqrt(2.);
 u(0,0)=u(1,0)=u(0,1)=c;u(1,1)=-c;
 for(int a=1;a<6;++a) {
  MatType k=a==1?MatType(.5*(basis[0]-basis[1])):MatType(c*basis[a]);
  MatType star=MatType::Zero(6,6);star(0,a)=star(a,0)=1.;
  require((u.transpose()*quaternionAction(k)*u-star).norm()<1e-12,"Quaternion2/Lorentz6 equivalence");
 }
 MatType star=MatType::Zero(6,6),leaf=star;star(0,1)=star(1,0)=1.;leaf(1,1)=1.;
 MatType a=numberMajorana(u*star*u.transpose()),b=numberMajorana(u*leaf*u.transpose());double t=std::log(3.);
 std::vector<MatType> bad={t*a,2.*t*b,-t*a,-2.*t*b};
 closeCheck(exactTrace(bad),-11200./81.,"Quaternion independent orbital negative Fock");
 r.record("quaternion_independent_orbital","negative",0,6,traceProduct(bad),-11200./81.,-1);
 // Odd real multiplicity excludes a common anti-linear symmetry squaring to -I.
 MatType k=basis[0]-basis[1];
 for(bool anti:{false,true}) {
  MatType h=quaternionAction(k,anti);h+=3.*eye(h.rows());
  Eigen::SelfAdjointEigenSolver<MatType> es(numberMajorana(h));int top=0;
  for(auto value:es.eigenvalues())if(std::abs(value-5.)<1e-10)++top;
  require(top==(anti?3:1),"Quaternion odd real multiplicity certificate");
 }
}
static void verifyCandidate(Recorder& r,rdGenerator& rng,const std::string& id,int samples,int sweeps) {
 auto m=makeSearchModel(id);auto& h=m.h;int slices=3;double beta=.8;
 for(int t=0;t<samples;++t){std::vector<double> f(h.fieldCount(slices));for(auto& x:f)x=(t%3==0 ? 1.5:t%3==1 ? .5:1.)*rng.rdNormal();
  auto p=h.pulses(slices,beta,f);auto z=traceProduct(p);closeCheck(z,determinantTrace(p),id+" independent determinant");
  if(h.modes<=6)closeCheck(z,exactTrace(p),id+" full Fock");
  r.record(id,"search_fields",t,h.modes,z,m.reference(slices,beta,f));
 }
 Walker walker(h,slices,beta,&rng);PfQMC q(&walker,5);
 for(int t=0;t<=sweeps;++t){auto f=walker.fields();auto p=h.pulses(slices,beta,f);auto z=traceProduct(p);
  closeCheck(q.sign,z/std::abs(z),id+" tracked sign");closeCheck(q.getSignRaw(),q.sign,id+" raw sign");
  MatType b=eye(q.nDim);for(auto& a:p)b=(matrixExp(a)*b).eval();MatType g=2.*(eye(q.nDim)+b).inverse();
  require((q.g-g).norm()<2e-7*std::max(1.,g.norm()),"Candidate sweep Green mismatch");
  r.record(id,"search_sweep",t,h.modes,z,m.reference(slices,beta,f),1,walker.attempts(),walker.accepted());
  if(t<sweeps){if(t%2==0)q.rightSweep();else q.leftSweep();}
 }
 require(walker.attempts()==sweeps*h.fieldCount(slices) && walker.accepted()>0,"Candidate walker never moved");
 std::cerr<<"PASS "<<id<<"\n";
}
static void controls(Recorder& r) {
 quaternionControls(r);
 // Traceless Lorentz Sym2: independent congruence and an exact rational negative control.
 auto lb=lorentzSymmetricBasis(5);MatType lk=MatType::Zero(5,5);lk(0,1)=lk(1,0)=1.;
 MatType ll=MatType::Zero(5,5);ll(0,2)=ll(2,0)=.17;
 MatType lg=matrixExp(ll)*matrixExp(.23*lk),ld(14,14);
 for(int a=0;a<14;++a)for(int b=0;b<14;++b)ld(a,b)=(lb[a]*lg*lb[b]*lg.transpose()).trace();
 require((ld-matrixExp(lorentzSymmetricAction(ll))*matrixExp(.23*lorentzSymmetricAction(lk))).norm()<1e-10,"Lorentz Sym2 congruence oracle");
 MatType lp=MatType::Zero(5,5);lp(0,0)=1.;
 MatType la=std::log(4.)*numberMajorana(lorentzSymmetricAction(lk));
 MatType lc=5.*std::log(1.25)*numberMajorana(lorentzSymmetricAction(lp));
 std::vector<MatType> lbad={la,lc,-la,-lc};double lref=-20.65029352655228;
 closeCheck(determinantTrace(lbad),lref,"Lorentz Sym2 exact rational negative control");
 r.record("lorentz_sym_projected_independent_center","negative",0,14,traceProduct(lbad),lref,-1);
 double a=std::log(3.);MatType b=MatType::Zero(4,4),p=b;b(0,1)=b(1,0)=1.;p(1,1)=1.;
 std::vector<MatType> bad={a*numberMajorana(b),2.*a*numberMajorana(p),-a*numberMajorana(b),-2.*a*numberMajorana(p)};
 closeCheck(exactTrace(bad),-2800./81.,"Lorentz broken-field exact negative");
 r.record("lorentz_independent_leaf","negative",0,4,traceProduct(bad),-2800./81.,-1);
 // The two-mode tensor and four-orbital Lorentz star are the SAME family.
 MatType sx=MatType::Zero(2,2),sy=sx,sz=sx;sx(0,1)=sx(1,0)=1.;sy(0,1)=-1.i;sy(1,0)=1.i;sz(0,0)=1.;sz(1,1)=-1.;
 std::vector<MatType> pauli={eye(2),sx,sy,sz};MatType u(4,4);
 for(int col=0;col<4;++col)for(int i=0;i<2;++i)for(int j=0;j<2;++j)u(2*i+j,col)=pauli[col](i,j)/std::sqrt(2.);
 require((u.adjoint()*u-eye(4)).norm()<1e-12,"Pauli basis unitarity");
 for(int i=1;i<4;++i){MatType star=MatType::Zero(4,4);star(0,i)=star(i,0)=1.;require((u.adjoint()*conjugateTensorAction(.5*pauli[i])*u-star).norm()<1e-12,"Tensor2/Lorentz4 equivalence");}
 // Uncorrelated tensor factors already give a complex weight: not sign=1.
 std::vector<MatType> broken;
 for(auto& k:std::vector<MatType>{.4*sx,.4*sy,.4*sz}){MatType h=MatType::Zero(4,4);for(int i=0;i<2;++i)for(int j=0;j<2;++j)h.block(2*i,2*j,2,2)=k(i,j)*eye(2);broken.push_back(numberMajorana(h));}
 auto phase=traceProduct(broken);require(std::abs(phase.imag())>.1,"Missing tensor correlation did not trigger phase control");
 r.record("tensor_unshared_fields","complex_control",0,4,phase,exactTrace(broken),2);
 // Exact zero at a real boost word, recovered by the global signed Pfaffian.
 double z=std::asinh(std::pow(.5,.25));std::vector<MatType> zero={numberMajorana(conjugateTensorAction(z*sx)),numberMajorana(conjugateTensorAction(z*sz)),numberMajorana(conjugateTensorAction(-z*sx)),numberMajorana(conjugateTensorAction(-z*sz))};
 closeCheck(exactTrace(zero),0.,"Lorentz/tensor zero");r.record("tensor_lorentz_zero","zero",0,4,traceProduct(zero,true),0.,0);
 // Independent orbital filters violate the Sym^2 shared-Cartan restriction.
 MatType x1=MatType::Zero(3,3),x2=x1;x1(0,1)=x1(1,0)=1.;x2(1,2)=x2(2,1)=1.;
 MatType kinetic=.5*(symmetricPowerAction(x1,2)+symmetricPowerAction(x2,2));
 int diag[6][6]={{0,1,4,0,-4,2},{-1,1,-2,-3,-4,2},{3,-2,-2,2,1,1},{-2,-4,-1,2,0,3},{-1,2,-4,0,-2,-2},{-1,-4,3,2,3,3}};
 std::vector<MatType> filters;for(int l=0;l<6;++l){MatType h=kinetic;for(int i=0;i<6;++i)h(i,i)=diag[l][i];filters.push_back(numberMajorana(h));}
 double oracle=-71463587094.3984974;closeCheck(exactTrace(filters),oracle,"Sym2 independent density Fock negative");
 r.record("symmetric_power_independent_orbitals","negative",0,6,traceProduct(filters),oracle,-1);
 // The spectral obstruction used in the Wei argument is present in the actual families.
 for(auto h:std::vector<MatType>{2.*eye(4)+b,3.*eye(4)+conjugateTensorAction(sz)}){
  Eigen::SelfAdjointEigenSolver<MatType> es(numberMajorana(h));auto ev=es.eigenvalues();require(ev(ev.size()-1)-ev(ev.size()-2)>.5,"Wei simple real eigenvalue certificate");}
 MatType cartan=MatType::Zero(3,3);cartan.diagonal()<<1.,3.,9.;Eigen::SelfAdjointEigenSolver<MatType> es(numberMajorana(symmetricPowerAction(cartan,2)));
 auto ev=es.eigenvalues();for(int i=1;i<ev.size();++i)require(ev(i)-ev(i-1)>1.,"Sym2 Wei simple spectrum certificate");
}
int main(int argc,char** argv) {
 try {
  mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=32,sweeps=8,seed=9112026;std::string id="all",csv;
  for(int i=1;i<argc;++i){std::string key=argv[i];require(i+1<argc,"Missing option value");std::string val=argv[++i];if(key=="--model")id=val;else if(key=="--samples")samples=std::stoi(val);else if(key=="--sweeps")sweeps=std::stoi(val);else if(key=="--seed")seed=std::stoi(val);else if(key=="--csv")csv=val;else throw std::runtime_error("Unknown option "+key);}
  require(samples>0 && sweeps>=2,"Need samples>0 and both sweep directions");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open CSV");}
  Recorder r(csv.empty()?std::cout:file);rdGenerator rng(seed);if(id=="all"||id=="controls")controls(r);
  auto names=id=="all"?searchModelNames():id=="controls"?std::vector<std::string>{}:std::vector<std::string>{id};
  for(auto& name:names)verifyCandidate(r,rng,name,samples,sweeps);
  std::cerr<<"PASS rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zeros<<" negative="<<r.negative<<" max_relative_error="<<r.maxError<<" max_phase_error="<<r.maxPhaseError<<" seed="<<seed<<"\n";
 }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}return 0;
}
