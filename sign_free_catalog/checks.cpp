#include "checks.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
namespace catalog {
void closeCheck(DataType a,DataType b,const std::string& label,double tol) {
    if(!std::isfinite(std::abs(a)) || !std::isfinite(std::abs(b)) || std::abs(a-b)>tol*std::max(1.,std::abs(b))) {
        std::cerr<<label<<": "<<a<<" != "<<b<<"\n"; throw std::runtime_error(label);
    }
}
Recorder::Recorder(std::ostream& stream):out(stream) {
 out<<std::setprecision(17)<<"model,kind,index,modes,weight_real,weight_imag,reference_real,reference_imag,relative_error,sign_real,sign_imag,attempts,accepted\n";
}
void Recorder::record(const std::string& id,const std::string& kind,int index,int n,DataType z,DataType ref,int expected,long attempts,long accepted) {
 closeCheck(z,ref,id+" "+kind);
 double err=std::abs(z-ref)/std::max(1.,std::abs(ref)); maxError=std::max(maxError,err);
 double nan=std::numeric_limits<double>::quiet_NaN(); DataType phase(nan,nan);
 if(expected==0) { require(std::abs(z)<1e-7*std::ldexp(1.,n),id+" is not zero"); ++zeros; }
 else {
  require(std::abs(z)>1e-12,id+" unresolved zero weight"); phase=z/std::abs(z);
  if(expected==1) { closeCheck(phase,1.,id+" positive phase");++positive;maxPhaseError=std::max(maxPhaseError,std::abs(phase-1.)); }
  if(expected==-1) { closeCheck(phase,-1.,id+" negative control");++negative; }
 }
 ++rows; out<<id<<","<<kind<<","<<index<<","<<n<<","<<z.real()<<","<<z.imag()<<","<<ref.real()<<","<<ref.imag()<<","<<err
   <<","<<phase.real()<<","<<phase.imag()<<","<<attempts<<","<<accepted<<"\n";
}
static MatType eye(int n) { return MatType::Identity(n,n); }
static MatType embed(const MatType& a,int n,int offset) { MatType b=MatType::Zero(n,n); b.block(offset,offset,a.rows(),a.cols())=a;return b; }
static DataType adjointReference(int d,int slices,double beta,const std::vector<double>& f) {
 MatType x=MatType::Zero(d,d),g=eye(d);
 for(int i=0;i<d;++i) x(i,i)=.08*(i-.5*(d-1));
 for(int i=0;i<d-1;++i) x(i,i+1)=x(i+1,i)=.3;
 double dt=beta/slices; int pos=0;
 for(int l=0;l<slices;++l) {
  g=(matrixExp(.5*dt*x)*g).eval();
  for(int k=1;k<d;++k) {
   MatType a=MatType::Zero(d,d); double c=1./std::sqrt(double(k*(k+1)));
   for(int i=0;i<k;++i) a(i,i)=c; a(k,k)=-k*c;
   g=(matrixExp(std::sqrt(.2*dt)*f[pos++]*a)*g).eval();
  }
  g=(matrixExp(.5*dt*x)*g).eval();
 }
 Eigen::ComplexEigenSolver<MatType> es(g); require(es.info()==Eigen::Success,"Adjoint reference eigenvalues");
 auto ev=es.eigenvalues(); DataType z=std::ldexp(1.,(d-1)/2);
 for(int i=0;i<d;++i) for(int j=i+1;j<d;++j) z*=ev(i)+ev(j);
 return z;
}
DataType independentTrace(const Hamiltonian& h,int slices,double beta,const std::vector<double>& f) {
 auto p=h.pulses(slices,beta,f);
 if(h.id.rfind("adjoint_d",0)==0) return adjointReference(std::stoi(h.id.substr(9)),slices,beta,f);
 if(h.id.rfind("pairing_lattice",0)==0) {
  DataType z=1.;
  for(int r=0;r<h.modes/4;++r) {std::vector<MatType> local;for(auto& a:p) local.push_back(a.block(8*r,8*r,8,8)); z*=exactTrace(local);}
  return z;
 }
 if(h.conserving) return determinantTrace(p);
 return exactTrace(p);
}
static void checkGreens(PfQMC& q,const std::vector<MatType>& p) {
 MatType b=eye(q.nDim);for(auto& a:p) b=(matrixExp(a)*b).eval();
 MatType g=2.*(eye(q.nDim)+b).inverse();
 require((g-q.g).norm()<2e-7*std::max(1.,g.norm()),"Sweep Green matrix disagrees with independent inverse");
}
void verifyModel(Recorder& r,rdGenerator& rng,const std::string& name,int samples,int sweeps) {
 auto h=makeModel(name); const int slices=4; const double beta=.9;
 for(int t=0;t<samples;++t) {
  std::vector<double> f(h.fieldCount(slices)); double scale=(t%3==0 ? 2.:t%3==1 ? .5:1.);
  for(auto& v:f) v=scale*rng.rdNormal();
  auto p=h.pulses(slices,beta,f); DataType z=traceProduct(p),ref=independentTrace(h,slices,beta,f);
  if(h.modes<=6) closeCheck(z,exactTrace(p),name+" Fock oracle");
  r.record(name,"random_fields",t,h.modes,z,ref);
 }
 Walker w(h,slices,beta,&rng); PfQMC q(&w,5);
 for(int t=0;t<=sweeps;++t) {
  auto p=h.pulses(slices,beta,w.fields()); auto z=traceProduct(p);
  closeCheck(q.sign,z/std::abs(z),name+" tracked sign");closeCheck(q.getSignRaw(),q.sign,name+" raw sign");
  checkGreens(q,p);
  r.record(name,"pfqmc_sweep",t,h.modes,z,independentTrace(h,slices,beta,w.fields()),1,w.attempts(),w.accepted());
  if(t<sweeps) {if(t%2==0) q.rightSweep();else q.leftSweep();}
 }
 require(w.attempts()==sweeps*h.fieldCount(slices) && w.accepted()>0,name+" did not sample fields");
}
static std::vector<MatType> tenPulses() {
 int n=6; double a=std::log(2.); MatType x=2.*a*embed(spinX(2),n,0),z=2.*a*embed(spinZ(2),n,0);
 MatType da=a*embed(eye(3),n,0),db=a*embed(eye(3),n,3),t=MatType::Zero(n,n);
 t.topRightCorner(3,3)=t.bottomLeftCorner(3,3)=a*eye(3);
 std::vector<MatType> p={t,-x,db,-z,t,x,z,da,z,-x};std::reverse(p.begin(),p.end());
 for(auto& m:p) m=numberMajorana(m);return p;
}
void verifyCore(Recorder& r,rdGenerator& rng) {
 // Unrestricted complex Gaussian traces exercise orientation and phases.
 for(int t=0;t<8;++t) {
  std::vector<MatType> p;
  for(int l=0;l<5;++l) {MatType a(6,6);for(int i=0;i<6;++i) for(int j=0;j<6;++j) a(i,j)=.17*DataType(rng.rdNormal(),rng.rdNormal());a=(a-a.transpose()).eval();p.push_back(a);}
  auto z=exactTrace(p);closeCheck(traceProduct(p,true),z,"complex spacetime");r.record("complex_oracle","control",t,3,traceProduct(p),z,2);
 }
 MatType a=numberMajorana(eye(1)); std::vector<MatType> lift(4,.5i*std::acos(-1.)*a);
 r.record("spin_lift_minus_identity","negative",0,1,traceProduct(lift),-2.,-1);
 auto spin=numberMajorana(spinZ(2)); std::vector<MatType> complexSpin={DataType(std::log(4.),std::acos(-1.))*spin};
 r.record("complex_spin1","negative",0,3,traceProduct(complexSpin),-4.5,-1);
 auto ten=tenPulses();double negative=-7772907825./67108864.;
 closeCheck(exactTrace(ten),negative,"ten pulse Fock");r.record("crosscell_ten_pulses","negative",0,6,traceProduct(ten),negative,-1);
 // Test the production sign tracker with negative total weight and real Hermitian slice traces.
 Spinless_tV w;w.nDim=12;for(auto& p:ten) w.op_array.push_back(new FieldOperator(p,&rng,false));
 MatType direction=.1*numberMajorana(embed(spinZ(2),6,3));
 w.op_array.push_back(new FieldOperator(direction,&rng,true));PfQMC q(&w,3);
 for(int t=0;t<6;++t) {
  auto p=ten;p.push_back(static_cast<FieldOperator*>(w.op_array.back())->value()*direction);
  auto ref=exactTrace(p);closeCheck(q.sign,ref/std::abs(ref),"negative model update sign");checkGreens(q,p);
  r.record("negative_walker","negative",t,6,traceProduct(p),ref,ref.real()<0 ? -1:1);
  if(t%2==0)q.rightSweep();else q.leftSweep();
 }
 double b=2.*std::asinh(std::pow(.5,.25));
 std::vector<MatType> zero={b*numberMajorana(spinX(2)),b*spin,-b*numberMajorana(spinX(2)),-b*spin};
 closeCheck(exactTrace(zero),0.,"spin1 exact zero");r.record("spin1_trace_zero","zero",0,3,traceProduct(zero,true),0.,0);
 // Explicit HS identity, including the centering constant in the interacting Hamiltonian.
 auto h=fourMode();int dim=1<<h.modes;MatType jac=MatType::Zero(24,24);
 for(int i=1;i<24;++i)jac(i-1,i)=jac(i,i-1)=std::sqrt(double(i));
 Eigen::SelfAdjointEigenSolver<MatType> es(jac);
 for(auto& c:h.channels) {
  auto o=fockMatrix(c.a);MatType integral=MatType::Zero(dim,dim);
  for(int k=0;k<24;++k) integral+=std::norm(es.eigenvectors()(0,k))*matrixExp(std::sqrt(.2*c.coupling)*es.eigenvalues()(k)*o);
  require((integral-matrixExp(.1*c.coupling*o*o)).norm()<1e-10,"Gaussian HS normalization");
 }
 auto kinetic=fockMatrix(h.q0);auto o1=fockMatrix(h.channels[0].a),o2=fockMatrix(h.channels[1].a);
 MatType direct=-kinetic-.35*(o1*o1-.5*eye(dim))-.55*(o2*o2-.5*eye(dim));
 require((direct-h.exactHamiltonian()).norm()<1e-12,"Four-mode density centering constant");
 // Fixed non-Gaussian unitary: compare conjugated physical Hamiltonian and every pulse.
 auto u=transformedHopping(4),base=openHopping(4);MatType U=eye(16);
 for(auto gate:u.controlledPhases){int i,j;double theta;std::tie(i,j,theta)=gate;for(int v=0;v<16;++v)if((v&(1<<i))&&(v&(1<<j)))U(v,v)*=std::exp(1.i*theta);}
 require((u.exactHamiltonian()-U*base.exactHamiltonian()*U.adjoint()).norm()<1e-12,"Local unitary Hamiltonian");
 std::vector<double> f(u.fieldCount(2));for(auto& x:f)x=rng.rdNormal();auto pulses=u.pulses(2,.5,f);MatType product=eye(16);
 for(auto& p:pulses) product=(matrixExp(U*fockMatrix(p)*U.adjoint())*product).eval();
 r.record("local_unitary_physical","Fock_conjugation",0,4,traceProduct(pulses),product.trace());
 // Reject invalid representation parameters instead of silently constructing them.
 auto rejects=[](const std::function<void()>& f){try{f();}catch(const std::runtime_error&){return;}throw std::runtime_error("Invalid model accepted");};
 rejects([]{exteriorSp(3,3);});rejects([]{exteriorSL(5,2);});rejects([]{spinCell(1);});rejects([]{sl2Orthogonal({0,2});});rejects([]{oddPairRing(5,2);});
 rejects([]{auto h=fourMode();h.channels[0].a=MatType::Zero(2,2);h.validate();});
 // Check independence-Metropolis decisions and the local Green update against Fock weights.
 rdGenerator actual(83910),oracle(83910);MatType dir=.6*fourMode().channels[0].a;
 FieldOperator op(dir,&actual,true);double phi=oracle.rdNormal();auto baseP=fourMode().pulses(2,.6,std::vector<double>(4,.4));
 int acceptCount=0,rejectCount=0;
 for(int t=0;t<64;++t) {
  double proposal=oracle.rdNormal(),uniform=oracle.rdUniform01();
  auto oldP=baseP,newP=baseP;oldP.insert(oldP.begin(),phi*dir);newP.insert(newP.begin(),proposal*dir);
  auto oldZ=exactTrace(oldP),newZ=exactTrace(newP);bool accepted=uniform<std::min(1.,std::abs(newZ/oldZ));
  MatType bmat=eye(8);for(auto& p:oldP)bmat=(matrixExp(p)*bmat).eval();MatType green=2.*(eye(8)+bmat).inverse();
  auto phase=op.update(green);if(accepted){phi=proposal;++acceptCount;}else ++rejectCount;
  closeCheck(op.value(),phi,"Metropolis decision");closeCheck(phase,accepted ? newZ/oldZ/std::abs(newZ/oldZ):DataType(1.),"Metropolis ratio phase");
  auto checkP=baseP;checkP.insert(checkP.begin(),phi*dir);bmat=eye(8);for(auto& p:checkP)bmat=(matrixExp(p)*bmat).eval();
  require((green-2.*(eye(8)+bmat).inverse()).norm()<1e-9,"Local Green update");
 }
 require(acceptCount>0 && rejectCount>0,"Metropolis test needs both outcomes");
 // Strong commutators have negative auxiliary eigenvalues. These are not tests
 // confined to a small neighborhood of the identity.
 for(int d:{4,6,8}) {
  MatType x=MatType::Zero(d,d),z=x;x(0,1)=x(1,0)=1.03;z(0,0)=1.03;z(1,1)=-1.03;
  MatType a=numberMajorana(wedgeAction(x,d/2)),b=numberMajorana(wedgeAction(z,d/2));std::vector<MatType> p={a,b,-a,-b};
  r.record("sl_exterior_"+std::to_string(d)+"_"+std::to_string(d/2),"strong_commutator",0,a.rows()/2,traceProduct(p),determinantTrace(p));
 }
 for(auto nk:std::vector<std::pair<int,int>>{{2,2},{3,3},{4,3},{5,5}}) {
  int n=nk.first,k=nk.second;MatType x=MatType::Zero(2*n,2*n),z=x;
  x(0,n)=x(n,0)=1.03;z(0,0)=1.03;z(n,n)=-1.03;auto p=primitiveBasis(n,k);
  MatType wx=wedgeAction(x,k),wz=wedgeAction(z,k);
  require(((eye(p.rows())-p*p.adjoint())*wx*p).norm()<1e-11,"Primitive subspace invariance X");
  require(((eye(p.rows())-p*p.adjoint())*wz*p).norm()<1e-11,"Primitive subspace invariance Z");
  MatType a=numberMajorana(p.adjoint()*wx*p),b=numberMajorana(p.adjoint()*wz*p);
  std::vector<MatType> pulses={a,b,-a,-b};
  r.record("sp_exterior_"+std::to_string(n)+"_"+std::to_string(k),"strong_commutator",0,p.cols(),traceProduct(pulses),determinantTrace(pulses),n==3 ? -1:1);
 }
}
}
