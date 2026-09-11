#include "checks.h"
#include <algorithm>
namespace catalog {
static MatType eye(int n){return MatType::Identity(n,n);}
static DataType staticReference(const StaticModel& m,const std::vector<int>& s,double beta) {
 auto a=m.generator(s);
 if(m.pairing!=0.) {
  Eigen::SelfAdjointEigenSolver<MatType> es(a);DataType z=std::exp(m.logScalar(s,beta));
  for(int i=0;i<a.rows();++i) z*=std::sqrt(2.*std::cosh(.5*beta*es.eigenvalues()(i)));
  // This positive thermal spectrum formula is ONLY an oracle for one Hermitian H.
  return z;
 }
 MatType h=-numberMatrix(a);Eigen::SelfAdjointEigenSolver<MatType> es(h);double z=std::exp(m.logScalar(s,beta));
 for(int i=0;i<h.rows();++i)z*=2.*std::cosh(.5*beta*es.eigenvalues()(i));return z;
}
void verifySpecials(Recorder& r,rdGenerator& rng,int samples) {
 const double beta=.9;
 for(int pairing=0;pairing<2;++pairing) {
  StaticModel m;m.pairing=pairing ? .17:0.;std::string name=pairing ? "static_bdg":"static_fk";
  const int n=m.lx*m.ly;
  // Enumerate every static background on a genuine 3x3 graph (no HS time histories).
  for(unsigned v=0;v<(1u<<n);++v) {
   std::vector<int> s(n);for(int i=0;i<n;++i)s[i]=(v&(1u<<i)) ? 1:-1;
   auto z=std::exp(m.logScalar(s,beta))*traceProduct({beta*m.generator(s)});
   r.record(name,"all_static_backgrounds",v,n,z,staticReference(m,s,beta));
  }
  std::vector<int> s(n,1);auto weight=[&](const std::vector<int>& a){return std::exp(m.logScalar(a,beta))*traceProduct({beta*m.generator(a)});};
  auto z=weight(s);int accepted=0;
  for(int t=0;t<samples*n;++t) {
   auto proposal=s;proposal[t%n]*=-1;auto next=weight(proposal);
   if(rng.rdUniform01()<std::min(1.,std::abs(next/z))){s=proposal;z=next;++accepted;}
   r.record(name,"static_metropolis",t,n,z,staticReference(m,s,beta),1,t+1,accepted);
  }
  require(accepted>0,"Static background sampler did not move");
  // Smaller sector oracle uses the full fermionic Fock exponential.
  m.lx=m.ly=2;
  for(unsigned v=0;v<16;++v) {std::vector<int>s(4);for(int i=0;i<4;++i)s[i]=(v&(1u<<i))?1:-1;
   r.record(name,"Fock_sector",v,4,traceProduct({beta*m.generator(s)}),exactTrace({beta*m.generator(s)}));}
 }
 for(auto& name:generatorNames()) {
  auto h=makeGenerator(name);
  for(int t=0;t<samples;++t) {
   std::vector<double> f(h.fieldCount(3));for(auto& v:f)v=1.5*rng.rdNormal();auto p=h.pulses(3,beta,f),reference=p;
   if(name.rfind("creation_only",0)==0) {
    auto base=name=="creation_only_chain" ? openHopping(4):spinCell(2);reference=base.pulses(3,beta,f);
   } else for(auto& a:reference) {
    MatType q=numberMatrix(a);
    for(int i=0;i<h.modes;++i)for(int j=0;j<h.modes;++j)if(i!=j && (name!="directed_scc" || i/2!=j/2))q(i,j)=0.;
    a=numberMajorana(q);
   }
   auto z=traceProduct(p);closeCheck(z,exactTrace(p),name+" Fock");
   r.record(name,"nonhermitian_fields",t,h.modes,z,determinantTrace(reference));
  }
 }
 // Independent quartic Hamiltonian oracle, not a positive Gaussian decoupling claim.
 PairHopping model;int n=2*model.lx*model.ly;MatType h=model.exactHamiltonian(),direct=MatType::Zero(1<<n,1<<n);
 for(int v=0;v<(1<<n);++v)direct(v,v)=model.diagonalEnergy(v);
 for(auto e:model.edges()) {
  MatType ci=fockMatrix(pairMajorana(n,2*e.first,2*e.first+1,1.,false));
  MatType cj=fockMatrix(pairMajorana(n,2*e.second,2*e.second+1,1.,false));
  direct-=model.hopping*(ci*cj.adjoint()+cj*ci.adjoint());
 }
 require((h-direct).norm()<1e-12,"Pair-hopping quartic Hamiltonian signs");
 // Every legal edge move, including spectator single-occupation blockers, closes positively.
 int index=0;
 for(unsigned v=0;v<(1u<<n);++v)for(auto e:model.edges())for(int flip=0;flip<2;++flip) {
  int i=flip?e.second:e.first,j=flip?e.first:e.second;unsigned a=3u<<(2*i),b=3u<<(2*j);
  if((v&a)==a && (v&b)==0) r.record("pair_hopping","all_closed_edge_hops",index++,n,model.closedHopPfaffian(v,v^a^b),1.);
 }
 for(int t=0;t<samples;++t) {
  // Closed paths include a pair moving around the 2D square, idle vertices, and blocked sites.
  std::vector<unsigned> path;
  if(t%3==0)path={3u,12u,192u,48u,3u};
  else if(t%3==1)path={3u|16u,12u|16u,12u|16u,3u|16u};
  else path={3u|192u,12u|192u,12u|48u,3u|48u,3u|192u};
  double dt=.11,ref=1.;
  for(size_t k=0;k+1<path.size();++k){ref*=std::exp(-dt*model.diagonalEnergy(path[k]));if(path[k]!=path[k+1])ref*=dt*model.hopping;}
  r.record("pair_hopping","closed_2d_path",t,n,model.pathPfaffian(path,dt),ref);
 }
 // Sum ALL length-3 worldlines of a two-site model; compare the transfer trace.
 model.lx=1;model.ly=2;n=4;h=model.exactHamiltonian();MatType transfer=eye(16);double dt=.07;
 for(int i=0;i<16;++i)for(int j=0;j<16;++j)if(i!=j)transfer(i,j)=-dt*h(i,j);
 for(int j=0;j<16;++j)transfer.col(j)*=std::exp(-dt*model.diagonalEnergy(j));
 DataType sum=0.;
 for(unsigned a=0;a<16;++a)for(unsigned b=0;b<16;++b)for(unsigned c=0;c<16;++c)
  if(std::abs(transfer(b,a)*transfer(c,b)*transfer(a,c))>0.)sum+=model.pathPfaffian({a,b,c,a},dt);
 r.record("pair_hopping","exhaustive_transfer_trace",0,n,sum,(transfer*transfer*transfer).trace());
 // T_dt=(I+dt*t*sum hops)exp(-dt*D) converges to exp(-beta*H).
 auto error=[&](int steps){double d=beta/steps;MatType t=eye(16),prod=eye(16);for(int i=0;i<16;++i)for(int j=0;j<16;++j)if(i!=j)t(i,j)=-d*h(i,j);for(int j=0;j<16;++j)t.col(j)*=std::exp(-d*model.diagonalEnergy(j));for(int k=0;k<steps;++k)prod=(t*prod).eval();return std::abs(prod.trace()-matrixExp(-beta*h).trace());};
 require(error(128)<.6*error(64),"Worldline discretization does not converge to specified H");
}
}
