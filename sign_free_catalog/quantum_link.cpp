#include "catalog/quantum_link.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
using namespace catalog;
struct Records {
 std::ostream& out;int rows=0,positive=0,zero=0,negative=0;double error=0.,phase=0.;
 Records(std::ostream& o):out(o){out<<std::setprecision(17)<<"case,index,modes,weight_real,weight_imag,reference_real,reference_imag,expected,relative_error,beta,path,times\n";}
 void add(std::string name,int index,int modes,DataType z,DataType ref,int expected=1,const LinkPath& p={}) {
  double e=std::abs(z-ref)/std::max(1.,std::abs(ref));require(std::isfinite(e)&&e<2e-7,"Quantum link Pfaffian/reference mismatch");error=std::max(error,e);
  if(!expected){require(std::abs(z)<1e-9,"Expected zero link path");++zero;}else{require(std::abs(z)>1e-14,"Unresolved link weight");double d=std::abs(z/std::abs(z)-double(expected));require(d<2e-7,"Wrong quantum link signed phase");if(expected==1){++positive;phase=std::max(phase,d);}else ++negative;}
  out<<name<<','<<index<<','<<modes<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<','<<p.beta<<',';for(auto s:p.states)out<<s<<':';out<<',';for(double t:p.times)out<<t<<':';out<<'\n';++rows;
 }
};
static void times(LinkPath& p){p.times.clear();for(int k=1;k<int(p.states.size());++k)p.times.push_back(p.beta*k/p.states.size());}
static DataType physicalWord(const QuantumLink& m,const LinkPath& p,const MatType& h) {
 DataType z=1.;double last=0.;for(int k=0;k<int(p.times.size());++k){z*=std::exp(-(p.times[k]-last)*h(p.states[k],p.states[k]).real())*(-h(p.states[k+1],p.states[k]));last=p.times[k];}return z*std::exp(-(p.beta-last)*h(p.states.back(),p.states.back()).real());
}
static void controls(Records& r) {
 // A closed word that is not a word followed by its reverse.
 QuantumLink square;square.sites=4;LinkPath diamond;std::uint64_t ds=5;diamond.states={ds};
 for(int edge:{0,2,0,2}){ds^=(std::uint64_t(1)<<edge)|(std::uint64_t(1)<<(edge+1))|(std::uint64_t(3)<<(4+2*edge));diamond.states.push_back(ds);}
 times(diamond);require(square.valid(diamond),"Non-retraced quantum-link cycle must close");
 r.add("nonretraced_closed_path",0,square.modes(),square.weight(diamond),square.reference(diamond),1,diamond);
 LinkPath boundary;boundary.states={10,3075,3269,204,10};times(boundary);
 require(square.valid(boundary),"Nonretraced boundary cycle must close");
 r.add("nonretraced_boundary_path",0,square.modes(),square.weight(boundary),square.reference(boundary),1,boundary);
 QuantumLink m;int n=m.modes();std::uint64_t state=3;LinkPath broken;broken.states={state};std::vector<Gaussian> factors;
 // A cyclic exchange of two c particles winds once and leaves all link pairs raised.
 for(int edge:{1,0,2}){int i=edge,j=(edge+1)%3;auto after=state^(std::uint64_t(1)<<i)^(std::uint64_t(1)<<j)^(std::uint64_t(3)<<(3+2*edge));factors.push_back(linkOccupationProjector(n,state));factors.push_back(gaussian(m.hopRotation(state,after)));state=after;broken.states.push_back(state);}
 times(broken);require(!m.valid(broken),"Winding incorrectly admitted with nonperiodic link flux");r.add("winding_nonclosed_flux",0,n,m.weight(broken),0.,0,broken);
 // Allowing independent link flips completes the same negative fermion exchange.
 for(int edge=0;edge<3;++edge){factors.push_back(linkOccupationProjector(n,state));auto c=pairMajorana(n,3+2*edge,4+2*edge,1.,false);factors.push_back(gaussian(-.5*std::acos(-1.)*(c-c.adjoint()).eval()));state^=std::uint64_t(3)<<(3+2*edge);}
 require(state==3,"Negative link-reset word not closed");factors.push_back(linkOccupationProjector(n,state));
 r.add("independent_link_flip_breaks_sign",0,n,traceGaussianFactors(factors),-1.,-1);
 // Odd auxiliary occupations block their bond, while remaining in full Fock space.
 LinkPath blocked;blocked.states={std::uint64_t(1)|(std::uint64_t(1)<<3),std::uint64_t(2)|(std::uint64_t(2)<<3),std::uint64_t(1)|(std::uint64_t(1)<<3)};times(blocked);r.add("odd_link_blocker",0,n,m.weight(blocked),0.,0,blocked);
 LinkPath invalid;invalid.states={0};invalid.beta=-1;r.add("invalid_time",0,n,m.weight(invalid),0.,0,invalid);
 // Zero edge current does not prevent exchange on a branching graph.
 // Star edges 0->1,0->2,0->3: two particles use the empty third arm to exchange.
 std::uint64_t ts=246;std::vector<Gaussian> tree;
 for(auto move:std::vector<std::pair<int,int>>{{1,0},{0,3},{2,0},{0,1},{3,0},{0,2}}) {
  int source=move.first,dest=move.second,edge=std::max(source,dest)-1;bool forward=source==0;
  MatType k=MatType::Zero(10,10);k(dest,source)=1.;k(source,dest)=-1.;auto pair=pairMajorana(10,4+2*edge,5+2*edge,1.,false);
  tree.push_back(linkOccupationProjector(10,ts));tree.push_back(gaussian(.5*std::acos(-1.)*(numberMajorana(k)+(forward?1.:-1.)*(pair-pair.adjoint()))));
  ts^=(std::uint64_t(1)<<source)|(std::uint64_t(1)<<dest)|(std::uint64_t(3)<<(4+2*edge));
 }
 require(ts==246,"Branching link exchange did not close");tree.push_back(linkOccupationProjector(10,ts));
 r.add("branching_graph_breaks_path_sign",0,10,traceGaussianFactors(tree),-1.,-1);
}
int main(int argc,char** argv) {
 try {
  mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=64;unsigned seed=230508917;std::string csv;
  for(int i=1;i<argc;++i){std::string key=argv[i];require(i+1<argc,"Missing quantum-link CLI value");std::string value=argv[++i];if(key=="--samples")samples=std::stoi(value);else if(key=="--seed")seed=std::stoul(value);else if(key=="--csv")csv=value;else throw std::runtime_error("Unknown option "+key);}
  require(samples>0,"Need positive samples");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open link CSV");}Records r(csv.empty()?std::cout:file);controls(r);std::mt19937_64 rng(seed);
  for(int geometry=0;geometry<3;++geometry){QuantumLink m;m.sites=geometry==0?2:geometry==1?3:4;m.periodic=geometry!=0;int n=m.modes();MatType h;
   if(n<=9){h=m.exactHamiltonian();require((h-h.adjoint()).norm()<1e-11,"Quantum link physical H not Hermitian");for(unsigned s=0;s<unsigned(h.rows());++s){require(std::abs(h(s,s)-m.diagonalEnergy(s))<1e-12,"Quantum link diagonal Hamiltonian mismatch");auto next=m.nextStates(s);for(unsigned t=0;t<unsigned(h.cols());++t)if(t!=s)require((std::abs(h(t,s))>1e-12)==(std::find(next.begin(),next.end(),t)!=next.end()),"Physical link transition mismatch");}}
   // Independent exhaustive phase-gauge certificate on every physical Fock edge.
   if(h.size())for(unsigned s=0;s<unsigned(h.rows());++s)for(auto t:m.nextStates(s)) {
    auto gauge=[&](std::uint64_t v){if(!m.periodic)return 1.;int particles=__builtin_popcountll(v&((std::uint64_t(1)<<m.sites)-1));int link=(v>>(m.sites+2*(int(m.edges().size())-1)))&1;return ((particles-1)*link)%2?-1.:1.;};
    require(std::abs(gauge(s)*gauge(t)*h(t,s)+m.hopping)<1e-11,"Quantum link full graph phase-gauge certificate");
   }
   int moved=0;
   for(int sample=0;sample<samples;++sample){LinkPath p;p.states={rng()%(std::uint64_t(1)<<n)};
    for(int k=0;k<2+sample%4;++k){auto next=m.nextStates(p.states.back());if(next.empty())break;p.states.push_back(next[rng()%next.size()]);}
    int length=p.states.size();for(int k=length-2;k>=0;--k)p.states.push_back(p.states[k]);times(p);if(length>1)++moved;
    auto z=m.weight(p),ref=m.reference(p);if(h.size())require(std::abs(ref-physicalWord(m,p,h))<1e-10,"Independent full Fock CT path oracle mismatch");r.add(m.periodic?"quantum_link_cycle":"quantum_link_open",sample,n,z,ref,1,p);
   }require(moved>0,"Quantum link random ensemble never hopped");
   // Exhaustively sum diagonal and one forward/backward pair at a fixed time:
   // Tr[D2 T D1 T D0] from full Fock must equal the sum over allowed paths.
   if(geometry==0){MatType kinetic=-h;for(int i=0;i<h.rows();++i)kinetic(i,i)=0.;MatType d0=MatType::Zero(h.rows(),h.rows()),d1=d0,d2=d0;for(int s=0;s<h.rows();++s){double e=h(s,s).real();d0(s,s)=std::exp(-.2*e);d1(s,s)=std::exp(-.3*e);d2(s,s)=std::exp(-.3*e);}DataType sum=0.;int index=0;for(int s=0;s<h.rows();++s)for(auto t:m.nextStates(s)){LinkPath p;p.states={std::uint64_t(s),t,std::uint64_t(s)};p.times={.2,.5};auto z=m.weight(p);r.add("complete_two_hop_word",index++,n,z,physicalWord(m,p,h),1,p);sum+=z;}require(std::abs(sum-(d2*kinetic*d1*kinetic*d0).trace())<1e-10,"Complete CT path sum mismatch");}
  }
  std::cerr<<"PASS quantum_link rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero<<" negative="<<r.negative<<" max_relative_error="<<r.error<<" max_phase_error="<<r.phase<<" seed="<<seed<<'\n';
 }catch(const std::exception& e){std::cerr<<"FAIL quantum_link: "<<e.what()<<'\n';return 1;}return 0;
}
