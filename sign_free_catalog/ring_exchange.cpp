#include "catalog/ring_exchange.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
using namespace catalog;
struct Records {
 std::ostream& out;int rows=0,positive=0,zero=0,negative=0;double error=0.,phase=0.;
 Records(std::ostream& o):out(o){out<<std::setprecision(17)<<"case,index,modes,weight_real,weight_imag,reference_real,reference_imag,expected,relative_error,beta,path,times\n";}
 void add(std::string name,int index,int modes,DataType z,DataType ref,int expected=1,const LinkPath& p={}){double e=std::abs(z-ref)/std::max(1.,std::abs(ref));require(std::isfinite(e)&&e<2e-7,"Ring signed Pfaffian/Fock mismatch");error=std::max(error,e);if(expected==0){require(std::abs(z)<1e-10,"Ring expected zero");++zero;}else{require(std::abs(z)>1e-14,"Unresolved ring weight");double d=std::abs(z/std::abs(z)-double(expected));require(d<2e-7,"Ring wrong sign");if(expected==1){++positive;phase=std::max(phase,d);}else ++negative;}out<<name<<','<<index<<','<<modes<<','<<z.real()<<','<<z.imag()<<','<<ref.real()<<','<<ref.imag()<<','<<expected<<','<<e<<','<<p.beta<<',';for(auto s:p.states)out<<s<<':';out<<',';for(double t:p.times)out<<t<<':';out<<'\n';++rows;}
};
static void times(LinkPath& p){p.times.clear();for(int k=1;k<int(p.states.size());++k)p.times.push_back(p.beta*k/p.states.size());}
static DataType physicalWord(const MatType& h,const LinkPath& p){DataType z=1.;double last=0.;for(int k=0;k<int(p.times.size());++k){z*=std::exp(-(p.times[k]-last)*h(p.states[k],p.states[k]).real())*(-h(p.states[k+1],p.states[k]));last=p.times[k];}return z*std::exp(-(p.beta-last)*h(p.states.back(),p.states.back()).real());}
int main(int argc,char** argv){try{
 mkl_set_num_threads(1);mkl_set_dynamic(0);int samples=64;unsigned seed=409112;std::string csv;for(int i=1;i<argc;++i){std::string k=argv[i];require(i+1<argc,"Missing ring CLI value");std::string v=argv[++i];if(k=="--samples")samples=std::stoi(v);else if(k=="--seed")seed=std::stoul(v);else if(k=="--csv")csv=v;else throw std::runtime_error("Unknown option "+k);}require(samples>0,"Need positive ring samples");std::ofstream file;if(!csv.empty()){file.open(csv);require(file.good(),"Cannot open ring CSV");}Records r(csv.empty()?std::cout:file);std::mt19937_64 rng(seed);
 RingExchange m;auto h=m.exactHamiltonian();require((h-h.adjoint()).norm()<1e-12,"Ring physical Hamiltonian not Hermitian");for(int s=0;s<h.rows();++s){require(std::abs(h(s,s)-m.diagonalEnergy(s))<1e-12,"Ring diagonal mismatch");auto next=m.nextStates(s);for(int t=0;t<h.rows();++t)if(t!=s){bool present=std::find(next.begin(),next.end(),t)!=next.end();require(present==(std::abs(h(t,s))>1e-12),"Ring CAR Hamiltonian transition mismatch");if(present)require(std::abs(h(t,s)+m.moveCoupling(s,t))<1e-12,"Pair ring CAR coefficient mismatch");}}
 LinkPath cycle;cycle.states={195,60,51,195};times(cycle);r.add("ring_and_two_bond_hops",0,8,m.weight(cycle),physicalWord(h,cycle),1,cycle);
 auto bad=m;bad.ring=-m.ring;require(!bad.signFreeParameters(),"Negative ring domain not identified");r.add("negative_ring_frustrated_cycle",0,8,bad.weight(cycle),physicalWord(bad.exactHamiltonian(),cycle),-1,cycle);
 auto zero=m;zero.ring=0.;r.add("zero_ring_coupling",0,8,zero.weight(cycle),0.,0,cycle);
 LinkPath blocked;blocked.states={1,2,1};times(blocked);r.add("odd_blocked_ring_path",0,8,m.weight(blocked),0.,0,blocked);
 // A single ring term has a matrix element between states differing on four cells.
 // Every sum of one/two-cell terms has identically zero at this position.
 auto noRing=m;noRing.ring=0.;require(std::abs(h(60,195)+m.ring)<1e-12 && std::abs(noRing.exactHamiltonian()(60,195))<1e-12,"Missing irreducible four-cell interaction");
 MatType T=-h;for(int s=0;s<T.rows();++s)T(s,s)=0.;MatType D=MatType::Zero(h.rows(),h.rows());for(int s=0;s<h.rows();++s)D(s,s)=std::exp(-.2*h(s,s).real());DataType sum=0.;int index=0;
 for(int s=0;s<h.rows();++s)for(auto a:m.nextStates(s))for(auto b:m.nextStates(a)){auto next=m.nextStates(b);if(std::find(next.begin(),next.end(),s)==next.end())continue;LinkPath p;p.states={std::uint64_t(s),a,b,std::uint64_t(s)};p.times={.2,.4,.6};auto z=m.weight(p);r.add("complete_three_vertex_paths",index++,8,z,physicalWord(h,p),1,p);sum+=z;}
 require(index>0 && std::abs(sum-(D*T*D*T*D*T*D).trace())<1e-10,"Complete ring CT word sum mismatch");
 for(int ly:{2,3}){RingExchange model;model.ly=ly;int moved=0;for(int sample=0;sample<samples;++sample){LinkPath p;p.states={rng()%(std::uint64_t(1)<<model.modes())};for(int k=0;k<2+sample%4;++k){auto next=model.nextStates(p.states.back());if(next.empty())break;p.states.push_back(next[rng()%next.size()]);}int length=p.states.size();for(int k=length-2;k>=0;--k)p.states.push_back(p.states[k]);times(p);if(length>1)++moved;auto ref=model.reference(p);if(ly==2)require(std::abs(ref-physicalWord(h,p))<1e-12,"Random ring full Fock oracle mismatch");r.add("random_pair_ring_path",sample,model.modes(),model.weight(p),ref,1,p);}require(moved>0,"Ring random ensemble never moved");}
 std::cerr<<"PASS ring_exchange rows="<<r.rows<<" positive="<<r.positive<<" zero="<<r.zero<<" negative="<<r.negative<<" max_relative_error="<<r.error<<" max_phase_error="<<r.phase<<" seed="<<seed<<'\n';
 }catch(const std::exception& e){std::cerr<<"FAIL ring_exchange: "<<e.what()<<'\n';return 1;}return 0;}
