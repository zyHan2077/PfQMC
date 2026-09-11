#include "checks.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
using namespace catalog;
int main(int argc,char** argv) {
 try {
  mkl_set_num_threads(1); mkl_set_dynamic(0);
  std::string model="all",csv; int samples=16,sweeps=6,seed=91011;
  for(int i=1;i<argc;++i) {
   std::string a=argv[i];
   if(a=="--help") { std::cout<<"sign_free_models --model all|checks|MODEL --samples N --sweeps N --seed N --csv PATH\nUse --list for Gaussian HS Hamiltonians. checks includes static FK/BdG, non-Hermitian generators, worldlines and controls.\n"; return 0; }
   if(a=="--list") { for(auto& n:modelNames()) {auto h=makeModel(n); std::cout<<n<<" modes="<<h.modes<<" fields/slice="<<h.channels.size()<<"\n";} return 0; }
   require(i+1<argc,"Missing value for "+a); std::string v=argv[++i];
   if(a=="--model") model=v; else if(a=="--csv") csv=v;
   else if(a=="--samples") samples=std::stoi(v); else if(a=="--sweeps") sweeps=std::stoi(v);
   else if(a=="--seed") seed=std::stoi(v); else throw std::runtime_error("Unknown option "+a);
  }
  require(samples>0 && sweeps>=2,"Use samples>0 and sweeps>=2 (both sweep directions)");
  std::ofstream file; if(!csv.empty()) {file.open(csv);require(file.good(),"Cannot write CSV: "+csv);}
  Recorder recorder(csv.empty() ? std::cout:file); rdGenerator rng(seed);
  if(model=="all" || model=="checks") {verifyCore(recorder,rng);verifySpecials(recorder,rng,samples);}
  auto names=model=="all" ? modelNames():model=="checks" ? std::vector<std::string>{}:std::vector<std::string>{model};
  for(auto& name:names) {
   auto start=std::chrono::steady_clock::now(); verifyModel(recorder,rng,name,samples,sweeps);
   std::cerr<<"PASS "<<name<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
  }
  std::cerr<<std::setprecision(12)<<"PASS rows="<<recorder.rows<<" positive="<<recorder.positive<<" zero="<<recorder.zeros
    <<" negative_controls="<<recorder.negative<<" max_relative_error="<<recorder.maxError<<" max_positive_phase_error="<<recorder.maxPhaseError
    <<" seed="<<seed<<" samples="<<samples<<" sweeps="<<sweeps<<"\n";
 } catch(const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<"\n"; return 1; }
 return 0;
}
