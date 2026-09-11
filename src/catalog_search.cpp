#include "catalog/search_models.h"
#include <map>
#include <numeric>
namespace catalog {
static MatType eye(int n) {return MatType::Identity(n,n);}
static MatType kron(const MatType& a,const MatType& b) {
 MatType c(a.rows()*b.rows(),a.cols()*b.cols());
 for(int i=0;i<a.rows();++i)for(int j=0;j<a.cols();++j)c.block(i*b.rows(),j*b.cols(),b.rows(),b.cols())=a(i,j)*b;
 return c;
}
MatType conjugateTensorAction(const MatType& k) {return kron(k,eye(k.rows()))+kron(eye(k.rows()),k.conjugate());}
std::vector<std::vector<int>> symmetricStates(int d,int q) {
 require(d>=2 && d<=8 && q>=1 && q<=8,"Symmetric power reference limits: 2<=d<=8, 1<=q<=8");
 std::vector<std::vector<int>> out;std::vector<int> v(d);
 std::function<void(int,int)> visit=[&](int pos,int left) {
  if(pos==d-1){v[pos]=left;out.push_back(v);return;}
  for(int k=left;k>=0;--k){v[pos]=k;visit(pos+1,left-k);}
 };visit(0,q);require(out.size()<=512,"Symmetric power dense orbital limit: 512");return out;
}
MatType symmetricPowerAction(const MatType& x,int q) {
 int d=x.rows();auto states=symmetricStates(d,q);std::map<std::vector<int>,int> index;
 for(int i=0;i<int(states.size());++i)index[states[i]]=i;
 MatType h=MatType::Zero(states.size(),states.size());
 for(int col=0;col<int(states.size());++col)for(int j=0;j<d;++j) {
  h(col,col)+=double(states[col][j])*x(j,j);
  if(states[col][j])for(int i=0;i<d;++i)if(i!=j) {
   auto next=states[col];--next[j];++next[i];
   h(index.at(next),col)+=std::sqrt(double((states[col][i]+1)*states[col][j]))*x(i,j);
  }
 }
 return h;
}
static std::vector<MatType> tracelessHermitian(int d) {
 std::vector<MatType> a;
 for(int k=1;k<d;++k){MatType x=MatType::Zero(d,d);double c=1./std::sqrt(double(k*(k+1)));for(int i=0;i<k;++i)x(i,i)=c;x(k,k)=-k*c;a.push_back(x);}
 for(int i=0;i<d;++i)for(int j=i+1;j<d;++j){MatType x=MatType::Zero(d,d);x(i,j)=x(j,i)=1./std::sqrt(2.);a.push_back(x);x(i,j)=1.i/std::sqrt(2.);x(j,i)=-1.i/std::sqrt(2.);a.push_back(x);}
 return a;
}
static MatType physical(const SearchModel& m,const MatType& a,double mu) {
 MatType h;
 if(m.family==SearchModel::Lorentz)h=a;
 else if(m.family==SearchModel::LorentzSymmetric)h=lorentzSymmetricAction(a);
 else if(m.family==SearchModel::ConjugateTensor)h=conjugateTensorAction(a);
 else if(m.family==SearchModel::QuaternionHermitian || m.family==SearchModel::QuaternionAntiHermitian)h=quaternionAction(a,m.family==SearchModel::QuaternionAntiHermitian);
 else h=symmetricPowerAction(a,m.degree);
 h+=mu*eye(m.cellModes);return h;
}
static SearchModel assemble(SearchModel m,const MatType& a0,double mu0,const std::vector<MatType>& directions,
                            const std::vector<double>& chemical,int lx,int ly) {
 require(lx>0 && ly>0 && directions.size()==chemical.size(),"Invalid candidate geometry");
 m.cells=lx*ly;m.h.modes=m.cells*m.cellModes;m.h.q0=MatType::Zero(2*m.h.modes,2*m.h.modes);
 m.auxiliary0.resize(m.cells,a0);m.chemical0.resize(m.cells,mu0);
 for(int r=0;r<m.cells;++r)m.h.q0.block(2*r*m.cellModes,2*r*m.cellModes,2*m.cellModes,2*m.cellModes)=numberMajorana(physical(m,a0,mu0));
 auto add=[&](int r,int a,int s,int b,double eta,const std::string& name,double coupling) {
  std::vector<MatType> aux(m.cells,MatType::Zero(m.auxiliaryDim,m.auxiliaryDim));std::vector<double> mus(m.cells,0.);
  aux[r]+=directions[a];mus[r]+=chemical[a];if(s>=0){aux[s]+=eta*directions[b];mus[s]+=eta*chemical[b];}
  MatType full=MatType::Zero(2*m.h.modes,2*m.h.modes);
  for(int cell=0;cell<m.cells;++cell)full.block(2*cell*m.cellModes,2*cell*m.cellModes,2*m.cellModes,2*m.cellModes)=numberMajorana(physical(m,aux[cell],mus[cell]));
  m.h.channels.push_back({name,full,coupling});m.auxiliaryChannels.push_back(aux);m.chemicalChannels.push_back(mus);
 };
 for(int r=0;r<m.cells;++r)for(int a=0;a<int(directions.size());++a)add(r,a,-1,0,0.,"onsite_"+std::to_string(r)+"_"+std::to_string(a),.12);
 int edge=0;
 for(int x=0;x<lx;++x)for(int y=0;y<ly;++y) {
  int r=x*ly+y;std::vector<int> neighbors;if(x+1<lx)neighbors.push_back(r+ly);if(y+1<ly)neighbors.push_back(r+1);
  for(int s:neighbors){add(r,0,s,1,edge%2 ? -1.:1.,"mixed_edge_"+std::to_string(edge),.16);++edge;}
 }
 m.h.validate();return m;
}
SearchModel lorentzModel(int boosts,int lx,int ly) {
 require(boosts>=3 && boosts<=15,"New Lorentz series requires 3<=boosts<=15; boosts=2 repeats spin1");
 SearchModel m;m.family=SearchModel::Lorentz;m.auxiliaryDim=m.cellModes=boosts+1;
 MatType a0=MatType::Zero(m.cellModes,m.cellModes);std::vector<MatType> dirs={a0};std::vector<double> mus={1.};
 for(int i=1;i<=boosts;++i){MatType b=MatType::Zero(m.cellModes,m.cellModes);b(0,i)=b(i,0)=1.;dirs.push_back(b);mus.push_back(0.);a0+=.16*std::cos(double(i))*b;}
 m.h.description="SO+(1,m) boosts plus independent common density; finite-range shared negative squares";
 return assemble(m,a0,.13,dirs,mus,lx,ly);
}
SearchModel conjugateTensorModel(int d,int lx,int ly) {
 require(d>=2 && d<=6,"Conjugate tensor reference limits: 2<=d<=6");SearchModel m;
 m.family=SearchModel::ConjugateTensor;m.auxiliaryDim=d;m.cellModes=d*d;
 auto basis=tracelessHermitian(d);std::vector<MatType> dirs={MatType::Zero(d,d)};std::vector<double> mus={1.};
 MatType a0=MatType::Zero(d,d);
 for(int i=0;i<int(basis.size());++i){dirs.push_back(basis[i]);mus.push_back(0.);a0+=.08*std::sin(i+1.)*basis[i];}
 m.h.description="Conjugate tensor (not direct sum), including imaginary Hermitian auxiliary directions";
 return assemble(m,a0,.11,dirs,mus,lx,ly);
}
SearchModel symmetricPowerModel(int d,int q,int lx,int ly) {
 require(d>=3 && q>=2,"d=2 repeats integer-spin representations; q=1 repeats an open chain");SearchModel m;
 m.family=SearchModel::SymmetricPower;m.auxiliaryDim=d;m.degree=q;m.cellModes=symmetricStates(d,q).size();
 MatType a0=MatType::Zero(d,d);std::vector<MatType> dirs;std::vector<double> mus;
 for(int i=0;i<d;++i){MatType a=MatType::Zero(d,d);a(i,i)=1.;dirs.push_back(a);mus.push_back(0.);a0(i,i)=.09*(i-.5*(d-1));}
 for(int i=0;i<d-1;++i)a0(i,i+1)=a0(i+1,i)=.3;
 m.h.description="Symmetric-power orbitals of a totally nonnegative Jacobi semigroup; only correlated diagonal HS";
 return assemble(m,a0,0.,dirs,mus,lx,ly);
}
std::vector<MatType> quaternionBasis(int n,bool anti) {
 require(n>=2 && n<=6,"Quaternion reference limits: 2<=n<=6");
 std::vector<MatType> units(4,MatType::Zero(2,2));units[0]=eye(2);
 units[1](0,0)=1.i;units[1](1,1)=-1.i;units[2](0,1)=1.;units[2](1,0)=-1.;units[3](0,1)=units[3](1,0)=1.i;
 std::vector<MatType> basis;
 for(int i=0;i<n;++i)for(int a=anti?1:0;a<(anti?4:1);++a){MatType t=MatType::Zero(2*n,2*n);t.block(2*i,2*i,2,2)=units[a];basis.push_back(t);}
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)for(int a=0;a<4;++a){MatType t=MatType::Zero(2*n,2*n);t.block(2*i,2*j,2,2)=units[a]/std::sqrt(2.);t.block(2*j,2*i,2,2)=(anti?-1.:1.)*units[a].adjoint()/std::sqrt(2.);basis.push_back(t);}
 return basis;
}
MatType quaternionAction(const MatType& k,bool anti) {
 require(k.rows()==k.cols() && k.rows()%2==0,"Quaternion matrix dimensions");int n=k.rows()/2;
 MatType j=MatType::Zero(2*n,2*n);for(int i=0;i<n;++i){j(2*i,2*i+1)=1.;j(2*i+1,2*i)=-1.;}
 require(k.allFinite() && (k-k.adjoint()).norm()<1e-10 && (k*j-j*k.conjugate()).norm()<1e-10,"Expected quaternion-Hermitian generator");
 auto basis=quaternionBasis(n,anti);MatType h(basis.size(),basis.size());
 for(int b=0;b<int(basis.size());++b){MatType mapped=k*basis[b]+basis[b]*k;for(int a=0;a<int(basis.size());++a)h(a,b)=.5*(basis[a].adjoint()*mapped).trace().real();}
 require((h-h.transpose()).norm()<1e-10,"Quaternion orbital action is not symmetric");return h;
}
SearchModel quaternionModel(int n,bool anti,int lx,int ly) {
 require(n>=2 && n<=6,"Quaternion reference limits: 2<=n<=6");SearchModel m;
 m.family=anti?SearchModel::QuaternionAntiHermitian:SearchModel::QuaternionHermitian;m.auxiliaryDim=2*n;m.cellModes=n*(2*n+(anti?1:-1));
 auto basis=quaternionBasis(n,false);std::vector<MatType> dirs={MatType::Zero(2*n,2*n)};std::vector<double> mus={1.};
 for(int k=1;k<n;++k){MatType a=MatType::Zero(2*n,2*n);double c=1./std::sqrt(double(k*(k+1)));for(int i=0;i<k;++i)a+=c*basis[i];a-=k*c*basis[k];dirs.push_back(a);mus.push_back(0.);}
 for(int k=n;k<int(basis.size());++k){dirs.push_back(basis[k]);mus.push_back(0.);}
 MatType a0=MatType::Zero(2*n,2*n);for(int k=1;k<int(dirs.size());++k)a0+=.06*std::sin(double(k))*dirs[k];
 m.h.description=anti?"Quaternion anti-Hermitian matrix orbitals (complexified Sym2)":"Quaternion Hermitian matrix orbitals (complexified exterior2)";
 return assemble(m,a0,.12,dirs,mus,lx,ly);
}
std::vector<MatType> lorentzSymmetricBasis(int n) {
 require(n>=5 && n<=9,"Traceless symmetric Lorentz family requires 5<=n<=9");
 MatType j=-eye(n);j(0,0)=1.;MatType v=MatType::Zero(n,1);v(0,0)=1.;v-=j.diagonal()/std::sqrt(double(n));
 MatType u=eye(n)-2.*v*v.transpose()/(v.transpose()*v)(0,0);
 std::vector<MatType> basis;
 for(int a=1;a<n;++a){MatType t=MatType::Zero(n,n);t.diagonal()=u.col(a);basis.push_back(t);}
 for(int a=0;a<n;++a)for(int b=a+1;b<n;++b){MatType t=MatType::Zero(n,n);t(a,b)=t(b,a)=1./std::sqrt(2.);basis.push_back(t);}
 return basis;
}
MatType lorentzSymmetricAction(const MatType& k) {
 require(k.rows()==k.cols() && (k-k.transpose()).norm()<1e-10 && k.imag().norm()<1e-10,"Real symmetric Lorentz generator required");
 auto basis=lorentzSymmetricBasis(k.rows());MatType h(basis.size(),basis.size());
 // Projection is explicit; callers must separately enforce kJ+Jk=0 for positivity.
 for(int b=0;b<int(basis.size());++b)for(int a=0;a<int(basis.size());++a)h(a,b)=(basis[a]*(k*basis[b]+basis[b]*k)).trace();
 return h;
}
SearchModel lorentzSymmetricModel(int n,int lx,int ly) {
 SearchModel m;m.family=SearchModel::LorentzSymmetric;m.auxiliaryDim=n;m.cellModes=lorentzSymmetricBasis(n).size();
 MatType a0=MatType::Zero(n,n);std::vector<MatType> dirs={a0};std::vector<double> mus={1.};
 for(int i=1;i<n;++i){MatType b=MatType::Zero(n,n);b(0,i)=b(i,0)=1.;dirs.push_back(b);mus.push_back(0.);a0+=.12*std::sin(double(i))*b;}
 m.h.description="Lorentz traceless symmetric matrix orbitals, one polynomial representation construction";
 return assemble(m,a0,.09,dirs,mus,lx,ly);
}
DataType SearchModel::reference(int slices,double beta,const std::vector<double>& fields) const {
 require(slices>0 && int(fields.size())==h.fieldCount(slices),"Candidate reference field dimensions");
 std::vector<MatType> g(cells,eye(auxiliaryDim));std::vector<double> sums(cells,0.);
 double dt=beta/slices;int pos=0;
 auto multiply=[&](const std::vector<MatType>& a,const std::vector<double>& mu,double scale){for(int r=0;r<cells;++r){g[r]=(matrixExp(scale*a[r])*g[r]).eval();sums[r]+=scale*mu[r];}};
 for(int l=0;l<slices;++l){multiply(auxiliary0,chemical0,.5*dt);for(int c=0;c<int(h.channels.size());++c)multiply(auxiliaryChannels[c],chemicalChannels[c],std::sqrt(dt*h.channels[c].coupling)*fields[pos++]);multiply(auxiliary0,chemical0,.5*dt);}
 DataType z=1.;
 for(int r=0;r<cells;++r) {
  Eigen::ComplexEigenSolver<MatType> es(g[r]);require(es.info()==Eigen::Success,"Candidate reference spectrum");auto ev=es.eigenvalues();double R=std::exp(sums[r]);DataType w=std::exp(-.5*cellModes*sums[r]);
  if(family==Lorentz){for(int i=0;i<auxiliaryDim;++i)w*=1.+R*ev(i);}
  else if(family==LorentzSymmetric){for(int i=0;i<auxiliaryDim;++i)for(int j=i;j<auxiliaryDim;++j)w*=1.+R*ev(i)*ev(j);w/=1.+R;}
  else if(family==ConjugateTensor){for(int i=0;i<auxiliaryDim;++i)w*=1.+R*std::norm(ev(i));for(int i=0;i<auxiliaryDim;++i)for(int j=i+1;j<auxiliaryDim;++j)w*=std::norm(1.+R*ev(i)*std::conj(ev(j)));}
  else if(family==SymmetricPower) {
   w=1.;for(int i=0;i<auxiliaryDim;++i)require(ev(i).real()>0 && std::abs(ev(i).imag())<1e-8*ev(i).real(),"TN auxiliary eigenvalue not positive real");
   for(auto& state:symmetricStates(auxiliaryDim,degree)){double energy=sums[r];for(int i=0;i<auxiliaryDim;++i)energy+=state[i]*std::log(ev(i).real());w*=2.*std::cosh(.5*energy);}
  } else {
   // Quaternion complex eigenvalues occur in conjugate pairs, including even real multiplicities.
   std::vector<DataType> remaining(ev.data(),ev.data()+ev.size()),half;
   while(!remaining.empty()) {
    auto x=remaining.back();remaining.pop_back();require(!remaining.empty(),"Unpaired quaternion eigenvalue");
    int match=0;for(int i=1;i<int(remaining.size());++i)if(std::abs(remaining[i]-std::conj(x))<std::abs(remaining[match]-std::conj(x)))match=i;
    require(std::abs(remaining[match]-std::conj(x))<1e-7*std::max(1.,std::abs(x)),"Quaternion spectral pairing failed");
    half.push_back(.5*(x+std::conj(remaining[match])));remaining.erase(remaining.begin()+match);
   }
   for(auto x:half)w*=1.+R*std::norm(x);
   for(int i=0;i<int(half.size());++i)for(int j=i+1;j<int(half.size());++j)w*=std::norm(1.+R*half[i]*half[j])*std::norm(1.+R*half[i]*std::conj(half[j]));
   if(family==QuaternionAntiHermitian)for(auto x:half)w*=std::norm(1.+R*x*x);
  }
  z*=w;
 }
 return z;
}
std::vector<std::string> searchModelNames(){return {"lorentz4","lorentz5","lorentz8","lorentz4_lattice","conjugate_tensor2","conjugate_tensor3","conjugate_tensor3_lattice","symmetric_power3_2","symmetric_power3_3","symmetric_power4_2","symmetric_power3_2_lattice","quaternion_herm2","quaternion_herm3","quaternion_herm3_lattice","quaternion_anti2","quaternion_anti3","quaternion_anti2_lattice","lorentz_sym5","lorentz_sym6","lorentz_sym5_lattice"};}
SearchModel makeSearchModel(const std::string& id) {
 SearchModel m;
 if(id=="lorentz4")m=lorentzModel();else if(id=="lorentz5")m=lorentzModel(4);else if(id=="lorentz8")m=lorentzModel(7);
 else if(id=="lorentz4_lattice")m=lorentzModel(3,2,2);
 else if(id=="conjugate_tensor2")m=conjugateTensorModel(2);else if(id=="conjugate_tensor3")m=conjugateTensorModel(3);
 else if(id=="conjugate_tensor3_lattice")m=conjugateTensorModel(3,2,1);
 else if(id=="symmetric_power3_2")m=symmetricPowerModel();else if(id=="symmetric_power3_3")m=symmetricPowerModel(3,3);
 else if(id=="symmetric_power4_2")m=symmetricPowerModel(4,2);else if(id=="symmetric_power3_2_lattice")m=symmetricPowerModel(3,2,2,2);
 else if(id=="quaternion_herm2")m=quaternionModel(2);else if(id=="quaternion_herm3")m=quaternionModel(3);
 else if(id=="quaternion_herm3_lattice")m=quaternionModel(3,false,2,1);
 else if(id=="quaternion_anti2")m=quaternionModel(2,true);else if(id=="quaternion_anti3")m=quaternionModel(3,true);
 else if(id=="quaternion_anti2_lattice")m=quaternionModel(2,true,2,1);
 else if(id=="lorentz_sym5")m=lorentzSymmetricModel(5);
 else if(id=="lorentz_sym6")m=lorentzSymmetricModel(6);
 else if(id=="lorentz_sym5_lattice")m=lorentzSymmetricModel(5,2,1);
 else throw std::runtime_error("Unknown search model: "+id);
 m.h.id=id;return m;
}
}
