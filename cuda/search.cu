#include "common.h"
#ifdef KEYHUNT_CPU_TEST
// Test-only runtime shim: production builds always require the CUDA runtime.
#include <cstring>
#include <cstdlib>
using cudaError_t = int;
static const int cudaSuccess=0, cudaMemcpyDeviceToHost=0;
static const char *cudaGetErrorString(int) { return "test allocation failure"; }
static int cudaMalloc(void **p, size_t n) { *p=std::malloc(n); return *p?0:1; }
static int cudaFree(void *p) { std::free(p); return 0; }
static int cudaMemcpy(void *a,const void *b,size_t n,int) { std::memcpy(a,b,n); return 0; }
static int cudaGetLastError() { return 0; }
static int cudaSetDevice(int) { return 0; }
static int cudaDeviceSynchronize() { return 0; }
#else
#include <cuda_runtime.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

static void check(cudaError_t e, const char *where) { if(e!=cudaSuccess) throw std::runtime_error(std::string(where)+": "+cudaGetErrorString(e)); }
#define CUDA(x) check((x),#x)

// Each thread computes one point independently. This is deliberately simple,
// not an optimized GPU implementation; CPU code performs the exact matching.
HD inline Affine calculate_point(U256 base, uint64_t step, unsigned i, Affine target, bool giant) {
  U256 k=add_raw(base,from64((uint64_t)i*step));
  Point p=multiply(k);
  if(giant) {
    // target - kG: negate the Jacobian y coordinate before mixed addition.
    if(!iszero(p.z) && !iszero(p.y)) p.y=sub_raw(prime(),p.y);
    p=add_mixed(p,target);
  }
  return affine(p);
}
#ifndef KEYHUNT_CPU_TEST
__global__ void point_kernel(U256 base, uint64_t step, unsigned count, Affine *output, Affine target, int giant) {
  unsigned i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<count) output[i]=calculate_point(base,step,i,target,giant!=0);
}
#endif

struct DeviceBuffer {
  Affine *p=nullptr;
  explicit DeviceBuffer(size_t n) { CUDA(cudaMalloc((void **)&p,n*sizeof(Affine))); }
  ~DeviceBuffer() { if(p) cudaFree(p); }
  DeviceBuffer(const DeviceBuffer&)=delete;
  DeviceBuffer& operator=(const DeviceBuffer&)=delete;
};
static void points(DeviceBuffer &device, std::vector<Affine> &out, U256 base, uint64_t step, unsigned count, Affine target, bool giant) {
#ifdef KEYHUNT_CPU_TEST
  for(unsigned i=0;i<count;++i) device.p[i]=calculate_point(base,step,i,target,giant);
#else
  point_kernel<<<(count+127)/128,128>>>(base,step,count,device.p,target,giant?1:0);
#endif
  CUDA(cudaGetLastError());
  CUDA(cudaMemcpy(out.data(),device.p,count*sizeof(Affine),cudaMemcpyDeviceToHost));
}
static std::vector<Affine> load_targets(const std::string &path) {
  std::ifstream in(path); if(!in) throw std::runtime_error("cannot open target file: "+path);
  std::vector<Affine> result; std::string line; unsigned n=0;
  while(std::getline(in,line)) { ++n; auto pos=line.find('#'); if(pos!=std::string::npos) line.resize(pos); std::istringstream s(line); std::string token, extra; if(!(s>>token)) continue; if(s>>extra) throw std::runtime_error("unexpected extra data at line "+std::to_string(n)); try { result.push_back(parse_public(token)); } catch(const std::exception &e) { throw std::runtime_error("line "+std::to_string(n)+": "+e.what()); } }
  if(in.bad()) throw std::runtime_error("failed reading target file");
  if(result.empty()) throw std::runtime_error("target file is empty");
  return result;
}
static void usage() { std::cout<<"Usage: keyhunt-cuda -m xpoint|bsgs -f publickeys.txt -r START:END [--batch N]\n"
  "START and END are inclusive hexadecimal private keys; one SEC public key per line.\n"; }
static unsigned parse_batch(const std::string &s) { if(s.empty() || s.find_first_not_of("0123456789")!=std::string::npos) throw std::runtime_error("invalid batch size"); unsigned long long n=std::stoull(s); if(!n || n>1048576) throw std::runtime_error("batch size must be 1..1048576"); return (unsigned)n; }
static uint64_t small_difference(const U256 &a, const U256 &b, uint64_t limit) {
  U256 d=sub_raw(a,b); for(int i=2;i<8;++i) if(d.v[i]) throw std::runtime_error("BSGS range exceeds maximum of 10^12 keys");
  uint64_t value=((uint64_t)d.v[1]<<32)|d.v[0]; if(value>=limit) throw std::runtime_error("BSGS range exceeds maximum of 10^12 keys"); return value;
}
static unsigned count_chunk(const U256 &cur,const U256 &end,unsigned batch) {
  U256 d=sub_raw(end,cur); for(int i=2;i<8;++i) if(d.v[i]) return batch;
  uint64_t n=((uint64_t)d.v[1]<<32)|d.v[0]; return n>=batch-1?batch:(unsigned)(n+1);
}
static void print_hit(const U256 &k,const Affine &p) { std::cout<<"Hit! Private Key: "<<hex256(k)<<"\nPublic key: 04"<<hex256(p.x)<<hex256(p.y)<<"\n"<<std::flush; }
static void xpoint(const std::vector<Affine> &targets,U256 start,const U256 &end,unsigned batch) {
  std::vector<U256> xs; for(const auto &t:targets) xs.push_back(t.x);
  auto less=[](const U256 &a,const U256 &b){return cmp(a,b)<0;}; std::sort(xs.begin(),xs.end(),less);
  DeviceBuffer device(batch); std::vector<Affine> out(batch); Affine unused={}; U256 cur=start;
  while(true) {
    unsigned n=count_chunk(cur,end,batch); points(device,out,cur,1,n,unused,false);
    for(unsigned i=0;i<n;++i) if(!out[i].infinity && std::binary_search(xs.begin(),xs.end(),out[i].x,less)) {
      U256 k=add_raw(cur,from64(i)); Affine verified=affine(multiply(k));
      if(point_cmp(verified,out[i])!=0) throw std::runtime_error("GPU/CPU verification failed");
      print_hit(k,verified);
    }
    U256 last=add_raw(cur,from64(n-1)); if(cmp(last,end)==0) break; cur=add_raw(last,one());
  }
}
struct Baby { Affine p; uint64_t j; };
static void bsgs(const std::vector<Affine> &targets,U256 start,const U256 &end,unsigned batch) {
  const uint64_t length=small_difference(end,start,1000000000000ULL)+1;
  uint64_t m=(uint64_t)std::sqrt((long double)length); while(m*m<length) ++m; while(m>1 && (m-1)*(m-1)>=length) --m;
  std::cerr<<"BSGS: "<<length<<" keys, "<<m<<" baby points, "<<targets.size()<<" target(s)\n";
  DeviceBuffer device(batch); std::vector<Affine> out(batch); std::vector<Baby> babies; babies.reserve((size_t)m); Affine unused={};
  for(uint64_t j=0;j<m;) { unsigned n=(unsigned)std::min<uint64_t>(batch,m-j); points(device,out,from64(j),1,n,unused,false); for(unsigned i=0;i<n;++i) babies.push_back({out[i],j+i}); j+=n; }
  std::sort(babies.begin(),babies.end(),[](const Baby &a,const Baby &b){return point_cmp(a.p,b.p)<0;});
  uint64_t giants=(length+m-1)/m;
  for(const Affine &target:targets) {
    bool found=false;
    for(uint64_t i=0;i<giants && !found;) {
      unsigned n=(unsigned)std::min<uint64_t>(batch,giants-i);
      points(device,out,add_raw(start,from64(i*m)),m,n,target,true);
      for(unsigned t=0;t<n;++t) {
        auto it=std::lower_bound(babies.begin(),babies.end(),out[t],[](const Baby &b,const Affine &p){return point_cmp(b.p,p)<0;});
        if(it!=babies.end() && point_cmp(it->p,out[t])==0) {
          uint64_t offset=(i+t)*m+it->j; if(offset>=length) continue;
          U256 key=add_raw(start,from64(offset)); Affine verified=affine(multiply(key));
          if(point_cmp(verified,target)!=0) throw std::runtime_error("GPU/CPU verification failed");
          print_hit(key,verified); found=true; break;
        }
      }
      i+=n;
    }
    if(!found) std::cout<<"No match for public key: "<<(target.y.v[0]&1U?"03":"02")<<hex256(target.x)<<"\n";
  }
}
int main(int argc,char **argv) {
  try {
    std::string mode,file,range; unsigned batch=256;
    for(int i=1;i<argc;++i) { std::string a=argv[i]; if(a=="-h"||a=="--help") {usage(); return 0;} if(a!="-m"&&a!="-f"&&a!="-r"&&a!="--batch") throw std::runtime_error("unknown option: "+a); if(++i>=argc) throw std::runtime_error("missing argument for "+a); if(a=="-m") mode=argv[i]; else if(a=="-f") file=argv[i]; else if(a=="-r") range=argv[i]; else batch=parse_batch(argv[i]); }
    if((mode!="xpoint"&&mode!="bsgs")||file.empty()) throw std::runtime_error("specify -m xpoint|bsgs and -f");
    auto p=range.find(':'); if(p==std::string::npos || range.find(':',p+1)!=std::string::npos) throw std::runtime_error("range must be START:END");
    U256 start=parse_hex(range.substr(0,p)), end=parse_hex(range.substr(p+1));
    if(iszero(start)||cmp(start,end)>0||cmp(end,order())>=0) throw std::runtime_error("range must satisfy 1 <= START <= END < secp256k1 order");
    auto targets=load_targets(file); CUDA(cudaSetDevice(0));
    if(mode=="xpoint") xpoint(targets,start,end,batch); else bsgs(targets,start,end,batch);
    CUDA(cudaDeviceSynchronize()); return 0;
  } catch(const std::exception &e) { std::cerr<<"error: "<<e.what()<<"\n"; return 1; }
}
