// Correctness-oriented secp256k1 primitives shared by the CUDA kernels and CPU tests.
// Not constant time. Not suitable for signing or handling secret keys.
#pragma once
#include <stdint.h>
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

struct U256 { uint32_t v[8]; }; // little-endian limbs
struct Point { U256 x, y, z; }; // Jacobian; z == 0 denotes infinity
struct Affine { U256 x, y; bool infinity; };

HD inline U256 zero() { U256 a = {}; return a; }
HD inline U256 one() { U256 a = {}; a.v[0]=1; return a; }
HD inline U256 prime() { U256 a={{0xfffffc2f,0xfffffffe,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff}}; return a; }
HD inline U256 order() { U256 a={{0xd0364141,0xbfd25e8c,0xaf48a03b,0xbaaedce6,0xfffffffe,0xffffffff,0xffffffff,0xffffffff}}; return a; }
HD inline Affine generator() {
  Affine g={{{0x16f81798,0x59f2815b,0x2dce28d9,0x029bfcdb,0xce870b07,0x55a06295,0xf9dcbbac,0x79be667e}},
            {{0xfb10d4b8,0x9c47d08f,0xa6855419,0xfd17b448,0x0e1108a8,0x5da4fbfc,0x26a3c465,0x483ada77}},false}; return g;
}
HD inline int cmp(const U256 &a,const U256 &b) { for(int i=7;i>=0;--i) if(a.v[i]!=b.v[i]) return a.v[i]<b.v[i]?-1:1; return 0; }
HD inline bool iszero(const U256 &a) { for(int i=0;i<8;++i) if(a.v[i]) return false; return true; }
HD inline U256 sub_raw(const U256 &a,const U256 &b) { U256 r; uint64_t borrow=0; for(int i=0;i<8;++i) { uint64_t x=(uint64_t)b.v[i]+borrow; r.v[i]=(uint32_t)((uint64_t)a.v[i]-x); borrow=(uint64_t)a.v[i]<x; } return r; }
HD inline U256 add_raw(const U256 &a,const U256 &b) { U256 r; uint64_t carry=0; for(int i=0;i<8;++i) { uint64_t x=(uint64_t)a.v[i]+b.v[i]+carry; r.v[i]=(uint32_t)x; carry=x>>32; } return r; }
HD inline U256 from64(uint64_t v) { U256 a=zero(); a.v[0]=(uint32_t)v; a.v[1]=(uint32_t)(v>>32); return a; }
// Inputs are reduced. a + b mod p without overflowing 256 bits.
HD inline U256 fadd(const U256 &a,const U256 &b) { U256 t=sub_raw(prime(),b); return cmp(a,t)>=0 ? sub_raw(a,t) : add_raw(a,b); }
HD inline U256 fsub(const U256 &a,const U256 &b) { return cmp(a,b)>=0 ? sub_raw(a,b) : sub_raw(prime(),sub_raw(b,a)); }
HD inline U256 fmul(U256 a,const U256 &b) { U256 r=zero(); for(int i=0;i<256;++i) { if((b.v[i/32]>>(i%32))&1U) r=fadd(r,a); a=fadd(a,a); } return r; }
HD inline U256 fpow(const U256 &a,const U256 &e) { U256 r=one(); for(int i=255;i>=0;--i) { r=fmul(r,r); if((e.v[i/32]>>(i%32))&1U) r=fmul(r,a); } return r; }
HD inline U256 inv(const U256 &a) { U256 e=prime(); e.v[0]-=2; return fpow(a,e); }
HD inline Point infinity() { Point p={}; return p; }
HD inline Point from_affine(const Affine &a) { if(a.infinity) return infinity(); Point p={a.x,a.y,one()}; return p; }
HD inline Point twice(const Point &p) {
  if(iszero(p.z)||iszero(p.y)) return infinity();
  U256 a=fmul(p.x,p.x), b=fmul(p.y,p.y), c=fmul(b,b);
  U256 d=fsub(fsub(fmul(fadd(p.x,b),fadd(p.x,b)),a),c); d=fadd(d,d);
  U256 e=fadd(fadd(a,a),a), f=fmul(e,e);
  Point r; r.x=fsub(f,fadd(d,d));
  U256 eight=fadd(c,c); eight=fadd(eight,eight); eight=fadd(eight,eight);
  r.y=fsub(fmul(e,fsub(d,r.x)),eight); r.z=fmul(fadd(p.y,p.y),p.z); return r;
}
HD inline Point add_mixed(const Point &p,const Affine &q) {
  if(q.infinity) return p;
  if(iszero(p.z)) return from_affine(q);
  U256 zz=fmul(p.z,p.z), u=fmul(q.x,zz), s=fmul(q.y,fmul(zz,p.z));
  U256 h=fsub(u,p.x), rr=fsub(s,p.y);
  if(iszero(h)) return iszero(rr)?twice(p):infinity();
  U256 hh=fmul(h,h), hhh=fmul(h,hh), v=fmul(p.x,hh);
  Point r; r.x=fsub(fsub(fmul(rr,rr),hhh),fadd(v,v));
  r.y=fsub(fmul(rr,fsub(v,r.x)),fmul(p.y,hhh)); r.z=fmul(p.z,h); return r;
}
HD inline Point multiply(const U256 &k) {
  Point p=infinity(); Affine g=generator();
  for(int i=255;i>=0;--i) { p=twice(p); if((k.v[i/32]>>(i%32))&1U) p=add_mixed(p,g); }
  return p;
}
HD inline Affine affine(const Point &p) {
  Affine a={}; if(iszero(p.z)) { a.infinity=true; return a; }
  U256 zi=inv(p.z), zz=fmul(zi,zi); a.x=fmul(p.x,zz); a.y=fmul(p.y,fmul(zz,zi)); return a;
}
HD inline Affine negate(Affine p) { if(!p.infinity && !iszero(p.y)) p.y=sub_raw(prime(),p.y); return p; }
HD inline bool oncurve(const Affine &p) { if(p.infinity || cmp(p.x,prime())>=0 || cmp(p.y,prime())>=0) return false; return cmp(fmul(p.y,p.y),fadd(fmul(fmul(p.x,p.x),p.x),from64(7)))==0; }
// p = 3 (mod 4), so square roots are a^((p+1)/4).
HD inline bool decompress(const U256 &x, bool odd, Affine &out) {
  if(cmp(x,prime())>=0) return false;
  U256 e={{0xbfffff0c,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0x3fffffff}};
  U256 rhs=fadd(fmul(fmul(x,x),x),from64(7)); U256 y=fpow(rhs,e);
  if(cmp(fmul(y,y),rhs)!=0) return false;
  if((y.v[0]&1U)!=(unsigned)odd) { if(iszero(y)) return false; y=sub_raw(prime(),y); }
  out={x,y,false}; return true;
}
