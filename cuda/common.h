#pragma once
#include "curve.h"
#include <string>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cctype>

inline U256 parse_hex(const std::string &s) {
  if(s.empty() || s.size()>64) throw std::runtime_error("expected 1 to 64 hexadecimal digits");
  U256 a=zero();
  for(char c:s) { int d=c>='0'&&c<='9'?c-'0':c>='a'&&c<='f'?c-'a'+10:c>='A'&&c<='F'?c-'A'+10:-1; if(d<0) throw std::runtime_error("invalid hexadecimal value: "+s); for(int i=7;i>0;--i) a.v[i]=(a.v[i]<<4)|(a.v[i-1]>>28); a.v[0]=(a.v[0]<<4)|d; }
  return a;
}
inline std::string hex256(const U256 &a) { std::ostringstream s; s<<std::hex<<std::setfill('0'); for(int i=7;i>=0;--i) s<<std::setw(8)<<a.v[i]; return s.str(); }
inline Affine parse_public(const std::string &s) {
  if(s.size()==66 && (s.substr(0,2)=="02" || s.substr(0,2)=="03")) { Affine a; if(!decompress(parse_hex(s.substr(2)),s[1]=='3',a)) throw std::runtime_error("invalid compressed public key"); return a; }
  if(s.size()==130 && s.substr(0,2)=="04") { Affine a={parse_hex(s.substr(2,64)),parse_hex(s.substr(66,64)),false}; if(!oncurve(a)) throw std::runtime_error("invalid uncompressed public key"); return a; }
  throw std::runtime_error("expected a compressed or uncompressed SEC public key");
}
inline int point_cmp(const Affine &a,const Affine &b) { if(a.infinity!=b.infinity) return a.infinity?-1:1; if(a.infinity) return 0; int c=cmp(a.x,b.x); return c?c:cmp(a.y,b.y); }
