// CPU harness for the exact field and curve code used in the CUDA kernels.
#include "common.h"
#include <iostream>
int main() {
  std::string op,a,b;
  while(std::cin>>op>>a) {
    try {
      if(op=="point") { Affine p=affine(multiply(parse_hex(a))); if(p.infinity) std::cout<<"infinity\n"; else std::cout<<hex256(p.x)<<" "<<hex256(p.y)<<"\n"; }
      else if(op=="pub") { Affine p=parse_public(a); std::cout<<hex256(p.x)<<" "<<hex256(p.y)<<"\n"; }
      else if(op=="field") { std::cin>>b; U256 x=parse_hex(a), y=parse_hex(b); std::cout<<hex256(fadd(x,y))<<" "<<hex256(fsub(x,y))<<" "<<hex256(fmul(x,y))<<"\n"; }
      else return 2;
    } catch(const std::exception &) { std::cout<<"invalid\n"; }
  }
}
