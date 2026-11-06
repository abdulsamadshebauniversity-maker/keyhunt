#!/usr/bin/env python3
"""Differential tests for code shared with the kernels; optional real GPU smoke tests."""
import os
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
P = 2**256 - 2**32 - 977
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
G = (0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
     0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8)

def add(a, b):
    if a is None: return b
    if b is None: return a
    x, y = a
    u, v = b
    if x == u and (y + v) % P == 0: return None
    slope = ((3*x*x if a == b else v-y) * pow((2*y if a == b else u-x) % P, -1, P)) % P
    z = (slope*slope-x-u) % P
    return z, (slope*(x-z)-y) % P

def mul(k):
    r, a = None, G
    while k:
        if k & 1: r = add(r, a)
        a = add(a, a)
        k >>= 1
    return r

def sec(p, compressed=True):
    x, y = p
    return (f'{2+(y&1):02x}{x:064x}' if compressed else f'04{x:064x}{y:064x}')

def main():
    rng = random.Random(49611)
    commands, expected = [], []
    for a, b in [(0, 0), (P-1, P-1), (P-1, 1), (1, P-1)] + [(rng.randrange(P),rng.randrange(P)) for _ in range(30)]:
        commands.append(f'field {a:x} {b:x}')
        expected.append(' '.join(f'{x:064x}' for x in ((a+b)%P, (a-b)%P, a*b%P)))
    for k in [0,1,2,3,7,17,255,256,65537,N-1,N] + [rng.randrange(1,N) for _ in range(8)]:
        commands.append(f'point {k:x}')
        p=mul(k)
        expected.append('infinity' if p is None else f'{p[0]:064x} {p[1]:064x}')
    for k in [1,2,3,N-1] + [rng.randrange(1,N) for _ in range(4)]:
        p=mul(k)
        for compressed in [True,False]:
            commands.append('pub '+sec(p,compressed)); expected.append(f'{p[0]:064x} {p[1]:064x}')
    commands.extend(['pub 02'+'ff'*32, 'pub 04'+'00'*64, 'pub 05'+'00'*32, 'pub 02'+'00'*32])
    expected.extend(['invalid']*4)
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d)/'curve-test'
        subprocess.run(['g++','-O2','-std=c++14','-Wall','-Wextra','-o',str(exe),str(ROOT/'cuda/test_driver.cpp')],check=True)
        run = subprocess.run([str(exe)],input='\n'.join(commands)+'\n',text=True,capture_output=True,check=True)
        actual=run.stdout.splitlines()
        assert len(actual)==len(expected), (len(actual),len(expected),run.stderr)
        for cmd, got, want in zip(commands,actual,expected):
            assert got==want, f'{cmd}: got {got}, expected {want}'
        print(f'CPU differential tests passed: {len(expected)} cases')
        def integration(binary, label):
            targets=Path(d)/'targets.txt'
            targets.write_text(sec(mul(1))+'\n'+sec(mul(3),False)+' # comment\n')
            def run_search(mode, start, end, batch=2, file=targets):
                return subprocess.run([str(binary),'-m',mode,'-f',str(file),'-r',f'{start:x}:{end:x}','--batch',str(batch)],text=True,capture_output=True)
            for mode in ('xpoint','bsgs'):
                r=run_search(mode,1,5)
                assert r.returncode==0, r.stderr
                assert r.stdout.count('Hit! Private Key:')==2, r.stdout
                for k in (1,3): assert f'{k:064x}' in r.stdout, r.stdout
                r=run_search(mode,4,5)
                assert r.returncode==0 and 'Hit!' not in r.stdout, (r.stdout,r.stderr)
                r=run_search(mode,0,1)
                assert r.returncode!=0, 'invalid range accepted'
                r=run_search(mode,2,1)
                assert r.returncode!=0, 'reversed range accepted'
            targets.write_text(sec(mul(N-1))+'\n')
            for mode in ('xpoint','bsgs'):
                r=run_search(mode,N-2,N-1,1)
                assert r.returncode==0 and f'{N-1:064x}' in r.stdout, (r.stdout,r.stderr)
            # XPOINT matches either sign; BSGS matches the full public key.
            r=run_search('xpoint',1,1); assert r.returncode==0 and 'Hit!' in r.stdout, (r.stdout,r.stderr)
            r=run_search('bsgs',1,1); assert r.returncode==0 and 'Hit!' not in r.stdout, (r.stdout,r.stderr)
            # A padded BSGS giant step must not return a key outside the interval.
            targets.write_text(sec(mul(6))+'\n')
            r=run_search('bsgs',1,5); assert r.returncode==0 and 'Hit!' not in r.stdout, (r.stdout,r.stderr)
            targets.write_text(sec(mul(5))+'\n')
            r=run_search('bsgs',1,5); assert r.returncode==0 and f'{5:064x}' in r.stdout, (r.stdout,r.stderr)
            # Ranges with large high limbs and a non-square number of keys.
            start=2**100+123
            targets.write_text(sec(mul(start+10))+'\n')
            for mode in ('xpoint','bsgs'):
                r=run_search(mode,start,start+12,3)
                assert r.returncode==0 and f'{start+10:064x}' in r.stdout, (r.stdout,r.stderr)
            r=run_search('bsgs',1,10**12+1)
            assert r.returncode!=0, 'oversized BSGS range accepted'
            # Exercise each BSGS remainder class and the infinity baby point.
            for low, high, key in [(10,10,10),(10,18,10),(10,18,13),(10,18,18),(10,18,19)]:
                targets.write_text(sec(mul(key))+'\n')
                r=run_search('bsgs',low,high,2)
                assert r.returncode==0, r.stderr
                assert ('Hit!' in r.stdout)==(key<=high), (low,high,key,r.stdout)
                if key<=high: assert f'{key:064x}' in r.stdout, r.stdout
            targets.write_text('not-a-key\n')
            r=run_search('xpoint',1,1); assert r.returncode!=0, 'invalid target accepted'
            print(label+' integration tests passed')
        host=Path(d)/'host-search'
        subprocess.run(['g++','-x','c++','-DKEYHUNT_CPU_TEST','-O2','-std=c++14','-Wall','-Wextra','-o',str(host),str(ROOT/'cuda/search.cu')],check=True)
        integration(host, 'CPU-simulated kernel')
        if os.environ.get('CUDA_TEST_GPU') == '1':
            integration(ROOT/'keyhunt-cuda', 'GPU')
        else:
            print('GPU integration tests skipped (set CUDA_TEST_GPU=1 after building on an NVIDIA host)')

if __name__=='__main__': main()
