# Experimental CUDA XPOINT and BSGS

`keyhunt-cuda` is a **separate, experimental** program. It does not replace
`keyhunt` and implements only sequential XPOINT and BSGS searches. It is not
command-line or feature compatible with all options in the CPU program.
Use it only for puzzles or keys you are authorized to recover.

## Requirements and build

Linux, an NVIDIA GPU and driver compatible with the CUDA toolkit, `nvcc`, a
compatible C++ compiler, `make`, and Python 3 for the tests. Build with:

```sh
make cuda
# Example for a particular GPU architecture, if needed:
make cuda NVCCFLAGS='-O2 -std=c++14 -arch=sm_80'
```

This builds `keyhunt-cuda`; the normal `make` continues to build the existing
CPU program. To rebuild with different flags, remove `keyhunt-cuda` first.

## Input and usage

Supply one valid compressed (02/03) or uncompressed (04) SEC secp256k1 public
key per line. Blank lines and `#` comments are allowed. The whole file is
validated before searching. For example:

```sh
printf '%s\n' 0279be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798 > /tmp/targets.txt
./keyhunt-cuda -m xpoint -f /tmp/targets.txt -r 1:10 --batch 256
./keyhunt-cuda -m bsgs   -f /tmp/targets.txt -r 1:10 --batch 256
```

`-r START:END` is required: **both endpoints are inclusive hexadecimal**
integers, between 1 and the secp256k1 group order minus 1. This is intentionally
explicit and should not be assumed to have the same endpoint semantics as
CPU keyhunt. There is no random search, bit-range shortcut, endomorphism,
stride, checkpoint/resume, multi-GPU selection, bloom filter, file caching,
statistics, or automatic result file. `--batch` is the number of GPU points
per launch (1 to 1,048,576; default 256). GPU device 0 is used. Hits are
printed to standard output, including private keys: protect that output.
Exit status is 0 for a completed search, even when no key matches; invalid
input, CUDA failures and failed hit verification exit nonzero. Interrupted
runs must be restarted.

XPOINT compares only the x-coordinate, so a key and its negation have the same
x-coordinate; public-key parity in the input does not restrict the match.
All matching keys in the interval are printed once, even if targets repeat.
BSGS uses the complete public key, so parity matters. It searches each target
independently using a shared baby-step table and prints at most one hit per
input target; duplicates may therefore produce duplicate hits. A missing BSGS
target is reported. Both modes validate hits again on the CPU.

## Design and limitations

Each CUDA thread independently multiplies a scalar by the generator; the CPU
sorts and matches the resulting points. BSGS builds `ceil(sqrt(length))` baby
points and evaluates `ceil(length / m)` giant points per target, using exact
full-point comparison and rejecting padded results outside the range. The
BSGS interval is limited to 10^12 keys, requiring up to one million baby
points in host memory. XPOINT has no such interval-size limit, but a large
search may take an impractical amount of time. This reference implementation
uses slow, simple finite-field arithmetic, transfers every generated point to
the host, and performs no GPU tuning. It is expected to be **much slower**
than optimized keyhunt implementations; GPU watchdog timeouts are possible.
It has not been benchmarked or validated on real GPU hardware and should not
be relied on for production searches without further independent verification.
CPU verification of hits does not prove that a GPU run missed no matches.

## Tests

```sh
make cuda-test                      # CPU differential + simulated integration tests
make cuda
CUDA_TEST_GPU=1 make cuda-test      # additionally exercise real CUDA kernels
```

The tests compare field arithmetic and scalar multiplication against an
independent Python big-integer implementation and exercise ranges, full-point
and x-only semantics, input validation, batch boundaries, group-order edges,
and BSGS padding. The CPU simulation exercises the same point calculations and
host search logic, but is **not** a substitute for CUDA execution on a GPU.
Real GPU tests must be run on a compatible NVIDIA system before relying on the
program. The arithmetic is not constant-time and must not be used for signing
or other secret-key cryptography.
