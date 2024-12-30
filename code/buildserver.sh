 #!/bin/bash

 rm -Rf kyber
 git clone https://github.com/pq-crystals/kyber
 cd kyber/ref
 git checkout 10b478fc3cc4ff6215eb0b6a11bd758bf0929cbd
 cd ../
 git apply ../patch.patch
 cd ref
 gcc -shared -fPIC -DKYBER_K=2 randombytes.c fips202.c symmetric-shake.c indcpa.c polyvec.c poly.c ntt.c cbd.c reduce.c verify.c kem.c -o libpqcrystals_kyber512_ref_ch.so
 cd ../../
 cp kyber/ref/libpqcrystals_kyber512_ref_ch.so libpqcrystals_kyber512_ref_ch.so
 gcc server_timing.c -o server_timing -L. -l:libpqcrystals_kyber512_ref_ch.so -Wl,-R. inverse.c
 # replace the above build command with the following one if you are running a Macos
 # gcc server_timing.c -o server_timing -L. -lpqcrystals_kyber512_ref_ch inverse.c
 rm -Rf kyber
