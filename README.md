# ML-KEM-768

This is a from-scratch Rust implementation of ML-KEM-768 (NIST's post-quantum
key encapsulation standard) in ~1,900 lines.

The layers of the implementation are:

```
Zq (Montgomery arithmetic)
  → Poly (ring Zq[X]/(X²⁵⁶+1) with NTT)
    → Vec3 / Mat3x3 (module lattice algebra)
      → KPKE (Kyber public-key encryption)
        → ML-KEM-768 (IND-CCA2 KEM)
```

## License

Licensed under either of

 * Apache License, Version 2.0, (http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license (http://opensource.org/licenses/MIT)
