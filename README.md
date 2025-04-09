<h1 align="center">RHPA (Rust Hadamard Products Arguments)</h1>

<p align="center">
    <a href="https://travis-ci.org/scipr-lab/ripp"><img src="https://travis-ci.org/scipr-lab/ripp.svg?branch=master"></a>
    <a href="https://github.com/scipr-lab/ripp/blob/master/AUTHORS"><img src="https://img.shields.io/badge/authors-SCIPR%20Lab-orange.svg"></a>
    <a href="https://github.com/scipr-lab/ripp/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/scipr-lab/ripp/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/scipr-lab/ripp"><img src="https://deps.rs/repo/github/scipr-lab/ripp/status.svg"></a>
</p>

## HPA

___HPA___ is a Rust library implementing Hadamard Product Arguments (HPA) based on Inner Pairing Products (IPP). The underlying HPA protocols and the HPA application built atop these are described in our paper *"[Hadamard Product Arguments and Their Applications][hpa]"*.

The library currently contains an implementation of our HPA mechanism utilizing IPP, enabling the aggregation of certain types of proofs. In the future, we intend to expand support for aggregating other proof systems (such as Garuda and Pari), implement related polynomial commitment schemes, and potentially other protocols leveraging HPA, as discussed in our [paper][hpa].

This library is released under the MIT License and the Apache v2 License (see [License](#license)).

[hpa]: https://eprint.iacr.org/2024/981

## License


This project implements HPA features by modifying and extending the ___RIPP___ library (see [https://github.com/arkworks-rs/ripp]). Please refer to the original ___RIPP___ project for its copyright and license information.

This ___HPA___ library itself is distributed under the terms of both the MIT License and the Apache License (Version 2.0). See LICENSE-MIT and LICENSE-APACHE files in the repository for details.

This library is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/5Donghwan/HPA.git
cd HPA
cargo build --release
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test
``` 

Lastly, the library comes with benchmarks.
```bash
cargo bench
cargo run --release --example groth16_aggregation
cargo run --release --example scaling-ipp
```

## License

RIPP is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in RIPP by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[ripp]: https://eprint.iacr.org/2019/1177

## Reference paper

[_Proofs for Inner Pairing Products and Applications_][ripp]    
[Benedikt Bünz](https://www.github.com/bbuenz), Mary Maller, [Pratyush Mishra](https://www.github.com/pratyush), [Psi Vesely](https://www.github.com/psivesely)    
*IACR ePrint Report 2019/1177*
