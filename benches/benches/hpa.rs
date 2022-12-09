use ark_bls12_381::Bls12_381;
use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
    // pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};
use ark_ec::{PairingEngine};
use ark_ff::{UniformRand};
use ark_inner_products::{
    ExtensionFieldElement, InnerProduct, PairingInnerProduct,
};
use ark_hpa::hpa::HPA;

use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
use blake2::Blake2b;
use digest::Digest;

use std::{ops::MulAssign, time::Instant};

fn bench_hpa<IP, LMC, RMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective, Message = P::G1Projective>,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar, Key = P::G1Projective, Message = P::G2Projective>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    LMC::Message: MulAssign<P::Fr>,
    RMC::Message: MulAssign<P::Fr>,
    IPC::Message: MulAssign<P::Fr>,
    IPC::Key: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
    RMC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<LMC::Scalar>,
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
    LMC::Output: MulAssign<LMC::Scalar>,
    // IPC::Message: AddAssign<LMC::Output>,
    // IPC::Message: AddAssign<RMC::Output>,
    // RMC::Output: AddAssign<LMC::Output>,
{
    let mut l = Vec::new(); 
    let mut r = Vec::new();
  
    for _ in 0..len {
        l.push(<IP::LeftMessage>::rand(rng));
        r.push(<IP::RightMessage>::rand(rng));
    }

    let (gamma2, gamma1) = HPA::<IP,LMC,RMC,IPC, D>::setup(rng, len).unwrap();
    let d1 = IP::inner_product(&l, &gamma2).unwrap();
    let d2 = IP::inner_product(&gamma1, &r).unwrap();
    let c = IP::inner_product(&l, &r).unwrap();

    let mut hpa_srs = HPA::<IP, LMC, RMC, IPC, D>::precompute((&(gamma1.clone()), &(gamma2.clone()))).unwrap();
    let mut start = Instant::now();
    let mut proof =
        HPA::<IP, LMC, RMC, IPC, D>::prove((&(l.clone()), &(r.clone())),
        //  (&(gamma1.clone()), &(gamma2.clone())), 
         (&(gamma1.clone()), &(gamma2.clone())), 
         (&(d1.clone()), &(d2.clone()), &(c.clone()))
        ).unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);
    start = Instant::now();
    let result = HPA::<IP, LMC, RMC, IPC, D>::verify(&mut hpa_srs, (&(gamma1.clone()), &(gamma2.clone())),
         (&(d1.clone()), &(d2.clone()), &(c.clone())), &mut proof)
        .unwrap();
    bench = start.elapsed().as_millis();
    println!("\t verification time: {} ms", bench);
    println!("result : {}", result);
}


fn main() {
    const LEN: usize = 16;
    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    let mut rng = StdRng::seed_from_u64(0u64);

    println!("Benchmarking HPA with vector length: {}", LEN);

    println!("1) Pairing inner product...");
    bench_hpa::<
        PairingInnerProduct<Bls12_381>,
        GC1,
        GC2,
        IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, LEN);

}
