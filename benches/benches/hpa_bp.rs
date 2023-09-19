use ark_bls12_381::Bls12_381;
use ark_dh_commitments::{
    // afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    // identity::IdentityCommitment,
    pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};
use ark_ec::PairingEngine;
use ark_ff::UniformRand;
use ark_inner_products::{
    InnerProduct, ScalarInnerProduct//, PairingInnerProduct, //ExtensionFieldElement
};
use ark_hpa_bp::hpa_bp::HPA;

use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
use blake2::Blake2b;
use digest::Digest;

use std::{time::Instant, env}; //, env,ops::MulAssign};


// fn bench_hpa_bp<IP, LMC, RMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
fn bench_hpa_bp<IP, CM, P, D, R: Rng>(rng: &mut R, len: usize)
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = CM::Message,
        RightMessage = CM::Message,
        Output = CM::Message,
    >,
    CM: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G1Projective, Message = P::Fr>,
    // LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G1Projective, Message = P::Fr>,
    // RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar, Key = P::G1Projective, Message = P::Fr>,
    // IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // LMC::Message: MulAssign<P::Fr>,
    // RMC::Message: MulAssign<P::Fr>,
    // IPC::Message: MulAssign<P::Fr>,
    // IPC::Key: MulAssign<P::Fr>,
    // LMC::Output: MulAssign<P::Fr>,
    // RMC::Output: MulAssign<P::Fr>,
    // IPC::Output: MulAssign<P::Fr>,
    // IPC::Output: MulAssign<LMC::Scalar>,
    // IP::LeftMessage: UniformRand,
    // IP::RightMessage: UniformRand,
    // LMC::Output: MulAssign<LMC::Scalar>,
{
    let mut v1 = Vec::new(); 
    let mut v2 = Vec::new();

    for _ in 0..len {
        v1.push(<CM as DoublyHomomorphicCommitment>::Scalar::rand(rng));
        v2.push(<CM as DoublyHomomorphicCommitment>::Scalar::rand(rng));
    }
    // let generator_g1 = <IP::LeftMessage>::rand(rng);
    // let generator_g2 = <IP::RightMessage>::rand(rng);
    
    let (u1, u2) = HPA::<IP, CM, D>::set_values(&v1, &v2).unwrap();

    let (gamma1, gamma2)= HPA::<IP,CM, D>::setup(rng, len).unwrap();
    
    
    // let mut h1 = Vec::new();
    // let mut h2 = Vec::new();
    // h1.push(<IP::LeftMessage>::rand(rng));
    // h2.push(<IP::RightMessage>::rand(rng));

    let (c, d1, d2, 
        x, d3, 
        gm, gm_vec, 
        r_c, r_x, w_vec)
        = HPA::<IP,CM, D>::init_commit(&v1, &v2, &gamma1, &gamma2, rng).unwrap();

    let (c_, d1_, d2_, 
        x_, d3_, r_c_, w_vec_)
        = HPA::<IP,CM, D>::init_commit2(&u1, &u2, &gamma1, &gamma2,  &gm_vec, &r_x, rng).unwrap();

    // X ?= X'
    let bool_x = x == x_;
    println!("X == X' : {}", bool_x);   


    // let mut hpa_srs = HPA::<IP, LMC, RMC, IPC, D>::precompute((&(gamma1.clone()), &(gamma2.clone())), &h1, &h2).unwrap();
    // let mut hpa_srs_ = HPA::<IP, LMC, RMC, IPC, D>::precompute((&(gamma1.clone()), &(gamma2.clone())), &h1, &h2).unwrap();


    let mut start = Instant::now();
    let mut proof =
        HPA::<IP, CM, D>::prove((&(v1.clone()), &(v2.clone()), &(w_vec.clone())),
        //  &hpa_srs, 
         (&(gamma1.clone()), &(gamma2.clone())), 
        //  (&(d1.clone()), &(d2.clone()), &(c.clone())),
         (&r_c, &r_x),// &r_y, &r_d1, &r_d2, &r_d3, &r_d4),
         &gm,
         rng
        ).unwrap();

    // let mut start = Instant::now();
    let mut proof_ =
        HPA::<IP, CM, D>::prove((&(u1.clone()), &(u2.clone()), &(w_vec_.clone())),// &(k_vec_.clone())),
            // &hpa_srs_, 
            (&(gamma1.clone()), &(gamma2.clone())), 
            //  (&(d1.clone()), &(d2.clone()), &(c.clone())),
            (&r_c_, &r_x), // &r_y, &r_d1, &r_d2, &r_d3, &r_d4),
            &gm,
            rng
        ).unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);


    start = Instant::now();
    let result = HPA::<IP,CM, D>::verify( (gamma1.clone(), gamma2.clone()),
        (&c, &x, &d1, &d2, &d3), &mut proof, &gm)
        .unwrap();
    let result2 = HPA::<IP, CM, D>::verify( (gamma1.clone(), gamma2.clone()),
    (&c_, &x_, &d1_, &d2_, &d3_), &mut proof_, &gm)
    .unwrap();
    bench = start.elapsed().as_millis();
    println!("\t verification time: {} ms", bench);

    println!("v1, v2 - result : {}", result);
    println!("u1, u2 - result : {}", result2);
}


fn main() { 
    let arg = env::args().nth(1).unwrap();
    let power: u32 =arg.parse().unwrap();
    
    let len = usize::pow(2,power);
    // const LEN: usize = 32;
    type CM = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
    // type GC1 = AFGHOCommitmentG1<Bls12_381>;
    // type GC2 = AFGHOCommitmentG2<Bls12_381>;
    let mut rng = StdRng::seed_from_u64(0u64);

    println!("Benchmarking HPA_with_zk with vector length: {}", len);

    println!("1) Pairing hadamard product...");
    bench_hpa_bp::<
        // InnerProduct<Bls12_381>,
        ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>,
        // GC1,
        // GC2,
        CM,
        // IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, len);

}
