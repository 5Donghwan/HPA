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
use ark_dory::dory::{
    DORY,
};
use ark_mv_product::mv_product::MVP;

use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
use blake2::Blake2b;
use digest::Digest;

use std::{ops::MulAssign, time::Instant};

fn bench_mvp<IP, LMC, RMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
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
    let (gamma2, gamma1) = DORY::<IP,LMC,RMC,IPC, D>::setup(rng, len).unwrap();

    // set matrix_a_trans : n*n length vector...
    /////// let matrix_a = MVP::<IP,LMC, RMC, IPC, D>::set_sparse_matrix(rng, len).unwrap();

    let mut matrix_a = Vec::new();
    for _ in 0..len{
        matrix_a.push(<LMC::Scalar>::rand(rng));
    }


    // for _ in 0..len*len{
    //     matrix_a.push(<LMC::Scalar>::rand(rng));
    // }
    
    // set z vector
    let mut z = Vec::new();
    for _ in 0..len{
        z.push(<LMC::Scalar>::rand(rng));
    }

    let generator_g1 = <IP::LeftMessage>::rand(rng);
    let _generator_g2 = <IP::RightMessage>::rand(rng);

    let mut start = Instant::now();
    let v_a = MVP::<IP, LMC, RMC, IPC, D>::collapse_matrix2(&matrix_a, &gamma2, len).unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t matrix collapsing time: {} ms", bench);

    let z_vec = MVP::<IP, LMC, RMC, IPC, D>::set_vector(&z, &generator_g1).unwrap();
    let a = MVP::<IP, LMC, RMC, IPC, D>::compute_az(&matrix_a, &z).unwrap();
    let a_vec = MVP::<IP, LMC, RMC, IPC, D>::set_vector(&a, &generator_g1).unwrap();


    let d1 = IP::inner_product(&z_vec, &gamma2).unwrap();
    let d2 = IP::inner_product(&gamma1, &v_a).unwrap();
    let c = IP::inner_product(&z_vec, &v_a).unwrap();
    
    let d1_ = IP::inner_product(&a_vec, &gamma2).unwrap();
    let d2_ = IP::inner_product(&gamma1, &gamma2).unwrap();
    let c_ = IP::inner_product(&a_vec, &gamma2).unwrap();

    start = Instant::now();
    let (x, bat_l, bat_r, delta) 
        = DORY::<IP, LMC, RMC, IPC, D>::batch_commit(
            &z_vec.clone(), &v_a.clone(),
            &a_vec.clone(), &gamma2.clone(), 
            // rng
        ).unwrap();
    bench = start.elapsed().as_millis();
    println!("\t batching time: {} ms", bench);

    let (bat_c, bat_d1, bat_d2)
        = DORY::<IP, LMC, RMC, IPC, D>::batch_verify (
            &c, &c_,
            &d1, &d1_,
            &d2, &d2_,
            &x,
            &delta
        ).unwrap();

    let dory_srs = DORY::<IP, LMC, RMC, IPC, D>::precompute((&(gamma1.clone()), &(gamma2.clone()))).unwrap();
    start = Instant::now();
    let mut proof =
        DORY::<IP, LMC, RMC, IPC, D>::prove((&(bat_l.clone()), &(bat_r.clone())),
        //  (&(gamma1.clone()), &(gamma2.clone())), 
         (&(gamma1.clone()), &(gamma2.clone())), 
         (&(bat_d1.clone()), &(bat_d2.clone()), &(bat_c.clone()))
        ).unwrap();
    bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);
    // let mut proof_ =
    //     DORY::<IP, LMC, RMC, IPC, D>::prove((&(a_vec.clone()), &(gamma2.clone())),
    //     //  (&(gamma1.clone()), &(gamma2.clone())), 
    //      (&(gamma1.clone()), &(gamma2.clone())), 
    //      (&(d1_.clone()), &(d2_.clone()), &(c_.clone()))
    //     ).unwrap();
    
    start = Instant::now();
    let eq_c = c == c_;
    let result = DORY::<IP, LMC, RMC, IPC, D>::verify(&mut dory_srs.clone(), (&(gamma1.clone()), &(gamma2.clone())),
         (&(bat_d1.clone()), &(bat_d2.clone()), &(bat_c.clone())), &mut proof)
        .unwrap();
    bench = start.elapsed().as_millis();
    println!("\t verification time: {} ms", bench);
    // let result_ = DORY::<IP, LMC, RMC, IPC, D>::verify(&mut dory_srs.clone(), (&(gamma1.clone()), &(gamma2.clone())),
    //      (&(d1_.clone()), &(d2_.clone()), &(c_.clone())), &mut proof_)
    //     .unwrap();
    
    println!("C == C' : {}", eq_c);
    println!("result : {}", result);
}


fn main() {
    const LEN: usize = 2048;
    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    let mut rng = StdRng::seed_from_u64(0u64);

    println!("Benchmarking MV-product with matrix size: {} * {}, vector length: {}", LEN, LEN, LEN);

    println!("1) Pairing inner product...");
    bench_mvp::<
        PairingInnerProduct<Bls12_381>,
        GC1,
        GC2,
        IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, LEN);

}
