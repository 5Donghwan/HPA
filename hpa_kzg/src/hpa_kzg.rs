extern crate ark_ff;
use self::ark_ff::{to_bytes, Field, One, PrimeField, Zero, UniformRand};
extern crate ark_ec;
use self::ark_ec::{msm::FixedBaseMSM,PairingEngine, ProjectiveCurve};
extern crate ark_poly;
use self::ark_poly::polynomial::{univariate::DensePolynomial, UVPolynomial};
extern crate ark_serialize;
use self::ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};
extern crate ark_std;
use self::ark_std::rand::Rng;
use self::ark_std::{end_timer, start_timer};
extern crate digest;
use self::digest::Digest;
use std::{convert::TryInto, marker::PhantomData, ops::MulAssign};
extern crate itertools;
use self::itertools::Itertools;
use crate::{mul_helper, Error, InnerProductArgumentError};
extern crate ark_dh_commitments;
use self::ark_dh_commitments::{
    // afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    // pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};
extern crate ark_inner_products;
use self::ark_inner_products::{InnerProduct, MultiexponentiationInnerProduct};
use self::ark_std::cfg_iter;

use std::{fmt};  //clone};

#[cfg(feature = "parallel")]
extern crate rayon;
use self::rayon::prelude::*;

struct List(Vec<u8>);
impl fmt::Display for List {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Extract the value using tuple indexing,
        // and create a reference to 'vec'.
        let vec = &self.0;
        write!(f, "[")?;
        // Iterate over 'v' in 'vec' while enumerating the iteration
        // count in 'count'.
        for (count, v) in vec.iter().enumerate() {
            // For every element except the first, add a comma.
            // Use the ? operator to return on errors.
            if count != 0 {
                write!(f, ", ")?;
            }
            write!(f, "{}", v)?;
        }
        // Close the opened bracket and return a fmt::Result value.
        write!(f, "]")
    }
}

pub struct HPA<IP, CM, P, D> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<CM>,
    _right_commitment: PhantomData<CM>,
    _inner_product_commitment: PhantomData<CM>,
    _pair: PhantomData<P>,
    _digest: PhantomData<D>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct HPAProof<IP, CM, P, D>
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = CM::Message,
        RightMessage = CM::Message,
        Output = CM::Message,
    >,
    CM: DoublyHomomorphicCommitment,
    // RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // RMC::Message: MulAssign<LMC::Scalar>,
    // IPC::Message: MulAssign<LMC::Scalar>,
    // RMC::Key: MulAssign<LMC::Scalar>,
    // IPC::Key: MulAssign<LMC::Scalar>,
    // RMC::Output: MulAssign<LMC::Scalar>,
    // IPC::Output: MulAssign<LMC::Scalar>,
{
    pub(crate) r_commitment_steps: Vec<(
        (IP::Output, IP::Output, IP::Output, IP::Output),
        (IP::Output, IP::Output, IP::Output, IP::Output),
        (CM::Output, CM::Output, CM::Output, CM::Output),
        (CM::Output, CM::Output, CM::Output, CM::Output),
    )>,
    pub(crate) r_d1_x: Vec<CM::Output>,
    pub(crate) r_d2_x: Vec<CM::Output>,

    pub(crate) e1: Vec<IP::LeftMessage>,
    pub(crate) e2: Vec<IP::RightMessage>,
    pub(crate) q1: CM::Message,
    pub(crate) q2: CM::Message,
    pub(crate) q3: CM::Message,
    pub(crate) q4: CM::Message,
    pub(crate) p1: CM::Output,
    pub(crate) p2: CM::Output,
    pub(crate) r: CM::Message,
    pub(crate) alpha_transcript: Vec<CM::Scalar>,
    final_ck:(CM::Key, CM::Key),
    final_ck_proof: (P::G1Projective, P::G1Projective),
    _hpa: PhantomData<HPA<IP, CM, P, D>>,
}

// impl<IP, CM, P, D> Clone for HPAProof<IP, CM, P, D>
// where
//     D: Digest,
//     P: PairingEngine,
//     IP: InnerProduct<
//         LeftMessage = CM::Message,
//         RightMessage = CM::Message,
//         Output = CM::Message,
//     >,
//     CM: DoublyHomomorphicCommitment,
// {
//     fn clone(&self) -> Self {
//         Self {
//             pub(crate) r_commitment_steps: self.r_commitment_steps.clone(),
//             pub(crate) r_d1_x: self.r_d1_x.clone(),
//             pub(crate) r_d2_x: self.r_d2_x.clone(),
        
//             pub(crate) e1: self.e1.clone(),
//             pub(crate) e2: self.e1.clone(),
//             pub(crate) q1: self.e1.clone(),
//             pub(crate) q2: self.e1.clone(),
//             pub(crate) q3: self.e1.clone(),
//             pub(crate) q4: self.e1.clone(),
//             pub(crate) p1: self.e1.clone(),
//             pub(crate) p2: self.e1.clone(),
//             pub(crate) r: self.e1.clone(),
//             final_ck: self. final_ck.clone(),
//             final_ck_proof: self.final_ck_proof.clone(),
//             _hpa: self._hpa.clone(),
//         }
//     }
// }

#[derive(Clone)]
pub struct SRS<P: PairingEngine> {
    pub g_alpha_powers: Vec<P::G1Projective>,
    pub g_beta_powers: Vec<P::G1Projective>,
    pub h_beta_powers : Vec<P::G2Projective>,
    pub h_beta: P::G2Projective,
    pub h_alpha: P::G2Projective,
}

#[derive(Clone)]
pub struct VerifierSRS<P: PairingEngine> {
    pub g: P::G1Projective,
    pub h: P::G2Projective,
    pub h_beta: P::G2Projective,
    pub h_alpha: P::G2Projective,
}

impl<P: PairingEngine> SRS<P> {
    pub fn get_commitment_keys(&self) -> (Vec<P::G1Projective>, Vec<P::G1Projective>) {
        let ck_1 = self.g_beta_powers.iter().step_by(2).cloned().collect();
        let ck_2 = self.g_alpha_powers.iter().step_by(2).cloned().collect();
        (ck_1, ck_2)
    }

    pub fn get_verifier_key(&self) -> VerifierSRS<P> {
        VerifierSRS { 
            g: self.g_alpha_powers[0].clone(),
            h: self.h_beta_powers[0].clone(),
            h_beta: self.h_beta.clone(),
            h_alpha: self.h_alpha.clone(),
        }
    }
}

// #[derive(Clone)]
// pub struct HPASRS<IP, CM, D>
// where
//     D: Digest,
//     IP: InnerProduct<
//         LeftMessage = LMC::Message,
//         RightMessage = RMC::Message,
//         Output = IPC::Message,
//     >,
//     CM: DoublyHomomorphicCommitment,
//     // RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     // IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     // RMC::Message: MulAssign<LMC::Scalar>,
//     // IPC::Message: MulAssign<LMC::Scalar>,
//     // RMC::Key: MulAssign<LMC::Scalar>,
//     // IPC::Key: MulAssign<LMC::Scalar>,
//     // RMC::Output: MulAssign<LMC::Scalar>,
//     // IPC::Output: MulAssign<LMC::Scalar>,
// {
//     pub(crate) delta1_l: Vec<IP::Output>,
//     pub(crate) delta1_r: Vec<IP::Output>,
//     pub(crate) delta2_l: Vec<IP::Output>,
//     pub(crate) delta2_r: Vec<IP::Output>,
//     pub(crate) kai: Vec<IP::Output>,
//     pub(crate) ht: IP::Output,

//     _hpa: PhantomData<HPA<IP, CM, D>>,
// }

#[derive(Clone)]
pub struct HPAAux<IP, CM, P, D>
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = CM::Message,
        RightMessage = CM::Message,
        Output = CM::Message,
    >,
    CM: DoublyHomomorphicCommitment,
    // RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // RMC::Message: MulAssign<LMC::Scalar>,
    // IPC::Message: MulAssign<LMC::Scalar>,
    // RMC::Key: MulAssign<LMC::Scalar>,
    // IPC::Key: MulAssign<LMC::Scalar>,
    // RMC::Output: MulAssign<LMC::Scalar>,
    // IPC::Output: MulAssign<LMC::Scalar>,
{
    //  pub(crate) r_transcript: Vec<(CM::Scalar, CM::Scalar, CM::Scalar, LMCMC::Scalar)>,
    // pub(crate) ck_base: (LMC::Key, RMC::Key),
    _hpa: PhantomData<HPA<IP, CM, P, D>>,
}

//TODO: Can extend HPA to support "identity commitments" in addition to "compact commitments", i.e. for SIPP

impl<IP, CM, P, D> HPA<IP, CM, P, D>
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = CM::Message,
        RightMessage = CM::Message,
        Output = CM::Message,
    >,
    CM: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G1Projective>,
    // LMC: DoublyHomomorphicCommitment,
    // RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // RMC::Message: MulAssign<LMC::Scalar>,
    // IPC::Message: MulAssign<LMC::Scalar>,
    // RMC::Key: MulAssign<LMC::Scalar>,
    // LMC::Key: MulAssign<LMC::Scalar>,
    // IPC::Key: MulAssign<LMC::Scalar>,
    // RMC::Output: MulAssign<LMC::Scalar>,
    // IPC::Output: MulAssign<LMC::Scalar>,
    CM::Message: MulAssign<CM::Message>,
    CM::Message: MulAssign<CM::Scalar>,
    CM::Key: MulAssign<CM::Scalar>,
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
    IP::Output: MulAssign<CM::Scalar>,
{

    pub fn set_values(
        l: &[<CM as DoublyHomomorphicCommitment>::Scalar],
        r: &[<CM as DoublyHomomorphicCommitment>::Scalar],
        // generator_g1: &IP::LeftMessage,
        // generator_g2: &IP::RightMessage,
    ) -> Result<
        (
            Vec<<CM as DoublyHomomorphicCommitment>::Scalar>,
            Vec<<CM as DoublyHomomorphicCommitment>::Scalar>,
            // Vec<IP::LeftMessage>,
            // Vec<IP::RightMessage>
        ), Error >
    {
        // let mut v1 = Vec::new();
        // let mut v2 = Vec::new();
        // let mut u1 = Vec::new();
        let mut u2 = Vec::new();

        let len = l.len();
        let one = <CM as DoublyHomomorphicCommitment>::Scalar::one();
        // let g2_one = mul_helper(&generator_g2.clone(), &one);
        for _i in 0..len {
            // v1.push(mul_helper(&generator_g1.clone(), &l[i]));
            // v2.push(mul_helper(&generator_g2.clone(), &r[i]));
            u2.push(one.clone());
        }
        let u1 = cfg_iter!(l)
                    .zip(r)
                    .map(|(b_1, b_2)| *b_1 * b_2.clone())
                    .collect::<Vec<CM::Scalar>>();
        
        Ok((u1, u2))
    }

    pub fn init_commit<R: Rng>(
        left_value: &Vec<IP::LeftMessage>,
        right_value: &Vec<IP::RightMessage>,
        gamma1: &Vec<CM::Key>,
        gamma2: &Vec<CM::Key>,
        // h1: &Vec<IP::LeftMessage>,
        // h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            CM::Output,
            CM::Output,
            IP::Output,
            CM::Output,
            CM::Scalar,
            Vec<CM::Scalar>,
            CM::Message,
            CM::Message,
            Vec<IP::LeftMessage>,
        ),
        Error,
    > {
        let l = left_value.clone();
        let r = right_value.clone();
        let gamma1 = gamma1.clone();
        let gamma2 = gamma2.clone();
        

        let r_c = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
        let r_x = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);

        let c = IP::inner_product(&l, &r)? + r_c.clone();
        let d1 = CM::commit(&gamma1, &l).unwrap();
        let d2 = CM::commit(&gamma2, &r)?;

        // Fiat-Schamir challenge
        // let gm = 'challenge: loop {
        //     let mut hash_input = Vec::new();
        //     //TODO: Should use CanonicalSerialize instead of ToBytes
        //     hash_input.extend_from_slice(&to_bytes![c, d1, d2]?);
        //     let gm: CM::Scalar =
        //         u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
        //             .into();
        //     break 'challenge gm;
        // };
        let gm = <CM as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let mut gm_vec = Vec::new();
        gm_vec.push(<CM as DoublyHomomorphicCommitment>::Scalar::one());
        for i in 1..l.len() {
            gm_vec.push(gm_vec[i - 1] * gm);
        }
        let mut w_vec = Vec::new();
        for i in 0..l.len() {
            w_vec.push( mul_helper(&l[i], &gm_vec[i]) );
        }
        
        let x = IP::inner_product(&w_vec, &r).unwrap() + r_x.clone();
        // let y = IP::inner_product(&k_vec, &r).unwrap() + mul_helper(&ht, &r_y);
        let d3 = CM::commit(&gamma1, &w_vec).unwrap();//IP::inner_product(&w_vec, &gamma2).unwrap();
        // let d4 = IP::inner_product(&k_vec, &gamma2).unwrap() + mul_helper(&ht, &r_d4);

        Ok((
            c, d1, d2, x, d3, gm, gm_vec, r_c, r_x, w_vec,
        ))
    }

    pub fn init_commit2<R: Rng>(
        left_value: &Vec<IP::LeftMessage>,
        right_value: &Vec<IP::RightMessage>,
        gamma1: &Vec<<CM as DoublyHomomorphicCommitment>::Key>,
        gamma2: &Vec<<CM as DoublyHomomorphicCommitment>::Key>,
        // h1: &Vec<IP::LeftMessage>,
        // h2: &Vec<IP::RightMessage>,
        // gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        gm_vec: &Vec<<CM as DoublyHomomorphicCommitment>::Scalar>,
        // r_c: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        // r_d1: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        // r_d2: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_x: &<CM as DoublyHomomorphicCommitment>::Message,
        // r_y: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        // r_d3: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        // r_d4: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            CM::Output,
            CM::Output,
            IP::Output,
            CM::Output,
            CM::Message,
            // IP::Output,
            // IP::Output,
            // LMC::Scalar,
            // Vec<LMC::Scalar>,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>,
        ),
        Error,
    > {
        let l = left_value.clone();
        let r = right_value.clone();
        let gamma1 = gamma1.clone();
        let gamma2 = gamma2.clone();


        let r_c = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);


        let c = IP::inner_product(&l, &r)? + r_c.clone();
        let d1 = CM::commit(&gamma1, &l).unwrap();//IP::inner_product(&l, &gamma2)?;
        let d2 = CM::commit(&gamma2, &r).unwrap();//IP::inner_product(&gamma1, &r)?;

        // Fiat-Schamir challenge
        // let gm = 'challenge: loop {
        //     let mut hash_input = Vec::new();
        //     //TODO: Should use CanonicalSerialize instead of ToBytes
        //     hash_input.extend_from_slice(&to_bytes![c, d1, d2]?);
        //     let gm: LMC::Scalar =
        //         u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
        //             .into();
        //     break 'challenge gm;
        // };
        // let mut gm_vec = Vec::new();
        // gm_vec.push(<LMC as DoublyHomomorphicCommitment>::Scalar::one());
        // for i in 1..l.len() {
        //     gm_vec.push(gm_vec[i - 1] * gm);
        // }
        let mut w_vec = Vec::new();
        for i in 0..l.len() {
            w_vec.push( mul_helper(&l[i],&gm_vec[i]));
        }

        let x = IP::inner_product(&w_vec, &r).unwrap() + r_x.clone();
        // let y = IP::inner_product(&k_vec, &r).unwrap() + mul_helper(&ht, &r_y);
        let d3 = CM::commit(&gamma1, &w_vec).unwrap();//IP::inner_product(&w_vec, &gamma2).unwrap();
        // let d4 = IP::inner_product(&k_vec, &gamma2).unwrap() + mul_helper(&ht, &r_d4);

        Ok((
            c, d1, d2, x, d3, r_c, w_vec,
        ))
    }


    pub fn batch_commit<R: Rng>(
        r_c: &CM::Message, r_c_: &CM::Message,
        r_x: &CM::Message,
        v1: &Vec<CM::Message>, u1: &Vec<CM::Message>,
        v2: &Vec<CM::Message>, u2: &Vec<CM::Message>,
        w_vec: &Vec<CM::Message>, w_vec_: &Vec<CM::Message>,
        // h1: &Vec<IP::LeftMessage>, h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            Vec<CM::Message>,
            Vec<CM::Message>,
            CM::Message,
            CM::Message,
            Vec<CM::Message>,
            CM::Scalar,
            // z_c, z_x, v1, v2, r_c, r_x, w_vec, delta
        ),
        Error,
    > {
        // let ht = IP::inner_product(&h1, &h2).unwrap();

        let r_zc = CM::Message::rand(rng);
        let r_zx = CM::Message::rand(rng);
        
        let z_c = IP::inner_product(&v1, &u2).unwrap() + IP::inner_product(&u1, &v2).unwrap() + r_zc.clone();
        let z_x = IP::inner_product(&w_vec, &u2).unwrap() + IP::inner_product(&w_vec_, &v2).unwrap() + r_zx.clone();

        let delta = 'challenge: loop {
            let mut hash_input = Vec::new();
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![z_c, z_x]?);
            let delta: CM::Scalar =
                u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
                    .into();
            break 'challenge delta;
        };
        let delta_sqr = delta * delta;

        let mut bat_v1 = Vec::new();
        let mut bat_v2 = Vec::new();
        let mut bat_w_vec = Vec::new();

        for i in 0..v1.len(){
            bat_v1.push(mul_helper(&v1[i], &delta) + u1[i].clone());
        }
        for i in 0..v1.len(){
            bat_v2.push(mul_helper(&v2[i], &delta) + u2[i].clone());
        }
        for i in 0..v1.len(){
            bat_w_vec.push(mul_helper(&w_vec[i], &delta) + w_vec_[i].clone());
        }

        let bat_r_c = mul_helper(r_c, &delta_sqr) + mul_helper(&r_zc, &delta) + r_c_.clone();
        // let bat_r_c = delta_sqr * r_c.clone() + delta * r_zc.clone() + r_c_;
        // let bat_r_x = delta_sqr * r_x.clone() + delta * r_zx + r_x;
        let bat_r_x = mul_helper(r_x, &delta_sqr) + mul_helper(&r_zx, &delta) + r_x.clone();


        Ok((
            z_c, z_x, bat_v1, bat_v2, bat_r_c, bat_r_x, bat_w_vec, delta
        ))
    }

    pub fn setup<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<SRS<P>, Error> {
        // println!("size : {}", size);
        let alpha  = <P::Fr>::rand(rng);
        let beta = <P::Fr>::rand(rng);
        let g = <P::G1Projective>::prime_subgroup_generator();
        let h = <P::G2Projective>::prime_subgroup_generator();
        Ok(
            SRS {
                g_alpha_powers: structured_generators_scalar_power(2 * size - 1, &g, &alpha),
                g_beta_powers: structured_generators_scalar_power(2 * size - 1, &g, &beta),
                h_beta_powers: structured_generators_scalar_power(2 * size - 1, &h, &beta),
                h_beta: h.mul(beta.into_repr()),
                h_alpha: h.mul(alpha.into_repr()),
            }
        )
    }
    // pub fn setup<R: Rng>(
    //     rng: &mut R,
    //     size: usize,
    // ) -> Result<(Vec<CM::Key>, Vec<CM::Key>), Error> {
    //     // println!("size : {}", size);
    //     let alpha  = <P::Fr>::rand(rng);
    //     let beta = <P::Fr>::rand(rng);
    //     let g = <P::G1Projective>::prime_subgroup_generator();
    //     let h = <P::G2Projective>::prime_subgroup_generator();
        
        
        
    //     //let gamma1 = CM::setup(rng, size)?;
    //     //let gamma2 = CM::setup(rng, size)?;
    //     // println!("gamma1_len : {}", gamma1.len());
    //     // println!("gamma2_len : {}", gamma2.len());

    //     Ok((gamma1, gamma2))
    // }

    // pub fn precompute(
    //     ck_message: (&[CM::Message], &[CM::Message]),
    //     // h1: &[LMC::Message],
    //     // h2: &[RMC::Message],
    // ) -> Result<HPASRS<IP, CM, D>, Error> {
    //     // loop : until ck.len() >= 1
    //     let (mut gamma1, mut gamma2) = ck_message.clone();
    //     let h1 = h1.clone();
    //     let h2 = h2.clone();
    //     // let mut i = ck_message.0.len();
    //     let mut delta1_l = Vec::new();
    //     let mut delta1_r = Vec::new();
    //     let mut delta2_r = Vec::new();
    //     let mut kai = Vec::new();
    //     let mut split = ck_message.0.len() / 2;
    //     while split >= 1 {
    //         kai.push(IP::inner_product(gamma1, gamma2)?);
    //         // Generate gamma1, gamma2, gamma1_prime, gamma2_prime
    //         let gamma1_prime = &gamma1[..split];
    //         let gamma2_prime = &gamma2[..split];
    //         // Generate gamma1_R, gamma2_R (replace gamma1_L to gamma1_prime)
    //         let gamma1_r = &gamma1[split..];
    //         let gamma2_r = &gamma2[split..];

    //         // Compute delta1_L, delta1_R, delta2_L, delta2_R, kai
    //         delta1_l.push(IP::inner_product(gamma1_prime, gamma2_prime)?);
    //         delta1_r.push(IP::inner_product(gamma1_r, gamma2_prime)?);
    //         delta2_r.push(IP::inner_product(gamma1_prime, gamma2_r)?);

    //         split = split / 2;
    //         gamma1 = gamma1_prime;
    //         gamma2 = gamma2_prime;
    //         // i = i/2;
    //     }
    //     let mut delta2_l = delta1_l.clone();
    //     delta1_l.reverse();
    //     delta1_r.reverse();
    //     delta2_l.reverse();
    //     delta2_r.reverse();
    //     kai.reverse();

    //     let ht = IP::inner_product(&h1, &h2)?;

    //     Ok(HPASRS {
    //         delta1_l: delta1_l,
    //         delta1_r: delta1_r,
    //         delta2_l: delta2_l,
    //         delta2_r: delta2_r,
    //         kai: kai,
    //         ht: ht,
    //         _hpa: PhantomData,
    //     })
    // }

    pub fn prove<R: Rng>(
        srs: &SRS<P>,
        values: (
            &[IP::LeftMessage],
            &[IP::RightMessage],
            &[IP::LeftMessage],
            // &[IP::LeftMessage],
        ),
        // srs: &HPASRS<IP, CM, D>,
        ck_message: (&[CM::Key], &[CM::Key]),
        // com: (&IP::Output, &IP::Output, &IP::Output),
        witness: (
            &<CM as DoublyHomomorphicCommitment>::Message,
            &<CM as DoublyHomomorphicCommitment>::Message,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<CM as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<HPAProof<IP, CM, P, D>, Error> {
        // if IP::inner_product(values.0, values.1)? != com.2.clone() {
        //     return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        // }
        // if values.0.len().count_ones() != 1 {
        //     // Power of 2 length
        //     return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
        //         values.0.len(),
        //         values.1.len(),
        //     )));
        // }
        // if !(LMC::verify(ck.1, values.0, com.0)?
        //     && RMC::verify(ck.0, values.1, com.1)?
        //     )
        // {
        //     return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        // }

        // TODO : compare ck and ck_message

        let (proof, _) = Self::prove_with_aux(
            srs,
            (values.0, values.1, values.2),// values.3),
            // srs,
            (ck_message.0, ck_message.1),
            witness,
            gm,
            rng,
        )?;
        Ok(proof)
    }

    pub fn verify(
        v_srs: &VerifierSRS<P>,
        ck_message: (Vec<CM::Key>, Vec<CM::Key>),
        com: (&IP::Output, &IP::Output, &CM::Output, &CM::Output, &CM::Output), // com ( c, x, d1, d2, d3 )
        proof: &mut HPAProof<IP, CM, P, D>,
        gm: &<CM as DoublyHomomorphicCommitment>::Scalar,
        // rng: &mut R
    ) -> Result<bool, Error> {
        if ck_message.0.len().count_ones() != 1 || ck_message.0.len() != ck_message.1.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                ck_message.0.len(),
                ck_message.1.len(),
            )));
        }
        // Calculate transcript
        
        let (mut transcript, ch_c) = Self::_compute_recursive_challenges(proof, gm)?;

        //let ( gamma1, gamma2) = ck_message.clone();
        let ( gamma1, gamma2) = proof.final_ck.clone();

        let round = transcript.len();
        // let mut c_prime : &IP::Output;
        // let mut d1_prime : &IP::Output;
        // let mut d2_prime : &IP::Output;

        // println!("round : {}", round);
        // println!("kai len : {}", srs.kai.len());
        // println!("r_commitment_steps len : {}", proof.r_commitment_steps.len());

        // let c = com.2.clone();
        // let d1 = com.0.clone();
        // let d2 = com.1.clone();

        let mut c_prime = com.0.clone();
        let mut x_prime = com.1.clone();
        let mut d1_prime = com.2.clone();
        let mut d2_prime = com.3.clone();
        let mut d3_prime = com.4.clone();
        let mut result = false;

        if round > 0 {
            for i in 0..round {
                //println!("check");
                // Verifier's work in reduce
                //let split = gamma1.len() / 2;

                let last_commitment = proof.r_commitment_steps.pop().unwrap();
                let last_transcript = transcript.pop().unwrap();
                let alpha_transcript = proof.alpha_transcript.clone();
                let last_d1_x = proof.r_d1_x.pop().unwrap();
                let last_d2_x = proof.r_d2_x.pop().unwrap();

                let c_l = last_commitment.0.0.clone();
                let c_r = last_commitment.0.1.clone();
                let c_x = last_commitment.0.2.clone();
                let x_l = last_commitment.1.0.clone();
                let x_r = last_commitment.1.1.clone();
                let x_plus = last_commitment.1.2.clone();
                let x_minus = last_commitment.1.3.clone();
                let d1_x = last_d1_x.clone();
                let d1_l = last_commitment.2.0.clone();
                let d1_r = last_commitment.2.1.clone();
                let d2_x = last_d2_x.clone();
                let d2_l = last_commitment.2.2.clone();
                let d2_r = last_commitment.2.3.clone();
                // let d1_plus = last_commitment.2.0.clone();
                // let d1_minus = last_commitment.2.1.clone();
                // let d2_plus = last_commitment.2.2.clone();
                // let d2_minus = last_commitment.2.3.clone();
                let d3_l = last_commitment.3.0.clone();
                let d3_r = last_commitment.3.1.clone();
                let d3_plus = last_commitment.3.2.clone();
                let d3_minus = last_commitment.3.3.clone();

                // let gamma1_l = &gamma1[..split];
                // let gamma1_r = &gamma1[split..];
                // let gamma2_l = &gamma2[..split];
                // let gamma2_r = &gamma2[split..];
                

                let alpha = last_transcript.0;
                // println!("ver al : {}", alpha);
                // let alpha_inv = last_transcript.1;
                let gm_inv = last_transcript.2;


                // gamma1 = cfg_iter!(gamma1_l)
                //     .map(|b| mul_helper(b, &alpha))
                //     .zip(gamma1_r)
                //     .map(|(b_1, b_2)| b_1 + b_2.clone())
                //     .collect::<Vec<CM::Key>>();
                // gamma2 = cfg_iter!(gamma2_l)
                //     .map(|b| mul_helper(b, &alpha))
                //     .zip(gamma2_r)
                //     .map(|(b_1, b_2)| b_1 + b_2.clone())
                //     .collect::<Vec<CM::Key>>();

                assert!(c_prime == c_l.clone() + c_r.clone());
                assert!(x_prime == x_l.clone() + x_r.clone());
                assert!(d1_prime == d1_l.clone() + d1_r.clone());
                assert!(d2_prime == d2_l.clone() + d2_r.clone());
                assert!(d3_prime == d3_l.clone() + d3_r.clone());

                let alpha_sqr = alpha * alpha;
                let alpha_gm_inv = alpha * gm_inv;
                // let alpha_inv_gm_inv = alpha_inv * gm_inv;

                c_prime = mul_helper(&c_l, &alpha_sqr) + c_r + mul_helper(&c_x, &alpha);
                x_prime = mul_helper(&x_l, &alpha_sqr) + mul_helper(&x_r, &gm_inv) + mul_helper(&x_plus, &alpha) + mul_helper(&x_minus, &alpha_gm_inv);
                d1_prime = mul_helper(&d1_l, &alpha_sqr) + d1_r + mul_helper(&d1_x, &alpha);
                // d1_prime = d1_prime + mul_helper(&d1_plus, &alpha) + mul_helper(&d1_minus, &alpha_inv);
                d2_prime = mul_helper(&d2_l, &alpha_sqr) + d2_r + mul_helper(&d2_x, &alpha);
                // d2_prime = d2_prime + mul_helper(&d2_plus, &alpha) + mul_helper(&d2_minus, &alpha_inv);
                d3_prime = mul_helper(&d3_l, &alpha_sqr) + mul_helper(&d3_r, &gm_inv) + mul_helper(&d3_plus, &alpha) + mul_helper(&d3_minus, &alpha_gm_inv);

              

                // Scalar product
                if i == round - 1 {
                      // Verify commitment keys wellformed
                let (ck_a_final, ck_b_final) = &proof.final_ck;
                let (ck_a_proof, ck_b_proof) = &proof.final_ck_proof;
                //alpha_transcript.reverse();
                // KZG challenge point
                let mut counter_nonce: usize = 0;
                let c = loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input.extend_from_slice(&to_bytes![
                        alpha_transcript.first().unwrap(),
                        ck_a_final,
                        ck_b_final
                    ]?);
                    if let Some(c) = CM::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                        break c;
                    };
                    counter_nonce += 1;
                };

                let _ck_a_valid = verify_commitment_key_g2_kzg_opening(
                    v_srs,
                    &ck_a_final,
                    &ck_a_proof,
                    &alpha_transcript,
                    &<P::Fr>::one(),
                    &c,
                )?;
                let _ck_b_valid = verify_commitment_key_g1_kzg_opening(
                    v_srs,
                    &ck_b_final,
                    &ck_b_proof,
                    &alpha_transcript,
                    &<P::Fr>::one(),
                    &c,
                )?;
                    let e1 = proof.e1.clone();
                    let e2 = proof.e2.clone();


                    let mut ch_c_vec = Vec::new();
                    ch_c_vec.push(ch_c.clone()); // ch_c_vec = {c^1, c^2, c^3, ..., c^7}
                    for i in 1..7{
                        ch_c_vec.push(ch_c_vec[i-1] * ch_c);
                    }

                    let mut result1 = false;
                    let mut result2 = false;
                    let mut result3 = false;

                    let temp_left = IP::inner_product(&e1, &e2).unwrap();
                    let mut temp_right = mul_helper(&c_prime, &(ch_c_vec[4].clone())) + mul_helper(&x_prime, &(ch_c_vec[5] + ch_c_vec[6]).clone()) 
                    + proof.q1.clone() + mul_helper(&proof.q2, &(ch_c_vec[3].clone())) + mul_helper(&proof.q3, &((ch_c_vec[0] + ch_c_vec[1]).clone())) + mul_helper(&proof.q4, &(ch_c_vec[2]));
                    
                    let zero = <CM as DoublyHomomorphicCommitment>::Scalar::zero();
                    let minus_one = zero - <CM as DoublyHomomorphicCommitment>::Scalar::one();

                    temp_right = temp_right + mul_helper(&proof.r, &minus_one);
                    if temp_left == temp_right {
                        result1 = true;
                    }
                    // println!("g : {}, e: {}", gamma1.len(), e1.len());

                    let temp_left = CM::commit(&[gamma1], &e1).unwrap();
                    let temp_right = mul_helper(&d1_prime, &ch_c_vec[0]) + mul_helper(&d3_prime, &(ch_c_vec[1] + ch_c_vec[2])) + proof.p1.clone();
                    if temp_left == temp_right {
                        result2 = true;
                    }

                    let temp_left = CM::commit(&[gamma2], &e2).unwrap();
                    let temp_right = mul_helper(&d2_prime, &ch_c_vec[3]) + proof.p2.clone();
                    if temp_left == temp_right {
                        result3 = true;
                    }
                    
                    result = result1 && result2 && result3;
                    //println!("ck_a_valid : {} , _ck_b_valid : {}", _ck_a_valid, _ck_b_valid);
                    result = result && _ck_a_valid && _ck_b_valid;
                    //println!("result1 : {}, result2 : {}, result3 : {}", result1, result2, result3);
                }
                
            }
            Ok(result)
        } else {
            Ok(result)
        }
    }

    pub fn prove_with_aux<R: Rng>(
        srs: &SRS<P>,
        values: (
            &[IP::LeftMessage],
            &[IP::RightMessage],
            &[IP::LeftMessage],
            // &[IP::LeftMessage],
        ),
        // srs: &HPASRS<IP, CM, D>,
        ck_message: (&[CM::Key], &[CM::Key]),
        witness: (
            &<CM as DoublyHomomorphicCommitment>::Message,
            &<CM as DoublyHomomorphicCommitment>::Message,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<CM as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<(HPAProof<IP, CM, P, D>, HPAAux<IP, CM, P, D>), Error> {
        let (v1, v2, w_vec) = values;
        // let (gamma1, gamma2) = ck;
        let (gamma1_message, gamma2_message) = ck_message;
        Self::_prove(
            srs,
            &(v1.to_vec(), v2.to_vec(), w_vec.to_vec()),// k_vec.to_vec()),
            // srs,
            (gamma1_message.to_vec(), gamma2_message.to_vec()),
            witness,
            gm,
            rng,
        )
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove<R: Rng>(
        srs: &SRS<P>,
        values: &(
            Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>,
            Vec<IP::LeftMessage>,
        ),
        ck_message: (Vec<CM::Key>, Vec<CM::Key>),
        witness: (
            &<CM as DoublyHomomorphicCommitment>::Message,
            &<CM as DoublyHomomorphicCommitment>::Message,
        ),
        gm: &<CM as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<(HPAProof<IP, CM, P, D>, HPAAux<IP, CM, P, D>), Error> {
        let (mut v1, mut v2, mut w_vec) = values.clone();
        let (mut gamma1_message, mut gamma2_message) = ck_message.clone();
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        let mut alpha_transcript = Vec::new();
        let mut r_d1_x = Vec::new();
        let mut r_d2_x = Vec::new();
        assert!(v1.len().is_power_of_two());

        let mut r_c = witness.0.clone();
        let mut r_x = witness.1.clone();

        let (_m_base, _ck_base) = 'recurse: loop {
            let recurse = start_timer!(|| format!("Recurse round size {}", m_a.len()));
            if v1.len() == 1 {
                // base case
                break 'recurse (
                    (v1[0].clone(), v2[0].clone()),
                    (gamma1_message[0].clone(), gamma2_message[0].clone()),
                );
            } else {
                // recursive step
                // Recurse with problem of half size
                let split = v1.len() / 2;

                let r_cl = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
                let r_cx = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
                let r_xl = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
                let r_x_plus = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
                let r_x_minus = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);

                let zero = <CM as DoublyHomomorphicCommitment>::Scalar::zero();
                let minus_one = zero - <CM as DoublyHomomorphicCommitment>::Scalar::one();

                let r_cr = r_c + mul_helper(&r_cl, &minus_one);
                let r_xr = r_x + mul_helper(&r_xl, &minus_one);
                // let r_d1r = r_d1 - r_d1l;
                // let r_d3r = r_d3 - r_d3l;

                let v1_l = &v1[..split];
                let v1_r = &v1[split..];
                let gamma1_l = &gamma1_message[..split];
                let gamma1_r = &gamma1_message[split..];

                let v2_l = &v2[..split];
                let v2_r = &v2[split..];
                let gamma2_l = &gamma2_message[..split];
                let gamma2_r = &gamma2_message[split..];

                let w_vec_l = &w_vec[..split];
                let w_vec_r = &w_vec[split..];


                // println!("{} {}",gamma1_r.len(), v1_l.len());
                let cl = start_timer!(|| "Compute D");

                let c_l = IP::inner_product(&v1_l, &v2_l).unwrap() + r_cl.clone();
                let c_r = IP::inner_product(&v1_r, &v2_r).unwrap() + r_cr.clone();
                let c_x = IP::inner_product(&v1_l, &v2_r).unwrap() + IP::inner_product(&v1_r, &v2_l).unwrap() + r_cx.clone();
                let x_l = IP::inner_product(&w_vec_l, &v2_l).unwrap() + r_xl.clone();
                let x_r = IP::inner_product(&w_vec_r, &v2_r).unwrap() + r_xr.clone();
                let x_plus = IP::inner_product(&w_vec_l, &v2_r).unwrap() + r_x_plus.clone();
                let x_minus = IP::inner_product(&w_vec_r, &v2_l).unwrap() + r_x_minus.clone();
                let d1_l = CM::commit(&gamma1_l, &v1_l).unwrap();
                let d1_r = CM::commit(&gamma1_r, &v1_r).unwrap();
                let d1_x = CM::commit(&gamma1_l, &v1_r).unwrap() + CM::commit(&gamma1_r, &v1_l).unwrap();
                // let d1_plus = CM::commit(&gamma1_r, &v1_l).unwrap();
                // let d1_minus = CM::commit(&gamma1_l, &v1_r).unwrap();
                let d2_l = CM::commit(&gamma2_l, &v2_l).unwrap();
                let d2_r = CM::commit(&gamma2_r, &v2_r).unwrap();
                let d2_x = CM::commit(&gamma2_l, &v2_r).unwrap() + CM::commit(&gamma2_r, &v2_l).unwrap();
                // let d2_plus = CM::commit(&gamma2_r, &v2_l).unwrap();
                // let d2_minus = CM::commit(&gamma2_l, &v2_r).unwrap();
                let d3_l = CM::commit(&gamma1_l, &w_vec_l).unwrap();
                let d3_r = CM::commit(&gamma1_r, &w_vec_r).unwrap();
                let d3_plus = CM::commit(&gamma1_r, &w_vec_l).unwrap();
                let d3_minus = CM::commit(&gamma1_l, &w_vec_r).unwrap();


                // Fiat-Shamir challenge
                let mut counter_nonce: usize = 0;
                //  let default_transcript = (Default::default(),Default::default(),Default::default(),Default::default());
                //  let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (alpha, alpha_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input.extend_from_slice(&to_bytes![
                        c_l, c_r, c_x, 
                        x_l, x_r, x_plus, x_minus,
                        d1_l, d1_r, d2_l, d2_r,
                        d3_l, d3_r, d3_plus, d3_minus
                    ]?);
                    let alpha: CM::Scalar = u128::from_be_bytes(
                        D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                    )
                    .into();

                    //  let list_beta_input = List(hash_input.clone());
                    //  println!("prove - hash beta - {}", list_beta_input);

                    if let Some(alpha_inv) = alpha.inverse() {
                        // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                        // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                        break 'challenge (alpha_inv, alpha);
                    }
                    counter_nonce += 1;
                };
                // println!("pro al : {}", alpha);

                end_timer!(cl);

                let rescale_v1 = start_timer!(|| "Rescale V1");
                v1 = cfg_iter!(v1_l)
                    .map(|a| mul_helper(a, &alpha))
                    .zip(v1_r)
                    .map(|(a_1, a_2)| a_1 + a_2.clone())
                    .collect::<Vec<CM::Message>>();
                end_timer!(rescale_v1);

                let rescale_v2 = start_timer!(|| "Rescale V2");
                v2 = cfg_iter!(v2_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(v2_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<CM::Message>>();
                end_timer!(rescale_v2);

                let rescale_w_vec = start_timer!(|| "Rescale W");
                let mut gm_inv = gm.inverse().unwrap();
                let w_len = w_vec.len();
                let exp_m = f32::log2(w_len as f32) as usize;
                for _ in 0..exp_m - 1 {
                    gm_inv = gm_inv * gm_inv;
                }

                w_vec = cfg_iter!(w_vec_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(w_vec_r)
                    .map(|(b_1, b_2)| b_1 + mul_helper(b_2, &gm_inv))
                    .collect::<Vec<CM::Message>>();
                end_timer!(rescale_w_vec);

                let alpha_sqr = alpha * alpha;
                let gm_alpha = alpha * gm_inv;
                r_c = mul_helper(&r_cl, &alpha_sqr)
                    + r_cr
                    + mul_helper(&r_cx, &alpha);
                r_x = mul_helper(&r_xl, &alpha_sqr)
                    + mul_helper(&r_xr, &gm_inv)
                    + mul_helper(&r_x_plus, &alpha)
                    + mul_helper(&r_x_minus, &gm_alpha);


                gamma1_message = cfg_iter!(gamma1_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(gamma1_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<CM::Key>>();
                gamma2_message = cfg_iter!(gamma2_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(gamma2_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<CM::Key>>();
                

                let com1 = (c_l, c_r, c_x.clone(), c_x.clone());
                let com2 = (x_l, x_r, x_plus, x_minus);
                let com3 = (d1_l, d1_r, d2_l, d2_r);
                let com4 = (d3_l, d3_r, d3_plus, d3_minus);

                r_commitment_steps.push((com1, com2, com3, com4));
                r_transcript.push((alpha, alpha_inv, gm_inv));
                alpha_transcript.push(alpha);
                r_d1_x.push(d1_x);
                r_d2_x.push(d2_x);

                end_timer!(recurse);
            }
        };

        
        let v1_val = v1.pop().unwrap();
        let v2_val = v2.pop().unwrap();
        let w_val = w_vec.pop().unwrap();

        // println!("v1 ?= w : {}", v1_val == w_val);

        let r_d1 = <IP::LeftMessage>::rand(rng);
        let r_d2 = <IP::RightMessage>::rand(rng);
        let r_q1 = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
        let r_q2 = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
        let r_q3 = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);
        let r_q4 = <CM as DoublyHomomorphicCommitment>::Message::rand(rng);

         let mut d1 = Vec::new();
         let mut d2 = Vec::new();
         d1.push(r_d1.clone());
         d2.push(r_d2.clone());

        let q1 = mul_helper(&r_d1, &r_d2) + r_q1.clone();
        let q2 = mul_helper(&r_d1, &v2_val) + r_q2.clone();
        let q3 = mul_helper(&v1_val, &r_d2) + r_q3.clone();
        let q4 = mul_helper(&w_val, &r_d2) + r_q4.clone();
        let p1 = CM::commit(&gamma1_message, &d1).unwrap();
        let p2 = CM::commit(&gamma2_message, &d2).unwrap();


        let ch_c = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&to_bytes![q1, q2, q3, q4, p1, p2]?);
            let ch_c: CM::Scalar =
                u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
                    .into();

            break 'challenge ch_c;
        };

        let ch_c_2 = ch_c * ch_c;
        let ch_c_3 = ch_c_2 * ch_c;
        let ch_c_4 = ch_c_3 * ch_c;
        let ch_c_5 = ch_c_4 * ch_c;
        let ch_c_6 = ch_c_5 * ch_c;
        let ch_c_7 = ch_c_6 * ch_c;


        let e1 = r_d1 + mul_helper(&v1_val, &ch_c) + mul_helper(&v1_val, &ch_c_2) + mul_helper(&w_val, &ch_c_3);
        let e2 = r_d2 + mul_helper(&v2_val, &ch_c_4);
        let r = r_q1 + mul_helper(&r_q2, &ch_c_4) + mul_helper(&r_q3, &(ch_c + ch_c_2)) + mul_helper(&r_q4, &ch_c_3)
            + mul_helper(&r_c, &ch_c_5) + mul_helper(&r_x, &(ch_c_6 + ch_c_7));

        let mut e1_vec = Vec::new();
        let mut e2_vec = Vec::new();
        e1_vec.push(e1);
        e2_vec.push(e2);
        
        // r_transcript.reverse();
        r_commitment_steps.reverse();
        alpha_transcript.reverse();
        r_d1_x.reverse();
        r_d2_x.reverse();

        // Prove final commitment keys are wellformed
        let (ck_a_final, ck_b_final) = _ck_base.clone();
        let transcript = alpha_transcript.clone();
        //let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();
        //let r_inverse = r_shift.inverse().unwrap();
        //println!("ck_a_final: {} ck_b_final: {}", ck_a_final, ck_b_final);
        // KZG challenge point
        let mut counter_nonce: usize = 0;
        let c = loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![
                transcript.first().unwrap(),
                ck_a_final,
                ck_b_final
            ]?);
            if let Some(c) = CM::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                break c;
            };
            counter_nonce += 1;
        };
        // Complete KZG proofs
        let ck_a_kzg_opening = prove_commitment_key_kzg_opening(
            &srs.g_beta_powers,
            &transcript,
            &<P::Fr>::one(),
            &c,
        )?;
        let ck_b_kzg_opening = prove_commitment_key_kzg_opening(
            &srs.g_alpha_powers,
            &transcript,
            &<P::Fr>::one(),
            &c,
        )?;

        Ok((
            HPAProof {
                r_commitment_steps,
                r_d1_x,
                r_d2_x,
                e1: e1_vec, 
                e2: e2_vec,
                q1, q2, q3, q4,
                p1, p2, 
                r,
                alpha_transcript,
                final_ck: (ck_a_final, ck_b_final),
                final_ck_proof: (ck_a_kzg_opening, ck_b_kzg_opening),
                _hpa: PhantomData,
            },
            HPAAux { _hpa: PhantomData },
        ))
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        proof: &HPAProof<IP, CM, P, D>,
        gm: &<CM as DoublyHomomorphicCommitment>::Scalar
    ) -> Result<
        (
            Vec<(CM::Scalar, CM::Scalar, CM::Scalar)>,
            CM::Scalar,
        ),
        Error,
    > {
        Self::_compute_recursive_challenges(proof, gm)
    }

    fn _compute_recursive_challenges(
        proof: &HPAProof<IP, CM, P, D>,
        gm: &<CM as DoublyHomomorphicCommitment>::Scalar
    ) -> Result<
        (
            Vec<(CM::Scalar, CM::Scalar,  CM::Scalar)>,
            CM::Scalar,
        ),
        Error,
    > {

        let mut r_transcript = Vec::new();

        for (com_1, com_2, com_3, com_4) in proof.r_commitment_steps.iter().rev() {
            // First Fiat-Shamir challenge
            let mut counter_nonce: usize = 0;
            let (alpha, alpha_inv) = 'challenge: loop {
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                hash_input.extend_from_slice(&to_bytes![
                    com_1.0, com_1.1, com_1.2,
                    com_2.0, com_2.1, com_2.2, com_2.3,
                    com_3.0, com_3.1, com_3.2, com_3.3,
                    com_4.0, com_4.1, com_4.2, com_4.3
                ]?);
                let alpha: CM::Scalar = u128::from_be_bytes(
                    D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                )
                .into();

                // let list_beta_input = List(hash_input.clone());
                // println!("prove - hash beta - {}", list_beta_input);

                if let Some(alpha_inv) = alpha.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (alpha_inv, alpha);
                }
                counter_nonce += 1;
            };


            let mut gm_inv = gm.inverse().unwrap();
            let len = proof.r_commitment_steps.len() - r_transcript.len() - 1;
            // println!("len : {}", len);
            for _ in 0..len {
                gm_inv = gm_inv * gm_inv;
            }

            r_transcript.push((alpha, alpha_inv, gm_inv));
        }
        r_transcript.reverse();
        //println!("r_transcript len : {}", r_transcript.len());

        let ch_c = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&to_bytes![proof.q1, proof.q2, proof.q3, proof.q4, proof.p1, proof.p2]?);
            let ch_c: CM::Scalar =
                u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
                    .into();

            break 'challenge ch_c;
        };

        Ok((r_transcript, ch_c))
    }

    pub fn batch_verify(
        c: &IP::Output,
        c_: &IP::Output,
        x: &IP::Output,
        x_: &IP::Output,
        d1: &CM::Output,
        d1_: &CM::Output,
        d2: &CM::Output,
        d2_: &CM::Output,
        d3: &CM::Output,
        d3_: &CM::Output,
        z_c: &IP::Output,
        z_x: &IP::Output,
        delta: &CM::Scalar,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            CM::Output,
            CM::Output,
            CM::Output,
            // c, x, d1, d2, d3
        ),
        Error,
    > {

        let delta = delta.clone();
        let delta_sqr = delta * delta;
        let bat_c = mul_helper(c, &delta_sqr) + mul_helper(z_c, &delta) + c_.clone();
        let bat_x = mul_helper(x, &delta_sqr) + mul_helper(z_x, &delta) + x_.clone();
        let bat_d1 = mul_helper(d1, &delta) + d1_.clone();
        let bat_d2 = mul_helper(d2, &delta) + d2_.clone();
        let bat_d3 = mul_helper(d3, &delta) + d3_.clone();

        Ok((bat_c, bat_x, bat_d1, bat_d2, bat_d3))
    }

}

// pub(crate) fn _compute_final_commitment_keys(
//     ck: (&[LMC::Key], &[RMC::Key]),
//     transcript: &Vec<LMC::Scalar>,
// ) -> Result<(LMC::Key, RMC::Key), Error> {
//     // Calculate base commitment keys
//     let (ck_a, ck_b) = ck;
//     assert!(ck_a.len().is_power_of_two());

//     let mut ck_a_agg_challenge_exponents = vec![LMC::Scalar::one()];
//     let mut ck_b_agg_challenge_exponents = vec![LMC::Scalar::one()];
//     for (i, c) in transcript.iter().enumerate() {
//         let c_inv = c.inverse().unwrap();
//         for j in 0..(2_usize).pow(i as u32) {
//             ck_a_agg_challenge_exponents.push(ck_a_agg_challenge_exponents[j] * &c_inv);
//             ck_b_agg_challenge_exponents.push(ck_b_agg_challenge_exponents[j] * c);
//         }
//     }
//     assert_eq!(ck_a_agg_challenge_exponents.len(), ck_a.len());
//     //TODO: Optimization: Use VariableMSM multiexponentiation
//     let ck_a_base_init = mul_helper(&ck_a[0], &ck_a_agg_challenge_exponents[0]);
//     let ck_a_base = ck_a[1..]
//         .iter()
//         .zip(&ck_a_agg_challenge_exponents[1..])
//         .map(|(g, x)| mul_helper(g, &x))
//         .fold(ck_a_base_init, |sum, x| sum + x);
//     //.reduce(|| ck_a_base_init.clone(), |sum, x| sum + x);
//     let ck_b_base_init = mul_helper(&ck_b[0], &ck_b_agg_challenge_exponents[0]);
//     let ck_b_base = ck_b[1..]
//         .iter()
//         .zip(&ck_b_agg_challenge_exponents[1..])
//         .map(|(g, x)| mul_helper(g, &x))
//         .fold(ck_b_base_init, |sum, x| sum + x);
//     //.reduce(|| ck_b_base_init.clone(), |sum, x| sum + x);
//     Ok((ck_a_base, ck_b_base))
// }

// pub(crate) fn _verify_base_commitment(
//     base_ck: (&LMC::Key, &RMC::Key, &Vec<IPC::Key>),
//     base_com: (LMC::Output, RMC::Output, IPC::Output),
//     proof: &HPAProof<IP, LMC, RMC, IPC, D>,
// ) -> Result<bool, Error> {
//     let (com_a, com_b, com_t) = base_com;
//     let (ck_a_base, ck_b_base, ck_t) = base_ck;
//     let a_base = vec![proof.r_base.0.clone()];
//     let b_base = vec![proof.r_base.1.clone()];
//     let t_base = vec![IP::inner_product(&a_base, &b_base)?];

//     Ok(LMC::verify(&vec![ck_a_base.clone()], &a_base, &com_a)?
//         && RMC::verify(&vec![ck_b_base.clone()], &b_base, &com_b)?
//         && IPC::verify(&ck_t, &t_base, &com_t)?)
// }

impl<IP, CM, P, D> Clone for HPAProof<IP, CM, P, D>
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = CM::Message,
        RightMessage = CM::Message,
        Output = CM::Message,
    >,
    CM: DoublyHomomorphicCommitment,
    // RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    // RMC::Message: MulAssign<LMC::Scalar>,
    // IPC::Message: MulAssign<LMC::Scalar>,
    // RMC::Key: MulAssign<LMC::Scalar>,
    // IPC::Key: MulAssign<LMC::Scalar>,
    // RMC::Output: MulAssign<LMC::Scalar>,
    // IPC::Output: MulAssign<LMC::Scalar>,
{
    fn clone(&self) -> Self {
        HPAProof {
            r_commitment_steps: self.r_commitment_steps.clone(),
            r_d1_x: self.r_d1_x.clone(),
            r_d2_x: self.r_d2_x.clone(),
            // c_x: self.c_x.clone(),
            // x_plus: self.x_plus.clone(),
            // x_minus: self.x_minus.clone(),
            // y_plus: self.y_plus.clone(),
            // y_minus: self.y_minus.clone(),
            e1: self.e1.clone(),
            e2: self.e2.clone(),
            q1: self.q1.clone(),
            q2: self.q2.clone(),
            q3: self.q3.clone(),
            q4: self.q4.clone(),
            p1: self.p1.clone(),
            p2: self.p2.clone(),
            r: self.r.clone(),
            alpha_transcript: self.alpha_transcript.clone(),
            final_ck: self.final_ck.clone(),
            final_ck_proof: self.final_ck_proof.clone(),
            // r1: self.r1.clone(),
            // r2: self.r2.clone(),
            // r3: self.r3.clone(),
            // r_base: self.r_base.clone(),
            _hpa: PhantomData,
        }
    }
}

pub fn prove_commitment_key_kzg_opening<G: ProjectiveCurve>(
    srs_powers: &Vec<G>,
    transcript: &Vec<G::ScalarField>,
    r_shift: &G::ScalarField,
    kzg_challenge: &G::ScalarField,
) -> Result<G, Error> {
    let ck_polynomial = DensePolynomial::from_coefficients_slice(
        &polynomial_coefficients_from_transcript(transcript, r_shift),
    );
    assert_eq!(srs_powers.len(), ck_polynomial.coeffs.len());

    let eval = start_timer!(|| "polynomial eval");
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(&transcript, kzg_challenge, &r_shift);
    end_timer!(eval);
    let quotient = start_timer!(|| "polynomial quotient");
    let quotient_polynomial = &(&ck_polynomial
        - &DensePolynomial::from_coefficients_vec(vec![ck_polynomial_c_eval]))
        / &(DensePolynomial::from_coefficients_vec(vec![
            -kzg_challenge.clone(),
            <G::ScalarField>::one(),
        ]));
    end_timer!(quotient);

    let mut quotient_polynomial_coeffs = quotient_polynomial.coeffs;
    quotient_polynomial_coeffs.resize(srs_powers.len(), <G::ScalarField>::zero());

    let multiexp = start_timer!(|| "opening multiexp");
    let opening =
        MultiexponentiationInnerProduct::inner_product(srs_powers, &quotient_polynomial_coeffs);
    end_timer!(multiexp);
    opening
}

//TODO: Figure out how to avoid needing two separate methods for verification of opposite groups
// pub fn verify_commitment_key_g2_kzg_opening<P: PairingEngine>(
//     v_srs: &VerifierSRS<P>,
//     ck_final: &P::G1Projective,
//     ck_opening: &P::G1Projective,
//     transcript: &Vec<P::Fr>,
//     r_shift: &P::Fr,
//     kzg_challenge: &P::Fr,
// ) -> Result<bool, Error> {
//     let ck_polynomial_c_eval =
//         polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
//     Ok(P::pairing(
//         v_srs.g,
//         *ck_final - &v_srs.h.mul(ck_polynomial_c_eval.into_repr()),
//     ) == P::pairing(
//         v_srs.h_beta - &v_srs.h_alpha.mul(kzg_challenge.into_repr()),
//         *ck_opening,
//     ))
// }
pub fn verify_commitment_key_g2_kzg_opening<P: PairingEngine>(
    v_srs: &VerifierSRS<P>,
    ck_final: &P::G1Projective,
    ck_opening: &P::G1Projective,
    transcript: &Vec<P::Fr>,
    r_shift: &P::Fr,
    kzg_challenge: &P::Fr,
) -> Result<bool, Error> {
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
    Ok(P::pairing(
        *ck_final - &v_srs.g.mul(ck_polynomial_c_eval.into_repr()),
        v_srs.h,
    ) == P::pairing(
        *ck_opening,
        v_srs.h_beta - &v_srs.h.mul(kzg_challenge.into_repr()),
    ))
}
pub fn verify_commitment_key_g1_kzg_opening<P: PairingEngine>(
    v_srs: &VerifierSRS<P>,
    ck_final: &P::G1Projective,
    ck_opening: &P::G1Projective,
    transcript: &Vec<P::Fr>,
    r_shift: &P::Fr,
    kzg_challenge: &P::Fr,
) -> Result<bool, Error> {
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
    Ok(P::pairing(
        *ck_final - &v_srs.g.mul(ck_polynomial_c_eval.into_repr()),
        v_srs.h,
    ) == P::pairing(
        *ck_opening,
        v_srs.h_alpha - &v_srs.h.mul(kzg_challenge.into_repr()),
    ))
}

pub fn structured_generators_scalar_power<G: ProjectiveCurve>(
    num: usize,
    g: &G,
    s: &G::ScalarField,
) -> Vec<G> {
    assert!(num > 0);
    let mut powers_of_scalar = vec![];
    let mut pow_s = G::ScalarField::one();
    for _ in 0..num {
        powers_of_scalar.push(pow_s);
        pow_s *= s;
    }

    let window_size = FixedBaseMSM::get_mul_window_size(num);

    let scalar_bits = G::ScalarField::size_in_bits();
    let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g.clone());
    let powers_of_g =
        FixedBaseMSM::multi_scalar_mul::<G>(scalar_bits, window_size, &g_table, &powers_of_scalar);
    powers_of_g
}

fn polynomial_evaluation_product_form_from_transcript<F: Field>(
    transcript: &Vec<F>,
    z: &F,
    r_shift: &F,
) -> F {
    let mut power_2_zr = (z.clone() * z) * r_shift;
    let mut product_form = Vec::new();
    for x in transcript.iter() {
        product_form.push(F::one() + (x.clone() * &power_2_zr));
        power_2_zr *= power_2_zr;
    }
    product_form.iter().product()
}

fn polynomial_coefficients_from_transcript<F: Field>(transcript: &Vec<F>, r_shift: &F) -> Vec<F> {
    let mut coefficients = vec![F::one()];
    let mut power_2_r = r_shift.clone();
    for (i, x) in transcript.iter().enumerate() {
        for j in 0..(2_usize).pow(i as u32) {
            coefficients.push(coefficients[j] * &(x.clone() * &power_2_r));
        }
        power_2_r *= power_2_r;
    }
    // Interleave with 0 coefficients
    coefficients
        .iter()
        .interleave(vec![F::zero()].iter().cycle().take(coefficients.len() - 1))
        .cloned()
        .collect()
}