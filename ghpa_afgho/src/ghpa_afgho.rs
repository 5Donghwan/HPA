extern crate ark_ff;
// extern crate ark_ec;
// use self::ark_ec::PairingEngine;

use self::ark_ff::{to_bytes, Field, One, Zero, UniformRand};
extern crate ark_serialize;
use self::ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};
extern crate ark_std;
use self::ark_std::rand::Rng;
use self::ark_std::{end_timer, start_timer};
extern crate digest;
use self::digest::Digest;
use std::{convert::TryInto, f32, marker::PhantomData, ops::MulAssign};

use crate::{mul_helper, Error, InnerProductArgumentError};
extern crate ark_dh_commitments;
use self::ark_dh_commitments::DoublyHomomorphicCommitment;
extern crate ark_inner_products;
use self::ark_inner_products::InnerProduct;
use self::ark_std::cfg_iter;

use std::fmt;

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

pub struct HPA<IP, LMC, RMC, IPC, D> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _right_commitment: PhantomData<RMC>,
    _inner_product_commitment: PhantomData<IPC>,
    _digest: PhantomData<D>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct HPAProof<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    // CM: DoublyHomomorphicCommitment,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    pub(crate) r_commitment_steps: Vec<(
        (IP::Output, IP::Output, IP::Output, IP::Output),
        (IP::Output, IP::Output, IP::Output, IP::Output),
        (IP::Output, IP::Output, IP::Output, IP::Output),
        (IP::Output, IP::Output, IP::Output, IP::Output),
        // (LMC::Output, LMC::Output, RMC::Output, RMC::Output),
        // (LMC::Output, LMC::Output, LMC::Output, LMC::Output),
    )>,
    pub(crate) r_d1_x: Vec<IP::Output>,
    pub(crate) r_d2_x: Vec<IP::Output>,
    // pub(crate) r_d1_x: Vec<LMC::Output>,
    // pub(crate) r_d2_x: Vec<RMC::Output>,

    pub(crate) e1: Vec<IP::LeftMessage>,
    pub(crate) e2: Vec<IP::RightMessage>,
    pub(crate) q1: IP::Output,
    pub(crate) q2: IP::Output,
    pub(crate) q3: IP::Output,
    pub(crate) q4: IP::Output,
    pub(crate) p1: IP::Output,
    pub(crate) p2: IP::Output,
    pub(crate) r: LMC::Scalar,

    _hpa: PhantomData<HPA<IP, LMC, RMC, IPC, D>>,
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
pub struct HPAAux<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    // CM: DoublyHomomorphicCommitment,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    //  pub(crate) r_transcript: Vec<(CM::Scalar, CM::Scalar, CM::Scalar, LMCMC::Scalar)>,
    // pub(crate) ck_base: (LMC::Key, RMC::Key),
    _hpa: PhantomData<HPA<IP, LMC, RMC, IPC, D>>,
}

//TODO: Can extend HPA to support "identity commitments" in addition to "compact commitments", i.e. for SIPP

impl<IP, LMC, RMC, IPC, D> HPA<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    // CM: DoublyHomomorphicCommitment,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    LMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
    // CM::Message: MulAssign<CM::Message>,
    // CM::Message: MulAssign<CM::Scalar>,
    // CM::Key: MulAssign<CM::Scalar>,
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
{

    pub fn set_values(
        l: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        r: &[<RMC as DoublyHomomorphicCommitment>::Scalar],
        generator_g1: &IP::LeftMessage,
        generator_g2: &IP::RightMessage,
    ) -> Result<
        (
            Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>,
            Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>
        ), Error >
    {
        let mut v1 = Vec::new();
        let mut v2 = Vec::new();
        let mut u1 = Vec::new();
        let mut u2 = Vec::new();

        let len = l.len();
        let one = <LMC as DoublyHomomorphicCommitment>::Scalar::one();
        let g2_one = mul_helper(&generator_g2.clone(), &one);
        for i in 0..len {
            v1.push(mul_helper(&generator_g1.clone(), &l[i]));
            v2.push(mul_helper(&generator_g2.clone(), &r[i]));
            u1.push(mul_helper(&generator_g1.clone(), &(l[i] * r[i])));
            u2.push(g2_one.clone());
        }
        // let u1 = cfg_iter!(l)
        //             .zip(r)
        //             .map(|(b_1, b_2)| *b_1 * b_2.clone())
        //             .collect::<Vec<CM::Scalar>>();
        
        Ok((v1, v2, u1, u2))
    }

    pub fn init_commit<R: Rng>(
        left_value: &Vec<IP::LeftMessage>,
        right_value: &Vec<IP::RightMessage>,
        gamma1: &Vec<IP::RightMessage>,
        gamma2: &Vec<IP::LeftMessage>,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            LMC::Scalar,
            Vec<LMC::Scalar>,
            LMC::Scalar,
            LMC::Scalar,
            Vec<IP::LeftMessage>,
            // c, d1, d2, x, d3, gm, gm_vec, r_c, r_x, w_vec,
        ),
        Error,
    > {
        let l = left_value.clone();
        let r = right_value.clone();
        let gamma1 = gamma1.clone();
        let gamma2 = gamma2.clone();
        let h1 = h1.clone();
        let h2 = h2.clone();
        let ht = IP::inner_product(&h1, &h2).unwrap();

        let r_c = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_x = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

        let c = IP::inner_product(&l, &r)? + mul_helper(&ht, &r_c);
        let d1 = IP::inner_product(&l, &gamma1).unwrap();
        let d2 = IP::inner_product(&gamma2, &r)?;

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
        // let gm = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let gm = 'challenge: loop {
            let mut hash_input = Vec::new();
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![c, d1, d2]?);
            let gm: LMC::Scalar =
                u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
                    .into();
            break 'challenge gm;
        };
        let mut gm_vec = Vec::new();
        gm_vec.push(<LMC as DoublyHomomorphicCommitment>::Scalar::one());
        for i in 1..l.len() {
            gm_vec.push(gm_vec[i - 1] * gm);
        }
        let mut w_vec = Vec::new();
        for i in 0..l.len() {
            w_vec.push( mul_helper(&l[i], &gm_vec[i]) );
        }
        
        let x = IP::inner_product(&w_vec, &r).unwrap() + mul_helper(&ht, &r_x);
        let d3 = IP::inner_product(&w_vec, &gamma1).unwrap();

        Ok((
            c, d1, d2, x, d3, gm, gm_vec, r_c, r_x, w_vec,
        ))
    }

    pub fn init_commit2<R: Rng>(
        left_value: &Vec<IP::LeftMessage>,
        right_value: &Vec<IP::RightMessage>,
        gamma1: &Vec<IP::RightMessage>,
        gamma2: &Vec<IP::LeftMessage>,
        gm_vec: &Vec<<LMC as DoublyHomomorphicCommitment>::Scalar>,
        r_x: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            LMC::Scalar,
            Vec<IP::LeftMessage>,
            // c, d1, d2, x, d3, r_c, w_vec,
        ),
        Error,
    > {
        let l = left_value.clone();
        let r = right_value.clone();
        let gamma1 = gamma1.clone();
        let gamma2 = gamma2.clone();
        let h1 = h1.clone();
        let h2 = h2.clone();
        let ht = IP::inner_product(&h1, &h2).unwrap();


        let r_c = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);


        let c = IP::inner_product(&l, &r)? + mul_helper(&ht, &r_c);
        let d1 = IP::inner_product(&l, &gamma1).unwrap();
        let d2 = IP::inner_product(&gamma2, &r).unwrap();

        let mut w_vec = Vec::new();
        for i in 0..l.len() {
            w_vec.push( mul_helper(&l[i],&gm_vec[i]));
        }

        let x = IP::inner_product(&w_vec, &r).unwrap() + mul_helper(&ht, &r_x);
        let d3 = IP::inner_product(&w_vec, &gamma1).unwrap();

        Ok((
            c, d1, d2, x, d3, r_c, w_vec,
        ))
    }

    pub fn setup<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<(Vec<LMC::Key>, Vec<RMC::Key>), Error> {
        // println!("size : {}", size);
        let gamma1 = LMC::setup(rng, size)?;
        let gamma2 = RMC::setup(rng, size)?;
        // println!("gamma1_len : {}", gamma1.len());
        // println!("gamma2_len : {}", gamma2.len());

        Ok((gamma1, gamma2))
    }

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
        values: (
            &[IP::LeftMessage],
            &[IP::RightMessage],
            &[IP::LeftMessage],
        ),
        // srs: &HPASRS<IP, CM, D>,
        ck_message: (&[RMC::Message], &[LMC::Message]),
        // com: (&IP::Output, &IP::Output, &IP::Output),
        witness: (
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<HPAProof<IP, LMC, RMC, IPC, D>, Error> {
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
            (values.0, values.1, values.2),// values.3),
            // srs,
            (ck_message.0, ck_message.1),
            witness,
            gm,
            h1, h2,
            rng,
        )?;
        Ok(proof)
    }

    pub fn verify(
        ck_message: (Vec<RMC::Message>, Vec<LMC::Message>),
        com: (&IP::Output, &IP::Output, &IP::Output, &IP::Output, &IP::Output), // com ( c, x, d1, d2, d3 )
        proof: &mut HPAProof<IP, LMC, RMC, IPC, D>,
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
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
        let h1 = h1.clone();
        let h2 = h2.clone();
        let ht = IP::inner_product(&h1, &h2).unwrap();
        
        let (mut transcript, ch_c) = Self::_compute_recursive_challenges(proof, gm)?;

        let (mut gamma1, mut gamma2) = ck_message.clone();

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
                // println!("check");
                // Verifier's work in reduce
                let split = gamma1.len() / 2;

                let last_commitment = proof.r_commitment_steps.pop().unwrap();
                let last_transcript = transcript.pop().unwrap();
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

                let gamma1_l = &gamma1[..split];
                let gamma1_r = &gamma1[split..];
                let gamma2_l = &gamma2[..split];
                let gamma2_r = &gamma2[split..];
                

                let alpha = last_transcript.0;
                // println!("ver al : {}", alpha);
                // let alpha_inv = last_transcript.1;
                let gm_inv = last_transcript.2;


                gamma1 = cfg_iter!(gamma1_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(gamma1_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<RMC::Message>>();
                gamma2 = cfg_iter!(gamma2_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(gamma2_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<LMC::Message>>();

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
                    
                    let zero = <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
                    let minus_one = zero - <LMC as DoublyHomomorphicCommitment>::Scalar::one();

                    let r = mul_helper(&ht, &proof.r);
                    temp_right = temp_right + mul_helper(&r, &minus_one);
                    if temp_left == temp_right {
                        result1 = true;
                    }
                    // println!("g : {}, e: {}", gamma1.len(), e1.len());

                    let temp_left = IP::inner_product(&e1, &gamma1).unwrap();
                    let temp_right = mul_helper(&d1_prime, &ch_c_vec[0]) + mul_helper(&d3_prime, &(ch_c_vec[1] + ch_c_vec[2])) + proof.p1.clone();
                    if temp_left == temp_right {
                        result2 = true;
                    }

                    let temp_left = IP::inner_product(&gamma2, &e2).unwrap();
                    let temp_right = mul_helper(&d2_prime, &ch_c_vec[3]) + proof.p2.clone();
                    if temp_left == temp_right {
                        result3 = true;
                    }

                    // let ch_d = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                    // let ch_d_inv = ch_d.inverse().unwrap();

                    // check pairing equation
                    // let kai_scalar =
                    //     IP::inner_product(&(gamma1[..1].to_vec()), &(gamma2[..1].to_vec()))?;

                    // let temp1 = mul_helper(&(gamma1[0]), &ch_d);
                    // e1[0] = e1[0].clone() + temp1;
                    // e2[0] = e2[0].clone() + mul_helper(&(gamma2[0]), &(ch_d_inv));

                    // let left = IP::inner_product(&e1, &e2)?;
                    
                    // let one = <LMC as DoublyHomomorphicCommitment>::Scalar::one();
                    // let zero: <LMC as DoublyHomomorphicCommitment>::Scalar =
                    //     <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
                    // let minus_one = zero - one;

                    // let mut ch_c_vec = Vec::new();
                    // ch_c_vec.push(ch_c.clone()); // ch_c_vec = {c^1, c^2, c^3, ..., c^11}
                    // for i in 1..11{
                    //     ch_c_vec.push(ch_c_vec[i-1] * ch_c);
                    // }

                    // let mut temp2 = c_prime.clone() + mul_helper(&x_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&y_prime, &ch_c_vec[8]);
                    // temp2 = mul_helper(&temp2, &ch_c_vec[1]);
                    // let mut temp3 = d1_prime.clone() + mul_helper(&d3_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&d4_prime, &ch_c_vec[8]);
                    // temp3 = mul_helper(&temp3, &(ch_c * ch_d_inv));
                    // let mut temp4 = proof.q1.clone() + mul_helper(&proof.q2, &ch_c_vec[2]) + mul_helper(&proof.q3, &ch_c_vec[5]) + mul_helper(&proof.p3, &ch_c_vec[8]);
                    // temp4 = mul_helper(&temp4, &ch_c);
                    // let temp5 = ch_d_inv * proof.r1 + ch_d * proof.r2 + proof.r3;
                    // let mut temp6 = mul_helper(&srs.ht, &temp5);
                    // temp6 = mul_helper(&temp6, &minus_one);

                    // let right = kai_scalar + temp2 + mul_helper(&d2_prime, &(ch_c * ch_d)) + temp3 + temp4 + mul_helper(&proof.p1, &ch_d_inv)
                    //     + mul_helper(&proof.p2, &ch_d) + proof.r.clone() + temp6;
                    
                    result = result1 && result2 && result3;
                }
            }
            Ok(result)
        } else {
            
            // let e1 = proof.e1.clone();
            // let e2 = proof.e2.clone();

            // // let ch_d = <CM as DoublyHomomorphicCommitment>::Scalar::rand(rng);
            // // let ch_d_inv = ch_d.inverse().unwrap();

            // // check pairing equation
            // // let kai_scalar = IP::inner_product(&(gamma1[..1].to_vec()), &(gamma2[..1].to_vec()))?;

            // // e1[0] = e1[0].clone() + mul_helper(&(gamma1[0]), &ch_d);
            // // e2[0] = e2[0].clone() + mul_helper(&(gamma2[0]), &(ch_d_inv));
            // // let left = IP::inner_product(&e1, &e2)?;

            // // let one = <LMC as DoublyHomomorphicCommitment>::Scalar::one();
            // // let zero: <LMC as DoublyHomomorphicCommitment>::Scalar =
            // //     <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
            // // let minus_one = zero - one;

            // let mut ch_c_vec = Vec::new();
            // ch_c_vec.push(ch_c.clone()); // ch_c_vec = {c^1, c^2, c^3, ..., c^11}
            // for i in 1..7{
            //     ch_c_vec.push(ch_c_vec[i-1] * ch_c);
            // }



            // let mut temp2 = c_prime.clone() + mul_helper(&x_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&y_prime, &ch_c_vec[8]);
            // temp2 = mul_helper(&temp2, &ch_c_vec[1]);
            // let mut temp3 = d1_prime.clone() + mul_helper(&d3_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&d4_prime, &ch_c_vec[8]);
            // temp3 = mul_helper(&temp3, &(ch_c * ch_d_inv));
            // let mut temp4 = proof.q1.clone() + mul_helper(&proof.q2, &ch_c_vec[2]) + mul_helper(&proof.q3, &ch_c_vec[5]) + mul_helper(&proof.p3, &ch_c_vec[8]);
            // temp4 = mul_helper(&temp4, &ch_c);
            // let temp5 = ch_d_inv * proof.r1 + ch_d * proof.r2 + proof.r3;
            // let mut temp6 = mul_helper(&srs.ht, &temp5);
            // temp6 = mul_helper(&temp6, &minus_one);

            // let right = kai_scalar + temp2 + mul_helper(&d2_prime, &(ch_c * ch_d)) + temp3 + temp4 + mul_helper(&proof.p1, &ch_d_inv)
            //     + mul_helper(&proof.p2, &ch_d) + proof.r.clone() + temp6;
            
            // result = left == right;
            Ok(result)
        }
    }

    pub fn prove_with_aux<R: Rng>(
        values: (
            &[IP::LeftMessage],
            &[IP::RightMessage],
            &[IP::LeftMessage],
            // &[IP::LeftMessage],
        ),
        // srs: &HPASRS<IP, CM, D>,
        ck_message: (&[RMC::Message], &[LMC::Message]),
        witness: (
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
            // &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<(HPAProof<IP, LMC, RMC, IPC, D>, HPAAux<IP, LMC, RMC, IPC, D>), Error> {
        let (v1, v2, w_vec) = values;
        // let (gamma1, gamma2) = ck;
        let (gamma1_message, gamma2_message) = ck_message;
        Self::_prove(
            &(v1.to_vec(), v2.to_vec(), w_vec.to_vec()),// k_vec.to_vec()),
            // srs,
            (gamma1_message.to_vec(), gamma2_message.to_vec()),
            witness,
            gm,
            h1, h2,
            rng,
        )
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove<R: Rng>(
        values: &(
            Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>,
            Vec<IP::LeftMessage>,
        ),
        ck_message: (Vec<RMC::Message>, Vec<LMC::Message>),
        witness: (
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<(HPAProof<IP, LMC, RMC, IPC, D>, HPAAux<IP, LMC, RMC, IPC, D>), Error> {
        let (mut v1, mut v2, mut w_vec) = values.clone();
        let (mut gamma1_message, mut gamma2_message) = ck_message.clone();
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        let mut r_d1_x = Vec::new();
        let mut r_d2_x = Vec::new();
        assert!(v1.len().is_power_of_two());

        let mut r_c = witness.0.clone();
        let mut r_x = witness.1.clone();
        let h1 = h1.clone();
        let h2 = h2.clone();
        let ht = IP::inner_product(&h1, &h2).unwrap();

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

                let r_cl = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_cx = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_xl = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_x_plus = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_x_minus = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

                let zero = <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
                let minus_one = zero - <LMC as DoublyHomomorphicCommitment>::Scalar::one();

                let r_cr = r_c + mul_helper(&r_cl, &minus_one);
                let r_xr = r_x + mul_helper(&r_xl, &minus_one);

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

                let c_l = IP::inner_product(&v1_l, &v2_l).unwrap() + mul_helper(&ht, &r_cl);
                let c_r = IP::inner_product(&v1_r, &v2_r).unwrap() + mul_helper(&ht, &r_cr);
                let c_x = IP::inner_product(&v1_l, &v2_r).unwrap() + IP::inner_product(&v1_r, &v2_l).unwrap() + mul_helper(&ht, &r_cx);
                let x_l = IP::inner_product(&w_vec_l, &v2_l).unwrap() + mul_helper(&ht, &r_xl);
                let x_r = IP::inner_product(&w_vec_r, &v2_r).unwrap() + mul_helper(&ht, &r_xr);
                let x_plus = IP::inner_product(&w_vec_l, &v2_r).unwrap() + mul_helper(&ht, &r_x_plus);
                let x_minus = IP::inner_product(&w_vec_r, &v2_l).unwrap() + mul_helper(&ht, &r_x_minus);
                let d1_l = IP::inner_product(&v1_l, &gamma1_l).unwrap();
                let d1_r = IP::inner_product(&v1_r, &gamma1_r).unwrap();
                let d1_x = IP::inner_product(&v1_r, &gamma1_l).unwrap() + IP::inner_product(&v1_l, &gamma1_r).unwrap();
                // let d1_plus = CM::commit(&gamma1_r, &v1_l).unwrap();
                // let d1_minus = CM::commit(&gamma1_l, &v1_r).unwrap();
                let d2_l = IP::inner_product(&gamma2_l, &v2_l).unwrap();
                let d2_r = IP::inner_product(&gamma2_r, &v2_r).unwrap();
                let d2_x = IP::inner_product(&gamma2_l, &v2_r).unwrap() + IP::inner_product(&gamma2_r, &v2_l).unwrap();
                // let d2_plus = CM::commit(&gamma2_r, &v2_l).unwrap();
                // let d2_minus = CM::commit(&gamma2_l, &v2_r).unwrap();
                let d3_l = IP::inner_product(&w_vec_l, &gamma1_l).unwrap();
                let d3_r = IP::inner_product(&w_vec_r, &gamma1_r).unwrap();
                let d3_plus = IP::inner_product(&w_vec_l, &gamma1_r).unwrap();
                let d3_minus = IP::inner_product(&w_vec_r, &gamma1_l).unwrap();


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
                    let alpha: LMC::Scalar = u128::from_be_bytes(
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
                    .collect::<Vec<LMC::Message>>();
                end_timer!(rescale_v1);

                let rescale_v2 = start_timer!(|| "Rescale V2");
                v2 = cfg_iter!(v2_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(v2_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<RMC::Message>>();
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
                    .collect::<Vec<LMC::Message>>();
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
                    .collect::<Vec<RMC::Message>>();
                gamma2_message = cfg_iter!(gamma2_l)
                    .map(|b| mul_helper(b, &alpha))
                    .zip(gamma2_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<LMC::Message>>();
                

                let com1 = (c_l, c_r, c_x.clone(), c_x.clone());
                let com2 = (x_l, x_r, x_plus, x_minus);
                let com3 = (d1_l, d1_r, d2_l, d2_r);
                let com4 = (d3_l, d3_r, d3_plus, d3_minus);

                r_commitment_steps.push((com1, com2, com3, com4));
                r_transcript.push((alpha, alpha_inv, gm_inv));
                r_d1_x.push(d1_x);
                r_d2_x.push(d2_x);

                end_timer!(recurse);
            }
        };

        
        let v1_val = v1[0].clone();
        let v2_val = v2[0].clone();
        let w_val = w_vec[0].clone();

        // println!("v1 ?= w : {}", v1_val == w_val);

        let r_d1 = <LMC::Scalar>::rand(rng);
        let r_d2 = <LMC::Scalar>::rand(rng);
        let r_q1 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_q2 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_q3 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_q4 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

         let mut d1 = Vec::new();
         let mut d2 = Vec::new();
         d1.push(mul_helper(&h1[0], &r_d1));
         d2.push(mul_helper(&h2[0], &r_d2));

        let q1 = IP::inner_product(&d1, &d2).unwrap() + mul_helper(&ht, &r_q1);
        let q2 = IP::inner_product(&d1, &v2).unwrap() + mul_helper(&ht, &r_q2);
        let q3 = IP::inner_product(&v1, &d2).unwrap() + mul_helper(&ht, &r_q3);
        let q4 = IP::inner_product(&w_vec, &d2).unwrap() + mul_helper(&ht, &r_q4);
        let p1 = IP::inner_product(&d1, &gamma1_message).unwrap();
        let p2 = IP::inner_product(&gamma2_message, &d2).unwrap();


        let ch_c = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&to_bytes![q1, q2, q3, q4, p1, p2]?);
            let ch_c: LMC::Scalar =
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


        let e1 = d1[0].clone() + mul_helper(&v1_val, &ch_c) + mul_helper(&v1_val, &ch_c_2) + mul_helper(&w_val, &ch_c_3);
        let e2 = d2[0].clone() + mul_helper(&v2_val, &ch_c_4);
        let r = r_q1 + mul_helper(&r_q2, &ch_c_4) + mul_helper(&r_q3, &(ch_c + ch_c_2)) + mul_helper(&r_q4, &ch_c_3)
            + mul_helper(&r_c, &ch_c_5) + mul_helper(&r_x, &(ch_c_6 + ch_c_7));

        let mut e1_vec = Vec::new();
        let mut e2_vec = Vec::new();
        e1_vec.push(e1);
        e2_vec.push(e2);
        
        // r_transcript.reverse();
        r_commitment_steps.reverse();
        r_d1_x.reverse();
        r_d2_x.reverse();

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

                _hpa: PhantomData,
            },
            HPAAux { _hpa: PhantomData },
        ))
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        proof: &HPAProof<IP, LMC, RMC, IPC, D>,
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar
    ) -> Result<
        (
            Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar)>,
            LMC::Scalar,
        ),
        Error,
    > {
        Self::_compute_recursive_challenges(proof, gm)
    }

    fn _compute_recursive_challenges(
        proof: &HPAProof<IP, LMC, RMC, IPC, D>,
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar
    ) -> Result<
        (
            Vec<(LMC::Scalar, LMC::Scalar,  LMC::Scalar)>,
            LMC::Scalar,
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
                let alpha: LMC::Scalar = u128::from_be_bytes(
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
        // println!("r_transcript len : {}", r_transcript.len());

        let ch_c = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&to_bytes![proof.q1, proof.q2, proof.q3, proof.q4, proof.p1, proof.p2]?);
            let ch_c: LMC::Scalar =
                u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
                    .into();

            break 'challenge ch_c;
        };

        Ok((r_transcript, ch_c))
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

impl<IP, LMC, RMC, IPC, D> Clone for HPAProof<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
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
            // r1: self.r1.clone(),
            // r2: self.r2.clone(),
            // r3: self.r3.clone(),
            // r_base: self.r_base.clone(),
            _hpa: PhantomData,
        }
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use ark_bls12_381::Bls12_381;
//     use ark_ec::PairingEngine;
//     use ark_ff::UniformRand;
//     use ark_std::rand::{rngs::StdRng, SeedableRng};
//     use blake2::Blake2b;

//     use ark_dh_commitments::{
//         afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
//         identity::IdentityCommitment,
//         pedersen::PedersenCommitment,
//         random_generators,
//     };
//     use ark_inner_products::{
//         ExtensionFieldElement, InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct,
//         ScalarInnerProduct,
//     };

//     type GC1 = AFGHOCommitmentG1<Bls12_381>;
//     type GC2 = AFGHOCommitmentG2<Bls12_381>;
//     type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
//     type SC2 = PedersenCommitment<<Bls12_381 as PairingEngine>::G2Projective>;
//     const TEST_SIZE: usize = 8;

//     #[test]
//     fn pairing_inner_product_test() {
//         type IP = PairingInnerProduct<Bls12_381>;
//         type IPC =
//             IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>;
//         type PairingDORY = DORY<IP, GC1, GC2, IPC, Blake2b>;

//         let mut rng = StdRng::seed_from_u64(0u64);
//         let (ck_a, ck_b, ck_t) = PairingDORY::setup(&mut rng, TEST_SIZE).unwrap();
//         let m_a = random_generators(&mut rng, TEST_SIZE);
//         let m_b = random_generators(&mut rng, TEST_SIZE);
//         let com_a = GC1::commit(&ck_a, &m_a).unwrap();
//         let com_b = GC2::commit(&ck_b, &m_b).unwrap();
//         let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
//         let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

//         let proof = PairingDORY::prove(
//             (&m_a, &m_b, &t[0]),
//             (&ck_a, &ck_b, &ck_t),
//             (&com_a, &com_b, &com_t),
//         )
//         .unwrap();

//         assert!(
//             PairingDORY::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
//         );
//     }

//     #[test]
//     fn multiexponentiation_inner_product_test() {
//         type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
//         type IPC = IdentityCommitment<
//             <Bls12_381 as PairingEngine>::G1Projective,
//             <Bls12_381 as PairingEngine>::Fr,
//         >;
//         type MultiExpDORY = DORY<IP, GC1, SC1, IPC, Blake2b>;

//         let mut rng = StdRng::seed_from_u64(0u64);
//         let (ck_a, ck_b, ck_t) = MultiExpDORY::setup(&mut rng, TEST_SIZE).unwrap();
//         let m_a = random_generators(&mut rng, TEST_SIZE);
//         let mut m_b = Vec::new();
//         for _ in 0..TEST_SIZE {
//             m_b.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
//         }
//         let com_a = GC1::commit(&ck_a, &m_a).unwrap();
//         let com_b = SC1::commit(&ck_b, &m_b).unwrap();
//         let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
//         let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

//         let proof = MultiExpDORY::prove(
//             (&m_a, &m_b, &t[0]),
//             (&ck_a, &ck_b, &ck_t),
//             (&com_a, &com_b, &com_t),
//         )
//         .unwrap();

//         assert!(
//             MultiExpDORY::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
//         );
//     }

//     #[test]
//     fn scalar_inner_product_test() {
//         type IP = ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>;
//         type IPC =
//             IdentityCommitment<<Bls12_381 as PairingEngine>::Fr, <Bls12_381 as PairingEngine>::Fr>;
//         type ScalarDORY = DORY<IP, SC2, SC2, IPC, Blake2b>;

//         let mut rng = StdRng::seed_from_u64(0u64);
//         let (ck_a, ck_b, ck_t) = ScalarDORY::setup(&mut rng, TEST_SIZE).unwrap();
//         let mut m_a = Vec::new();
//         let mut m_b = Vec::new();
//         for _ in 0..TEST_SIZE {
//             m_a.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
//             m_b.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
//         }
//         let com_a = SC2::commit(&ck_a, &m_a).unwrap();
//         let com_b = SC2::commit(&ck_b, &m_b).unwrap();
//         let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
//         let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

//         let proof = ScalarDORY::prove(
//             (&m_a, &m_b, &t[0]),
//             (&ck_a, &ck_b, &ck_t),
//             (&com_a, &com_b, &com_t),
//         )
//         .unwrap();

//         assert!(
//             ScalarDORY::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
//         );
//     }
// }
