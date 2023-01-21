extern crate ark_ff;
// use ark_ec::PairingEngine;

use self::ark_ff::{to_bytes, Field, One, UniformRand, Zero};
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
    )>,

    pub(crate) c_x: Vec<IP::Output>,
    pub(crate) x_plus: Vec<IP::Output>,
    pub(crate) x_minus: Vec<IP::Output>,
    pub(crate) y_plus: Vec<IP::Output>,
    pub(crate) y_minus: Vec<IP::Output>,

    pub(crate) e1: Vec<IP::LeftMessage>,
    pub(crate) e2: Vec<IP::RightMessage>,
    pub(crate) p1: IP::Output,
    pub(crate) p2: IP::Output,
    pub(crate) p3: IP::Output,
    pub(crate) q1: IP::Output,
    pub(crate) q2: IP::Output,
    pub(crate) q3: IP::Output,
    pub(crate) r: IP::Output,
    pub(crate) r1: LMC::Scalar,
    pub(crate) r2: LMC::Scalar,
    pub(crate) r3: LMC::Scalar,
    // pub(crate) r_base: (LMC::Message, RMC::Message),
    _hpa: PhantomData<HPA<IP, LMC, RMC, IPC, D>>,
}

#[derive(Clone)]
pub struct HPASRS<IP, LMC, RMC, IPC, D>
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
    pub(crate) delta1_l: Vec<IP::Output>,
    pub(crate) delta1_r: Vec<IP::Output>,
    pub(crate) delta2_l: Vec<IP::Output>,
    pub(crate) delta2_r: Vec<IP::Output>,
    pub(crate) kai: Vec<IP::Output>,
    pub(crate) ht: IP::Output,

    _hpa: PhantomData<HPA<IP, LMC, RMC, IPC, D>>,
}

#[derive(Clone)]
pub struct HPAAux<IP, LMC, RMC, IPC, D>
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
    // pub(crate) r_transcript: Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar)>,
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
    LMC::Output: MulAssign<LMC::Scalar>,
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
{

    pub fn set_values(
        l: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        r: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
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

        Ok((v1, v2, u1, u2))
    }

    pub fn init_commit<R: Rng>(
        left_value: &Vec<IP::LeftMessage>,
        right_value: &Vec<IP::RightMessage>,
        gamma1: &Vec<IP::LeftMessage>,
        gamma2: &Vec<IP::RightMessage>,
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
            IP::Output,
            IP::Output,
            LMC::Scalar,
            Vec<LMC::Scalar>,
            LMC::Scalar,
            LMC::Scalar,
            LMC::Scalar,
            LMC::Scalar,
            LMC::Scalar,
            LMC::Scalar,
            LMC::Scalar,
            Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>,
        ),
        Error,
    > {
        let l = left_value.clone();
        let r = right_value.clone();
        let gamma1 = gamma1.clone();
        let gamma2 = gamma2.clone();
        let h1 = h1.clone();
        let h2 = h2.clone();
        let ht = IP::inner_product(&h1, &h2)?;

        let r_c = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_d1 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_d2 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

        let c = IP::inner_product(&l, &r)? + mul_helper(&ht, &r_c);
        let d1 = IP::inner_product(&l, &gamma2)? + mul_helper(&ht, &r_d1);
        let d2 = IP::inner_product(&gamma1, &r)? + mul_helper(&ht, &r_d2);

        // Fiat-Schamir challenge
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
            w_vec.push(mul_helper(&l[i], &gm_vec[i]));
        }
        let mut k_vec = Vec::new();
        let zero = <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
        let temp = <IP::LeftMessage>::rand(rng);
        let zero_g1 = mul_helper(&temp, &zero);
        for _ in 0..l.len() {
            k_vec.push(zero_g1.clone());
        }

        let r_x = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_y = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_d3 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_d4 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

        let x = IP::inner_product(&w_vec, &r).unwrap() + mul_helper(&ht, &r_x);
        let y = IP::inner_product(&k_vec, &r).unwrap() + mul_helper(&ht, &r_y);
        let d3 = IP::inner_product(&w_vec, &gamma2).unwrap() + mul_helper(&ht, &r_d3);
        let d4 = IP::inner_product(&k_vec, &gamma2).unwrap() + mul_helper(&ht, &r_d4);

        Ok((
            c, d1, d2, x, y, d3, d4, gm, gm_vec, r_c, r_d1, r_d2, r_x, r_y, r_d3, r_d4, w_vec,
            k_vec,
        ))
    }

    pub fn init_commit2<R: Rng>(
        left_value: &Vec<IP::LeftMessage>,
        right_value: &Vec<IP::RightMessage>,
        gamma1: &Vec<IP::LeftMessage>,
        gamma2: &Vec<IP::RightMessage>,
        h1: &Vec<IP::LeftMessage>,
        h2: &Vec<IP::RightMessage>,
        // gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        gm_vec: &Vec<<LMC as DoublyHomomorphicCommitment>::Scalar>,
        r_c: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_d1: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_d2: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_x: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_y: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_d3: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        r_d4: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            // LMC::Scalar,
            // Vec<LMC::Scalar>,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            // LMC::Scalar,
            Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>,
        ),
        Error,
    > {
        let l = left_value.clone();
        let r = right_value.clone();
        let gamma1 = gamma1.clone();
        let gamma2 = gamma2.clone();
        let h1 = h1.clone();
        let h2 = h2.clone();
        let ht = IP::inner_product(&h1, &h2)?;

        // let r_c = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        // let r_d1 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        // let r_d2 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

        let c = IP::inner_product(&l, &r)? + mul_helper(&ht, &r_c);
        let d1 = IP::inner_product(&l, &gamma2)? + mul_helper(&ht, &r_d1);
        let d2 = IP::inner_product(&gamma1, &r)? + mul_helper(&ht, &r_d2);

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
            w_vec.push(mul_helper(&l[i], &gm_vec[i]));
        }
        let mut k_vec = Vec::new();
        let zero = <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
        let temp = <IP::LeftMessage>::rand(rng);
        let zero_g1 = mul_helper(&temp, &zero);
        for _ in 0..l.len() {
            k_vec.push(zero_g1.clone());
        }

        // let r_x = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        // let r_y = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        // let r_d3 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        // let r_d4 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

        let x = IP::inner_product(&w_vec, &r).unwrap() + mul_helper(&ht, &r_x);
        let y = IP::inner_product(&k_vec, &r).unwrap() + mul_helper(&ht, &r_y);
        let d3 = IP::inner_product(&w_vec, &gamma2).unwrap() + mul_helper(&ht, &r_d3);
        let d4 = IP::inner_product(&k_vec, &gamma2).unwrap() + mul_helper(&ht, &r_d4);

        Ok((
            c, d1, d2, x, y, d3, d4, w_vec,
            k_vec,
        ))
    }

pub fn batch_commit<R: Rng>(
        r_c: &LMC::Scalar, r_x: &LMC::Scalar, r_y:&LMC::Scalar,
        r_d1: &LMC::Scalar, r_d2: &LMC::Scalar, r_d3: &LMC::Scalar, r_d4: &LMC::Scalar,
        v1: &Vec<IP::LeftMessage>, u1: &Vec<IP::LeftMessage>,
        v2: &Vec<IP::RightMessage>, u2: &Vec<IP::RightMessage>,
        w_vec: &Vec<IP::LeftMessage>, w_vec_: &Vec<IP::LeftMessage>,
        k_vec: &Vec<IP::LeftMessage>, k_vec_: &Vec<IP::LeftMessage>,
        h1: &Vec<IP::LeftMessage>, h2: &Vec<IP::RightMessage>,
        rng: &mut R,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            IP::Output,
            Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>,
            LMC::Scalar, LMC::Scalar, LMC::Scalar,
            LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar,
            Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>,
            LMC::Scalar,
            // z_c, z_x, z_y, v1, v2, r_c, r_x, r_y, r_d1, r_d2, r_d3, r_d4, w_vec, k_vec, delta
        ),
        Error,
    > {
        let ht = IP::inner_product(&h1, &h2).unwrap();

        let r_zc = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_zx = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_zy = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        
        let z_c = IP::inner_product(&v1, &u2).unwrap() + IP::inner_product(&u1, &v2).unwrap() + mul_helper(&ht, &r_zc);
        let z_x = IP::inner_product(&w_vec, &u2).unwrap() + IP::inner_product(&w_vec_, &v2).unwrap() + mul_helper(&ht, &r_zx);
        let z_y = IP::inner_product(&k_vec, &u2).unwrap() + IP::inner_product(&k_vec_, &v2).unwrap() + mul_helper(&ht, &r_zy);


        let delta = 'challenge: loop {
            let mut hash_input = Vec::new();
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![z_c, z_x]?);
            let delta: LMC::Scalar =
                u128::from_be_bytes(D::digest(&hash_input).as_slice()[0..16].try_into().unwrap())
                    .into();
            break 'challenge delta;
        };
        let delta_sqr = delta * delta;

        let mut bat_v1 = Vec::new();
        let mut bat_v2 = Vec::new();
        let mut bat_w_vec = Vec::new();
        let mut bat_k_vec = Vec::new();

        for i in 0..v1.len(){
            bat_v1.push(mul_helper(&v1[i], &delta) + u1[i].clone());
        }
        for i in 0..v1.len(){
            bat_v2.push(mul_helper(&v2[i], &delta) + u2[i].clone());
        }
        for i in 0..v1.len(){
            bat_w_vec.push(mul_helper(&w_vec[i], &delta) + w_vec_[i].clone());
        }
        for i in 0..v1.len(){
            bat_k_vec.push(mul_helper(&k_vec[i], &delta) + k_vec_[i].clone());
        }

        let bat_r_c = delta_sqr * r_c + delta * r_zc + r_c;
        let bat_r_x = delta_sqr * r_x.clone() + delta * r_zx + r_x;
        let bat_r_y = delta_sqr * r_y.clone() + delta * r_zx + r_y;
        let bat_r_d1 = delta * r_d1 + r_d1;
        let bat_r_d2 = delta * r_d2 + r_d2;
        let bat_r_d3 = delta * r_d3 + r_d3;
        let bat_r_d4 = delta * r_d4 + r_d4;


        Ok((
            z_c, z_x, z_y, 
            bat_v1, bat_v2, 
            bat_r_c, bat_r_x, bat_r_y, 
            bat_r_d1, bat_r_d2, bat_r_d3, bat_r_d4, 
            bat_w_vec, bat_k_vec, 
            delta
        ))
    }


    pub fn setup<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<(Vec<LMC::Key>, Vec<RMC::Key>), Error> {
        let gamma1 = RMC::setup(rng, size)?;
        let gamma2 = LMC::setup(rng, size)?;

        Ok((gamma2, gamma1))
    }

    pub fn precompute(
        ck_message: (&[LMC::Message], &[RMC::Message]),
        h1: &[LMC::Message],
        h2: &[RMC::Message],
    ) -> Result<HPASRS<IP, LMC, RMC, IPC, D>, Error> {
        // loop : until ck.len() >= 1
        let (mut gamma1, mut gamma2) = ck_message.clone();
        let h1 = h1.clone();
        let h2 = h2.clone();
        // let mut i = ck_message.0.len();
        let mut delta1_l = Vec::new();
        let mut delta1_r = Vec::new();
        let mut delta2_r = Vec::new();
        let mut kai = Vec::new();
        let mut split = ck_message.0.len() / 2;
        while split >= 1 {
            kai.push(IP::inner_product(gamma1, gamma2)?);
            // Generate gamma1, gamma2, gamma1_prime, gamma2_prime
            let gamma1_prime = &gamma1[..split];
            let gamma2_prime = &gamma2[..split];
            // Generate gamma1_R, gamma2_R (replace gamma1_L to gamma1_prime)
            let gamma1_r = &gamma1[split..];
            let gamma2_r = &gamma2[split..];

            // Compute delta1_L, delta1_R, delta2_L, delta2_R, kai
            delta1_l.push(IP::inner_product(gamma1_prime, gamma2_prime)?);
            delta1_r.push(IP::inner_product(gamma1_r, gamma2_prime)?);
            delta2_r.push(IP::inner_product(gamma1_prime, gamma2_r)?);

            split = split / 2;
            gamma1 = gamma1_prime;
            gamma2 = gamma2_prime;
            // i = i/2;
        }
        let mut delta2_l = delta1_l.clone();
        delta1_l.reverse();
        delta1_r.reverse();
        delta2_l.reverse();
        delta2_r.reverse();
        kai.reverse();

        let ht = IP::inner_product(&h1, &h2)?;

        Ok(HPASRS {
            delta1_l: delta1_l,
            delta1_r: delta1_r,
            delta2_l: delta2_l,
            delta2_r: delta2_r,
            kai: kai,
            ht: ht,
            _hpa: PhantomData,
        })
    }

    pub fn prove<R: Rng>(
        values: (
            &[IP::LeftMessage],
            &[IP::RightMessage],
            &[IP::LeftMessage],
            &[IP::LeftMessage],
        ),
        srs: &HPASRS<IP, LMC, RMC, IPC, D>,
        ck_message: (&[LMC::Message], &[RMC::Message]),
        // com: (&IP::Output, &IP::Output, &IP::Output),
        witness: (
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
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
            (values.0, values.1, values.2, values.3),
            srs,
            (ck_message.0, ck_message.1),
            witness,
            gm,
            rng,
        )?;
        Ok(proof)
    }

    pub fn verify<R: Rng>(
        srs: &mut HPASRS<IP, LMC, RMC, IPC, D>, //
        // ck: (&[LMC::Key], &[RMC::Key]),
        ck_message: (&[LMC::Message], &[RMC::Message]),
        com: (&IP::Output, &IP::Output, &IP::Output, &IP::Output, &IP::Output, &IP::Output, &IP::Output), // com ( c, x, y, d1, d2, d3, d4 )
        proof: &mut HPAProof<IP, LMC, RMC, IPC, D>,
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R
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

        let gamma1 = ck_message.0.clone();
        let gamma2 = ck_message.1.clone();

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
        let mut y_prime = com.2.clone();
        let mut d1_prime = com.3.clone();
        let mut d2_prime = com.4.clone();
        let mut d3_prime = com.5.clone();
        let mut d4_prime = com.6.clone();
        let mut result = false;

        if round > 0 {
            for i in 0..round {
                // println!("check");
                // Verifier's work in reduce
                let last_commitment = proof.r_commitment_steps.pop().unwrap();
                let last_transcript = transcript.pop().unwrap();

                let c_l = last_commitment.0.0.clone();
                let c_r = last_commitment.0.1.clone();
                let x_l = last_commitment.0.2.clone();
                let x_r = last_commitment.0.3.clone();
                let d1_l = last_commitment.1.0.clone();
                let d1_r = last_commitment.1.1.clone();
                let d3_l = last_commitment.1.2.clone();
                let d3_r = last_commitment.1.3.clone();

                assert!(c_prime == c_l + c_r);
                assert!(x_prime == x_l + x_r);
                assert!(d1_prime == d1_l + d1_r);
                assert!(d3_prime == d3_l + d3_r);

                let alpha = last_transcript.0;
                let alpha_inv = last_transcript.1;
                let beta = last_transcript.2;
                let beta_inv = last_transcript.3;
                let gm_inv = last_transcript.4;
                let alpha_sqr = alpha * alpha;
                let alpha_sqr_beta_inv = alpha_sqr * beta_inv;
                let gm_inv_beta_inv = gm_inv * beta_inv;
                let alpha_gm_inv = alpha * gm_inv;
                let alpha_beta_inv = alpha * beta_inv;
                let alpha_inv_beta = alpha_inv * beta;

                let last_c_x = proof.c_x.pop().unwrap();
                let last_x_plus = proof.x_plus.pop().unwrap();
                let last_x_minus = proof.x_minus.pop().unwrap();
                let last_y_plus = proof.y_plus.pop().unwrap();
                let last_y_minus = proof.y_minus.pop().unwrap();

                c_prime = mul_helper(&last_commitment.0.0, &alpha_sqr) + last_commitment.0.1 + mul_helper(&last_commitment.1.0, &alpha_sqr_beta_inv) 
                    + mul_helper(&last_commitment.1.1, &beta_inv) + mul_helper(&last_c_x, &alpha);
                x_prime = mul_helper(&last_commitment.0.2, &alpha_sqr) + mul_helper(&last_commitment.0.3, &gm_inv) + mul_helper(&last_commitment.1.2, &alpha_sqr_beta_inv)
                    + mul_helper(&last_commitment.1.3, &gm_inv_beta_inv) + mul_helper(&last_x_plus, &alpha) + mul_helper(&last_x_minus, &alpha_gm_inv);
                y_prime = y_prime + srs.kai.pop().unwrap() + mul_helper(&d2_prime, &beta) + mul_helper(&d4_prime, &beta_inv) + mul_helper(&last_y_minus, &alpha) + mul_helper(&last_y_plus, &alpha_inv);
                d1_prime = mul_helper(&last_commitment.2.0, &alpha) + last_commitment.2.1;
                d2_prime = mul_helper(&last_commitment.2.2, &alpha) + last_commitment.2.3 + mul_helper(&srs.delta2_l.pop().unwrap(), &alpha_beta_inv) + mul_helper(&srs.delta2_r.pop().unwrap(), &beta_inv);
                d3_prime = mul_helper(&last_commitment.3.0, &alpha) + mul_helper(&last_commitment.3.1, &gm_inv);
                d4_prime = mul_helper(&last_commitment.3.2, &alpha_inv) + last_commitment.3.3 + mul_helper(&srs.delta1_l.pop().unwrap(), &alpha_inv_beta) + mul_helper(&srs.delta1_r.pop().unwrap(), &beta);

                // Scalar product
                if i == round - 1 {
                    let mut e1 = proof.e1.clone();
                    let mut e2 = proof.e2.clone();

                    let ch_d = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                    let ch_d_inv = ch_d.inverse().unwrap();

                    // check pairing equation
                    let kai_scalar =
                        IP::inner_product(&(gamma1[..1].to_vec()), &(gamma2[..1].to_vec()))?;

                    let temp1 = mul_helper(&(gamma1[0]), &ch_d);
                    e1[0] = e1[0].clone() + temp1;
                    e2[0] = e2[0].clone() + mul_helper(&(gamma2[0]), &(ch_d_inv));

                    let left = IP::inner_product(&e1, &e2)?;
                    
                    let one = <LMC as DoublyHomomorphicCommitment>::Scalar::one();
                    let zero: <LMC as DoublyHomomorphicCommitment>::Scalar =
                        <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
                    let minus_one = zero - one;

                    let mut ch_c_vec = Vec::new();
                    ch_c_vec.push(ch_c.clone()); // ch_c_vec = {c^1, c^2, c^3, ..., c^11}
                    for i in 1..11{
                        ch_c_vec.push(ch_c_vec[i-1] * ch_c);
                    }

                    let mut temp2 = c_prime.clone() + mul_helper(&x_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&y_prime, &ch_c_vec[8]);
                    temp2 = mul_helper(&temp2, &ch_c_vec[1]);
                    let mut temp3 = d1_prime.clone() + mul_helper(&d3_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&d4_prime, &ch_c_vec[8]);
                    temp3 = mul_helper(&temp3, &(ch_c * ch_d_inv));
                    let mut temp4 = proof.q1.clone() + mul_helper(&proof.q2, &ch_c_vec[2]) + mul_helper(&proof.q3, &ch_c_vec[5]) + mul_helper(&proof.p3, &ch_c_vec[8]);
                    temp4 = mul_helper(&temp4, &ch_c);
                    let temp5 = ch_d_inv * proof.r1 + ch_d * proof.r2 + proof.r3;
                    let mut temp6 = mul_helper(&srs.ht, &temp5);
                    temp6 = mul_helper(&temp6, &minus_one);

                    let right = kai_scalar + temp2 + mul_helper(&d2_prime, &(ch_c * ch_d)) + temp3 + temp4 + mul_helper(&proof.p1, &ch_d_inv)
                        + mul_helper(&proof.p2, &ch_d) + proof.r.clone() + temp6;
                    
                    result = left == right;
                }
            }
            Ok(result)
        } else {
            
            let mut e1 = proof.e1.clone();
            let mut e2 = proof.e2.clone();

            let ch_d = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
            let ch_d_inv = ch_d.inverse().unwrap();

            // check pairing equation
            let kai_scalar = IP::inner_product(&(gamma1[..1].to_vec()), &(gamma2[..1].to_vec()))?;

            e1[0] = e1[0].clone() + mul_helper(&(gamma1[0]), &ch_d);
            e2[0] = e2[0].clone() + mul_helper(&(gamma2[0]), &(ch_d_inv));
            let left = IP::inner_product(&e1, &e2)?;

            let one = <LMC as DoublyHomomorphicCommitment>::Scalar::one();
            let zero: <LMC as DoublyHomomorphicCommitment>::Scalar =
                <LMC as DoublyHomomorphicCommitment>::Scalar::zero();
            let minus_one = zero - one;

            let mut ch_c_vec = Vec::new();
            ch_c_vec.push(ch_c.clone()); // ch_c_vec = {c^1, c^2, c^3, ..., c^11}
            for i in 1..11{
                ch_c_vec.push(ch_c_vec[i-1] * ch_c);
            }

            let mut temp2 = c_prime.clone() + mul_helper(&x_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&y_prime, &ch_c_vec[8]);
            temp2 = mul_helper(&temp2, &ch_c_vec[1]);
            let mut temp3 = d1_prime.clone() + mul_helper(&d3_prime, &(ch_c_vec[2] + ch_c_vec[5])) + mul_helper(&d4_prime, &ch_c_vec[8]);
            temp3 = mul_helper(&temp3, &(ch_c * ch_d_inv));
            let mut temp4 = proof.q1.clone() + mul_helper(&proof.q2, &ch_c_vec[2]) + mul_helper(&proof.q3, &ch_c_vec[5]) + mul_helper(&proof.p3, &ch_c_vec[8]);
            temp4 = mul_helper(&temp4, &ch_c);
            let temp5 = ch_d_inv * proof.r1 + ch_d * proof.r2 + proof.r3;
            let mut temp6 = mul_helper(&srs.ht, &temp5);
            temp6 = mul_helper(&temp6, &minus_one);

            let right = kai_scalar + temp2 + mul_helper(&d2_prime, &(ch_c * ch_d)) + temp3 + temp4 + mul_helper(&proof.p1, &ch_d_inv)
                + mul_helper(&proof.p2, &ch_d) + proof.r.clone() + temp6;
            
            result = left == right;
            Ok(result)
        }
    }

    pub fn prove_with_aux<R: Rng>(
        values: (
            &[IP::LeftMessage],
            &[IP::RightMessage],
            &[IP::LeftMessage],
            &[IP::LeftMessage],
        ),
        srs: &HPASRS<IP, LMC, RMC, IPC, D>,
        ck_message: (&[LMC::Message], &[RMC::Message]),
        witness: (
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<(HPAProof<IP, LMC, RMC, IPC, D>, HPAAux<IP, LMC, RMC, IPC, D>), Error> {
        let (v1, v2, w_vec, k_vec) = values;
        // let (gamma1, gamma2) = ck;
        let (gamma1_message, gamma2_message) = ck_message;
        Self::_prove(
            &(v1.to_vec(), v2.to_vec(), w_vec.to_vec(), k_vec.to_vec()),
            srs,
            (gamma1_message.to_vec(), gamma2_message.to_vec()),
            witness,
            gm,
            rng,
        )
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove<R: Rng>(
        values: &(
            Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>,
            Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>,
        ),
        srs: &HPASRS<IP, LMC, RMC, IPC, D>,
        ck_message: (Vec<LMC::Message>, Vec<RMC::Message>),
        witness: (
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
            &<LMC as DoublyHomomorphicCommitment>::Scalar,
        ),
        gm: &<LMC as DoublyHomomorphicCommitment>::Scalar,
        rng: &mut R,
    ) -> Result<(HPAProof<IP, LMC, RMC, IPC, D>, HPAAux<IP, LMC, RMC, IPC, D>), Error> {
        let (mut v1, mut v2, mut w_vec, mut k_vec) = values.clone();
        // let (gamma1, gamma2) = ck.clone();
        let (mut gamma1_message, mut gamma2_message) = ck_message.clone();
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        assert!(v1.len().is_power_of_two());

        let ht = srs.ht.clone();
        let mut r_c = witness.0.clone();
        let mut r_x = witness.1.clone();
        let mut r_y = witness.2.clone();
        let mut r_d1 = witness.3.clone();
        let mut r_d2 = witness.4.clone();
        let mut r_d3 = witness.5.clone();
        let mut r_d4 = witness.6.clone();
        let mut c_x_vec = Vec::new();
        let mut x_plus_vec = Vec::new();
        let mut x_minus_vec = Vec::new();
        let mut y_plus_vec = Vec::new();
        let mut y_minus_vec = Vec::new();

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
                let r_xl = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d1l = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d3l = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d1l_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d1r_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d2l_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d2r_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d3l_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d3r_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d4l_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_d4r_prime = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

                let r_cr = r_c - r_cl;
                let r_xr = r_x - r_xl;
                let r_d1r = r_d1 - r_d1l;
                let r_d3r = r_d3 - r_d3l;

                let v1_l = &v1[..split];
                let v1_r = &v1[split..];
                let gamma1_prime = &gamma1_message[..split];
                // let gamma1_r = &gamma1_message[split..];

                let v2_l = &v2[..split];
                let v2_r = &v2[split..];
                let gamma2_prime = &gamma2_message[..split];
                let gamma2_r = &gamma2_message[split..];

                let w_vec_l = &w_vec[..split];
                let w_vec_r = &w_vec[split..];

                let k_vec_l = &k_vec[..split];
                let k_vec_r = &k_vec[split..];

                let cl = start_timer!(|| "Compute D");
                let c_l = IP::inner_product(&v1_l, &v2_l).unwrap() + mul_helper(&ht, &r_cl);
                let c_r = IP::inner_product(&v1_r, &v2_r).unwrap() + mul_helper(&ht, &r_cr);
                let x_l = IP::inner_product(&w_vec_l, &v2_l).unwrap() + mul_helper(&ht, &r_xl);
                let x_r = IP::inner_product(&w_vec_r, &v2_r).unwrap() + mul_helper(&ht, &r_xr);
                let d1_l = IP::inner_product(&v1_l, &gamma2_prime)? + mul_helper(&ht, &r_d1l);
                let d1_r = IP::inner_product(&v1_r, &gamma2_r)? + mul_helper(&ht, &r_d1r);
                let d3_l = IP::inner_product(&w_vec_l, &gamma2_prime)? + mul_helper(&ht, &r_d3l);
                let d3_r = IP::inner_product(&w_vec_r, &gamma2_r)? + mul_helper(&ht, &r_d3r);
                let d1_l_prime = IP::inner_product(&v1_l, &gamma2_prime).unwrap()
                    + mul_helper(&ht, &r_d1l_prime);
                let d1_r_prime = IP::inner_product(&v1_r, &gamma2_prime).unwrap()
                    + mul_helper(&ht, &r_d1r_prime);
                let d2_l_prime = IP::inner_product(&gamma1_prime, &v2_l).unwrap()
                    + mul_helper(&ht, &r_d2l_prime);
                let d2_r_prime = IP::inner_product(&gamma1_prime, &v2_r).unwrap()
                    + mul_helper(&ht, &r_d2r_prime);
                let d3_l_prime = IP::inner_product(&w_vec_l, &gamma2_prime).unwrap()
                    + mul_helper(&ht, &r_d3l_prime);
                let d3_r_prime = IP::inner_product(&w_vec_r, &gamma2_prime).unwrap()
                    + mul_helper(&ht, &r_d3r_prime);
                let d4_l_prime = IP::inner_product(&k_vec_l, &gamma2_prime).unwrap()
                    + mul_helper(&ht, &r_d4l_prime);
                let d4_r_prime = IP::inner_product(&k_vec_r, &gamma2_prime).unwrap()
                    + mul_helper(&ht, &r_d4r_prime);

                // Fiat-Shamir challenge
                let mut counter_nonce: usize = 0;
                //  let default_transcript = (Default::default(),Default::default(),Default::default(),Default::default());
                //  let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (beta, beta_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input.extend_from_slice(&to_bytes![
                        c_l, c_r, x_l, x_r, d1_l, d1_r, d3_l, d3_r, d1_l_prime, d1_r_prime,
                        d2_l_prime, d2_r_prime, d3_l_prime, d3_r_prime, d4_l_prime, d4_r_prime
                    ]?);
                    let beta: LMC::Scalar = u128::from_be_bytes(
                        D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                    )
                    .into();

                    //  let list_beta_input = List(hash_input.clone());
                    //  println!("prove - hash beta - {}", list_beta_input);

                    if let Some(beta_inv) = beta.inverse() {
                        // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                        // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                        break 'challenge (beta_inv, beta);
                    }
                    counter_nonce += 1;
                };

                end_timer!(cl);
                let gamma1_message_temp = gamma1_message.clone();
                let gamma1_beta = cfg_iter!(gamma1_message_temp)
                    .map(|gamma| mul_helper(gamma, &beta))
                    .collect::<Vec<LMC::Message>>();
                let gamma2_message_temp = gamma2_message.clone();
                let gamma2_beta_inv = cfg_iter!(gamma2_message_temp)
                    .map(|gamma| mul_helper(gamma, &beta_inv))
                    .collect::<Vec<RMC::Message>>();

                // println!("v1 len : {}", v1.len());

                for i in 0..v1.len() {
                    k_vec[i] = k_vec[i].clone() + gamma1_beta[i].clone();
                    v2[i] = v2[i].clone() + gamma2_beta_inv[i].clone();
                }

                r_y = r_y.clone() + beta * r_d2.clone() + beta_inv * r_d4.clone();

                let r_c_x = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_x_plus = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_x_minus = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_y_plus = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
                let r_y_minus = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

                // compute C and message rescale

                let cr = start_timer!(|| "Compute C");

                let v2_l = v2[..split].to_vec();
                let v2_r = v2[split..].to_vec();
                let k_vec_l = k_vec[..split].to_vec();
                let k_vec_r = k_vec[split..].to_vec();

                let c_x = IP::inner_product(&v1_l, &v2_r).unwrap()
                    + IP::inner_product(&v1_r, &v2_l).unwrap()
                    + mul_helper(&ht, &r_c_x);
                let x_plus =
                    IP::inner_product(&w_vec_l, &v2_r).unwrap() + mul_helper(&ht, &r_x_plus);
                let x_minus =
                    IP::inner_product(&w_vec_r, &v2_l).unwrap() + mul_helper(&ht, &r_x_minus);
                let y_plus =
                    IP::inner_product(&k_vec_l, &v2_r).unwrap() + mul_helper(&ht, &r_y_plus);
                let y_minus =
                    IP::inner_product(&k_vec_r, &v2_l).unwrap() + mul_helper(&ht, &r_y_minus);

                end_timer!(cr);

                // Second Fiat-Shamir challenge
                let (alpha, alpha_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input
                        .extend_from_slice(&to_bytes![c_x, x_plus, x_minus, y_plus, y_minus]?);
                    let alpha: LMC::Scalar = u128::from_be_bytes(
                        D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                    )
                    .into();

                    // let list_alpha_input = List(hash_input.clone());
                    //  println!("prove - hash alpha - {}", list_alpha_input);

                    if let Some(alpha_inv) = alpha.inverse() {
                        // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                        // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                        break 'challenge (alpha_inv, alpha);
                    }
                    counter_nonce += 1;
                };

                // Set up values for next step of recursion
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

                let rescale_k_vec = start_timer!(|| "Rescale K");
                k_vec = cfg_iter!(k_vec_l)
                    .map(|b| mul_helper(b, &alpha_inv))
                    .zip(k_vec_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<LMC::Message>>();
                end_timer!(rescale_k_vec);

                let alpha_sqr = alpha * alpha;
                let alpha_sqr_beta_inv = alpha_sqr * beta_inv;
                let gm_beta_inv = gm_inv * beta_inv;
                let gm_alpha = gm_inv * alpha;

                r_c = alpha_sqr * r_cl
                    + r_cr
                    + alpha_sqr_beta_inv * r_d1l
                    + beta_inv * r_d1r
                    + alpha * r_c_x;
                r_x = alpha_sqr * r_xl
                    + alpha_sqr_beta_inv * r_d3l
                    + alpha * r_x_plus
                    + gm_inv * r_xr
                    + gm_beta_inv * r_d3r
                    + gm_alpha * r_x_minus;
                r_y = r_y + alpha * r_y_minus + alpha_inv * r_y_plus;
                r_d1 = alpha * r_d1l_prime + r_d1r_prime;
                r_d2 = alpha * r_d2l_prime + r_d2r_prime;
                r_d3 = alpha * r_d3l_prime + gm_inv * r_d3r_prime;
                r_d4 = alpha_inv * r_d4l_prime + r_d4r_prime;

                // println!("v1 len : {}", v1.len());

                // gamma1 = gamma1_prime.to_vec();
                // gamma2 = gamma2_prime.to_vec();
                gamma1_message = gamma1_message[..split].to_vec();
                gamma2_message = gamma2_message[..split].to_vec();

                let com1 = (c_l, c_r, x_l, x_r);
                let com2 = (d1_l, d1_r, d3_l, d3_r);
                let com3 = (d1_l_prime, d1_r_prime, d2_l_prime, d2_r_prime);
                let com4 = (d3_l_prime, d3_r_prime, d4_l_prime, d4_r_prime);

                r_commitment_steps.push((com1, com2, com3, com4));
                r_transcript.push((alpha, alpha_inv, beta, beta_inv, gm_inv));
                c_x_vec.push(c_x);
                x_plus_vec.push(x_plus);
                x_minus_vec.push(x_minus);
                y_plus_vec.push(y_plus);
                y_minus_vec.push(y_minus);
                end_timer!(recurse);
            }
        };

        let r_p1 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_p2 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_p3 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_q1 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_q2 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_q3 = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);
        let r_r = <LMC as DoublyHomomorphicCommitment>::Scalar::rand(rng);

        let mut d1 = Vec::new();
        let mut d2 = Vec::new();
        d1.push(<IP::LeftMessage>::rand(rng));
        d2.push(<IP::RightMessage>::rand(rng));

        let p1 = IP::inner_product(&d1, &gamma2_message).unwrap() + mul_helper(&ht, &r_p1);
        let p2 = IP::inner_product(&gamma1_message, &d2).unwrap() + mul_helper(&ht, &r_p2);
        let p3 = IP::inner_product(&k_vec, &d2).unwrap() + mul_helper(&ht, &r_p3);
        let q1 = IP::inner_product(&d1, &v2).unwrap()
            + IP::inner_product(&v1, &d2).unwrap()
            + mul_helper(&ht, &r_q1);
        let q2 = IP::inner_product(&v1, &d2).unwrap() + mul_helper(&ht, &r_q2);
        let q3 = IP::inner_product(&w_vec, &d2).unwrap() + mul_helper(&ht, &r_q3);
        let r = IP::inner_product(&d1, &d2).unwrap() + mul_helper(&ht, &r_r);

        let ch_c = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&to_bytes![p1, p2, p3, q1, q2, q3, r]?);
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
        let ch_c_8 = ch_c_7 * ch_c;
        let ch_c_9 = ch_c_8 * ch_c;
        let ch_c_10 = ch_c_9 * ch_c;
        let ch_c_11 = ch_c_10 * ch_c;

        let mut e1 = Vec::new();
        let mut e2 = Vec::new();
        e1.push(d1[0].clone() + mul_helper(&v1[0], &(ch_c + ch_c_4)) + mul_helper(&w_vec[0], &ch_c_7) + mul_helper(&k_vec[0], &ch_c_10));
        e2.push(d2[0].clone() + mul_helper(&v2[0], &ch_c));
        let r1 = r_p1 + ch_c * r_d1 + (ch_c_4 + ch_c_7) * r_d3 + ch_c_10 * r_d4;
        let r2 = r_p2 + ch_c * r_d2;
        let r3 = r_r + ch_c_2 * r_c + (ch_c_5 + ch_c_8) * r_x + ch_c_11 * r_y + ch_c * r_q1 + ch_c_4 * r_q2 + ch_c_7 * r_q3 + ch_c_10 * r_p3;

        // r_transcript.reverse();
        r_commitment_steps.reverse();
        c_x_vec.reverse();
        x_plus_vec.reverse();
        x_minus_vec.reverse();
        y_plus_vec.reverse();
        y_minus_vec.reverse();
        Ok((
            HPAProof {
                r_commitment_steps,
                c_x:c_x_vec, 
                x_plus:x_plus_vec, 
                x_minus:x_minus_vec, 
                y_plus:y_plus_vec, 
                y_minus:y_minus_vec,

                e1, e2,
                p1, p2, p3,
                q1, q2, q3,
                r,
                r1,
                r2,
                r3,
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
            Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar)>,
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
            Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar)>,
            LMC::Scalar,
        ),
        Error,
    > {

        let mut r_transcript = Vec::new();
        let mut c_x = proof.c_x.clone();
        let mut x_plus = proof.x_plus.clone();
        let mut x_minus = proof.x_minus.clone();
        let mut y_plus = proof.y_plus.clone();
        let mut y_minus = proof.y_minus.clone();

        for (com_1, com_2, com_3, com_4) in proof.r_commitment_steps.iter().rev() {
            // First Fiat-Shamir challenge
            let mut counter_nonce: usize = 0;
            // let default_transcript = (Default::default(), Default::default(),Default::default(),Default::default());
            // let transcript = r_transcript.last().unwrap_or(&default_transcript);
            let (beta, beta_inv) = 'challenge: loop {
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                hash_input.extend_from_slice(&to_bytes![com_1.0, com_1.1, com_1.2, com_1.3, com_2.0, com_2.1, com_2.2, com_2.3,
                    com_3.0, com_3.1, com_3.2, com_3.3, com_4.0, com_4.1, com_4.2, com_4.3]?);
                let beta: LMC::Scalar = u128::from_be_bytes(
                    D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                )
                .into();

                // let list_beta_input = List(hash_input.clone());
                // println!("prove - hash beta - {}", list_beta_input);

                if let Some(beta_inv) = beta.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (beta_inv, beta);
                }
                counter_nonce += 1;
            };

            // Second Fiat-Shamir challenge
            let (alpha, alpha_inv) = 'challenge: loop {
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                hash_input.extend_from_slice(&to_bytes![c_x.pop().unwrap(), x_plus.pop().unwrap(), x_minus.pop().unwrap(), y_plus.pop().unwrap(), y_minus.pop().unwrap()]?);        ////////////////
                let alpha: LMC::Scalar = u128::from_be_bytes(
                    D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                )
                .into();

                // let list_alpha_input = List(hash_input.clone());
                // println!("prove - hash alpha - {}", list_alpha_input);

                if let Some(alpha_inv) = alpha.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (alpha_inv, alpha);
                }
                counter_nonce += 1;
            };


            // TODO: compute gm, and push it also
            let mut gm_inv = gm.inverse().unwrap();
            let len = c_x.len();
            for _ in 0..len {
                gm_inv = gm_inv * gm_inv;
            }

            r_transcript.push((alpha, alpha_inv, beta, beta_inv, gm_inv));
        }
        r_transcript.reverse();

        let ch_c = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&to_bytes![proof.p1, proof.p2, proof.p3, proof.q1, proof.q2, proof.q3, proof.r]?);
            let ch_c: LMC::Scalar =
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
        y: &IP::Output,
        y_: &IP::Output,
        d1: &IP::Output,
        d1_: &IP::Output,
        d2: &IP::Output,
        d2_: &IP::Output,
        d3: &IP::Output,
        d3_: &IP::Output,
        d4: &IP::Output,
        d4_: &IP::Output,
        z_c: &IP::Output,
        z_x: &IP::Output,
        z_y: &IP::Output,
        delta: &LMC::Scalar,
    ) -> Result<
        (
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            IP::Output,
            // c, x, y, d1, d2, d3, d4
        ),
        Error,
    > {

        let delta = delta.clone();
        let delta_sqr = delta * delta;
        let bat_c = mul_helper(c, &delta_sqr) + mul_helper(z_c, &delta) + c_.clone();
        let bat_x = mul_helper(x, &delta_sqr) + mul_helper(z_x, &delta) + x_.clone();
        let bat_y = mul_helper(y, &delta_sqr) + mul_helper(z_y, &delta) + y_.clone();
        let bat_d1 = mul_helper(d1, &delta) + d1_.clone();
        let bat_d2 = mul_helper(d2, &delta) + d2_.clone();
        let bat_d3 = mul_helper(d3, &delta) + d3_.clone();
        let bat_d4 = mul_helper(d4, &delta) + d4_.clone();

        Ok((bat_c, bat_x, bat_y, bat_d1, bat_d2, bat_d3, bat_d4))
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
            c_x: self.c_x.clone(),
            x_plus: self.x_plus.clone(),
            x_minus: self.x_minus.clone(),
            y_plus: self.y_plus.clone(),
            y_minus: self.y_minus.clone(),
            e1: self.e1.clone(),
            e2: self.e2.clone(),
            p1: self.p1.clone(),
            p2: self.p2.clone(),
            p3: self.p3.clone(),
            q1: self.q1.clone(),
            q2: self.q2.clone(),
            q3: self.q3.clone(),
            r: self.r.clone(),
            r1: self.r1.clone(),
            r2: self.r2.clone(),
            r3: self.r3.clone(),
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
