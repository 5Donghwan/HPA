extern crate ark_ff;
use self::ark_ff::{to_bytes, Field};
extern crate ark_serialize;
use self::ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
extern crate ark_std;
use self::ark_std::rand::Rng;
use self::ark_std::{end_timer, start_timer};
extern  crate digest;
use self::digest::Digest;
use std::{convert::TryInto, marker::PhantomData, ops::MulAssign};

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
            if count != 0 { write!(f, ", ")?; }
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
        (IP::Output, IP::Output, IPC::Message), 
        (IP::Output, IP::Output, IPC::Message),
    )>,
    pub(crate) e1: Vec<IP::LeftMessage>,
    pub(crate) e2: Vec<IP::RightMessage>,
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
    // IPC::Message: AddAssign<LMC::Output>,
    // // IPC::Message: AddAssign<IPC::Message>,
    // IPC::Message: AddAssign<RMC::Output>,
    // IPC::Message: AddAssign<LMC::Output>,
    // RMC::Output: AddAssign<LMC::Output>,
    // IPC::Output: AddAssign<LMC::Output>,
{
    pub fn setup<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<(Vec<LMC::Key>, Vec<RMC::Key>), Error> {
        Ok((
            LMC::setup(rng, size)?,
            RMC::setup(rng, size)?,
        ))
    }

    pub fn precompute(
        ck_message: (&[LMC::Message], &[RMC::Message]),
    ) -> Result<HPASRS<IP,LMC,RMC,IPC,D>, Error> {
        // loop : until ck.len() >= 1
        let (mut gamma1, mut gamma2) = ck_message.clone();
        // let mut i = ck_message.0.len();
        let mut delta1_l = Vec::new();
        let mut delta1_r = Vec::new();
        let mut delta2_r = Vec::new();
        let mut kai = Vec::new();
        let mut split = ck_message.0.len()/2;
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
            
            split = split/2;
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

        Ok(HPASRS { delta1_l: delta1_l, delta1_r: delta1_r, delta2_l: delta2_l, delta2_r: delta2_r, kai: kai, _hpa: PhantomData })
    }

    pub fn prove(
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        // ck: (&[RMC::Key], &[LMC::Key]),
        ck_message: (&[LMC::Message], &[RMC::Message]),
        com: (&IP::Output, &IP::Output, &IP::Output),
    ) -> Result<HPAProof<IP, LMC, RMC, IPC, D>, Error> {
        if IP::inner_product(values.0, values.1)? != com.2.clone() {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }
        if values.0.len().count_ones() != 1 {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                values.0.len(),
                values.1.len(),
            )));
        }
        // if !(LMC::verify(ck.1, values.0, com.0)?
        //     && RMC::verify(ck.0, values.1, com.1)?
        //     )
        // {
        //     return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        // }

        // TODO : compare ck and ck_message

        let (proof, _) =
            Self::prove_with_aux((values.0, values.1), 
            // (ck.0, ck.1), 
            (ck_message.0, ck_message.1))?;
        Ok(proof)
    }

    pub fn verify(
        srs: &mut HPASRS<IP, LMC, RMC, IPC, D>,  //
        // ck: (&[LMC::Key], &[RMC::Key]),
        ck_message: (&[LMC::Message], &[RMC::Message]),
        com: (&IP::Output, &IP::Output, &IP::Output), // com ( d1, d2, c )
        proof: &mut HPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<bool, Error> {
        if ck_message.0.len().count_ones() != 1 || ck_message.0.len() != ck_message.1.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                ck_message.0.len(),
                ck_message.1.len(),
            )));
        }
        // Calculate transcript
        let mut transcript = Self::_compute_recursive_challenges(
            proof,
        )?;

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

        let mut c_prime = com.2.clone();
        let mut d1_prime = com.0.clone();
        let mut d2_prime = com.1.clone();
        let mut result = false;

        if round > 0 {
            for i in 0..round {
                // println!("check");
                // Verifier's work in reduce
                let last_commitment = proof.r_commitment_steps.pop().unwrap();
                let last_transcript = transcript.pop().unwrap();
                let temp2 = mul_helper(&d1_prime, &(last_transcript.3));
                let temp = mul_helper(&d2_prime, &(last_transcript.2));
                let temp = temp + temp2;//add_helper(&temp, &temp2);
                let last_kai = srs.kai.pop().unwrap();
                let temp = last_kai + temp;//add_helper(&last_kai, &temp);

                c_prime = c_prime + temp + mul_helper(&(last_commitment.0.2), &(last_transcript.0)) + mul_helper(&(last_commitment.1.2), &(last_transcript.1)); 
                let temp = mul_helper(&(last_commitment.0.0.clone()), &(last_transcript.0)) + last_commitment.1.0;
                d1_prime = mul_helper(&(srs.delta1_l.pop().unwrap()), &(last_transcript.0 * last_transcript.2)) + mul_helper(&(srs.delta1_r.pop().unwrap()), &(last_transcript.2));
                d1_prime = d1_prime + temp; //add_helper(&d1_prime, &temp);
                let temp2 =  mul_helper(&(last_commitment.0.1), &(last_transcript.1)) + last_commitment.1.1;
                d2_prime = mul_helper(&(srs.delta2_l.pop().unwrap()), &(last_transcript.1 * last_transcript.3)) + mul_helper(&(srs.delta2_r.pop().unwrap()), &(last_transcript.3));
                d2_prime = d2_prime + temp2;

                // Scalar product
                if i == round-1 {
                    let mut e1 = proof.e1.clone();
                    let mut e2 = proof.e2.clone();


                    // Fiat-Schamir challenge
                    let (d, d_inv) = 'challenge: loop {
                        let mut hash_input = Vec::new();
                        //TODO: Should use CanonicalSerialize instead of ToBytes
                        hash_input.extend_from_slice(&to_bytes![
                            e1[0],e2[0]
                        ]?);
                        let d:LMC::Scalar = u128::from_be_bytes(
                            D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                        )
                        .into();
                        if let Some(d_inv) = d.inverse() {
                            // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                            // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                            break 'challenge (d_inv, d);
                        }
                    };

                    // check pairing equation
                    let kai_scalar = IP::inner_product(&(gamma1[..1].to_vec()), &(gamma2[..1].to_vec()))?;


                    let temp3 = mul_helper(&(gamma1[0]), &d);
                    e1[0] = e1[0].clone() + temp3;
                    e2[0] = e2[0].clone() + mul_helper(&(gamma2[0]), &(d_inv));

                    let left = IP::inner_product(&e1, &e2)?;
                    let temp1 = c_prime.clone() + kai_scalar;
                    let temp2 = mul_helper(&d2_prime, &d);
                    let temp4 = mul_helper(&d1_prime, &d_inv);
                    let temp5 = temp2 + temp4;//add_helper(&temp2, &temp4);
                    let right = temp1 + temp5;//add_helper(&temp1,&temp5);
                    // let left = IP::inner_product(&e1, &e2)?;
                    // let temp1 = c_prime + kai_scalar;
                    // let temp2 = mul_helper(&d2_prime, &d) + mul_helper(&d1_prime, &d_inv);
                    // let right = temp1 +temp2;
                    result = left == right;
                }

                // c = c_prime;
                // d1 = d1_prime;
                // d2 = d2_prime;
            }
            Ok(result)
        }
        else{
            let mut e1 = proof.e1.clone();
            let mut e2 = proof.e2.clone();


            // Fiat-Schamir challenge
            let (d, d_inv) = 'challenge: loop {
            let mut hash_input = Vec::new();
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![
                e1[0],e2[0]
                ]?);
            let d:LMC::Scalar = u128::from_be_bytes(
                D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                )
                .into();
                if let Some(d_inv) = d.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (d_inv, d);
                }
            };

            // check pairing equation
            let kai_scalar = IP::inner_product(&(gamma1[..1].to_vec()), &(gamma2[..1].to_vec()))?;


            let temp3 = mul_helper(&(gamma1[0]), &d);
            e1[0] = e1[0].clone() + temp3;
            e2[0] = e2[0].clone() + mul_helper(&(gamma2[0]), &(d_inv));
            let left = IP::inner_product(&e1, &e2)?;
            let temp1 = c_prime.clone() + kai_scalar;
            let temp2 = mul_helper(&d2_prime, &d);
            let temp4 = mul_helper(&d1_prime, &d_inv);
            let temp5 = temp2 + temp4;//add_helper(&temp2, &temp4);
            let right = temp1 + temp5;//add_helper(&temp1,&temp5);
            result = left == right;
            Ok(result)
        }
    }

    pub fn prove_with_aux(
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        // ck: (&[RMC::Key], &[LMC::Key]),
        ck_message: (&[LMC::Message], &[RMC::Message]),
    ) -> Result<
        (
            HPAProof<IP, LMC, RMC, IPC, D>,
            HPAAux<IP, LMC, RMC, IPC, D>,
        ),
        Error,
    > {
        let (v1, v2) = values;
        // let (gamma1, gamma2) = ck;
        let (gamma1_message, gamma2_message) = ck_message;
        Self::_prove(
            &(v1.to_vec(), v2.to_vec()),
            // (gamma1.to_vec(), gamma2.to_vec()),
            (gamma1_message.to_vec(), gamma2_message.to_vec()),
        )
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove(
        values: &(Vec<IP::LeftMessage>, Vec<IP::RightMessage>),
        // ck: (Vec<RMC::Key>, Vec<LMC::Key>),
        ck_message: (Vec<LMC::Message>, Vec<RMC::Message>),
    ) -> Result<
        (
            HPAProof<IP, LMC, RMC, IPC, D>,
            HPAAux<IP, LMC, RMC, IPC, D>,
        ),
        Error,
    > {
        let ( mut v1,   mut v2) = values.clone() ;
        // let (gamma1, gamma2) = ck.clone();
        let (mut gamma1_message, mut gamma2_message) = ck_message.clone();
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        assert!(v1.len().is_power_of_two());

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

                let v1_l = &v1[..split];
                let v1_r = &v1[split..];
                let gamma1_prime = &gamma1_message[..split];
                
                let v2_l = &v2[..split];
                let v2_r = &v2[split..];
                let gamma2_prime = &gamma2_message[..split];

                let cl = start_timer!(|| "Compute D");
                let d1_l = IP::inner_product(&v1_l, &gamma2_prime)?;
                let d1_r = IP::inner_product(&v1_r, &gamma2_prime)?;
                let d2_l = IP::inner_product(&gamma1_prime, &v2_l)?;
                let d2_r = IP::inner_product(&gamma1_prime, &v2_r)?;
                
                 // Fiat-Shamir challenge
                 let mut counter_nonce: usize = 0;
                //  let default_transcript = (Default::default(),Default::default(),Default::default(),Default::default());
                //  let transcript = r_transcript.last().unwrap_or(&default_transcript);
                 let (beta, beta_inv) = 'challenge: loop {
                     let mut hash_input = Vec::new();
                     hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                     //TODO: Should use CanonicalSerialize instead of ToBytes
                     hash_input.extend_from_slice(&to_bytes![
                        d1_l, d1_r, d2_l, d2_r
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
                    .map(|gamma| mul_helper(gamma, &beta)).collect::<Vec<LMC::Message>>();
                let gamma2_message_temp = gamma2_message.clone();
                let gamma2_beta_inv = cfg_iter!(gamma2_message_temp)
                    .map(|gamma| mul_helper(gamma, &beta_inv)).collect::<Vec<RMC::Message>>();

                // println!("v1 len : {}", v1.len());
                for i in 0..v1.len() {
                    v1[i] = v1[i].clone() + gamma1_beta[i].clone();
                    v2[i] = v2[i].clone() + gamma2_beta_inv[i].clone();
                }

                // compute C and message rescale

                let cr = start_timer!(|| "Compute C");

                let v1_l = v1[..split].to_vec();
                let v1_r = v1[split..].to_vec();
                let v2_l = v2[..split].to_vec();
                let v2_r = v2[split..].to_vec();

                let c_plus = IP::inner_product(&v1_l, &v2_r).unwrap();
                let c_minus = IP::inner_product(&v1_r, &v2_l).unwrap();
           
                end_timer!(cr);

                // Second Fiat-Shamir challenge
                let (alpha, alpha_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input.extend_from_slice(&to_bytes![
                        c_plus, c_minus
                    ]?);
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
                    .map(|b| mul_helper(b, &alpha_inv))
                    .zip(v2_r)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<RMC::Message>>();
                end_timer!(rescale_v2);

                // gamma1 = gamma1_prime.to_vec();
                // gamma2 = gamma2_prime.to_vec();
                gamma1_message = gamma1_message[..split].to_vec();
                gamma2_message = gamma2_message[..split].to_vec();



                let com1 = (d1_l, d2_l, c_plus);
                let com2 = (d1_r, d2_r, c_minus);

                r_commitment_steps.push((com1, com2));
                r_transcript.push((alpha, alpha_inv, beta, beta_inv));
                end_timer!(recurse);
            }

        };

        // r_transcript.reverse();
        r_commitment_steps.reverse();
        Ok((
            HPAProof {
                r_commitment_steps,
                e1: v1,
                e2: v2,
                // r_base: m_base,
                _hpa: PhantomData,
            },
            HPAAux {
                // r_transcript,
                // ck_base,
                _hpa: PhantomData,
            },
        ))
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        proof: &HPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar)>, Error> {
        Self::_compute_recursive_challenges(proof)
    }

    fn _compute_recursive_challenges(
        proof: &HPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar)>, Error> {
        // let (mut com1, mut com2) = (proof.r_commitment_steps[0], proof.r_commitment_steps[1]);

        let mut r_transcript = Vec::new();
        for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
            // First Fiat-Shamir challenge
            let mut counter_nonce: usize = 0;
            // let default_transcript = (Default::default(), Default::default(),Default::default(),Default::default());
            // let transcript = r_transcript.last().unwrap_or(&default_transcript);
            let (beta, beta_inv) = 'challenge: loop {
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                hash_input.extend_from_slice(&to_bytes![
                    com_1.0, com_2.0, com_1.1, com_2.1
                ]?);
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
                hash_input.extend_from_slice(&to_bytes![
                    com_1.2, com_2.2
                ]?);
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

        
            r_transcript.push((alpha, alpha_inv, beta, beta_inv));
        }
        r_transcript.reverse();
        
        Ok(r_transcript)
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
    //     proof: &DORYProof<IP, LMC, RMC, IPC, D>,
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
            e1: self.e1.clone(),
            e2: self.e2.clone(),
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