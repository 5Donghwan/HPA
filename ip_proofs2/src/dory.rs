use ark_ec::PairingEngine;
use ark_ff::{to_bytes, Field, One};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer};
use digest::Digest;
use std::{convert::TryInto, marker::PhantomData, ops::MulAssign};

use crate::{mul_helper, Error, InnerProductArgumentError, add_helper};
use ark_dh_commitments::DoublyHomomorphicCommitment;
use ark_inner_products::InnerProduct;
use ark_std::cfg_iter;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct DORY<IP, LMC, RMC, IPC, D> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _right_commitment: PhantomData<RMC>,
    _inner_product_commitment: PhantomData<IPC>,
    _digest: PhantomData<D>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct DORYProof<IP, LMC, RMC, IPC, D>
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
        (LMC::Output, RMC::Output, IPC::Output),
        (LMC::Output, RMC::Output, IPC::Output),
    )>,
    pub(crate) r_base: (LMC::Message, RMC::Message),
    _dory: PhantomData<DORY<IP, LMC, RMC, IPC, D>>,
}

#[derive(Clone)]
pub struct DORYAux<IP, LMC, RMC, IPC, D>
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
    pub(crate) r_transcript: Vec<LMC::Scalar>,
    pub(crate) ck_base: (LMC::Key, RMC::Key),
    _dory: PhantomData<DORY<IP, LMC, RMC, IPC, D>>,
}

//TODO: Can extend DORY to support "identity commitments" in addition to "compact commitments", i.e. for SIPP

impl<IP, LMC, RMC, IPC, D> DORY<IP, LMC, RMC, IPC, D>
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
    pub fn setup<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<(Vec<LMC::Key>, Vec<RMC::Key>, IPC::Key), Error> {
        Ok((
            LMC::setup(rng, size)?,
            RMC::setup(rng, size)?,
            IPC::setup(rng, 1)?.pop().unwrap(),
        ))
    }

    pub fn prove(
        values: (&[IP::LeftMessage], &[IP::RightMessage], &IP::Output),
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
    ) -> Result<DORYProof<IP, LMC, RMC, IPC, D>, Error> {
        if IP::inner_product(values.0, values.1)? != values.2.clone() {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }
        if values.0.len().count_ones() != 1 {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                values.0.len(),
                values.1.len(),
            )));
        }
        if !(LMC::verify(ck.0, values.0, com.0)?
            && RMC::verify(ck.1, values.1, com.1)?
            && IPC::verify(&vec![ck.2.clone()], &vec![values.2.clone()], com.2)?)
        {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }

        let (proof, _) =
            Self::prove_with_aux((values.0, values.1), (ck.0, ck.1, &vec![ck.2.clone()]))?;
        Ok(proof)
    }

    pub fn verify(
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &DORYProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<bool, Error> {
        if ck.0.len().count_ones() != 1 || ck.0.len() != ck.1.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                ck.0.len(),
                ck.1.len(),
            )));
        }
        // Calculate base commitment and transcript
        let (base_com, transcript) = Self::_compute_recursive_challenges(
            (com.0.clone(), com.1.clone(), com.2.clone()),
            proof,
        )?;
        // Calculate base commitment keys
        let (ck_a_base, ck_b_base) = Self::_compute_final_commitment_keys(ck, &transcript)?;
        // Verify base commitment
        Self::_verify_base_commitment(
            (&ck_a_base, &ck_b_base, &vec![ck.2.clone()]),
            base_com,
            proof,
        )
    }

    pub fn prove_with_aux(
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &[RMC::Key], &[IPC::Key]),
    ) -> Result<
        (
            DORYProof<IP, LMC, RMC, IPC, D>,
            DORYAux<IP, LMC, RMC, IPC, D>,
        ),
        Error,
    > {
        let (m_a, m_b) = values;
        let (ck_a, ck_b, ck_t) = ck;
        
        Self::_prove(
            (m_a.to_vec(), m_b.to_vec()),
            (ck_a.to_vec(), ck_b.to_vec(), ck_t.to_vec()),
        )
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove(
        values: (Vec<IP::LeftMessage>, Vec<IP::RightMessage>),
        ck: (Vec<RMC::Key>, Vec<LMC::Key>),
        ck_message: (Vec<LMC::Message>, Vec<RMC::Message>),
    ) -> Result<
        (
            DORYProof<IP, LMC, RMC, IPC, D>,
            DORYAux<IP, LMC, RMC, IPC, D>,
        ),
        Error,
    > {
        let (mut v1, mut v2) = values;
        let (mut gamma1, mut gamma2) = ck;
        let (mut gamma1_message, mut gamma2_message) = ck_message;
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        assert!(v1.len().is_power_of_two());
        let (m_base, ck_base) = 'recurse: loop {
            let recurse = start_timer!(|| format!("Recurse round size {}", m_a.len()));
            if v1.len() == 1 {
                // base case
                break 'recurse (
                    (v1[0].clone(), v2[0].clone()),
                    (gamma1[0].clone(), gamma2[0].clone()),
                );
            } else {
                // recursive step
                // Recurse with problem of half size
                let split = v1.len() / 2;

                let v1_L = &v1[split..];
                let v1_R = &v1[..split];
                let gamma1_prime = &gamma1[..split];
                let gamma1_message_prime = &gamma1_message[..split];

                let v2_L = &v2[..split];
                let v2_R = &v2[split..];
                let gamma2_prime = &gamma2[split..];
                let gamma2_message_prime = &gamma2_message[split..];

                let cl = start_timer!(|| "Compute D");
                let d1_L = LMC::commit(gamma2_prime, v1_L)?;
                let d1_R = LMC::commit(gamma2_prime, v1_R)?;
                let d2_L = RMC::commit(gamma1_prime, v2_L)?;
                let d2_R = RMC::commit(gamma1_prime, v2_R)?;
                
                 // Fiat-Shamir challenge
                 let mut counter_nonce: usize = 0;
                 let default_transcript = Default::default();
                 let transcript = r_transcript.last().unwrap_or(&default_transcript);
                 let (beta, beta_inv) = 'challenge: loop {
                     let mut hash_input = Vec::new();
                     hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                     //TODO: Should use CanonicalSerialize instead of ToBytes
                     hash_input.extend_from_slice(&to_bytes![
                         transcript, d1_L, d1_R, d2_L, d2_R
                     ]?);
                     let beta: LMC::Scalar = u128::from_be_bytes(
                         D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                     )
                     .into();
                     if let Some(beta_inv) = beta.inverse() {
                         // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                         // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                         break 'challenge (beta_inv, beta);
                     }
                     counter_nonce += 1;
                 };
                 
                end_timer!(cl);

                let gamma1_message_temp = cfg_iter!(gamma1_message)
                    .map(|gamma| mul_helper(gamma, &beta)).collect::<Vec<LMC::Message>>();


                let gamma2_message_temp = cfg_iter!(gamma2_message)
                    .map(|gamma| mul_helper(gamma, &beta_inv)).collect::<Vec<RMC::Message>>();

                for i in 0..v1.len() {
                    v1[i] = v1[i] + gamma1_message_temp[i];
                    v2[i] = v2[i] + gamma2_message_temp[i];
                }

                // compute C and message rescale


                let cr = start_timer!(|| "Compute C");

                let c_plus = vec![IP::inner_product(&v1_L, &v2_R).unwrap()];
                let c_minus = vec![IP::inner_product(&v1_R, &v2_L).unwrap()];
           
                end_timer!(cr);

                // Second Fiat-Shamir challenge
                let (alpha, alpha_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input.extend_from_slice(&to_bytes![
                        transcript, c_plus, c_minus
                    ]?);
                    let alpha: LMC::Scalar = u128::from_be_bytes(
                        D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                    )
                    .into();
                    if let Some(alpha_inv) = alpha.inverse() {
                        // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                        // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                        break 'challenge (alpha_inv, alpha);
                    }
                    counter_nonce += 1;
                };

                // Set up values for next step of recursion
                let rescale_v1 = start_timer!(|| "Rescale V1");
                v1 = cfg_iter!(v1_L)
                    .map(|a| mul_helper(a, &alpha))
                    .zip(v1_R)
                    .map(|(a_1, a_2)| a_1 + a_2.clone())
                    .collect::<Vec<LMC::Message>>();
                end_timer!(rescale_v1);

                let rescale_v2 = start_timer!(|| "Rescale V2");
                v2 = cfg_iter!(v2_L)
                    .map(|b| mul_helper(b, &alpha_inv))
                    .zip(v2_R)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<RMC::Message>>();
                end_timer!(rescale_v2);

                r_commitment_steps.push((d1_L, d1_R, d2_L, d2_R, c_plus, c_minus));
                r_transcript.push(alpha, beta);
                end_timer!(recurse);
            }
        };
        r_transcript.reverse();
        r_commitment_steps.reverse();
        Ok((
            DORYProof {
                r_commitment_steps,
                r_base: m_base,
                _dory: PhantomData,
            },
            DORYAux {
                r_transcript,
                ck_base,
                _dory: PhantomData,
            },
        ))
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &DORYProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<((LMC::Output, RMC::Output, IPC::Output), Vec<LMC::Scalar>), Error> {
        Self::_compute_recursive_challenges((com.0.clone(), com.1.clone(), com.2.clone()), proof)
    }

    fn _compute_recursive_challenges(
        com: (LMC::Output, RMC::Output, IPC::Output),
        proof: &DORYProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<((LMC::Output, RMC::Output, IPC::Output), Vec<LMC::Scalar>), Error> {
        let (mut com_a, mut com_b, mut com_t) = com;
        let mut r_transcript = Vec::new();
        for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
            // Fiat-Shamir challenge
            let mut counter_nonce: usize = 0;
            let default_transcript = Default::default();
            let transcript = r_transcript.last().unwrap_or(&default_transcript);
            let (c, c_inv) = 'challenge: loop {
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                hash_input.extend_from_slice(&to_bytes![
                    transcript, com_1.0, com_1.1, com_1.2, com_2.0, com_2.1, com_2.2
                ]?);
                let c: LMC::Scalar = u128::from_be_bytes(
                    D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                )
                .into();
                if let Some(c_inv) = c.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (c_inv, c);
                }
                counter_nonce += 1;
            };

            com_a = mul_helper(&com_1.0, &c) + com_a.clone() + mul_helper(&com_2.0, &c_inv);
            com_b = mul_helper(&com_1.1, &c) + com_b.clone() + mul_helper(&com_2.1, &c_inv);
            com_t = mul_helper(&com_1.2, &c) + com_t.clone() + mul_helper(&com_2.2, &c_inv);

            r_transcript.push(c);
        }
        r_transcript.reverse();
        Ok(((com_a, com_b, com_t), r_transcript))
    }

    pub(crate) fn _compute_final_commitment_keys(
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        transcript: &Vec<LMC::Scalar>,
    ) -> Result<(LMC::Key, RMC::Key), Error> {
        // Calculate base commitment keys
        let (ck_a, ck_b, _) = ck;
        assert!(ck_a.len().is_power_of_two());

        let mut ck_a_agg_challenge_exponents = vec![LMC::Scalar::one()];
        let mut ck_b_agg_challenge_exponents = vec![LMC::Scalar::one()];
        for (i, c) in transcript.iter().enumerate() {
            let c_inv = c.inverse().unwrap();
            for j in 0..(2_usize).pow(i as u32) {
                ck_a_agg_challenge_exponents.push(ck_a_agg_challenge_exponents[j] * &c_inv);
                ck_b_agg_challenge_exponents.push(ck_b_agg_challenge_exponents[j] * c);
            }
        }
        assert_eq!(ck_a_agg_challenge_exponents.len(), ck_a.len());
        //TODO: Optimization: Use VariableMSM multiexponentiation
        let ck_a_base_init = mul_helper(&ck_a[0], &ck_a_agg_challenge_exponents[0]);
        let ck_a_base = ck_a[1..]
            .iter()
            .zip(&ck_a_agg_challenge_exponents[1..])
            .map(|(g, x)| mul_helper(g, &x))
            .fold(ck_a_base_init, |sum, x| sum + x);
        //.reduce(|| ck_a_base_init.clone(), |sum, x| sum + x);
        let ck_b_base_init = mul_helper(&ck_b[0], &ck_b_agg_challenge_exponents[0]);
        let ck_b_base = ck_b[1..]
            .iter()
            .zip(&ck_b_agg_challenge_exponents[1..])
            .map(|(g, x)| mul_helper(g, &x))
            .fold(ck_b_base_init, |sum, x| sum + x);
        //.reduce(|| ck_b_base_init.clone(), |sum, x| sum + x);
        Ok((ck_a_base, ck_b_base))
    }

    pub(crate) fn _verify_base_commitment(
        base_ck: (&LMC::Key, &RMC::Key, &Vec<IPC::Key>),
        base_com: (LMC::Output, RMC::Output, IPC::Output),
        proof: &DORYProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<bool, Error> {
        let (com_a, com_b, com_t) = base_com;
        let (ck_a_base, ck_b_base, ck_t) = base_ck;
        let a_base = vec![proof.r_base.0.clone()];
        let b_base = vec![proof.r_base.1.clone()];
        let t_base = vec![IP::inner_product(&a_base, &b_base)?];

        Ok(LMC::verify(&vec![ck_a_base.clone()], &a_base, &com_a)?
            && RMC::verify(&vec![ck_b_base.clone()], &b_base, &com_b)?
            && IPC::verify(&ck_t, &t_base, &com_t)?)
    }
}

impl<IP, LMC, RMC, IPC, D> Clone for DORYProof<IP, LMC, RMC, IPC, D>
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
        DORYProof {
            r_commitment_steps: self.r_commitment_steps.clone(),
            r_base: self.r_base.clone(),
            _dory: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_ff::UniformRand;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use blake2::Blake2b;

    use ark_dh_commitments::{
        afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
        identity::IdentityCommitment,
        pedersen::PedersenCommitment,
        random_generators,
    };
    use ark_inner_products::{
        ExtensionFieldElement, InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct,
        ScalarInnerProduct,
    };

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
    type SC2 = PedersenCommitment<<Bls12_381 as PairingEngine>::G2Projective>;
    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC =
            IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>;
        type PairingDORY = DORY<IP, GC1, GC2, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = PairingDORY::setup(&mut rng, TEST_SIZE).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = GC2::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = PairingDORY::prove(
            (&m_a, &m_b, &t[0]),
            (&ck_a, &ck_b, &ck_t),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            PairingDORY::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
        );
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
        type IPC = IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >;
        type MultiExpDORY = DORY<IP, GC1, SC1, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = MultiExpDORY::setup(&mut rng, TEST_SIZE).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_b.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
        }
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = SC1::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = MultiExpDORY::prove(
            (&m_a, &m_b, &t[0]),
            (&ck_a, &ck_b, &ck_t),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            MultiExpDORY::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
        );
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>;
        type IPC =
            IdentityCommitment<<Bls12_381 as PairingEngine>::Fr, <Bls12_381 as PairingEngine>::Fr>;
        type ScalarDORY = DORY<IP, SC2, SC2, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = ScalarDORY::setup(&mut rng, TEST_SIZE).unwrap();
        let mut m_a = Vec::new();
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
            m_b.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
        }
        let com_a = SC2::commit(&ck_a, &m_a).unwrap();
        let com_b = SC2::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = ScalarDORY::prove(
            (&m_a, &m_b, &t[0]),
            (&ck_a, &ck_b, &ck_t),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            ScalarDORY::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
        );
    }
}
