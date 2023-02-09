extern crate ark_ff;
use self::ark_ff::{Zero, UniformRand};
// extern crate ark_serialize;
// use self::ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
extern crate ark_std;
use self::ark_std::rand::Rng;
// use self::ark_std::{end_timer, start_timer};
extern  crate digest;
use self::digest::Digest;
use std::{marker::PhantomData, ops::MulAssign}; //convert::TryInto

use crate::{mul_helper, Error};//, InnerProductArgumentError};
extern crate ark_dh_commitments;
use self::ark_dh_commitments::DoublyHomomorphicCommitment;
extern crate ark_inner_products;
use self::ark_inner_products::InnerProduct;
// use self::ark_std::cfg_iter;

use std::fmt;

#[cfg(feature = "parallel")]
// extern crate rayon;
// use self::rayon::prelude::*;

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

pub struct MVP<IP, LMC, RMC, IPC, D> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _right_commitment: PhantomData<RMC>,
    _inner_product_commitment: PhantomData<IPC>,
    _digest: PhantomData<D>,
}

// #[derive(CanonicalSerialize, CanonicalDeserialize)]
// pub struct DORYProof<IP, LMC, RMC, IPC, D>
// where
//     D: Digest,
//     IP: InnerProduct<
//         LeftMessage = LMC::Message,
//         RightMessage = RMC::Message,
//         Output = IPC::Message,
//     >,
//     LMC: DoublyHomomorphicCommitment,
//     RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     RMC::Message: MulAssign<LMC::Scalar>,
//     IPC::Message: MulAssign<LMC::Scalar>,
//     RMC::Key: MulAssign<LMC::Scalar>,
//     IPC::Key: MulAssign<LMC::Scalar>,
//     RMC::Output: MulAssign<LMC::Scalar>,
//     IPC::Output: MulAssign<LMC::Scalar>,
// {
//     pub(crate) r_commitment_steps: Vec<(
//         (IP::Output, IP::Output, IPC::Message), 
//         (IP::Output, IP::Output, IPC::Message),
//     )>,
//     pub(crate) e1: Vec<IP::LeftMessage>,
//     pub(crate) e2: Vec<IP::RightMessage>,
//     // pub(crate) r_base: (LMC::Message, RMC::Message),
//     _dory: PhantomData<DORY<IP, LMC, RMC, IPC, D>>,
// }


// #[derive(Clone)]
// pub struct DORYSRS<IP, LMC, RMC, IPC, D>
// where
//     D: Digest,
//     IP: InnerProduct<
//         LeftMessage = LMC::Message,
//         RightMessage = RMC::Message,
//         Output = IPC::Message,
//     >,
//     LMC: DoublyHomomorphicCommitment,
//     RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     RMC::Message: MulAssign<LMC::Scalar>,
//     IPC::Message: MulAssign<LMC::Scalar>,
//     RMC::Key: MulAssign<LMC::Scalar>,
//     IPC::Key: MulAssign<LMC::Scalar>,
//     RMC::Output: MulAssign<LMC::Scalar>,
//     IPC::Output: MulAssign<LMC::Scalar>,
// {
//     pub(crate) delta1_l: Vec<IP::Output>,
//     pub(crate) delta1_r: Vec<IP::Output>,
//     pub(crate) delta2_l: Vec<IP::Output>,
//     pub(crate) delta2_r: Vec<IP::Output>,
//     pub(crate) kai: Vec<IP::Output>,
    
//     _dory: PhantomData<DORY<IP, LMC, RMC, IPC, D>>,
// }

// #[derive(Clone)]
// pub struct DORYAux<IP, LMC, RMC, IPC, D>
// where
//     D: Digest,
//     IP: InnerProduct<
//         LeftMessage = LMC::Message,
//         RightMessage = RMC::Message,
//         Output = IPC::Message,
//     >,
//     LMC: DoublyHomomorphicCommitment,
//     RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
//     RMC::Message: MulAssign<LMC::Scalar>,
//     IPC::Message: MulAssign<LMC::Scalar>,
//     RMC::Key: MulAssign<LMC::Scalar>,
//     IPC::Key: MulAssign<LMC::Scalar>,
//     RMC::Output: MulAssign<LMC::Scalar>,
//     IPC::Output: MulAssign<LMC::Scalar>,
// {
//     // pub(crate) r_transcript: Vec<(LMC::Scalar, LMC::Scalar, LMC::Scalar, LMC::Scalar)>,
//     // pub(crate) ck_base: (LMC::Key, RMC::Key),
//     _dory: PhantomData<DORY<IP, LMC, RMC, IPC, D>>,
// }

//TODO: Can extend DORY to support "identity commitments" in addition to "compact commitments", i.e. for SIPP

impl<IP, LMC, RMC, IPC, D> MVP<IP, LMC, RMC, IPC, D>
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
    // collapse_matrix(&matrix_a, &gamma2)
    pub fn collapse_matrix2(
        matrix_a: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // r: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        gamma2: &[IP::RightMessage],
        // generator_g2: &IP::RightMessage,
        size: usize,
    ) -> Result<
            // Vec<IP::LeftMessage>,
            // Vec<IP::RightMessage>,
            // Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>
        , Error >
    {
        let matrix_a = matrix_a.clone();
        let gamma2 = gamma2.clone();
        // let generator_g2 = generator_g2.clone();
        let mut v_a = Vec::new();
        // let zero = <LMC::Scalar>::zero();

        // for i in 0..size{
        //     let mut temp = mul_helper(&generator_g2, &zero);
        //     for j in 0..size{
        //         if matrix_a[j*size + i] != zero{
        //             let temp2 = mul_helper(&gamma2[j], &matrix_a[j*size+i]);
        //             temp = temp + temp2;
        //         }
        //     }
        //     v_a.push(temp);
        // }

        for i in 0..size{
            let temp = mul_helper(&gamma2[i], &matrix_a[i]);
            v_a.push(temp);
        }

        Ok(v_a)
    }

    // collapse_matrix(&matrix_a, &gamma2)
    pub fn collapse_matrix1(
        matrix_a: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // r: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        gamma1: &[IP::LeftMessage],
        // generator_g2: &IP::RightMessage,
        size: usize,
    ) -> Result<
            // Vec<IP::LeftMessage>,
            // Vec<IP::RightMessage>,
            // Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>
        , Error >
    {
        let matrix_a = matrix_a.clone();
        let gamma1 = gamma1.clone();
        // let generator_g2 = generator_g2.clone();
        let mut v_a = Vec::new();
        // let zero = <LMC::Scalar>::zero();

        // for i in 0..size{
        //     let mut temp = mul_helper(&generator_g2, &zero);
        //     for j in 0..size{
        //         if matrix_a[j*size + i] != zero{
        //             let temp2 = mul_helper(&gamma2[j], &matrix_a[j*size+i]);
        //             temp = temp + temp2;
        //         }
        //     }
        //     v_a.push(temp);
        // }

        for i in 0..size{
            let temp = mul_helper(&gamma1[i], &matrix_a[i]);
            v_a.push(temp);
        }

        Ok(v_a)
    }

    // set sparse matrix
    pub fn set_sparse_matrix<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<Vec<LMC::Scalar>, Error> {
    let mut matrix_a = Vec::new();
    for i in 0..size{
        for j in 0..size{
            if i == j {
                matrix_a.push(<LMC::Scalar>::rand(rng));
            }
            else{
                matrix_a.push(<LMC::Scalar>::zero());
            }
        }
    }

    Ok(
        matrix_a
    )
    }

    // generator_g^vec_a
    pub fn set_vector(
        vec: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // r: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // generator_g1: &IP::LeftMessage,
        generator_g1: &IP::LeftMessage,
        // size: usize,
    ) -> Result<
            // Vec<IP::LeftMessage>,
            // Vec<IP::RightMessage>,
            // Vec<IP::LeftMessage>,
            Vec<IP::LeftMessage>
        , Error >
    {
        let vec = vec.clone();
        let generator_g1 = generator_g1.clone();

        // let matrix_a = matrix_a.clone();
        // let generator_g2 = generator_g2.clone();
        let mut a_vec = Vec::new();
        // let zero = <LMC as DoublyHomomorphicCommitment::Scalar>::zero();

        for i in 0..vec.len(){
            a_vec.push(mul_helper(&generator_g1, &vec[i]));
        }

        Ok(a_vec)
    }

    pub fn set_vector2(
        vec: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // r: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // generator_g1: &IP::LeftMessage,
        generator_g2: &IP::RightMessage,
        // size: usize,
    ) -> Result<
            // Vec<IP::LeftMessage>,
            // Vec<IP::RightMessage>,
            // Vec<IP::LeftMessage>,
            Vec<IP::RightMessage>
        , Error >
    {
        let vec = vec.clone();
        let generator_g2 = generator_g2.clone();

        // let matrix_a = matrix_a.clone();
        // let generator_g2 = generator_g2.clone();
        let mut a_vec = Vec::new();
        // let zero = <LMC as DoublyHomomorphicCommitment::Scalar>::zero();

        for i in 0..vec.len(){
            a_vec.push(mul_helper(&generator_g2, &vec[i]));
        }

        Ok(a_vec)
    }

    // A * z = a ::: compute a
    pub fn compute_az(
        matrix_a: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        z: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // r: &[<LMC as DoublyHomomorphicCommitment>::Scalar],
        // generator_g1: &IP::LeftMessage,
        // gamma1: &[IP::LeftMessage],
        // size: usize,
    ) -> Result<
            // Vec<IP::LeftMessage>,
            // Vec<IP::RightMessage>,
            // Vec<IP::LeftMessage>,
            Vec<LMC::Scalar>
        , Error >
    {
        let matrix_a = matrix_a.clone();
        let z = z.clone();
        // let gamma1 = gamma1.clone();
        let mut a_vec = Vec::new();
        let size = z.len();
        // let zero = <LMC as DoublyHomomorphicCommitment::Scalar>::zero();

        // for i in 0..size{
        //     let mut temp = <LMC::Scalar>::zero();
        //     for j in 0..size{
        //         temp = temp + matrix_a[i*size + j] * z[j];
        //     }
        //     a_vec.push(temp);
        // }

        for i in 0..size{
            let temp = matrix_a[i] * z[i];
                    
            a_vec.push(temp);
        }

        Ok(a_vec)
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