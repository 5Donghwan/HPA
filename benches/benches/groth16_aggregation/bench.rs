use ark_bls12_377::{
    constraints::PairingVar as BLS12PairingVar, Bls12_377, Fr as BLS12Fr,
    FrParameters as BLS12FrParameters,
};
use ark_bw6_761::Fr as BW6Fr;
use ark_crypto_primitives::{
    snark::*,
    prf::{
        blake2s::{constraints::Blake2sGadget,Blake2s},
        constraints::PRFGadget,
        PRF,
    },
};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{
    biginteger::BigInteger, FftParameters, One, PrimeField, ToConstraintField, UniformRand,
};
use ark_groth16::{constraints::*, Groth16, PreparedVerifyingKey, Proof, VerifyingKey};
use ark_r1cs_std::{fields::fp::FpVar, prelude::*};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};

use ark_bls12_381::{Bls12_381,Fr, G1Projective};
use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
    // pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};

use ark_inner_products::{
    ExtensionFieldElement, InnerProduct, PairingInnerProduct,
};
use ark_hpa::hpa::HPA;
// use ark_ip_proofs::applications::groth16_aggregation::{
//     aggregate_proofs, setup_inner_product, verify_aggregate_proof,
// };

use ark_std::rand::{rngs::StdRng, SeedableRng,Rng};
use blake2::Blake2b;
use digest::Digest;
use std::{ops::MulAssign, time::Instant,env};

#[derive(Clone)]
// struct TestCircuit {
//     public_inputs: Vec<Fr>,
//     witness_input: Fr,
//     public_sum: Fr,
// }
// impl <F: PrimeField>ConstraintSynthesizer<F> for TestCircuit {
//     fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
//         let input_variables =
//             Vec::<F>::new_input_vec(cs.clone(), || Ok(self.public_inputs.clone()))?;
//         let sum = FpVar::new_input(cs.clone(), || Ok(&self.public_sum))?;
//         let witness = FpVar::new_witness(cs.clone(), || Ok(&self.witness_input))?;

//         let mut computed_sum = witness;
//         for x in &input_variables {
//             computed_sum += x;
//         }

//         sum.enforce_equal(&computed_sum)?;

//         Ok(())
//     }
// }

#[derive(Default)]
struct SingleBlake2SCircuit {
    input: [u8; 32],
    output: [u8; 32],
}

impl<F: PrimeField> ConstraintSynthesizer<F> for SingleBlake2SCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let seed = UInt8::constant_vec(&[0; 32]);
        let input = UInt8::new_witness_vec(cs.clone(), &self.input)?;
        let hash = <Blake2sGadget as PRFGadget<_, F>>::OutputVar::new_variable(
            cs.clone(),
            || Ok(self.output.clone()),
            AllocationMode::Input,
        )?;
        hash.enforce_equal(&<Blake2sGadget as PRFGadget<_, F>>::evaluate(
            &seed, &input,
        )?)?;
        Ok(())
    }
}

#[derive(Clone)]
struct ManyBlake2SCircuit {
    input: Vec<[u8; 32]>,
    output: Vec<[u8; 32]>,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for ManyBlake2SCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let seed = UInt8::constant_vec(&[0; 32]);

        for (hash_input, hash_output) in self.input.iter().zip(&self.output) {
            let input = UInt8::new_witness_vec(cs.clone(), hash_input)?;
            let hash = <Blake2sGadget as PRFGadget<_, F>>::OutputVar::new_variable(
                cs.clone(),
                || Ok(hash_output.clone()),
                AllocationMode::Input,
            )?;
            hash.enforce_equal(&<Blake2sGadget as PRFGadget<_, F>>::evaluate(
                &seed, &input,
            )?)?;
        }
        Ok(())
    }
}

#[derive(Clone)]
struct AggregateBlake2SCircuitVerificationCircuit {
    hash_outputs: Vec<[u8; 32]>,
    proofs: Vec<Proof<Bls12_377>>,
    vk: VerifyingKey<Bls12_377>,
}

impl ConstraintSynthesizer<BW6Fr> for AggregateBlake2SCircuitVerificationCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<BW6Fr>) -> Result<(), SynthesisError> {
        let input_gadgets = self
            .hash_outputs
            .iter()
            .map(|h| h.to_field_elements())
            .collect::<Option<Vec<Vec<BLS12Fr>>>>()
            .unwrap()
            .iter()
            .map(|h_as_bls_fr| {
                let h_as_bls_fr_bytes = h_as_bls_fr
                    .iter()
                    .map(|bls_fr| {
                        bls_fr
                            .into_repr()
                            .as_ref()
                            .iter()
                            .map(|bls_fr_int| bls_fr_int.to_le_bytes().to_vec())
                            .collect::<Vec<Vec<u8>>>()
                    })
                    .collect::<Vec<Vec<Vec<u8>>>>()
                    .iter()
                    .flatten()
                    .flatten()
                    .cloned()
                    .collect::<Vec<u8>>();
                UInt8::new_input_vec(cs.clone(), &h_as_bls_fr_bytes)
            })
            .collect::<Result<Vec<Vec<UInt8<BW6Fr>>>, SynthesisError>>()?
            // Allocated as BLS12-377 field elements byte representation packed together to BW6-761 field elements
            // Now split BW6-761 byte representation back to iterator over BLS12-377 field element byte representations
            .iter()
            .map(|h_as_bls_fr_bytes| {
                let bls_field_element_size_in_bytes =
                    <BLS12FrParameters as FftParameters>::BigInt::NUM_LIMBS * 8;
                h_as_bls_fr_bytes
                    .chunks(bls_field_element_size_in_bytes)
                    .map(|bls_field_element_chunk| bls_field_element_chunk.to_vec())
                    .collect::<Vec<Vec<UInt8<BW6Fr>>>>()
            })
            .collect::<Vec<Vec<Vec<UInt8<BW6Fr>>>>>();

        let vk_gadget =
            VerifyingKeyVar::<Bls12_377, BLS12PairingVar>::new_constant(cs.clone(), &self.vk)?;

        let proof_gadgets =
            self.proofs
                .iter()
                .map(|proof| {
                    ProofVar::<Bls12_377, BLS12PairingVar>::new_witness(cs.clone(), || {
                        Ok(proof.clone())
                    })
                })
                .collect::<Result<Vec<ProofVar<Bls12_377, BLS12PairingVar>>, SynthesisError>>()?;

        assert_eq!(input_gadgets.len(), proof_gadgets.len());

        for (input_gadget, proof_gadget) in input_gadgets.iter().zip(&proof_gadgets) {
            let input = input_gadget
                .iter()
                .map(|bytes| {
                    bytes
                        .iter()
                        .flat_map(|b| b.to_bits_le().unwrap())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let input = BooleanInputVar::new(input);
            <Groth16VerifierGadget<Bls12_377, BLS12PairingVar>>::verify(
                &vk_gadget,
                &input,
                proof_gadget,
            )?
            .enforce_equal(&Boolean::TRUE)?;
        }

        Ok(())
    }
}

struct AggregateBlake2SCircuitVerificationCircuitInput {
    hash_outputs: Vec<[u8; 32]>,
}

impl ToConstraintField<BW6Fr> for AggregateBlake2SCircuitVerificationCircuitInput {
    fn to_field_elements(&self) -> Option<Vec<BW6Fr>> {
        let mut fr_elements: Vec<BW6Fr> = vec![];

        for h_as_bls_fr_bytes in self
            .hash_outputs
            .iter()
            .map(|h| h.to_field_elements())
            .collect::<Option<Vec<Vec<BLS12Fr>>>>()?
            .iter()
            .map(|h_as_bls_fr| {
                h_as_bls_fr
                    .iter()
                    .map(|bls_fr| {
                        bls_fr
                            .into_repr()
                            .as_ref()
                            .iter()
                            .map(|bls_fr_int| bls_fr_int.to_le_bytes().to_vec())
                            .collect::<Vec<Vec<u8>>>()
                    })
                    .collect::<Vec<Vec<Vec<u8>>>>()
                    .iter()
                    .flatten()
                    .flatten()
                    .cloned()
                    .collect::<Vec<u8>>()
            })
        {
            fr_elements.extend_from_slice(&h_as_bls_fr_bytes.to_field_elements()?);
        }

        Some(fr_elements)
    }
}


fn bench_hpa_aggregate<IP, LMC, RMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
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
    let mut start = Instant::now();
    let mut time;
    let mut rng_: StdRng = StdRng::seed_from_u64(0u64);
    let generate_all_proofs = true;
    const NUM_PUBLIC_INPUTS: usize = 16;


    let mut hash_inputs = vec![];
    let mut hash_outputs = vec![];
    for _ in 0..len {
        let hash_input = UniformRand::rand(&mut rng_);
        hash_outputs.push(Blake2s::evaluate(&[0; 32], &hash_input).unwrap());
        hash_inputs.push(hash_input);
    }
    //let parameters = Groth16::<P>::setup(test_circuit, &mut rng_).unwrap();
    let hash_circuit_parameters =
            Groth16::<P>::setup(SingleBlake2SCircuit::default(), &mut rng_).unwrap();
        time = start.elapsed().as_millis();
    // Generate parameters for inner product aggregation
    //let srs = setup_inner_product::<_, Blake2b, _>(&mut rng, len).unwrap();

    // Generate proofs
    println!("Generating {} Groth16 proofs...", len);
    
   // Prove individual proofs
        start = Instant::now();
        let hash_circuit_parameters =
            Groth16::<P>::setup(SingleBlake2SCircuit::default(), &mut rng_).unwrap();
        time = start.elapsed().as_millis();
        if generate_all_proofs {
            // csv_writer
            //     .write_record(&[
            //         1.to_string(),
            //         num_proofs.to_string(),
            //         "single_circuit".to_string(),
            //         "setup".to_string(),
            //         time.to_string(),
            //     ])
            //     .unwrap();
            // csv_writer.flush().unwrap();
        }

        start = Instant::now();
        let mut proofs: Vec<Proof<P>> = vec![];
        for i in 0..len {
            proofs.push(
                Groth16::<P>::prove(
                    &hash_circuit_parameters.0,
                    SingleBlake2SCircuit {
                        input: hash_inputs[i].clone(),
                        output: hash_outputs[i].clone(),
                    },
                    &mut rng_,
                )
                .unwrap(),
            );
            if !generate_all_proofs {
                break;
            }
            //assert!(Groth16::<Bls12_377, SingleBlake2SCircuit, [u8; 32]>::verify(&hash_circuit_parameters.1, &hash_outputs[i], &proofs[i]).unwrap());
        }
        time = start.elapsed().as_millis();

        if !generate_all_proofs {
            proofs = vec![proofs[0].clone(); len];
            hash_inputs = vec![hash_inputs[0]; len];
            hash_outputs = vec![hash_outputs[0]; len];
        } else {
            // csv_writer
            //     .write_record(&[
            //         1.to_string(),
            //         len.to_string(),
            //         "single_circuit".to_string(),
            //         "prove".to_string(),
            //         time.to_string(),
            //     ])
            //     .unwrap();
            // csv_writer.flush().unwrap();
        


        //let result = Groth16::<Bls12_381, TestCircuit, [Fr]>::verify(&parameters.1, &statement, &proof).unwrap();
        //assert!(result);
    }
    
    //let generator_g1 = <IP::LeftMessage>::rand(rng);
    //let generator_g2 = <IP::RightMessage>::rand(rng);
    
    //let (v1, v2, u1, u2) = HPA::<IP, LMC, RMC, IPC, D>::set_values(&l, &r, &generator_g1, &generator_g2).unwrap();
    let a = proofs
        .iter()
        .map(|proof| proof.a.into_projective())
        // .collect::<Vec<P::G1Projective>>();
        .collect::<Vec<<P as PairingEngine>::G1Projective>>();
    let b = proofs
        .iter()
        .map(|proof| proof.b.into_projective())
        // .collect::<Vec<P::G2Projective>>();
        .collect::<Vec<<P as PairingEngine>::G2Projective>>();
    let cc: Vec<<P as PairingEngine>::G1Projective> = proofs
        .iter()
        .map(|proof| proof.c.into_projective())
        // .collect::<Vec<P::G1Projective>>();
        .collect::<Vec<<P as PairingEngine>::G1Projective>>();
    let mut _d = vec![hash_circuit_parameters.0.vk.delta_g2.into_projective().clone(); len];

    let mut _h = vec![hash_circuit_parameters.0.vk.gamma_g2.into_projective().clone(); len];

    let (gamma2, gamma1) = HPA::<IP,LMC,RMC,IPC, D>::setup(rng, len).unwrap();
    
    
    let mut h1 = Vec::new();
    let mut h2 = Vec::new();
    h1.push(<IP::LeftMessage>::rand(rng));
    h2.push(<IP::RightMessage>::rand(rng));

    //precomputing commitment keys are done before sending commitments
    let mut hpa_srs = HPA::<IP, LMC, RMC, IPC, D>::precompute((&(gamma1.clone()), &(gamma2.clone())), &h1, &h2).unwrap();
    let mut hpa_srs_ = HPA::<IP, LMC, RMC, IPC, D>::precompute((&(gamma1.clone()), &(gamma2.clone())), &h1, &h2).unwrap();

    let (c, d1, d2,
        x, y, d3, d4,
        gm, gm_vec, r_c, r_d1, r_d2, r_x, r_y, r_d3, r_d4,
        w_vec, k_vec)
        = HPA::<IP,LMC,RMC,IPC,D>::init_commit(&a, &b, &gamma1, &gamma2, &h1, &h2, rng).unwrap();

    let (c_, d1_, d2_,
        x_, y_, d3_, d4_,
        w_vec_, k_vec_)
            = HPA::<IP,LMC,RMC,IPC,D>::init_commit2(&cc, &_d, &gamma1, &gamma2, &h1, &h2, &gm_vec,
            &r_c, &r_d1, &r_d2, &r_x, &r_y, &r_d3, &r_d4, rng).unwrap();

    // X ?= X'

    let bool_x = x == x_;
    println!("X == X' : {}", bool_x);   



    let mut start = Instant::now();
    let mut proof =
        HPA::<IP, LMC, RMC, IPC, D>::prove((&(a.clone()), &(b.clone()), &(w_vec.clone()), &(k_vec.clone())),
            &hpa_srs, 
            (&(gamma1.clone()), &(gamma2.clone())), 
            //  (&(d1.clone()), &(d2.clone()), &(c.clone())),
            (&r_c, &r_x, &r_y, &r_d1, &r_d2, &r_d3, &r_d4),
            &gm,
            rng
        ).unwrap();

    // let mut start = Instant::now();
    let mut proof_ =
        HPA::<IP, LMC, RMC, IPC, D>::prove((&(cc.clone()), &(_d.clone()), &(w_vec_.clone()), &(k_vec_.clone())),
            &hpa_srs_, 
            (&(gamma1.clone()), &(gamma2.clone())), 
            //  (&(d1.clone()), &(d2.clone()), &(c.clone())),
            (&r_c, &r_x, &r_y, &r_d1, &r_d2, &r_d3, &r_d4),
            &gm,
            rng
        ).unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);


    start = Instant::now();
    let result = HPA::<IP, LMC, RMC, IPC, D>::verify(&mut hpa_srs, (&(gamma1.clone()), &(gamma2.clone())),
        (&c, &x, &y, &d1, &d2, &d3, &d4), &mut proof, &gm, rng)
        .unwrap();
    let result2 = HPA::<IP, LMC, RMC, IPC, D>::verify(&mut hpa_srs_, (&(gamma1.clone()), &(gamma2.clone())),
    (&c_, &x_, &y_, &d1_, &d2_, &d3_, &d4_), &mut proof_, &gm, rng)
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
    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    let mut rng = StdRng::seed_from_u64(0u64);

    println!("Benchmarking HPA_Groth16_proof_aggregation: {}", len);

    println!("1) Pairing hadamard product...");
    bench_hpa_aggregate::<
        PairingInnerProduct<Bls12_381>,
        GC1,
        GC2,
        IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, len);
}


pub fn batch_verify_proof<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    public_inputs: &[Vec<E::Fr>],
    proofs: &[Proof<E>],
) -> Result<bool, SynthesisError> {
    let mut rng = StdRng::seed_from_u64(0u64);
    let mut r_powers = Vec::with_capacity(proofs.len());
    for _ in 0..proofs.len() {
        let challenge: E::Fr = u128::rand(&mut rng).into();
        r_powers.push(challenge);
    }

    let combined_inputs = public_inputs
        .iter()
        .zip(&r_powers)
        .map(|(input, r)| {
            let mut g_ic = pvk.vk.gamma_abc_g1[0].into_projective();
            for (i, b) in input.iter().zip(pvk.vk.gamma_abc_g1.iter().skip(1)) {
                g_ic += &b.mul(i.into_repr());
            }
            g_ic.mul(r.into_repr())
        })
        .sum::<E::G1Projective>()
        .into_affine();

    let combined_proof_a_s = proofs
        .iter()
        .zip(&r_powers)
        .map(|(proof, r)| proof.a.mul(*r))
        .collect::<Vec<_>>();
    let combined_proof_a_s = E::G1Projective::batch_normalization_into_affine(&combined_proof_a_s);
    let ml_inputs = proofs
        .iter()
        .zip(&combined_proof_a_s)
        .map(|(proof, a)| ((*a).into(), proof.b.into()))
        .collect::<Vec<_>>();
    let a_r_times_b = E::miller_loop(ml_inputs.iter());

    let combined_c_s = proofs
        .iter()
        .zip(&r_powers)
        .map(|(proof, r)| proof.c.mul(*r))
        .sum::<E::G1Projective>()
        .into_affine();

    let sum_of_rs = (&r_powers).iter().copied().sum::<E::Fr>();
    let combined_alpha = (-pvk.vk.alpha_g1.mul(sum_of_rs)).into_affine();
    let qap = E::miller_loop(
        [
            (combined_alpha.into(), pvk.vk.beta_g2.into()),
            (combined_inputs.into(), pvk.gamma_g2_neg_pc.clone()),
            (combined_c_s.into(), pvk.delta_g2_neg_pc.clone()),
        ]
        .iter(),
    );

    let test =
        E::final_exponentiation(&(qap * &a_r_times_b)).ok_or(SynthesisError::UnexpectedIdentity)?;

    Ok(test == E::Fqk::one())
}