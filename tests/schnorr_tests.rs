//! BIP 340 Schnorr sign/verify tests.

use blvm_secp256k1::schnorr::{
    schnorr_sign, schnorr_verify, schnorr_verify_batch, xonly_pubkey_from_secret,
};

fn hex32(s: &str) -> [u8; 32] {
    let v = hex::decode(s).unwrap();
    v.try_into().unwrap()
}

fn hex64(s: &str) -> [u8; 64] {
    let v = hex::decode(s).unwrap();
    v.try_into().unwrap()
}

#[test]
fn test_schnorr_bip340_vectors() {
    let vectors: &[(&str, &str, &str, &str, &str, bool)] = &[
        (
            "0000000000000000000000000000000000000000000000000000000000000003",
            "F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9",
            "0000000000000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000000000000000000000000000000000000000000000",
            "E907831F80848D1069A5371B402410364BDF1C5F8307B0084C55F1CE2DCA821525F66A4A85EA8B71E482A74F382D2CE5EBEEE8FDB2172F477DF4900D310536C0",
            true,
        ),
        (
            "B7E151628AED2A6ABF7158809CF4F3C762E7160F38B4DA56A784D9045190CFEF",
            "DFF1D77F2A671C5F36183726DB2341BE58FEAE1DA2DECED843240F7B502BA659",
            "0000000000000000000000000000000000000000000000000000000000000001",
            "243F6A8885A308D313198A2E03707344A4093822299F31D0082EFA98EC4E6C89",
            "6896BD60EEAE296DB48A229FF71DFE071BDE413E6D43F917DC8DCF8C78DE33418906D11AC976ABCCB20B091292BFF4EA897EFCB639EA871CFA95F6DE339E4B0A",
            true,
        ),
        (
            "C90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B14E5C9",
            "DD308AFEC5777E13121FA72B9CC1B7CC0139715309B086C960E18FD969774EB8",
            "C87AA53824B4D7AE2EB035A2B5BBBCCC080E76CDC6D1692C4B0B62D798E6D906",
            "7E2D58D8B3BCDF1ABADEC7829054F90DDA9805AAB56C77333024B9D0A508B75C",
            "5831AAEED7B44BB74E5EAB94BA9D4294C49BCF2A60728D8B4C200F50DD313C1BAB745879A5AD954A72C45A91C3A51D3C7ADEA98D82F8481E0E1E03674A6F3FB7",
            true,
        ),
        (
            "0B432B2677937381AEF05BB02A66ECD012773062CF3FA2549E44F58ED2401710",
            "25D1DFF95105F5253C4022F628A996AD3A0D95FBF21D468A1B33F8C160D8F517",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "7EB0509757E246F19449885651611CB965ECC1A187DD51B64FDA1EDC9637D5EC97582B9CB13DB3933705B32BA982AF5AF25FD78881EBB32771FC5922EFC66EA3",
            true,
        ),
    ];

    for (seckey_hex, pk_hex, aux_hex, msg_hex, sig_hex, expect_ok) in vectors {
        let pk = hex32(pk_hex);
        let msg = hex::decode(msg_hex).unwrap();
        let sig = hex64(sig_hex);

        assert_eq!(
            schnorr_verify(&sig, &msg, &pk),
            *expect_ok,
            "verify mismatch for pk={}",
            pk_hex
        );

        if !seckey_hex.is_empty() {
            let seckey = hex32(seckey_hex);
            let aux = hex32(aux_hex);
            let computed_pk = xonly_pubkey_from_secret(&seckey).expect("xonly_pubkey");
            assert_eq!(computed_pk, pk, "pubkey mismatch");

            let our_sig = schnorr_sign(&seckey, &msg, &aux).expect("sign");
            assert!(schnorr_verify(&our_sig, &msg, &pk), "our sig should verify");
        }
    }
}

#[test]
fn test_schnorr_verify_fail_vectors() {
    let pk = hex32("DFF1D77F2A671C5F36183726DB2341BE58FEAE1DA2DECED843240F7B502BA659");
    let msg =
        hex::decode("243F6A8885A308D313198A2E03707344A4093822299F31D0082EFA98EC4E6C89").unwrap();

    let fail_sig = "FFF97BD5755EEEA420453A14355235D382F6472F8568A18B2F057A14602975563CC27944640AC607CD107AE10923D9EF7A73C643E166BE5EBEAFA34B1AC553E2";
    let sig = hex64(fail_sig);
    assert!(
        !schnorr_verify(&sig, &msg, &pk),
        "should fail: has_even_y(R) false"
    );
}

#[test]
fn test_schnorr_our_sign_matches_bip340_vector0() {
    let seckey = hex32("0000000000000000000000000000000000000000000000000000000000000003");
    let aux = hex32("0000000000000000000000000000000000000000000000000000000000000000");
    let msg =
        hex::decode("0000000000000000000000000000000000000000000000000000000000000000").unwrap();
    let expected_sig = hex64("E907831F80848D1069A5371B402410364BDF1C5F8307B0084C55F1CE2DCA821525F66A4A85EA8B71E482A74F382D2CE5EBEEE8FDB2172F477DF4900D310536C0");

    let our_sig = schnorr_sign(&seckey, &msg, &aux).expect("sign");
    assert_eq!(
        our_sig, expected_sig,
        "our sign should match BIP 340 vector 0"
    );
}

#[test]
fn test_schnorr_sign_verify_roundtrip() {
    let seckey = hex32("B7E151628AED2A6ABF7158809CF4F3C762E7160F38B4DA56A784D9045190CFEF");
    let aux = [0u8; 32];
    let msg = b"test message";

    let pk = xonly_pubkey_from_secret(&seckey).unwrap();
    let sig = schnorr_sign(&seckey, msg, &aux).unwrap();
    assert!(schnorr_verify(&sig, msg, &pk));
}

#[test]
fn test_schnorr_verify_batch() {
    let seckey = hex32("0000000000000000000000000000000000000000000000000000000000000001");
    let pk = xonly_pubkey_from_secret(&seckey).unwrap();
    let msg1 = b"hello";
    let msg2 = b"world";
    let aux = [0u8; 32];
    let sig1 = schnorr_sign(&seckey, msg1, &aux).unwrap();
    let sig2 = schnorr_sign(&seckey, msg2, &aux).unwrap();

    assert!(schnorr_verify(&sig1, msg1, &pk), "individual verify sig1");
    assert!(schnorr_verify(&sig2, msg2, &pk), "individual verify sig2");

    let sigs = vec![sig1, sig2];
    let msgs: Vec<&[u8]> = vec![msg1, msg2];
    let pubkeys = vec![pk, pk];

    assert!(
        schnorr_verify_batch(&sigs, &msgs, &pubkeys),
        "batch should verify"
    );

    // Large batch (90 sigs) exercises Pippenger path
    let mut large_sigs = Vec::with_capacity(90);
    let mut large_msgs: Vec<Vec<u8>> = Vec::with_capacity(90);
    let mut large_pubkeys = Vec::with_capacity(90);
    for i in 0..90 {
        let msg = format!("batch message {}", i);
        let sig = schnorr_sign(&seckey, msg.as_bytes(), &aux).unwrap();
        large_sigs.push(sig);
        large_msgs.push(msg.into_bytes());
        large_pubkeys.push(pk);
    }
    let large_msgs_refs: Vec<&[u8]> = large_msgs.iter().map(|m| m.as_slice()).collect();
    assert!(
        schnorr_verify_batch(&large_sigs, &large_msgs_refs, &large_pubkeys),
        "large batch (90 sigs, Pippenger path) should verify"
    );
}
