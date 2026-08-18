//! BIP327 NonceGen vectors (subset matching our API shape).

use blvm_secp256k1::musig::nonce_gen;

fn hex32(s: &str) -> [u8; 32] {
    hex::decode(s).unwrap().try_into().unwrap()
}

fn hex33(s: &str) -> [u8; 33] {
    hex::decode(s).unwrap().try_into().unwrap()
}

fn hex64(s: &str) -> [u8; 64] {
    hex::decode(s).unwrap().try_into().unwrap()
}

fn hex66(s: &str) -> [u8; 66] {
    hex::decode(s).unwrap().try_into().unwrap()
}

/// BIP327 nonce_gen_vectors.json case 3: no sk / aggpk / msg / extra_in.
#[test]
fn bip327_nonce_gen_case_no_sk() {
    let mut rand_ = hex32("0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F");
    let pk = hex33("02F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9");
    let expected_secnonce = hex::decode(
        "89BDD787D0284E5E4D5FC572E49E316BAB7E21E3B1830DE37DFE80156FA41A6D\
         0B17AE8D024C53679699A6FD7944D9C4A366B514BAF43088E0708B1023DD2897\
         02F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9",
    )
    .unwrap();
    let expected_pubnonce = hex66(
        "02C96E7CB1E8AA5DAC64D872947914198F607D90ECDE5200DE52978AD5DED63C00\
         0299EC5117C2D29EDEE8A2092587C3909BE694D5CFF0667D6C02EA4059F7CD9786",
    );

    let ((sec_k, sec_pk), pubnonce) =
        nonce_gen(&mut rand_, None, &pk, None, None, None).expect("nonce_gen");
    assert_eq!(rand_, [0u8; 32], "session_secrand must be wiped");
    assert_eq!(sec_pk, pk);
    let mut sec = [0u8; 97];
    sec[..64].copy_from_slice(&sec_k);
    sec[64..].copy_from_slice(&sec_pk);
    assert_eq!(sec.as_slice(), expected_secnonce.as_slice());
    assert_eq!(pubnonce, expected_pubnonce);
}

/// Exercises BIP327 `rand = sk XOR hash_MuSig/aux(rand')` (no aggpk/extra).
#[test]
fn bip327_nonce_gen_sk_xor_aux_no_aggpk() {
    let mut rand_ = hex32("0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F");
    let sk = hex32("0202020202020202020202020202020202020202020202020202020202020202");
    let pk = hex33("024D4B6CD1361032CA9BD2AEB9D900AA4D45D9EAD80AC9423374C451A7254D0766");
    let msg = hex32("0101010101010101010101010101010101010101010101010101010101010101");
    // Expected from BIP327 reference.py nonce_gen_internal(sk, pk, None, msg, None).
    let expected_sec_k = hex64(
        "2cd082e87affdf4ce950f9a405661bb817d6c20b6540ebd16a3f244c6bb2ee02\
         a532083c1c33d81cae81001cb6f8d74900e9b70f561c5152b4c5cb48062ad78e",
    );
    let expected_pubnonce = hex66(
        "03f015b1ad984d3972f3575963e890acc026770eb1622ca8742310dcdd60acd334\
         02a5809feecd8ee3e7f2c671cf27f5cde34c6e6d2ed66c16266483abde99568070",
    );

    let ((sec_k, sec_pk), pubnonce) =
        nonce_gen(&mut rand_, Some(&sk), &pk, Some(&msg), None, None).expect("nonce_gen");
    assert_eq!(sec_pk, pk);
    assert_eq!(sec_k, expected_sec_k);
    assert_eq!(pubnonce, expected_pubnonce);
}
