//! Every committed row of `bip340_test_vectors.csv`.

use blvm_secp256k1::schnorr::{schnorr_sign, schnorr_verify, xonly_pubkey_from_secret};

const PATH: &str = "tests/data/bip340_test_vectors.csv";
const ROWS: usize = 19;

fn hex_array<const N: usize>(text: &str) -> Option<[u8; N]> {
    let bytes = hex::decode(text).ok()?;
    bytes.try_into().ok()
}

#[test]
fn committed_bip340_vectors_match_schnorr_verify() {
    let text = std::fs::read_to_string(PATH).expect(PATH);
    let mut rows = 0usize;
    for (lineno, line) in text.lines().enumerate() {
        if lineno == 0 || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        assert!(cols.len() >= 7, "short row {lineno}");
        let secret = cols[1];
        let pubkey = cols[2];
        let aux = cols[3];
        let message = hex::decode(cols[4]).unwrap_or_else(|_| panic!("row {lineno} message"));
        let signature = cols[5];
        let expect_ok = cols[6] == "TRUE";
        let got_ok = match (hex_array::<32>(pubkey), hex_array::<64>(signature)) {
            (Some(pk), Some(sig)) => schnorr_verify(&sig, &message, &pk),
            _ => false,
        };
        assert_eq!(got_ok, expect_ok, "row {lineno} {}", cols[7]);

        if !secret.is_empty() {
            let secret = hex_array::<32>(secret).expect("secret");
            let aux = hex_array::<32>(aux).expect("aux");
            let pubkey = hex_array::<32>(pubkey).expect("pubkey");
            assert_eq!(
                xonly_pubkey_from_secret(&secret).expect("xonly"),
                pubkey,
                "row {lineno} pubkey"
            );
            let signed = schnorr_sign(&secret, &message, &aux).expect("sign");
            assert_eq!(
                signed,
                hex_array::<64>(signature).expect("sig"),
                "row {lineno} sign"
            );
        }
        rows += 1;
    }
    assert_eq!(rows, ROWS);
}
