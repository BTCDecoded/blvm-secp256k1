//! Rejects the parsers already return for a bad DER tag and a bad pubkey prefix.

use blvm_secp256k1::ecdsa::{ecdsa_sig_parse_der, ge_from_compressed, ge_from_uncompressed};

#[test]
fn ecdsa_sig_parse_der_rejects_an_empty_buffer() {
    match ecdsa_sig_parse_der(&[]) {
        None => {}
        Some(_) => panic!("empty DER parsed as a signature"),
    }
}

#[test]
fn ecdsa_sig_parse_der_rejects_a_non_sequence_tag() {
    match ecdsa_sig_parse_der(&[0x02, 0x00]) {
        None => {}
        Some(_) => panic!("integer tag parsed as a DER signature"),
    }
}

#[test]
fn ge_from_compressed_rejects_a_zero_prefix() {
    let mut bytes = [0u8; 33];
    bytes[0] = 0x00;
    match ge_from_compressed(&bytes) {
        None => {}
        Some(_) => panic!("prefix 0x00 parsed as a compressed public key"),
    }
}

#[test]
fn ge_from_uncompressed_rejects_a_compressed_prefix() {
    let mut bytes = [0u8; 65];
    bytes[0] = 0x02;
    match ge_from_uncompressed(&bytes) {
        None => {}
        Some(_) => panic!("prefix 0x02 parsed as an uncompressed public key"),
    }
}
