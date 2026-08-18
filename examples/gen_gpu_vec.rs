
fn main() {
  use blvm_secp256k1::ecdsa::{ecdsa_sign_compact_rfc6979, ge_to_compressed, pubkey_from_secret, ecdsa_verify_one_compact};
  use blvm_secp256k1::scalar::Scalar;
  let mut sk=[0u8;32]; sk[31]=1;
  let mut msg=[0u8;32]; msg[31]=7;
  let mut s=Scalar::zero(); s.set_b32(&sk);
  let pk=ge_to_compressed(&pubkey_from_secret(&s));
  let sig=ecdsa_sign_compact_rfc6979(&msg,&sk).unwrap();
  assert!(ecdsa_verify_one_compact(&sig,&msg,&pk));
  print!("MSG="); for b in &msg { print!("{:02x}", b); } println!();
  print!("PK="); for b in &pk { print!("{:02x}", b); } println!();
  print!("SIG="); for b in &sig { print!("{:02x}", b); } println!();
}
