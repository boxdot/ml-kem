use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};

pub(crate) fn shake256(input: &[u8], out: &mut [u8]) {
    let mut hasher = Shake256::default();
    hasher.update(input);
    hasher.finalize_xof().read(out);
}
