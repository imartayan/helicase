use helicase::dna_format::BitMask;

pub fn bits(mask: &BitMask) -> Vec<bool> {
    (0..mask.len()).map(|i| mask.get(i)).collect()
}

// Reports the index of the first differing event instead of the whole outputs.
pub fn assert_same<T: PartialEq>(expected: &[T], got: &[T], id: &str) {
    let len = expected.len().max(got.len());
    if let Some(i) = (0..len).find(|&i| expected.get(i) != got.get(i)) {
        panic!(
            "{id}: event {i} differs ({} vs {} events)",
            expected.len(),
            got.len()
        );
    }
}
