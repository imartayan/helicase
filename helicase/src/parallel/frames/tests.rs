use super::zstd_frame_sizes;
use zstd::zstd_safe::{find_frame_compressed_size, get_frame_content_size};

fn compress(data: &[u8], level: i32, content_size: bool, checksum: bool) -> Vec<u8> {
    use zstd::stream::raw::CParameter;
    let mut compressor = zstd::bulk::Compressor::new(level).unwrap();
    compressor
        .set_parameter(CParameter::ContentSizeFlag(content_size))
        .unwrap();
    compressor
        .set_parameter(CParameter::ChecksumFlag(checksum))
        .unwrap();
    compressor.compress(data).unwrap()
}

#[test]
fn zstd_frame_sizes_match_libzstd() {
    let mut random = vec![0u8; 700_000];
    let mut state = 0x9e3779b97f4a7c15u64;
    for b in random.iter_mut() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        *b = state as u8;
    }
    let dna: Vec<u8> = random.iter().map(|&b| b"ACGT"[(b & 3) as usize]).collect();
    let inputs: [&[u8]; 7] = [
        b"",
        b"A",
        &[b'N'; 300],      // RLE block, 2-byte content size
        &[b'N'; 1 << 20],  // several RLE blocks
        &random[..70_000], // raw block, 4-byte content size
        &random,           // several raw blocks
        &dna,              // compressed blocks
    ];
    for data in inputs {
        for (level, content_size, checksum) in
            [(1, true, false), (3, false, true), (19, false, false)]
        {
            let frame = compress(data, level, content_size, checksum);
            let (comp_len, uncomp_len) = zstd_frame_sizes(&frame).unwrap();
            assert_eq!(comp_len, find_frame_compressed_size(&frame).unwrap());
            match get_frame_content_size(&frame).unwrap() {
                Some(len) => assert_eq!(uncomp_len, len as usize),
                None => assert!(uncomp_len >= data.len()),
            }
            let mut out = vec![0; uncomp_len];
            let len = zstd::bulk::Decompressor::new()
                .unwrap()
                .decompress_to_buffer(&frame, &mut out[..])
                .unwrap();
            assert_eq!(&out[..len], data);
            for cut in 0..frame.len() {
                assert_eq!(zstd_frame_sizes(&frame[..cut]), None, "cut {cut}");
            }
        }
    }
}
