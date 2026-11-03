use super::*;

#[test]
fn test_reader_input_fused_on_error() {
    struct FailingReader {
        data: Vec<u8>,
        pos: usize,
    }
    impl std::io::Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.pos < self.data.len() {
                let n = (self.data.len() - self.pos).min(buf.len());
                buf[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
                self.pos += n;
                Ok(n)
            } else {
                Err(std::io::Error::other("simulated I/O failure"))
            }
        }
    }

    let mut input = ReaderInput::new(FailingReader {
        data: b">test\nACGT\n".to_vec(),
        pos: 0,
    })
    .unwrap();
    // First block: Ok
    let first = input.next();
    assert!(first.is_some() && first.unwrap().is_ok());

    // Next block: Err
    let err = input.next();
    assert!(err.is_some() && err.unwrap().is_err());

    // Fused: all subsequent calls must be None, without attempting read again
    assert!(input.next().is_none());
    assert!(input.next().is_none());
}
