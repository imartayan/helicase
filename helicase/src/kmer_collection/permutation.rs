use std::cmp::min;
use std::mem;

/// returns the min of each non trivial orbits.
/// The permutation must be mutable as elements will be annoted
/// during the algorithm.
/// At the end of the algorithm, the permutation will be restored at its initial value.
///
/// To permut should have the last bit not used free (if usize is 64bits, then permut should
/// be less that 2**63 elements, so less than 73 exabyte, according to a quick computation)
pub fn compute_orbit(permut: &mut [usize]) -> Vec<usize> {
    let mut orbits_min = vec![];
    let mut offset: usize = 0;
    let n = permut.len();
    let mut min_orb: Option<usize> = None;
    let last_bit = 1usize << (mem::size_of::<usize>() * 8 - 1);
    assert!(n < last_bit);
    while offset < n {
        if permut[offset] >= last_bit || permut[offset] == offset {
            if let Some(value) = min_orb {
                orbits_min.push(value)
            }
            offset += 1;
            min_orb = None;
            continue;
        }
        let noffset: usize = permut[offset];
        min_orb = Some(match min_orb {
            Some(value) => min(value, noffset),
            None => min(offset, noffset),
        });
        permut[offset] |= last_bit;
        offset = noffset;
    }
    if let Some(value) = min_orb {
        orbits_min.push(value)
    }
    let mask = !last_bit;
    println!("orbit:first part");
    for el in permut.iter_mut() {
        *el &= mask;
    }
    println!("orbit:second part");
    orbits_min
}

pub fn onto_map_copy<T: Copy>(
    to_permut: &mut [T],
    permutation: &[usize],
    out_size: usize,
) -> Vec<T> {
    let mut out = vec![to_permut[0]; out_size];
    for i in 0..out_size {
        out[i] = to_permut[permutation[i]];
    }
    out
}

pub fn permut_slice_orbit<'a, T: Copy>(
    to_permut: &'a mut [T],
    permutation: &[usize],
    orbits_min: &'a [usize],
) {
    for el in orbits_min.iter() {
        let mut offset = *el;
        let mut noffset = permutation[*el];
        while noffset != *el {
            to_permut.swap(offset, noffset);
            let buff_offset = noffset;
            noffset = permutation[noffset];
            offset = buff_offset;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_compute_orbit() {
        let permut: &mut [usize] = &mut [1, 2, 0, 4, 3, 5, 8, 7, 6];
        let orbits = compute_orbit(permut);
        assert_eq!(orbits.len(), 3, "orbit size");
        assert_eq!(orbits[0], 0);
        assert_eq!(orbits[1], 3);
        assert_eq!(orbits[2], 6);

        let permut2: &mut [usize] = &mut [1, 3, 0, 2];
        let orbits2 = compute_orbit(permut2);
        assert_eq!(orbits2.len(), 1, "orbit size");
        assert_eq!(orbits2[0], 0);
    }

    #[test]
    fn test_permut() {
        let permut: &mut [usize] = &mut [1, 2, 0, 4, 3, 5, 8, 7, 6];
        let to_permut: &mut [usize] = &mut [0, 1, 2, 3, 4, 5, 6, 7, 8];
        let orbits = compute_orbit(permut);
        permut_slice_orbit(to_permut, permut, &orbits);
        assert_eq!(format!("{:?}", to_permut), format!("{:?}", permut));

        let permut2: &mut [usize] = &mut [1, 3, 0, 2];
        let to_permut2: &mut [usize] = &mut [0, 1, 2, 3];
        let orbits2 = compute_orbit(permut2);
        permut_slice_orbit(to_permut2, permut2, &orbits2);
        assert_eq!(format!("{:?}", to_permut2), format!("{:?}", permut2));
    }

    #[test]
    fn test_sort() {
        let data: &mut [char] = &mut ['b', 'a', 'b', 'c'];
        let mut permut = vec![];
        for i in 0..data.len() {
            permut.push(i);
        }
        permut.sort_by_key(|i| data[*i]);
        let orbits = compute_orbit(&mut permut);
        permut_slice_orbit(data, &permut, &orbits);
        assert_eq!(format!("{:?}", data), format!("{:?}", ['a', 'b', 'b', 'c']));
    }
}
