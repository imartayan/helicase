use arrow::array::{Array, ArrayRef, BooleanArray};
use arrow::compute::filter;
use arrow::compute::kernels::boolean::or_kleene;
use arrow::compute::kernels::cmp::neq;
use arrow::compute::kernels::concat::concat;

use eyre::{Result, eyre};

fn downcast_bool(array: &ArrayRef) -> Result<&BooleanArray> {
    array
        .as_any()
        .downcast_ref::<BooleanArray>()
        .ok_or_else(|| eyre!("expected BooleanArray"))
}

fn boundary_mask_for_sorted_columns(columns: &[ArrayRef]) -> Result<BooleanArray> {
    if columns.is_empty() {
        return Ok(BooleanArray::from(Vec::<bool>::new()));
    }

    let n = columns[0].len();
    if columns.iter().any(|c| c.len() != n) {
        return Err(eyre!("all columns must have the same length",));
    }

    if n <= 1 {
        return Ok(BooleanArray::from(vec![true; n]));
    }

    let mut boundary = neq(&columns[0].slice(0, n - 1), &columns[0].slice(1, n - 1))?;
    for col in &columns[1..] {
        let diff = neq(&col.slice(0, n - 1), &col.slice(1, n - 1))?;
        boundary = or_kleene(&boundary, &diff)?;
    }

    Ok(boundary)
}

/// Deduplicate already-sorted key columns.
///
/// Returns:
/// - the deduplicated columns
/// - optionally, a per-run count column (one count per deduplicated row)
pub fn dedup_sorted_columns(columns: &[ArrayRef]) -> Result<Vec<ArrayRef>> {
    if columns.is_empty() {
        return Ok(Vec::new());
    }

    let n = columns[0].len();
    if columns.iter().any(|c| c.len() != n) {
        return Err(eyre!("all columns must have the same length",));
    }

    if n == 0 {
        return Ok(columns.to_vec());
    }

    // boundary.len() == n - 1
    let boundary = boundary_mask_for_sorted_columns(columns)?;

    // keep_mask.len() == n, keep_mask[0] = true, keep_mask[i] = boundary[i-1]
    let first = BooleanArray::from(vec![true]);
    let keep_mask_ref = concat(&[&first, &boundary])?;
    let keep_mask = downcast_bool(&keep_mask_ref)?;

    // Deduplicate columns with one kernel pass per column.
    let mut deduped = Vec::with_capacity(columns.len());
    for col in columns {
        let filtered = filter(col.as_ref(), keep_mask)?;
        deduped.push(filtered);
    }

    Ok(deduped)
}
