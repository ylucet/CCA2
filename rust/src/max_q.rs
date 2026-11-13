//! max_q -- the pointwise maximum of two QuaCons, exactly (port of max_q.py).

use crate::cells::{assemble_qua_con_cells, con_row, Cell};
use crate::errors::R;
use crate::qua_con::QuaCon;
use crate::rat_q::*;
use num_traits::Signed;

fn rows_of(g: &QuaCon, k: usize) -> Vec<Vec<Z>> {
    // MATLAB keeps [conic, side] and DROPS the branch column here; kept verbatim.
    g.fc[k]
        .iter()
        .map(|r| {
            let mut v = g.ecq[r[0] as usize - 1].clone();
            v.push(zi(r[1]));
            v
        })
        .collect()
}

pub fn max_q(g1: &QuaCon, g2: &QuaCon) -> R<QuaCon> {
    let mut cells = vec![];
    for i in 0..g1.nf {
        let ci = rows_of(g1, i);
        for j in 0..g2.nf {
            let cj = rows_of(g2, j);
            let mut base = ci.clone();
            base.extend(cj.iter().cloned());
            let (n1, d1) = (&g1.fn_[i], &g1.fd[i]);
            let (n2, d2) = (&g2.fn_[j], &g2.fd[j]);
            let (dn, _) = sub(n1, d1, n2, d2)?;
            if all_zero(&dn) {
                cells.push(Cell { num: n1.clone(), den: d1.clone(), con: base });
                continue;
            }
            let r = face_conic(&dn);
            let rc = conic(&r)?;
            let nz = r.iter().find(|x| !num_traits::Zero::is_zero(*x)).unwrap();
            let flip = if nz.is_negative() { -1 } else { 1 };
            let mut c1 = base.clone();
            c1.push({
                let mut v = rc.clone();
                v.push(zi(flip));
                v
            });
            let mut c2 = base;
            c2.push({
                let mut v = rc;
                v.push(zi(-flip));
                v
            });
            cells.push(Cell { num: n1.clone(), den: d1.clone(), con: c1 });
            cells.push(Cell { num: n2.clone(), den: d2.clone(), con: c2 });
        }
    }
    let _ = con_row;
    assemble_qua_con_cells(cells)
}
