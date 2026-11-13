//! Named refusals (port of errors.py). Python raises `CCA2Error(identifier, message)`; the port
//! returns `Err(CErr)` carrying the SAME identifier string, so `try/except CCA2Error` becomes a
//! `match` on the identifier.

#[derive(Debug, Clone)]
pub struct CErr {
    pub id: String,
    pub msg: String,
}

impl std::fmt::Display for CErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.msg.is_empty() {
            write!(f, "{}", self.id)
        } else {
            write!(f, "{}: {}", self.id, self.msg)
        }
    }
}

impl std::error::Error for CErr {}

pub type R<T> = Result<T, CErr>;

pub fn cerr(id: &str, msg: impl Into<String>) -> CErr {
    CErr { id: id.to_string(), msg: msg.into() }
}

/// `error(id, msg)`: return early with a named refusal.
#[macro_export]
macro_rules! bail {
    ($id:expr) => { return Err($crate::errors::cerr($id, "")) };
    ($id:expr, $($arg:tt)*) => { return Err($crate::errors::cerr($id, format!($($arg)*))) };
}
