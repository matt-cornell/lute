use crate::prelude::*;
use rand::prelude::*;

#[test]
fn simple() {
    assert_eq!(smiles!("C").empirical(), empirical!("CH4"));
    assert_eq!(smiles!("C=C").empirical(), empirical!("C2H4"));
    assert_eq!(smiles!("CC#C").empirical(), empirical!("C3H4"));
}

#[test]
fn format() {
    assert_eq!(&empirical!("CH4").to_string(), "CH4");
    assert_eq!(&empirical!("C6H12O6").to_string(), "C6H12O6");
}

#[cfg(feature = "rand")]
#[test]
fn roundtrip() {
    for _ in 0..100 {
        let e: EmpiricalFormula = rand::thread_rng().sample(crate::rand::Salt::default());
        let fmt = e.to_string();
        assert_eq!(fmt.parse::<EmpiricalFormula>(), Ok(e));
    }
}
