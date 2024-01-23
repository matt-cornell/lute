use chem_sim::atomic::SmilesParser;
use petgraph::dot::Dot;

fn main() {
    for arg in std::env::args().skip(1) {
        let parser = SmilesParser::new(&arg);
        match parser.parse() {
            Ok(graph) => println!("{}", Dot::new(&graph)),
            Err(err) => eprintln!("{err}")
        }
    }
}
