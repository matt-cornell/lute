# Lute
This project is going to be a general chemistry library, though only a fraction of my goals are implemented right now.  
Molecules here are represented as undirected graphs from `petgraph`, with many methods being generic over the exact type of graph used. One big feature is the `Arena`, which can store and deduplicate molecules, which should allow for quick substructure checks.  
Molecules can be input and output using SMILES (output canonicalization isn't quite working yet, though), and can be rendered using `coordgen` to generate atom coordinates, or drawn to a bitmap format using `resvg` (with features of the same name).  
Eventually, I plan to compute more complex properties and add the ability to simulate reactions, but that still seems to be a long way off. For now, only features can be found, look over the `Molecule` trait to see what exactly can be done.