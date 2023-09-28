use torx::Graph;

use super::Topology;

impl From<Topology> for Graph {
    fn from(value: Topology) -> Self {
        let mut g = Graph::new();
        for atom in value.atoms() {
            g.add_node(atom.index);
        }
        for bond in value.bonds() {
            g.add_edge(bond.0.index, bond.1.index);
        }
        g
    }
}
