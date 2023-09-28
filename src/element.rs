use std::{collections::HashMap, hash::Hash, sync::LazyLock};

use crate::quantity::Quantity;

#[derive(Clone, Debug, PartialEq)]
#[allow(unused)]
pub(crate) struct Element {
    atomic_number: usize,
    name: &'static str,
    symbol: &'static str,
    mass: Quantity<f64>,
}

impl Hash for Element {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.atomic_number.hash(state);
        self.name.hash(state);
        self.symbol.hash(state);
    }
}

impl Eq for Element {}

macro_rules! e {
    ($a:expr, $n:expr, $s:expr, $m:expr) => {
        Element::new($a, $n, $s, $m)
    };
}

impl Element {
    pub(crate) fn new(
        atomic_number: usize,
        name: &'static str,
        symbol: &'static str,
        mass: Quantity<f64>,
    ) -> Self {
        Self {
            atomic_number,
            name,
            symbol,
            mass,
        }
    }
}

pub fn element(symbol: &'static str) -> Element {
    BY_SYMBOL[symbol].clone()
}

pub(crate) static EP: LazyLock<Element> =
    LazyLock::new(|| Element::new(0, "EP", "EP", Quantity::new(0.0, "EP")));

pub(crate) static BY_SYMBOL: LazyLock<HashMap<&str, Element>> =
    LazyLock::new(|| {
        HashMap::from([
            ("H", e!(1, "hydrogen", "H", q!(1.007947, "daltons"))),
            ("D", e!(1, "deuterium", "D", q!(2.01355321270, "daltons"))),
            ("He", e!(2, "helium", "He", q!(4.003, "daltons"))),
            ("Li", e!(3, "lithium", "Li", q!(6.9412, "daltons"))),
            ("Be", e!(4, "beryllium", "Be", q!(9.0121823, "daltons"))),
            ("B", e!(5, "boron", "B", q!(10.8117, "daltons"))),
            ("C", e!(6, "carbon", "C", q!(12.01078, "daltons"))),
            ("N", e!(7, "nitrogen", "N", q!(14.00672, "daltons"))),
            ("O", e!(8, "oxygen", "O", q!(15.99943, "daltons"))),
            ("F", e!(9, "fluorine", "F", q!(18.99840325, "daltons"))),
            ("Ne", e!(10, "neon", "Ne", q!(20.17976, "daltons"))),
            ("Na", e!(11, "sodium", "Na", q!(22.989769282, "daltons"))),
            ("Mg", e!(12, "magnesium", "Mg", q!(24.30506, "daltons"))),
            ("Al", e!(13, "aluminum", "Al", q!(26.98153868, "daltons"))),
            ("Si", e!(14, "silicon", "Si", q!(28.08553, "daltons"))),
            ("P", e!(15, "phosphorus", "P", q!(30.9737622, "daltons"))),
            ("S", e!(16, "sulfur", "S", q!(32.0655, "daltons"))),
            ("Cl", e!(17, "chlorine", "Cl", q!(35.4532, "daltons"))),
            ("Ar", e!(18, "argon", "Ar", q!(39.9481, "daltons"))),
            ("K", e!(19, "potassium", "K", q!(39.09831, "daltons"))),
            ("Ca", e!(20, "calcium", "Ca", q!(40.0784, "daltons"))),
            ("Sc", e!(21, "scandium", "Sc", q!(44.9559126, "daltons"))),
            ("Ti", e!(22, "titanium", "Ti", q!(47.8671, "daltons"))),
            ("V", e!(23, "vanadium", "V", q!(50.94151, "daltons"))),
            ("Cr", e!(24, "chromium", "Cr", q!(51.99616, "daltons"))),
            ("Mn", e!(25, "manganese", "Mn", q!(54.9380455, "daltons"))),
            ("Fe", e!(26, "iron", "Fe", q!(55.8452, "daltons"))),
            ("Co", e!(27, "cobalt", "Co", q!(58.9331955, "daltons"))),
            ("Ni", e!(28, "nickel", "Ni", q!(58.69342, "daltons"))),
            ("Cu", e!(29, "copper", "Cu", q!(63.5463, "daltons"))),
            ("Zn", e!(30, "zinc", "Zn", q!(65.4094, "daltons"))),
            ("Ga", e!(31, "gallium", "Ga", q!(69.7231, "daltons"))),
            ("Ge", e!(32, "germanium", "Ge", q!(72.641, "daltons"))),
            ("As", e!(33, "arsenic", "As", q!(74.921602, "daltons"))),
            ("Se", e!(34, "selenium", "Se", q!(78.963, "daltons"))),
            ("Br", e!(35, "bromine", "Br", q!(79.9041, "daltons"))),
            ("Kr", e!(36, "krypton", "Kr", q!(83.7982, "daltons"))),
            ("Rb", e!(37, "rubidium", "Rb", q!(85.46783, "daltons"))),
            ("Sr", e!(38, "strontium", "Sr", q!(87.621, "daltons"))),
            ("Y", e!(39, "yttrium", "Y", q!(88.905852, "daltons"))),
            ("Zr", e!(40, "zirconium", "Zr", q!(91.2242, "daltons"))),
            ("Nb", e!(41, "niobium", "Nb", q!(92.906382, "daltons"))),
            ("Mo", e!(42, "molybdenum", "Mo", q!(95.942, "daltons"))),
            ("Tc", e!(43, "technetium", "Tc", q!(98.0, "daltons"))),
            ("Ru", e!(44, "ruthenium", "Ru", q!(101.072, "daltons"))),
            ("Rh", e!(45, "rhodium", "Rh", q!(102.905502, "daltons"))),
            ("Pd", e!(46, "palladium", "Pd", q!(106.421, "daltons"))),
            ("Ag", e!(47, "silver", "Ag", q!(107.86822, "daltons"))),
            ("Cd", e!(48, "cadmium", "Cd", q!(112.4118, "daltons"))),
            ("In", e!(49, "indium", "In", q!(114.8183, "daltons"))),
            ("Sn", e!(50, "tin", "Sn", q!(118.7107, "daltons"))),
            ("Sb", e!(51, "antimony", "Sb", q!(121.7601, "daltons"))),
            ("Te", e!(52, "tellurium", "Te", q!(127.603, "daltons"))),
            ("I", e!(53, "iodine", "I", q!(126.904473, "daltons"))),
            ("Xe", e!(54, "xenon", "Xe", q!(131.2936, "daltons"))),
            ("Cs", e!(55, "cesium", "Cs", q!(132.90545192, "daltons"))),
            ("Ba", e!(56, "barium", "Ba", q!(137.3277, "daltons"))),
            ("La", e!(57, "lanthanum", "La", q!(138.905477, "daltons"))),
            ("Ce", e!(58, "cerium", "Ce", q!(140.1161, "daltons"))),
            (
                "Pr",
                e!(59, "praseodymium", "Pr", q!(140.907652, "daltons")),
            ),
            ("Nd", e!(60, "neodymium", "Nd", q!(144.2423, "daltons"))),
            ("Pm", e!(61, "promethium", "Pm", q!(145.0, "daltons"))),
            ("Sm", e!(62, "samarium", "Sm", q!(150.362, "daltons"))),
            ("Eu", e!(63, "europium", "Eu", q!(151.9641, "daltons"))),
            ("Gd", e!(64, "gadolinium", "Gd", q!(157.253, "daltons"))),
            ("Tb", e!(65, "terbium", "Tb", q!(158.925352, "daltons"))),
            ("Dy", e!(66, "dysprosium", "Dy", q!(162.5001, "daltons"))),
            ("Ho", e!(67, "holmium", "Ho", q!(164.930322, "daltons"))),
            ("Er", e!(68, "erbium", "Er", q!(167.2593, "daltons"))),
            ("Tm", e!(69, "thulium", "Tm", q!(168.934212, "daltons"))),
            ("Yb", e!(70, "ytterbium", "Yb", q!(173.043, "daltons"))),
            ("Lu", e!(71, "lutetium", "Lu", q!(174.9671, "daltons"))),
            ("Hf", e!(72, "hafnium", "Hf", q!(178.492, "daltons"))),
            ("Ta", e!(73, "tantalum", "Ta", q!(180.947882, "daltons"))),
            ("W", e!(74, "tungsten", "W", q!(183.841, "daltons"))),
            ("Re", e!(75, "rhenium", "Re", q!(186.2071, "daltons"))),
            ("Os", e!(76, "osmium", "Os", q!(190.233, "daltons"))),
            ("Ir", e!(77, "iridium", "Ir", q!(192.2173, "daltons"))),
            ("Pt", e!(78, "platinum", "Pt", q!(195.0849, "daltons"))),
            ("Au", e!(79, "gold", "Au", q!(196.9665694, "daltons"))),
            ("Hg", e!(80, "mercury", "Hg", q!(200.592, "daltons"))),
            ("Tl", e!(81, "thallium", "Tl", q!(204.38332, "daltons"))),
            ("Pb", e!(82, "lead", "Pb", q!(207.21, "daltons"))),
            ("Bi", e!(83, "bismuth", "Bi", q!(208.980401, "daltons"))),
            ("Po", e!(84, "polonium", "Po", q!(209.0, "daltons"))),
            ("At", e!(85, "astatine", "At", q!(210.0, "daltons"))),
            ("Rn", e!(86, "radon", "Rn", q!(222.018, "daltons"))),
            ("Fr", e!(87, "francium", "Fr", q!(223.0, "daltons"))),
            ("Ra", e!(88, "radium", "Ra", q!(226.0, "daltons"))),
            ("Ac", e!(89, "actinium", "Ac", q!(227.0, "daltons"))),
            ("Th", e!(90, "thorium", "Th", q!(232.038062, "daltons"))),
            (
                "Pa",
                e!(91, "protactinium", "Pa", q!(231.035882, "daltons")),
            ),
            ("U", e!(92, "uranium", "U", q!(238.028913, "daltons"))),
            ("Np", e!(93, "neptunium", "Np", q!(237.0, "daltons"))),
            ("Pu", e!(94, "plutonium", "Pu", q!(244.0, "daltons"))),
            ("Am", e!(95, "americium", "Am", q!(243.0, "daltons"))),
            ("Cm", e!(96, "curium", "Cm", q!(247.0, "daltons"))),
            ("Bk", e!(97, "berkelium", "Bk", q!(247.0, "daltons"))),
            ("Cf", e!(98, "californium", "Cf", q!(251.0, "daltons"))),
            ("Es", e!(99, "einsteinium", "Es", q!(252.0, "daltons"))),
            ("Fm", e!(100, "fermium", "Fm", q!(257.0, "daltons"))),
            ("Md", e!(101, "mendelevium", "Md", q!(258.0, "daltons"))),
            ("No", e!(102, "nobelium", "No", q!(259.0, "daltons"))),
            ("Lr", e!(103, "lawrencium", "Lr", q!(262.0, "daltons"))),
            ("Rf", e!(104, "rutherfordium", "Rf", q!(261.0, "daltons"))),
            ("Db", e!(105, "dubnium", "Db", q!(262.0, "daltons"))),
            ("Sg", e!(106, "seaborgium", "Sg", q!(266.0, "daltons"))),
            ("Bh", e!(107, "bohrium", "Bh", q!(264.0, "daltons"))),
            ("Hs", e!(108, "hassium", "Hs", q!(269.0, "daltons"))),
            ("Mt", e!(109, "meitnerium", "Mt", q!(268.0, "daltons"))),
            ("Ds", e!(110, "darmstadtium", "Ds", q!(281.0, "daltons"))),
            ("Rg", e!(111, "roentgenium", "Rg", q!(272.0, "daltons"))),
            ("Uub", e!(112, "ununbium", "Uub", q!(285.0, "daltons"))),
            ("Uut", e!(113, "ununtrium", "Uut", q!(284.0, "daltons"))),
            ("Uuq", e!(114, "ununquadium", "Uuq", q!(289.0, "daltons"))),
            ("Uup", e!(115, "ununpentium", "Uup", q!(288.0, "daltons"))),
            ("Uuh", e!(116, "ununhexium", "Uuh", q!(292.0, "daltons"))),
        ])
    });
