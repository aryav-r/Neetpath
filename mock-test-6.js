// mock-test-6.js - NEET Mock Test 6 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_6 = {
    id: "neet-006",
    title: "Full Syllabus Mock 6", 
    level: "hard",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A body of mass 10 kg moving with velocity 20 m/s collides head-on with a body of mass 20 kg moving with velocity 10 m/s in opposite direction. If the collision is perfectly inelastic, the final velocity is:",
                    options: ["Zero", "3.33 m/s", "6.67 m/s", "10 m/s"],
                    answer: 0,
                    explain: "Conservation of momentum: m₁v₁ - m₂v₂ = (m₁ + m₂)v_final. 10×20 - 20×10 = 30×v_final ⟹ v_final = 0"
                },
                {
                    id: "p2",
                    text: "A uniform disc of radius R rolls without slipping down an inclined plane of angle θ. The acceleration of its center of mass is:",
                    options: ["g sin θ", "(2/3)g sin θ", "(3/4)g sin θ", "(5/7)g sin θ"],
                    answer: 1,
                    explain: "For rolling disc: a = g sin θ/(1 + I/mR²) = g sin θ/(1 + 1/2) = (2/3)g sin θ"
                },
                {
                    id: "p3",
                    text: "The electric field due to an infinite sheet of charge with surface charge density σ is:",
                    options: ["σ/(4πε₀)", "σ/(2ε₀)", "σ/ε₀", "2σ/ε₀"],
                    answer: 1,
                    explain: "Using Gauss's law with cylindrical surface: E = σ/(2ε₀) on both sides of sheet"
                },
                {
                    id: "p4",
                    text: "A simple pendulum has period T on Earth. On Moon (g_moon = g_earth/6), its period becomes:",
                    options: ["T/√6", "T√6", "6T", "T/6"],
                    answer: 1,
                    explain: "T = 2π√(L/g). When g becomes g/6: T' = 2π√(L/(g/6)) = T√6"
                },
                {
                    id: "p5",
                    text: "In a p-n junction diode, the depletion layer thickness increases when:",
                    options: ["Forward biased", "Reverse biased", "No bias", "Temperature increases"],
                    answer: 1,
                    explain: "Reverse bias widens depletion layer as more charge carriers are pulled away from junction"
                },
                {
                    id: "p6",
                    text: "A wave travels with speed v = fλ. If frequency doubles and speed remains same, wavelength:",
                    options: ["Doubles", "Halves", "Remains same", "Becomes 1/4"],
                    answer: 1,
                    explain: "v = fλ. If v constant and f doubles, then λ must halve to maintain equation"
                },
                {
                    id: "p7",
                    text: "The magnetic field inside a toroidal solenoid is:",
                    options: ["Uniform everywhere", "Zero outside, uniform inside", "Maximum at center", "Minimum at center"],
                    answer: 1,
                    explain: "Toroidal field is confined inside torus, zero outside, uniform inside due to symmetry"
                },
                {
                    id: "p8",
                    text: "A charge +q is placed at center of a hollow conducting sphere. The electric field just outside the surface is:",
                    options: ["Zero", "kq/R²", "kq/(4R²)", "Depends on thickness"],
                    answer: 1,
                    explain: "Conductor redistributes charge on outer surface. Field just outside = kq/R² by Gauss's law"
                },
                {
                    id: "p9",
                    text: "The peak value of AC voltage is 220V. Its RMS value is:",
                    options: ["220V", "156V", "110V", "311V"],
                    answer: 1,
                    explain: "V_RMS = V_peak/√2 = 220/√2 = 156V for sinusoidal AC"
                },
                {
                    id: "p10",
                    text: "A particle moves in circular path under centripetal force F = mv²/r. If speed doubles, force becomes:",
                    options: ["2F", "4F", "F/2", "F/4"],
                    answer: 1,
                    explain: "F ∝ v². If v → 2v, then F → 4F"
                },
                {
                    id: "p11",
                    text: "The work done in moving a charge q through potential difference V is:",
                    options: ["qV", "q/V", "V/q", "qV²"],
                    answer: 0,
                    explain: "Work done W = qΔV = qV when moving charge through potential difference"
                },
                {
                    id: "p12",
                    text: "A convex mirror has focal length 20 cm. An object at 30 cm gives image at:",
                    options: ["12 cm", "-12 cm", "60 cm", "-60 cm"],
                    answer: 1,
                    explain: "1/f = 1/u + 1/v. For convex mirror f = -20 cm: 1/(-20) = 1/(-30) + 1/v ⟹ v = -12 cm"
                },
                {
                    id: "p13",
                    text: "The time constant of RC circuit is 2 seconds. After 6 seconds, charge on capacitor becomes:",
                    options: ["1/e times", "1/e² times", "1/e³ times", "Zero"],
                    answer: 2,
                    explain: "q(t) = q₀e^(-t/RC). At t = 6s, τ = 2s: q = q₀e^(-6/2) = q₀e^(-3) = q₀/e³"
                },
                {
                    id: "p14",
                    text: "A satellite orbits Earth at height equal to Earth's radius. Its orbital period compared to geostationary satellite is:",
                    options: ["Less", "More", "Same", "Cannot determine"],
                    answer: 0,
                    explain: "T² ∝ r³. At height R, orbital radius = 2R < geostationary radius, so period is less"
                },
                {
                    id: "p15",
                    text: "The de Broglie wavelength of 1 keV electron is approximately:",
                    options: ["0.39 nm", "3.9 nm", "39 nm", "0.039 nm"],
                    answer: 0,
                    explain: "λ = h/p = h/√(2mE) = 6.63×10⁻³⁴/√(2×9.1×10⁻³¹×1.6×10⁻¹⁶) ≈ 0.39 nm"
                },
                {
                    id: "p16",
                    text: "Two coherent sources produce interference pattern. If path difference is λ/4, the phase difference is:",
                    options: ["π/4", "π/2", "π", "2π"],
                    answer: 1,
                    explain: "Phase difference = (2π/λ) × path difference = (2π/λ) × (λ/4) = π/2"
                },
                {
                    id: "p17",
                    text: "A gas undergoes isothermal expansion from volume V to 2V. Work done by gas is:",
                    options: ["nRT ln 2", "nRT", "2nRT", "nRT/2"],
                    answer: 0,
                    explain: "For isothermal process: W = nRT ln(V₂/V₁) = nRT ln(2V/V) = nRT ln 2"
                },
                {
                    id: "p18",
                    text: "The energy of X-ray photon with wavelength 1 Å is:",
                    options: ["1.24 keV", "12.4 keV", "124 keV", "0.124 keV"],
                    answer: 1,
                    explain: "E = hc/λ = 1240 eV·nm / 0.1 nm = 12.4 keV"
                },
                {
                    id: "p19",
                    text: "A spring-mass system oscillates with amplitude A. At displacement A/2, ratio of kinetic to total energy is:",
                    options: ["1:4", "3:4", "1:2", "1:3"],
                    answer: 1,
                    explain: "Total energy = ½kA². At x = A/2: PE = ½k(A/2)² = ½kA²/4. KE = Total - PE = 3kA²/8. Ratio KE:Total = 3:4"
                },
                {
                    id: "p20",
                    text: "The magnetic flux through a coil changes from 0.1 Wb to 0.3 Wb in 0.2 s. Induced EMF is:",
                    options: ["1 V", "2 V", "0.5 V", "4 V"],
                    answer: 0,
                    explain: "EMF = -dΦ/dt = -(0.3-0.1)/0.2 = -0.2/0.2 = -1 V. Magnitude = 1 V"
                },
                {
                    id: "p21",
                    text: "A uniform rod of mass m and length L is pivoted at distance L/3 from one end. Its moment of inertia about pivot is:",
                    options: ["mL²/9", "mL²/12", "4mL²/27", "2mL²/9"],
                    answer: 2,
                    explain: "I = ICM + md². ICM = mL²/12, d = L/6. I = mL²/12 + m(L/6)² = mL²/12 + mL²/36 = 4mL²/27"
                },
                {
                    id: "p22",
                    text: "In photoelectric effect, maximum kinetic energy depends on:",
                    options: ["Intensity only", "Frequency only", "Both intensity and frequency", "Work function only"],
                    answer: 1,
                    explain: "KEmax = hf - φ depends only on frequency f, independent of intensity"
                },
                {
                    id: "p23",
                    text: "A wire of resistance 10Ω is cut into 10 equal parts and connected in parallel. Net resistance is:",
                    options: ["1Ω", "0.1Ω", "10Ω", "100Ω"],
                    answer: 1,
                    explain: "Each part has resistance 1Ω. Ten 1Ω resistors in parallel: 1/R = 10/1 = 10, so R = 0.1Ω"
                },
                {
                    id: "p24",
                    text: "The impedance of series LR circuit with L = 3H, R = 4Ω at ω = 1 rad/s is:",
                    options: ["5Ω", "7Ω", "3Ω", "4Ω"],
                    answer: 0,
                    explain: "Z = √(R² + (ωL)²) = √(4² + (1×3)²) = √(16 + 9) = 5Ω"
                },
                {
                    id: "p25",
                    text: "A particle in SHM has displacement x = A sin(ωt). At t = T/4, its velocity is:",
                    options: ["0", "ωA", "-ωA", "ωA/2"],
                    answer: 0,
                    explain: "v = dx/dt = ωA cos(ωt). At t = T/4, ωt = π/2, so v = ωA cos(π/2) = 0"
                },
                {
                    id: "p26",
                    text: "The electric potential at distance r from point charge q is V. Electric field at same point is:",
                    options: ["V/r", "V²/r", "Vr", "V/r²"],
                    answer: 0,
                    explain: "V = kq/r and E = kq/r². Therefore E = V/r"
                },
                {
                    id: "p27",
                    text: "A body falls from height h. Time to fall last h/4 distance is:",
                    options: ["t/2", "t(1-√3/2)", "t/4", "t√3/2"],
                    answer: 1,
                    explain: "Total time t = √(2h/g). Time to fall 3h/4 = t√3/2. Last h/4 takes t - t√3/2 = t(1-√3/2)"
                },
                {
                    id: "p28",
                    text: "The threshold frequency for photoelectric effect is f₀. For frequency 2f₀, stopping potential is:",
                    options: ["hf₀/e", "2hf₀/e", "hf₀/(2e)", "3hf₀/e"],
                    answer: 0,
                    explain: "eV₀ = h(2f₀) - hf₀ = hf₀. Therefore V₀ = hf₀/e"
                },
                {
                    id: "p29",
                    text: "A transformer has 100 primary turns and 500 secondary turns. If primary voltage is 220V, secondary voltage is:",
                    options: ["44V", "220V", "1100V", "2750V"],
                    answer: 2,
                    explain: "Vs/Vp = Ns/Np = 500/100 = 5. Therefore Vs = 5 × 220 = 1100V"
                },
                {
                    id: "p30",
                    text: "The magnetic field at center of current-carrying circular loop is B₀. At distance equal to radius on axis, field is:",
                    options: ["B₀/2", "B₀/√2", "B₀/(2√2)", "B₀/4"],
                    answer: 2,
                    explain: "On axis at distance R: B = μ₀IR²/(2(R²+R²)^(3/2)) = μ₀I/(4√2R) = B₀/(2√2)"
                },
                {
                    id: "p31",
                    text: "A gas at 27°C is heated to 327°C at constant pressure. Volume ratio V₂/V₁ is:",
                    options: ["2", "12.1", "1.5", "327/27"],
                    answer: 0,
                    explain: "Charles's law: V₁/T₁ = V₂/T₂. T₁ = 300K, T₂ = 600K. V₂/V₁ = T₂/T₁ = 600/300 = 2"
                },
                {
                    id: "p32",
                    text: "The work function of cesium is 2.1 eV. Threshold wavelength is:",
                    options: ["590 nm", "295 nm", "1180 nm", "495 nm"],
                    answer: 0,
                    explain: "λ₀ = hc/φ = 1240 eV·nm / 2.1 eV = 590 nm"
                },
                {
                    id: "p33",
                    text: "A capacitor C discharges through resistor R. After time t = RC, charge becomes:",
                    options: ["C/e", "CV/e", "Q₀/e", "Q₀e"],
                    answer: 2,
                    explain: "q(t) = q₀e^(-t/RC). At t = RC: q = q₀e^(-1) = q₀/e"
                },
                {
                    id: "p34",
                    text: "Two masses m₁ and m₂ are connected by string over pulley. The tension in string is:",
                    options: ["(m₁ + m₂)g/2", "2m₁m₂g/(m₁ + m₂)", "m₁m₂g/(m₁ + m₂)", "(m₁ - m₂)g"],
                    answer: 1,
                    explain: "For Atwood machine: T = 2m₁m₂g/(m₁ + m₂) assuming m₁ ≠ m₂"
                },
                {
                    id: "p35",
                    text: "A plane mirror rotates through angle θ. The reflected ray rotates through:",
                    options: ["θ", "2θ", "θ/2", "4θ"],
                    answer: 1,
                    explain: "When mirror rotates by θ, normal rotates by θ, reflected ray deviates by 2θ"
                },
                {
                    id: "p36",
                    text: "The binding energy per nucleon is maximum around mass number:",
                    options: ["16", "56", "120", "238"],
                    answer: 1,
                    explain: "Iron-56 has maximum binding energy per nucleon, making it most stable nucleus"
                },
                {
                    id: "p37",
                    text: "A charge q moves in magnetic field B with velocity v. Force experienced is:",
                    options: ["qvB sinθ", "qvB cosθ", "qvB", "qB/v"],
                    answer: 0,
                    explain: "Magnetic force F = qv × B = qvB sinθ where θ is angle between v and B"
                },
                {
                    id: "p38",
                    text: "The distance of closest approach in Rutherford scattering depends on:",
                    options: ["Kinetic energy only", "Charge only", "Both kinetic energy and charge", "Mass only"],
                    answer: 2,
                    explain: "d = kZe²/KE where Z is atomic number, e is elementary charge"
                },
                {
                    id: "p39",
                    text: "A mass attached to spring oscillates with period T. If mass is quadrupled, new period is:",
                    options: ["T", "2T", "4T", "T/2"],
                    answer: 1,
                    explain: "T = 2π√(m/k). When m → 4m: T' = 2π√(4m/k) = 2T"
                },
                {
                    id: "p40",
                    text: "The coefficient of restitution for perfectly inelastic collision is:",
                    options: ["0", "1", "0.5", "-1"],
                    answer: 0,
                    explain: "For perfectly inelastic collision, objects stick together, e = 0"
                },
                {
                    id: "p41",
                    text: "A wire carrying current I is placed in magnetic field B. Force per unit length is:",
                    options: ["BI", "B/I", "I/B", "BI²"],
                    answer: 0,
                    explain: "Force per unit length on current-carrying conductor: F/L = BI (when perpendicular)"
                },
                {
                    id: "p42",
                    text: "The power dissipated in resistor R carrying current I is:",
                    options: ["I²R", "IR", "I²/R", "R/I²"],
                    answer: 0,
                    explain: "Power P = I²R from Joule heating in resistor"
                },
                {
                    id: "p43",
                    text: "A lens has power +5D. Its focal length is:",
                    options: ["0.2 m", "5 m", "0.5 m", "2 m"],
                    answer: 0,
                    explain: "Power P = 1/f. For P = +5D: f = 1/5 = 0.2 m"
                },
                {
                    id: "p44",
                    text: "The energy density of electromagnetic wave is proportional to:",
                    options: ["E", "E²", "E³", "√E"],
                    answer: 1,
                    explain: "Energy density u = ½ε₀E² + B²/(2μ₀) ∝ E² in EM wave"
                },
                {
                    id: "p45",
                    text: "A radioactive sample has activity 1000 Bq. After 2 half-lives, activity becomes:",
                    options: ["500 Bq", "250 Bq", "125 Bq", "750 Bq"],
                    answer: 1,
                    explain: "After n half-lives, activity = A₀/2ⁿ. After 2 half-lives: 1000/4 = 250 Bq"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which of the following has maximum number of lone pairs on central atom?",
                    options: ["NH₃", "H₂O", "XeF₄", "ClF₃"],
                    answer: 2,
                    explain: "XeF₄ has 2 lone pairs on Xe atom in square planar geometry, maximum among given options."
                },
                {
                    id: "c2",
                    text: "The correct order of second ionization energy is:",
                    options: ["Li > Be > B > C", "C > B > Be > Li", "Be > C > B > Li", "Li > C > B > Be"],
                    answer: 2,
                    explain: "Second IE order: Be (removing from 1s) > C > B > Li. Be has highest due to removing electron from inner shell."
                },
                {
                    id: "c3",
                    text: "Which of the following is not a π-acid ligand?",
                    options: ["CO", "NO", "NH₃", "PF₃"],
                    answer: 2,
                    explain: "NH₃ is σ-donor only. CO, NO, PF₃ are π-acid ligands that accept electron density into empty π* orbitals."
                },
                {
                    id: "c4",
                    text: "The number of stereoisomers of [CoCl₂(en)₂]⁺ is:",
                    options: ["2", "3", "4", "6"],
                    answer: 0,
                    explain: "This complex has 2 stereoisomers: cis and trans forms due to arrangement of Cl ligands."
                },
                {
                    id: "c5",
                    text: "Which alcohol shows fastest rate of dehydration with conc. H₂SO₄?",
                    options: ["Primary", "Secondary", "Tertiary", "All equal"],
                    answer: 2,
                    explain: "Tertiary alcohols form most stable carbocations, hence fastest dehydration via E1 mechanism."
                },
                {
                    id: "c6",
                    text: "The hybridization of phosphorus in PF₅ is:",
                    options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
                    answer: 1,
                    explain: "PF₅ has 5 bonding pairs around P requiring sp³d hybridization (trigonal bipyramidal)."
                },
                {
                    id: "c7",
                    text: "Which of the following is most stable alkene?",
                    options: ["1-butene", "2-butene", "2-methylpropene", "Cyclobutene"],
                    answer: 2,
                    explain: "2-methylpropene is most substituted alkene, hence most stable due to hyperconjugation."
                },
                {
                    id: "c8",
                    text: "The correct order of acidic strength is:",
                    options: ["HClO₄ > HClO₃ > HClO₂ > HClO", "HClO > HClO₂ > HClO₃ > HClO₄", "HClO₃ > HClO₄ > HClO₂ > HClO", "HClO₂ > HClO₃ > HClO₄ > HClO"],
                    answer: 0,
                    explain: "Acidic strength increases with number of oxygen atoms: HClO₄ > HClO₃ > HClO₂ > HClO."
                },
                {
                    id: "c9",
                    text: "Which of the following shows maximum covalent character?",
                    options: ["LiCl", "BeF₂", "BF₃", "CF₄"],
                    answer: 3,
                    explain: "CF₄ has maximum covalent character due to high electronegativity of both C and F."
                },
                {
                    id: "c10",
                    text: "The number of unpaired electrons in Mn³⁺ ion is:",
                    options: ["2", "3", "4", "5"],
                    answer: 2,
                    explain: "Mn³⁺ has d⁴ configuration. In octahedral field (high spin): t₂g³ eg¹ = 4 unpaired electrons."
                },
                {
                    id: "c11",
                    text: "Which of the following undergoes Cannizzaro reaction?",
                    options: ["Acetaldehyde", "Formaldehyde", "Propionaldehyde", "Butyraldehyde"],
                    answer: 1,
                    explain: "Formaldehyde lacks α-hydrogen, undergoes Cannizzaro reaction (disproportionation)."
                },
                {
                    id: "c12",
                    text: "The bond angle in water molecule is approximately:",
                    options: ["109.5°", "107°", "105°", "104.5°"],
                    answer: 3,
                    explain: "H₂O has bent structure with bond angle ~104.5° due to two lone pairs on oxygen."
                },
                {
                    id: "c13",
                    text: "Which of the following is strongest reducing agent?",
                    options: ["Li", "Na", "K", "Cs"],
                    answer: 0,
                    explain: "Li has most negative standard reduction potential, hence strongest reducing agent."
                },
                {
                    id: "c14",
                    text: "The correct order of stability of halides is:",
                    options: ["MF > MCl > MBr > MI", "MI > MBr > MCl > MF", "MCl > MF > MBr > MI", "MBr > MI > MCl > MF"],
                    answer: 0,
                    explain: "Fluorides are most stable due to high lattice energy and strong ionic bonding."
                },
                {
                    id: "c15",
                    text: "Which of the following has highest electron affinity?",
                    options: ["F", "Cl", "Br", "I"],
                    answer: 1,
                    explain: "Cl has highest electron affinity due to optimal size. F has lower EA due to electron-electron repulsion."
                },
                {
                    id: "c16",
                    text: "The number of σ and π bonds in acetylene (C₂H₂) are:",
                    options: ["3σ, 2π", "2σ, 3π", "5σ, 0π", "4σ, 1π"],
                    answer: 0,
                    explain: "C₂H₂ has 3 σ bonds (2 C-H + 1 C-C) and 2 π bonds in C≡C triple bond."
                },
                {
                    id: "c17",
                    text: "Which of the following is not aromatic?",
                    options: ["Benzene", "Pyridine", "Pyrrole", "Cyclooctatetraene"],
                    answer: 3,
                    explain: "Cyclooctatetraene has 8π electrons, violates Hückel's 4n+2 rule, hence not aromatic."
                },
                {
                    id: "c18",
                    text: "The hybridization of carbon atoms in graphite is:",
                    options: ["sp", "sp²", "sp³", "sp³d"],
                    answer: 1,
                    explain: "Carbon atoms in graphite are sp² hybridized forming planar hexagonal layers."
                },
                {
                    id: "c19",
                    text: "Which of the following shows maximum hydrogen bonding?",
                    options: ["HF", "H₂O", "NH₃", "H₂S"],
                    answer: 0,
                    explain: "HF shows strongest hydrogen bonding due to high electronegativity of fluorine."
                },
                {
                    id: "c20",
                    text: "The coordination number of Na⁺ in rock salt structure is:",
                    options: ["4", "6", "8", "12"],
                    answer: 1,
                    explain: "In NaCl (rock salt) structure, Na⁺ is surrounded by 6 Cl⁻ ions octahedrally."
                },
                {
                    id: "c21",
                    text: "Which of the following is most basic?",
                    options: ["CH₃NH₂", "(CH₃)₂NH", "(CH₃)₃N", "NH₃"],
                    answer: 1,
                    explain: "Secondary amine (CH₃)₂NH is most basic due to optimal +I effect and steric factors."
                },
                {
                    id: "c22",
                    text: "The oxidation state of Fe in K₄[Fe(CN)₆] is:",
                    options: ["+2", "+3", "+4", "+6"],
                    answer: 0,
                    explain: "In K₄[Fe(CN)₆]: 4(+1) + Fe + 6(-1) = 0, so Fe = +2."
                },
                {
                    id: "c23",
                    text: "Which of the following shows maximum lattice energy?",
                    options: ["NaF", "NaCl", "MgO", "CaO"],
                    answer: 2,
                    explain: "MgO has highest lattice energy due to +2,-2 charges and small ionic radii."
                },
                {
                    id: "c24",
                    text: "The shape of ammonia molecule is:",
                    options: ["Tetrahedral", "Trigonal planar", "Pyramidal", "Linear"],
                    answer: 2,
                    explain: "NH₃ has pyramidal shape due to one lone pair on nitrogen (VSEPR theory)."
                },
                {
                    id: "c25",
                    text: "Which of the following is strongest oxidizing agent?",
                    options: ["F₂", "Cl₂", "Br₂", "I₂"],
                    answer: 0,
                    explain: "F₂ has highest reduction potential, making it strongest oxidizing agent."
                },
                {
                    id: "c26",
                    text: "The number of atoms per unit cell in face-centered cubic structure is:",
                    options: ["1", "2", "4", "8"],
                    answer: 2,
                    explain: "FCC has 8 corner atoms (8×1/8) + 6 face atoms (6×1/2) = 4 atoms per unit cell."
                },
                {
                    id: "c27",
                    text: "Which of the following undergoes nucleophilic substitution fastest?",
                    options: ["CH₃Cl", "CH₃Br", "CH₃I", "CH₃F"],
                    answer: 2,
                    explain: "CH₃I undergoes SN2 fastest due to best leaving group (I⁻) among halogens."
                },
                {
                    id: "c28",
                    text: "The correct order of thermal stability of carbonates is:",
                    options: ["BeCO₃ > MgCO₃ > CaCO₃", "CaCO₃ > MgCO₃ > BeCO₃", "MgCO₃ > BeCO₃ > CaCO₃", "All equal"],
                    answer: 1,
                    explain: "Thermal stability of carbonates increases down the group: CaCO₃ > MgCO₃ > BeCO₃."
                },
                {
                    id: "c29",
                    text: "Which of the following is paramagnetic?",
                    options: ["N₂", "O₂", "F₂", "Ne₂"],
                    answer: 1,
                    explain: "O₂ has two unpaired electrons in π* orbitals, making it paramagnetic."
                },
                {
                    id: "c30",
                    text: "The number of lone pairs in XeF₂ is:",
                    options: ["1", "2", "3", "4"],
                    answer: 2,
                    explain: "Xe has 8 valence electrons, uses 2 for bonding with F atoms, leaving 6 electrons = 3 lone pairs."
                },
                {
                    id: "c31",
                    text: "Which of the following shows maximum dipole moment?",
                    options: ["BF₃", "NH₃", "CO₂", "CH₄"],
                    answer: 1,
                    explain: "NH₃ has pyramidal structure with lone pair, creating significant dipole moment."
                },
                {
                    id: "c32",
                    text: "The bond order in CO molecule is:",
                    options: ["2", "2.5", "3", "3.5"],
                    answer: 2,
                    explain: "CO has triple bond character with bond order 3 from molecular orbital theory."
                },
                {
                    id: "c33",
                    text: "Which of the following is most stable carbocation?",
                    options: ["CH₃⁺", "(CH₃)₂CH⁺", "(CH₃)₃C⁺", "C₆H₅CH₂⁺"],
                    answer: 3,
                    explain: "Benzyl cation is most stable due to resonance with benzene ring."
                },
                {
                    id: "c34",
                    text: "The hybridization of central atom in SF₄ is:",
                    options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
                    answer: 1,
                    explain: "SF₄ has 5 electron pairs (4 bonding + 1 lone pair) requiring sp³d hybridization."
                },
                {
                    id: "c35",
                    text: "Which of the following has maximum boiling point?",
                    options: ["He", "Ne", "Ar", "Kr"],
                    answer: 3,
                    explain: "Kr has highest atomic mass and strongest van der Waals forces, highest boiling point."
                },
                {
                    id: "c36",
                    text: "The correct order of ionic radii is:",
                    options: ["Li⁺ > Na⁺ > K⁺", "K⁺ > Na⁺ > Li⁺", "Na⁺ > K⁺ > Li⁺", "All equal"],
                    answer: 1,
                    explain: "Ionic radii increase down the group: K⁺ > Na⁺ > Li⁺ due to additional electron shells."
                },
                {
                    id: "c37",
                    text: "Which of the following undergoes electrophilic substitution most readily?",
                    options: ["Benzene", "Toluene", "Nitrobenzene", "Chlorobenzene"],
                    answer: 1,
                    explain: "Toluene has electron-donating methyl group that activates benzene ring."
                },
                {
                    id: "c38",
                    text: "The number of geometrical isomers of [Ma₂b₂] (square planar) is:",
                    options: ["2", "3", "4", "5"],
                    answer: 0,
                    explain: "Square planar [Ma₂b₂] has 2 geometrical isomers: cis and trans."
                },
                {
                    id: "c39",
                    text: "Which of the following is strongest acid?",
                    options: ["CH₃COOH", "CHCl₂COOH", "CCl₃COOH", "CF₃COOH"],
                    answer: 3,
                    explain: "CF₃COOH is strongest due to maximum -I effect of three fluorine atoms."
                },
                {
                    id: "c40",
                    text: "The magnetic moment of Cu²⁺ ion is:",
                    options: ["1.73 BM", "2.83 BM", "3.87 BM", "0 BM"],
                    answer: 0,
                    explain: "Cu²⁺ has d⁹ configuration with 1 unpaired electron. μ = √[n(n+2)] = √3 = 1.73 BM."
                },
                {
                    id: "c41",
                    text: "Which of the following exhibits optical isomerism?",
                    options: ["CH₄", "CHCl₃", "CHClBrI", "C₂H₆"],
                    answer: 2,
                    explain: "CHClBrI has 4 different groups on carbon, making it chiral and optically active."
                },
                {
                    id: "c42",
                    text: "The rate constant units for second order reaction are:",
                    options: ["s⁻¹", "mol⁻¹ L s⁻¹", "mol L⁻¹ s⁻¹", "mol² L⁻² s⁻¹"],
                    answer: 1,
                    explain: "Second order rate constant has units mol⁻¹ L s⁻¹ or M⁻¹ s⁻¹."
                },
                {
                    id: "c43",
                    text: "Which of the following has maximum melting point?",
                    options: ["Diamond", "Graphite", "Fullerene", "Graphene"],
                    answer: 0,
                    explain: "Diamond has highest melting point due to strong 3D network of C-C covalent bonds."
                },
                {
                    id: "c44",
                    text: "The number of bridging CO ligands in [Fe₂(CO)₉] is:",
                    options: ["0", "2", "3", "6"],
                    answer: 2,
                    explain: "Fe₂(CO)₉ has 3 bridging CO groups connecting two Fe atoms."
                },
                {
                    id: "c45",
                    text: "Which of the following is Lewis acid?",
                    options: ["NH₃", "BF₃", "H₂O", "OH⁻"],
                    answer: 1,
                    explain: "BF₃ has empty p-orbital and can accept electron pair, making it Lewis acid."
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "Which of the following is the correct sequence of Calvin cycle?",
                    options: ["Carboxylation → Reduction → Regeneration", "Reduction → Carboxylation → Regeneration", "Regeneration → Carboxylation → Reduction", "Carboxylation → Regeneration → Reduction"],
                    answer: 0,
                    explain: "Calvin cycle has three phases: carboxylation (CO₂ fixation), reduction (NADPH use), and regeneration (RuBP formation)."
                },
                {
                    id: "b2",
                    text: "Which tissue is responsible for transport of organic food in plants?",
                    options: ["Xylem", "Phloem", "Cambium", "Epidermis"],
                    answer: 1,
                    explain: "Phloem transports organic nutrients, mainly sucrose, from source to sink tissues."
                },
                {
                    id: "b3",
                    text: "The maximum number of electrons that can be excited in photosystem II is:",
                    options: ["1", "2", "4", "8"],
                    answer: 1,
                    explain: "PSII reaction center P680 can lose maximum 2 electrons which are replaced by water splitting."
                },
                {
                    id: "b4",
                    text: "Which hormone is known as stress hormone in plants?",
                    options: ["Auxin", "Gibberellin", "Abscisic acid", "Cytokinin"],
                    answer: 2,
                    explain: "Abscisic acid helps plants cope with stress conditions like drought, salinity, and cold."
                },
                {
                    id: "b5",
                    text: "The arrangement of sepals and petals in flower bud is called:",
                    options: ["Phyllotaxy", "Vernation", "Aestivation", "Placentation"],
                    answer: 2,
                    explain: "Aestivation refers to arrangement of sepals and petals in flower bud before opening."
                },
                {
                    id: "b6",
                    text: "Which of the following is an example of aggregate fruit?",
                    options: ["Orange", "Apple", "Blackberry", "Pineapple"],
                    answer: 2,
                    explain: "Blackberry is aggregate fruit formed from multiple carpels of single flower."
                },
                {
                    id: "b7",
                    text: "The water potential of pure water at standard conditions is:",
                    options: ["Zero", "Positive", "Negative", "Infinite"],
                    answer: 0,
                    explain: "Pure water has water potential of zero, used as reference for measuring water potential."
                },
                {
                    id: "b8",
                    text: "Which of the following is involved in opening of stomata?",
                    options: ["Decrease in turgor pressure", "Increase in CO₂", "Accumulation of K⁺ ions", "Darkness"],
                    answer: 2,
                    explain: "K⁺ ion accumulation in guard cells increases turgor pressure, causing stomatal opening."
                },
                {
                    id: "b9",
                    text: "The enzyme that joins CO₂ to RuBP in Calvin cycle is:",
                    options: ["PEP carboxylase", "RuBisCO", "Carbonic anhydrase", "Pyruvate kinase"],
                    answer: 1,
                    explain: "RuBisCO (Ribulose bisphosphate carboxylase oxygenase) catalyzes CO₂ fixation in Calvin cycle."
                },
                {
                    id: "b10",
                    text: "Which of the following is characteristic of C₄ plants?",
                    options: ["Low CO₂ compensation point", "High photorespiration", "Mesophyll cell specialization", "All of above"],
                    answer: 0,
                    explain: "C₄ plants have low CO₂ compensation point due to efficient CO₂ concentration mechanism."
                },
                {
                    id: "b11",
                    text: "The protective tissue in woody stems is:",
                    options: ["Epidermis", "Cork", "Collenchyma", "Sclerenchyma"],
                    answer: 1,
                    explain: "Cork tissue replaces epidermis in woody stems, providing protection against water loss."
                },
                {
                    id: "b12",
                    text: "Which of the following shows circinate vernation?",
                    options: ["Palm leaves", "Fern fronds", "Grass leaves", "Rose leaves"],
                    answer: 1,
                    explain: "Fern fronds show circinate vernation where young leaves are coiled like watch spring."
                },
                {
                    id: "b13",
                    text: "The term 'polyembryony' refers to:",
                    options: ["Multiple ovules", "Multiple embryos in single seed", "Multiple seeds", "Multiple cotyledons"],
                    answer: 1,
                    explain: "Polyembryony is occurrence of more than one embryo in single seed, common in citrus fruits."
                },
                {
                    id: "b14",
                    text: "Which of the following is involved in nitrogen fixation?",
                    options: ["Nitrogenase", "Nitrate reductase", "Glutamine synthetase", "Urease"],
                    answer: 0,
                    explain: "Nitrogenase enzyme complex converts atmospheric nitrogen to ammonia in biological nitrogen fixation."
                },
                {
                    id: "b15",
                    text: "The region of root where most water absorption occurs is:",
                    options: ["Root cap", "Meristematic zone", "Root hair zone", "Maturation zone"],
                    answer: 2,
                    explain: "Root hair zone has maximum surface area for water and mineral absorption."
                },
                {
                    id: "b16",
                    text: "Which of the following is a quantitative short day plant?",
                    options: ["Rice", "Wheat", "Spinach", "Cucumber"],
                    answer: 0,
                    explain: "Rice is facultative short day plant that flowers faster with short days but can flower with long days."
                },
                {
                    id: "b17",
                    text: "The phenomenon of seed formation without fertilization is:",
                    options: ["Parthenocarpy", "Apomixis", "Parthenogenesis", "Apogamy"],
                    answer: 1,
                    explain: "Apomixis is asexual reproduction through seed formation without fertilization."
                },
                {
                    id: "b18",
                    text: "Which of the following is essential for photophosphorylation?",
                    options: ["CO₂", "NADP⁺", "ADP + Pi", "RuBP"],
                    answer: 2,
                    explain: "Photophosphorylation requires ADP and inorganic phosphate to synthesize ATP."
                },
                {
                    id: "b19",
                    text: "The vascular bundles in monocot stems are:",
                    options: ["Open and arranged in ring", "Closed and scattered", "Open and scattered", "Closed and arranged in ring"],
                    answer: 1,
                    explain: "Monocot stems have closed vascular bundles (no cambium) scattered in ground tissue."
                },
                {
                    id: "b20",
                    text: "Which of the following is involved in gibberellin biosynthesis?",
                    options: ["Tryptophan", "Tyrosine", "Mevalonic acid", "Shikimic acid"],
                    answer: 2,
                    explain: "Gibberellins are synthesized from mevalonic acid pathway like other terpenoids."
                },
                {
                    id: "b21",
                    text: "The loss of water from aerial parts of plant is called:",
                    options: ["Transpiration", "Guttation", "Bleeding", "Exudation"],
                    answer: 0,
                    explain: "Transpiration is loss of water vapor from aerial parts, mainly through stomata."
                },
                {
                    id: "b22",
                    text: "Which of the following is characteristic of wind-pollinated flowers?",
                    options: ["Bright colors", "Strong fragrance", "Feathery stigma", "Nectar production"],
                    answer: 2,
                    explain: "Wind-pollinated flowers have feathery stigma to capture airborne pollen effectively."
                },
                {
                    id: "b23",
                    text: "The primary acceptor of CO₂ in C₃ plants is:",
                    options: ["PEP", "RuBP", "Oxaloacetic acid", "Pyruvic acid"],
                    answer: 1,
                    explain: "RuBP (Ribulose bisphosphate) is primary CO₂ acceptor in C₃ plants during Calvin cycle."
                },
                {
                    id: "b24",
                    text: "Which of the following shows thigmotropism?",
                    options: ["Root tips", "Shoot tips", "Tendrils", "Leaves"],
                    answer: 2,
                    explain: "Tendrils show thigmotropism (growth response to touch) by coiling around support."
                },
                {
                    id: "b25",
                    text: "The site of protein synthesis in chloroplasts is:",
                    options: ["Thylakoid", "Stroma", "Inner membrane", "Outer membrane"],
                    answer: 1,
                    explain: "Chloroplast stroma contains 70S ribosomes that synthesize chloroplast proteins."
                },
                {
                    id: "b26",
                    text: "Which of the following is involved in photoperiodism?",
                    options: ["Cryptochrome", "Phytochrome", "Phototropin", "Rhodopsin"],
                    answer: 1,
                    explain: "Phytochrome system detects day length and mediates photoperiodic flowering responses."
                },
                {
                    id: "b27",
                    text: "The cambial ring in dicot stem is formed by:",
                    options: ["Vascular cambium only", "Cork cambium only", "Both vascular and cork cambium", "Fascicular and interfascicular cambium"],
                    answer: 3,
                    explain: "Cambial ring is formed by fascicular cambium in bundles and interfascicular cambium between bundles."
                },
                {
                    id: "b28",
                    text: "Which of the following is major site of photosynthesis?",
                    options: ["Palisade parenchyma", "Spongy parenchyma", "Bundle sheath", "Epidermis"],
                    answer: 0,
                    explain: "Palisade parenchyma has densely packed chloroplasts and is major site of photosynthesis."
                },
                {
                    id: "b29",
                    text: "The phenomenon where some plants require cold treatment for flowering is:",
                    options: ["Photoperiodism", "Vernalization", "Stratification", "Etiolation"],
                    answer: 1,
                    explain: "Vernalization is cold treatment required by some plants to induce flowering."
                },
                {
                    id: "b30",
                    text: "Which of the following is example of lateral meristem?",
                    options: ["Apical meristem", "Intercalary meristem", "Vascular cambium", "Primary meristem"],
                    answer: 2,
                    explain: "Vascular cambium is lateral meristem responsible for secondary growth in woody plants."
                },
                {
                    id: "b31",
                    text: "The endosperm in coconut is:",
                    options: ["Solid", "Liquid", "Both solid and liquid", "Absent"],
                    answer: 2,
                    explain: "Coconut has both solid endosperm (coconut meat) and liquid endosperm (coconut water)."
                },
                {
                    id: "b32",
                    text: "Which of the following is essential for seed germination?",
                    options: ["Light", "Oxygen", "Soil", "Fertilizer"],
                    answer: 1,
                    explain: "Oxygen is essential for aerobic respiration during seed germination process."
                },
                {
                    id: "b33",
                    text: "The pressure flow hypothesis explains:",
                    options: ["Water transport", "Mineral transport", "Sugar transport", "Hormone transport"],
                    answer: 2,
                    explain: "Pressure flow hypothesis explains translocation of sugars in phloem from source to sink."
                },
                {
                    id: "b34",
                    text: "Which of the following is involved in cell wall synthesis?",
                    options: ["Ribosomes", "Golgi apparatus", "Mitochondria", "Chloroplasts"],
                    answer: 1,
                    explain: "Golgi apparatus modifies and packages cellulose and other cell wall components."
                },
                {
                    id: "b35",
                    text: "The opening and closing of flowers in response to light is called:",
                    options: ["Phototropism", "Photoperiodism", "Photonasty", "Phototaxis"],
                    answer: 2,
                    explain: "Photonastic movements are opening/closing responses of flowers to light intensity changes."
                },
                {
                    id: "b36",
                    text: "Which of the following is stored in aleurone layer of seeds?",
                    options: ["Starch", "Proteins", "Lipids", "Cellulose"],
                    answer: 1,
                    explain: "Aleurone layer in cereal seeds stores proteins and enzymes needed for germination."
                },
                {
                    id: "b37",
                    text: "The term 'bolting' refers to:",
                    options: ["Rapid flowering", "Rapid stem elongation", "Rapid root growth", "Rapid leaf growth"],
                    answer: 1,
                    explain: "Bolting is rapid stem elongation that typically occurs before flowering, often induced by gibberellins."
                },
                {
                    id: "b38",
                    text: "Which of the following is involved in abscission?",
                    options: ["Auxin", "Cytokinin", "Ethylene", "Gibberellin"],
                    answer: 2,
                    explain: "Ethylene promotes abscission by activating hydrolytic enzymes in abscission zone."
                },
                {
                    id: "b39",
                    text: "The primary function of root hairs is:",
                    options: ["Protection", "Support", "Absorption", "Photosynthesis"],
                    answer: 2,
                    explain: "Root hairs increase surface area for water and mineral absorption from soil."
                },
                {
                    id: "b40",
                    text: "Which of the following shows epigeal germination?",
                    options: ["Pea", "Bean", "Maize", "Wheat"],
                    answer: 1,
                    explain: "Bean shows epigeal germination where cotyledons emerge above ground and become photosynthetic."
                },
                {
                    id: "b41",
                    text: "The region of stem from which leaves arise is called:",
                    options: ["Internode", "Node", "Axil", "Apex"],
                    answer: 1,
                    explain: "Node is region of stem where leaves are attached and axillary buds are present."
                },
                {
                    id: "b42",
                    text: "Which of the following is involved in stomatal regulation?",
                    options: ["Light intensity", "CO₂ concentration", "Water availability", "All of above"],
                    answer: 3,
                    explain: "Stomatal opening/closing is regulated by light, CO₂ levels, water status, and temperature."
                },
                {
                    id: "b43",
                    text: "The waxy coating on leaf surface is called:",
                    options: ["Cutin", "Suberin", "Lignin", "Pectin"],
                    answer: 0,
                    explain: "Cutin forms waxy cuticle on leaf epidermis to prevent water loss."
                },
                {
                    id: "b44",
                    text: "Which of the following is found in sieve tube elements?",
                    options: ["Nucleus", "P-protein", "Ribosomes", "Vacuole"],
                    answer: 1,
                    explain: "Sieve tube elements contain P-protein (phloem protein) but lack nucleus and ribosomes."
                },
                {
                    id: "b45",
                    text: "The arrangement of leaves in which each leaf is at different height is:",
                    options: ["Opposite", "Whorled", "Alternate", "Clustered"],
                    answer: 2,
                    explain: "Alternate phyllotaxy has single leaf at each node at different heights on stem."
                },

                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Which of the following is the largest cell in human body?",
                    options: ["Neuron", "Ovum", "Muscle cell", "Hepatocyte"],
                    answer: 1,
                    explain: "Human ovum (egg cell) is largest cell with diameter of about 0.1 mm."
                },
                {
                    id: "z2",
                    text: "The process by which glucose is converted to lactate in absence of oxygen is:",
                    options: ["Aerobic respiration", "Anaerobic respiration", "Fermentation", "Glycolysis"],
                    answer: 2,
                    explain: "Lactic acid fermentation converts glucose to lactate without oxygen in muscle cells."
                },
                {
                    id: "z3",
                    text: "Which hormone regulates sleep-wake cycle?",
                    options: ["Cortisol", "Melatonin", "Thyroxine", "Insulin"],
                    answer: 1,
                    explain: "Melatonin from pineal gland regulates circadian rhythms and sleep-wake cycles."
                },
                {
                    id: "z4",
                    text: "The gap between two neurons is called:",
                    options: ["Node of Ranvier", "Synapse", "Neuromuscular junction", "Axon hillock"],
                    answer: 1,
                    explain: "Synapse is junction between two neurons where neurotransmitters are released."
                },
                {
                    id: "z5",
                    text: "Which of the following is not a function of liver?",
                    options: ["Protein synthesis", "Glycogen storage", "Bile production", "Insulin production"],
                    answer: 3,
                    explain: "Insulin is produced by pancreatic beta cells, not liver. Liver performs other listed functions."
                },
                {
                    id: "z6",
                    text: "The functional unit of muscle contraction is:",
                    options: ["Myofibril", "Sarcomere", "Actin filament", "Myosin filament"],
                    answer: 1,
                    explain: "Sarcomere is basic contractile unit between two Z-lines containing actin and myosin."
                },
                {
                    id: "z7",
                    text: "Which part of nephron is involved in counter-current mechanism?",
                    options: ["Proximal tubule", "Loop of Henle", "Distal tubule", "Collecting duct"],
                    answer: 1,
                    explain: "Loop of Henle creates counter-current flow essential for urine concentration."
                },
                {
                    id: "z8",
                    text: "The hormone that stimulates uterine contractions is:",
                    options: ["Progesterone", "Estrogen", "Oxytocin", "Prolactin"],
                    answer: 2,
                    explain: "Oxytocin from posterior pituitary stimulates uterine contractions during labor."
                },
                {
                    id: "z9",
                    text: "Which of the following is involved in adaptive immunity?",
                    options: ["Neutrophils", "Lymphocytes", "Macrophages", "Eosinophils"],
                    answer: 1,
                    explain: "B and T lymphocytes provide specific adaptive immune responses."
                },
                {
                    id: "z10",
                    text: "The pH of normal human blood is approximately:",
                    options: ["6.8", "7.0", "7.4", "8.0"],
                    answer: 2,
                    explain: "Normal blood pH is maintained around 7.4 (slightly alkaline) by buffer systems."
                },
                {
                    id: "z11",
                    text: "Which enzyme breaks down starch in mouth?",
                    options: ["Pepsin", "Amylase", "Lipase", "Trypsin"],
                    answer: 1,
                    explain: "Salivary amylase (ptyalin) begins starch digestion in oral cavity."
                },
                {
                    id: "z12",
                    text: "The respiratory center is located in:",
                    options: ["Cerebrum", "Cerebellum", "Medulla oblongata", "Spinal cord"],
                    answer: 2,
                    explain: "Medulla oblongata contains respiratory control center that regulates breathing rhythm."
                },
                {
                    id: "z13",
                    text: "Which of the following is characteristic of cardiac muscle?",
                    options: ["Voluntary control", "Multinucleated", "Intercalated discs", "Non-striated"],
                    answer: 2,
                    explain: "Cardiac muscle has intercalated discs with gap junctions for synchronized contraction."
                },
                {
                    id: "z14",
                    text: "The hormone that increases metabolic rate is:",
                    options: ["Insulin", "Glucagon", "Thyroxine", "Cortisol"],
                    answer: 2,
                    explain: "Thyroxine (T₄) and T₃ from thyroid increase basal metabolic rate."
                },
                {
                    id: "z15",
                    text: "Which part of ear detects sound vibrations?",
                    options: ["Cochlea", "Vestibule", "Semicircular canals", "Eustachian tube"],
                    answer: 0,
                    explain: "Cochlea contains organ of Corti with hair cells that detect sound vibrations."
                },
                {
                    id: "z16",
                    text: "The yellowish pigment in urine is:",
                    options: ["Bilirubin", "Urochrome", "Hemoglobin", "Creatinine"],
                    answer: 1,
                    explain: "Urochrome gives urine its characteristic yellow color, derived from hemoglobin breakdown."
                },
                {
                    id: "z17",
                    text: "Which of the following crosses placenta?",
                    options: ["Antibodies", "Nutrients", "Oxygen", "All of above"],
                    answer: 3,
                    explain: "Placenta allows passage of antibodies (IgG), nutrients, oxygen, and waste products."
                },
                {
                    id: "z18",
                    text: "The period of ventricular relaxation is called:",
                    options: ["Systole", "Diastole", "Cardiac cycle", "Bradycardia"],
                    answer: 1,
                    explain: "Diastole is relaxation phase when ventricles fill with blood."
                },
                {
                    id: "z19",
                    text: "Which vitamin deficiency causes rickets?",
                    options: ["Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D"],
                    answer: 3,
                    explain: "Vitamin D deficiency impairs calcium absorption, causing rickets in children."
                },
                {
                    id: "z20",
                    text: "The cells that produce myelin in central nervous system are:",
                    options: ["Schwann cells", "Oligodendrocytes", "Astrocytes", "Microglia"],
                    answer: 1,
                    explain: "Oligodendrocytes produce myelin sheaths around axons in CNS."
                },
                {
                    id: "z21",
                    text: "Which of the following is not a component of blood?",
                    options: ["Plasma", "RBCs", "Platelets", "Lymph"],
                    answer: 3,
                    explain: "Lymph is tissue fluid, not blood component. Blood contains plasma, RBCs, WBCs, and platelets."
                },
                {
                    id: "z22",
                    text: "The hormone that regulates water reabsorption in kidneys is:",
                    options: ["Aldosterone", "ADH", "Renin", "Angiotensin"],
                    answer: 1,
                    explain: "ADH (antidiuretic hormone) increases water reabsorption in collecting duct."
                },
                {
                    id: "z23",
                    text: "Which of the following has highest oxygen carrying capacity?",
                    options: ["Hemoglobin", "Myoglobin", "Plasma", "Cytochrome"],
                    answer: 1,
                    explain: "Myoglobin has higher oxygen affinity than hemoglobin, stores oxygen in muscles."
                },
                {
                    id: "z24",
                    text: "The structure that prevents backflow in veins is:",
                    options: ["Muscle layer", "Valves", "Elastic fibers", "Smooth muscle"],
                    answer: 1,
                    explain: "Veins have one-way valves to prevent backflow and maintain blood flow toward heart."
                },
                {
                    id: "z25",
                    text: "Which of the following is involved in blood clotting?",
                    options: ["Fibrinogen", "Prothrombin", "Platelets", "All of above"],
                    answer: 3,
                    explain: "Blood clotting involves platelets aggregation, fibrinogen to fibrin conversion, and prothrombin activation."
                },
                {
                    id: "z26",
                    text: "The minimum amount of urine that must be produced daily is:",
                    options: ["500 mL", "1000 mL", "1500 mL", "2000 mL"],
                    answer: 0,
                    explain: "Minimum obligatory urine volume is ~500 mL/day to eliminate metabolic waste."
                },
                {
                    id: "z27",
                    text: "Which part of brain coordinates voluntary movements?",
                    options: ["Cerebrum", "Cerebellum", "Medulla", "Hypothalamus"],
                    answer: 1,
                    explain: "Cerebellum coordinates voluntary movements, maintains balance and motor learning."
                },
                {
                    id: "z28",
                    text: "The hormone that stimulates red blood cell production is:",
                    options: ["Insulin", "Growth hormone", "Erythropoietin", "Thyroxine"],
                    answer: 2,
                    explain: "Erythropoietin from kidneys stimulates RBC production in bone marrow."
                },
                {
                    id: "z29",
                    text: "Which of the following is essential amino acid?",
                    options: ["Glycine", "Alanine", "Lysine", "Serine"],
                    answer: 2,
                    explain: "Lysine is essential amino acid that must be obtained from diet as body cannot synthesize it."
                },
                {
                    id: "z30",
                    text: "The process of removing metabolic waste from blood is called:",
                    options: ["Excretion", "Secretion", "Absorption", "Filtration"],
                    answer: 0,
                    explain: "Excretion is elimination of metabolic wastes from body, primarily by kidneys."
                },
                {
                    id: "z31",
                    text: "Which hormone is produced during pregnancy?",
                    options: ["FSH", "LH", "hCG", "GnRH"],
                    answer: 2,
                    explain: "Human chorionic gonadotropin (hCG) is produced by placenta during pregnancy."
                },
                {
                    id: "z32",
                    text: "The site of gas exchange in lungs is:",
                    options: ["Bronchi", "Bronchioles", "Alveoli", "Pleura"],
                    answer: 2,
                    explain: "Alveoli are tiny air sacs where oxygen and carbon dioxide exchange occurs."
                },
                {
                    id: "z33",
                    text: "Which of the following regulates blood calcium levels?",
                    options: ["Insulin", "Parathyroid hormone", "Growth hormone", "Thyroxine"],
                    answer: 1,
                    explain: "Parathyroid hormone increases blood calcium by promoting bone resorption."
                },
                {
                    id: "z34",
                    text: "The longest phase of cardiac cycle is:",
                    options: ["Atrial systole", "Ventricular systole", "Ventricular diastole", "Isovolumic contraction"],
                    answer: 2,
                    explain: "Ventricular diastole (relaxation and filling) is longest phase of cardiac cycle."
                },
                {
                    id: "z35",
                    text: "Which of the following is not a granulocyte?",
                    options: ["Neutrophil", "Eosinophil", "Basophil", "Monocyte"],
                    answer: 3,
                    explain: "Monocyte is agranulocyte lacking visible cytoplasmic granules."
                },
                {
                    id: "z36",
                    text: "The normal range of human body temperature is:",
                    options: ["35-36°C", "36-37°C", "37-38°C", "38-39°C"],
                    answer: 1,
                    explain: "Normal human body temperature ranges from 36-37°C with average of 37°C."
                },
                {
                    id: "z37",
                    text: "Which part of stomach connects to small intestine?",
                    options: ["Fundus", "Body", "Antrum", "Pylorus"],
                    answer: 3,
                    explain: "Pyloric region connects stomach to duodenum through pyloric sphincter."
                },
                {
                    id: "z38",
                    text: "The condition of excessive glucose in urine is called:",
                    options: ["Glycemia", "Glycosuria", "Gluconeogenesis", "Glycogenolysis"],
                    answer: 1,
                    explain: "Glycosuria is presence of glucose in urine, typically occurs when blood glucose exceeds renal threshold."
                },
                {
                    id: "z39",
                    text: "Which of the following produces growth hormone?",
                    options: ["Posterior pituitary", "Anterior pituitary", "Hypothalamus", "Thyroid"],
                    answer: 1,
                    explain: "Growth hormone is secreted by anterior pituitary gland (adenohypophysis)."
                },
                {
                    id: "z40",
                    text: "The primary function of large intestine is:",
                    options: ["Digestion", "Absorption of nutrients", "Water absorption", "Enzyme secretion"],
                    answer: 2,
                    explain: "Large intestine primarily absorbs water and electrolytes from indigestible food matter."
                },
                {
                    id: "z41",
                    text: "Which cranial nerve controls facial expressions?",
                    options: ["Trigeminal", "Facial", "Glossopharyngeal", "Vagus"],
                    answer: 1,
                    explain: "Facial nerve (VII cranial nerve) controls muscles of facial expression."
                },
                {
                    id: "z42",
                    text: "The lifespan of human RBCs is approximately:",
                    options: ["60 days", "90 days", "120 days", "150 days"],
                    answer: 2,
                    explain: "Human red blood cells have average lifespan of 120 days before being removed by spleen."
                },
                {
                    id: "z43",
                    text: "Which of the following is involved in temperature regulation?",
                    options: ["Cerebrum", "Cerebellum", "Hypothalamus", "Medulla"],
                    answer: 2,
                    explain: "Hypothalamus acts as thermostat, regulating body temperature through various mechanisms."
                },
                {
                    id: "z44",
                    text: "The process of formation of female gametes is called:",
                    options: ["Spermatogenesis", "Oogenesis", "Gametogenesis", "Fertilization"],
                    answer: 1,
                    explain: "Oogenesis is process of ovum formation in ovaries from primordial germ cells."
                },
                {
                    id: "z45",
                    text: "Which of the following prevents entry of food into windpipe?",
                    options: ["Epiglottis", "Uvula", "Soft palate", "Tonsils"],
                    answer: 0,
                    explain: "Epiglottis covers laryngeal opening during swallowing to prevent food aspiration."
                }
            ]
        }
    ]
};
