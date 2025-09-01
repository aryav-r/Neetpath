// mock-test-2.js - NEET Mock Test 2 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// Follows exact NEET pattern and difficulty

window.MOCK_TEST_2 = {
    id: "neet-002",
    title: "Full Syllabus Mock 2", 
    level: "medium",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A body is thrown vertically upward with velocity u. The time to reach maximum height and the time to return to the ground are in the ratio:",
                    options: ["1:1", "1:2", "2:1", "1:√2"],
                    answer: 0,
                    explain: "Time to reach max height = u/g, time to return = 2u/g. Ratio = (u/g):(2u/g) = 1:2. But time up equals time down, so 1:1."
                },
                {
                    id: "p2",
                    text: "A uniform rod of length L and mass M is pivoted at L/4 from one end. The moment of inertia about the pivot is:",
                    options: ["ML²/12", "7ML²/48", "3ML²/16", "5ML²/48"],
                    answer: 1,
                    explain: "I = ICM + Md² where ICM = ML²/12, d = L/4. I = ML²/12 + M(L/4)² = ML²/12 + ML²/16 = 7ML²/48"
                },
                {
                    id: "p3",
                    text: "In a series RLC circuit, the phase difference between voltage across R and voltage across LC combination is:",
                    options: ["0°", "90°", "180°", "270°"],
                    answer: 1,
                    explain: "Voltage across R is in phase with current, voltage across LC is 90° out of phase with current, so phase difference is 90°."
                },
                {
                    id: "p4",
                    text: "A satellite orbiting Earth has kinetic energy K. Its total energy is:",
                    options: ["K", "-K", "K/2", "-K/2"],
                    answer: 3,
                    explain: "For circular orbits: KE = K, PE = -2K, Total energy = K + (-2K) = -K. So E = -K/2 when K is given."
                },
                {
                    id: "p5",
                    text: "The ratio of de Broglie wavelengths of electron and proton accelerated through same potential is:",
                    options: ["1:1", "√(mp/me)", "√(me/mp)", "mp/me"],
                    answer: 2,
                    explain: "λ = h/√(2mV). λe/λp = √(mp/me) since both accelerated through same V."
                },
                {
                    id: "p6",
                    text: "A convex mirror with focal length 20 cm forms image at 12 cm. The object distance is:",
                    options: ["-30 cm", "30 cm", "-7.5 cm", "7.5 cm"],
                    answer: 0,
                    explain: "1/f = 1/u + 1/v ⟹ 1/(-20) = 1/u + 1/(-12) ⟹ u = -30 cm"
                },
                {
                    id: "p7",
                    text: "In photoelectric effect, if frequency of incident light is doubled, the stopping potential becomes:",
                    options: ["Double", "Four times", "Increases by hf/e", "Unchanged"],
                    answer: 2,
                    explain: "eV₀ = hf - φ. When frequency doubles: eV₀' = 2hf - φ = (hf - φ) + hf = eV₀ + hf. So V₀' increases by hf/e."
                },
                {
                    id: "p8",
                    text: "A spring-block system has period T. If both mass and spring constant are doubled, new period is:",
                    options: ["T", "T/2", "2T", "T√2"],
                    answer: 0,
                    explain: "T = 2π√(m/k). When both m and k are doubled: T' = 2π√(2m/2k) = 2π√(m/k) = T"
                },
                {
                    id: "p9",
                    text: "Two coherent sources with phase difference π interfere. If individual amplitudes are A, resultant amplitude is:",
                    options: ["0", "A", "2A", "A√2"],
                    answer: 0,
                    explain: "Resultant amplitude = √(A² + A² + 2A²cosπ) = √(2A² - 2A²) = 0 (destructive interference)"
                },
                {
                    id: "p10",
                    text: "A charged particle moves in a uniform magnetic field in a circular path. Its kinetic energy:",
                    options: ["Increases", "Decreases", "Remains constant", "First increases then decreases"],
                    answer: 2,
                    explain: "Magnetic force is always perpendicular to velocity, does no work, so kinetic energy remains constant."
                },
                {
                    id: "p11",
                    text: "The magnetic field inside a long solenoid with n turns per unit length carrying current I is:",
                    options: ["μ₀I", "μ₀nI", "μ₀nI/2", "2μ₀nI"],
                    answer: 1,
                    explain: "Magnetic field inside long solenoid: B = μ₀nI where n is turns per unit length."
                },
                {
                    id: "p12",
                    text: "In Young's double slit experiment, fringe width is 0.5 mm. If screen distance is doubled, fringe width becomes:",
                    options: ["0.25 mm", "0.5 mm", "1.0 mm", "2.0 mm"],
                    answer: 2,
                    explain: "Fringe width β = λD/d. When D is doubled, β becomes 2 × 0.5 = 1.0 mm"
                },
                {
                    id: "p13",
                    text: "A gas expands adiabatically to twice its volume. If γ = 1.4, pressure becomes:",
                    options: ["P/2", "P/2^1.4", "P/2^0.4", "P√2"],
                    answer: 1,
                    explain: "For adiabatic process: PVᵞ = constant. P₁V₁ᵞ = P₂V₂ᵞ ⟹ P₂ = P₁(V₁/V₂)ᵞ = P₁(1/2)^1.4 = P₁/2^1.4"
                },
                {
                    id: "p14",
                    text: "The half-life of a radioactive sample is 10 years. What fraction decays in 30 years?",
                    options: ["1/8", "7/8", "1/4", "3/4"],
                    answer: 1,
                    explain: "In 30 years = 3 half-lives, remaining fraction = (1/2)³ = 1/8. Decayed fraction = 1 - 1/8 = 7/8"
                },
                {
                    id: "p15",
                    text: "A wire carrying current I is bent into a circular loop. The magnetic field at center is B. If same wire is bent into 2 turns, field at center becomes:",
                    options: ["B/2", "B", "2B", "4B"],
                    answer: 3,
                    explain: "For 2 turns with half radius: B ∝ nI/r where n = 2, r = r₀/2. So B' = 4B₀"
                },
                {
                    id: "p16",
                    text: "A body moves with constant speed in a circular path. Which statement is correct?",
                    options: [
                        "Acceleration is zero",
                        "Velocity is constant", 
                        "Acceleration is centripetal",
                        "No force acts on body"
                    ],
                    answer: 2,
                    explain: "Speed constant but velocity changes direction, requiring centripetal acceleration toward center."
                },
                {
                    id: "p17",
                    text: "The binding energy per nucleon is maximum for:",
                    options: ["Light nuclei", "Heavy nuclei", "Medium mass nuclei", "All nuclei"],
                    answer: 2,
                    explain: "Binding energy per nucleon peaks around Fe-56 (medium mass), making these most stable."
                },
                {
                    id: "p18",
                    text: "In a step-up transformer, if primary has 100 turns and secondary has 500 turns, voltage ratio is:",
                    options: ["1:5", "5:1", "1:25", "25:1"],
                    answer: 0,
                    explain: "Voltage ratio = turns ratio = Ns/Np = 500/100 = 5:1. Secondary voltage is 5 times primary."
                },
                {
                    id: "p19",
                    text: "The coefficient of restitution for perfectly elastic collision is:",
                    options: ["0", "1", "∞", "-1"],
                    answer: 1,
                    explain: "Coefficient of restitution e = 1 for perfectly elastic collision (kinetic energy conserved)."
                },
                {
                    id: "p20",
                    text: "A lens has focal length 20 cm in air. When immersed in water (μ = 4/3), focal length becomes:",
                    options: ["15 cm", "26.7 cm", "80 cm", "60 cm"],
                    answer: 2,
                    explain: "1/f' = (μ_lens/μ_medium - 1)(1/R₁ - 1/R₂). f' = f × μ_medium/(μ_lens - μ_medium) ≈ 80 cm"
                },
                {
                    id: "p21",
                    text: "A particle executes SHM with amplitude 4 cm. At what displacement is kinetic energy 3 times potential energy?",
                    options: ["1 cm", "2 cm", "3 cm", "√3 cm"],
                    answer: 1,
                    explain: "KE = 3PE and KE + PE = Total energy. So 3PE + PE = 4PE = Total = ½kA². At x: PE = ½kx², so ½kx² = ¼(½kA²) ⟹ x = A/2 = 2 cm"
                },
                {
                    id: "p22",
                    text: "Two resistors of 4Ω and 12Ω are connected in parallel. Their equivalent resistance is:",
                    options: ["3Ω", "4Ω", "8Ω", "16Ω"],
                    answer: 0,
                    explain: "1/R = 1/4 + 1/12 = 3/12 + 1/12 = 4/12 = 1/3. So R = 3Ω"
                },
                {
                    id: "p23",
                    text: "The work function of a metal is 2.5 eV. The threshold frequency is:",
                    options: ["6.05 × 10¹⁴ Hz", "6.05 × 10¹³ Hz", "3.8 × 10¹⁵ Hz", "2.5 × 10¹⁵ Hz"],
                    answer: 0,
                    explain: "hf₀ = φ ⟹ f₀ = φ/h = 2.5 × 1.6 × 10⁻¹⁹/(6.63 × 10⁻³⁴) = 6.05 × 10¹⁴ Hz"
                },
                {
                    id: "p24",
                    text: "A ball is thrown horizontally from height h. Time of flight depends on:",
                    options: ["Initial velocity only", "Height only", "Both height and velocity", "Mass of ball"],
                    answer: 1,
                    explain: "For projectile motion: t = √(2h/g). Time depends only on height, not horizontal velocity."
                },
                {
                    id: "p25",
                    text: "The electric potential at distance r from point charge q is V. The electric field at same point is:",
                    options: ["V/r", "V/r²", "Vr", "V²/r"],
                    answer: 0,
                    explain: "V = kq/r and E = kq/r². Therefore E = V/r"
                },
                {
                    id: "p26",
                    text: "A conducting rod moves in magnetic field. The induced EMF depends on:",
                    options: ["Velocity only", "Field strength only", "Length of rod only", "All of the above"],
                    answer: 3,
                    explain: "Motional EMF = BLv where B is field, L is length, v is velocity. Depends on all three."
                },
                {
                    id: "p27",
                    text: "The ratio of wavelengths of first Balmer and first Lyman lines in hydrogen spectrum is:",
                    options: ["27:5", "5:27", "9:4", "4:9"],
                    answer: 0,
                    explain: "First Balmer (n=3→2): 1/λ₁ ∝ (1/4 - 1/9) = 5/36. First Lyman (n=2→1): 1/λ₂ ∝ (1/1 - 1/4) = 3/4. λ₁/λ₂ = (3/4)/(5/36) = 27/5"
                },
                {
                    id: "p28",
                    text: "A carnot engine operates between 400K and 300K. Its efficiency is:",
                    options: ["25%", "33%", "75%", "67%"],
                    answer: 0,
                    explain: "Efficiency = 1 - T₂/T₁ = 1 - 300/400 = 1 - 0.75 = 0.25 = 25%"
                },
                {
                    id: "p29",
                    text: "The momentum of photon of wavelength λ is:",
                    options: ["hλ", "h/λ", "hλ/c", "hc/λ"],
                    answer: 1,
                    explain: "Photon momentum p = E/c = hf/c = hc/λc = h/λ"
                },
                {
                    id: "p30",
                    text: "In a uniform electric field, equipotential surfaces are:",
                    options: ["Spherical", "Cylindrical", "Plane surfaces", "Irregular"],
                    answer: 2,
                    explain: "In uniform field, equipotential surfaces are planes perpendicular to field direction."
                },
                {
                    id: "p31",
                    text: "Two identical capacitors are first connected in series, then in parallel. The ratio of energies stored is:",
                    options: ["1:2", "1:4", "2:1", "4:1"],
                    answer: 1,
                    explain: "Series: C_eq = C/2, E₁ = ½(C/2)V² = CV²/4. Parallel: C_eq = 2C, E₂ = ½(2C)V² = CV². Ratio = 1:4"
                },
                {
                    id: "p32",
                    text: "A mass attached to spring undergoes SHM. If mass is doubled, frequency becomes:",
                    options: ["f/2", "f/√2", "f√2", "2f"],
                    answer: 1,
                    explain: "f = (1/2π)√(k/m). When mass doubles: f' = (1/2π)√(k/2m) = f/√2"
                },
                {
                    id: "p33",
                    text: "The magnetic flux through a loop is maximum when angle between field and normal is:",
                    options: ["0°", "30°", "60°", "90°"],
                    answer: 0,
                    explain: "Flux Φ = BA cos θ is maximum when θ = 0° (field parallel to normal)."
                },
                {
                    id: "p34",
                    text: "A wire of resistance R is stretched to double its length. The new resistance becomes:",
                    options: ["R/4", "R/2", "2R", "4R"],
                    answer: 3,
                    explain: "R ∝ L/A. When length doubles, area halves: R' = R × 2 × 2 = 4R"
                },
                {
                    id: "p35",
                    text: "The dimensional formula of Planck's constant is:",
                    options: ["ML²T⁻¹", "ML²T⁻²", "MLT⁻¹", "ML⁻¹T⁻¹"],
                    answer: 0,
                    explain: "h = E/f has dimensions [ML²T⁻²]/[T⁻¹] = [ML²T⁻¹]"
                },
                {
                    id: "p36",
                    text: "In an AC circuit with only capacitor, current leads voltage by:",
                    options: ["0°", "45°", "90°", "180°"],
                    answer: 2,
                    explain: "In purely capacitive circuit, current leads voltage by 90°."
                },
                {
                    id: "p37",
                    text: "The gravitational potential energy of system of two masses separated by distance r is:",
                    options: ["-Gm₁m₂/r", "Gm₁m₂/r", "-Gm₁m₂/r²", "Gm₁m₂/r²"],
                    answer: 0,
                    explain: "Gravitational potential energy U = -Gm₁m₂/r (negative because it's attractive)."
                },
                {
                    id: "p38",
                    text: "The speed of electromagnetic waves in vacuum is:",
                    options: ["√(μ₀/ε₀)", "√(ε₀/μ₀)", "1/√(μ₀ε₀)", "√(μ₀ε₀)"],
                    answer: 2,
                    explain: "Speed of EM waves c = 1/√(μ₀ε₀) where μ₀ is permeability and ε₀ is permittivity of free space."
                },
                {
                    id: "p39",
                    text: "A particle moving in circular path has centripetal acceleration 4 m/s². If radius is doubled, acceleration becomes:",
                    options: ["1 m/s²", "2 m/s²", "8 m/s²", "16 m/s²"],
                    answer: 1,
                    explain: "a = v²/r. If radius doubles and speed constant: a' = v²/(2r) = a/2 = 2 m/s²"
                },
                {
                    id: "p40",
                    text: "The electric field inside a hollow charged conductor is:",
                    options: ["Maximum", "Minimum", "Zero", "Variable"],
                    answer: 2,
                    explain: "Electric field inside hollow conductor is zero due to electrostatic shielding."
                },
                {
                    id: "p41",
                    text: "In Compton scattering, the wavelength shift depends on:",
                    options: ["Initial wavelength", "Final wavelength", "Scattering angle", "Energy of photon"],
                    answer: 2,
                    explain: "Compton shift Δλ = (h/mc)(1 - cos θ) depends only on scattering angle θ."
                },
                {
                    id: "p42",
                    text: "The resolving power of microscope can be increased by:",
                    options: ["Using light of longer wavelength", "Using light of shorter wavelength", "Decreasing aperture", "Increasing distance"],
                    answer: 1,
                    explain: "Resolving power ∝ 1/λ. Shorter wavelength gives better resolution."
                },
                {
                    id: "p43",
                    text: "A liquid drop has surface energy E. When it breaks into 8 droplets, surface energy becomes:",
                    options: ["E/2", "2E", "4E", "8E"],
                    answer: 1,
                    explain: "Surface area increases by factor 2 when drop breaks into 8 equal droplets, so surface energy becomes 2E."
                },
                {
                    id: "p44",
                    text: "The magnetic moment of current loop depends on:",
                    options: ["Current only", "Area only", "Number of turns only", "All of the above"],
                    answer: 3,
                    explain: "Magnetic moment m = nIA where n = turns, I = current, A = area. Depends on all three."
                },
                {
                    id: "p45",
                    text: "In X-ray spectrum, Kα line arises due to transition:",
                    options: ["L→K", "M→L", "N→M", "K→L"],
                    answer: 0,
                    explain: "Kα line results from electron transition from L shell to K shell (L→K)."
                }
            ]
        },
        {
            name: "Chemistry", 
            questions: [
                {
                    id: "c1",
                    text: "Which of the following species has maximum number of unpaired electrons?",
                    options: ["Mn²⁺", "Fe³⁺", "Co²⁺", "Ni²⁺"],
                    answer: 0,
                    explain: "Mn²⁺: [Ar] 3d⁵ has 5 unpaired electrons (high spin). Fe³⁺: [Ar] 3d⁵ also has 5, but Mn²⁺ is d⁵ system."
                },
                {
                    id: "c2",
                    text: "The correct increasing order of ionic radii is:",
                    options: [
                        "Li⁺ < Na⁺ < K⁺ < Rb⁺",
                        "Rb⁺ < K⁺ < Na⁺ < Li⁺", 
                        "Na⁺ < Li⁺ < K⁺ < Rb⁺",
                        "K⁺ < Na⁺ < Rb⁺ < Li⁺"
                    ],
                    answer: 0,
                    explain: "Ionic radii increase down the group due to addition of electron shells: Li⁺ < Na⁺ < K⁺ < Rb⁺"
                },
                {
                    id: "c3",
                    text: "The molecular geometry of IF₅ is:",
                    options: ["Trigonal bipyramidal", "Square pyramidal", "Pentagonal planar", "Octahedral"],
                    answer: 1,
                    explain: "IF₅ has 6 electron pairs (5 bonding + 1 lone pair) around I, giving square pyramidal geometry."
                },
                {
                    id: "c4",
                    text: "Which of the following shows maximum catenation?",
                    options: ["C", "Si", "Ge", "Sn"],
                    answer: 0,
                    explain: "Carbon shows maximum catenation due to small size, high bond strength, and ability to form multiple bonds."
                },
                {
                    id: "c5",
                    text: "The number of stereoisomers possible for CH₃CHBrCHBrCH₃ is:",
                    options: ["2", "3", "4", "6"],
                    answer: 1,
                    explain: "Two chiral centers give 2ⁿ = 4 stereoisomers, but meso form reduces it to 3 (2 enantiomers + 1 meso)."
                },
                {
                    id: "c6",
                    text: "Which of the following has highest boiling point?",
                    options: ["NH₃", "PH₃", "AsH₃", "SbH₃"],
                    answer: 0,
                    explain: "NH₃ has highest boiling point due to strong hydrogen bonding, despite being smallest molecule."
                },
                {
                    id: "c7",
                    text: "The oxidation state of phosphorus in H₃PO₂ is:",
                    options: ["+1", "+2", "+3", "+5"],
                    answer: 0,
                    explain: "In H₃PO₂: 3(+1) + P + 2(-2) = 0 ⟹ P = +1. Also, it's a monobasic acid."
                },
                {
                    id: "c8",
                    text: "Which reagent is used for Markovnikov addition of HBr to alkenes?",
                    options: ["HBr + peroxides", "HBr in dark", "HBr + UV light", "HBr + heat"],
                    answer: 1,
                    explain: "HBr in dark (absence of peroxides) follows Markovnikov rule. Peroxides cause anti-Markovnikov addition."
                },
                {
                    id: "c9",
                    text: "The hybridization of central atom in PCl₅ is:",
                    options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
                    answer: 1,
                    explain: "PCl₅ has 5 bonding pairs around P, requiring sp³d hybridization (trigonal bipyramidal)."
                },
                {
                    id: "c10",
                    text: "Which of the following is most stable carbocation?",
                    options: ["CH₃⁺", "(CH₃)₂CH⁺", "(CH₃)₃C⁺", "C₆H₅CH₂⁺"],
                    answer: 3,
                    explain: "Benzyl carbocation is most stable due to resonance delocalization with benzene ring."
                },
                {
                    id: "c11",
                    text: "The correct order of acidic strength is:",
                    options: [
                        "HF > HCl > HBr > HI",
                        "HI > HBr > HCl > HF", 
                        "HCl > HF > HBr > HI",
                        "HBr > HI > HCl > HF"
                    ],
                    answer: 1,
                    explain: "Acidic strength increases down the group due to decreasing bond strength: HI > HBr > HCl > HF"
                },
                {
                    id: "c12",
                    text: "Which of the following exhibits optical isomerism?",
                    options: ["CHCl₃", "CH₂Cl₂", "CHClBrI", "CCl₄"],
                    answer: 2,
                    explain: "CHClBrI has four different groups attached to carbon, making it chiral and optically active."
                },
                {
                    id: "c13",
                    text: "The rate law for reaction A + B → C is rate = k[A][B]. The order of reaction is:",
                    options: ["Zero", "First", "Second", "Third"],
                    answer: 2,
                    explain: "Overall order = sum of individual orders = 1 + 1 = 2 (second order)"
                },
                {
                    id: "c14",
                    text: "Which of the following is not a buffer solution?",
                    options: [
                        "CH₃COOH + CH₃COONa",
                        "NH₄OH + NH₄Cl", 
                        "HCl + NaCl",
                        "H₂CO₃ + NaHCO₃"
                    ],
                    answer: 2,
                    explain: "HCl + NaCl doesn't form buffer as HCl is strong acid. Buffers need weak acid/base with their salts."
                },
                {
                    id: "c15",
                    text: "The IUPAC name of (CH₃)₂CHCH₂COOH is:",
                    options: [
                        "3-methylbutanoic acid",
                        "2-methylbutanoic acid", 
                        "Isovaleric acid",
                        "3-methylpentanoic acid"
                    ],
                    answer: 0,
                    explain: "5-carbon chain with COOH at C1, methyl branch at C3: 3-methylbutanoic acid"
                },
                {
                    id: "c16",
                    text: "Which of the following has sp² hybridization?",
                    options: ["BeCl₂", "BCl₃", "CCl₄", "PCl₅"],
                    answer: 1,
                    explain: "BCl₃ has 3 bonding pairs around B with trigonal planar geometry requiring sp² hybridization."
                },
                {
                    id: "c17",
                    text: "The complex [Ni(CO)₄] has:",
                    options: ["Square planar geometry", "Tetrahedral geometry", "Octahedral geometry", "Linear geometry"],
                    answer: 1,
                    explain: "[Ni(CO)₄] has tetrahedral geometry. Ni is in 0 oxidation state with sp³ hybridization."
                },
                {
                    id: "c18",
                    text: "Which of the following is strongest reducing agent?",
                    options: ["H₂", "CO", "C", "Al"],
                    answer: 3,
                    explain: "Aluminum is strongest reducing agent with highest tendency to lose electrons (lowest reduction potential)."
                },
                {
                    id: "c19",
                    text: "The number of α-hydrogen atoms in propanoic acid is:",
                    options: ["2", "3", "5", "8"],
                    answer: 0,
                    explain: "CH₃CH₂COOH: α-hydrogens are on carbon adjacent to COOH group. CH₂ has 2 α-hydrogens."
                },
                {
                    id: "c20",
                    text: "Which of the following undergoes nucleophilic addition reaction?",
                    options: ["Benzene", "Ethene", "Ethyne", "Ethanal"],
                    answer: 3,
                    explain: "Ethanal (acetaldehyde) has C=O group which undergoes nucleophilic addition due to electrophilic carbon."
                },
                {
                    id: "c21",
                    text: "The electronegativity difference in HF molecule is:",
                    options: ["1.9", "2.1", "1.7", "1.5"],
                    answer: 0,
                    explain: "Electronegativity: F = 4.0, H = 2.1. Difference = 4.0 - 2.1 = 1.9"
                },
                {
                    id: "c22",
                    text: "Which of the following has maximum dipole moment?",
                    options: ["CO₂", "BF₃", "NH₃", "CCl₄"],
                    answer: 2,
                    explain: "NH₃ has pyramidal shape with lone pair, creating maximum dipole moment among given options."
                },
                {
                    id: "c23",
                    text: "The coordination number of Ni²⁺ in [Ni(en)₃]²⁺ is:",
                    options: ["3", "4", "6", "8"],
                    answer: 2,
                    explain: "Ethylenediamine (en) is bidentate ligand. 3 en ligands × 2 coordination sites = 6"
                },
                {
                    id: "c24",
                    text: "Which of the following shows maximum hydrogen bonding?",
                    options: ["H₂O", "HF", "NH₃", "H₂S"],
                    answer: 1,
                    explain: "HF shows strongest hydrogen bonding due to highest electronegativity of fluorine."
                },
                {
                    id: "c25",
                    text: "The number of π bonds in benzene is:",
                    options: ["3", "6", "9", "12"],
                    answer: 0,
                    explain: "Benzene has 3 delocalized π bonds (6 π electrons forming 3 π molecular orbitals)."
                },
                {
                    id: "c26",
                    text: "Which of the following is paramagnetic?",
                    options: ["N₂", "O₂", "F₂", "Ne₂"],
                    answer: 1,
                    explain: "O₂ has 2 unpaired electrons in antibonding π* orbitals, making it paramagnetic."
                },
                {
                    id: "c27",
                    text: "The bond order in CO molecule is:",
                    options: ["1", "2", "3", "2.5"],
                    answer: 2,
                    explain: "CO has triple bond character: C≡O with bond order 3 (2σ + 1π from molecular orbital theory)."
                },
                {
                    id: "c28",
                    text: "Which of the following is most basic amine?",
                    options: ["Aniline", "Methylamine", "Dimethylamine", "Ammonia"],
                    answer: 2,
                    explain: "Dimethylamine is most basic due to +I effect of two methyl groups and optimal steric factors."
                },
                {
                    id: "c29",
                    text: "The shape of ClO₄⁻ ion is:",
                    options: ["Tetrahedral", "Square planar", "Pyramidal", "Bent"],
                    answer: 0,
                    explain: "ClO₄⁻ has 4 bonding pairs around Cl with no lone pairs, giving tetrahedral geometry."
                },
                {
                    id: "c30",
                    text: "Which of the following has highest lattice energy?",
                    options: ["LiF", "NaCl", "KBr", "CsI"],
                    answer: 0,
                    explain: "LiF has highest lattice energy due to small size and high charges of Li⁺ and F⁻ ions."
                },
                {
                    id: "c31",
                    text: "The number of moles of KMnO₄ required to oxidize 1 mole of FeSO₄ in acidic medium is:",
                    options: ["1", "0.2", "5", "0.6"],
                    answer: 1,
                    explain: "MnO₄⁻ + 5Fe²⁺ + 8H⁺ → Mn²⁺ + 5Fe³⁺ + 4H₂O. 1 mol KMnO₄ oxidizes 5 mol FeSO₄, so need 1/5 = 0.2 mol"
                },
                {
                    id: "c32",
                    text: "Which of the following exhibits tautomerism?",
                    options: ["Phenol", "Acetone", "Ethanol", "Diethyl ether"],
                    answer: 1,
                    explain: "Acetone shows keto-enol tautomerism: CH₃COCH₃ ⇌ CH₃C(OH)=CH₂"
                },
                {
                    id: "c33",
                    text: "The first ionization energy is lowest for:",
                    options: ["Li", "Na", "K", "Cs"],
                    answer: 3,
                    explain: "Ionization energy decreases down the group. Cs has lowest IE due to largest atomic size."
                },
                {
                    id: "c34",
                    text: "Which of the following is not aromatic?",
                    options: ["Benzene", "Pyridine", "Furan", "Cyclooctatetraene"],
                    answer: 3,
                    explain: "Cyclooctatetraene is not planar and doesn't follow Hückel's rule, hence not aromatic."
                },
                {
                    id: "c35",
                    text: "The unit of rate constant for zero order reaction is:",
                    options: ["mol L⁻¹ s⁻¹", "s⁻¹", "mol⁻¹ L s⁻¹", "mol⁻² L² s⁻¹"],
                    answer: 0,
                    explain: "For zero order: rate = k[A]⁰ = k. So k has same units as rate: mol L⁻¹ s⁻¹"
                },
                {
                    id: "c36",
                    text: "Which of the following is strongest acid?",
                    options: ["HClO", "HClO₂", "HClO₃", "HClO₄"],
                    answer: 3,
                    explain: "HClO₄ is strongest due to maximum +I effect of oxygen atoms stabilizing conjugate base."
                },
                {
                    id: "c37",
                    text: "The hybridization state of carbon in diamond is:",
                    options: ["sp", "sp²", "sp³", "sp³d"],
                    answer: 2,
                    explain: "Each carbon in diamond is bonded to 4 other carbons tetrahedrally, requiring sp³ hybridization."
                },
                {
                    id: "c38",
                    text: "Which of the following shows linkage isomerism?",
                    options: ["[Co(NH₃)₆]Cl₃", "[Co(NH₃)₅NO₂]Cl₂", "[Co(NH₃)₄Cl₂]Cl", "[Co(en)₃]Cl₃"],
                    answer: 1,
                    explain: "NO₂⁻ can bind through N or O, showing linkage isomerism in [Co(NH₃)₅NO₂]Cl₂"
                },
                {
                    id: "c39",
                    text: "The oxidation number of Cr in CrO₅ is:",
                    options: ["+6", "+8", "+10", "+5"],
                    answer: 0,
                    explain: "In CrO₅ (chromium peroxide): Cr + 5(-2) = -4. But structure is CrO(O₂)₂, so Cr = +6"
                },
                {
                    id: "c40",
                    text: "Which of the following is Lewis base but not Brønsted base?",
                    options: ["NH₃", "OH⁻", "BF₃", "CO"],
                    answer: 3,
                    explain: "CO can donate electron pair (Lewis base) but cannot accept proton (not Brønsted base)."
                },
                {
                    id: "c41",
                    text: "The correct order of stability of carbocations is:",
                    options: [
                        "1° > 2° > 3°",
                        "3° > 2° > 1°", 
                        "2° > 1° > 3°",
                        "1° = 2° = 3°"
                    ],
                    answer: 1,
                    explain: "Carbocation stability: 3° > 2° > 1° due to hyperconjugation and inductive effects."
                },
                {
                    id: "c42",
                    text: "Which of the following has maximum melting point?",
                    options: ["NaCl", "MgO", "AlN", "SiC"],
                    answer: 3,
                    explain: "SiC has maximum melting point due to strong covalent bonding in giant covalent structure."
                },
                {
                    id: "c43",
                    text: "The number of unpaired electrons in Fe³⁺ (high spin) is:",
                    options: ["3", "4", "5", "6"],
                    answer: 2,
                    explain: "Fe³⁺: [Ar] 3d⁵. In high spin state, all 5 d-electrons are unpaired in separate orbitals."
                },
                {
                    id: "c44",
                    text: "Which of the following undergoes aldol condensation?",
                    options: ["Formaldehyde", "Acetaldehyde", "Benzaldehyde", "Acetone"],
                    answer: 1,
                    explain: "Acetaldehyde has α-hydrogen atoms required for aldol condensation reaction."
                },
                {
                    id: "c45",
                    text: "The entropy change is maximum for:",
                    options: ["Melting", "Vaporization", "Sublimation", "Crystallization"],
                    answer: 1,
                    explain: "Vaporization shows maximum entropy increase due to largest change in molecular disorder."
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY Questions (45)
                {
                    id: "b1",
                    text: "Which of the following is characteristic of C4 photosynthesis?",
                    options: [
                        "Spatial separation of initial CO₂ fixation and Calvin cycle",
                        "Temporal separation of CO₂ fixation", 
                        "No photorespiration occurs",
                        "Both A and C"
                    ],
                    answer: 3,
                    explain: "C4 plants have spatial separation (mesophyll and bundle sheath cells) and reduced photorespiration due to CO₂ concentrating mechanism."
                },
                {
                    id: "b2",
                    text: "The enzyme responsible for CO₂ fixation in CAM plants is:",
                    options: ["RuBisCO", "PEP carboxylase", "Carbonic anhydrase", "Both A and B"],
                    answer: 3,
                    explain: "CAM plants use PEP carboxylase at night and RuBisCO during day for temporal separation of CO₂ fixation."
                },
                {
                    id: "b3",
                    text: "Which plant hormone is synthesized in root tips and transported upwards?",
                    options: ["Auxin", "Cytokinin", "Gibberellin", "Abscisic acid"],
                    answer: 1,
                    explain: "Cytokinins are synthesized in root apical meristems and transported through xylem to shoots."
                },
                {
                    id: "b4",
                    text: "The phenomenon of bending of plant parts towards light is called:",
                    options: ["Geotropism", "Phototropism", "Chemotropism", "Thigmotropism"],
                    answer: 1,
                    explain: "Phototropism is growth movement in response to unidirectional light stimulus, mediated by auxin."
                },
                {
                    id: "b5",
                    text: "Which of the following is NOT involved in opening of stomata?",
                    options: ["K⁺ influx", "Malate synthesis", "Starch to sugar conversion", "Ca²⁺ influx"],
                    answer: 3,
                    explain: "Ca²⁺ influx is involved in stomatal closure, not opening. K⁺ influx, malate synthesis, and starch breakdown promote opening."
                },
                {
                    id: "b6",
                    text: "The site of protein synthesis in chloroplasts is:",
                    options: ["Stroma", "Thylakoid membrane", "Inner membrane", "Intermembrane space"],
                    answer: 0,
                    explain: "70S ribosomes in chloroplast stroma synthesize proteins encoded by chloroplast DNA."
                },
                {
                    id: "b7",
                    text: "Which of the following represents gametophytic generation in angiosperms?",
                    options: ["Entire plant body", "Pollen grain and embryo sac", "Only embryo sac", "Seed and fruit"],
                    answer: 1,
                    explain: "Both pollen grain (male gametophyte) and embryo sac (female gametophyte) represent gametophytic generation."
                },
                {
                    id: "b8",
                    text: "The term 'polyembryony' refers to:",
                    options: [
                        "Multiple ovules in ovary",
                        "Multiple embryos in single seed", 
                        "Multiple seeds in fruit",
                        "Multiple flowers in inflorescence"
                    ],
                    answer: 1,
                    explain: "Polyembryony is occurrence of more than one embryo in a single seed, common in citrus."
                },
                {
                    id: "b9",
                    text: "Which of the following is essential for synthesis of chlorophyll?",
                    options: ["Iron", "Magnesium", "Calcium", "Potassium"],
                    answer: 1,
                    explain: "Magnesium is central atom in chlorophyll molecule, essential for its synthesis and function."
                },
                {
                    id: "b10",
                    text: "The pressure flow hypothesis explains:",
                    options: ["Water transport in xylem", "Translocation in phloem", "Gas exchange", "Ion uptake"],
                    answer: 1,
                    explain: "Pressure flow hypothesis explains translocation of sugars in phloem from source to sink."
                },
                {
                    id: "b11",
                    text: "Which type of fruit is formed without fertilization?",
                    options: ["True fruit", "False fruit", "Parthenocarpic fruit", "Aggregate fruit"],
                    answer: 2,
                    explain: "Parthenocarpic fruits develop without fertilization, seedless fruits like bananas and grapes."
                },
                {
                    id: "b12",
                    text: "The vascular bundles in monocot stem are:",
                    options: ["Open and scattered", "Closed and scattered", "Open and arranged in ring", "Closed and arranged in ring"],
                    answer: 1,
                    explain: "Monocot stems have closed vascular bundles (no cambium) scattered throughout the ground tissue."
                },
                {
                    id: "b13",
                    text: "Which of the following is a day-neutral plant?",
                    options: ["Chrysanthemum", "Spinach", "Cucumber", "Poinsettia"],
                    answer: 2,
                    explain: "Cucumber is day-neutral plant, flowering is not affected by photoperiod length."
                },
                {
                    id: "b14",
                    text: "The process by which water is lost from plant surfaces is called:",
                    options: ["Transpiration", "Guttation", "Bleeding", "Imbibition"],
                    answer: 0,
                    explain: "Transpiration is loss of water vapor from plant surfaces, mainly through stomata."
                },
                {
                    id: "b15",
                    text: "Which of the following shows circadian rhythm?",
                    options: ["Stomatal opening", "Leaf movement", "Flower opening", "All of above"],
                    answer: 3,
                    explain: "All these processes show circadian rhythms - endogenous biological clocks with ~24 hour cycles."
                },
                {
                    id: "b16",
                    text: "The phenomenon where pollen from same flower pollinates the stigma is:",
                    options: ["Autogamy", "Geitonogamy", "Xenogamy", "Allogamy"],
                    answer: 0,
                    explain: "Autogamy is self-pollination within the same flower, ensuring reproductive success."
                },
                {
                    id: "b17",
                    text: "Which plant tissue is responsible for radial conduction of water?",
                    options: ["Xylem vessels", "Xylem parenchyma", "Phloem parenchyma", "Cambium"],
                    answer: 1,
                    explain: "Xylem parenchyma cells facilitate radial movement of water and minerals from xylem to surrounding tissues."
                },
                {
                    id: "b18",
                    text: "The term 'vernalization' refers to:",
                    options: [
                        "Cold treatment for flowering",
                        "Light treatment for flowering", 
                        "Hormone treatment",
                        "Nutrient treatment"
                    ],
                    answer: 0,
                    explain: "Vernalization is cold treatment required by some plants for flowering, particularly winter varieties."
                },
                {
                    id: "b19",
                    text: "Which of the following is NOT a function of abscisic acid?",
                    options: ["Seed dormancy", "Stomatal closure", "Stress tolerance", "Cell division"],
                    answer: 3,
                    explain: "ABA inhibits cell division. Its functions include maintaining dormancy, closing stomata, and stress responses."
                },
                {
                    id: "b20",
                    text: "The arrangement of ovules in which micropyle faces the hilum is:",
                    options: ["Orthotropous", "Anatropous", "Campylotropous", "Amphitropous"],
                    answer: 0,
                    explain: "In orthotropous ovules, micropyle is at the tip opposite to hilum (straight ovule)."
                },
                {
                    id: "b21",
                    text: "Which of the following is a quantitative long day plant?",
                    options: ["Wheat", "Rice", "Soybean", "Cockscomb"],
                    answer: 0,
                    explain: "Wheat is facultative long day plant - flowers faster with long days but can flower with short days too."
                },
                {
                    id: "b22",
                    text: "The guard cells of stomata in grasses are:",
                    options: ["Kidney shaped", "Dumbbell shaped", "Spherical", "Rectangular"],
                    answer: 1,
                    explain: "Guard cells in grass family are dumbbell or bone-shaped, different from kidney-shaped in dicots."
                },
                {
                    id: "b23",
                    text: "Which of the following shows apical dominance?",
                    options: ["Auxiliary bud growth", "Terminal bud growth", "Root growth", "Flower development"],
                    answer: 1,
                    explain: "Apical dominance is suppression of lateral bud growth by terminal/apical bud through auxin."
                },
                {
                    id: "b24",
                    text: "The first product of dark reaction of photosynthesis is:",
                    options: ["3-PGA", "RuBP", "DHAP", "G3P"],
                    answer: 0,
                    explain: "3-phosphoglycerate (3-PGA) is first stable product when CO₂ combines with RuBP in Calvin cycle."
                },
                {
                    id: "b25",
                    text: "Which of the following is synthesized during G₁ phase of cell cycle?",
                    options: ["DNA", "Histone proteins", "rRNA", "Both B and C"],
                    answer: 3,
                    explain: "During G₁ phase, cells synthesize proteins including histones and rRNA, but not DNA."
                },
                {
                    id: "b26",
                    text: "The water potential of pure water is:",
                    options: ["Zero", "Positive", "Negative", "Infinite"],
                    answer: 0,
                    explain: "Pure water has water potential of zero at standard conditions, used as reference point."
                },
                {
                    id: "b27",
                    text: "Which of the following is involved in nitrogen fixation?",
                    options: ["Nitrogenase", "Nitrate reductase", "Glutamine synthetase", "Urease"],
                    answer: 0,
                    explain: "Nitrogenase enzyme complex converts atmospheric nitrogen to ammonia in nitrogen fixation."
                },
                {
                    id: "b28",
                    text: "The site of respiration in plant cells is:",
                    options: ["Chloroplasts", "Mitochondria", "Ribosomes", "Peroxisomes"],
                    answer: 1,
                    explain: "Mitochondria are sites of cellular respiration where glucose is oxidized to produce ATP."
                },
                {
                    id: "b29",
                    text: "Which of the following is characteristic of xerophytes?",
                    options: ["Reduced leaf surface", "Waxy cuticle", "Succulent stems", "All of above"],
                    answer: 3,
                    explain: "Xerophytes show multiple adaptations: reduced leaf area, thick waxy cuticle, water storage tissues."
                },
                {
                    id: "b30",
                    text: "The process of formation of wood is called:",
                    options: ["Lignification", "Suberization", "Cutinization", "Periderm formation"],
                    answer: 0,
                    explain: "Lignification is deposition of lignin in cell walls, characteristic of wood formation."
                },
                {
                    id: "b31",
                    text: "Which of the following shows determinate growth?",
                    options: ["Root", "Stem", "Leaf", "All of above"],
                    answer: 2,
                    explain: "Leaves show determinate growth - they stop growing after reaching mature size, unlike roots and stems."
                },
                {
                    id: "b32",
                    text: "The phenomenon of seed germination stimulated by red light is mediated by:",
                    options: ["Cryptochrome", "Phytochrome", "Phototropin", "UVR8"],
                    answer: 1,
                    explain: "Phytochrome system mediates red/far-red light responses including seed germination."
                },
                {
                    id: "b33",
                    text: "Which of the following is a proteolytic enzyme?",
                    options: ["Amylase", "Pepsin", "Lipase", "Catalase"],
                    answer: 1,
                    explain: "Pepsin is proteolytic enzyme that breaks down proteins in acidic conditions of stomach."
                },
                {
                    id: "b34",
                    text: "The vascular cambium is a:",
                    options: ["Primary meristem", "Secondary meristem", "Intercalary meristem", "Apical meristem"],
                    answer: 1,
                    explain: "Vascular cambium is secondary meristem responsible for secondary growth in woody plants."
                },
                {
                    id: "b35",
                    text: "Which of the following is NOT a component of phloem?",
                    options: ["Sieve tubes", "Companion cells", "Phloem parenchyma", "Vessel elements"],
                    answer: 3,
                    explain: "Vessel elements are components of xylem, not phloem. Phloem consists of sieve elements and associated cells."
                },
                {
                    id: "b36",
                    text: "The loss of turgor in guard cells leads to:",
                    options: ["Stomatal opening", "Stomatal closure", "Increased transpiration", "Photosynthesis"],
                    answer: 1,
                    explain: "Loss of turgor pressure in guard cells causes stomatal closure, reducing water loss."
                },
                {
                    id: "b37",
                    text: "Which of the following is a long day plant?",
                    options: ["Rice", "Soybean", "Spinach", "Chrysanthemum"],
                    answer: 2,
                    explain: "Spinach is a long day plant requiring long photoperiods for flowering."
                },
                {
                    id: "b38",
                    text: "The primary function of lenticels is:",
                    options: ["Water absorption", "Gas exchange", "Support", "Transport"],
                    answer: 1,
                    explain: "Lenticels are pores in bark that allow gas exchange between internal tissues and atmosphere."
                },
                {
                    id: "b39",
                    text: "Which of the following shows epigeal germination?",
                    options: ["Pea", "Bean", "Maize", "Wheat"],
                    answer: 1,
                    explain: "Bean shows epigeal germination where cotyledons emerge above ground during seedling development."
                },
                {
                    id: "b40",
                    text: "The hormone responsible for fruit ripening is:",
                    options: ["Auxin", "Gibberellin", "Ethylene", "Cytokinin"],
                    answer: 2,
                    explain: "Ethylene is gaseous hormone that triggers fruit ripening and senescence processes."
                },
                {
                    id: "b41",
                    text: "Which of the following is involved in photoperiodism?",
                    options: ["Chlorophyll", "Phytochrome", "Carotenoids", "Anthocyanins"],
                    answer: 1,
                    explain: "Phytochrome pigment system detects day length and mediates photoperiodic responses."
                },
                {
                    id: "b42",
                    text: "The cork cambium produces:",
                    options: ["Cork towards outside", "Secondary cortex towards inside", "Both A and B", "Xylem and phloem"],
                    answer: 2,
                    explain: "Cork cambium produces cork (phellem) outwardly and secondary cortex (phelloderm) inwardly."
                },
                {
                    id: "b43",
                    text: "Which of the following is NOT a micronutrient?",
                    options: ["Iron", "Boron", "Nitrogen", "Zinc"],
                    answer: 2,
                    explain: "Nitrogen is a macronutrient required in large quantities. Fe, B, and Zn are micronutrients."
                },
                {
                    id: "b44",
                    text: "The term 'bolting' refers to:",
                    options: ["Rapid stem elongation", "Flower development", "Fruit setting", "Seed dispersal"],
                    answer: 0,
                    explain: "Bolting is rapid elongation of stem prior to flowering, often induced by gibberellins."
                },
                {
                    id: "b45",
                    text: "Which of the following shows alternate phyllotaxy?",
                    options: ["Calotropis", "Guava", "China rose", "Alstonia"],
                    answer: 2,
                    explain: "China rose (Hibiscus) shows alternate phyllotaxy with single leaf at each node."
                },

                // ZOOLOGY Questions (45)
                {
                    id: "z1",
                    text: "The juxtaglomerular apparatus is involved in regulation of:",
                    options: ["Blood glucose", "Blood pressure", "Body temperature", "pH balance"],
                    answer: 1,
                    explain: "JGA secretes renin which regulates blood pressure through renin-angiotensin-aldosterone system."
                },
                {
                    id: "z2",
                    text: "Which of the following hormones is NOT produced by anterior pituitary?",
                    options: ["Growth hormone", "Prolactin", "ADH", "ACTH"],
                    answer: 2,
                    explain: "ADH (antidiuretic hormone) is produced by hypothalamus and stored in posterior pituitary."
                },
                {
                    id: "z3",
                    text: "The region of sarcomere that contains only thin filaments is:",
                    options: ["A band", "I band", "H zone", "M line"],
                    answer: 1,
                    explain: "I band contains only actin (thin) filaments and appears light under microscope."
                },
                {
                    id: "z4",
                    text: "Which of the following is NOT a function of liver?",
                    options: ["Bile production", "Glycogen storage", "Insulin production", "Detoxification"],
                    answer: 2,
                    explain: "Insulin is produced by pancreatic beta cells, not liver. Liver produces bile, stores glycogen, and detoxifies."
                },
                {
                    id: "z5",
                    text: "The corpus luteum secretes:",
                    options: ["Estrogen only", "Progesterone only", "Both estrogen and progesterone", "FSH and LH"],
                    answer: 2,
                    explain: "Corpus luteum secretes both progesterone (mainly) and estrogen to maintain pregnancy."
                },
                {
                    id: "z6",
                    text: "Which of the following is the correct sequence of cardiac cycle?",
                    options: [
                        "Systole → Diastole",
                        "Atrial systole → Ventricular systole → Diastole", 
                        "Ventricular systole → Atrial systole → Diastole",
                        "Diastole → Systole"
                    ],
                    answer: 1,
                    explain: "Cardiac cycle: Atrial contraction → Ventricular contraction → Complete heart relaxation"
                },
                {
                    id: "z7",
                    text: "The immune cells that provide immunological memory are:",
                    options: ["Plasma cells", "Memory cells", "Helper T cells", "Suppressor T cells"],
                    answer: 1,
                    explain: "Memory B and T cells provide long-term immunity and rapid response upon re-exposure."
                },
                {
                    id: "z8",
                    text: "Which part of nephron is impermeable to urea?",
                    options: ["Proximal tubule", "Loop of Henle", "Distal tubule", "Collecting duct"],
                    answer: 2,
                    explain: "Distal convoluted tubule is impermeable to urea, while collecting duct has variable permeability."
                },
                {
                    id: "z9",
                    text: "The hormone that stimulates milk production is:",
                    options: ["Oxytocin", "Prolactin", "Estrogen", "Progesterone"],
                    answer: 1,
                    explain: "Prolactin from anterior pituitary stimulates milk production in mammary glands."
                },
                {
                    id: "z10",
                    text: "Which of the following vitamins is water-soluble?",
                    options: ["Vitamin A", "Vitamin D", "Vitamin K", "Vitamin C"],
                    answer: 3,
                    explain: "Vitamin C (ascorbic acid) is water-soluble, while A, D, K are fat-soluble vitamins."
                },
                {
                    id: "z11",
                    text: "The pacemaker of heart is located in:",
                    options: ["Right atrium", "Left atrium", "Right ventricle", "Left ventricle"],
                    answer: 0,
                    explain: "Sinoatrial (SA) node, the natural pacemaker, is located in wall of right atrium."
                },
                {
                    id: "z12",
                    text: "Which enzyme converts pepsinogen to pepsin?",
                    options: ["Trypsin", "HCl", "Pepsin", "Chymotrypsin"],
                    answer: 2,
                    explain: "Pepsin itself converts pepsinogen to pepsin (autocatalytic activation) in acidic environment."
                },
                {
                    id: "z13",
                    text: "The minimum amount of urine that must be excreted per day is:",
                    options: ["500 ml", "1000 ml", "1500 ml", "2000 ml"],
                    answer: 0,
                    explain: "Obligatory urine volume is ~500 ml/day to excrete metabolic wastes at maximum concentration."
                },
                {
                    id: "z14",
                    text: "Which of the following is NOT involved in blood clotting?",
                    options: ["Fibrinogen", "Prothrombin", "Heparin", "Calcium"],
                    answer: 2,
                    explain: "Heparin is anticoagulant that prevents clotting. Fibrinogen, prothrombin, and Ca²⁺ promote clotting."
                },
                {
                    id: "z15",
                    text: "The respiratory center is located in:",
                    options: ["Cerebrum", "Cerebellum", "Medulla oblongata", "Pons"],
                    answer: 2,
                    explain: "Primary respiratory center controlling automatic breathing is in medulla oblongata."
                },
                {
                    id: "z16",
                    text: "Which of the following crosses blood-brain barrier easily?",
                    options: ["Glucose", "Insulin", "Large proteins", "Most drugs"],
                    answer: 0,
                    explain: "Glucose crosses blood-brain barrier via specific transporters to supply brain's energy needs."
                },
                {
                    id: "z17",
                    text: "The acrosome of sperm contains:",
                    options: ["Mitochondria", "Hydrolytic enzymes", "Nucleus", "Centrioles"],
                    answer: 1,
                    explain: "Acrosome contains enzymes like hyaluronidase and acrosin that help penetrate egg coats."
                },
                {
                    id: "z18",
                    text: "Which of the following hormones has shortest half-life?",
                    options: ["Insulin", "Thyroxine", "Cortisol", "Growth hormone"],
                    answer: 0,
                    explain: "Insulin has very short half-life (~5-10 minutes) allowing rapid glucose regulation."
                },
                {
                    id: "z19",
                    text: "The bicuspid valve is located between:",
                    options: [
                        "Right atrium and right ventricle",
                        "Left atrium and left ventricle", 
                        "Right ventricle and pulmonary artery",
                        "Left ventricle and aorta"
                    ],
                    answer: 1,
                    explain: "Bicuspid (mitral) valve prevents backflow from left ventricle to left atrium."
                },
                {
                    id: "z20",
                    text: "Which of the following is NOT a granulocyte?",
                    options: ["Neutrophil", "Eosinophil", "Basophil", "Monocyte"],
                    answer: 3,
                    explain: "Monocyte is agranulocyte (lacks visible granules). Neutrophils, eosinophils, basophils are granulocytes."
                },
                {
                    id: "z21",
                    text: "The hormone that regulates calcium homeostasis is:",
                    options: ["Insulin", "Calcitonin", "Parathyroid hormone", "Both B and C"],
                    answer: 3,
                    explain: "Both calcitonin (lowers Ca²⁺) and PTH (raises Ca²⁺) regulate calcium homeostasis."
                },
                {
                    id: "z22",
                    text: "Which part of brain controls body temperature?",
                    options: ["Cerebrum", "Hypothalamus", "Medulla", "Pons"],
                    answer: 1,
                    explain: "Hypothalamus contains thermoregulatory center that maintains body temperature."
                },
                {
                    id: "z23",
                    text: "The sliding filament theory explains:",
                    options: ["Muscle contraction", "Nerve conduction", "Blood clotting", "Hormone action"],
                    answer: 0,
                    explain: "Sliding filament theory explains muscle contraction through actin-myosin interaction."
                },
                {
                    id: "z24",
                    text: "Which of the following has highest concentration in urine?",
                    options: ["Glucose", "Urea", "Proteins", "RBCs"],
                    answer: 1,
                    explain: "Urea is major nitrogenous waste with highest concentration in normal urine."
                },
                {
                    id: "z25",
                    text: "The cells that secrete insulin are:",
                    options: ["Alpha cells", "Beta cells", "Delta cells", "Gamma cells"],
                    answer: 1,
                    explain: "Beta cells in pancreatic islets of Langerhans secrete insulin hormone."
                },
                {
                    id: "z26",
                    text: "Which cranial nerve controls heart rate?",
                    options: ["Vagus", "Trigeminal", "Facial", "Glossopharyngeal"],
                    answer: 0,
                    explain: "Vagus nerve (X cranial nerve) provides parasympathetic innervation to heart, slowing heart rate."
                },
                {
                    id: "z27",
                    text: "The functional unit of muscle contraction is:",
                    options: ["Myofibril", "Sarcomere", "Muscle fiber", "Fascicle"],
                    answer: 1,
                    explain: "Sarcomere is the basic contractile unit between two Z-lines containing actin and myosin."
                },
                {
                    id: "z28",
                    text: "Which of the following is reabsorbed maximally in PCT?",
                    options: ["Water", "Glucose", "Sodium", "All of above"],
                    answer: 3,
                    explain: "Proximal convoluted tubule reabsorbs ~65% of filtered water, sodium, and 100% of glucose."
                },
                {
                    id: "z29",
                    text: "The period of cardiac cycle when all valves are closed is:",
                    options: ["Atrial systole", "Ventricular systole", "Isovolumetric contraction", "Complete diastole"],
                    answer: 2,
                    explain: "During isovolumetric contraction, all valves are closed and ventricular pressure builds up."
                },
                {
                    id: "z30",
                    text: "Which of the following is NOT an accessory reproductive gland in males?",
                    options: ["Prostate", "Seminal vesicles", "Bulbourethral gland", "Epididymis"],
                    answer: 3,
                    explain: "Epididymis is part of male reproductive tract, not an accessory gland. It stores and matures sperm."
                },
                {
                    id: "z31",
                    text: "The oxygen-carrying capacity of blood is determined by:",
                    options: ["Plasma proteins", "Hemoglobin concentration", "WBC count", "Platelet count"],
                    answer: 1,
                    explain: "Hemoglobin concentration determines oxygen-carrying capacity as each Hb molecule binds 4 O₂."
                },
                {
                    id: "z32",
                    text: "Which of the following occurs during inspiration?",
                    options: ["Diaphragm relaxes", "Ribs move downward", "Lung volume increases", "Intrapulmonary pressure increases"],
                    answer: 2,
                    explain: "During inspiration, diaphragm contracts, ribs move up and out, increasing lung volume."
                },
                {
                    id: "z33",
                    text: "The condition where lens loses its elasticity is:",
                    options: ["Myopia", "Hyperopia", "Presbyopia", "Astigmatism"],
                    answer: 2,
                    explain: "Presbyopia is age-related loss of lens elasticity, reducing accommodation power."
                },
                {
                    id: "z34",
                    text: "Which hormone is called 'fight or flight' hormone?",
                    options: ["Insulin", "Adrenaline", "Thyroxine", "Cortisol"],
                    answer: 1,
                    explain: "Adrenaline (epinephrine) prepares body for fight or flight response during stress."
                },
                {
                    id: "z35",
                    text: "The glomerular filtration rate in humans is approximately:",
                    options: ["125 ml/min", "180 ml/min", "200 ml/min", "150 ml/min"],
                    answer: 0,
                    explain: "Normal GFR is ~125 ml/min or ~180 L/day, indicating kidney function."
                },
                {
                    id: "z36",
                    text: "Which of the following is involved in acquired immunity?",
                    options: ["Skin", "Mucus", "Antibodies", "Lysozyme"],
                    answer: 2,
                    explain: "Antibodies provide specific acquired immunity, while others are part of innate immunity."
                },
                {
                    id: "z37",
                    text: "The site of fertilization in human female is:",
                    options: ["Ovary", "Cervix", "Ampulla of fallopian tube", "Uterus"],
                    answer: 2,
                    explain: "Fertilization typically occurs in ampulla region of fallopian tube."
                },
                {
                    id: "z38",
                    text: "Which of the following is NOT a function of kidney?",
                    options: ["Osmoregulation", "pH regulation", "RBC production", "Glucose production"],
                    answer: 2,
                    explain: "RBCs are produced in bone marrow. Kidneys produce erythropoietin hormone that stimulates RBC production."
                },
                {
                    id: "z39",
                    text: "The neurotransmitter released at all autonomic ganglia is:",
                    options: ["Adrenaline", "Noradrenaline", "Acetylcholine", "Dopamine"],
                    answer: 2,
                    explain: "Acetylcholine is released at all autonomic ganglia (both sympathetic and parasympathetic)."
                },
                {
                    id: "z40",
                    text: "Which of the following increases during exercise?",
                    options: ["Stroke volume", "Heart rate", "Cardiac output", "All of above"],
                    answer: 3,
                    explain: "During exercise, heart rate, stroke volume, and cardiac output all increase to meet oxygen demand."
                },
                {
                    id: "z41",
                    text: "The hormone responsible for milk ejection reflex is:",
                    options: ["Prolactin", "Oxytocin", "Estrogen", "Progesterone"],
                    answer: 1,
                    explain: "Oxytocin causes contraction of mammary alveoli, leading to milk ejection (let-down reflex)."
                },
                {
                    id: "z42",
                    text: "Which of the following is characteristic of cardiac muscle?",
                    options: ["Voluntary", "Multinucleated", "Intercalated discs", "Spindle shaped"],
                    answer: 2,
                    explain: "Cardiac muscle has intercalated discs containing gap junctions for synchronized contraction."
                },
                {
                    id: "z43",
                    text: "The blind spot in retina is due to:",
                    options: ["Absence of blood vessels", "Presence of optic nerve", "Absence of photoreceptors", "Both B and C"],
                    answer: 3,
                    explain: "Blind spot is where optic nerve exits retina - no photoreceptors present there."
                },
                {
                    id: "z44",
                    text: "Which of the following is NOT a sexually transmitted infection?",
                    options: ["Gonorrhea", "Syphilis", "Hepatitis B", "Tuberculosis"],
                    answer: 3,
                    explain: "Tuberculosis is airborne disease transmitted through respiratory droplets, not sexual contact."
                },
                {
                    id: "z45",
                    text: "The gestation period in humans is:",
                    options: ["36 weeks", "38 weeks", "40 weeks", "42 weeks"],
                    answer: 2,
                    explain: "Human gestation period is approximately 40 weeks (280 days) from last menstrual period."
                }
            ]
        }
    ]
};
