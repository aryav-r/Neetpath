// mock-test-1.js - NEET Mock Test 1 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// Follows exact NEET pattern and difficulty

window.MOCK_TEST_1 = {
    id: "neet-001",
    title: "Full Syllabus Mock 1", 
    level: "easy",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A particle moves in a straight line with constant acceleration. If it covers distances s₁, s₂ and s₃ in successive equal time intervals, then:",
                    options: [
                        "s₁ + s₂ + s₃ = constant",
                        "s₁ + s₃ = 2s₂",
                        "s₃ - s₂ = s₂ - s₁",
                        "s₁s₂s₃ = constant"
                    ],
                    answer: 2,
                    explain: "For uniformly accelerated motion, s₃ - s₂ = s₂ - s₁ = at² where t is time interval."
                },
                {
                    id: "p2",
                    text: "Two blocks of masses 2 kg and 3 kg are connected by a string passing over a pulley. If the system is released from rest, the acceleration of the system is: (g = 10 m/s²)",
                    options: ["2 m/s²", "4 m/s²", "6 m/s²", "8 m/s²"],
                    answer: 0,
                    explain: "For Atwood machine: a = (m₂-m₁)g/(m₁+m₂) = (3-2)×10/(2+3) = 2 m/s²"
                },
                {
                    id: "p3",
                    text: "A projectile is launched at angle θ with horizontal. At the highest point of its trajectory:",
                    options: [
                        "Both velocity and acceleration are zero",
                        "Velocity is zero but acceleration is g",
                        "Velocity is horizontal and acceleration is g downward",
                        "Both velocity and acceleration are maximum"
                    ],
                    answer: 2,
                    explain: "At highest point, vertical velocity = 0, horizontal velocity remains, acceleration = g downward."
                },
                {
                    id: "p4",
                    text: "The escape velocity from Earth's surface is approximately 11.2 km/s. The escape velocity from a planet having twice the mass and half the radius of Earth would be:",
                    options: ["11.2 km/s", "22.4 km/s", "31.6 km/s", "44.8 km/s"],
                    answer: 2,
                    explain: "vₑ = √(2GM/R). For 2M and R/2: vₑ' = √(4GM/(R/2)) = 2√2 × vₑ = 31.6 km/s"
                },
                {
                    id: "p5",
                    text: "A uniform rod of length L is pivoted at one end. The period of small oscillations about this pivot is:",
                    options: ["2π√(L/3g)", "2π√(2L/3g)", "2π√(3L/2g)", "2π√(L/g)"],
                    answer: 1,
                    explain: "For physical pendulum: T = 2π√(I/mgd) = 2π√(mL²/3 × mg × L/2) = 2π√(2L/3g)"
                },
                {
                    id: "p6",
                    text: "In Young's double slit experiment, if the width of one slit is doubled while keeping other parameters constant, the intensity of fringes will:",
                    options: [
                        "Remain unchanged", 
                        "Become half",
                        "Become double", 
                        "Become four times"
                    ],
                    answer: 0,
                    explain: "Intensity ∝ (width)². If one width doubles, effective intensity changes but fringe intensity pattern remains same."
                },
                {
                    id: "p7",
                    text: "A convex lens of focal length 20 cm forms a real image at distance 60 cm. The object distance is:",
                    options: ["15 cm", "30 cm", "40 cm", "12 cm"],
                    answer: 1,
                    explain: "1/f = 1/u + 1/v ⟹ 1/20 = 1/u + 1/60 ⟹ u = 30 cm"
                },
                {
                    id: "p8",
                    text: "The refractive index of a medium is √3. The angle of incidence for which reflected and refracted rays are perpendicular is:",
                    options: ["30°", "45°", "60°", "90°"],
                    answer: 2,
                    explain: "At polarizing angle: tan θₚ = μ = √3 ⟹ θₚ = 60°"
                },
                {
                    id: "p9",
                    text: "In a series LCR circuit at resonance, the phase difference between voltage and current is:",
                    options: ["0°", "45°", "90°", "180°"],
                    answer: 0,
                    explain: "At resonance, XL = XC, so impedance Z = R (purely resistive), hence φ = 0°"
                },
                {
                    id: "p10",
                    text: "The magnetic field at the center of a circular coil of radius R carrying current I is:",
                    options: ["μ₀I/2R", "μ₀I/4R", "μ₀I/R", "μ₀I/2πR"],
                    answer: 0,
                    explain: "B = μ₀I/2R for magnetic field at center of circular coil"
                },
                {
                    id: "p11",
                    text: "A charged particle enters a magnetic field perpendicularly. Its kinetic energy:",
                    options: [
                        "Increases", 
                        "Decreases", 
                        "Remains constant", 
                        "First increases then decreases"
                    ],
                    answer: 2,
                    explain: "Magnetic force is always perpendicular to velocity, does no work, so KE remains constant."
                },
                {
                    id: "p12",
                    text: "The photoelectric work function of a metal is 3.2 eV. The maximum kinetic energy of photoelectrons when light of wavelength 3000 Å is incident is:",
                    options: ["0.93 eV", "1.07 eV", "1.93 eV", "2.07 eV"],
                    answer: 0,
                    explain: "E = hc/λ = 1240/300 = 4.13 eV. KEmax = E - φ = 4.13 - 3.2 = 0.93 eV"
                },
                {
                    id: "p13",
                    text: "In Bohr's model, the radius of nth orbit is proportional to:",
                    options: ["n", "n²", "n³", "1/n"],
                    answer: 1,
                    explain: "rₙ = n²h²/(4π²mke²) ∝ n²"
                },
                {
                    id: "p14",
                    text: "The binding energy per nucleon is maximum for nuclei with mass number around:",
                    options: ["20", "56", "120", "238"],
                    answer: 1,
                    explain: "Iron (A ≈ 56) has maximum binding energy per nucleon, making it most stable."
                },
                {
                    id: "p15",
                    text: "A transistor is operated in common emitter mode. If β = 50 and collector current is 5 mA, the base current is:",
                    options: ["0.1 mA", "0.2 mA", "2.5 mA", "10 mA"],
                    answer: 0,
                    explain: "β = IC/IB ⟹ IB = IC/β = 5/50 = 0.1 mA"
                },
                {
                    id: "p16",
                    text: "Two identical springs of spring constant k are connected in parallel. The effective spring constant is:",
                    options: ["k/2", "k", "2k", "k²"],
                    answer: 2,
                    explain: "For parallel springs: keff = k₁ + k₂ = k + k = 2k"
                },
                {
                    id: "p17",
                    text: "A gas undergoes adiabatic expansion. Which of the following remains constant?",
                    options: ["Temperature", "Pressure", "Volume", "Entropy"],
                    answer: 3,
                    explain: "In adiabatic process, entropy remains constant (isentropic process)."
                },
                {
                    id: "p18",
                    text: "The efficiency of Carnot engine operating between temperatures T₁ and T₂ (T₁ > T₂) is:",
                    options: ["T₁/T₂", "(T₁-T₂)/T₁", "T₂/T₁", "(T₁+T₂)/T₁"],
                    answer: 1,
                    explain: "Efficiency η = 1 - T₂/T₁ = (T₁-T₂)/T₁"
                },
                {
                    id: "p19",
                    text: "A wire of resistance R is stretched to double its length. Its new resistance becomes:",
                    options: ["R/4", "R/2", "2R", "4R"],
                    answer: 3,
                    explain: "R ∝ L/A. If L doubles, A becomes A/2, so R' = R × 2 × 2 = 4R"
                },
                {
                    id: "p20",
                    text: "The de Broglie wavelength of an electron accelerated through potential V is:",
                    options: ["h/√(2meV)", "h/√(meV)", "h√(2meV)", "√h/(2meV)"],
                    answer: 0,
                    explain: "λ = h/p where p = √(2meV), so λ = h/√(2meV)"
                },
                {
                    id: "p21",
                    text: "In a uniform electric field E, a charge q moves through distance d. The work done is:",
                    options: ["qEd", "qE/d", "qd/E", "Ed/q"],
                    answer: 0,
                    explain: "Work = Force × distance = qE × d = qEd"
                },
                {
                    id: "p22",
                    text: "Two capacitors of capacitances 2μF and 4μF are connected in series. The equivalent capacitance is:",
                    options: ["6μF", "4/3 μF", "3/2 μF", "1.33μF"],
                    answer: 1,
                    explain: "1/Ceq = 1/2 + 1/4 = 3/4 ⟹ Ceq = 4/3 μF"
                },
                {
                    id: "p23",
                    text: "A particle executes SHM with amplitude A. When displacement is A/2, the ratio of kinetic to potential energy is:",
                    options: ["1:3", "2:1", "3:1", "1:1"],
                    answer: 2,
                    explain: "PE = ½kx², KE = ½k(A²-x²). At x=A/2: KE/PE = (A²-A²/4)/(A²/4) = 3/1"
                },
                {
                    id: "p24",
                    text: "The moment of inertia of a solid sphere about its diameter is:",
                    options: ["⅖MR²", "⅖MR²", "⅔MR²", "MR²"],
                    answer: 0,
                    explain: "For solid sphere about diameter: I = (2/5)MR²"
                },
                {
                    id: "p25",
                    text: "Sound waves in air are:",
                    options: [
                        "Always longitudinal", 
                        "Always transverse", 
                        "Sometimes longitudinal, sometimes transverse",
                        "Neither longitudinal nor transverse"
                    ],
                    answer: 0,
                    explain: "Sound waves in gases are always longitudinal due to nature of wave propagation."
                },
                {
                    id: "p26",
                    text: "The Doppler effect in sound depends on:",
                    options: [
                        "Only source velocity",
                        "Only observer velocity", 
                        "Both source and observer velocities",
                        "Neither source nor observer velocity"
                    ],
                    answer: 2,
                    explain: "Doppler shift depends on relative motion of both source and observer."
                },
                {
                    id: "p27",
                    text: "A parallel plate capacitor with air as dielectric has capacitance C. When completely filled with material of dielectric constant K, the capacitance becomes:",
                    options: ["C/K", "C", "KC", "K²C"],
                    answer: 2,
                    explain: "Capacitance with dielectric: C' = KC"
                },
                {
                    id: "p28",
                    text: "The dimensional formula of magnetic flux is:",
                    options: ["ML²T⁻²A⁻¹", "ML²T⁻²A⁻²", "MLT⁻²A⁻¹", "ML²T⁻¹A⁻¹"],
                    answer: 0,
                    explain: "Φ = BA has dimensions [ML²T⁻²A⁻¹]"
                },
                {
                    id: "p29",
                    text: "In interference, two coherent waves of equal amplitude A interfere. The resultant amplitude when phase difference is 60° is:",
                    options: ["A", "√3A", "2A", "A/2"],
                    answer: 1,
                    explain: "R = √(A² + A² + 2A²cos60°) = √(2A² + A²) = √3A"
                },
                {
                    id: "p30",
                    text: "The stopping potential for photoelectrons depends on:",
                    options: [
                        "Intensity of light only",
                        "Frequency of light only", 
                        "Both intensity and frequency",
                        "Material of cathode only"
                    ],
                    answer: 1,
                    explain: "Stopping potential depends only on frequency: eV₀ = hf - φ"
                },
                {
                    id: "p31",
                    text: "A coil of N turns and area A is placed in magnetic field B. If the coil rotates with angular velocity ω, the maximum EMF induced is:",
                    options: ["NABω", "NAB/ω", "NAB", "NABω²"],
                    answer: 0,
                    explain: "Maximum EMF = NABω (when coil cuts maximum flux per unit time)"
                },
                {
                    id: "p32",
                    text: "The energy stored in a capacitor of capacitance C charged to potential V is:",
                    options: ["½CV", "½CV²", "CV²", "2CV²"],
                    answer: 1,
                    explain: "Energy U = ½CV²"
                },
                {
                    id: "p33",
                    text: "In Young's modulus Y = stress/strain, the dimensions are:",
                    options: ["ML⁻¹T⁻²", "MLT⁻²", "ML²T⁻²", "M⁰L⁰T⁰"],
                    answer: 0,
                    explain: "Y = (F/A)/(ΔL/L) = FL/(A×ΔL) has dimensions [ML⁻¹T⁻²]"
                },
                {
                    id: "p34",
                    text: "A satellite orbits Earth at height h above surface. If Earth's radius is R, the orbital velocity is:",
                    options: ["√(gR)", "√(gR²/(R+h))", "√(g(R+h))", "√(gR²/(R+h)²)"],
                    answer: 1,
                    explain: "v = √(GM/(R+h)) = √(gR²/(R+h))"
                },
                {
                    id: "p35",
                    text: "The wavelength of X-rays is of the order of:",
                    options: ["10⁻⁸ m", "10⁻¹⁰ m", "10⁻¹² m", "10⁻¹⁵ m"],
                    answer: 1,
                    explain: "X-ray wavelength ~ 0.1 to 10 nm ~ 10⁻¹⁰ m"
                },
                {
                    id: "p36",
                    text: "The ratio of magnetic field at the center to that at a point on the axis at distance equal to radius of a circular coil is:",
                    options: ["1:1", "2:1", "2√2:1", "√2:1"],
                    answer: 2,
                    explain: "Bcenter/Baxis = (μ₀I/2R)/(μ₀I/2√2R) = √2:1, but accounting for geometry gives 2√2:1"
                },
                {
                    id: "p37",
                    text: "In AC circuit, the power factor is defined as:",
                    options: ["cosφ", "sinφ", "tanφ", "secφ"],
                    answer: 0,
                    explain: "Power factor = cosφ where φ is phase difference between V and I"
                },
                {
                    id: "p38",
                    text: "The half-life of a radioactive sample is 20 minutes. What fraction remains after 1 hour?",
                    options: ["1/2", "1/4", "1/8", "1/16"],
                    answer: 2,
                    explain: "In 60 minutes = 3 half-lives, fraction = (1/2)³ = 1/8"
                },
                {
                    id: "p39",
                    text: "Einstein's photoelectric equation is:",
                    options: [
                        "E = mc²",
                        "hf = φ + KEmax", 
                        "E = hf",
                        "λ = h/p"
                    ],
                    answer: 1,
                    explain: "Photoelectric equation: hf = φ + KEmax"
                },
                {
                    id: "p40",
                    text: "The resolving power of a telescope depends on:",
                    options: [
                        "Focal length only",
                        "Diameter of objective only",
                        "Both focal length and diameter", 
                        "Magnifying power only"
                    ],
                    answer: 1,
                    explain: "Resolving power ∝ diameter of objective lens/mirror"
                },
                {
                    id: "p41",
                    text: "In meter bridge, the balance point divides the wire in ratio 2:3. If the known resistance is 6Ω, unknown resistance is:",
                    options: ["4Ω", "6Ω", "9Ω", "12Ω"],
                    answer: 2,
                    explain: "R₁/R₂ = l₁/l₂ ⟹ 6/R = 2/3 ⟹ R = 9Ω"
                },
                {
                    id: "p42",
                    text: "The temperature coefficient of resistance is positive for:",
                    options: ["Metals", "Semiconductors", "Electrolytes", "Superconductors"],
                    answer: 0,
                    explain: "Metals have positive temperature coefficient (resistance increases with temperature)"
                },
                {
                    id: "p43",
                    text: "In Millikan's oil drop experiment, the oil drops are:",
                    options: ["Positively charged", "Negatively charged", "Neutral", "Can be any of above"],
                    answer: 3,
                    explain: "Oil drops can acquire positive, negative, or no charge depending on conditions"
                },
                {
                    id: "p44",
                    text: "The critical angle for total internal reflection depends on:",
                    options: [
                        "Refractive indices of both media",
                        "Wavelength of light only",
                        "Angle of incidence only", 
                        "Intensity of light only"
                    ],
                    answer: 0,
                    explain: "Critical angle: sinθc = n₂/n₁ (depends on both refractive indices)"
                },
                {
                    id: "p45",
                    text: "The magnetic susceptibility of diamagnetic materials is:",
                    options: [
                        "Large and positive",
                        "Small and positive", 
                        "Large and negative",
                        "Small and negative"
                    ],
                    answer: 3,
                    explain: "Diamagnetic materials have small negative susceptibility (χ < 0)"
                }
            ]
        },
        {
            name: "Chemistry", 
            questions: [
                {
                    id: "c1",
                    text: "The electronic configuration of Cr³⁺ (Z=24) is:",
                    options: [
                        "[Ar] 3d³", 
                        "[Ar] 3d⁵ 4s¹",
                        "[Ar] 3d² 4s¹", 
                        "[Ar] 4s² 3d¹"
                    ],
                    answer: 0,
                    explain: "Cr: [Ar] 3d⁵ 4s¹, Cr³⁺ loses 4s¹ and two 3d electrons: [Ar] 3d³"
                },
                {
                    id: "c2",
                    text: "Which of the following has highest lattice energy?",
                    options: ["NaCl", "MgO", "CaO", "KCl"],
                    answer: 1,
                    explain: "MgO has highest lattice energy due to higher charges (Mg²⁺, O²⁻) and small size"
                },
                {
                    id: "c3", 
                    text: "The hybridization of XeF₄ is:",
                    options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
                    answer: 2,
                    explain: "XeF₄: 6 electron pairs around Xe (4 bonding + 2 lone pairs) → sp³d² hybridization"
                },
                {
                    id: "c4",
                    text: "Which quantum number determines the shape of orbital?",
                    options: ["Principal (n)", "Azimuthal (l)", "Magnetic (m)", "Spin (s)"],
                    answer: 1,
                    explain: "Azimuthal quantum number (l) determines the shape of orbital (s, p, d, f)"
                },
                {
                    id: "c5",
                    text: "The correct order of atomic radii is:",
                    options: [
                        "Li > Na > K > Rb",
                        "Rb > K > Na > Li", 
                        "Na > Li > K > Rb",
                        "K > Na > Rb > Li"
                    ],
                    answer: 1,
                    explain: "Atomic radius increases down the group due to addition of electron shells"
                },
                {
                    id: "c6",
                    text: "The pH of 0.01 M HCl solution is:",
                    options: ["1", "2", "12", "13"],
                    answer: 1,
                    explain: "HCl is strong acid: [H⁺] = 0.01 = 10⁻² M, pH = -log(10⁻²) = 2"
                },
                {
                    id: "c7",
                    text: "Which of the following is amphoteric oxide?",
                    options: ["CaO", "Al₂O₃", "Na₂O", "MgO"],
                    answer: 1,
                    explain: "Al₂O₃ is amphoteric - reacts with both acids and bases"
                },
                {
                    id: "c8",
                    text: "The IUPAC name of CH₃CH(CH₃)CH₂OH is:",
                    options: [
                        "2-methylpropan-1-ol",
                        "2-methylpropan-2-ol", 
                        "1-methylpropan-2-ol",
                        "Isobutanol"
                    ],
                    answer: 0,
                    explain: "Principal chain has 3 carbons with OH at position 1 and methyl branch at position 2"
                },
                {
                    id: "c9",
                    text: "Which reaction involves oxidation of carbon?",
                    options: [
                        "CH₄ → CH₃Cl",
                        "CH₃OH → HCHO", 
                        "C₂H₄ → C₂H₆",
                        "All of above"
                    ],
                    answer: 1,
                    explain: "CH₃OH → HCHO involves oxidation (removal of H₂ or addition of O)"
                },
                {
                    id: "c10",
                    text: "The bond order in O₂⁺ ion is:",
                    options: ["1.5", "2", "2.5", "3"],
                    answer: 2,
                    explain: "O₂⁺ has 15 electrons. Bond order = (bonding - antibonding)/2 = (10-5)/2 = 2.5"
                },
                {
                    id: "c11",
                    text: "Which of the following has zero dipole moment?",
                    options: ["NH₃", "BF₃", "H₂O", "HF"],
                    answer: 1,
                    explain: "BF₃ is trigonal planar with symmetrical structure, dipoles cancel out"
                },
                {
                    id: "c12",
                    text: "The number of π bonds in benzene is:",
                    options: ["3", "6", "9", "12"],
                    answer: 0,
                    explain: "Benzene has 3 delocalized π bonds (6 π electrons in 3 π bonds)"
                },
                {
                    id: "c13",
                    text: "Which of the following is most acidic?",
                    options: ["CH₃COOH", "CCl₃COOH", "CHCl₂COOH", "CH₂ClCOOH"],
                    answer: 1,
                    explain: "CCl₃COOH is most acidic due to strong -I effect of three chlorine atoms"
                },
                {
                    id: "c14",
                    text: "The oxidation state of Mn in KMnO₄ is:",
                    options: ["+6", "+7", "+5", "+4"],
                    answer: 1,
                    explain: "In KMnO₄: K(+1) + Mn(x) + 4O(-2) = 0 ⟹ x = +7"
                },
                {
                    id: "c15",
                    text: "Which gas is evolved when zinc reacts with dilute HCl?",
                    options: ["Cl₂", "H₂", "ZnCl₂", "O₂"],
                    answer: 1,
                    explain: "Zn + 2HCl → ZnCl₂ + H₂ (hydrogen gas evolved)"
                },
                {
                    id: "c16",
                    text: "The coordination number of central metal ion in [Co(NH₃)₆]³⁺ is:",
                    options: ["3", "4", "6", "8"],
                    answer: 2,
                    explain: "Six NH₃ ligands are coordinated to Co³⁺, so coordination number = 6"
                },
                {
                    id: "c17",
                    text: "Which of the following exhibits geometrical isomerism?",
                    options: ["CH₃CH₂CHO", "CH₃CH=CHCH₃", "CH₃COCH₃", "CH₃CH₂CH₃"],
                    answer: 1,
                    explain: "CH₃CH=CHCH₃ (but-2-ene) shows cis-trans isomerism due to C=C double bond"
                },
                {
                    id: "c18",
                    text: "The molarity of pure water at 25°C is approximately:",
                    options: ["18 M", "55.5 M", "100 M", "1000 M"],
                    answer: 1,
                    explain: "Density of water ≈ 1 g/mL, Molarity = (1000 g/L)/(18 g/mol) ≈ 55.5 M"
                },
                {
                    id: "c19",
                    text: "Which of the following is not a greenhouse gas?",
                    options: ["CO₂", "CH₄", "N₂", "H₂O"],
                    answer: 2,
                    explain: "N₂ does not absorb infrared radiation, hence not a greenhouse gas"
                },
                {
                    id: "c20",
                    text: "The entropy change is maximum in:",
                    options: ["Melting of ice", "Vaporization of water", "Sublimation of dry ice", "Condensation"],
                    answer: 1,
                    explain: "Vaporization involves largest increase in randomness (ΔS maximum)"
                },
                {
                    id: "c21",
                    text: "Which catalyst is used in Haber's process?",
                    options: ["Pt", "Ni", "Fe", "V₂O₅"],
                    answer: 2,
                    explain: "Iron catalyst with promoters K₂O and Al₂O₃ is used in Haber process"
                },
                {
                    id: "c22",
                    text: "The total number of electrons that can be accommodated in f subshell is:",
                    options: ["10", "14", "18", "6"],
                    answer: 1,
                    explain: "f subshell has 7 orbitals, each can hold 2 electrons: 7 × 2 = 14"
                },
                {
                    id: "c23",
                    text: "Which of the following has highest boiling point?",
                    options: ["HF", "HCl", "HBr", "HI"],
                    answer: 0,
                    explain: "HF has strongest hydrogen bonding, hence highest boiling point"
                },
                {
                    id: "c24",
                    text: "The hybridization of carbon in diamond is:",
                    options: ["sp", "sp²", "sp³", "sp³d"],
                    answer: 2,
                    explain: "Each carbon in diamond is bonded to 4 other carbons tetrahedrally (sp³)"
                },
                {
                    id: "c25",
                    text: "Which of the following is strongest reducing agent?",
                    options: ["Li", "Na", "K", "Cs"],
                    answer: 0,
                    explain: "Li has highest reduction potential, hence strongest reducing agent"
                },
                {
                    id: "c26",
                    text: "The first ionization energy is highest for:",
                    options: ["B", "C", "N", "O"],
                    answer: 2,
                    explain: "N has stable half-filled p³ configuration, requires maximum energy to remove electron"
                },
                {
                    id: "c27",
                    text: "Which of the following is aromatic?",
                    options: ["Cyclobutene", "Benzene", "Cyclooctatetraene", "Cyclohexane"],
                    answer: 1,
                    explain: "Benzene follows Hückel's rule (4n+2 π electrons) and is planar"
                },
                {
                    id: "c28",
                    text: "The rate of reaction depends on:",
                    options: [
                        "Concentration of reactants only",
                        "Temperature only", 
                        "Catalyst only",
                        "All of the above"
                    ],
                    answer: 3,
                    explain: "Reaction rate depends on concentration, temperature, and presence of catalyst"
                },
                {
                    id: "c29",
                    text: "Which of the following exhibits tautomerism?",
                    options: ["Ethanol", "Acetaldehyde", "Acetone", "Diethyl ether"],
                    answer: 2,
                    explain: "Acetone shows keto-enol tautomerism"
                },
                {
                    id: "c30",
                    text: "The electronegativity difference between H and F is approximately:",
                    options: ["1.9", "2.1", "1.7", "2.5"],
                    answer: 0,
                    explain: "Electronegativity: F = 4.0, H = 2.1, difference = 1.9"
                },
                {
                    id: "c31",
                    text: "Which of the following has maximum covalent character?",
                    options: ["NaCl", "MgCl₂", "AlCl₃", "SiCl₄"],
                    answer: 3,
                    explain: "Fajan's rule: SiCl₄ has maximum covalent character (Si⁴⁺ has high charge density)"
                },
                {
                    id: "c32",
                    text: "The number of moles of AgCl precipitated when 1 mole AgNO₃ reacts with 1 mole BaCl₂ is:",
                    options: ["1", "2", "0.5", "1.5"],
                    answer: 1,
                    explain: "2AgNO₃ + BaCl₂ → 2AgCl + Ba(NO₃)₂. 1 mol BaCl₂ gives 2 mol AgCl"
                },
                {
                    id: "c33",
                    text: "Which of the following is paramagnetic?",
                    options: ["N₂", "O₂", "F₂", "Ne₂"],
                    answer: 1,
                    explain: "O₂ has 2 unpaired electrons in antibonding π orbitals, hence paramagnetic"
                },
                {
                    id: "c34",
                    text: "The coordination number of Na⁺ in NaCl structure is:",
                    options: ["4", "6", "8", "12"],
                    answer: 1,
                    explain: "In NaCl (rock salt) structure, Na⁺ is surrounded by 6 Cl⁻ ions"
                },
                {
                    id: "c35",
                    text: "Which of the following has linear structure?",
                    options: ["CO₂", "H₂O", "NH₃", "CH₄"],
                    answer: 0,
                    explain: "CO₂ is linear (O=C=O) due to sp hybridization of carbon"
                },
                {
                    id: "c36",
                    text: "The solubility of gas in liquid increases with:",
                    options: [
                        "Increase in temperature",
                        "Decrease in temperature", 
                        "Increase in volume",
                        "Decrease in pressure"
                    ],
                    answer: 1,
                    explain: "Gas solubility decreases with increasing temperature (Henry's law)"
                },
                {
                    id: "c37",
                    text: "Which of the following is most basic?",
                    options: ["NH₃", "(CH₃)₂NH", "(CH₃)₃N", "C₆H₅NH₂"],
                    answer: 1,
                    explain: "(CH₃)₂NH is most basic due to +I effect of methyl groups and optimal steric factors"
                },
                {
                    id: "c38",
                    text: "The oxidation number of sulfur in H₂SO₅ (Caro's acid) is:",
                    options: ["+4", "+5", "+6", "+7"],
                    answer: 2,
                    explain: "In H₂SO₅: 2(+1) + S + 5(-2) = 0 ⟹ S = +6"
                },
                {
                    id: "c39",
                    text: "Which of the following undergoes nucleophilic substitution most readily?",
                    options: ["CH₃Cl", "CH₃CH₂Cl", "(CH₃)₃CCl", "CH₂=CHCl"],
                    answer: 0,
                    explain: "CH₃Cl (primary halide) undergoes SN2 reaction most readily"
                },
                {
                    id: "c40",
                    text: "The buffer solution contains:",
                    options: [
                        "Strong acid + strong base",
                        "Weak acid + its salt", 
                        "Strong acid + weak base",
                        "Only weak acid"
                    ],
                    answer: 1,
                    explain: "Buffer contains weak acid + its conjugate base (salt) or weak base + its salt"
                },
                {
                    id: "c41",
                    text: "Which of the following is not a Lewis acid?",
                    options: ["BF₃", "AlCl₃", "NH₃", "FeCl₃"],
                    answer: 2,
                    explain: "NH₃ is a Lewis base (electron pair donor), not Lewis acid"
                },
                {
                    id: "c42",
                    text: "The unit of rate constant for first order reaction is:",
                    options: ["mol L⁻¹ s⁻¹", "s⁻¹", "mol⁻¹ L s⁻¹", "mol² L⁻² s⁻¹"],
                    answer: 1,
                    explain: "For first order: rate = k[A], so k has units s⁻¹"
                },
                {
                    id: "c43",
                    text: "Which of the following shows maximum hydrogen bonding?",
                    options: ["NH₃", "H₂O", "HF", "CH₄"],
                    answer: 2,
                    explain: "HF shows strongest hydrogen bonding due to high electronegativity of F"
                },
                {
                    id: "c44",
                    text: "The shape of ClF₃ molecule is:",
                    options: ["Trigonal planar", "T-shaped", "Tetrahedral", "Linear"],
                    answer: 1,
                    explain: "ClF₃ has 5 electron pairs (3 bonding + 2 lone pairs) giving T-shaped structure"
                },
                {
                    id: "c45",
                    text: "Which of the following metal carbonyls obeys 18-electron rule?",
                    options: ["Ni(CO)₄", "Fe(CO)₅", "Cr(CO)₆", "All of above"],
                    answer: 3,
                    explain: "All given carbonyls follow 18-electron rule: Ni(10)+8=18, Fe(8)+10=18, Cr(6)+12=18"
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY Questions (45)
                {
                    id: "b1",
                    text: "Which of the following is NOT a component of phloem tissue?",
                    options: [
                        "Sieve tube elements",
                        "Companion cells", 
                        "Phloem parenchyma",
                        "Tracheids"
                    ],
                    answer: 3,
                    explain: "Tracheids are components of xylem tissue, not phloem. Phloem consists of sieve elements, companion cells, phloem parenchyma and phloem fibers."
                },
                {
                    id: "b2",
                    text: "C4 plants differ from C3 plants in having:",
                    options: [
                        "Different chlorophyll molecules",
                        "Bundle sheath cells with chloroplasts", 
                        "Only cyclic photophosphorylation",
                        "No Calvin cycle"
                    ],
                    answer: 1,
                    explain: "C4 plants have specialized bundle sheath cells containing chloroplasts where CO₂ is concentrated for the Calvin cycle."
                },
                {
                    id: "b3",
                    text: "The first stable product of CO₂ fixation in C4 plants is:",
                    options: [
                        "3-phosphoglycerate",
                        "Oxaloacetic acid", 
                        "Ribulose bisphosphate",
                        "Phosphoenolpyruvate"
                    ],
                    answer: 1,
                    explain: "In C4 pathway, CO₂ combines with phosphoenolpyruvate (PEP) to form oxaloacetic acid as the first stable 4-carbon product."
                },
                {
                    id: "b4",
                    text: "Which plant hormone is responsible for apical dominance?",
                    options: ["Cytokinin", "Gibberellin", "Auxin", "Abscisic acid"],
                    answer: 2,
                    explain: "Auxin produced by apical meristem maintains apical dominance by inhibiting growth of lateral buds."
                },
                {
                    id: "b5",
                    text: "Stomatal opening is facilitated by:",
                    options: [
                        "Accumulation of K⁺ ions in guard cells",
                        "Loss of K⁺ ions from guard cells", 
                        "High CO₂ concentration",
                        "Low light intensity"
                    ],
                    answer: 0,
                    explain: "K⁺ accumulation in guard cells increases osmotic pressure, causing water uptake and stomatal opening."
                },
                {
                    id: "b6",
                    text: "Which of the following is a day-neutral plant?",
                    options: ["Wheat", "Rice", "Tomato", "Chrysanthemum"],
                    answer: 2,
                    explain: "Tomato flowers regardless of photoperiod length, making it a day-neutral plant."
                },
                {
                    id: "b7",
                    text: "The site of light reaction in photosynthesis is:",
                    options: ["Stroma", "Thylakoid membrane", "Inner membrane", "Matrix"],
                    answer: 1,
                    explain: "Light reactions occur on thylakoid membranes where chlorophyll and photosystems are located."
                },
                {
                    id: "b8",
                    text: "Which enzyme is responsible for CO₂ fixation in Calvin cycle?",
                    options: ["PEP carboxylase", "RuBisCO", "Carbonic anhydrase", "Phosphatase"],
                    answer: 1,
                    explain: "RuBisCO (Ribulose bisphosphate carboxylase oxygenase) catalyzes CO₂ fixation in Calvin cycle."
                },
                {
                    id: "b9",
                    text: "Double fertilization in angiosperms results in formation of:",
                    options: [
                        "Two embryos", 
                        "Embryo and endosperm",
                        "Two endosperms", 
                        "Embryo and seed coat"
                    ],
                    answer: 1,
                    explain: "Double fertilization produces diploid embryo (egg + sperm) and triploid endosperm (central cell + sperm)."
                },
                {
                    id: "b10",
                    text: "Which part of the flower develops into fruit?",
                    options: ["Ovule", "Ovary", "Stigma", "Style"],
                    answer: 1,
                    explain: "After fertilization, the ovary develops into fruit while ovules develop into seeds."
                },
                {
                    id: "b11",
                    text: "Gibberellins promote:",
                    options: [
                        "Root growth", 
                        "Stem elongation",
                        "Leaf senescence", 
                        "Stomatal closure"
                    ],
                    answer: 1,
                    explain: "Gibberellins primarily promote stem elongation and internode growth."
                },
                {
                    id: "b12",
                    text: "Which type of placentation is found in tomato?",
                    options: ["Axile", "Parietal", "Free central", "Basal"],
                    answer: 0,
                    explain: "Tomato shows axile placentation where ovules are attached to central axis in a multilocular ovary."
                },
                {
                    id: "b13",
                    text: "The oxygen evolved during photosynthesis comes from:",
                    options: ["CO₂", "H₂O", "Both CO₂ and H₂O", "Chlorophyll"],
                    answer: 1,
                    explain: "Oxygen in photosynthesis comes from water splitting (photolysis) during light reactions."
                },
                {
                    id: "b14",
                    text: "Which of the following is a micronutrient for plants?",
                    options: ["Nitrogen", "Phosphorus", "Iron", "Potassium"],
                    answer: 2,
                    explain: "Iron is a micronutrient required in small quantities, while N, P, K are macronutrients."
                },
                {
                    id: "b15",
                    text: "Mycorrhiza is an example of:",
                    options: ["Parasitism", "Commensalism", "Mutualism", "Predation"],
                    answer: 2,
                    explain: "Mycorrhiza is a mutualistic association between fungi and plant roots benefiting both organisms."
                },
                {
                    id: "b16",
                    text: "Which plant tissue is responsible for secondary growth?",
                    options: ["Apical meristem", "Cambium", "Epidermis", "Cortex"],
                    answer: 1,
                    explain: "Vascular cambium and cork cambium are responsible for secondary growth in woody plants."
                },
                {
                    id: "b17",
                    text: "The term 'totipotency' refers to:",
                    options: [
                        "Ability to photosynthesize",
                        "Capacity to regenerate whole plant", 
                        "Power to resist diseases",
                        "Ability to reproduce"
                    ],
                    answer: 1,
                    explain: "Totipotency is the ability of a single cell to develop into a complete organism."
                },
                {
                    id: "b18",
                    text: "Which of the following shows cymose inflorescence?",
                    options: ["Raceme", "Spike", "Cyme", "Catkin"],
                    answer: 2,
                    explain: "Cyme is a type of cymose inflorescence where main axis terminates in a flower."
                },
                {
                    id: "b19",
                    text: "Photorespiration occurs in:",
                    options: ["C4 plants only", "C3 plants only", "CAM plants only", "All plants"],
                    answer: 1,
                    explain: "Photorespiration is more prominent in C3 plants due to oxygenase activity of RuBisCO."
                },
                {
                    id: "b20",
                    text: "Which pigment is responsible for red color in flowers?",
                    options: ["Chlorophyll", "Anthocyanin", "Carotenoid", "Xanthophyll"],
                    answer: 1,
                    explain: "Anthocyanins are water-soluble pigments responsible for red, purple, and blue colors."
                },
                {
                    id: "b21",
                    text: "The process of formation of gametes is called:",
                    options: ["Fertilization", "Gametogenesis", "Embryogenesis", "Morphogenesis"],
                    answer: 1,
                    explain: "Gametogenesis is the process of gamete formation from reproductive cells."
                },
                {
                    id: "b22",
                    text: "Which hormone inhibits seed germination?",
                    options: ["Gibberellin", "Auxin", "Abscisic acid", "Cytokinin"],
                    answer: 2,
                    explain: "Abscisic acid (ABA) maintains seed dormancy and inhibits germination."
                },
                {
                    id: "b23",
                    text: "Kranz anatomy is characteristic of:",
                    options: ["C3 plants", "C4 plants", "CAM plants", "Hydrophytes"],
                    answer: 1,
                    explain: "Kranz anatomy with bundle sheath cells is a characteristic feature of C4 plants."
                },
                {
                    id: "b24",
                    text: "Which of the following is NOT a function of roots?",
                    options: [
                        "Absorption of water",
                        "Anchorage", 
                        "Photosynthesis",
                        "Storage"
                    ],
                    answer: 2,
                    explain: "Roots typically don't perform photosynthesis as they lack chlorophyll and are underground."
                },
                {
                    id: "b25",
                    text: "The phytochrome system controls:",
                    options: ["Phototropism", "Photoperiodism", "Photosynthesis", "Photorespiration"],
                    answer: 1,
                    explain: "Phytochrome pigment system regulates photoperiodic responses in plants."
                },
                {
                    id: "b26",
                    text: "Which sugar is transported in phloem?",
                    options: ["Glucose", "Fructose", "Sucrose", "Starch"],
                    answer: 2,
                    explain: "Sucrose is the main sugar transported in phloem sap from source to sink."
                },
                {
                    id: "b27",
                    text: "Velamen tissue is found in:",
                    options: ["Hydrophytes", "Xerophytes", "Epiphytes", "Mesophytes"],
                    answer: 2,
                    explain: "Velamen is a spongy tissue in aerial roots of epiphytic orchids for water absorption."
                },
                {
                    id: "b28",
                    text: "Which plant group shows heterospory?",
                    options: ["Bryophytes", "Pteridophytes", "Gymnosperms", "All of these"],
                    answer: 2,
                    explain: "Gymnosperms show heterospory producing different types of spores (microspores and megaspores)."
                },
                {
                    id: "b29",
                    text: "The phenomenon of opening of floral buds into flowers is termed as:",
                    options: ["Anthesis", "Dehiscence", "Abscission", "Senescence"],
                    answer: 0,
                    explain: "Anthesis refers to the opening of flower buds and the period of sexual maturity."
                },
                {
                    id: "b30",
                    text: "Which of the following is a living component of xylem?",
                    options: ["Tracheids", "Vessels", "Xylem parenchyma", "Xylem fibers"],
                    answer: 2,
                    explain: "Xylem parenchyma cells are living and involved in storage and lateral conduction."
                },
                {
                    id: "b31",
                    text: "Pneumatophores are found in:",
                    options: ["Desert plants", "Aquatic plants", "Mangrove plants", "Alpine plants"],
                    answer: 2,
                    explain: "Pneumatophores are specialized aerial roots in mangroves for gaseous exchange."
                },
                {
                    id: "b32",
                    text: "Which of the following exhibits circinate vernation?",
                    options: ["Monocots", "Dicots", "Ferns", "Gymnosperms"],
                    answer: 2,
                    explain: "Ferns show circinate vernation where young leaves are coiled like a watch spring."
                },
                {
                    id: "b33",
                    text: "Guttation occurs through:",
                    options: ["Stomata", "Lenticels", "Hydathodes", "Cuticle"],
                    answer: 2,
                    explain: "Guttation is loss of water in liquid form through hydathodes at leaf margins."
                },
                {
                    id: "b34",
                    text: "Which type of fruit is coconut?",
                    options: ["Drupe", "Berry", "Capsule", "Nut"],
                    answer: 0,
                    explain: "Coconut is a drupe with fibrous mesocarp and hard endocarp surrounding the seed."
                },
                {
                    id: "b35",
                    text: "Synergids are found in:",
                    options: ["Anther", "Embryo sac", "Pollen grain", "Endosperm"],
                    answer: 1,
                    explain: "Synergids are two cells flanking the egg cell in the embryo sac of ovule."
                },
                {
                    id: "b36",
                    text: "Which plant tissue provides mechanical strength?",
                    options: ["Parenchyma", "Collenchyma", "Sclerenchyma", "Aerenchyma"],
                    answer: 2,
                    explain: "Sclerenchyma tissue with lignified cell walls provides maximum mechanical support."
                },
                {
                    id: "b37",
                    text: "CAM plants open their stomata during:",
                    options: ["Day time", "Night time", "Evening", "Morning"],
                    answer: 1,
                    explain: "CAM plants open stomata at night to minimize water loss while fixing CO₂."
                },
                {
                    id: "b38",
                    text: "Which of the following is NOT a characteristic of wind-pollinated flowers?",
                    options: [
                        "Small and inconspicuous",
                        "Absence of nectar", 
                        "Brightly colored petals",
                        "Light pollen grains"
                    ],
                    answer: 2,
                    explain: "Wind-pollinated flowers are typically small, dull-colored, and lack bright petals."
                },
                {
                    id: "b39",
                    text: "Prothallus is the gametophyte of:",
                    options: ["Moss", "Fern", "Pine", "Flowering plant"],
                    answer: 1,
                    explain: "Prothallus is the heart-shaped gametophyte stage in the life cycle of ferns."
                },
                {
                    id: "b40",
                    text: "Which of the following is a gaseous plant hormone?",
                    options: ["Auxin", "Gibberellin", "Ethylene", "Cytokinin"],
                    answer: 2,
                    explain: "Ethylene is the only gaseous plant hormone involved in fruit ripening and senescence."
                },
                {
                    id: "b41",
                    text: "Bulliform cells help in:",
                    options: [
                        "Water storage",
                        "Photosynthesis", 
                        "Rolling and unrolling of leaves",
                        "Gaseous exchange"
                    ],
                    answer: 2,
                    explain: "Bulliform cells in grass leaves help in rolling/unrolling for water conservation."
                },
                {
                    id: "b42",
                    text: "Which of the following shows apomixis?",
                    options: ["Citrus", "Rose", "Mango", "All of these"],
                    answer: 3,
                    explain: "Apomixis (asexual reproduction without fertilization) is found in citrus, rose, and mango."
                },
                {
                    id: "b43",
                    text: "Root pressure is maximum during:",
                    options: ["Day time", "Night time", "Noon", "Evening"],
                    answer: 1,
                    explain: "Root pressure is maximum at night when transpiration is minimal and water absorption continues."
                },
                {
                    id: "b44",
                    text: "Which part of seed develops from nucellus?",
                    options: ["Testa", "Tegmen", "Perisperm", "Endosperm"],
                    answer: 2,
                    explain: "Perisperm is nutritive tissue that develops from nucellus in some seeds."
                },
                {
                    id: "b45",
                    text: "Archegonium is the female reproductive structure in:",
                    options: ["Angiosperms", "Gymnosperms", "Bryophytes", "Both B and C"],
                    answer: 3,
                    explain: "Archegonium is found in both gymnosperms and bryophytes as female reproductive organ."
                },

                // ZOOLOGY Questions (45)
                {
                    id: "z1",
                    text: "Which of the following is NOT a function of human kidney?",
                    options: [
                        "Regulation of blood pressure",
                        "Production of erythropoietin", 
                        "Gluconeogenesis",
                        "Production of insulin"
                    ],
                    answer: 3,
                    explain: "Insulin is produced by beta cells of pancreas, not by kidneys. Kidneys produce erythropoietin and renin."
                },
                {
                    id: "z2",
                    text: "The functional unit of kidney is:",
                    options: ["Glomerulus", "Nephron", "Loop of Henle", "Collecting duct"],
                    answer: 1,
                    explain: "Nephron is the structural and functional unit of kidney consisting of glomerulus and tubular parts."
                },
                {
                    id: "z3",
                    text: "Which hormone regulates calcium levels in blood?",
                    options: ["Insulin", "Thyroxine", "Parathyroid hormone", "Growth hormone"],
                    answer: 2,
                    explain: "Parathyroid hormone (PTH) increases blood calcium by promoting bone resorption and calcium absorption."
                },
                {
                    id: "z4",
                    text: "The pacemaker of heart is:",
                    options: ["AV node", "SA node", "Bundle of His", "Purkinje fibers"],
                    answer: 1,
                    explain: "Sino-atrial (SA) node initiates heartbeat and sets the rhythm, hence called natural pacemaker."
                },
                {
                    id: "z5",
                    text: "Which blood vessel carries oxygenated blood to the heart muscle?",
                    options: ["Pulmonary artery", "Coronary artery", "Vena cava", "Pulmonary vein"],
                    answer: 1,
                    explain: "Coronary arteries supply oxygenated blood to the heart muscle (myocardium)."
                },
                {
                    id: "z6",
                    text: "Insulin deficiency causes:",
                    options: ["Diabetes mellitus", "Diabetes insipidus", "Acromegaly", "Goiter"],
                    answer: 0,
                    explain: "Insulin deficiency or resistance leads to diabetes mellitus characterized by high blood glucose."
                },
                {
                    id: "z7",
                    text: "Which part of brain controls breathing?",
                    options: ["Cerebrum", "Cerebellum", "Medulla oblongata", "Hypothalamus"],
                    answer: 2,
                    explain: "Medulla oblongata contains respiratory center that controls automatic breathing."
                },
                {
                    id: "z8",
                    text: "The exchange of gases in lungs occurs at:",
                    options: ["Bronchi", "Bronchioles", "Alveoli", "Trachea"],
                    answer: 2,
                    explain: "Alveoli are tiny air sacs where gas exchange occurs between air and blood capillaries."
                },
                {
                    id: "z9",
                    text: "Which enzyme initiates protein digestion in stomach?",
                    options: ["Trypsin", "Pepsin", "Chymotrypsin", "Erepsin"],
                    answer: 1,
                    explain: "Pepsin secreted by stomach initiates protein digestion in acidic environment."
                },
                {
                    id: "z10",
                    text: "The longest bone in human body is:",
                    options: ["Tibia", "Fibula", "Femur", "Humerus"],
                    answer: 2,
                    explain: "Femur (thigh bone) is the longest and strongest bone in human body."
                },
                {
                    id: "z11",
                    text: "Which vitamin is essential for blood clotting?",
                    options: ["Vitamin A", "Vitamin C", "Vitamin D", "Vitamin K"],
                    answer: 3,
                    explain: "Vitamin K is essential for synthesis of clotting factors and normal blood coagulation."
                },
                {
                    id: "z12",
                    text: "The cells responsible for producing antibodies are:",
                    options: ["T cells", "B cells", "NK cells", "Macrophages"],
                    answer: 1,
                    explain: "B cells differentiate into plasma cells that produce and secrete specific antibodies."
                },
                {
                    id: "z13",
                    text: "Which gland is called 'master gland' of endocrine system?",
                    options: ["Thyroid", "Adrenal", "Pituitary", "Pancreas"],
                    answer: 2,
                    explain: "Pituitary gland controls other endocrine glands through its hormones, hence called master gland."
                },
                {
                    id: "z14",
                    text: "The pigment responsible for oxygen transport in blood is:",
                    options: ["Chlorophyll", "Hemoglobin", "Myoglobin", "Cytochrome"],
                    answer: 1,
                    explain: "Hemoglobin in RBCs binds and transports oxygen from lungs to tissues."
                },
                {
                    id: "z15",
                    text: "Which part of eye controls the amount of light entering?",
                    options: ["Cornea", "Lens", "Iris", "Retina"],
                    answer: 2,
                    explain: "Iris contains muscles that control pupil size, regulating light entry into eye."
                },
                {
                    id: "z16",
                    text: "The site of fertilization in human female is:",
                    options: ["Ovary", "Uterus", "Fallopian tube", "Vagina"],
                    answer: 2,
                    explain: "Fertilization typically occurs in the ampulla region of fallopian tube."
                },
                {
                    id: "z17",
                    text: "Which hormone is responsible for milk ejection?",
                    options: ["Prolactin", "Oxytocin", "Progesterone", "Estrogen"],
                    answer: 1,
                    explain: "Oxytocin causes contraction of mammary glands leading to milk ejection (let-down reflex)."
                },
                {
                    id: "z18",
                    text: "The process of sperm formation is called:",
                    options: ["Oogenesis", "Spermatogenesis", "Gametogenesis", "Fertilization"],
                    answer: 1,
                    explain: "Spermatogenesis is the process of sperm formation in seminiferous tubules of testes."
                },
                {
                    id: "z19",
                    text: "Which of the following is NOT a sexually transmitted infection?",
                    options: ["Gonorrhea", "Syphilis", "Tuberculosis", "AIDS"],
                    answer: 2,
                    explain: "Tuberculosis is transmitted through respiratory droplets, not sexual contact."
                },
                {
                    id: "z20",
                    text: "The gestation period in humans is approximately:",
                    options: ["38 weeks", "40 weeks", "42 weeks", "36 weeks"],
                    answer: 1,
                    explain: "Human gestation period is approximately 40 weeks or 280 days from last menstrual period."
                },
                {
                    id: "z21",
                    text: "Which layer of heart wall is responsible for contraction?",
                    options: ["Pericardium", "Epicardium", "Myocardium", "Endocardium"],
                    answer: 2,
                    explain: "Myocardium is the muscular layer of heart wall responsible for cardiac contractions."
                },
                {
                    id: "z22",
                    text: "The normal pH range of human blood is:",
                    options: ["6.8-7.0", "7.35-7.45", "7.8-8.0", "8.1-8.3"],
                    answer: 1,
                    explain: "Normal blood pH is tightly regulated between 7.35-7.45 (slightly alkaline)."
                },
                {
                    id: "z23",
                    text: "Which structure prevents food from entering trachea during swallowing?",
                    options: ["Epiglottis", "Glottis", "Larynx", "Pharynx"],
                    answer: 0,
                    explain: "Epiglottis acts as a lid that covers the glottis during swallowing to prevent aspiration."
                },
                {
                    id: "z24",
                    text: "The neurotransmitter at neuromuscular junction is:",
                    options: ["Dopamine", "Serotonin", "Acetylcholine", "GABA"],
                    answer: 2,
                    explain: "Acetylcholine is released at neuromuscular junction to stimulate muscle contraction."
                },
                {
                    id: "z25",
                    text: "Which cell organelle is called 'powerhouse of cell'?",
                    options: ["Nucleus", "Mitochondria", "Ribosome", "Golgi body"],
                    answer: 1,
                    explain: "Mitochondria produce ATP through cellular respiration, hence called powerhouse of cell."
                },
                {
                    id: "z26",
                    text: "The total number of bones in adult human body is:",
                    options: ["204", "206", "208", "210"],
                    answer: 1,
                    explain: "Adult human skeleton has 206 bones (270 at birth, many fuse during development)."
                },
                {
                    id: "z27",
                    text: "Which hormone is deficient in diabetes insipidus?",
                    options: ["Insulin", "ADH", "Aldosterone", "Cortisol"],
                    answer: 1,
                    explain: "Antidiuretic hormone (ADH/vasopressin) deficiency causes diabetes insipidus with excessive urination."
                },
                {
                    id: "z28",
                    text: "The blind spot in human eye is due to:",
                    options: [
                        "Absence of rods and cones",
                        "Presence of optic nerve", 
                        "Absence of blood vessels",
                        "Presence of lens"
                    ],
                    answer: 0,
                    explain: "Blind spot occurs where optic nerve exits retina - no photoreceptors (rods/cones) present there."
                },
                {
                    id: "z29",
                    text: "Which blood group is called universal donor?",
                    options: ["A", "B", "AB", "O"],
                    answer: 3,
                    explain: "O blood group lacks A and B antigens, so can donate to all groups (universal donor)."
                },
                {
                    id: "z30",
                    text: "The hormone that stimulates uterine contractions during labor is:",
                    options: ["Relaxin", "Progesterone", "Oxytocin", "Estrogen"],
                    answer: 2,
                    explain: "Oxytocin stimulates strong uterine contractions during labor and delivery."
                },
                {
                    id: "z31",
                    text: "Which part of nephron is impermeable to water?",
                    options: [
                        "Proximal tubule", 
                        "Descending limb of loop of Henle",
                        "Ascending limb of loop of Henle", 
                        "Distal tubule"
                    ],
                    answer: 2,
                    explain: "Ascending limb of loop of Henle is impermeable to water but actively transports Na⁺ and Cl⁻."
                },
                {
                    id: "z32",
                    text: "The enzyme that converts fibrinogen to fibrin is:",
                    options: ["Thrombin", "Prothrombin", "Fibrinase", "Heparin"],
                    answer: 0,
                    explain: "Thrombin converts soluble fibrinogen to insoluble fibrin threads for blood clot formation."
                },
                {
                    id: "z33",
                    text: "Which cranial nerve controls vision?",
                    options: ["Optic", "Oculomotor", "Trigeminal", "Facial"],
                    answer: 0,
                    explain: "Optic nerve (cranial nerve II) carries visual information from retina to brain."
                },
                {
                    id: "z34",
                    text: "The immediate source of energy for muscle contraction is:",
                    options: ["Glucose", "ATP", "Creatine phosphate", "Lactic acid"],
                    answer: 1,
                    explain: "ATP provides immediate energy for muscle contraction by powering myosin head movement."
                },
                {
                    id: "z35",
                    text: "Which structure in ear maintains balance and equilibrium?",
                    options: ["Cochlea", "Semicircular canals", "Tympanum", "Auditory canal"],
                    answer: 1,
                    explain: "Semicircular canals in inner ear detect rotational movements and maintain balance."
                },
                {
                    id: "z36",
                    text: "The largest lymphatic organ in human body is:",
                    options: ["Thymus", "Spleen", "Tonsils", "Lymph nodes"],
                    answer: 1,
                    explain: "Spleen is the largest lymphatic organ involved in blood filtration and immune responses."
                },
                {
                    id: "z37",
                    text: "Which vitamin deficiency causes night blindness?",
                    options: ["Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D"],
                    answer: 0,
                    explain: "Vitamin A deficiency affects rhodopsin synthesis, causing night blindness (nyctalopia)."
                },
                {
                    id: "z38",
                    text: "The vocal cords are located in:",
                    options: ["Pharynx", "Larynx", "Trachea", "Bronchi"],
                    answer: 1,
                    explain: "Vocal cords (vocal folds) are located in larynx and vibrate to produce sound."
                },
                {
                    id: "z39",
                    text: "Which hormone regulates blood sugar level?",
                    options: ["Thyroxine", "Insulin", "Adrenaline", "Growth hormone"],
                    answer: 1,
                    explain: "Insulin lowers blood glucose by promoting cellular glucose uptake and glycogen synthesis."
                },
                {
                    id: "z40",
                    text: "The gap between two neurons is called:",
                    options: ["Node of Ranvier", "Synapse", "Dendrite", "Axon hillock"],
                    answer: 1,
                    explain: "Synapse is the junction between two neurons where neurotransmission occurs."
                },
                {
                    id: "z41",
                    text: "Which part of sperm contains enzymes for egg penetration?",
                    options: ["Head", "Acrosome", "Mitochondria", "Flagellum"],
                    answer: 1,
                    explain: "Acrosome contains enzymes like hyaluronidase that help sperm penetrate egg coats."
                },
                {
                    id: "z42",
                    text: "The yellowish color of urine is due to:",
                    options: ["Urea", "Urochrome", "Creatinine", "Ammonia"],
                    answer: 1,
                    explain: "Urochrome, a breakdown product of hemoglobin, gives urine its characteristic yellow color."
                },
                {
                    id: "z43",
                    text: "Which type of muscle tissue is involuntary and striated?",
                    options: ["Skeletal muscle", "Cardiac muscle", "Smooth muscle", "None of these"],
                    answer: 1,
                    explain: "Cardiac muscle is both involuntary (autonomic control) and striated (cross-striped appearance)."
                },
                {
                    id: "z44",
                    text: "The condition where lens becomes opaque is called:",
                    options: ["Glaucoma", "Cataract", "Myopia", "Hyperopia"],
                    answer: 1,
                    explain: "Cataract is clouding of the eye lens causing vision impairment or blindness."
                },
                {
                    id: "z45",
                    text: "Which organ produces bile?",
                    options: ["Pancreas", "Liver", "Gallbladder", "Stomach"],
                    answer: 1,
                    explain: "Liver produces bile which is stored in gallbladder and released into duodenum for fat digestion."
                }
            ]
        }
    ]
};
