// mock-test-10.js - NEET Mock Test 10 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_10 = {
    id: "neet-010",
    title: "Full Syllabus Mock 10", 
    level: "medium",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A body of mass 5 kg is acted upon by forces 3î + 4ĵ N and -2î + 7ĵ N. The magnitude of acceleration is:",
                    options: ["1 m/s²", "√26 m/s²", "2.2 m/s²", "3 m/s²"],
                    answer: 2,
                    explain: "Net force = (3-2)î + (4+7)ĵ = î + 11ĵ N. |F| = √(1² + 11²) = √122 N. Acceleration = √122/5 = 2.2 m/s²"
                },
                {
                    id: "p2",
                    text: "A wheel of radius R rolls without slipping. The velocity of point at top of wheel is:",
                    options: ["v", "2v", "v/2", "√2v"],
                    answer: 1,
                    explain: "For rolling motion, velocity of center = v. Velocity of top point = v(translation) + v(rotation) = v + v = 2v"
                },
                {
                    id: "p3",
                    text: "Two capacitors 4μF and 6μF are connected in series across 100V. Voltage across 4μF capacitor is:",
                    options: ["40V", "60V", "50V", "25V"],
                    answer: 1,
                    explain: "In series, voltage divides inversely as capacitance. V₄ = V × C₆/(C₄+C₆) = 100 × 6/(4+6) = 60V"
                },
                {
                    id: "p4",
                    text: "In single slit diffraction, if slit width is doubled, the angular width of central maximum becomes:",
                    options: ["Double", "Half", "Four times", "One-fourth"],
                    answer: 1,
                    explain: "Angular width of central maximum = 2λ/a. When slit width 'a' doubles, angular width becomes half"
                },
                {
                    id: "p5",
                    text: "Energy of photon with wavelength 6000 Å is approximately:",
                    options: ["2.07 eV", "3.1 eV", "1.24 eV", "4.14 eV"],
                    answer: 0,
                    explain: "E = hc/λ = (6.626×10⁻³⁴ × 3×10⁸)/(6000×10⁻¹⁰ × 1.6×10⁻¹⁹) = 2.07 eV"
                },
                {
                    id: "p6",
                    text: "A conducting rod moves in magnetic field. If both velocity and field are doubled, induced EMF becomes:",
                    options: ["Same", "Double", "Four times", "Half"],
                    answer: 2,
                    explain: "EMF = BLv. When both B and v are doubled: EMF' = (2B)L(2v) = 4BLv = 4 times original EMF"
                },
                {
                    id: "p7",
                    text: "Two discs of same mass and radius, one solid and other hollow, roll down incline. Which reaches bottom first?",
                    options: ["Solid disc", "Hollow disc", "Both together", "Cannot determine"],
                    answer: 0,
                    explain: "For rolling: a = g sinθ/(1 + I/mR²). Solid disc has smaller I/mR² ratio, hence larger acceleration"
                },
                {
                    id: "p8",
                    text: "In RLC circuit, power factor is 0.6. If resistance is 3Ω, impedance is:",
                    options: ["5Ω", "4Ω", "3.6Ω", "1.8Ω"],
                    answer: 0,
                    explain: "Power factor = R/Z = 0.6. Given R = 3Ω, so Z = R/0.6 = 3/0.6 = 5Ω"
                },
                {
                    id: "p9",
                    text: "In reversible adiabatic expansion, internal energy of ideal gas:",
                    options: ["Increases", "Decreases", "Remains constant", "May increase or decrease"],
                    answer: 1,
                    explain: "In adiabatic expansion, gas does work against external pressure. Since Q = 0, ΔU = -W < 0 (decreases)"
                },
                {
                    id: "p10",
                    text: "Current loop in uniform magnetic field experiences:",
                    options: ["Force only", "Torque only", "Both force and torque", "Neither force nor torque"],
                    answer: 1,
                    explain: "In uniform magnetic field, net force on current loop is zero but there can be torque if plane of loop is not perpendicular to field"
                },
                {
                    id: "p11",
                    text: "Two SHM of same amplitude but frequencies in ratio 1:2 are superposed. Resultant motion is:",
                    options: ["SHM", "Periodic but not SHM", "Non-periodic", "Uniform"],
                    answer: 1,
                    explain: "Superposition of two SHM with frequency ratio 1:2 gives periodic motion but not simple harmonic"
                },
                {
                    id: "p12",
                    text: "Compton shift depends on:",
                    options: ["Incident photon energy", "Scattering angle", "Target material", "Both energy and angle"],
                    answer: 1,
                    explain: "Compton shift Δλ = (h/mₑc)(1 - cosθ) depends only on scattering angle θ, not on incident energy or target"
                },
                {
                    id: "p13",
                    text: "A battery of EMF 12V and internal resistance 2Ω is connected to external resistance 4Ω. Terminal voltage is:",
                    options: ["6V", "8V", "10V", "12V"],
                    answer: 1,
                    explain: "Current I = E/(r+R) = 12/(2+4) = 2A. Terminal voltage V = E - Ir = 12 - 2×2 = 8V"
                },
                {
                    id: "p14",
                    text: "Range of projectile is maximum when angle of projection is:",
                    options: ["30°", "45°", "60°", "90°"],
                    answer: 1,
                    explain: "Range R = u²sin2θ/g is maximum when sin2θ = 1, i.e., 2θ = 90°, so θ = 45°"
                },
                {
                    id: "p15",
                    text: "Self-inductance of solenoid is proportional to:",
                    options: ["n", "n²", "1/n", "√n"],
                    answer: 1,
                    explain: "L = μ₀n²Al where n is turns per unit length. Self-inductance is proportional to square of turn density"
                },
                {
                    id: "p16",
                    text: "A ball dropped from 20m height loses 25% energy on first bounce. Height after bounce is:",
                    options: ["15m", "10m", "5m", "7.5m"],
                    answer: 0,
                    explain: "After losing 25% energy, remaining energy = 75% = 0.75. Height ∝ Energy, so h = 0.75 × 20 = 15m"
                },
                {
                    id: "p17",
                    text: "Quality factor Q of LCR series circuit is:",
                    options: ["ωL/R", "R/ωL", "1/ωRC", "ωRC"],
                    answer: 0,
                    explain: "Quality factor Q = ωL/R = 1/(ωRC) at resonance where ωL = 1/(ωC)"
                },
                {
                    id: "p18",
                    text: "Wave function Ψ represents:",
                    options: ["Probability", "Probability density", "Probability amplitude", "Energy"],
                    answer: 2,
                    explain: "Wave function Ψ is probability amplitude. |Ψ|² gives probability density of finding particle"
                },
                {
                    id: "p19",
                    text: "Time period of compound pendulum is minimum when:",
                    options: ["k = l", "k² = l²", "k² + l² is minimum", "l = 0"],
                    answer: 0,
                    explain: "T = 2π√[(k² + l²)/(gl)]. For minimum T, differentiate w.r.t. l and set = 0, gives k = l"
                },
                {
                    id: "p20",
                    text: "Electric field just outside conductor surface having surface charge density σ is:",
                    options: ["σ/ε₀", "σ/2ε₀", "2σ/ε₀", "0"],
                    answer: 0,
                    explain: "Using Gauss's law for conductor surface: E = σ/ε₀ just outside the surface"
                },
                {
                    id: "p21",
                    text: "Minimum deviation by prism occurs when:",
                    options: ["i₁ = i₂", "r₁ = r₂", "Both i₁ = i₂ and r₁ = r₂", "i₁ = 0"],
                    answer: 2,
                    explain: "At minimum deviation, ray passes symmetrically through prism: i₁ = i₂ and r₁ = r₂"
                },
                {
                    id: "p22",
                    text: "In wave interference, if path difference is (n + 1/2)λ, interference is:",
                    options: ["Constructive", "Destructive", "Partially constructive", "No interference"],
                    answer: 1,
                    explain: "Path difference = (n + 1/2)λ corresponds to odd multiple of λ/2, giving destructive interference"
                },
                {
                    id: "p23",
                    text: "In cyclic process, change in internal energy is:",
                    options: ["Positive", "Negative", "Zero", "Cannot determine"],
                    answer: 2,
                    explain: "In cyclic process, system returns to initial state, so ΔU = 0"
                },
                {
                    id: "p24",
                    text: "Activity of radioactive sample is 8 Curie. After 3 half-lives, activity becomes:",
                    options: ["1 Curie", "2 Curie", "4 Curie", "0.5 Curie"],
                    answer: 0,
                    explain: "After n half-lives: A = A₀/2ⁿ. After 3 half-lives: A = 8/2³ = 8/8 = 1 Curie"
                },
                {
                    id: "p25",
                    text: "Magnetic energy stored in inductor is:",
                    options: ["½LI²", "LI²", "2LI²", "LI"],
                    answer: 0,
                    explain: "Magnetic energy stored in inductor: U = ½LI² where L is inductance and I is current"
                },
                {
                    id: "p26",
                    text: "Angular momentum of electron in hydrogen atom is quantized as:",
                    options: ["nh/2π", "nh", "h/2πn", "nh/4π"],
                    answer: 0,
                    explain: "Bohr's quantization condition: mvr = nh/2π where n is principal quantum number"
                },
                {
                    id: "p27",
                    text: "Brewster's angle for glass (μ = 1.5) is approximately:",
                    options: ["33°", "48°", "56°", "63°"],
                    answer: 2,
                    explain: "Brewster's angle: tanθ = μ = 1.5, so θ = arctan(1.5) ≈ 56°"
                },
                {
                    id: "p28",
                    text: "Standing wave has nodes separated by distance:",
                    options: ["λ/4", "λ/2", "λ", "2λ"],
                    answer: 1,
                    explain: "In standing wave, adjacent nodes are separated by λ/2 where λ is wavelength"
                },
                {
                    id: "p29",
                    text: "Otto cycle efficiency depends on:",
                    options: ["Compression ratio only", "Heat capacities only", "Both compression ratio and heat capacities", "Working substance"],
                    answer: 2,
                    explain: "Otto cycle efficiency η = 1 - 1/r^(γ-1) depends on compression ratio r and γ = Cp/Cv"
                },
                {
                    id: "p30",
                    text: "Lorentz force on charged particle is maximum when velocity is:",
                    options: ["Parallel to B", "Antiparallel to B", "Perpendicular to B", "At 45° to B"],
                    answer: 2,
                    explain: "Lorentz force F = qv×B is maximum when sinθ = 1, i.e., when v ⊥ B (θ = 90°)"
                },
                {
                    id: "p31",
                    text: "Doppler shift in light is due to:",
                    options: ["Relative motion only", "Medium properties", "Wave nature", "Relativistic effects"],
                    answer: 3,
                    explain: "Doppler shift in light requires relativistic treatment as light speed is independent of source motion"
                },
                {
                    id: "p32",
                    text: "Radius of gyration of solid cylinder about its axis is:",
                    options: ["R/2", "R/√2", "R/√3", "R"],
                    answer: 1,
                    explain: "For solid cylinder: I = ½MR², so k = √(I/M) = √(½R²) = R/√2"
                },
                {
                    id: "p33",
                    text: "RMS value of AC voltage V = V₀sinωt is:",
                    options: ["V₀", "V₀/2", "V₀/√2", "V₀√2"],
                    answer: 2,
                    explain: "For sinusoidal AC: Vrms = V₀/√2 where V₀ is peak voltage"
                },
                {
                    id: "p34",
                    text: "Power of convex lens is +5D. When immersed in water (μ = 4/3), power becomes:",
                    options: ["+3.75D", "+1.25D", "+2.5D", "+6.67D"],
                    answer: 1,
                    explain: "P₍water₎/P₍air₎ = (μ_lens/μ_water - 1)/(μ_lens/μ_air - 1). For glass in water: P = P₀/4, so P = 5/4 = 1.25D"
                },
                {
                    id: "p35",
                    text: "Nuclear binding energy curve has maximum at:",
                    options: ["A = 4", "A = 56", "A = 238", "A = 12"],
                    answer: 1,
                    explain: "Binding energy per nucleon is maximum around mass number A = 56 (Iron-56), making it most stable nucleus"
                },
                {
                    id: "p36",
                    text: "For damped oscillations, amplitude decreases:",
                    options: ["Linearly", "Exponentially", "Quadratically", "As inverse square"],
                    answer: 1,
                    explain: "In damped harmonic motion, amplitude decreases exponentially: A(t) = A₀e^(-γt/2m)"
                },
                {
                    id: "p37",
                    text: "Ampere's circuital law is:",
                    options: ["∮B⋅dl = μ₀I", "∮E⋅dl = 0", "∮D⋅dA = Q", "∮B⋅dA = 0"],
                    answer: 0,
                    explain: "Ampere's law: ∮B⋅dl = μ₀Ienc where Ienc is current enclosed by the path"
                },
                {
                    id: "p38",
                    text: "Turns ratio of step-down transformer is 10:1. If primary voltage is 220V, secondary voltage is:",
                    options: ["22V", "44V", "110V", "2200V"],
                    answer: 0,
                    explain: "For transformer: Vs/Vp = Ns/Np = 1/10. So Vs = 220/10 = 22V"
                },
                {
                    id: "p39",
                    text: "Growth of current in RL circuit follows:",
                    options: ["Linear law", "Exponential law", "Square law", "Inverse law"],
                    answer: 1,
                    explain: "Current growth in RL circuit: I(t) = I₀(1 - e^(-Rt/L)) follows exponential law"
                },
                {
                    id: "p40",
                    text: "Equipotential surfaces and electric field lines are:",
                    options: ["Parallel", "Antiparallel", "Perpendicular", "At 45°"],
                    answer: 2,
                    explain: "Electric field lines are always perpendicular to equipotential surfaces"
                },
                {
                    id: "p41",
                    text: "Bulk modulus of incompressible liquid is:",
                    options: ["Zero", "Infinite", "Unity", "Negative"],
                    answer: 1,
                    explain: "For incompressible liquid, volume change is zero, so bulk modulus K = -V(dP/dV) = ∞"
                },
                {
                    id: "p42",
                    text: "In AC circuit with pure inductor, current lags voltage by:",
                    options: ["0°", "45°", "90°", "180°"],
                    answer: 2,
                    explain: "In pure inductive circuit, current lags voltage by 90° (π/2 radians)"
                },
                {
                    id: "p43",
                    text: "Coherent sources have:",
                    options: ["Same amplitude", "Same frequency", "Constant phase difference", "Same intensity"],
                    answer: 2,
                    explain: "Coherent sources maintain constant phase relationship for observable interference pattern"
                },
                {
                    id: "p44",
                    text: "Geostationary satellite orbits at height approximately:",
                    options: ["200 km", "2000 km", "36000 km", "200000 km"],
                    answer: 2,
                    explain: "Geostationary satellite orbits at height ~36,000 km above Earth's surface with 24-hour period"
                },
                {
                    id: "p45",
                    text: "Common emitter amplifier has voltage gain 100. If input impedance is 1kΩ and load resistance is 10kΩ, current gain is:",
                    options: ["10", "100", "1000", "10000"],
                    answer: 0,
                    explain: "Voltage gain = Current gain × (Load resistance/Input impedance). 100 = β × (10k/1k), so β = 10"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which has the highest melting point?",
                    options: ["NaF", "NaCl", "NaBr", "NaI"],
                    answer: 0,
                    explain: "NaF has highest lattice energy due to small F⁻ ion size, resulting in highest melting point among sodium halides"
                },
                {
                    id: "c2",
                    text: "The molecular geometry of XeOF₄ is:",
                    options: ["Tetrahedral", "Square pyramidal", "Octahedral", "Square planar"],
                    answer: 1,
                    explain: "XeOF₄ has 6 electron pairs (5 bonding + 1 lone pair) around Xe, giving square pyramidal geometry"
                },
                {
                    id: "c3",
                    text: "Which reaction follows SN1 mechanism?",
                    options: ["CH₃Br + OH⁻", "(CH₃)₃CBr + H₂O", "CH₃CH₂Br + CN⁻", "CH₃I + NH₃"],
                    answer: 1,
                    explain: "Tertiary halide (CH₃)₃CBr with polar protic solvent H₂O favors SN1 mechanism via carbocation intermediate"
                },
                {
                    id: "c4",
                    text: "The oxidation state of chromium in Cr₂O₇²⁻ is:",
                    options: ["+3", "+4", "+5", "+6"],
                    answer: 3,
                    explain: "In dichromate ion Cr₂O₇²⁻: 2(Cr) + 7(-2) = -2, solving gives Cr = +6"
                },
                {
                    id: "c5",
                    text: "How many stereoisomers does CHBrClF have?",
                    options: ["0", "1", "2", "4"],
                    answer: 2,
                    explain: "CHBrClF has one chiral carbon with 4 different substituents, giving 2 enantiomers (stereoisomers)"
                },
                {
                    id: "c6",
                    text: "The strongest Brønsted base is:",
                    options: ["OH⁻", "NH₂⁻", "CH₃⁻", "H⁻"],
                    answer: 3,
                    explain: "H⁻ (hydride ion) is strongest Brønsted base as its conjugate acid H₂ is extremely weak"
                },
                {
                    id: "c7",
                    text: "Electronic configuration of Mn²⁺ is:",
                    options: ["[Ar] 3d⁵", "[Ar] 3d³", "[Ar] 3d⁴ 4s¹", "[Ar] 3d⁶"],
                    answer: 0,
                    explain: "Mn (25): [Ar] 3d⁵ 4s². Mn²⁺ loses 4s² electrons: [Ar] 3d⁵"
                },
                {
                    id: "c8",
                    text: "Which molecule has see-saw geometry?",
                    options: ["SF₄", "CF₄", "XeF₄", "BF₄⁻"],
                    answer: 0,
                    explain: "SF₄ has 5 electron pairs (4 bonding + 1 lone pair) giving see-saw or distorted tetrahedral geometry"
                },
                {
                    id: "c9",
                    text: "Which cannot exist according to molecular orbital theory?",
                    options: ["He₂⁺", "He₂", "H₂⁺", "Li₂"],
                    answer: 1,
                    explain: "He₂ has bond order = (2-2)/2 = 0, so it cannot exist as stable molecule"
                },
                {
                    id: "c10",
                    text: "Entropy change is negative for:",
                    options: ["Melting", "Vaporization", "Crystallization", "Sublimation"],
                    answer: 2,
                    explain: "Crystallization involves ordering of particles, decreasing randomness, hence negative entropy change"
                },
                {
                    id: "c11",
                    text: "Best reducing agent among alkali metals in aqueous solution is:",
                    options: ["Li", "Na", "K", "Cs"],
                    answer: 0,
                    explain: "Despite standard potential, Li is best reducing agent in aqueous solution due to high hydration enthalpy"
                },
                {
                    id: "c12",
                    text: "Number of unpaired electrons in [Fe(H₂O)₆]³⁺ is:",
                    options: ["1", "3", "5", "6"],
                    answer: 2,
                    explain: "Fe³⁺ (d⁵) with weak field H₂O ligands forms high-spin complex with 5 unpaired electrons"
                },
                {
                    id: "c13",
                    text: "Which undergoes Friedel-Crafts acylation most easily?",
                    options: ["Benzene", "Toluene", "Anisole", "Nitrobenzene"],
                    answer: 2,
                    explain: "Anisole has strongly activating -OCH₃ group that enhances electrophilic aromatic substitution"
                },
                {
                    id: "c14",
                    text: "The correct order of atomic radii is:",
                    options: ["F < Cl < Br < I", "I < Br < Cl < F", "Cl < F < Br < I", "F < I < Br < Cl"],
                    answer: 0,
                    explain: "Atomic radius increases down the group: F < Cl < Br < I due to additional electron shells"
                },
                {
                    id: "c15",
                    text: "Henderson-Hasselbalch equation is:",
                    options: ["pH = pKa + log[A⁻]/[HA]", "pH = pKa - log[A⁻]/[HA]", "pH = pKw - pOH", "pH = -log[H⁺]"],
                    answer: 0,
                    explain: "Henderson-Hasselbalch equation: pH = pKa + log([A⁻]/[HA]) for buffer calculations"
                },
                {
                    id: "c16",
                    text: "Spin-only magnetic moment of [CoF₆]³⁻ is:",
                    options: ["0 BM", "2.83 BM", "4.90 BM", "5.92 BM"],
                    answer: 2,
                    explain: "Co³⁺ (d⁶) with weak field F⁻ ligands forms high-spin complex with 4 unpaired electrons. μ = √[4×6] = 4.90 BM"
                },
                {
                    id: "c17",
                    text: "Which has maximum covalent character in bonding?",
                    options: ["LiF", "LiCl", "LiBr", "LiI"],
                    answer: 3,
                    explain: "By Fajan's rules, LiI has maximum covalent character due to large, polarizable I⁻ ion"
                },
                {
                    id: "c18",
                    text: "Arrhenius equation is:",
                    options: ["k = Ae^(-Ea/RT)", "k = Ae^(Ea/RT)", "k = A + Ea/RT", "k = A - Ea/RT"],
                    answer: 0,
                    explain: "Arrhenius equation: k = Ae^(-Ea/RT) relates rate constant to activation energy and temperature"
                },
                {
                    id: "c19",
                    text: "Which shows maximum hydrogen bonding?",
                    options: ["CH₃OH", "C₂H₅OH", "C₃H₇OH", "C₄H₉OH"],
                    answer: 0,
                    explain: "CH₃OH has highest OH/molecular weight ratio, giving maximum hydrogen bonding per unit mass"
                },
                {
                    id: "c20",
                    text: "Which is anti-aromatic?",
                    options: ["Benzene", "Naphthalene", "Cyclobutadiene", "Pyridine"],
                    answer: 2,
                    explain: "Cyclobutadiene has 4π electrons (4n, n=1), making it anti-aromatic according to Hückel's rule"
                },
                {
                    id: "c21",
                    text: "Which has zero dipole moment?",
                    options: ["H₂S", "SO₂", "CO₂", "H₂O"],
                    answer: 2,
                    explain: "CO₂ has linear geometry with equal and opposite dipoles that cancel out, giving zero net dipole moment"
                },
                {
                    id: "c22",
                    text: "Denticity of EDTA is:",
                    options: ["4", "5", "6", "8"],
                    answer: 2,
                    explain: "EDTA has 6 donor sites (2 nitrogens + 4 carboxyl oxygens), making it hexadentate ligand"
                },
                {
                    id: "c23",
                    text: "Strongest acid among these is:",
                    options: ["CH₃COOH", "ClCH₂COOH", "Cl₂CHCOOH", "Cl₃CCOOH"],
                    answer: 3,
                    explain: "Trichloroacetic acid (Cl₃CCOOH) is strongest due to maximum -I effect of three chlorine atoms"
                },
                {
                    id: "c24",
                    text: "Bond order of NO is:",
                    options: ["2", "2.5", "3", "1.5"],
                    answer: 1,
                    explain: "NO has 15 electrons. Bond order = (bonding - antibonding)/2 = (10-5)/2 = 2.5"
                },
                {
                    id: "c25",
                    text: "Which undergoes nucleophilic addition reaction?",
                    options: ["Ethene", "Benzene", "Acetone", "Methane"],
                    answer: 2,
                    explain: "Acetone has electrophilic carbonyl carbon that readily undergoes nucleophilic addition reactions"
                },
                {
                    id: "c26",
                    text: "Number of electrons in the outermost shell of Xe in XeF₂ is:",
                    options: ["8", "10", "12", "14"],
                    answer: 2,
                    explain: "Xe in XeF₂ has 8 original + 2 from F atoms + 2 lone pairs = 12 electrons in expanded octet"
                },
                {
                    id: "c27",
                    text: "Which is most paramagnetic?",
                    options: ["O₂", "O₂⁺", "O₂⁻", "O₂²⁻"],
                    answer: 1,
                    explain: "O₂⁺ has 15 electrons with 3 unpaired electrons, making it most paramagnetic"
                },
                {
                    id: "c28",
                    text: "Oxidation number of phosphorus in H₃PO₂ is:",
                    options: ["+1", "+2", "+3", "+5"],
                    answer: 0,
                    explain: "In hypophosphorous acid H₃PO₂: 3(+1) + P + 2(-2) = 0, solving gives P = +1"
                },
                {
                    id: "c29",
                    text: "Coordination isomerism is shown by:",
                    options: ["[Co(NH₃)₆][Cr(CN)₆]", "[Co(NH₃)₅Cl]SO₄", "[Co(en)₃]Cl₃", "[Pt(NH₃)₂Cl₂]"],
                    answer: 0,
                    explain: "[Co(NH₃)₆][Cr(CN)₆] can show coordination isomerism by interchange of ligands between metal centers"
                },
                {
                    id: "c30",
                    text: "Coordination number in CsCl structure is:",
                    options: ["6", "8", "4", "12"],
                    answer: 1,
                    explain: "In CsCl structure, each Cs⁺ is surrounded by 8 Cl⁻ ions and vice versa, giving coordination number 8"
                },
                {
                    id: "c31",
                    text: "Which has lowest boiling point?",
                    options: ["CH₄", "SiH₄", "GeH₄", "SnH₄"],
                    answer: 0,
                    explain: "CH₄ has weakest van der Waals forces due to smallest size and lowest polarizability"
                },
                {
                    id: "c32",
                    text: "Half-life of first-order reaction is:",
                    options: ["Inversely proportional to k", "Directly proportional to k", "Independent of k", "Proportional to k²"],
                    answer: 0,
                    explain: "For first-order reaction: t₁/₂ = 0.693/k, showing inverse proportionality to rate constant"
                },
                {
                    id: "c33",
                    text: "Which is linear molecule?",
                    options: ["SO₂", "H₂S", "BeF₂", "H₂O"],
                    answer: 2,
                    explain: "BeF₂ has sp hybridization with no lone pairs on Be, giving linear geometry"
                },
                {
                    id: "c34",
                    text: "Most electronegative element is:",
                    options: ["F", "O", "N", "Cl"],
                    answer: 0,
                    explain: "Fluorine is most electronegative element with Pauling electronegativity value of 4.0"
                },
                {
                    id: "c35",
                    text: "Aldol condensation product of acetaldehyde is:",
                    options: ["Crotonaldehyde", "Butanal", "Acetone", "Propanol"],
                    answer: 0,
                    explain: "Acetaldehyde undergoes aldol condensation to form β-hydroxybutanal, which dehydrates to crotonaldehyde"
                },
                {
                    id: "c36",
                    text: "Bond angle in ammonia is less than tetrahedral due to:",
                    options: ["Lone pair-bond pair repulsion", "Bond pair-bond pair repulsion", "Small size of N", "sp² hybridization"],
                    answer: 0,
                    explain: "Lone pair on nitrogen exerts greater repulsion than bond pairs, compressing H-N-H bond angle to ~107°"
                },
                {
                    id: "c37",
                    text: "Which element has highest first ionization energy in period 3?",
                    options: ["Na", "Mg", "Al", "Ar"],
                    answer: 3,
                    explain: "Argon has highest ionization energy in period 3 due to complete octet and highest nuclear charge"
                },
                {
                    id: "c38",
                    text: "Maximum number of isomers for [MA₃B₃] octahedral complex is:",
                    options: ["1", "2", "3", "4"],
                    answer: 1,
                    explain: "Octahedral [MA₃B₃] complex shows 2 geometrical isomers: facial (fac) and meridional (mer)"
                },
                {
                    id: "c39",
                    text: "Which metal has highest reduction potential?",
                    options: ["Cu", "Ag", "Au", "Pt"],
                    answer: 2,
                    explain: "Gold (Au³⁺/Au) has highest reduction potential (+1.50 V), making it most noble metal"
                },
                {
                    id: "c40",
                    text: "Hybridization of carbon in fullerene (C₆₀) is:",
                    options: ["sp", "sp²", "sp³", "sp³d"],
                    answer: 1,
                    explain: "Carbon atoms in fullerene C₆₀ are sp² hybridized forming curved surface with delocalized π electrons"
                },
                {
                    id: "c41",
                    text: "Which shows cis-trans isomerism?",
                    options: ["CH₂=CHCl", "CHCl=CHCl", "CH₃CH=CHCl", "CH₃CH=CH₂"],
                    answer: 1,
                    explain: "1,2-dichloroethene (CHCl=CHCl) can exist as cis and trans isomers around C=C double bond"
                },
                {
                    id: "c42",
                    text: "EAN of central metal in [Ni(CO)₄] is:",
                    options: ["18", "28", "36", "54"],
                    answer: 2,
                    explain: "Ni (28) + 4×2 (from CO) = 36 electrons, achieving krypton configuration"
                },
                {
                    id: "c43",
                    text: "Which is amphoteric hydroxide?",
                    options: ["Mg(OH)₂", "Al(OH)₃", "Ca(OH)₂", "Ba(OH)₂"],
                    answer: 1,
                    explain: "Al(OH)₃ can act as both acid and base, dissolving in both HCl and NaOH solutions"
                },
                {
                    id: "c44",
                    text: "Ethyne has how many σ bonds?",
                    options: ["2", "3", "4", "5"],
                    answer: 1,
                    explain: "C₂H₂ has 3 σ bonds: 2 C-H bonds and 1 C-C σ bond. The triple bond C≡C contains 1σ + 2π"
                },
                {
                    id: "c45",
                    text: "Which bond is strongest?",
                    options: ["C-C", "C=C", "C≡C", "C-H"],
                    answer: 2,
                    explain: "Carbon-carbon triple bond (C≡C) has highest bond energy (~835 kJ/mol) among carbon bonds"
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "In C₄ plants, CO₂ acceptor in mesophyll cells is:",
                    options: ["RuBP", "PEP", "Pyruvate", "Malate"],
                    answer: 1,
                    explain: "In C₄ plants, phosphoenolpyruvate (PEP) acts as CO₂ acceptor in mesophyll cells, forming oxaloacetic acid"
                },
                {
                    id: "b2",
                    text: "Aleurone layer in seeds is rich in:",
                    options: ["Carbohydrates", "Proteins", "Lipids", "Nucleic acids"],
                    answer: 1,
                    explain: "Aleurone layer contains protein bodies (aleurone grains) rich in storage proteins and hydrolytic enzymes"
                },
                {
                    id: "b3",
                    text: "Which hormone promotes root hair formation?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Ethylene"],
                    answer: 0,
                    explain: "Auxin promotes root hair elongation and development, enhancing water and mineral absorption capacity"
                },
                {
                    id: "b4",
                    text: "Cork cambium is also called:",
                    options: ["Fascicular cambium", "Phellogen", "Interfascicular cambium", "Pericycle"],
                    answer: 1,
                    explain: "Cork cambium (phellogen) produces cork cells (phellem) externally and secondary cortex (phelloderm) internally"
                },
                {
                    id: "b5",
                    text: "Which shows both isogamy and anisogamy?",
                    options: ["Spirogyra", "Chlamydomonas", "Ulothrix", "Zygnema"],
                    answer: 1,
                    explain: "Chlamydomonas exhibits both isogamous reproduction (morphologically similar gametes) and anisogamous reproduction"
                },
                {
                    id: "b6",
                    text: "Self-incompatibility prevents:",
                    options: ["Cross-pollination", "Self-pollination", "Fertilization", "Fruit formation"],
                    answer: 1,
                    explain: "Self-incompatibility mechanisms prevent self-fertilization, promoting outcrossing and genetic diversity"
                },
                {
                    id: "b7",
                    text: "CAM photosynthesis was first discovered in:",
                    options: ["Cactus", "Pineapple", "Sedum", "Crassulaceae"],
                    answer: 3,
                    explain: "Crassulacean Acid Metabolism (CAM) was first discovered in Crassulaceae family plants like Sedum"
                },
                {
                    id: "b8",
                    text: "Bulliform cells help in:",
                    options: ["Water storage", "Leaf rolling", "Photosynthesis", "Transpiration"],
                    answer: 1,
                    explain: "Bulliform cells in grass leaves lose water and shrink, causing leaf rolling to reduce transpiration"
                },
                {
                    id: "b9",
                    text: "Cyathium inflorescence is characteristic of:",
                    options: ["Euphorbia", "Ficus", "Dorstenia", "Artocarpus"],
                    answer: 0,
                    explain: "Cyathium is specialized inflorescence found in Euphorbia with involucral bracts and reduced flowers"
                },
                {
                    id: "b10",
                    text: "Polyembryony in Citrus is due to:",
                    options: ["Multiple fertilization", "Nucellar budding", "Cleavage of embryo", "Accessory embryo sacs"],
                    answer: 1,
                    explain: "Adventive polyembryony in Citrus results from nucellar cells developing into embryos without fertilization"
                },
                {
                    id: "b11",
                    text: "Skoog and Miller discovered:",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Ethylene"],
                    answer: 2,
                    explain: "Skoog and Miller discovered cytokinins (kinetin) while working on tissue culture of tobacco callus"
                },
                {
                    id: "b12",
                    text: "Photosystem I absorbs light maximally at:",
                    options: ["680 nm", "700 nm", "660 nm", "720 nm"],
                    answer: 1,
                    explain: "Photosystem I (PSI) has absorption maximum at 700 nm, hence called P700"
                },
                {
                    id: "b13",
                    text: "Short-day plants are actually:",
                    options: ["Day-length sensitive", "Night-length sensitive", "Light-quality sensitive", "Temperature sensitive"],
                    answer: 1,
                    explain: "Short-day plants are actually long-night plants, requiring dark period longer than critical duration"
                },
                {
                    id: "b14",
                    text: "Plasmolysis occurs when cell is placed in:",
                    options: ["Isotonic solution", "Hypotonic solution", "Hypertonic solution", "Pure water"],
                    answer: 2,
                    explain: "In hypertonic solution, water moves out of cell, causing protoplasm to shrink away from cell wall (plasmolysis)"
                },
                {
                    id: "b15",
                    text: "Pressure flow theory explains:",
                    options: ["Water transport", "Mineral transport", "Food translocation", "Gas exchange"],
                    answer: 2,
                    explain: "Münch's pressure flow theory explains translocation of organic solutes (mainly sucrose) in phloem"
                },
                {
                    id: "b16",
                    text: "Synergids help in:",
                    options: ["Nutrition of embryo sac", "Guidance of pollen tube", "Formation of endosperm", "Seed coat development"],
                    answer: 1,
                    explain: "Synergids secrete chemotropic substances that guide pollen tube toward embryo sac for fertilization"
                },
                {
                    id: "b17",
                    text: "Isobilateral leaves have:",
                    options: ["Stomata on upper surface only", "Stomata on lower surface only", "Stomata on both surfaces", "No stomata"],
                    answer: 2,
                    explain: "Isobilateral leaves have similar structure on both sides with stomata present on both upper and lower surfaces"
                },
                {
                    id: "b18",
                    text: "Tyloses are found in:",
                    options: ["Phloem", "Xylem vessels", "Cambium", "Cork"],
                    answer: 1,
                    explain: "Tyloses are balloon-like outgrowths of parenchyma cells that block xylem vessels, forming heartwood"
                },
                {
                    id: "b19",
                    text: "Carotenoids protect chlorophyll from:",
                    options: ["High temperature", "Photo-oxidation", "UV rays", "Dehydration"],
                    answer: 1,
                    explain: "Carotenoids act as accessory pigments and protect chlorophyll from photo-oxidative damage"
                },
                {
                    id: "b20",
                    text: "Porogamy refers to:",
                    options: ["Pollen tube entry through micropyle", "Pollen tube entry through chalaza", "Pollen tube entry through integuments", "Direct fertilization"],
                    answer: 0,
                    explain: "Porogamy is normal mode of pollen tube entry into ovule through micropyle opening"
                },
                {
                    id: "b21",
                    text: "Photomorphogenesis is controlled by:",
                    options: ["Chlorophyll", "Phytochrome", "Carotenoids", "Anthocyanins"],
                    answer: 1,
                    explain: "Phytochrome photoreceptor controls various photomorphogenic responses like seed germination, flowering"
                },
                {
                    id: "b22",
                    text: "Geitonogamy involves:",
                    options: ["Self-pollination within same flower", "Cross-pollination between different plants", "Pollination between different flowers of same plant", "Wind pollination"],
                    answer: 2,
                    explain: "Geitonogamy is transfer of pollen from anther of one flower to stigma of another flower on same plant"
                },
                {
                    id: "b23",
                    text: "Calvin cycle occurs in:",
                    options: ["Thylakoid membrane", "Thylakoid lumen", "Stroma", "Outer chloroplast membrane"],
                    answer: 2,
                    explain: "Calvin cycle (dark reaction) takes place in chloroplast stroma where CO₂ fixation occurs"
                },
                {
                    id: "b24",
                    text: "Hormone that breaks seed dormancy:",
                    options: ["ABA", "Gibberellin", "Auxin", "Cytokinin"],
                    answer: 1,
                    explain: "Gibberellins overcome seed dormancy by activating hydrolytic enzymes that mobilize stored food materials"
                },
                {
                    id: "b25",
                    text: "Vessels are absent in:",
                    options: ["Angiosperms", "Most gymnosperms", "Pteridophytes", "All of these"],
                    answer: 1,
                    explain: "Most gymnosperms lack vessels in their xylem and have only tracheids for water conduction"
                },
                {
                    id: "b26",
                    text: "Photorespiration involves:",
                    options: ["Only chloroplasts", "Only mitochondria", "Only peroxisomes", "All three organelles"],
                    answer: 3,
                    explain: "Photorespiration requires cooperation of chloroplasts, peroxisomes, and mitochondria for complete cycle"
                },
                {
                    id: "b27",
                    text: "Apical dominance is due to:",
                    options: ["Cytokinin", "Gibberellin", "Auxin", "Ethylene"],
                    answer: 2,
                    explain: "Auxin produced by apical bud inhibits growth of lateral buds, maintaining apical dominance"
                },
                {
                    id: "b28",
                    text: "Antipodals are:",
                    options: ["3 cells at micropylar end", "3 cells at chalazal end", "2 polar nuclei", "1 egg cell"],
                    answer: 1,
                    explain: "Three antipodal cells are located at chalazal end of embryo sac, opposite to egg apparatus"
                },
                {
                    id: "b29",
                    text: "Nitrogenase enzyme is sensitive to:",
                    options: ["Light", "Temperature", "Oxygen", "pH"],
                    answer: 2,
                    explain: "Nitrogenase is irreversibly inactivated by oxygen, requiring anaerobic conditions for nitrogen fixation"
                },
                {
                    id: "b30",
                    text: "Secondary wall in plant cells is composed of:",
                    options: ["Cellulose only", "Lignin only", "Cellulose and lignin", "Pectin only"],
                    answer: 2,
                    explain: "Secondary cell wall contains cellulose fibers embedded in lignin matrix, providing mechanical strength"
                },
                {
                    id: "b31",
                    text: "Sensitive plant (Mimosa pudica) shows:",
                    options: ["Phototropism", "Geotropism", "Thigmonasty", "Hydrotropism"],
                    answer: 2,
                    explain: "Mimosa pudica exhibits thigmonasty (touch-sensitive movement) through rapid turgor changes"
                },
                {
                    id: "b32",
                    text: "Annual rings are formed due to:",
                    options: ["Primary growth", "Seasonal activity of cambium", "Cork formation", "Leaf fall"],
                    answer: 1,
                    explain: "Annual rings result from seasonal variation in cambial activity producing different types of secondary xylem"
                },
                {
                    id: "b33",
                    text: "Abscission zone develops in:",
                    options: ["Stem", "Leaf base", "Root", "Flower"],
                    answer: 1,
                    explain: "Abscission zone forms at leaf base where separation layer develops for leaf fall"
                },
                {
                    id: "b34",
                    text: "Active absorption of water involves:",
                    options: ["Transpiration", "Root pressure", "Cohesion", "Adhesion"],
                    answer: 1,
                    explain: "Active absorption creates root pressure through active salt accumulation and osmotic water uptake"
                },
                {
                    id: "b35",
                    text: "CO₂ concentration point in C₄ plants:",
                    options: ["Mesophyll cells", "Bundle sheath cells", "Epidermal cells", "Guard cells"],
                    answer: 1,
                    explain: "Bundle sheath cells concentrate CO₂ around RuBisCO, minimizing photorespiration in C₄ plants"
                },
                {
                    id: "b36",
                    text: "Hydathodes are involved in:",
                    options: ["Transpiration", "Guttation", "Respiration", "Photosynthesis"],
                    answer: 1,
                    explain: "Hydathodes are water-secreting structures involved in guttation, especially during high humidity"
                },
                {
                    id: "b37",
                    text: "Opposite phyllotaxy is seen in:",
                    options: ["China rose", "Calotropis", "Alstonia", "Nerium"],
                    answer: 1,
                    explain: "Calotropis shows opposite phyllotaxy where two leaves arise from each node at opposite sides"
                },
                {
                    id: "b38",
                    text: "Bolting in plants is caused by:",
                    options: ["Auxin", "Cytokinin", "Gibberellin", "ABA"],
                    answer: 2,
                    explain: "Gibberellins cause bolting (rapid stem elongation) especially in rosette plants before flowering"
                },
                {
                    id: "b39",
                    text: "Succulent plants show:",
                    options: ["C₃ photosynthesis", "C₄ photosynthesis", "CAM photosynthesis", "Normal photosynthesis"],
                    answer: 2,
                    explain: "Most succulents use CAM photosynthesis to conserve water in arid environments"
                },
                {
                    id: "b40",
                    text: "Transpiration is mainly controlled by:",
                    options: ["Cuticle", "Stomata", "Lenticels", "Hydathodes"],
                    answer: 1,
                    explain: "Stomata regulate about 90% of transpiration through opening and closing mechanisms"
                },
                {
                    id: "b41",
                    text: "Totipotency in plants means:",
                    options: ["Unlimited growth", "Ability to form whole plant", "Unlimited division", "Ability to differentiate"],
                    answer: 1,
                    explain: "Totipotency is capacity of plant cell to develop into complete organism under suitable conditions"
                },
                {
                    id: "b42",
                    text: "Emerson enhancement effect demonstrates:",
                    options: ["Role of chlorophyll-a", "Cooperation of two photosystems", "Light saturation", "CO₂ fixation"],
                    answer: 1,
                    explain: "Emerson effect shows that PSI and PSII work together more efficiently than individually"
                },
                {
                    id: "b43",
                    text: "Double fertilization results in formation of:",
                    options: ["Embryo only", "Endosperm only", "Both embryo and endosperm", "Seed coat"],
                    answer: 2,
                    explain: "Double fertilization produces diploid embryo (egg + sperm) and triploid endosperm (polar nuclei + sperm)"
                },
                {
                    id: "b44",
                    text: "Heartwood is formed due to:",
                    options: ["Secondary growth", "Deposition of metabolites", "Loss of water", "All of these"],
                    answer: 3,
                    explain: "Heartwood formation involves all processes: secondary growth, metabolite deposition, and water loss"
                },
                {
                    id: "b45",
                    text: "2,4-D is a synthetic:",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "ABA"],
                    answer: 0,
                    explain: "2,4-Dichlorophenoxyacetic acid (2,4-D) is synthetic auxin used as herbicide and in tissue culture"
                },
                
                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Which bone forms the bridge of nose?",
                    options: ["Frontal", "Nasal", "Maxilla", "Ethmoid"],
                    answer: 1,
                    explain: "Nasal bones form the bony bridge of nose and support the upper part of nasal cavity"
                },
                {
                    id: "z2",
                    text: "Somatostatin is secreted by:",
                    options: ["Alpha cells", "Beta cells", "Delta cells", "PP cells"],
                    answer: 2,
                    explain: "Delta cells in pancreatic islets secrete somatostatin, which inhibits insulin and glucagon release"
                },
                {
                    id: "z3",
                    text: "Tissue thromboplastin is released by:",
                    options: ["Platelets", "Damaged tissues", "Plasma", "Liver"],
                    answer: 1,
                    explain: "Tissue thromboplastin (Factor III) is released from damaged blood vessels and tissues to initiate clotting"
                },
                {
                    id: "z4",
                    text: "Juxtaglomerular apparatus consists of:",
                    options: ["Glomerulus only", "Bowman's capsule only", "Macula densa and JG cells", "Entire nephron"],
                    answer: 2,
                    explain: "JGA includes macula densa (DCT cells) and juxtaglomerular cells (modified smooth muscle) that secrete renin"
                },
                {
                    id: "z5",
                    text: "Broca's area controls:",
                    options: ["Speech comprehension", "Speech production", "Memory", "Vision"],
                    answer: 1,
                    explain: "Broca's area in frontal lobe controls motor aspects of speech production and articulation"
                },
                {
                    id: "z6",
                    text: "Bilirubin is formed from breakdown of:",
                    options: ["Heme", "Globin", "Both heme and globin", "Iron"],
                    answer: 0,
                    explain: "Bilirubin is formed from breakdown of heme portion of hemoglobin by reticuloendothelial system"
                },
                {
                    id: "z7",
                    text: "Cholecystokinin stimulates release of:",
                    options: ["Pancreatic juice only", "Bile only", "Both pancreatic juice and bile", "Gastric juice"],
                    answer: 2,
                    explain: "CCK stimulates release of enzyme-rich pancreatic juice and contraction of gall bladder for bile release"
                },
                {
                    id: "z8",
                    text: "Bundle of His divides into:",
                    options: ["Two bundle branches", "Three bundle branches", "Four bundle branches", "Many Purkinje fibers"],
                    answer: 0,
                    explain: "Bundle of His divides into right and left bundle branches that conduct impulses to respective ventricles"
                },
                {
                    id: "z9",
                    text: "Conn's syndrome is due to excess of:",
                    options: ["Cortisol", "Aldosterone", "Adrenaline", "Growth hormone"],
                    answer: 1,
                    explain: "Conn's syndrome (primary hyperaldosteronism) results from excess aldosterone production by adrenal cortex"
                },
                {
                    id: "z10",
                    text: "Oocyte is arrested in:",
                    options: ["Prophase I", "Metaphase I", "Prophase II", "Metaphase II"],
                    answer: 3,
                    explain: "Secondary oocyte is arrested in metaphase II and completes meiosis only if fertilized by sperm"
                },
                {
                    id: "z11",
                    text: "Xerophthalmia is caused by deficiency of:",
                    options: ["Vitamin A", "Vitamin B₁", "Vitamin C", "Vitamin D"],
                    answer: 0,
                    explain: "Vitamin A deficiency causes xerophthalmia (dry eye) due to impaired mucus production by conjunctival cells"
                },
                {
                    id: "z12",
                    text: "Natural killer (NK) cells are part of:",
                    options: ["Adaptive immunity", "Innate immunity", "Humoral immunity", "Cell-mediated immunity only"],
                    answer: 1,
                    explain: "NK cells are lymphocytes that provide innate immunity by killing virus-infected and tumor cells"
                },
                {
                    id: "z13",
                    text: "Osmolality of urine is determined in:",
                    options: ["Proximal tubule", "Loop of Henle", "Distal tubule", "Collecting duct"],
                    answer: 3,
                    explain: "Final osmolality of urine is determined in collecting duct under ADH influence"
                },
                {
                    id: "z14",
                    text: "Relaxin is secreted during:",
                    options: ["Pregnancy", "Lactation", "Menstruation", "Ovulation"],
                    answer: 0,
                    explain: "Relaxin is secreted by corpus luteum and placenta during pregnancy to relax pelvic ligaments"
                },
                {
                    id: "z15",
                    text: "ABO blood group system was discovered by:",
                    options: ["Landsteiner", "Wiener", "Levine", "Stetson"],
                    answer: 0,
                    explain: "Karl Landsteiner discovered ABO blood group system in 1900, earning Nobel Prize in 1930"
                },
                {
                    id: "z16",
                    text: "Surfactant is produced by:",
                    options: ["Type I pneumocytes", "Type II pneumocytes", "Alveolar macrophages", "Clara cells"],
                    answer: 1,
                    explain: "Type II pneumocytes (septal cells) produce surfactant that reduces surface tension in alveoli"
                },
                {
                    id: "z17",
                    text: "Otoliths are present in:",
                    options: ["Cochlea", "Semicircular canals", "Utricle and saccule", "Middle ear"],
                    answer: 2,
                    explain: "Otoliths (calcium carbonate crystals) are present in utricle and saccule for detecting linear acceleration"
                },
                {
                    id: "z18",
                    text: "Clearance value is highest for:",
                    options: ["Inulin", "Creatinine", "PAH", "Glucose"],
                    answer: 2,
                    explain: "Para-aminohippuric acid (PAH) has highest clearance as it's both filtered and actively secreted"
                },
                {
                    id: "z19",
                    text: "Diabetes insipidus is caused by deficiency of:",
                    options: ["Insulin", "ADH", "Growth hormone", "Thyroid hormones"],
                    answer: 1,
                    explain: "Diabetes insipidus results from ADH deficiency or resistance, causing excessive water loss"
                },
                {
                    id: "z20",
                    text: "Capacitation of sperm occurs in:",
                    options: ["Epididymis", "Vas deferens", "Female reproductive tract", "Seminal vesicles"],
                    answer: 2,
                    explain: "Sperm capacitation occurs in female reproductive tract, enabling sperm to fertilize ovum"
                },
                {
                    id: "z21",
                    text: "Coronary circulation supplies:",
                    options: ["Heart muscle", "Brain", "Kidneys", "Liver"],
                    answer: 0,
                    explain: "Coronary arteries arise from aortic root and supply oxygenated blood to heart muscle (myocardium)"
                },
                {
                    id: "z22",
                    text: "Bohr effect relates to:",
                    options: ["Oxygen transport", "CO₂ transport", "Both oxygen and CO₂ transport", "Blood pH"],
                    answer: 0,
                    explain: "Bohr effect describes decreased oxygen affinity of hemoglobin in presence of increased CO₂ and H⁺"
                },
                {
                    id: "z23",
                    text: "Tetany is caused by deficiency of:",
                    options: ["Sodium", "Potassium", "Calcium", "Magnesium"],
                    answer: 2,
                    explain: "Hypocalcemia causes tetany due to increased neuromuscular excitability and muscle spasms"
                },
                {
                    id: "z24",
                    text: "Hering-Breuer reflex prevents:",
                    options: ["Hyperventilation", "Hypoventilation", "Over-inflation of lungs", "Apnea"],
                    answer: 2,
                    explain: "Hering-Breuer reflex inhibits inspiration when lungs are stretched, preventing over-inflation"
                },
                {
                    id: "z25",
                    text: "Satellite cells are associated with:",
                    options: ["Cardiac muscle", "Skeletal muscle", "Smooth muscle", "Nervous tissue"],
                    answer: 1,
                    explain: "Satellite cells are stem cells associated with skeletal muscle fibers for muscle repair and growth"
                },
                {
                    id: "z26",
                    text: "Ito cells are found in:",
                    options: ["Liver", "Kidney", "Lung", "Spleen"],
                    answer: 0,
                    explain: "Ito cells (hepatic stellate cells) store vitamin A and produce collagen in liver during fibrosis"
                },
                {
                    id: "z27",
                    text: "Pellagra is caused by deficiency of:",
                    options: ["Thiamine", "Riboflavin", "Niacin", "Pyridoxine"],
                    answer: 2,
                    explain: "Pellagra results from niacin (vitamin B₃) deficiency, characterized by 4 D's: dermatitis, diarrhea, dementia, death"
                },
                {
                    id: "z28",
                    text: "Calcium ions bind to:",
                    options: ["Actin", "Myosin", "Tropomyosin", "Troponin"],
                    answer: 3,
                    explain: "Calcium ions bind to troponin C, causing conformational change that moves tropomyosin and exposes myosin-binding sites"
                },
                {
                    id: "z29",
                    text: "Lipolysis is stimulated by:",
                    options: ["Insulin", "Glucagon", "Cortisol", "Both glucagon and cortisol"],
                    answer: 3,
                    explain: "Both glucagon and cortisol stimulate lipolysis (breakdown of fats) for energy during fasting states"
                },
                {
                    id: "z30",
                    text: "M cells are found in:",
                    options: ["Peyer's patches", "Spleen", "Thymus", "Lymph nodes"],
                    answer: 0,
                    explain: "M cells (microfold cells) in Peyer's patches transport antigens from gut lumen to immune cells"
                },
                {
                    id: "z31",
                    text: "Brown adipose tissue generates heat through:",
                    options: ["Shivering", "Non-shivering thermogenesis", "Muscle contraction", "Vasoconstriction"],
                    answer: 1,
                    explain: "Brown fat produces heat through non-shivering thermogenesis via uncoupling protein (thermogenin)"
                },
                {
                    id: "z32",
                    text: "Quickening refers to:",
                    options: ["Fertilization", "Implantation", "First fetal movements felt", "Birth"],
                    answer: 2,
                    explain: "Quickening is first perception of fetal movements by mother, typically around 18-20 weeks"
                },
                {
                    id: "z33",
                    text: "Erythropoietin is produced mainly by:",
                    options: ["Bone marrow", "Liver", "Kidneys", "Spleen"],
                    answer: 2,
                    explain: "About 90% of erythropoietin is produced by peritubular interstitial cells in kidneys"
                },
                {
                    id: "z34",
                    text: "Suprachiasmatic nucleus regulates:",
                    options: ["Body temperature", "Circadian rhythms", "Blood pressure", "Appetite"],
                    answer: 1,
                    explain: "Suprachiasmatic nucleus in hypothalamus is master biological clock controlling circadian rhythms"
                },
                {
                    id: "z35",
                    text: "Lens is suspended by:",
                    options: ["Ciliary muscle", "Suspensory ligaments", "Iris", "Choroid"],
                    answer: 1,
                    explain: "Lens is suspended by suspensory ligaments (zonule of Zinn) attached to ciliary body"
                },
                {
                    id: "z36",
                    text: "Obligatory water reabsorption occurs in:",
                    options: ["Proximal tubule", "Loop of Henle", "Distal tubule", "Collecting duct"],
                    answer: 0,
                    explain: "About 65% of filtered water is obligatorily reabsorbed in proximal tubule along with solutes"
                },
                {
                    id: "z37",
                    text: "Umami taste is detected by:",
                    options: ["Sweet receptors", "Glutamate receptors", "Salt receptors", "Bitter receptors"],
                    answer: 1,
                    explain: "Umami (savory) taste is detected by glutamate receptors responsive to amino acids like glutamate"
                },
                {
                    id: "z38",
                    text: "Dawn phenomenon in diabetes is due to:",
                    options: ["Evening insulin", "Morning cortisol surge", "Late dinner", "Exercise"],
                    answer: 1,
                    explain: "Dawn phenomenon results from early morning surge of cortisol and growth hormone causing hyperglycemia"
                },
                {
                    id: "z39",
                    text: "Decidual reaction occurs in:",
                    options: ["Ovary", "Fallopian tube", "Endometrium", "Cervix"],
                    answer: 2,
                    explain: "Decidual reaction involves transformation of endometrial stromal cells during implantation"
                },
                {
                    id: "z40",
                    text: "True ribs directly attach to:",
                    options: ["Vertebrae", "Sternum", "Other ribs", "Clavicle"],
                    answer: 1,
                    explain: "First seven pairs of ribs (true ribs) directly attach to sternum through their own costal cartilages"
                },
                {
                    id: "z41",
                    text: "Castle's intrinsic factor is necessary for absorption of:",
                    options: ["Vitamin B₁", "Vitamin B₆", "Vitamin B₁₂", "Folic acid"],
                    answer: 2,
                    explain: "Intrinsic factor from parietal cells is essential for vitamin B₁₂ absorption in terminal ileum"
                },
                {
                    id: "z42",
                    text: "Synovial fluid is produced by:",
                    options: ["Cartilage", "Synovial membrane", "Bone", "Ligaments"],
                    answer: 1,
                    explain: "Synovial membrane (synovium) lining joint capsule produces synovial fluid for lubrication"
                },
                {
                    id: "z43",
                    text: "TSH releasing hormone is produced by:",
                    options: ["Anterior pituitary", "Posterior pituitary", "Hypothalamus", "Thyroid"],
                    answer: 2,
                    explain: "TRH (thyrotropin-releasing hormone) is produced by hypothalamus to stimulate TSH release from anterior pituitary"
                },
                {
                    id: "z44",
                    text: "Papillary muscles prevent:",
                    options: ["Valve stenosis", "Valve regurgitation", "Arrhythmias", "Heart failure"],
                    answer: 1,
                    explain: "Papillary muscles contract during systole to prevent AV valve regurgitation via chordae tendineae"
                },
                {
                    id: "z45",
                    text: "Zona reaction prevents:",
                    options: ["Fertilization", "Polyspermy", "Implantation", "Ovulation"],
                    answer: 1,
                    explain: "Zona reaction involves hardening of zona pellucida after sperm penetration to prevent polyspermy"
                }
            ]
        }
    ]
};
