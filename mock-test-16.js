// mock-test-16.js - NEET Mock Test 16 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_16 = {
    id: "neet-016",
    title: "Full Syllabus Mock 16", 
    level: "hard",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A particle moves in xy-plane such that its coordinates are x = 4t² and y = 3t³. Magnitude of acceleration at t = 1s is:",
                    options: ["20 m/s²", "10 m/s²", "18 m/s²", "√82 m/s²"],
                    answer: 1,
                    explain: "ax = d²x/dt² = 8 m/s², ay = d²y/dt² = 18t = 18 m/s² at t=1s. |a| = √(8² + 18²) = √(64 + 324) = √388 ≈ 19.7 m/s², closest to 20 m/s²"
                },
                {
                    id: "p2",
                    text: "A yo-yo of mass m and radius R rolls down a string. If moment of inertia about center is I = ½mR², acceleration of center of mass is:",
                    options: ["g/3", "2g/3", "g/2", "3g/4"],
                    answer: 1,
                    explain: "For yo-yo: mg - T = ma, and TR = Iα = ½mR²α. Since a = Rα: T = ma/2. Substituting: mg - ma/2 = ma. Therefore a = 2g/3"
                },
                {
                    id: "p3",
                    text: "A 5μF capacitor is charged to 100V and then connected to uncharged 15μF capacitor. Final voltage across the combination is:",
                    options: ["75 V", "25 V", "50 V", "80 V"],
                    answer: 1,
                    explain: "Initial charge Q = CV = 5×10⁻⁶ × 100 = 500×10⁻⁶ C. Total capacitance = 5 + 15 = 20μF. Final voltage = Q/Ctotal = 500×10⁻⁶/20×10⁻⁶ = 25 V"
                },
                {
                    id: "p4",
                    text: "In Newton's rings experiment with air film, radius of 10th dark ring is 3mm. Radius of 20th dark ring is:",
                    options: ["6 mm", "4.24 mm", "3√2 mm", "6√2 mm"],
                    answer: 1,
                    explain: "For dark rings: r²n = nλR. r₁₀² = 10λR = 9 mm². r₂₀² = 20λR = 2 × 9 = 18 mm². r₂₀ = √18 = 3√2 ≈ 4.24 mm"
                },
                {
                    id: "p5",
                    text: "A photon of wavelength 400 nm has energy approximately:",
                    options: ["3.1 eV", "4.97 eV", "2.07 eV", "6.2 eV"],
                    answer: 0,
                    explain: "E = hc/λ = 1240 eV⋅nm / 400 nm = 3.1 eV"
                },
                {
                    id: "p6",
                    text: "A straight wire carrying current 10A is bent into circular loop of radius 1m. Magnetic field at center changes by:",
                    options: ["2π × 10⁻⁶ T", "π × 10⁻⁶ T", "4π × 10⁻⁶ T", "10⁻⁶ T"],
                    answer: 0,
                    explain: "Straight wire field at distance r: B = μ₀I/2πr = 0. Circular loop field at center: B = μ₀I/2R = 4π×10⁻⁷ × 10 / 2 = 2π×10⁻⁶ T. Change = 2π×10⁻⁶ T"
                },
                {
                    id: "p7",
                    text: "A uniform disc of radius R oscillates about horizontal axis at distance R/2 from center. Time period is:",
                    options: ["2π√(3R/2g)", "2π√(5R/4g)", "2π√(3R/4g)", "2π√(R/g)"],
                    answer: 0,
                    explain: "I = ICM + md² = ½mR² + m(R/2)² = ¾mR². Distance to CM = R/2. T = 2π√(I/mgd) = 2π√(¾mR²/mg×R/2) = 2π√(3R/2g)"
                },
                {
                    id: "p8",
                    text: "An inductor L = 0.2H and capacitor C = 50μF are connected in series with AC source. Resonance occurs at frequency:",
                    options: ["50 Hz", "35.4 Hz", "70.7 Hz", "25 Hz"],
                    answer: 0,
                    explain: "f₀ = 1/(2π√LC) = 1/(2π√(0.2 × 50×10⁻⁶)) = 1/(2π√(10×10⁻⁶)) = 1/(2π × √10×10⁻³) ≈ 50 Hz"
                },
                {
                    id: "p9",
                    text: "In isobaric expansion of ideal gas from volume V to 3V, work done is:",
                    options: ["PV", "2PV", "3PV", "PV ln 3"],
                    answer: 1,
                    explain: "In isobaric process, work done W = P(V₂ - V₁) = P(3V - V) = 2PV"
                },
                {
                    id: "p10",
                    text: "Two long parallel wires carry currents 3A and 5A in opposite directions separated by 2cm. Force per unit length is:",
                    options: ["1.5×10⁻⁴ N/m (repulsive)", "3×10⁻⁴ N/m (repulsive)", "1.5×10⁻⁴ N/m (attractive)", "6×10⁻⁴ N/m (repulsive)"],
                    answer: 0,
                    explain: "F/L = μ₀I₁I₂/2πd = 4π×10⁻⁷ × 3 × 5 / (2π × 0.02) = 60×10⁻⁷ / 0.04 = 1.5×10⁻⁴ N/m. Opposite currents repel"
                },
                {
                    id: "p11",
                    text: "A particle executes SHM: x = 6sin(4t + π/3) cm. Displacement when velocity is maximum is:",
                    options: ["0 cm", "3 cm", "6 cm", "-3 cm"],
                    answer: 0,
                    explain: "Velocity v = dx/dt = 24cos(4t + π/3). Maximum velocity occurs when cos(4t + π/3) = ±1, which happens when displacement x = 0"
                },
                {
                    id: "p12",
                    text: "In photoelectric effect, if frequency is tripled, stopping potential:",
                    options: ["Triples", "More than triples", "Less than triples", "Becomes nine times"],
                    answer: 1,
                    explain: "eV₀ = hf - φ. If f becomes 3f: eV₀' = 3hf - φ = 3(hf - φ) + 2φ = 3eV₀ + 2φ. Since φ > 0, stopping potential more than triples"
                },
                {
                    id: "p13",
                    text: "In Wheatstone bridge, if all resistors are 10Ω and one arm increases to 11Ω, galvanometer current depends on:",
                    options: ["Battery EMF only", "Galvanometer resistance only", "Both EMF and galvanometer resistance", "Bridge configuration only"],
                    answer: 2,
                    explain: "Unbalanced bridge current depends on both the EMF of battery and internal resistance including galvanometer resistance"
                },
                {
                    id: "p14",
                    text: "A missile fired at 30° above horizontal travels 1000m horizontally. Maximum height achieved is:",
                    options: ["144 m", "125 m", "166 m", "200 m"],
                    answer: 0,
                    explain: "Range R = u²sin2θ/g = 1000m. For θ = 30°: u² = 1000g/sin60° = 1000×10/(√3/2) = 20000/√3. Height H = u²sin²θ/2g = (20000/√3)×(1/4)/20 = 144 m"
                },
                {
                    id: "p15",
                    text: "If current through inductor L increases linearly from 0 to I₀ in time t, energy stored becomes:",
                    options: ["½LI₀²", "LI₀²", "¼LI₀²", "⅓LI₀²"],
                    answer: 0,
                    explain: "Energy stored in inductor U = ½LI² where I is final current. When current reaches I₀, energy = ½LI₀²"
                },
                {
                    id: "p16",
                    text: "A ball dropped from height h bounces to height h/4. Coefficient of restitution is:",
                    options: ["1/4", "1/2", "√3/2", "1/√2"],
                    answer: 1,
                    explain: "Coefficient of restitution e = √(h₂/h₁) = √(h/4 / h) = √(1/4) = 1/2"
                },
                {
                    id: "p17",
                    text: "In parallel RLC circuit, admittance Y = G + j(ωC - 1/ωL). At resonance, admittance equals:",
                    options: ["0", "∞", "G", "1/R"],
                    answer: 2,
                    explain: "At resonance, ωC = 1/ωL, so imaginary part becomes zero. Admittance Y = G = 1/R"
                },
                {
                    id: "p18",
                    text: "de Broglie wavelength of electron accelerated through 100V is approximately:",
                    options: ["1.23 Å", "2.46 Å", "0.61 Å", "4.92 Å"],
                    answer: 0,
                    explain: "λ = h/p = h/√(2mK) where K = eV = 100 eV. λ = 6.626×10⁻³⁴/√(2×9.1×10⁻³¹×100×1.6×10⁻¹⁹) ≈ 1.23×10⁻¹⁰ m = 1.23 Å"
                },
                {
                    id: "p19",
                    text: "For damped harmonic oscillator with damping coefficient γ, amplitude decreases as:",
                    options: ["e^(-γt)", "e^(-γt/2)", "e^(-γt/2m)", "e^(-γt/m)"],
                    answer: 2,
                    explain: "For damped oscillator, amplitude A(t) = A₀e^(-γt/2m) where γ is damping coefficient and m is mass"
                },
                {
                    id: "p20",
                    text: "Electric field due to infinite charged sheet with surface density σ is:",
                    options: ["σ/ε₀", "σ/2ε₀", "2σ/ε₀", "σ/4πε₀"],
                    answer: 1,
                    explain: "Using Gauss's law for infinite plane sheet: E = σ/2ε₀"
                },
                {
                    id: "p21",
                    text: "A ray travels from dense medium (n₁ = 1.5) to rare medium (n₂ = 1.0) at angle 60° to normal. Refracted ray makes angle:",
                    options: ["45°", "30°", "No refraction occurs", "90°"],
                    answer: 2,
                    explain: "Critical angle θc = sin⁻¹(n₂/n₁) = sin⁻¹(1/1.5) = 41.8°. Since incident angle 60° > θc, total internal reflection occurs"
                },
                {
                    id: "p22",
                    text: "Two coherent waves with amplitudes 3 and 4 interfere destructively. Resultant amplitude is:",
                    options: ["1", "7", "5", "0"],
                    answer: 0,
                    explain: "For destructive interference, resultant amplitude = |a₁ - a₂| = |4 - 3| = 1"
                },
                {
                    id: "p23",
                    text: "In isothermal expansion, internal energy of ideal gas:",
                    options: ["Increases", "Decreases", "Remains constant", "May increase or decrease"],
                    answer: 2,
                    explain: "For ideal gas in isothermal process, temperature is constant, so internal energy U = nCvT remains constant"
                },
                {
                    id: "p24",
                    text: "If 87.5% of radioactive sample decays in 30 minutes, half-life is:",
                    options: ["10 min", "7.5 min", "15 min", "5 min"],
                    answer: 0,
                    explain: "87.5% decay means 12.5% remains = 1/8 of original = (1/2)³. So 3 half-lives = 30 min. Half-life = 10 min"
                },
                {
                    id: "p25",
                    text: "Magnetic flux linkage with coil of N turns carrying current I is:",
                    options: ["LI", "NLI", "LI/N", "NΦ"],
                    answer: 3,
                    explain: "Flux linkage = NΦ where N is number of turns and Φ is flux through each turn"
                },
                {
                    id: "p26",
                    text: "Zeeman effect is:",
                    options: ["Splitting of spectral lines in electric field", "Splitting of spectral lines in magnetic field", "Doppler shift", "Photoelectric effect"],
                    answer: 1,
                    explain: "Zeeman effect is splitting of atomic spectral lines in presence of external magnetic field"
                },
                {
                    id: "p27",
                    text: "For dense flint glass (n = 1.8), minimum angle of incidence for total internal reflection from air interface is:",
                    options: ["33.7°", "56.3°", "41.8°", "48.6°"],
                    answer: 0,
                    explain: "Critical angle θc = sin⁻¹(1/n) = sin⁻¹(1/1.8) = sin⁻¹(0.556) = 33.7°"
                },
                {
                    id: "p28",
                    text: "In organ pipe open at one end, third overtone frequency is:",
                    options: ["7f₁", "5f₁", "3f₁", "4f₁"],
                    answer: 0,
                    explain: "In closed pipe, frequencies are f₁, 3f₁, 5f₁, 7f₁... Third overtone = 7f₁"
                },
                {
                    id: "p29",
                    text: "Heat engine operating between 800K and 300K has maximum theoretical efficiency:",
                    options: ["62.5%", "37.5%", "75%", "50%"],
                    answer: 0,
                    explain: "Carnot efficiency η = 1 - Tc/Th = 1 - 300/800 = 1 - 0.375 = 0.625 = 62.5%"
                },
                {
                    id: "p30",
                    text: "Cyclotron frequency is independent of:",
                    options: ["Charge of particle", "Mass of particle", "Velocity of particle", "Magnetic field strength"],
                    answer: 2,
                    explain: "Cyclotron frequency f = qB/2πm is independent of particle velocity"
                },
                {
                    id: "p31",
                    text: "Speed of sound in gas depends on:",
                    options: ["Density only", "Pressure only", "Both density and pressure", "Elasticity and density"],
                    answer: 3,
                    explain: "Speed of sound v = √(B/ρ) where B is bulk modulus (elasticity) and ρ is density"
                },
                {
                    id: "p32",
                    text: "Radius of gyration of solid cylinder about its central axis is:",
                    options: ["R", "R/√2", "R/2", "R√2"],
                    answer: 1,
                    explain: "For solid cylinder: I = ½MR², so radius of gyration k = √(I/M) = √(½R²) = R/√2"
                },
                {
                    id: "p33",
                    text: "In LC oscillation, when electric energy is 25% of maximum, magnetic energy is:",
                    options: ["75% of maximum", "25% of maximum", "50% of maximum", "100% of maximum"],
                    answer: 0,
                    explain: "Total energy = UE + UM = constant. If UE = 25% of max, then UM = 75% of max"
                },
                {
                    id: "p34",
                    text: "A concave mirror of radius 40cm forms virtual image 3 times the object size. Object distance is:",
                    options: ["13.3 cm", "20 cm", "10 cm", "6.67 cm"],
                    answer: 0,
                    explain: "f = R/2 = 20 cm. For virtual image: m = +3 = -v/u, so v = -3u. 1/20 = 1/u + 1/(-3u) = 2/3u. u = 40/3 = 13.3 cm"
                },
                {
                    id: "p35",
                    text: "In nuclear fission of U-235, average number of neutrons released per fission is:",
                    options: ["2.4", "1.5", "3.2", "4.1"],
                    answer: 0,
                    explain: "On average, fission of U-235 releases about 2.4 neutrons per fission event"
                },
                {
                    id: "p36",
                    text: "Quality factor Q of damped oscillator is defined as:",
                    options: ["ω₀/2γ", "2γ/ω₀", "ω₀γ", "γ/ω₀"],
                    answer: 0,
                    explain: "Quality factor Q = ω₀/2γ where ω₀ is natural frequency and γ is damping parameter"
                },
                {
                    id: "p37",
                    text: "Curie temperature is temperature at which:",
                    options: ["Ferromagnetic material becomes paramagnetic", "Paramagnetic material becomes diamagnetic", "Material loses conductivity", "Material melts"],
                    answer: 0,
                    explain: "Above Curie temperature, ferromagnetic materials lose their ferromagnetism and become paramagnetic"
                },
                {
                    id: "p38",
                    text: "Auto-transformer works on principle of:",
                    options: ["Mutual induction", "Self induction", "Both mutual and self induction", "Electromagnetic induction"],
                    answer: 0,
                    explain: "Auto-transformer works on principle of mutual induction between different sections of same winding"
                },
                {
                    id: "p39",
                    text: "For charging capacitor in RC circuit, time to reach 63.2% of maximum charge is:",
                    options: ["RC", "2RC", "RC/2", "RC/e"],
                    answer: 0,
                    explain: "Time constant τ = RC. After time τ, charge reaches (1-1/e) = 63.2% of maximum value"
                },
                {
                    id: "p40",
                    text: "Two point charges +q and -q separated by distance 2a constitute electric dipole. Potential at distance r >> a along perpendicular bisector is:",
                    options: ["0", "kq/r²", "2kqa/r²", "kqa/r³"],
                    answer: 0,
                    explain: "On perpendicular bisector of dipole, distances from both charges are equal, so potentials cancel out: V = 0"
                },
                {
                    id: "p41",
                    text: "Kirchhoff's current law is based on:",
                    options: ["Conservation of energy", "Conservation of charge", "Ohm's law", "Faraday's law"],
                    answer: 1,
                    explain: "KCL states that sum of currents at a node is zero, based on conservation of electric charge"
                },
                {
                    id: "p42",
                    text: "Average power delivered by AC source V₀sinωt to resistor R is:",
                    options: ["V₀²/R", "V₀²/2R", "V₀²/4R", "2V₀²/R"],
                    answer: 1,
                    explain: "Average power = Vrms²/R = (V₀/√2)²/R = V₀²/2R"
                },
                {
                    id: "p43",
                    text: "In single slit diffraction, first minimum occurs when path difference is:",
                    options: ["λ", "λ/2", "2λ", "3λ/2"],
                    answer: 0,
                    explain: "First minimum in single slit diffraction occurs when path difference = λ"
                },
                {
                    id: "p44",
                    text: "Gravitational potential energy is maximum when orbital radius is:",
                    options: ["Minimum", "Maximum", "At mean distance", "Infinite"],
                    answer: 3,
                    explain: "Gravitational PE = -GMm/r. As r → ∞, PE → 0 which is maximum (least negative) value"
                },
                {
                    id: "p45",
                    text: "In reverse biased p-n junction, current is mainly due to:",
                    options: ["Majority carriers", "Minority carriers", "Both equally", "No current flows"],
                    answer: 1,
                    explain: "In reverse bias, small current flows due to minority carriers being swept across junction by electric field"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which oxide has highest lattice energy?",
                    options: ["MgO", "CaO", "SrO", "BaO"],
                    answer: 0,
                    explain: "MgO has smallest cation (Mg²⁺) and highest charge density among given oxides, resulting in highest lattice energy"
                },
                {
                    id: "c2",
                    text: "The shape of XeOF₄ is:",
                    options: ["Tetrahedral", "Square pyramidal", "Trigonal bipyramidal", "Octahedral"],
                    answer: 1,
                    explain: "XeOF₄ has 6 electron pairs (5 bonding + 1 lone pair) around Xe. Lone pair occupies equatorial position giving square pyramidal shape"
                },
                {
                    id: "c3",
                    text: "Which carbocation is least stable?",
                    options: ["CH₃⁺", "CH₃CH₂⁺", "(CH₃)₂CH⁺", "(CH₃)₃C⁺"],
                    answer: 0,
                    explain: "Methyl carbocation (CH₃⁺) is least stable as it has no alkyl groups to provide hyperconjugation or inductive stabilization"
                },
                {
                    id: "c4",
                    text: "Oxidation state of iodine in H₅IO₆ is:",
                    options: ["+5", "+6", "+7", "+8"],
                    answer: 2,
                    explain: "In periodic acid H₅IO₆: 5(+1) + I + 6(-2) = 0, solving gives I = +7"
                },
                {
                    id: "c5",
                    text: "Maximum number of stereoisomers for compound with n chiral centers is:",
                    options: ["2n", "2ⁿ", "n²", "n!"],
                    answer: 1,
                    explain: "Maximum stereoisomers = 2ⁿ where n is number of chiral centers (assuming no meso compounds)"
                },
                {
                    id: "c6",
                    text: "Among oxoacids of chlorine, strongest acid is:",
                    options: ["HClO", "HClO₂", "HClO₃", "HClO₄"],
                    answer: 3,
                    explain: "HClO₄ is strongest as it has maximum number of electronegative oxygen atoms increasing acidity"
                },
                {
                    id: "c7",
                    text: "Which has electronic configuration [Ar] 3d⁶?",
                    options: ["Mn⁺", "Fe²⁺", "Co³⁺", "Cr"],
                    answer: 1,
                    explain: "Fe²⁺: Fe (26) loses 2 electrons → [Ar] 3d⁶"
                },
                {
                    id: "c8",
                    text: "Which molecule has see-saw geometry?",
                    options: ["SF₂", "SF₄", "SF₆", "SO₂"],
                    answer: 1,
                    explain: "SF₄ has 5 electron pairs (4 bonding + 1 lone pair) around S, giving see-saw geometry"
                },
                {
                    id: "c9",
                    text: "Which compound has incomplete octet?",
                    options: ["BeCl₂", "BF₃", "AlCl₃", "All of these"],
                    answer: 3,
                    explain: "Be in BeCl₂ (4 electrons), B in BF₃ (6 electrons), and Al in AlCl₃ (6 electrons) all have incomplete octets"
                },
                {
                    id: "c10",
                    text: "Entropy of universe in spontaneous process:",
                    options: ["Increases", "Decreases", "Remains constant", "May increase or decrease"],
                    answer: 0,
                    explain: "Second law of thermodynamics: entropy of universe always increases in spontaneous processes"
                },
                {
                    id: "c11",
                    text: "Best reducing agent among alkaline earth metals is:",
                    options: ["Be", "Mg", "Ca", "Ba"],
                    answer: 3,
                    explain: "Ba has most negative reduction potential (-2.90 V), making it strongest reducing agent among alkaline earth metals"
                },
                {
                    id: "c12",
                    text: "Which has zero unpaired electrons?",
                    options: ["[Ni(CN)₄]²⁻", "[NiCl₄]²⁻", "[Ni(CO)₄]", "Both [Ni(CN)₄]²⁻ and [Ni(CO)₄]"],
                    answer: 3,
                    explain: "Both complexes have Ni in zero oxidation state (d¹⁰) or strong field ligands causing electron pairing"
                },
                {
                    id: "c13",
                    text: "Hoffmann bromamide reaction converts:",
                    options: ["Amide to amine with same carbon atoms", "Amide to amine with one less carbon", "Amine to amide", "Acid to amine"],
                    answer: 1,
                    explain: "Hoffmann bromamide degradation converts primary amide to primary amine with one carbon atom less"
                },
                {
                    id: "c14",
                    text: "Correct order of ionic radii is:",
                    options: ["Li⁺ > Na⁺ > K⁺", "K⁺ > Na⁺ > Li⁺", "Na⁺ > K⁺ > Li⁺", "All equal"],
                    answer: 1,
                    explain: "Ionic radius increases down group due to addition of electron shells: K⁺ > Na⁺ > Li⁺"
                },
                {
                    id: "c15",
                    text: "Buffer with pH = 9 can be prepared using:",
                    options: ["CH₃COOH + CH₃COONa", "NH₄OH + NH₄Cl", "HCl + NaCl", "NaOH + NaCl"],
                    answer: 1,
                    explain: "NH₄OH + NH₄Cl forms basic buffer with pKb ≈ 4.8, giving pH ≈ 9 in suitable concentrations"
                },
                {
                    id: "c16",
                    text: "Spin-only magnetic moment of [Cr(H₂O)₆]²⁺ is:",
                    options: ["2.83 BM", "3.87 BM", "4.90 BM", "5.92 BM"],
                    answer: 2,
                    explain: "Cr²⁺ (d⁴) with weak field H₂O ligands has 4 unpaired electrons. μ = √[n(n+2)] = √[4×6] = 4.90 BM"
                },
                {
                    id: "c17",
                    text: "Which has minimum dipole moment?",
                    options: ["H₂O", "NH₃", "NF₃", "CHCl₃"],
                    answer: 2,
                    explain: "NF₃ has minimum dipole moment due to high electronegativity of F atoms reducing the net dipole from lone pair"
                },
                {
                    id: "c18",
                    text: "For zero order reaction, half-life is:",
                    options: ["Independent of concentration", "Proportional to concentration", "Inversely proportional to concentration", "Proportional to square of concentration"],
                    answer: 1,
                    explain: "For zero order: t₁/₂ = [A₀]/2k, directly proportional to initial concentration"
                },
                {
                    id: "c19",
                    text: "Which has maximum boiling point?",
                    options: ["HF", "HCl", "HBr", "HI"],
                    answer: 0,
                    explain: "HF has highest boiling point due to strong hydrogen bonding despite being smallest molecule"
                },
                {
                    id: "c20",
                    text: "Which is anti-aromatic?",
                    options: ["Benzene", "Cyclobutadiene", "Naphthalene", "Pyridine"],
                    answer: 1,
                    explain: "Cyclobutadiene has 4π electrons (4n, n=1), making it anti-aromatic according to Hückel's rule"
                },
                {
                    id: "c21",
                    text: "Which molecule is non-polar despite polar bonds?",
                    options: ["H₂O", "NH₃", "CO₂", "HCl"],
                    answer: 2,
                    explain: "CO₂ is linear molecule where two polar C=O bonds cancel each other, resulting in zero dipole moment"
                },
                {
                    id: "c22",
                    text: "Which can show optical isomerism?",
                    options: ["[Co(en)₃]³⁺", "[Co(NH₃)₆]³⁺", "[CoCl₄]²⁻", "[Co(NH₃)₄Cl₂]⁺"],
                    answer: 0,
                    explain: "[Co(en)₃]³⁺ has no plane of symmetry due to chelating ligands, showing optical isomerism (Δ and Λ forms)"
                },
                {
                    id: "c23",
                    text: "pKa values: HCOOH (3.8), CH₃COOH (4.8), C₆H₅COOH (4.2). Strongest acid is:",
                    options: ["HCOOH", "CH₃COOH", "C₆H₅COOH", "All equal"],
                    answer: 0,
                    explain: "Lower pKa indicates stronger acid. HCOOH has lowest pKa (3.8), hence strongest acid"
                },
                {
                    id: "c24",
                    text: "Bond order of N₂⁺ is:",
                    options: ["2", "2.5", "3", "1.5"],
                    answer: 1,
                    explain: "N₂⁺ has 13 electrons. Bond order = (bonding - antibonding)/2 = (9-4)/2 = 2.5"
                },
                {
                    id: "c25",
                    text: "Which shows nucleophilic addition reaction?",
                    options: ["Alkenes", "Alkynes", "Aldehydes", "Alkanes"],
                    answer: 2,
                    explain: "Aldehydes have electrophilic carbonyl carbon that readily undergoes nucleophilic addition reactions"
                },
                {
                    id: "c26",
                    text: "Central atom hybridization in ICl₄⁻ is:",
                    options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
                    answer: 2,
                    explain: "ICl₄⁻ has 6 electron pairs (4 bonding + 2 lone pairs) around I, requiring sp³d² hybridization"
                },
                {
                    id: "c27",
                    text: "Which is diamagnetic?",
                    options: ["O₂", "NO", "N₂", "NO₂"],
                    answer: 2,
                    explain: "N₂ has all electrons paired in molecular orbitals, making it diamagnetic"
                },
                {
                    id: "c28",
                    text: "Oxidation number of Mn in K₂MnO₄ is:",
                    options: ["+4", "+6", "+7", "+3"],
                    answer: 1,
                    explain: "In manganate ion MnO₄²⁻: Mn + 4(-2) = -2, solving gives Mn = +6"
                },
                {
                    id: "c29",
                    text: "Coordination isomerism is possible in:",
                    options: ["[Co(NH₃)₆][Cr(CN)₆]", "[Co(NH₃)₅Cl]Cl₂", "[Co(en)₃]Cl₃", "[CoCl₄]²⁻"],
                    answer: 0,
                    explain: "Complex salts with different complex ions can show coordination isomerism by exchanging ligands"
                },
                {
                    id: "c30",
                    text: "In spinel structure AB₂O₄, coordination numbers are:",
                    options: ["A: 4, B: 6", "A: 6, B: 4", "A: 8, B: 4", "A: 6, B: 8"],
                    answer: 0,
                    explain: "In normal spinel structure, A²⁺ ions occupy tetrahedral sites (CN = 4) and B³⁺ ions occupy octahedral sites (CN = 6)"
                },
                {
                    id: "c31",
                    text: "Which has lowest surface tension?",
                    options: ["Water", "Mercury", "Ethanol", "Benzene"],
                    answer: 2,
                    explain: "Ethanol has lower surface tension than water due to weaker hydrogen bonding and hydrocarbon chain"
                },
                {
                    id: "c32",
                    text: "For reaction 2A + B → C, if doubling [A] increases rate 4 times, order with respect to A is:",
                    options: ["1", "2", "3", "0"],
                    answer: 1,
                    explain: "Rate = k[A]ⁿ[B]ᵐ. If [A] doubles and rate increases 4 times: 2ⁿ = 4, so n = 2"
                },
                {
                    id: "c33",
                    text: "Which has pyramidal structure?",
                    options: ["BF₃", "NH₃", "CO₂", "BeF₂"],
                    answer: 1,
                    explain: "NH₃ has trigonal pyramidal structure due to 3 bonding pairs and 1 lone pair around nitrogen"
                },
                {
                    id: "c34",
                    text: "Which element has highest first ionization energy?",
                    options: ["He", "Ne", "Ar", "Kr"],
                    answer: 0,
                    explain: "Helium has highest first ionization energy (2372 kJ/mol) among all elements due to small size and complete shell"
                },
                {
                    id: "c35",
                    text: "Wurtz reaction is used to prepare:",
                    options: ["Alkenes", "Alkanes", "Alkynes", "Alcohols"],
                    answer: 1,
                    explain: "Wurtz reaction involves coupling of alkyl halides using sodium metal to form alkanes"
                },
                {
                    id: "c36",
                    text: "VSEPR theory is based on:",
                    options: ["Minimizing electron pair repulsion", "Maximizing bonding", "Hybridization", "Electronegativity"],
                    answer: 0,
                    explain: "VSEPR theory predicts molecular geometry by minimizing repulsion between electron pairs around central atom"
                },
                {
                    id: "c37",
                    text: "Which halogen has maximum electron affinity?",
                    options: ["F", "Cl", "Br", "I"],
                    answer: 1,
                    explain: "Chlorine has maximum electron affinity due to optimal size - not too small like F (repulsion) or too large like Br, I"
                },
                {
                    id: "c38",
                    text: "Square planar complexes show which type of isomerism?",
                    options: ["Optical only", "Geometrical only", "Both optical and geometrical", "Neither"],
                    answer: 1,
                    explain: "Square planar complexes [MA₂B₂] can show cis-trans geometrical isomerism but rarely optical isomerism"
                },
                {
                    id: "c39",
                    text: "Which metal has highest density?",
                    options: ["Gold", "Platinum", "Osmium", "Iridium"],
                    answer: 2,
                    explain: "Osmium has highest density (22.59 g/cm³) among all elements"
                },
                {
                    id: "c40",
                    text: "Fullerene C₆₀ has how many double bonds?",
                    options: ["30", "60", "90", "120"],
                    answer: 0,
                    explain: "C₆₀ has 30 double bonds in its soccer ball structure with alternating single and double bonds"
                },
                {
                    id: "c41",
                    text: "Zaitsev rule predicts:",
                    options: ["Most substituted alkene as major product", "Least substituted alkene as major product", "Equal mixture of products", "No elimination"],
                    answer: 0,
                    explain: "Zaitsev rule states that elimination reactions prefer to form most substituted (thermodynamically stable) alkene"
                },
                {
                    id: "c42",
                    text: "EAN of chromium in [Cr(CO)₆] is:",
                    options: ["36", "54", "18", "24"],
                    answer: 0,
                    explain: "Cr (24) + 6×2 from CO ligands = 24 + 12 = 36 electrons (krypton configuration)"
                },
                {
                    id: "c43",
                    text: "Which can act as both oxidizing and reducing agent?",
                    options: ["HNO₃", "H₂SO₄", "H₂O₂", "HCl"],
                    answer: 2,
                    explain: "H₂O₂ can act as oxidizing agent (reduced to H₂O) or reducing agent (oxidized to O₂) depending on conditions"
                },
                {
                    id: "c44",
                    text: "Cyclopropane has how many sigma bonds?",
                    options: ["6", "9", "12", "3"],
                    answer: 1,
                    explain: "C₃H₆ has 9 σ bonds: 6 C-H bonds + 3 C-C bonds"
                },
                {
                    id: "c45",
                    text: "Which has strongest metallic bonding?",
                    options: ["Na", "Mg", "Al", "Si"],
                    answer: 2,
                    explain: "Aluminum has strongest metallic bonding due to 3 valence electrons and small size providing strong electrostatic attraction"
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "Kranz anatomy is characterized by:",
                    options: ["Scattered vascular bundles", "Bundle sheath with thick walls", "Absence of mesophyll", "Parallel venation"],
                    answer: 1,
                    explain: "Kranz anatomy in C₄ plants features bundle sheath cells with thick walls and no intercellular spaces, forming a wreath-like arrangement around vascular bundles"
                },
                {
                    id: "b2",
                    text: "In orchid seeds, food is stored in:",
                    options: ["Endosperm", "Perisperm", "Cotyledons", "No food storage"],
                    answer: 3,
                    explain: "Orchid seeds are dust-like with minimal food reserves, depending on mycorrhizal association for germination and early growth"
                },
                {
                    id: "b3",
                    text: "Which hormone promotes xylem differentiation?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Ethylene"],
                    answer: 0,
                    explain: "Auxin promotes differentiation of procambium into xylem tissue and maintains cambial activity"
                },
                {
                    id: "b4",
                    text: "Phelloderm is produced by:",
                    options: ["Cork cambium inwardly", "Cork cambium outwardly", "Vascular cambium", "Apical meristem"],
                    answer: 0,
                    explain: "Phelloderm (secondary cortex) is produced by cork cambium (phellogen) toward the inner side"
                },
                {
                    id: "b5",
                    text: "Stomata in moss capsule are:",
                    options: ["Functional", "Non-functional", "Absent", "Present only in young stage"],
                    answer: 1,
                    explain: "Stomata in moss capsules are non-functional as they lack guard cells and cannot regulate opening/closing"
                },
                {
                    id: "b6",
                    text: "Pollen-pistil recognition involves:",
                    options: ["Chemical signals only", "Physical barriers only", "Both chemical and physical factors", "Genetic factors only"],
                    answer: 2,
                    explain: "Pollen-pistil interaction involves both chemical recognition (proteins, enzymes) and physical compatibility"
                },
                {
                    id: "b7",
                    text: "Which step requires most ATP in Calvin cycle?",
                    options: ["Carboxylation", "Reduction", "Regeneration", "All require equal ATP"],
                    answer: 1,
                    explain: "Reduction phase of Calvin cycle requires 6 NADPH and 6 ATP to convert 6 molecules of 3-PGA to G3P"
                },
                {
                    id: "b8",
                    text: "Hydathodes are modified:",
                    options: ["Stomata", "Trichomes", "Lenticels", "Emergences"],
                    answer: 0,
                    explain: "Hydathodes are modified stomata at leaf margins that remain permanently open for guttation"
                },
                {
                    id: "b9",
                    text: "Coenobium is characteristic of:",
                    options: ["Spirogyra", "Volvox", "Chlamydomonas", "Ulothrix"],
                    answer: 1,
                    explain: "Volvox forms coenobium, a hollow sphere of cells working as integrated unit with coordinated movement"
                },
                {
                    id: "b10",
                    text: "Agamospermy results in:",
                    options: ["Hybrid vigor", "Genetic uniformity", "Increased variation", "Sterile offspring"],
                    answer: 1,
                    explain: "Agamospermy (asexual seed formation) produces genetically uniform offspring identical to parent"
                },
                {
                    id: "b11",
                    text: "Which hormone was first discovered?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "ABA"],
                    answer: 0,
                    explain: "Auxin was first plant hormone discovered by Charles Darwin and later isolated by Went"
                },
                {
                    id: "b12",
                    text: "Photosystem II contains:",
                    options: ["P680 only", "P680 and water-splitting complex", "P700 only", "Both P680 and P700"],
                    answer: 1,
                    explain: "PSII contains P680 reaction center and oxygen-evolving complex that splits water molecules"
                },
                {
                    id: "b13",
                    text: "Long day plants flower when:",
                    options: ["Day length exceeds critical period", "Night length is below critical period", "Both conditions met", "Temperature is optimal"],
                    answer: 2,
                    explain: "Long day plants require day length above and night length below their respective critical periods"
                },
                {
                    id: "b14",
                    text: "Soil water potential is usually:",
                    options: ["Positive", "Negative", "Zero", "Variable depending on soil type"],
                    answer: 1,
                    explain: "Soil water potential is typically negative due to matric forces and solute effects binding water to soil particles"
                },
                {
                    id: "b15",
                    text: "P-protein plugs in sieve tubes are formed during:",
                    options: ["Normal transport", "Injury response", "Seasonal changes", "Aging"],
                    answer: 1,
                    explain: "P-protein plugs form at sieve plates as response to injury or stress, blocking phloem transport"
                },
                {
                    id: "b16",
                    text: "Double fertilization was first observed in:",
                    options: ["Lilium", "Fritillaria", "Petunia", "Tobacco"],
                    answer: 0,
                    explain: "Double fertilization was first described by Nawaschin in Lilium and Fritillaria"
                },
                {
                    id: "b17",
                    text: "Amphivasal vascular bundles are found in:",
                    options: ["Monocot stems", "Dicot stems", "Monocot roots", "Some fern rhizomes"],
                    answer: 3,
                    explain: "Amphivasal bundles (xylem surrounding phloem) occur in some ferns and few monocots like Dracaena"
                },
                {
                    id: "b18",
                    text: "Reaction wood in angiosperms is:",
                    options: ["Compression wood", "Tension wood", "Both types", "Neither type"],
                    answer: 1,
                    explain: "Angiosperms form tension wood on upper side of bent branches, rich in cellulose"
                },
                {
                    id: "b19",
                    text: "Phycocyanin pigment is found in:",
                    options: ["Green algae", "Brown algae", "Red algae", "Blue-green algae"],
                    answer: 3,
                    explain: "Phycocyanin is blue pigment characteristic of cyanobacteria (blue-green algae)"
                },
                {
                    id: "b20",
                    text: "Amphitropous ovule shows:",
                    options: ["No curvature", "Complete inversion", "Horseshoe-shaped curvature", "Partial inversion"],
                    answer: 2,
                    explain: "Amphitropous ovule shows horseshoe-shaped curvature with micropyle at right angles to funicle"
                },
                {
                    id: "b21",
                    text: "Phototropism is mediated by:",
                    options: ["Phytochrome", "Cryptochrome", "Phototropin", "All photoreceptors"],
                    answer: 2,
                    explain: "Phototropin is blue light receptor specifically responsible for phototropic responses in plants"
                },
                {
                    id: "b22",
                    text: "Pseudogamy occurs in:",
                    options: ["Normal sexual reproduction", "Apomixis requiring pollination", "Vegetative reproduction", "Fragmentation"],
                    answer: 1,
                    explain: "Pseudogamy is type of apomixis where pollination is required to stimulate development but no fertilization occurs"
                },
                {
                    id: "b23",
                    text: "Bundle sheath cells in C₄ plants:",
                    options: ["Lack chloroplasts", "Have agranal chloroplasts", "Have normal chloroplasts", "Are non-photosynthetic"],
                    answer: 1,
                    explain: "Bundle sheath cells in C₄ plants contain agranal chloroplasts specialized for Calvin cycle"
                },
                {
                    id: "b24",
                    text: "Samara fruit is found in:",
                    options: ["Maple", "Mango", "Apple", "Orange"],
                    answer: 0,
                    explain: "Samara is winged dry fruit found in maple, ash, and elm trees for wind dispersal"
                },
                {
                    id: "b25",
                    text: "Strasburger cells are:",
                    options: ["Companion cells", "Albuminous cells", "Sieve cells", "Phloem parenchyma"],
                    answer: 1,
                    explain: "Strasburger cells (albuminous cells) are associated with sieve cells in gymnosperms"
                },
                {
                    id: "b26",
                    text: "Red drop phenomenon indicates:",
                    options: ["Photosystem saturation", "Two light reactions", "Chlorophyll breakdown", "Photorespiration"],
                    answer: 1,
                    explain: "Red drop phenomenon demonstrated existence of two photosystems working in series"
                },
                {
                    id: "b27",
                    text: "Cambial activity is highest during:",
                    options: ["Spring", "Summer", "Autumn", "Winter"],
                    answer: 0,
                    explain: "Cambial activity peaks in spring when favorable conditions promote rapid cell division"
                },
                {
                    id: "b28",
                    text: "Ubisch bodies are found in:",
                    options: ["Microspores", "Megaspores", "Tapetum", "Pollen tube"],
                    answer: 2,
                    explain: "Ubisch bodies are small spherical structures in tapetum cells containing sporopollenin"
                },
                {
                    id: "b29",
                    text: "Leghemoglobin in root nodules:",
                    options: ["Transports oxygen", "Maintains low oxygen levels", "Stores nitrogen", "Fixes carbon"],
                    answer: 1,
                    explain: "Leghemoglobin maintains optimal low oxygen levels around nitrogenase enzyme while allowing controlled oxygen supply"
                },
                {
                    id: "b30",
                    text: "Stone cells are type of:",
                    options: ["Parenchyma", "Collenchyma", "Sclerenchyma", "Aerenchyma"],
                    answer: 2,
                    explain: "Stone cells (sclereids) are type of sclerenchyma with thick lignified walls providing mechanical support"
                },
                {
                    id: "b31",
                    text: "Thigmonasty is shown by:",
                    options: ["Sunflower", "Venus flytrap", "Pea tendrils", "Root tips"],
                    answer: 1,
                    explain: "Venus flytrap shows thigmonasty (touch-sensitive rapid movement) to capture prey"
                },
                {
                    id: "b32",
                    text: "Metaxylem vessels have:",
                    options: ["Annular thickening", "Spiral thickening", "Reticulate or pitted thickening", "No thickening"],
                    answer: 2,
                    explain: "Metaxylem vessels show reticulate (net-like) or pitted thickening patterns"
                },
                {
                    id: "b33",
                    text: "Abscission layer formation involves:",
                    options: ["Cell wall strengthening", "Cell wall weakening", "Cell division", "Cell enlargement"],
                    answer: 1,
                    explain: "Abscission involves enzymatic weakening of cell walls in separation layer, facilitating organ detachment"
                },
                {
                    id: "b34",
                    text: "Symplastic transport involves movement through:",
                    options: ["Cell walls", "Intercellular spaces", "Plasmodesmata", "Xylem vessels"],
                    answer: 2,
                    explain: "Symplastic pathway involves movement through plasmodesmata connecting cytoplasm of adjacent cells"
                },
                {
                    id: "b35",
                    text: "NADP-malic enzyme C₄ plants are adapted to:",
                    options: ["Cool climates", "Hot dry climates", "Shade conditions", "Waterlogged conditions"],
                    answer: 1,
                    explain: "NADP-ME type C₄ plants like maize are typically adapted to hot, dry conditions"
                },
                {
                    id: "b36",
                    text: "Root cap is absent in:",
                    options: ["Terrestrial roots", "Aerial roots", "Aquatic roots", "Parasitic roots"],
                    answer: 2,
                    explain: "Root cap is typically absent in aquatic roots as it's not needed for soil penetration"
                },
                {
                    id: "b37",
                    text: "Tetramerous flowers are characteristic of:",
                    options: ["Monocots", "Dicots", "Cruciferae", "Leguminosae"],
                    answer: 2,
                    explain: "Four-merous (tetramerous) flowers are characteristic feature of Cruciferae (Brassicaceae) family"
                },
                {
                    id: "b38",
                    text: "Somatic hybridization involves:",
                    options: ["Gamete fusion", "Protoplast fusion", "Pollen fusion", "Tissue fusion"],
                    answer: 1,
                    explain: "Somatic hybridization involves fusion of protoplasts from different plant species or varieties"
                },
                {
                    id: "b39",
                    text: "Pneumatophores help in:",
                    options: ["Water absorption", "Nutrient storage", "Gaseous exchange", "Mechanical support"],
                    answer: 2,
                    explain: "Pneumatophores are aerial roots in mangroves that facilitate gaseous exchange in waterlogged conditions"
                },
                {
                    id: "b40",
                    text: "Pressure bomb technique measures:",
                    options: ["Osmotic pressure", "Turgor pressure", "Water potential", "Root pressure"],
                    answer: 2,
                    explain: "Pressure bomb (Scholander pressure chamber) measures water potential of plant tissues"
                },
                {
                    id: "b41",
                    text: "Somaclonal variation is useful for:",
                    options: ["Maintaining genetic uniformity", "Creating genetic diversity", "Preventing mutations", "Tissue preservation"],
                    answer: 1,
                    explain: "Somaclonal variation generates genetic diversity useful in plant breeding and crop improvement"
                },
                {
                    id: "b42",
                    text: "Cytochrome b6f complex connects:",
                    options: ["PSII to PSI", "PSI to NADP⁺", "ATP synthase to PSI", "Water to PSII"],
                    answer: 0,
                    explain: "Cytochrome b6f complex transfers electrons from PSII to PSI in photosynthetic electron transport chain"
                },
                {
                    id: "b43",
                    text: "Composite fruit develops from:",
                    options: ["Single flower", "Multiple flowers", "Inflorescence", "Single carpel"],
                    answer: 2,
                    explain: "Composite (multiple) fruits develop from entire inflorescence with many flowers contributing to single fruit structure"
                },
                {
                    id: "b44",
                    text: "Tyloses are formed by:",
                    options: ["Vessel elements", "Tracheids", "Adjacent parenchyma cells", "Fibers"],
                    answer: 2,
                    explain: "Tyloses are balloon-like outgrowths from adjacent parenchyma cells that intrude into vessel lumens"
                },
                {
                    id: "b45",
                    text: "Synthetic cytokinin commonly used is:",
                    options: ["BAP", "Kinetin", "Zeatin", "Both BAP and Kinetin"],
                    answer: 3,
                    explain: "Both BAP (benzylaminopurine) and kinetin are widely used synthetic cytokinins in tissue culture"
                },
                
                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Jugular foramen is passage for:",
                    options: ["Internal carotid artery", "Internal jugular vein", "Vertebral artery", "External carotid artery"],
                    answer: 1,
                    explain: "Jugular foramen in skull base allows passage of internal jugular vein, glossopharyngeal, vagus, and accessory nerves"
                },
                {
                    id: "z2",
                    text: "Urogastrone is:",
                    options: ["Gastrin", "Secretin", "EGF", "CCK"],
                    answer: 2,
                    explain: "Urogastrone is alternate name for epidermal growth factor (EGF) that inhibits gastric acid secretion"
                },
                {
                    id: "z3",
                    text: "Protein C deficiency causes:",
                    options: ["Bleeding disorder", "Thrombotic tendency", "Hemolytic anemia", "Platelet dysfunction"],
                    answer: 1,
                    explain: "Protein C is natural anticoagulant; its deficiency increases risk of venous thrombosis"
                },
                {
                    id: "z4",
                    text: "Liddle syndrome involves:",
                    options: ["Increased sodium reabsorption", "Decreased sodium reabsorption", "Increased potassium retention", "Decreased aldosterone"],
                    answer: 0,
                    explain: "Liddle syndrome is caused by overactive epithelial sodium channels causing excessive sodium retention"
                },
                {
                    id: "z5",
                    text: "Mirror neurons are involved in:",
                    options: ["Motor control", "Vision", "Imitation and empathy", "Memory"],
                    answer: 2,
                    explain: "Mirror neurons fire both when performing action and observing others perform same action, important for imitation"
                },
                {
                    id: "z6",
                    text: "Gilbert syndrome affects:",
                    options: ["Bilirubin conjugation", "Heme synthesis", "Iron metabolism", "Albumin synthesis"],
                    answer: 0,
                    explain: "Gilbert syndrome involves reduced UDP-glucuronosyltransferase activity affecting bilirubin conjugation"
                },
                {
                    id: "z7",
                    text: "GIP (glucose-dependent insulinotropic peptide) is secreted by:",
                    options: ["Alpha cells", "Beta cells", "K cells", "L cells"],
                    answer: 2,
                    explain: "GIP is incretin hormone secreted by K cells in duodenum and jejunum in response to glucose"
                },
                {
                    id: "z8",
                    text: "Pulsus paradoxus is:",
                    options: ["Stronger pulse during inspiration", "Weaker pulse during inspiration", "Irregular pulse", "Absent pulse"],
                    answer: 1,
                    explain: "Pulsus paradoxus is exaggerated decrease in systolic pressure during inspiration (>10 mmHg)"
                },
                {
                    id: "z9",
                    text: "Apparent mineralocorticoid excess involves:",
                    options: ["Increased aldosterone", "Decreased aldosterone", "11β-HSD2 deficiency", "21-hydroxylase deficiency"],
                    answer: 2,
                    explain: "AME syndrome involves 11β-hydroxysteroid dehydrogenase type 2 deficiency allowing cortisol to act on aldosterone receptors"
                },
                {
                    id: "z10",
                    text: "Cortical granules in oocyte contain:",
                    options: ["Enzymes for zona hardening", "Nutrients", "mRNA", "Mitochondria"],
                    answer: 0,
                    explain: "Cortical granules release enzymes that modify zona pellucida to prevent polyspermy after fertilization"
                },
                {
                    id: "z11",
                    text: "Bitot's spots appear as:",
                    options: ["Yellow patches", "Foamy white patches", "Red spots", "Dark areas"],
                    answer: 1,
                    explain: "Bitot's spots are foamy, whitish patches on conjunctiva characteristic of vitamin A deficiency"
                },
                {
                    id: "z12",
                    text: "Pentraxins include:",
                    options: ["C-reactive protein", "Serum amyloid P", "Both CRP and SAP", "Immunoglobulins"],
                    answer: 2,
                    explain: "Pentraxins are family of proteins including C-reactive protein and serum amyloid P component"
                },
                {
                    id: "z13",
                    text: "Fanconi syndrome affects:",
                    options: ["Glomerular filtration", "Proximal tubular reabsorption", "Distal tubular secretion", "Collecting duct function"],
                    answer: 1,
                    explain: "Fanconi syndrome involves generalized proximal tubular dysfunction affecting reabsorption of multiple substances"
                },
                {
                    id: "z14",
                    text: "Sperm undergo capacitation for approximately:",
                    options: ["30 minutes", "2-3 hours", "6-8 hours", "12-24 hours"],
                    answer: 2,
                    explain: "Sperm capacitation in female reproductive tract takes approximately 6-8 hours to complete"
                },
                {
                    id: "z15",
                    text: "Cartwright blood group is associated with:",
                    options: ["AChE activity", "Complement activation", "Hemolysis", "Drug metabolism"],
                    answer: 0,
                    explain: "Cartwright blood group antigens are located on acetylcholinesterase enzyme"
                },
                {
                    id: "z16",
                    text: "Alveolar ventilation equals:",
                    options: ["Tidal volume × respiratory rate", "(Tidal volume - dead space) × respiratory rate", "Vital capacity", "Total lung capacity"],
                    answer: 1,
                    explain: "Alveolar ventilation = (tidal volume - anatomical dead space) × respiratory rate"
                },
                {
                    id: "z17",
                    text: "Saccule responds to:",
                    options: ["Angular acceleration", "Linear acceleration in horizontal plane", "Linear acceleration in vertical plane", "Sound waves"],
                    answer: 2,
                    explain: "Saccule contains otoliths and responds to linear acceleration in vertical plane (up-down motion)"
                },
                {
                    id: "z18",
                    text: "Normal filtration fraction is:",
                    options: ["10%", "20%", "30%", "40%"],
                    answer: 1,
                    explain: "Filtration fraction (GFR/RPF) is normally about 20%, meaning 20% of renal plasma flow is filtered"
                },
                {
                    id: "z19",
                    text: "McCune-Albright syndrome involves:",
                    options: ["GH deficiency", "Precocious puberty", "Growth retardation", "Delayed puberty"],
                    answer: 1,
                    explain: "McCune-Albright syndrome can cause precocious puberty due to autonomous hormone production"
                },
                {
                    id: "z20",
                    text: "Gastrulation occurs during:",
                    options: ["1st week", "2nd week", "3rd week", "4th week"],
                    answer: 2,
                    explain: "Gastrulation occurs during 3rd week of human embryonic development, forming three germ layers"
                },
                {
                    id: "z21",
                    text: "Windkessel effect refers to:",
                    options: ["Elastic recoil of arteries", "Venous return", "Cardiac contractility", "Blood viscosity"],
                    answer: 0,
                    explain: "Windkessel effect describes elastic recoil of large arteries helping maintain diastolic pressure"
                },
                {
                    id: "z22",
                    text: "Sulfhemoglobin has characteristic:",
                    options: ["Bright red color", "Cherry red color", "Chocolate brown color", "Greenish color"],
                    answer: 3,
                    explain: "Sulfhemoglobin appears greenish and cannot carry oxygen effectively"
                },
                {
                    id: "z23",
                    text: "Romberg's test assesses:",
                                        options: ["Motor function", "Sensory function", "Proprioception and balance", "Cognitive function"],
                    answer: 2,
                    explain: "Romberg's test assesses proprioception and balance by testing ability to maintain posture with eyes closed"
                },
                {
                    id: "z24",
                    text: "Functional residual capacity is sum of:",
                    options: ["TV + IRV", "ERV + RV", "TV + ERV", "IRV + ERV"],
                    answer: 1,
                    explain: "FRC = Expiratory reserve volume (ERV) + Residual volume (RV), lung volume at resting expiratory level"
                },
                {
                    id: "z25",
                    text: "White muscle fibers are:",
                    options: ["Type I", "Type IIa", "Type IIx", "Type III"],
                    answer: 2,
                    explain: "Type IIx (fast-twitch glycolytic) fibers appear white due to low myoglobin content and fewer mitochondria"
                },
                {
                    id: "z26",
                    text: "Kupffer cells are:",
                    options: ["Hepatocytes", "Macrophages", "Stellate cells", "Endothelial cells"],
                    answer: 1,
                    explain: "Kupffer cells are specialized macrophages located in liver sinusoids that phagocytose pathogens and debris"
                },
                {
                    id: "z27",
                    text: "Pellagra is characterized by:",
                    options: ["3 Ds", "4 Ds", "5 Ds", "2 Ds"],
                    answer: 1,
                    explain: "Pellagra (niacin deficiency) is characterized by 4 Ds: dermatitis, diarrhea, dementia, and death"
                },
                {
                    id: "z28",
                    text: "α-actinin is located at:",
                    options: ["Z-disc", "M-line", "A-band", "I-band"],
                    answer: 0,
                    explain: "α-actinin is major component of Z-disc that anchors thin filaments and maintains sarcomere structure"
                },
                {
                    id: "z29",
                    text: "Apelin has what effect on appetite?",
                    options: ["Increases", "Decreases", "No effect", "Variable effect"],
                    answer: 1,
                    explain: "Apelin is adipokine that generally decreases appetite and has cardiovascular protective effects"
                },
                {
                    id: "z30",
                    text: "Peyer's patches are found in:",
                    options: ["Duodenum", "Jejunum", "Ileum", "All parts of small intestine"],
                    answer: 2,
                    explain: "Peyer's patches are lymphoid aggregations predominantly located in ileum for immune surveillance"
                },
                {
                    id: "z31",
                    text: "Non-shivering thermogenesis is mediated by:",
                    options: ["Skeletal muscle", "Brown adipose tissue", "Liver", "Brain"],
                    answer: 1,
                    explain: "Brown adipose tissue generates heat through uncoupling protein-1 without muscle contraction"
                },
                {
                    id: "z32",
                    text: "Bloody show occurs due to:",
                    options: ["Cervical dilatation", "Membrane rupture", "Placental separation", "Uterine contraction"],
                    answer: 0,
                    explain: "Bloody show is mucus plug discharge with blood during cervical effacement and early dilatation"
                },
                {
                    id: "z33",
                    text: "Target cells are seen in:",
                    options: ["Thalassemia", "Sickle cell disease", "Liver disease", "All of these"],
                    answer: 3,
                    explain: "Target cells (codocytes) appear in various conditions including thalassemia, liver disease, and hemoglobinopathies"
                },
                {
                    id: "z34",
                    text: "Ultradian rhythms have period:",
                    options: ["Less than 20 hours", "Exactly 24 hours", "More than 28 hours", "Variable"],
                    answer: 0,
                    explain: "Ultradian rhythms have periods shorter than circadian rhythms (less than 20 hours), like REM sleep cycles"
                },
                {
                    id: "z35",
                    text: "Mydriasis is caused by:",
                    options: ["Parasympathetic stimulation", "Sympathetic stimulation", "Both", "Neither"],
                    answer: 1,
                    explain: "Mydriasis (pupil dilation) results from sympathetic stimulation or parasympathetic blockade"
                },
                {
                    id: "z36",
                    text: "Creatinine clearance overestimates GFR because:",
                    options: ["Tubular reabsorption", "Tubular secretion", "Protein binding", "Metabolism"],
                    answer: 1,
                    explain: "Creatinine undergoes slight tubular secretion, causing creatinine clearance to slightly overestimate true GFR"
                },
                {
                    id: "z37",
                    text: "Phantom limb sensation involves:",
                    options: ["Peripheral nerve damage", "Spinal cord plasticity", "Cortical reorganization", "All of these"],
                    answer: 3,
                    explain: "Phantom limb sensations result from complex interactions at peripheral, spinal, and cortical levels"
                },
                {
                    id: "z38",
                    text: "Postprandial glucose peak occurs at:",
                    options: ["30 minutes", "60-90 minutes", "2 hours", "3 hours"],
                    answer: 1,
                    explain: "Blood glucose typically peaks 60-90 minutes after meal intake in healthy individuals"
                },
                {
                    id: "z39",
                    text: "Wolffian ducts develop into:",
                    options: ["Female reproductive tract", "Male reproductive tract", "Urinary system", "Both male tract and urinary system"],
                    answer: 3,
                    explain: "Wolffian (mesonephric) ducts form male reproductive tract and contribute to kidney/ureter development"
                },
                {
                    id: "z40",
                    text: "Ehlers-Danlos syndrome affects:",
                    options: ["Collagen synthesis", "Elastin production", "Keratin formation", "Fibrin polymerization"],
                    answer: 0,
                    explain: "Ehlers-Danlos syndrome involves defects in collagen synthesis causing joint hypermobility and skin hyperextensibility"
                },
                {
                    id: "z41",
                    text: "Gastrinoma most commonly occurs in:",
                    options: ["Stomach", "Duodenum", "Pancreas", "Both duodenum and pancreas"],
                    answer: 3,
                    explain: "Gastrinomas occur most commonly in duodenum and pancreas (gastrinoma triangle)"
                },
                {
                    id: "z42",
                    text: "Fibrous joints are classified as:",
                    options: ["Synarthroses", "Amphiarthroses", "Diarthroses", "Variable mobility"],
                    answer: 0,
                    explain: "Fibrous joints are synarthroses (immovable joints) like sutures in skull"
                },
                {
                    id: "z43",
                    text: "Reverse T3 (rT3) is:",
                    options: ["More active than T3", "Inactive metabolite", "Thyroid hormone precursor", "Storage form"],
                    answer: 1,
                    explain: "Reverse T3 is biologically inactive metabolite formed by 5-deiodinase acting on T4"
                },
                {
                    id: "z44",
                    text: "Laplace's law for heart states that wall tension is proportional to:",
                    options: ["Pressure × radius", "Pressure/radius", "Radius/pressure", "Pressure² × radius"],
                    answer: 0,
                    explain: "Laplace's law: Wall tension = Pressure × radius / wall thickness, explaining why dilated hearts require more energy"
                },
                {
                    id: "z45",
                    text: "Fast block to polyspermy involves:",
                    options: ["Membrane depolarization", "Cortical granule exocytosis", "Zona pellucida hardening", "Calcium release"],
                    answer: 0,
                    explain: "Fast block involves immediate membrane depolarization preventing additional sperm binding within seconds of fertilization"
                }
            ]
        }
    ]
};

