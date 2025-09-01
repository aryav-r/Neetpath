v// mock-test-14.js - NEET Mock Test 14 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_14 = {
    id: "neet-014",
    title: "Full Syllabus Mock 14", 
    level: "hard",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A particle moves such that its position vector is r⃗ = 3t²î + 4t³ĵ. The magnitude of acceleration at t = 2s is:",
                    options: ["6 m/s²", "48 m/s²", "50 m/s²", "52 m/s²"],
                    answer: 3,
                    explain: "Acceleration a⃗ = d²r⃗/dt² = 6î + 24tĵ. At t = 2s: a⃗ = 6î + 48ĵ. Magnitude |a⃗| = √(6² + 48²) = √(36 + 2304) = √2340 ≈ 48.4 m/s²"
                },
                {
                    id: "p2",
                    text: "A solid cylinder and solid sphere of equal mass and radius start from rest and roll down identical inclines. The ratio of their kinetic energies at the bottom is:",
                    options: ["6:5", "5:6", "7:5", "5:7"],
                    answer: 2,
                    explain: "For cylinder: KE = ½mv²(1 + I/mR²) = ½mv²(1 + ½) = ¾mv². For sphere: KE = ½mv²(1 + ⅖) = 7/10 mv². Since mgh is same, v² ∝ 1/(1 + I/mR²). Ratio = (7/10):(3/4) × (7/5):(10/6) = 7:5"
                },
                {
                    id: "p3",
                    text: "A 10μF capacitor charged to 100V is connected to an uncharged 20μF capacitor. The energy lost during redistribution is:",
                    options: ["16.7 mJ", "33.3 mJ", "50 mJ", "25 mJ"],
                    answer: 0,
                    explain: "Initial energy = ½ × 10×10⁻⁶ × (100)² = 50 mJ. Final voltage = 100 × 10/(10+20) = 33.3V. Final energy = ½ × 30×10⁻⁶ × (33.3)² = 16.7 mJ. Energy lost = 50 - 33.3 = 16.7 mJ"
                },
                {
                    id: "p4",
                    text: "In Lloyd's mirror experiment, if distance between source and screen is 2 m and height of source above mirror is 1 mm, fringe width is 0.5 mm. Wavelength of light is:",
                    options: ["500 nm", "250 nm", "750 nm", "600 nm"],
                    answer: 0,
                    explain: "In Lloyd's mirror: d = 2h = 2 × 1 mm = 2 mm. Fringe width β = λD/d = 0.5 mm. λ = βd/D = 0.5×10⁻³ × 2×10⁻³/2 = 0.5×10⁻⁶ m = 500 nm"
                },
                {
                    id: "p5",
                    text: "In photoelectric effect, if frequency of incident light is doubled, stopping potential:",
                    options: ["Doubles", "More than doubles", "Less than doubles", "Becomes four times"],
                    answer: 1,
                    explain: "KEmax = hf - φ. If f becomes 2f: KE'max = 2hf - φ = 2(hf - φ) + φ = 2KEmax + φ. Since eV₀ = KEmax, stopping potential more than doubles"
                },
                {
                    id: "p6",
                    text: "A conducting rod of length 0.5 m rotates about one end in uniform magnetic field 0.2 T. If angular velocity is 4 rad/s, EMF between ends is:",
                    options: ["0.2 V", "0.4 V", "0.1 V", "0.8 V"],
                    answer: 2,
                    explain: "EMF for rotating rod = ½BωL² = ½ × 0.2 × 4 × (0.5)² = ½ × 0.2 × 4 × 0.25 = 0.1 V"
                },
                {
                    id: "p7",
                    text: "A uniform disc is rolling on horizontal surface with velocity v. It then rolls up an incline of angle θ. Maximum distance along incline is:",
                    options: ["v²/g sin θ", "3v²/4g sin θ", "3v²/2g sin θ", "v²/2g sin θ"],
                    answer: 1,
                    explain: "Initial KE = ½mv² + ½Iω² = ½mv² + ¼mv² = ¾mv². At maximum height: PE = mgh = mg(s sin θ). Energy conservation: ¾mv² = mgs sin θ. s = 3v²/4g sin θ"
                },
                {
                    id: "p8",
                    text: "In series LCR circuit, if L = 0.1 H, C = 10μF, R = 100Ω, resonant frequency is:",
                    options: ["159 Hz", "318 Hz", "50 Hz", "100 Hz"],
                    answer: 0,
                    explain: "Resonant frequency f₀ = 1/(2π√LC) = 1/(2π√(0.1 × 10×10⁻⁶)) = 1/(2π × 10⁻³) = 1000/2π ≈ 159 Hz"
                },
                {
                    id: "p9",
                    text: "A monoatomic ideal gas undergoes process PV² = constant. Molar heat capacity for this process is:",
                    options: ["3R/2", "0", "-3R/2", "R"],
                    answer: 2,
                    explain: "For PVⁿ = constant, C = Cv + R/(1-n). Here n = 2, Cv = 3R/2 for monoatomic gas. C = 3R/2 + R/(1-2) = 3R/2 - R = -3R/2"
                },
                {
                    id: "p10",
                    text: "Two long parallel wires separated by 10 cm carry currents 5A and 8A in same direction. Force per unit length between them is:",
                    options: ["8×10⁻⁵ N/m", "4×10⁻⁵ N/m", "2×10⁻⁵ N/m", "6×10⁻⁵ N/m"],
                    answer: 0,
                    explain: "Force per unit length = μ₀I₁I₂/2πd = (4π×10⁻⁷ × 5 × 8)/(2π × 0.1) = (2×10⁻⁶ × 40)/0.2 = 8×10⁻⁵ N/m"
                },
                {
                    id: "p11",
                    text: "A particle executes SHM with amplitude 4 cm and period 2 s. Speed when displacement is 2 cm from mean position is:",
                    options: ["2π√3 cm/s", "4π√3 cm/s", "π√3 cm/s", "6π√3 cm/s"],
                    answer: 0,
                    explain: "ω = 2π/T = π rad/s. v = ω√(A² - x²) = π√(16 - 4) = π√12 = 2π√3 cm/s"
                },
                {
                    id: "p12",
                    text: "In Compton scattering, for scattering angle 90°, fractional change in wavelength is:",
                    options: ["h/mₑc", "2h/mₑc", "h/2mₑc", "3h/2mₑc"],
                    answer: 0,
                    explain: "Δλ/λ = (h/mₑcλ)(1 - cos θ). For θ = 90°: Δλ/λ = (h/mₑcλ)(1 - 0) = h/mₑcλ. But since we need fractional change, answer is h/mₑc when λ is normalized"
                },
                {
                    id: "p13",
                    text: "A cube of side 'a' has charges +q at four corners of one face and -q at four corners of opposite face. Electric dipole moment is:",
                    options: ["qa", "4qa", "2qa", "qa/2"],
                    answer: 1,
                    explain: "Net charge +4q at one face center, -4q at opposite face. Separation = a. Dipole moment = (4q) × a = 4qa"
                },
                {
                    id: "p14",
                    text: "A particle projected at 45° has range 100 m on horizontal ground. Maximum height reached is:",
                    options: ["50 m", "25 m", "75 m", "100 m"],
                    answer: 1,
                    explain: "Range R = u² sin 2θ/g = u²/g (for θ = 45°). Maximum height H = u² sin² θ/2g = u²/4g (for θ = 45°). H = R/4 = 100/4 = 25 m"
                },
                {
                    id: "p15",
                    text: "A toroidal solenoid has 500 turns wound over a ring of mean radius 20 cm. For current 2 A, magnetic field inside is:",
                    options: ["10⁻³ T", "2×10⁻³ T", "5×10⁻⁴ T", "1.5×10⁻³ T"],
                    answer: 0,
                    explain: "B = μ₀nI/2πr for toroidal solenoid. B = (4π×10⁻⁷ × 500 × 2)/(2π × 0.2) = (2×10⁻⁴ × 2×10³)/0.4 = 10⁻³ T"
                },
                {
                    id: "p16",
                    text: "A ball dropped from 20 m height bounces to 15 m. Coefficient of restitution is:",
                    options: ["0.75", "0.87", "0.6", "0.5"],
                    answer: 1,
                    explain: "e = √(h₂/h₁) = √(15/20) = √(3/4) = √3/2 ≈ 0.87"
                },
                {
                    id: "p17",
                    text: "In parallel RLC circuit, at frequency below resonance:",
                    options: ["Circuit is capacitive", "Circuit is inductive", "Circuit is resistive", "Current is maximum"],
                    answer: 1,
                    explain: "Below resonance frequency, XL > XC in parallel circuit. Net reactance is inductive, so circuit behaves inductively"
                },
                {
                    id: "p18",
                    text: "Energy levels of hydrogen atom are En = -13.6/n² eV. Energy of Lyman α line is:",
                    options: ["10.2 eV", "3.4 eV", "1.9 eV", "12.1 eV"],
                    answer: 0,
                    explain: "Lyman α corresponds to n = 2 → n = 1 transition. Energy = E₂ - E₁ = -13.6/4 - (-13.6/1) = -3.4 + 13.6 = 10.2 eV"
                },
                {
                    id: "p19",
                    text: "Two identical springs each with spring constant k support a mass m. If connected in parallel, frequency of oscillation is:",
                    options: ["√(k/m)/2π", "√(2k/m)/2π", "√(k/2m)/2π", "√(4k/m)/2π"],
                    answer: 1,
                    explain: "Parallel springs: keff = k₁ + k₂ = 2k. Frequency f = (1/2π)√(keff/m) = (1/2π)√(2k/m)"
                },
                {
                    id: "p20",
                    text: "A point charge q is at distance d from center of earthed conducting sphere of radius R (d > R). Force on charge is:",
                    options: ["kq²R/d²(d²-R²)", "kq²R/4d²(d²-R²)", "kq²R/2d(d²-R²)", "kq²/4d²"],
                    answer: 0,
                    explain: "Using method of images, image charge = -qR/d at distance R²/d from center. Force = kq(-qR/d)/(d-R²/d)² = -kq²R/d²(d²-R²)"
                },
                {
                    id: "p21",
                    text: "A spherical mirror has focal length 20 cm. For object at 30 cm, image distance is:",
                    options: ["60 cm", "12 cm", "-60 cm", "-12 cm"],
                    answer: 0,
                    explain: "Using mirror equation: 1/f = 1/u + 1/v. 1/20 = 1/30 + 1/v. 1/v = 1/20 - 1/30 = (3-2)/60 = 1/60. v = 60 cm"
                },
                {
                    id: "p22",
                    text: "Two sound waves of frequencies 254 Hz and 256 Hz interfere. Beat frequency heard is:",
                    options: ["2 Hz", "255 Hz", "510 Hz", "1 Hz"],
                    answer: 0,
                    explain: "Beat frequency = |f₁ - f₂| = |254 - 256| = 2 Hz"
                },
                {
                    id: "p23",
                    text: "For reversible isothermal expansion of ideal gas, entropy change is:",
                    options: ["nR ln(V₂/V₁)", "nR ln(V₁/V₂)", "Zero", "nCv ln(T₂/T₁)"],
                    answer: 0,
                    explain: "For isothermal process: ΔS = nR ln(V₂/V₁) = nR ln(P₁/P₂)"
                },
                {
                    id: "p24",
                    text: "A radioactive sample has initial activity 8000 Bq. After 6 hours, activity is 1000 Bq. Half-life is:",
                    options: ["2 hours", "3 hours", "1.5 hours", "4 hours"],
                    answer: 0,
                    explain: "Activity reduces from 8000 to 1000 in 6 hours = reduction by factor of 8 = 2³. So 3 half-lives = 6 hours. Half-life = 2 hours"
                },
                {
                    id: "p25",
                    text: "Magnetic energy stored in solenoid with inductance L carrying current I is:",
                    options: ["LI²", "½LI²", "2LI²", "LI"],
                    answer: 1,
                    explain: "Magnetic energy stored = ½LI² where L is inductance and I is current"
                },
                {
                    id: "p26",
                    text: "Principal quantum number determines:",
                    options: ["Shape of orbital", "Orientation of orbital", "Energy of electron", "Spin of electron"],
                    answer: 2,
                    explain: "Principal quantum number n primarily determines the energy of electron in hydrogen-like atoms: En ∝ -1/n²"
                },
                {
                    id: "p27",
                    text: "A ray entering glass slab (μ = 1.5) at 60° emerges at angle:",
                    options: ["60°", "35.3°", "40°", "45°"],
                    answer: 0,
                    explain: "For parallel-sided glass slab, emergent ray is parallel to incident ray regardless of thickness. Emergent angle = incident angle = 60°"
                },
                {
                    id: "p28",
                    text: "In organ pipe closed at one end, third harmonic has frequency:",
                    options: ["5f₁", "3f₁", "7f₁", "4f₁"],
                    answer: 0,
                    explain: "In closed pipe, only odd harmonics exist: f₁, 3f₁, 5f₁, 7f₁... Third harmonic = 5f₁"
                },
                {
                    id: "p29",
                    text: "COP of heat pump is related to COP of refrigerator as:",
                    options: ["COPhp = COPref + 1", "COPhp = COPref - 1", "COPhp = 2COPref", "COPhp = COPref/2"],
                    answer: 0,
                    explain: "COPhp = Q₁/W, COPref = Q₂/W. Since Q₁ = Q₂ + W: COPhp = (Q₂ + W)/W = Q₂/W + 1 = COPref + 1"
                },
                {
                    id: "p30",
                    text: "Cyclotron frequency depends on:",
                    options: ["Velocity of particle", "Mass and charge of particle", "Radius of path", "Kinetic energy"],
                    answer: 1,
                    explain: "Cyclotron frequency f = qB/2πm depends only on charge q, magnetic field B, and mass m of particle"
                },
                {
                    id: "p31",
                    text: "When sound travels from air to water, which property remains unchanged?",
                    options: ["Wavelength", "Frequency", "Speed", "Intensity"],
                    answer: 1,
                    explain: "Frequency is characteristic of source and remains unchanged when wave enters different medium"
                },
                {
                    id: "p32",
                    text: "A thin rod of mass M and length L has moment of inertia about perpendicular axis through center:",
                    options: ["ML²/12", "ML²/3", "ML²/6", "ML²/4"],
                    answer: 0,
                    explain: "For uniform thin rod about perpendicular axis through center: I = ML²/12"
                },
                {
                    id: "p33",
                    text: "In AC circuit, power factor is maximum when:",
                    options: ["XL = XC", "XL > XC", "XL < XC", "R = 0"],
                    answer: 0,
                    explain: "Power factor = cos φ = R/Z is maximum when Z is minimum, which occurs at resonance when XL = XC"
                },
                {
                    id: "p34",
                    text: "A convex lens of focal length 15 cm forms real image at 30 cm. Object is at:",
                    options: ["30 cm", "10 cm", "45 cm", "7.5 cm"],
                    answer: 0,
                    explain: "Using lens equation: 1/f = 1/u + 1/v. 1/15 = 1/u + 1/30. 1/u = 1/15 - 1/30 = 1/30. u = 30 cm"
                },
                {
                    id: "p35",
                    text: "Mass of deuteron is 2.014 u. Its binding energy is approximately:",
                    options: ["2.2 MeV", "1.1 MeV", "4.4 MeV", "0.55 MeV"],
                    answer: 0,
                    explain: "Mass defect = (1.007 + 1.008 - 2.014) u = 0.001 u. Binding energy = 0.001 × 931.5 MeV ≈ 2.2 MeV"
                },
                {
                    id: "p36",
                    text: "In driven oscillations, resonance occurs when driving frequency equals:",
                    options: ["Natural frequency", "Twice natural frequency", "Half natural frequency", "Angular frequency"],
                    answer: 0,
                    explain: "Maximum amplitude in forced oscillations occurs when driving frequency matches natural frequency of the system"
                },
                {
                    id: "p37",
                    text: "Magnetic susceptibility of diamagnetic material is:",
                    options: ["Large positive", "Small positive", "Small negative", "Large negative"],
                    answer: 2,
                    explain: "Diamagnetic materials have small negative magnetic susceptibility (χ ≈ -10⁻⁵)"
                },
                {
                    id: "p38",
                    text: "In ideal transformer, if primary has 100 turns and secondary has 400 turns, voltage ratio is:",
                    options: ["1:4", "4:1", "1:2", "2:1"],
                    answer: 0,
                    explain: "Voltage ratio = turns ratio = Ns/Np = 400/100 = 4:1 or 1:4 (secondary:primary)"
                },
                {
                    id: "p39",
                    text: "Time constant of CR circuit (capacitor discharging through resistor) is:",
                    options: ["RC", "R/C", "C/R", "1/RC"],
                    answer: 0,
                    explain: "Time constant τ = RC for both charging and discharging of capacitor through resistor"
                },
                {
                    id: "p40",
                    text: "Work done in rotating dipole from θ = 0° to θ = 180° in uniform electric field E is:",
                    options: ["pE", "2pE", "-pE", "-2pE"],
                    answer: 1,
                    explain: "Work done = U(180°) - U(0°) = (-pE cos 180°) - (-pE cos 0°) = pE - (-pE) = 2pE"
                },
                {
                    id: "p41",
                    text: "In adiabatic process for ideal gas, PVᵞ = constant. For monoatomic gas, γ equals:",
                    options: ["7/5", "5/3", "4/3", "3/2"],
                    answer: 1,
                    explain: "γ = Cp/Cv. For monoatomic gas: Cv = 3R/2, Cp = 5R/2. γ = (5R/2)/(3R/2) = 5/3"
                },
                {
                    id: "p42",
                    text: "In pure capacitive circuit, current leads voltage by:",
                    options: ["90°", "45°", "180°", "0°"],
                    answer: 0,
                    explain: "In pure capacitive circuit, current leads voltage by 90° (π/2 radians)"
                },
                {
                    id: "p43",
                    text: "For destructive interference in thin films, path difference should be:",
                    options: ["nλ", "(2n+1)λ/2", "nλ/2", "2nλ"],
                    answer: 1,
                    explain: "For destructive interference with phase change: path difference = (2n+1)λ/2 where n = 0,1,2..."
                },
                {
                    id: "p44",
                    text: "Binding energy per nucleon is maximum for nuclei with mass number around:",
                    options: ["20", "56", "120", "240"],
                    answer: 1,
                    explain: "Iron-56 has maximum binding energy per nucleon (~8.8 MeV), making it most stable nucleus"
                },
                {
                    id: "p45",
                    text: "In p-n junction, depletion layer width depends on:",
                    options: ["Forward bias voltage", "Reverse bias voltage", "Temperature", "All of these"],
                    answer: 3,
                    explain: "Depletion layer width increases with reverse bias, decreases with forward bias, and varies with temperature due to thermal generation of carriers"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which has maximum lattice energy?",
                    options: ["MgO", "CaO", "SrO", "BaO"],
                    answer: 0,
                    explain: "MgO has smallest cation (Mg²⁺) and highest charge density, resulting in maximum lattice energy"
                },
                {
                    id: "c2",
                    text: "The shape of ClF₄⁻ ion is:",
                    options: ["Tetrahedral", "Square planar", "See-saw", "Square pyramidal"],
                    answer: 1,
                    explain: "ClF₄⁻ has 6 electron pairs (4 bonding + 2 lone pairs) around Cl, giving square planar geometry"
                },
                {
                    id: "c3",
                    text: "Which shows maximum SN1 reactivity?",
                    options: ["C₆H₅CH₂Cl", "(C₆H₅)₂CHCl", "(C₆H₅)₃CCl", "C₆H₅CH₂CH₂Cl"],
                    answer: 2,
                    explain: "Trityl chloride (C₆H₅)₃CCl forms most stable carbocation due to resonance with three benzene rings"
                },
                {
                    id: "c4",
                    text: "Oxidation state of S in H₂S₂O₈ is:",
                    options: ["+6", "+7", "+5", "+4"],
                    answer: 0,
                    explain: "In peroxydisulfuric acid H₂S₂O₈: 2(+1) + 2S + 8(-2) = 0, solving gives S = +6"
                },
                {
                    id: "c5",
                    text: "Maximum number of stereoisomers for a compound with 2 chiral centers is:",
                    options: ["2", "4", "6", "8"],
                    answer: 1,
                    explain: "Maximum stereoisomers = 2ⁿ where n is number of chiral centers. For n=2: 2² = 4"
                },
                {
                    id: "c6",
                    text: "Strongest Brønsted base among these is:",
                    options: ["CH₃COO⁻", "C₆H₅O⁻", "OH⁻", "NH₂⁻"],
                    answer: 3,
                    explain: "NH₂⁻ (amide ion) is strongest base as NH₃ is weakest acid among conjugate acids"
                },
                {
                    id: "c7",
                    text: "Which transition metal shows maximum number of oxidation states?",
                    options: ["Sc", "Mn", "Fe", "Zn"],
                    answer: 1,
                    explain: "Manganese shows oxidation states from -3 to +7, maximum among transition metals"
                },
                {
                    id: "c8",
                    text: "Which has maximum bond angle?",
                    options: ["NH₃", "PH₃", "AsH₃", "SbH₃"],
                    answer: 0,
                    explain: "Bond angle decreases down the group due to decreasing electronegativity and increasing size. NH₃ has maximum bond angle (~107°)"
                },
                {
                    id: "c9",
                    text: "Which violates both octet rule and VSEPR theory?",
                    options: ["BF₃", "PF₅", "SF₆", "None of these"],
                    answer: 3,
                    explain: "All given compounds follow VSEPR theory for geometry prediction. BF₃ violates octet rule but follows VSEPR"
                },
                {
                    id: "c10",
                    text: "ΔG° = 0 indicates:",
                    options: ["Reaction is spontaneous", "Reaction is non-spontaneous", "System is at equilibrium", "Reaction is irreversible"],
                    answer: 2,
                    explain: "ΔG° = 0 means system is at equilibrium with no net change in free energy"
                },
                {
                    id: "c11",
                    text: "Which is strongest reducing agent in aqueous solution?",
                    options: ["Mg", "Al", "Zn", "Fe"],
                    answer: 1,
                    explain: "Aluminum has most negative reduction potential in aqueous solution, making it strongest reducing agent among given options"
                },
                {
                    id: "c12",
                    text: "Which complex is diamagnetic?",
                    options: ["[Fe(CN)₆]⁴⁻", "[Fe(H₂O)₆]²⁺", "[Co(H₂O)₆]²⁺", "[Ni(H₂O)₆]²⁺"],
                    answer: 0,
                    explain: "[Fe(CN)₆]⁴⁻ has Fe²⁺ (d⁶) with strong field CN⁻ ligands forming low-spin complex with all electrons paired"
                },
                {
                    id: "c13",
                    text: "Which gives carbylamine test?",
                    options: ["Primary amines", "Secondary amines", "Tertiary amines", "Quaternary ammonium salts"],
                    answer: 0,
                    explain: "Only primary amines give positive carbylamine (isocyanide) test with CHCl₃ and alcoholic KOH"
                },
                {
                    id: "c14",
                    text: "Correct order of solubility in water is:",
                    options: ["LiF > NaF > KF", "KF > NaF > LiF", "NaF > LiF > KF", "NaF > KF > LiF"],
                    answer: 1,
                    explain: "Solubility increases down group for fluorides due to decreasing lattice energy. KF > NaF > LiF"
                },
                {
                    id: "c15",
                    text: "Which buffer has maximum buffer capacity at pH 4.8?",
                    options: ["CH₃COOH/CH₃COONa (pKa = 4.8)", "HCOOH/HCOONa (pKa = 3.8)", "NH₄Cl/NH₄OH (pKb = 4.8)", "H₃PO₄/H₂PO₄⁻ (pKa = 2.1)"],
                    answer: 0,
                    explain: "Buffer capacity is maximum when pH = pKa. CH₃COOH/CH₃COONa has pKa = 4.8"
                },
                {
                    id: "c16",
                    text: "Magnetic moment of [V(H₂O)₆]³⁺ is:",
                    options: ["2.83 BM", "3.87 BM", "4.90 BM", "1.73 BM"],
                    answer: 0,
                    explain: "V³⁺ (d²) has 2 unpaired electrons. μ = √[n(n+2)] = √[2×4] = 2.83 BM"
                },
                {
                    id: "c17",
                    text: "Which bond is most polar?",
                    options: ["C-N", "C-O", "C-F", "C-Cl"],
                    answer: 2,
                    explain: "C-F bond has maximum electronegativity difference (4.0 - 2.5 = 1.5), hence most polar"
                },
                {
                    id: "c18",
                    text: "For consecutive reaction A → B → C, rate determining step is:",
                    options: ["First step always", "Slowest step", "Second step always", "Fastest step"],
                    answer: 1,
                    explain: "Rate determining step is always the slowest step in consecutive reactions"
                },
                {
                    id: "c19",
                    text: "Which shows maximum boiling point?",
                    options: ["C₂H₅OH", "C₂H₅SH", "C₂H₅NH₂", "(CH₃)₂O"],
                    answer: 0,
                    explain: "Ethanol forms strongest hydrogen bonds due to OH group and high electronegativity of oxygen"
                },
                {
                    id: "c20",
                    text: "Which is aromatic according to Hückel's rule?",
                    options: ["Cyclobutadiene", "[8]Annulene", "Tropylium ion", "All of these"],
                    answer: 2,
                    explain: "Tropylium ion (C₇H₇⁺) has 6π electrons (4n+2, n=1) and is aromatic"
                },
                {
                    id: "c21",
                    text: "Which has zero dipole moment due to symmetry?",
                    options: ["trans-1,2-dichloroethene", "cis-1,2-dichloroethene", "1,1-dichloroethene", "Chloroethene"],
                    answer: 0,
                    explain: "trans-1,2-dichloroethene has symmetrical structure with C-Cl dipoles canceling each other"
                },
                {
                    id: "c22",
                    text: "Which is strongest field ligand?",
                    options: ["H₂O", "NH₃", "CN⁻", "CO"],
                    answer: 3,
                    explain: "CO is strongest field ligand due to both σ-donation and π-back bonding capability"
                },
                {
                    id: "c23",
                    text: "pKa order for carboxylic acids is:",
                    options: ["CCl₃COOH > CHCl₂COOH > CH₂ClCOOH > CH₃COOH", "CH₃COOH > CH₂ClCOOH > CHCl₂COOH > CCl₃COOH", "Same for all", "Depends on temperature"],
                    answer: 1,
                    explain: "More electron-withdrawing groups decrease pKa (increase acidity). Order: CH₃COOH > CH₂ClCOOH > CHCl₂COOH > CCl₃COOH"
                },
                {
                    id: "c24",
                    text: "Bond order of NO⁺ is:",
                    options: ["2", "2.5", "3", "1.5"],
                    answer: 2,
                    explain: "NO⁺ has 14 electrons. Bond order = (bonding - antibonding)/2 = (10-4)/2 = 3"
                },
                {
                    id: "c25",
                    text: "Which shows electrophilic addition reaction?",
                    options: ["Benzene", "Ethene", "Ethyne", "Both ethene and ethyne"],
                    answer: 3,
                    explain: "Both alkenes (ethene) and alkynes (ethyne) undergo electrophilic addition reactions"
                },
                {
                    id: "c26",
                    text: "Number of hybrid orbitals in SF₄ is:",
                    options: ["4", "5", "6", "3"],
                    answer: 1,
                    explain: "SF₄ has 5 electron pairs (4 bonding + 1 lone pair), requiring 5 hybrid orbitals (sp³d)"
                },
                {
                    id: "c27",
                    text: "Which is ferromagnetic?",
                    options: ["O₂", "CrO₂", "MnO", "Both O₂ and CrO₂"],
                    answer: 1,
                    explain: "CrO₂ is ferromagnetic with parallel alignment of magnetic moments. O₂ is paramagnetic, not ferromagnetic"
                },
                {
                    id: "c28",
                    text: "Oxidation number of Cr in K₂Cr₂O₇ is:",
                    options: ["+6", "+7", "+12", "+3"],
                    answer: 0,
                    explain: "In dichromate K₂Cr₂O₇: 2(+1) + 2(Cr) + 7(-2) = 0, solving gives Cr = +6"
                },
                {
                    id: "c29",
                    text: "Which shows optical isomerism?",
                    options: ["[Pt(NH₃)₂Cl₂]", "[Co(en)₂Cl₂]⁺", "[Ni(NH₃)₆]²⁺", "[Cu(NH₃)₄]²⁺"],
                    answer: 1,
                    explain: "[Co(en)₂Cl₂]⁺ has asymmetric arrangement leading to non-superimposable mirror images (optical isomerism)"
                },
                {
                    id: "c30",
                    text: "In zinc blende structure, coordination number is:",
                    options: ["4:4", "6:6", "8:8", "6:3"],
                    answer: 0,
                    explain: "In zinc blende (ZnS) structure, both Zn²⁺ and S²⁻ have tetrahedral coordination (4:4)"
                },
                {
                    id: "c31",
                    text: "Which has lowest viscosity?",
                    options: ["Glycerol", "Ethylene glycol", "Ethanol", "Water"],
                    answer: 3,
                    explain: "Water has weakest intermolecular forces and least complex structure, hence lowest viscosity"
                },
                {
                    id: "c32",
                    text: "For pseudo first order reaction, rate constant has units:",
                    options: ["s⁻¹", "mol L⁻¹ s⁻¹", "mol⁻¹ L s⁻¹", "mol⁻² L² s⁻¹"],
                    answer: 0,
                    explain: "Pseudo first order reactions follow first order kinetics with rate constant units of s⁻¹"
                },
                {
                    id: "c33",
                    text: "Which has maximum bond angle?",
                    options: ["SO₂", "H₂O", "CO₂", "H₂S"],
                    answer: 2,
                    explain: "CO₂ has linear geometry with 180° bond angle, maximum among given options"
                },
                {
                    id: "c34",
                    text: "Ionization energy trend in a period is:",
                    options: ["Increases left to right", "Decreases left to right", "Remains constant", "First increases then decreases"],
                    answer: 0,
                    explain: "Ionization energy generally increases across period due to increasing effective nuclear charge"
                },
                {
                    id: "c35",
                    text: "Michael addition involves:",
                    options: ["α,β-unsaturated carbonyl compounds", "Simple alkenes", "Aromatic compounds", "Alkynes"],
                    answer: 0,
                    explain: "Michael addition is 1,4-addition of nucleophiles to α,β-unsaturated carbonyl compounds"
                },
                {
                    id: "c36",
                    text: "Which factor does not affect VSEPR geometry?",
                    options: ["Number of electron pairs", "Electronegativity of central atom", "Size of ligands", "Lone pair-bond pair repulsion"],
                    answer: 1,
                    explain: "VSEPR geometry depends on electron pair arrangement and repulsions, not on electronegativity of central atom"
                },
                {
                    id: "c37",
                    text: "Which group shows maximum +I effect?",
                    options: ["-CH₃", "-C₂H₅", "-(CH₃)₃", "-CH₂CH₃"],
                    answer: 2,
                    explain: "Tertiary butyl group -(CH₃)₃C shows maximum +I effect due to maximum electron donation"
                },
                {
                    id: "c38",
                    text: "Which can show both geometrical and optical isomerism?",
                    options: ["[Co(NH₃)₄Cl₂]⁺", "[Co(en)₂Cl₂]⁺", "[Pt(NH₃)₂Cl₂]", "[Co(NH₃)₆]³⁺"],
                    answer: 1,
                    explain: "[Co(en)₂Cl₂]⁺ shows cis-trans isomerism (geometrical) and each form shows optical isomerism"
                },
                {
                    id: "c39",
                    text: "Standard electrode potential is highest for:",
                    options: ["Au³⁺/Au", "Ag⁺/Ag", "Cu²⁺/Cu", "Zn²⁺/Zn"],
                    answer: 0,
                    explain: "Au³⁺/Au has highest standard reduction potential (+1.50 V), making Au least reactive"
                },
                {
                    id: "c40",
                    text: "Which carbon allotrope is hardest?",
                    options: ["Graphite", "Diamond", "Fullerene", "Graphene"],
                    answer: 1,
                    explain: "Diamond has 3D network of covalent bonds making it hardest carbon allotrope"
                },
                {
                    id: "c41",
                    text: "Markovnikov's rule is violated by:",
                    options: ["HBr + alkene in presence of peroxides", "HCl + alkene", "H₂SO₄ + alkene", "HI + alkene"],
                    answer: 0,
                    explain: "Anti-Markovnikov addition occurs when HBr adds to alkenes in presence of peroxides (free radical mechanism)"
                },
                {
                    id: "c42",
                    text: "EAN of central metal in [Fe(CO)₅] is:",
                    options: ["36", "54", "18", "86"],
                    answer: 0,
                    explain: "Fe (26) + 5×2 from CO ligands = 26 + 10 = 36 electrons (krypton configuration)"
                },
                {
                    id: "c43",
                    text: "Which can act as both Lewis acid and base?",
                    options: ["BF₃", "AlCl₃", "H₂O", "NH₃"],
                    answer: 2,
                    explain: "H₂O can donate lone pair (Lewis base) or accept lone pair via H⁺ (Lewis acid) - amphoteric behavior"
                },
                {
                    id: "c44",
                    text: "Benzene has how many pi bonds?",
                    options: ["3", "6", "9", "12"],
                    answer: 0,
                    explain: "Benzene has 3 delocalized π bonds forming aromatic sextet with 6 π electrons"
                },
                {
                    id: "c45",
                    text: "Which intermolecular force is strongest?",
                    options: ["van der Waals", "Dipole-dipole", "Hydrogen bonding", "Ion-dipole"],
                    answer: 3,
                    explain: "Ion-dipole interactions are strongest intermolecular forces, followed by hydrogen bonding"
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "In C₃ plants, RuBisCO shows maximum affinity for:",
                    options: ["CO₂ only", "O₂ only", "CO₂ but can bind O₂", "Neither CO₂ nor O₂"],
                    answer: 2,
                    explain: "RuBisCO has higher affinity for CO₂ but can also bind O₂, leading to photorespiration in C₃ plants"
                },
                {
                    id: "b2",
                    text: "Polyembryony is most commonly found in:",
                    options: ["Wheat", "Rice", "Citrus", "Maize"],
                    answer: 2,
                    explain: "Citrus fruits commonly show polyembryony with multiple embryos developing from nucellar cells"
                },
                {
                    id: "b3",
                    text: "Which hormone is responsible for parthenocarpy in tomatoes?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Both auxin and gibberellin"],
                    answer: 3,
                    explain: "Both auxins and gibberellins can induce parthenocarpy (fruit development without fertilization) in tomatoes"
                },
                {
                    id: "b4",
                    text: "Secondary phloem develops from:",
                    options: ["Cork cambium", "Vascular cambium", "Interfascicular cambium", "Fascicular cambium"],
                    answer: 1,
                    explain: "Secondary phloem is produced by vascular cambium toward the outer side during secondary growth"
                },
                {
                    id: "b5",
                    text: "Which structure is unique to bryophytes?",
                    options: ["Rhizoids", "Archegonia", "Antheridia", "Protonema"],
                    answer: 3,
                    explain: "Protonema is filamentous juvenile stage unique to bryophytes, developing from spore germination"
                },
                {
                    id: "b6",
                    text: "Chasmogamy refers to:",
                    options: ["Closed flower pollination", "Open flower pollination", "Self-incompatibility", "Cross-incompatibility"],
                    answer: 1,
                    explain: "Chasmogamy is pollination in open flowers where anthers and stigma are exposed"
                },
                {
                    id: "b7",
                    text: "Emerson effect demonstrates:",
                    options: ["Single light reaction", "Two light reactions working together", "Dark reaction", "Photorespiration"],
                    answer: 1,
                    explain: "Emerson enhancement effect shows that PSI and PSII work more efficiently together than separately"
                },
                {
                    id: "b8",
                    text: "Velamen tissue is found in:",
                    options: ["Hydrophytes", "Xerophytes", "Epiphytes", "Mesophytes"],
                    answer: 2,
                    explain: "Velamen is spongy tissue covering aerial roots of epiphytic orchids for water absorption from atmosphere"
                },
                {
                    id: "b9",
                    text: "Hypanthodium inflorescence is characteristic of:",
                    options: ["Ficus", "Helianthus", "Brassica", "Legumes"],
                    answer: 0,
                    explain: "Hypanthodium (specialized receptacle with flowers inside) is characteristic of Ficus species"
                },
                {
                    id: "b10",
                    text: "Adventive embryony is a form of:",
                    options: ["Sexual reproduction", "Apomixis", "Vegetative propagation", "Fragmentation"],
                    answer: 1,
                    explain: "Adventive embryony is type of apomixis where embryos develop from somatic cells of nucellus or integument"
                },
                {
                    id: "b11",
                    text: "α-amylase synthesis during seed germination is induced by:",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "ABA"],
                    answer: 1,
                    explain: "Gibberellic acid induces α-amylase synthesis in aleurone layer for starch breakdown during germination"
                },
                {
                    id: "b12",
                    text: "ATP synthase complex is located in:",
                    options: ["Thylakoid membrane", "Stroma", "Outer membrane", "Intermembrane space"],
                    answer: 0,
                    explain: "ATP synthase (CF₀-CF₁ complex) spans thylakoid membrane with CF₁ head in stroma"
                },
                {
                    id: "b13",
                    text: "Florigen is now identified as:",
                    options: ["GA₃", "Zeatin", "FT protein", "Abscisic acid"],
                    answer: 2,
                    explain: "Florigen, the hypothetical flowering hormone, is now identified as FT (Flowering Locus T) protein"
                },
                {
                    id: "b14",
                    text: "Pressure potential in plant cells is:",
                    options: ["Always negative", "Always positive", "Can be positive or negative", "Always zero"],
                    answer: 2,
                    explain: "Pressure potential can be positive (turgor pressure in turgid cells) or negative (tension in xylem)"
                },
                {
                    id: "b15",
                    text: "Callose is deposited in:",
                    options: ["Sieve plates", "Vessel walls", "Tracheid walls", "Companion cell walls"],
                    answer: 0,
                    explain: "Callose is β-1,3-glucan deposited in sieve plates to regulate phloem transport"
                },
                {
                    id: "b16",
                    text: "Perisperm develops from:",
                    options: ["Nucellus", "Integuments", "Endosperm", "Embryo sac"],
                    answer: 0,
                    explain: "Perisperm is nutritive tissue that develops from nucellus and persists in some seeds"
                },
                {
                    id: "b17",
                    text: "Which venation pattern is found in Calophyllum?",
                    options: ["Pinnate", "Palmate", "Parallel", "Reticulate with parallel secondaries"],
                    answer: 3,
                    explain: "Calophyllum shows reticulate venation with closely parallel secondary veins"
                },
                {
                    id: "b18",
                    text: "Libriform fibers are found in:",
                    options: ["Primary xylem", "Secondary xylem", "Phloem", "Cortex"],
                    answer: 1,
                    explain: "Libriform fibers are thick-walled, elongated cells providing mechanical support in secondary xylem"
                },
                {
                    id: "b19",
                    text: "Fucoxanthin is found in:",
                    options: ["Green algae", "Brown algae", "Red algae", "Blue-green algae"],
                    answer: 1,
                    explain: "Fucoxanthin is brown carotenoid pigment characteristic of brown algae (Phaeophyceae)"
                },
                {
                    id: "b20",
                    text: "Apotropous ovule orientation means:",
                    options: ["Micropyle towards hilum", "Micropyle away from hilum", "Curved embryo sac", "Straight ovule"],
                    answer: 1,
                    explain: "In apotropous ovule, micropyle is directed away from hilum, opposite to anatropous condition"
                },
                {
                    id: "b21",
                    text: "Shade avoidance response involves:",
                    options: ["Blue light", "Red light", "Far-red light", "Red/far-red light ratio"],
                    answer: 3,
                    explain: "Shade avoidance involves detection of low red/far-red light ratio by phytochrome system"
                },
                {
                    id: "b22",
                    text: "Cleistogamous flowers ensure:",
                    options: ["Cross-pollination", "Self-pollination", "No pollination", "Wind pollination"],
                    answer: 1,
                    explain: "Cleistogamous flowers remain closed and ensure automatic self-pollination without opening"
                },
                {
                    id: "b23",
                    text: "Glycolate is produced during:",
                    options: ["Calvin cycle", "Photorespiration", "Respiration", "Glycolysis"],
                    answer: 1,
                    explain: "Glycolate is 2-carbon compound produced when RuBisCO acts as oxygenase in photorespiration"
                },
                {
                    id: "b24",
                    text: "Which type of fruit is developed from inferior ovary?",
                    options: ["Drupe", "Berry", "Pome", "Capsule"],
                    answer: 2,
                    explain: "Pome fruits (apple, pear) develop from inferior ovary with fleshy thalamus surrounding true fruit"
                },
                {
                    id: "b25",
                    text: "Companion cells are absent in:",
                    options: ["Angiosperms", "Gymnosperms", "Pteridophytes", "Both gymnosperms and pteridophytes"],
                    answer: 3,
                    explain: "Companion cells are present only in angiosperms. Gymnosperms and pteridophytes have albuminous cells instead"
                },
                {
                    id: "b26",
                    text: "Oxygen-evolving complex contains:",
                    options: ["Manganese", "Iron", "Copper", "Zinc"],
                    answer: 0,
                    explain: "Oxygen-evolving complex in PSII contains manganese cluster that catalyzes water splitting"
                },
                {
                    id: "b27",
                    text: "Which meristem shows most rapid cell division?",
                    options: ["Apical meristem", "Intercalary meristem", "Lateral meristem", "All equal"],
                    answer: 1,
                    explain: "Intercalary meristem at base of internodes shows most rapid cell division, especially in grasses"
                },
                {
                    id: "b28",
                    text: "Endothecium is part of:",
                    options: ["Anther wall", "Ovary wall", "Seed coat", "Pollen grain"],
                    answer: 0,
                    explain: "Endothecium is inner layer of anther wall with fibrous thickenings that help in anther dehiscence"
                },
                {
                    id: "b29",
                    text: "Which bacteria can fix nitrogen in non-leguminous plants?",
                    options: ["Rhizobium", "Frankia", "Azotobacter", "Clostridium"],
                    answer: 1,
                    explain: "Frankia (actinomycete) forms symbiotic nitrogen-fixing nodules with non-legumes like Alnus, Casuarina"
                },
                {
                    id: "b30",
                    text: "Which tissue provides tensile strength to plant organs?",
                    options: ["Parenchyma", "Collenchyma", "Sclerenchyma", "Aerenchyma"],
                    answer: 2,
                    explain: "Sclerenchyma fibers provide tensile strength due to thick, lignified walls and elongated structure"
                },
                {
                    id: "b31",
                    text: "Circumnutation is:",
                    options: ["Response to gravity", "Spiral growth pattern", "Response to light", "Response to touch"],
                    answer: 1,
                    explain: "Circumnutation is inherent spiral or helical growth movement of plant organs during growth"
                },
                {
                    id: "b32",
                    text: "Which type of secondary wall thickening allows maximum stretching?",
                    options: ["Annular", "Spiral", "Reticulate", "Pitted"],
                    answer: 1,
                    explain: "Spiral thickening in protoxylem allows maximum stretching during organ elongation"
                },
                {
                    id: "b33",
                    text: "1-aminocyclopropane-1-carboxylic acid is precursor of:",
                    options: ["Auxin", "Cytokinin", "Ethylene", "Gibberellin"],
                    answer: 2,
                    explain: "ACC (1-aminocyclopropane-1-carboxylic acid) is immediate precursor of ethylene biosynthesis"
                },
                {
                    id: "b34",
                    text: "Active transport in plants requires:",
                    options: ["ATP only", "Carrier proteins only", "Both ATP and carrier proteins", "Neither ATP nor carriers"],
                    answer: 2,
                    explain: "Active transport requires both ATP for energy and specific carrier proteins for selective transport"
                },
                {
                    id: "b35",
                    text: "Stomatal opening in CAM plants occurs during:",
                    options: ["Day time", "Night time", "Dawn", "Dusk"],
                    answer: 1,
                    explain: "CAM plants open stomata at night to minimize water loss while fixing CO₂ as malic acid"
                },
                {
                    id: "b36",
                    text: "Casparian strips prevent:",
                    options: ["Apoplastic movement", "Symplastic movement", "Both movements", "Neither movement"],
                    answer: 0,
                    explain: "Casparian strips force water and solutes to move through symplastic pathway, preventing apoplastic bypass"
                },
                {
                    id: "b37",
                    text: "Isobilateral leaves have:",
                    options: ["Palisade on upper side only", "Palisade on both sides", "Spongy parenchyma only", "No differentiation"],
                    answer: 1,
                    explain: "Isobilateral leaves have palisade tissue on both upper and lower sides with stomata on both surfaces"
                },
                {
                    id: "b38",
                    text: "Mutation breeding involves:",
                    options: ["Cross-pollination", "Induced mutations", "Selection", "Hybridization"],
                    answer: 1,
                    explain: "Mutation breeding uses physical or chemical mutagens to induce mutations for crop improvement"
                },
                {
                    id: "b39",
                    text: "Buttress roots are found in:",
                    options: ["Desert plants", "Aquatic plants", "Tropical rainforest trees", "Alpine plants"],
                    answer: 2,
                    explain: "Buttress roots provide mechanical support to tall tropical rainforest trees in shallow soil"
                },
                {
                    id: "b40",
                    text: "Which process is driven by transpiration?",
                    options: ["Translocation", "Water absorption", "Photosynthesis", "Respiration"],
                    answer: 1,
                    explain: "Transpiration creates negative pressure that drives water absorption and transport in xylem"
                },
                {
                    id: "b41",
                    text: "Somatic embryogenesis involves:",
                    options: ["Fertilization", "Cell division only", "Dedifferentiation and redifferentiation", "Meiosis"],
                    answer: 2,
                    explain: "Somatic embryogenesis involves dedifferentiation of somatic cells followed by redifferentiation into embryos"
                },
                {
                    id: "b42",
                    text: "Which is the mobile electron carrier in photosynthesis?",
                    options: ["Cytochrome f", "Plastocyanin", "Ferredoxin", "All of these"],
                    answer: 3,
                    explain: "Plastocyanin, ferredoxin, and cytochrome f all act as mobile electron carriers in photosynthetic electron transport"
                },
                {
                    id: "b43",
                    text: "False fruit develops from:",
                    options: ["Ovary only", "Thalamus/receptacle", "Stamens", "Sepals"],
                    answer: 1,
                    explain: "False fruits (pseudocarps) develop from floral parts other than ovary, mainly thalamus/receptacle"
                },
                {
                    id: "b44",
                    text: "Spring wood has:",
                    options: ["Smaller vessels", "Larger vessels", "No vessels", "Same size vessels as autumn wood"],
                    answer: 1,
                    explain: "Spring wood has larger vessels and thinner walls to meet high water demand during active growth"
                },
                {
                    id: "b45",
                    text: "TIBA (tri-iodobenzoic acid) is:",
                    options: ["Auxin", "Anti-auxin", "Gibberellin", "Anti-gibberellin"],
                    answer: 1,
                    explain: "TIBA is auxin transport inhibitor that blocks polar auxin transport in plants"
                },
                
                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Maxilla forms part of:",
                    options: ["Hard palate", "Nasal septum", "Orbit", "All of these"],
                    answer: 3,
                    explain: "Maxilla contributes to hard palate, lateral nasal wall, floor of orbit, and forms part of nasal septum"
                },
                {
                    id: "z2",
                    text: "Enterocrinin stimulates:",
                    options: ["Gastric secretion", "Pancreatic secretion", "Intestinal secretion", "Bile secretion"],
                    answer: 2,
                    explain: "Enterocrinin stimulates intestinal juice secretion from crypts of Lieberkühn"
                },
                {
                    id: "z3",
                    text: "Hageman factor is:",
                    options: ["Factor X", "Factor XI", "Factor XII", "Factor XIII"],
                    answer: 2,
                    explain: "Hageman factor is Factor XII that initiates intrinsic pathway of coagulation cascade"
                },
                {
                    id: "z4",
                    text: "Renin is secreted by:",
                    options: ["Macula densa", "Juxtaglomerular cells", "Mesangial cells", "Podocytes"],
                    answer: 1,
                    explain: "Renin is secreted by juxtaglomerular cells in response to decreased blood pressure or sodium"
                },
                {
                    id: "z5",
                    text: "Declarative memory is processed in:",
                    options: ["Cerebellum", "Hippocampus", "Amygdala", "Thalamus"],
                    answer: 1,
                    explain: "Hippocampus is crucial for declarative (explicit) memory formation and consolidation"
                },
                {
                    id: "z6",
                    text: "Delta bilirubin is:",
                    options: ["Unconjugated", "Conjugated", "Covalently bound to albumin", "Free in plasma"],
                    answer: 2,
                    explain: "Delta bilirubin is conjugated bilirubin covalently bound to albumin, has longer half-life"
                },
                {
                    id: "z7",
                    text: "Villikinin stimulates:",
                    options: ["Villus contraction", "Villus elongation", "Crypt secretion", "Peristalsis"],
                    answer: 0,
                    explain: "Villikinin causes rhythmic contractions of intestinal villi to aid absorption"
                },
                {
                    id: "z8",
                    text: "Physiological splitting of S2 occurs during:",
                    options: ["Inspiration", "Expiration", "Both", "Neither"],
                    answer: 0,
                    explain: "Physiological splitting of S2 into A2-P2 occurs during inspiration due to delayed closure of pulmonary valve"
                },
                {
                    id: "z9",
                    text: "Aldosterone escape phenomenon prevents:",
                    options: ["Hyponatremia", "Hypernatremia", "Edema formation", "Dehydration"],
                    answer: 2,
                    explain: "Aldosterone escape prevents excessive sodium retention and edema despite continued aldosterone presence"
                },
                {
                    id: "z10",
                    text: "Cortical reaction is triggered by:",
                    options: ["Sperm binding", "Sperm penetration", "Calcium influx", "All of these"],
                    answer: 2,
                    explain: "Cortical reaction is triggered by calcium influx following sperm penetration, causing zona hardening"
                },
                {
                    id: "z11",
                    text: "Kwashiorkor is primarily due to:",
                    options: ["Protein deficiency", "Calorie deficiency", "Both protein and calorie deficiency", "Vitamin deficiency"],
                    answer: 0,
                    explain: "Kwashiorkor is caused primarily by protein deficiency with relatively adequate caloric intake"
                },
                {
                    id: "z12",
                    text: "Langerhans cells are:",
                    options: ["Pancreatic islet cells", "Epidermal dendritic cells", "Hepatic cells", "Neural cells"],
                    answer: 1,
                    explain: "Langerhans cells are dendritic cells in epidermis that present antigens and initiate immune responses"
                },
                {
                    id: "z13",
                    text: "Which hormone increases sodium excretion?",
                    options: ["Aldosterone", "ANF", "ADH", "Cortisol"],
                    answer: 1,
                    explain: "Atrial natriuretic factor (ANF) promotes sodium excretion and diuresis"
                },
                {
                    id: "z14",
                    text: "Leydig cells secrete:",
                    options: ["FSH", "LH", "Testosterone", "Inhibin"],
                    answer: 2,
                    explain: "Leydig cells (interstitial cells) in testes secrete testosterone in response to LH"
                },
                {
                    id: "z15",
                    text: "Lewis blood group antigens are found on:",
                    options: ["RBCs only", "Plasma proteins", "Secretions", "Both RBCs and secretions"],
                    answer: 3,
                    explain: "Lewis antigens are present on RBC surface and in body secretions like saliva and plasma"
                },
                {
                    id: "z16",
                    text: "Which lung capacity includes all four lung volumes?",
                    options: ["Vital capacity", "Functional residual capacity", "Total lung capacity", "Inspiratory capacity"],
                    answer: 2,
                    explain: "Total lung capacity = TV + IRV + ERV + RV, includes all four lung volumes"
                },
                {
                    id: "z17",
                    text: "Tectorial membrane is part of:",
                    options: ["Vestibular system", "Cochlear system", "Middle ear", "External ear"],
                                        answer: 1,
                    explain: "Tectorial membrane overlies hair cells in organ of Corti within cochlea for hearing mechanism"
                },
                {
                    id: "z18",
                    text: "Normal albumin/globulin ratio in plasma is:",
                    options: ["1:1", "1.5:1", "2:1", "3:1"],
                    answer: 1,
                    explain: "Normal A/G ratio is approximately 1.5:1 with albumin 3.5-5g/dL and globulin 2-3g/dL"
                },
                {
                    id: "z19",
                    text: "ACTH deficiency causes:",
                    options: ["Primary adrenal insufficiency", "Secondary adrenal insufficiency", "Tertiary adrenal insufficiency", "No adrenal effects"],
                    answer: 1,
                    explain: "ACTH deficiency from pituitary causes secondary adrenal insufficiency with low cortisol but normal aldosterone"
                },
                {
                    id: "z20",
                    text: "Sperm maturation occurs in:",
                    options: ["Seminiferous tubules", "Epididymis", "Vas deferens", "Seminal vesicles"],
                    answer: 1,
                    explain: "Sperm undergo final maturation and gain motility during transit through epididymis"
                },
                {
                    id: "z21",
                    text: "Afterload is increased by:",
                    options: ["Vasodilation", "Vasoconstriction", "Decreased blood volume", "Increased venous return"],
                    answer: 1,
                    explain: "Vasoconstriction increases arterial pressure, increasing afterload (resistance against which heart pumps)"
                },
                {
                    id: "z22",
                    text: "Oxygen debt occurs due to:",
                    options: ["Excess oxygen consumption", "Anaerobic metabolism", "Hyperventilation", "Carbon dioxide retention"],
                    answer: 1,
                    explain: "Oxygen debt results from anaerobic metabolism during intense exercise, requiring extra O₂ for recovery"
                },
                {
                    id: "z23",
                    text: "Positive Babinski sign indicates:",
                    options: ["Normal reflex", "Upper motor neuron lesion", "Lower motor neuron lesion", "Cerebellar lesion"],
                    answer: 1,
                    explain: "Positive Babinski sign (dorsiflexion of big toe) indicates upper motor neuron lesion"
                },
                {
                    id: "z24",
                    text: "Expiratory reserve volume is approximately:",
                    options: ["500 mL", "1000 mL", "1200 mL", "1500 mL"],
                    answer: 2,
                    explain: "ERV (maximum air that can be exhaled after normal expiration) is approximately 1200 mL"
                },
                {
                    id: "z25",
                    text: "Red muscle fibers are rich in:",
                    options: ["Glycogen", "Myoglobin", "Creatine phosphate", "Lactate"],
                    answer: 1,
                    explain: "Red muscle fibers (Type I) are rich in myoglobin and mitochondria for aerobic metabolism"
                },
                {
                    id: "z26",
                    text: "Central vein in liver lobule drains into:",
                    options: ["Portal vein", "Hepatic artery", "Hepatic vein", "Bile duct"],
                    answer: 2,
                    explain: "Central vein of hepatic lobule drains deoxygenated blood into hepatic vein, then to IVC"
                },
                {
                    id: "z27",
                    text: "Scurvy is caused by deficiency of:",
                    options: ["Vitamin B₁", "Vitamin B₁₂", "Vitamin C", "Vitamin D"],
                    answer: 2,
                    explain: "Scurvy results from vitamin C (ascorbic acid) deficiency affecting collagen synthesis"
                },
                {
                    id: "z28",
                    text: "Myomesin is located at:",
                    options: ["Z-disc", "M-line", "A-band", "I-band"],
                    answer: 1,
                    explain: "Myomesin is structural protein located at M-line that holds thick filaments in register"
                },
                {
                    id: "z29",
                    text: "Adiponectin has what effect on insulin sensitivity?",
                    options: ["Decreases", "Increases", "No effect", "Variable effect"],
                    answer: 1,
                    explain: "Adiponectin enhances insulin sensitivity and has anti-inflammatory effects"
                },
                {
                    id: "z30",
                    text: "Goblet cells secrete:",
                    options: ["Enzymes", "Hormones", "Mucus", "Acid"],
                    answer: 2,
                    explain: "Goblet cells in intestinal epithelium secrete mucus for lubrication and protection"
                },
                {
                    id: "z31",
                    text: "Which tissue has highest metabolic rate?",
                    options: ["Muscle", "Brain", "Liver", "Heart"],
                    answer: 1,
                    explain: "Brain tissue has highest metabolic rate, consuming ~20% of total body oxygen at rest"
                },
                {
                    id: "z32",
                    text: "Amniocentesis is usually performed at:",
                    options: ["8-12 weeks", "15-18 weeks", "24-28 weeks", "32-36 weeks"],
                    answer: 1,
                    explain: "Amniocentesis for genetic testing is typically performed at 15-18 weeks of pregnancy"
                },
                {
                    id: "z33",
                    text: "Heinz bodies are seen in:",
                    options: ["Normal RBCs", "G6PD deficiency", "Iron deficiency", "B12 deficiency"],
                    answer: 1,
                    explain: "Heinz bodies (denatured hemoglobin) appear in RBCs during oxidative stress, especially in G6PD deficiency"
                },
                {
                    id: "z34",
                    text: "Sleep spindles are generated by:",
                    options: ["Cortex", "Thalamus", "Hypothalamus", "Brainstem"],
                    answer: 1,
                    explain: "Sleep spindles are generated by thalamic reticular nucleus during stage 2 NREM sleep"
                },
                {
                    id: "z35",
                    text: "Miosis is caused by:",
                    options: ["Sympathetic stimulation", "Parasympathetic stimulation", "Both", "Neither"],
                    answer: 1,
                    explain: "Miosis (pupil constriction) is caused by parasympathetic stimulation via oculomotor nerve"
                },
                {
                    id: "z36",
                    text: "Juxtaglomerular apparatus regulates:",
                    options: ["Blood pressure", "Blood glucose", "Blood calcium", "Blood pH"],
                    answer: 0,
                    explain: "JGA secretes renin which regulates blood pressure through renin-angiotensin-aldosterone system"
                },
                {
                    id: "z37",
                    text: "Adaptation in sensory receptors involves:",
                    options: ["Increased sensitivity", "Decreased sensitivity", "Constant sensitivity", "Variable sensitivity"],
                    answer: 1,
                    explain: "Sensory adaptation involves decreased sensitivity to constant stimuli over time"
                },
                {
                    id: "z38",
                    text: "Gluconeogenesis is stimulated by:",
                    options: ["Insulin", "Glucagon", "IGF-1", "Leptin"],
                    answer: 1,
                    explain: "Glucagon stimulates gluconeogenesis (glucose synthesis from non-carbohydrate sources) in liver"
                },
                {
                    id: "z39",
                    text: "FSH receptor is present on:",
                    options: ["Leydig cells", "Sertoli cells", "Germ cells", "Interstitial cells"],
                    answer: 1,
                    explain: "FSH receptors are present on Sertoli cells in males and granulosa cells in females"
                },
                {
                    id: "z40",
                    text: "Pectus excavatum is:",
                    options: ["Protruding chest", "Sunken chest", "Barrel chest", "Normal chest"],
                    answer: 1,
                    explain: "Pectus excavatum (funnel chest) is congenital deformity with sunken/depressed sternum"
                },
                {
                    id: "z41",
                    text: "Chief cells secrete:",
                    options: ["HCl", "Pepsinogen", "Mucus", "Gastrin"],
                    answer: 1,
                    explain: "Chief cells in gastric glands secrete pepsinogen, the inactive precursor of pepsin"
                },
                {
                    id: "z42",
                    text: "Menisci in knee joint are made of:",
                    options: ["Hyaline cartilage", "Fibrocartilage", "Elastic cartilage", "Bone"],
                    answer: 1,
                    explain: "Menisci are C-shaped fibrocartilaginous structures that cushion and stabilize knee joint"
                },
                {
                    id: "z43",
                    text: "T4 to T3 conversion occurs mainly in:",
                    options: ["Thyroid", "Liver", "Kidney", "All peripheral tissues"],
                    answer: 3,
                    explain: "Peripheral conversion of T4 to active T3 occurs in liver, kidney, and other tissues via deiodinases"
                },
                {
                    id: "z44",
                    text: "Starling's law of heart states that:",
                    options: ["Force of contraction ∝ initial length", "Heart rate ∝ stroke volume", "Cardiac output ∝ blood pressure", "Contractility ∝ afterload"],
                    answer: 0,
                    explain: "Starling's law states that force of cardiac contraction is proportional to initial muscle fiber length (preload)"
                },
                {
                    id: "z45",
                    text: "Capacitation occurs in:",
                    options: ["Epididymis", "Vas deferens", "Female reproductive tract", "Seminiferous tubules"],
                    answer: 2,
                    explain: "Sperm capacitation (biochemical changes enabling fertilization) occurs in female reproductive tract"
                }
            ]
        }
    ]
};

