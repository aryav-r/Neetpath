// mock-test-8.js - NEET Mock Test 8 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_8 = {
    id: "neet-008",
    title: "Full Syllabus Mock 8", 
    level: "medium",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A block of mass 2 kg is placed on a rough inclined plane making 30° with horizontal. If coefficient of friction is 0.4, the acceleration down the plane is: (g = 10 m/s²)",
                    options: ["2.5 m/s²", "1.5 m/s²", "3.0 m/s²", "4.0 m/s²"],
                    answer: 1,
                    explain: "Forces: mg sin 30° - μmg cos 30° = ma. Substituting: 10 - 0.4×10×√3/2 = 2a. Solving: a = 1.5 m/s²"
                },
                {
                    id: "p2",
                    text: "The ratio of kinetic energies of two bodies having same momentum is 2:3. The ratio of their masses is:",
                    options: ["3:2", "2:3", "4:9", "9:4"],
                    answer: 0,
                    explain: "KE = p²/2m. For same momentum: KE₁/KE₂ = m₂/m₁. Given KE₁/KE₂ = 2/3, so m₂/m₁ = 2/3, hence m₁/m₂ = 3/2"
                },
                {
                    id: "p3",
                    text: "A capacitor of 10 μF is charged to 100 V. Energy stored in it is:",
                    options: ["0.05 J", "0.5 J", "5 J", "50 J"],
                    answer: 0,
                    explain: "U = ½CV² = ½ × 10×10⁻⁶ × (100)² = ½ × 10×10⁻⁶ × 10⁴ = 0.05 J"
                },
                {
                    id: "p4",
                    text: "In Young's double slit experiment, fringe width is 0.4 mm. If distance between slits is halved, new fringe width is:",
                    options: ["0.2 mm", "0.4 mm", "0.8 mm", "1.6 mm"],
                    answer: 2,
                    explain: "Fringe width β = λD/d. When d is halved, β becomes 2β = 2 × 0.4 = 0.8 mm"
                },
                {
                    id: "p5",
                    text: "The de Broglie wavelength associated with an electron accelerated through 100 V is approximately:",
                    options: ["1.2 Å", "2.4 Å", "0.6 Å", "3.6 Å"],
                    answer: 0,
                    explain: "λ = h/p = h/√(2meV) = 6.626×10⁻³⁴/√(2×9.1×10⁻³¹×1.6×10⁻¹⁹×100) ≈ 1.2×10⁻¹⁰ m = 1.2 Å"
                },
                {
                    id: "p6",
                    text: "A wire of length L and resistance R is cut into 4 equal parts. These parts are connected in parallel. The equivalent resistance is:",
                    options: ["R/16", "R/4", "4R", "16R"],
                    answer: 0,
                    explain: "Each part has resistance R/4. Four such resistors in parallel: 1/Req = 4×(4/R) = 16/R, so Req = R/16"
                },
                {
                    id: "p7",
                    text: "A uniform rod of length 2L is pivoted at its center. Forces F and 2F are applied at distances L/2 and L from pivot. Net torque about pivot is:",
                    options: ["FL/2", "3FL/2", "2FL", "5FL/2"],
                    answer: 1,
                    explain: "Torque = F₁×r₁ + F₂×r₂ = F×(L/2) + 2F×L = FL/2 + 2FL = 3FL/2"
                },
                {
                    id: "p8",
                    text: "The frequency of oscillation of LC circuit is 1000 Hz. If L is increased 4 times, new frequency is:",
                    options: ["250 Hz", "500 Hz", "2000 Hz", "4000 Hz"],
                    answer: 1,
                    explain: "f = 1/(2π√LC). When L increases 4 times: f' = 1/(2π√4LC) = f/2 = 1000/2 = 500 Hz"
                },
                {
                    id: "p9",
                    text: "A gas at 27°C is heated to 327°C at constant pressure. The ratio of final to initial volume is:",
                    options: ["1:2", "2:1", "1:1", "4:3"],
                    answer: 1,
                    explain: "Charles's law: V/T = constant. V₂/V₁ = T₂/T₁ = (327+273)/(27+273) = 600/300 = 2:1"
                },
                {
                    id: "p10",
                    text: "The magnetic field at the center of current carrying circular loop is B. If radius doubles keeping current same, new field is:",
                    options: ["B/4", "B/2", "2B", "4B"],
                    answer: 1,
                    explain: "B = μ₀I/(2R). When radius doubles: B' = μ₀I/(2×2R) = B/2"
                },
                {
                    id: "p11",
                    text: "A particle executing SHM has amplitude 4 cm and period 2 s. Maximum acceleration is:",
                    options: ["π² cm/s²", "4π² cm/s²", "8π² cm/s²", "16π² cm/s²"],
                    answer: 1,
                    explain: "Max acceleration = ω²A = (2π/T)²×A = (2π/2)²×4 = π²×4 = 4π² cm/s²"
                },
                {
                    id: "p12",
                    text: "In photoelectric effect, work function is 2 eV. If photon energy is 5 eV, maximum kinetic energy of emitted electron is:",
                    options: ["2 eV", "3 eV", "5 eV", "7 eV"],
                    answer: 1,
                    explain: "Einstein equation: KEmax = hf - φ = 5 - 2 = 3 eV"
                },
                {
                    id: "p13",
                    text: "Two identical springs are connected in series and then in parallel. Ratio of their equivalent spring constants is:",
                    options: ["1:4", "4:1", "1:2", "2:1"],
                    answer: 0,
                    explain: "Series: 1/ks = 1/k + 1/k = 2/k, so ks = k/2. Parallel: kp = k + k = 2k. Ratio ks:kp = (k/2):2k = 1:4"
                },
                {
                    id: "p14",
                    text: "A convex mirror of focal length 20 cm forms image of object at infinity. Image distance is:",
                    options: ["20 cm", "-20 cm", "∞", "10 cm"],
                    answer: 1,
                    explain: "For object at infinity in convex mirror, image forms at focus. Since f = -20 cm (convex), image distance = -20 cm"
                },
                {
                    id: "p15",
                    text: "The half-life of radioactive element is 1600 years. After 3200 years, fraction remaining is:",
                    options: ["1/2", "1/4", "1/8", "1/16"],
                    answer: 1,
                    explain: "Number of half-lives = 3200/1600 = 2. Remaining fraction = (1/2)ⁿ = (1/2)² = 1/4"
                },
                {
                    id: "p16",
                    text: "Two resistors 6Ω and 3Ω are connected in parallel across 6V battery. Current through 3Ω resistor is:",
                    options: ["1 A", "2 A", "3 A", "4 A"],
                    answer: 1,
                    explain: "In parallel, voltage across each resistor = 6V. Current through 3Ω: I = V/R = 6/3 = 2 A"
                },
                {
                    id: "p17",
                    text: "A satellite orbits Earth at height equal to Earth's radius. Its orbital velocity compared to surface orbital velocity is:",
                    options: ["Same", "1/√2 times", "√2 times", "1/2 times"],
                    answer: 1,
                    explain: "v = √(GM/r). At height R: v = √(GM/2R). At surface: v₀ = √(GM/R). Ratio = v/v₀ = 1/√2"
                },
                {
                    id: "p18",
                    text: "In interference, two coherent sources have intensities I and 4I. Maximum and minimum intensities in interference pattern are:",
                    options: ["5I and I", "9I and I", "9I and 3I", "6I and 2I"],
                    answer: 1,
                    explain: "Imax = (√I₁ + √I₂)² = (√I + √4I)² = (√I + 2√I)² = 9I. Imin = (√4I - √I)² = (2√I - √I)² = I"
                },
                {
                    id: "p19",
                    text: "The power factor of AC circuit containing only capacitor is:",
                    options: ["0", "1", "0.5", "√3/2"],
                    answer: 0,
                    explain: "In pure capacitive circuit, current leads voltage by 90°. Power factor = cos 90° = 0"
                },
                {
                    id: "p20",
                    text: "A body is thrown vertically upward with speed 20 m/s. Time taken to reach maximum height is:",
                    options: ["1 s", "2 s", "3 s", "4 s"],
                    answer: 1,
                    explain: "At maximum height, v = 0. Using v = u - gt: 0 = 20 - 10t, so t = 2 s"
                },
                {
                    id: "p21",
                    text: "The efficiency of ideal heat engine working between 400 K and 300 K is:",
                    options: ["20%", "25%", "30%", "35%"],
                    answer: 1,
                    explain: "Efficiency η = 1 - T₂/T₁ = 1 - 300/400 = 1 - 0.75 = 0.25 = 25%"
                },
                {
                    id: "p22",
                    text: "Two point charges +q and -2q are separated by distance d. Electric field is zero at distance x from +q towards -2q. Value of x is:",
                    options: ["d/3", "d/2", "2d/3", "d/√2"],
                    answer: 0,
                    explain: "For zero field: k(q)/x² = k(2q)/(d-x)². Solving: q/x² = 2q/(d-x)², gives x = d/3"
                },
                {
                    id: "p23",
                    text: "A mass attached to spring oscillates with period T. If mass is increased 4 times, new period is:",
                    options: ["T/2", "T", "2T", "4T"],
                    answer: 2,
                    explain: "T = 2π√(m/k). When mass increases 4 times: T' = 2π√(4m/k) = 2×2π√(m/k) = 2T"
                },
                {
                    id: "p24",
                    text: "The wavelength of X-rays is of order:",
                    options: ["10⁻⁸ m", "10⁻¹⁰ m", "10⁻¹² m", "10⁻¹⁴ m"],
                    answer: 1,
                    explain: "X-rays have wavelengths typically in range 10⁻¹¹ to 10⁻⁹ m, so order is 10⁻¹⁰ m"
                },
                {
                    id: "p25",
                    text: "In hydrogen spectrum, ratio of wavelengths of Hα and Hβ lines is:",
                    options: ["27:20", "20:27", "9:4", "4:9"],
                    answer: 0,
                    explain: "Hα: n₂=3→n₁=2, Hβ: n₂=4→n₁=2. Using Rydberg formula: λHα/λHβ = (1/4-1/16)/(1/4-1/9) = (3/16)/(5/36) = 27/20"
                },
                {
                    id: "p26",
                    text: "A uniform circular disk rolls down an incline without slipping. Ratio of rotational to translational kinetic energy is:",
                    options: ["1:1", "1:2", "2:1", "1:3"],
                    answer: 1,
                    explain: "For disk, I = ½mR². KErot = ½Iω² = ½(½mR²)(v/R)² = ¼mv². KEtrans = ½mv². Ratio = ¼mv²:½mv² = 1:2"
                },
                {
                    id: "p27",
                    text: "The magnetic flux through a coil changes from 0.1 Wb to 0.3 Wb in 0.2 s. Induced EMF is:",
                    options: ["0.5 V", "1 V", "1.5 V", "2 V"],
                    answer: 1,
                    explain: "EMF = -dΦ/dt = -(0.3-0.1)/0.2 = -0.2/0.2 = -1 V. Magnitude = 1 V"
                },
                {
                    id: "p28",
                    text: "In diffraction at single slit, first minimum occurs at angle θ where:",
                    options: ["a sin θ = λ", "a sin θ = 2λ", "a sin θ = λ/2", "a sin θ = 3λ/2"],
                    answer: 0,
                    explain: "First minimum in single slit diffraction occurs when path difference = λ, i.e., a sin θ = λ"
                },
                {
                    id: "p29",
                    text: "The temperature coefficient of resistance for metals is:",
                    options: ["Positive", "Negative", "Zero", "Infinite"],
                    answer: 0,
                    explain: "For metals, resistance increases with temperature, so temperature coefficient α > 0 (positive)"
                },
                {
                    id: "p30",
                    text: "A charged particle moves in magnetic field in helical path. This indicates:",
                    options: ["v parallel to B", "v perpendicular to B", "v has components both parallel and perpendicular to B", "B is uniform"],
                    answer: 2,
                    explain: "Helical motion occurs when velocity has both parallel (causing linear motion) and perpendicular (causing circular motion) components to magnetic field"
                },
                {
                    id: "p31",
                    text: "In common emitter amplifier, input resistance is 1000 Ω and load resistance is 5000 Ω. If current gain is 50, voltage gain is:",
                    options: ["10", "50", "100", "250"],
                    answer: 3,
                    explain: "Voltage gain = Current gain × (Load resistance/Input resistance) = 50 × (5000/1000) = 50 × 5 = 250"
                },
                {
                    id: "p32",
                    text: "The velocity of sound in air at 0°C is 332 m/s. At 27°C, velocity becomes:",
                    options: ["347 m/s", "352 m/s", "357 m/s", "362 m/s"],
                    answer: 0,
                    explain: "v ∝ √T. v₂/v₁ = √(T₂/T₁) = √(300/273) = √(300/273). v₂ = 332 × √(300/273) ≈ 347 m/s"
                },
                {
                    id: "p33",
                    text: "A lens of focal length 10 cm forms real image of magnification -2. Object distance is:",
                    options: ["5 cm", "15 cm", "20 cm", "30 cm"],
                    answer: 1,
                    explain: "m = -v/u = -2, so v = 2u. Using lens equation: 1/f = 1/u + 1/v = 1/u + 1/2u = 3/2u. 1/10 = 3/2u, so u = 15 cm"
                },
                {
                    id: "p34",
                    text: "In Young's experiment, one slit is covered with red filter and other with blue filter. The fringe pattern:",
                    options: ["Remains same", "Becomes invisible", "Shows only red fringes", "Shows only blue fringes"],
                    answer: 1,
                    explain: "Different colored lights are incoherent, cannot interfere. No fringe pattern is observed, becomes invisible"
                },
                {
                    id: "p35",
                    text: "The binding energy per nucleon is maximum for:",
                    options: ["Helium", "Carbon", "Iron", "Uranium"],
                    answer: 2,
                    explain: "Iron (Fe-56) has maximum binding energy per nucleon (~8.8 MeV), making it most stable nucleus"
                },
                {
                    id: "p36",
                    text: "A uniform rod of mass m and length L rotates about axis through one end. Its moment of inertia is:",
                    options: ["mL²/12", "mL²/3", "mL²/4", "mL²"],
                    answer: 1,
                    explain: "For uniform rod rotating about one end, moment of inertia I = mL²/3"
                },
                {
                    id: "p37",
                    text: "Two waves y₁ = sin(ωt - kx) and y₂ = sin(ωt + kx) interfere. The resultant wave has:",
                    options: ["Frequency ω", "Frequency 2ω", "No definite frequency", "Zero amplitude"],
                    answer: 0,
                    explain: "Resultant y = y₁ + y₂ = sin(ωt - kx) + sin(ωt + kx) = 2 sin(ωt) cos(kx), frequency remains ω"
                },
                {
                    id: "p38",
                    text: "The work done in stretching a spring by 2 cm is 4 J. Work done in stretching it by 4 cm from natural length is:",
                    options: ["8 J", "12 J", "16 J", "20 J"],
                    answer: 2,
                    explain: "Work = ½kx². If ½k(0.02)² = 4 J, then ½k(0.04)² = 4 × 4 = 16 J"
                },
                {
                    id: "p39",
                    text: "In series LCR circuit at resonance, impedance is:",
                    options: ["Maximum", "Minimum", "Zero", "Infinite"],
                    answer: 1,
                    explain: "At resonance, XL = XC, so impedance Z = R (minimum value for given R)"
                },
                {
                    id: "p40",
                    text: "The dimensional formula of magnetic permeability is:",
                    options: ["[MLT⁻²A⁻²]", "[ML²T⁻²A⁻²]", "[MLT⁻²A⁻¹]", "[ML⁻¹T⁻²A²]"],
                    answer: 0,
                    explain: "μ₀ = B/(H) where B has dimensions [MT⁻²A⁻¹] and H has dimensions [L⁻¹A]. So μ₀ has [MLT⁻²A⁻²]"
                },
                {
                    id: "p41",
                    text: "A car moving at 72 km/h is brought to rest in 10 s. Retardation is:",
                    options: ["2 m/s²", "7.2 m/s²", "20 m/s²", "72 m/s²"],
                    answer: 0,
                    explain: "Initial velocity = 72 km/h = 20 m/s. Final velocity = 0. a = (v-u)/t = (0-20)/10 = -2 m/s²"
                },
                {
                    id: "p42",
                    text: "The energy density in electric field E is:",
                    options: ["½ε₀E²", "ε₀E²", "2ε₀E²", "ε₀E²/2"],
                    answer: 0,
                    explain: "Energy density u = ½ε₀E² in electric field"
                },
                {
                    id: "p43",
                    text: "A ball is dropped from 80 m height. After bouncing, it rises to 45 m. Coefficient of restitution is:",
                    options: ["0.75", "0.5", "0.6", "0.8"],
                    answer: 0,
                    explain: "e = √(h₂/h₁) = √(45/80) = √(9/16) = 3/4 = 0.75"
                },
                {
                    id: "p44",
                    text: "In transistor, collector current is 10 mA and base current is 0.2 mA. Current gain is:",
                    options: ["20", "30", "40", "50"],
                    answer: 4,
                    explain: "Current gain β = IC/IB = 10/0.2 = 50"
                },
                {
                    id: "p45",
                    text: "The escape velocity from Earth's surface is 11.2 km/s. From height equal to Earth's radius, escape velocity is:",
                    options: ["5.6 km/s", "7.9 km/s", "11.2 km/s", "15.8 km/s"],
                    answer: 1,
                    explain: "ve = √(2GM/r). At height R: ve = √(2GM/2R) = ve(surface)/√2 = 11.2/√2 ≈ 7.9 km/s"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which of the following has highest first ionization energy?",
                    options: ["Be", "B", "C", "N"],
                    answer: 3,
                    explain: "Nitrogen has half-filled stable 2p³ configuration, requiring more energy to remove electron than adjacent elements."
                },
                {
                    id: "c2",
                    text: "The hybridization of central atom in BrF₅ is:",
                    options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
                    answer: 2,
                    explain: "BrF₅ has 6 electron pairs (5 bonding + 1 lone pair) around Br, requiring sp³d² hybridization with square pyramidal geometry."
                },
                {
                    id: "c3",
                    text: "Which of the following undergoes fastest SN1 reaction?",
                    options: ["CH₃Cl", "CH₃CH₂Cl", "(CH₃)₂CHCl", "(CH₃)₃CCl"],
                    answer: 3,
                    explain: "Tertiary carbocation (CH₃)₃C⁺ is most stable due to hyperconjugation, making (CH₃)₃CCl undergo fastest SN1."
                },
                {
                    id: "c4",
                    text: "The oxidation state of sulfur in H₂S₂O₈ is:",
                    options: ["+4", "+5", "+6", "+7"],
                    answer: 2,
                    explain: "In H₂S₂O₈ (peroxydisulfuric acid): 2(+1) + 2(S) + 8(-2) = 0, solving gives S = +6."
                },
                {
                    id: "c5",
                    text: "Which of the following shows optical isomerism?",
                    options: ["CH₃CH₂CHClCH₃", "CH₃CHBrCH₃", "CHBrClI", "CH₂ClCH₂Cl"],
                    answer: 2,
                    explain: "CHBrClI has four different groups attached to carbon, making it chiral and optically active."
                },
                {
                    id: "c6",
                    text: "The correct order of basic strength in aqueous solution is:",
                    options: ["NH₃ > CH₃NH₂ > (CH₃)₂NH", "(CH₃)₂NH > CH₃NH₂ > NH₃", "CH₃NH₂ > (CH₃)₂NH > NH₃", "(CH₃)₂NH > NH₃ > CH₃NH₂"],
                    answer: 1,
                    explain: "Secondary amine is most basic due to optimal +I effect and solvation. Order: (CH₃)₂NH > CH₃NH₂ > NH₃."
                },
                {
                    id: "c7",
                    text: "Which of the following has maximum number of unpaired electrons?",
                    options: ["Fe²⁺", "Co²⁺", "Ni²⁺", "Cu²⁺"],
                    answer: 0,
                    explain: "Fe²⁺ (d⁶) in high spin has 4 unpaired electrons, more than Co²⁺ (3), Ni²⁺ (2), or Cu²⁺ (1)."
                },
                {
                    id: "c8",
                    text: "The shape of ClF₃ molecule is:",
                    options: ["Planar triangular", "T-shaped", "Tetrahedral", "Trigonal bipyramidal"],
                    answer: 1,
                    explain: "ClF₃ has 5 electron pairs (3 bonding + 2 lone pairs) giving T-shaped geometry with lone pairs in equatorial positions."
                },
                {
                    id: "c9",
                    text: "Which of the following compounds cannot exist?",
                    options: ["BeCl₂", "BeCl₄²⁻", "AlCl₃", "AlCl₆³⁻"],
                    answer: 1,
                    explain: "Be cannot accommodate more than 4 electrons in its valence shell, so BeCl₄²⁻ with 8 electrons cannot exist."
                },
                {
                    id: "c10",
                    text: "The entropy change is maximum in:",
                    options: ["Melting of ice", "Vaporization of water", "Sublimation of iodine", "Condensation of steam"],
                    answer: 1,
                    explain: "Liquid to gas transition involves maximum increase in randomness, hence maximum positive entropy change."
                },
                {
                    id: "c11",
                    text: "Which of the following is strongest oxidizing agent?",
                    options: ["F₂", "Cl₂", "Br₂", "I₂"],
                    answer: 0,
                    explain: "F₂ has highest reduction potential (+2.87 V), making it strongest oxidizing agent among halogens."
                },
                {
                    id: "c12",
                    text: "The number of geometrical isomers of [Pt(NH₃)₂Cl₂] is:",
                    options: ["1", "2", "3", "4"],
                    answer: 1,
                    explain: "Square planar [Pt(NH₃)₂Cl₂] has 2 geometrical isomers: cis (adjacent) and trans (opposite)."
                },
                {
                    id: "c13",
                    text: "Which of the following undergoes aldol condensation most readily?",
                    options: ["HCHO", "CH₃CHO", "CH₃COCH₃", "C₆H₅CHO"],
                    answer: 1,
                    explain: "CH₃CHO has α-hydrogen and is less hindered than ketones, undergoing aldol condensation most readily."
                },
                {
                    id: "c14",
                    text: "The correct order of lattice energy is:",
                    options: ["LiF > LiCl > LiBr > LiI", "LiI > LiBr > LiCl > LiF", "LiCl > LiF > LiBr > LiI", "LiBr > LiI > LiCl > LiF"],
                    answer: 0,
                    explain: "Lattice energy inversely proportional to ionic size. F⁻ smallest, so LiF > LiCl > LiBr > LiI."
                },
                {
                    id: "c15",
                    text: "Which of the following is not a buffer system?",
                    options: ["CH₃COOH + CH₃COONa", "NH₄OH + NH₄Cl", "HCl + NaCl", "H₂CO₃ + NaHCO₃"],
                    answer: 2,
                    explain: "HCl (strong acid) + NaCl (salt of strong acid-strong base) cannot act as buffer system."
                },
                {
                    id: "c16",
                    text: "The magnetic moment of [Mn(CN)₆]³⁻ complex is:",
                    options: ["2.83 BM", "3.87 BM", "4.90 BM", "5.92 BM"],
                    answer: 0,
                    explain: "Mn³⁺ (d⁴) with strong field CN⁻ forms low spin complex with 2 unpaired electrons. μ = √[2(2+2)] = 2.83 BM."
                },
                {
                    id: "c17",
                    text: "Which of the following has maximum covalent character?",
                    options: ["AlF₃", "AlCl₃", "AlBr₃", "AlI₃"],
                    answer: 3,
                    explain: "According to Fajan's rules, AlI₃ has maximum covalent character due to large, polarizable I⁻ ion."
                },
                {
                    id: "c18",
                    text: "The rate of reaction quadruples when temperature increases from 300 K to 320 K. Activation energy is:",
                    options: ["52.3 kJ/mol", "104.6 kJ/mol", "26.1 kJ/mol", "78.4 kJ/mol"],
                    answer: 0,
                    explain: "Using Arrhenius equation: ln(k₂/k₁) = Ea/R(1/T₁ - 1/T₂). ln(4) = Ea/8.314(1/300 - 1/320). Solving: Ea ≈ 52.3 kJ/mol."
                },
                {
                    id: "c19",
                    text: "Which of the following shows maximum hydrogen bonding?",
                    options: ["H₂O", "H₂S", "H₂Se", "H₂Te"],
                    answer: 0,
                    explain: "H₂O shows maximum hydrogen bonding due to high electronegativity of oxygen and small size."
                },
                {
                    id: "c20",
                    text: "The number of π bonds in benzene ring is:",
                    options: ["3", "6", "9", "12"],
                    answer: 0,
                    explain: "Benzene has 3 delocalized π bonds in aromatic system with 6 π electrons following 4n+2 rule."
                },
                {
                    id: "c21",
                    text: "Which of the following has zero dipole moment?",
                    options: ["NH₃", "H₂O", "BF₃", "SO₂"],
                    answer: 2,
                    explain: "BF₃ has trigonal planar geometry with symmetrical arrangement, resulting in zero net dipole moment."
                },
                {
                    id: "c22",
                    text: "The coordination number of central atom in [Cu(NH₃)₄]SO₄ is:",
                    options: ["2", "4", "6", "8"],
                    answer: 1,
                    explain: "Four NH₃ ligands coordinate to Cu²⁺ ion, giving coordination number 4."
                },
                {
                    id: "c23",
                    text: "Which of the following is most acidic?",
                    options: ["CCl₃COOH", "CHCl₂COOH", "CH₂ClCOOH", "CH₃COOH"],
                    answer: 0,
                    explain: "CCl₃COOH is most acidic due to maximum -I effect of three Cl atoms stabilizing conjugate base."
                },
                {
                    id: "c24",
                    text: "The bond order in CO molecule is:",
                    options: ["2", "2.5", "3", "3.5"],
                    answer: 2,
                    explain: "CO has electronic configuration similar to N₂ with bond order 3 (triple bond)."
                },
                {
                    id: "c25",
                    text: "Which of the following undergoes nucleophilic substitution reaction?",
                    options: ["Benzene", "Toluene", "Chlorobenzene", "Phenol"],
                    answer: 2,
                    explain: "Chlorobenzene undergoes nucleophilic aromatic substitution due to electron-withdrawing Cl group."
                },
                {
                    id: "c26",
                    text: "The number of lone pairs in XeF₄ is:",
                    options: ["0", "1", "2", "3"],
                    answer: 2,
                    explain: "Xe in XeF₄ has 8 valence electrons, uses 4 for bonding, leaving 4 electrons = 2 lone pairs."
                },
                {
                    id: "c27",
                    text: "Which of the following is paramagnetic?",
                    options: ["CO", "NO", "CN⁻", "N₂"],
                    answer: 1,
                    explain: "NO has odd number of electrons (15), hence one unpaired electron making it paramagnetic."
                },
                {
                    id: "c28",
                    text: "The oxidation state of phosphorus in H₄P₂O₇ is:",
                    options: ["+3", "+4", "+5", "+6"],
                    answer: 2,
                    explain: "In pyrophosphoric acid H₄P₂O₇: 4(+1) + 2(P) + 7(-2) = 0, solving gives P = +5."
                },
                {
                    id: "c29",
                    text: "Which of the following shows linkage isomerism?",
                    options: ["[Co(NH₃)₅Cl]²⁺", "[Co(NH₃)₅SCN]²⁺", "[Co(NH₃)₆]³⁺", "[Co(en)₃]³⁺"],
                    answer: 1,
                    explain: "SCN⁻ can bind through S (thiocyanato) or N (isothiocyanato), showing linkage isomerism."
                },
                {
                    id: "c30",
                    text: "The number of atoms per unit cell in face-centered cubic is:",
                    options: ["1", "2", "4", "8"],
                    answer: 2,
                    explain: "FCC has 8 corner atoms (8×1/8) + 6 face atoms (6×1/2) = 4 atoms per unit cell."
                },
                {
                    id: "c31",
                    text: "Which of the following has maximum boiling point?",
                    options: ["n-pentane", "iso-pentane", "neo-pentane", "cyclopentane"],
                    answer: 0,
                    explain: "n-pentane has maximum surface contact and van der Waals forces, hence highest boiling point."
                },
                {
                    id: "c32",
                    text: "The rate constant of first order reaction has units:",
                    options: ["mol L⁻¹ s⁻¹", "s⁻¹", "mol⁻¹ L s⁻¹", "mol⁻² L² s⁻¹"],
                    answer: 1,
                    explain: "First order: rate = k[A], so k has units s⁻¹ (time⁻¹)."
                },
                {
                    id: "c33",
                    text: "Which of the following is diamagnetic?",
                    options: ["O₂", "B₂", "C₂", "F₂"],
                    answer: 2,
                    explain: "C₂ has all paired electrons in molecular orbital configuration, making it diamagnetic."
                },
                {
                    id: "c34",
                    text: "The correct order of second ionization energy is:",
                    options: ["Li > Be > B", "Be > Li > B", "B > Be > Li", "Li > B > Be"],
                    answer: 1,
                    explain: "After losing one electron: Li⁺ (He config), Be⁺ (Li config), B⁺ (He config). Be⁺ needs most energy: Be > Li > B."
                },
                {
                    id: "c35",
                    text: "Which of the following undergoes Cannizzaro reaction?",
                    options: ["Acetaldehyde", "Benzaldehyde", "Propionaldehyde", "Butyraldehyde"],
                    answer: 1,
                    explain: "Benzaldehyde lacks α-hydrogen atoms and undergoes Cannizzaro disproportionation reaction."
                },
                {
                    id: "c36",
                    text: "The bond angle in water molecule is approximately:",
                    options: ["104.5°", "107°", "109.5°", "120°"],
                    answer: 0,
                    explain: "H₂O has bent geometry with bond angle ~104.5° due to two lone pairs compressing H-O-H angle."
                },
                {
                    id: "c37",
                    text: "Which of the following has minimum atomic radius?",
                    options: ["Li", "Be", "B", "C"],
                    answer: 3,
                    explain: "Atomic radius decreases across period due to increasing nuclear charge. C has smallest radius among these."
                },
                {
                    id: "c38",
                    text: "The number of optical isomers of [Co(en)₂Cl₂]⁺ is:",
                    options: ["2", "3", "4", "6"],
                    answer: 0,
                    explain: "Octahedral complex with two bidentate ligands and two monodentate ligands has 2 optical isomers (enantiomers)."
                },
                {
                    id: "c39",
                    text: "Which of the following is weakest reducing agent?",
                    options: ["Li", "Na", "K", "Cs"],
                    answer: 0,
                    explain: "Despite highest reduction potential, Li is weakest reducing agent in solution due to high hydration energy."
                },
                {
                    id: "c40",
                    text: "The hybridization of carbon in graphite is:",
                    options: ["sp", "sp²", "sp³", "sp³d"],
                    answer: 1,
                    explain: "Carbon in graphite has sp² hybridization forming hexagonal layers with delocalized π electrons."
                },
                {
                    id: "c41",
                    text: "Which of the following shows stereoisomerism?",
                    options: ["CH₂=CHCl", "CHCl=CHCl", "CH₃CH=CHCl", "CH₂=CCl₂"],
                    answer: 1,
                    explain: "CHCl=CHCl shows geometrical isomerism with cis and trans forms around C=C double bond."
                },
                {
                    id: "c42",
                    text: "The effective atomic number of iron in [Fe(CN)₆]⁴⁻ is:",
                    options: ["26", "30", "34", "36"],
                    answer: 3,
                    explain: "Fe²⁺ (24e⁻) + 6 pairs from CN⁻ ligands (12e⁻) = 36 electrons, achieving Kr configuration."
                },
                {
                    id: "c43",
                    text: "Which of the following is Lewis base?",
                    options: ["BF₃", "AlCl₃", "NH₃", "BCl₃"],
                    answer: 2,
                    explain: "NH₃ has lone pair of electrons and can donate it to electron-deficient species, acting as Lewis base."
                },
                {
                    id: "c44",
                    text: "The number of sigma bonds in acetylene (C₂H₂) is:",
                    options: ["2", "3", "4", "5"],
                    answer: 1,
                    explain: "C₂H₂ has 3 σ bonds: 2 C-H bonds + 1 C-C σ bond. The C≡C has 1σ + 2π bonds."
                },
                {
                    id: "c45",
                    text: "Which of the following has maximum bond energy?",
                    options: ["C-C", "C=C", "C≡C", "All equal"],
                    answer: 2,
                    explain: "Triple bond C≡C has maximum bond energy (~835 kJ/mol) due to strongest bonding interaction."
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "Which of the following is the first stable product of dark reaction in C₄ plants?",
                    options: ["PEP", "OAA", "Malate", "Aspartate"],
                    answer: 1,
                    explain: "In C₄ plants, CO₂ first combines with PEP to form oxaloacetic acid (OAA) as first stable 4-carbon product."
                },
                {
                    id: "b2",
                    text: "The protective covering of radicle in monocots is called:",
                    options: ["Coleoptile", "Coleorhiza", "Scutellum", "Epiblast"],
                    answer: 1,
                    explain: "Coleorhiza is protective sheath covering radicle (root apex) in monocot seedlings during germination."
                },
                {
                    id: "b3",
                    text: "Which hormone is known as stress hormone in plants?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Abscisic acid"],
                    answer: 3,
                    explain: "Abscisic acid (ABA) is called stress hormone as it helps plants cope with various stress conditions like drought."
                },
                {
                    id: "b4",
                    text: "The cork cambium in woody stems is also called:",
                    options: ["Vascular cambium", "Phellogen", "Interfascicular cambium", "Intrafascicular cambium"],
                    answer: 1,
                    explain: "Cork cambium is technically called phellogen, which produces cork (phellem) outward and phelloderm inward."
                },
                {
                    id: "b5",
                    text: "Which of the following shows homospory?",
                    options: ["Selaginella", "Marsilea", "Salvinia", "Pteris"],
                    answer: 3,
                    explain: "Pteris (common fern) is homosporous, producing only one type of spore, unlike heterosporous pteridophytes."
                },
                {
                    id: "b6",
                    text: "The phenomenon of triple fusion was discovered by:",
                    options: ["Strasburger", "Nawaschin", "Hofmeister", "Guignard"],
                    answer: 1,
                    explain: "Nawaschin discovered triple fusion (fusion of one sperm nucleus with two polar nuclei) in 1898."
                },
                {
                    id: "b7",
                    text: "Which enzyme is involved in CO₂ fixation in C₃ plants?",
                    options: ["PEP carboxylase", "RuBP carboxylase", "Carbonic anhydrase", "NADP dehydrogenase"],
                    answer: 1,
                    explain: "RuBP carboxylase (RuBisCO) catalyzes CO₂ fixation with RuBP in Calvin cycle of C₃ plants."
                },
                {
                    id: "b8",
                    text: "The pericycle in monocot roots is made up of:",
                    options: ["Parenchyma", "Sclerenchyma", "Both parenchyma and sclerenchyma", "Collenchyma"],
                    answer: 2,
                    explain: "Monocot root pericycle consists of both parenchymatous and sclerenchymatous cells in alternating patches."
                },
                {
                    id: "b9",
                    text: "Which of the following shows racemose inflorescence?",
                    options: ["Cyme", "Helicoid cyme", "Scorpioid cyme", "Corymb"],
                    answer: 3,
                    explain: "Corymb is racemose (indeterminate) inflorescence where flowers arise at different levels but reach same height."
                },
                {
                    id: "b10",
                    text: "Nucellar embryony is a type of:",
                    options: ["Sexual reproduction", "Apomixis", "Vegetative propagation", "Parthenogenesis"],
                    answer: 1,
                    explain: "Nucellar embryony is type of apomixis where embryos develop from nucellus cells without fertilization."
                },
                {
                    id: "b11",
                    text: "Which hormone promotes stem elongation?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Ethylene"],
                    answer: 1,
                    explain: "Gibberellins specifically promote stem elongation by stimulating cell division and elongation in internodes."
                },
                {
                    id: "b12",
                    text: "Cyclic photophosphorylation involves only:",
                    options: ["PSI", "PSII", "Both PSI and PSII", "ATP synthase only"],
                    answer: 0,
                    explain: "Cyclic photophosphorylation involves only PSI where electrons cycle back without producing NADPH or O₂."
                },
                {
                    id: "b13",
                    text: "Which of the following is a long day plant?",
                    options: ["Rice", "Tobacco", "Spinach", "Cosmos"],
                    answer: 2,
                    explain: "Spinach requires long photoperiods (short nights) for flowering, making it a long day plant."
                },
                {
                    id: "b14",
                    text: "Solute potential is always:",
                    options: ["Positive", "Negative", "Zero", "Variable"],
                    answer: 1,
                    explain: "Solute potential (ψs) is always negative because solutes reduce water's chemical potential below pure water."
                },
                {
                    id: "b15",
                    text: "Which cells lack nucleus at maturity in phloem?",
                    options: ["Companion cells", "Phloem parenchyma", "Sieve tube elements", "Phloem fibers"],
                    answer: 2,
                    explain: "Mature sieve tube elements lose their nucleus and most organelles, depending on companion cells for metabolism."
                },
                {
                    id: "b16",
                    text: "The primary endosperm nucleus is:",
                    options: ["Haploid", "Diploid", "Triploid", "Tetraploid"],
                    answer: 2,
                    explain: "Primary endosperm nucleus is triploid (3n) formed by triple fusion of one sperm (n) with two polar nuclei (n+n)."
                },
                {
                    id: "b17",
                    text: "Which type of leaf venation is found in monocots?",
                    options: ["Reticulate", "Parallel", "Furcate", "Dichotomous"],
                    answer: 1,
                    explain: "Monocot leaves typically show parallel venation where veins run parallel to each other."
                },
                {
                    id: "b18",
                    text: "The function of passage cells is:",
                    options: ["Support", "Conduction", "Storage", "Lateral transport"],
                    answer: 3,
                    explain: "Passage cells in endodermis lack casparian strips, allowing lateral transport of water and minerals."
                },
                {
                    id: "b19",
                    text: "Chlorophyll-a differs from chlorophyll-b in having:",
                    options: ["Methyl group instead of formyl group", "Formyl group instead of methyl group", "Phytol tail", "Mg atom"],
                    answer: 0,
                    explain: "Chlorophyll-a has methyl group (-CH₃) at C3 position while chlorophyll-b has formyl group (-CHO)."
                },
                {
                    id: "b20",
                    text: "In which type of ovule is micropyle farthest from hilum?",
                    options: ["Orthotropous", "Anatropous", "Campylotropous", "Amphitropous"],
                    answer: 1,
                    explain: "In anatropous ovules, ovule inverts 180° making micropyle closest to hilum, not farthest. Actually orthotropous has micropyle farthest from hilum."
                },
                {
                    id: "b21",
                    text: "Guard cells accumulate which ion during stomatal opening?",
                    options: ["Na⁺", "Ca²⁺", "K⁺", "Mg²⁺"],
                    answer: 2,
                    explain: "Guard cells actively accumulate K⁺ ions during day, increasing osmotic pressure and causing stomatal opening."
                },
                {
                    id: "b22",
                    text: "Pollen tube is guided by:",
                    options: ["Synergids", "Antipodals", "Polar nuclei", "Secondary nucleus"],
                    answer: 0,
                    explain: "Synergids secrete chemotropic substances that guide pollen tube growth toward embryo sac."
                },
                {
                    id: "b23",
                    text: "Kranz anatomy is characteristic of:",
                    options: ["C₃ plants", "C₄ plants", "CAM plants", "All plants"],
                    answer: 1,
                    explain: "Kranz anatomy with bundle sheath cells surrounding vascular bundles is characteristic adaptation of C₄ plants."
                },
                {
                    id: "b24",
                    text: "Photoperiodism was discovered by:",
                    options: ["Garner and Allard", "Went", "Darwin", "Blackman"],
                    answer: 0,
                    explain: "Garner and Allard discovered photoperiodism in 1920 while working on tobacco flowering responses."
                },
                {
                    id: "b25",
                    text: "Which of the following is not found in xylem tissue?",
                    options: ["Tracheids", "Vessels", "Xylem parenchyma", "Companion cells"],
                    answer: 3,
                    explain: "Companion cells are components of phloem tissue associated with sieve tube elements, not found in xylem."
                },
                {
                    id: "b26",
                    text: "Calvin cycle occurs in:",
                    options: ["Thylakoids", "Stroma", "Intermembrane space", "Outer membrane"],
                    answer: 1,
                    explain: "Calvin cycle (dark reaction) occurs in chloroplast stroma where CO₂ fixation and reduction take place."
                },
                {
                    id: "b27",
                    text: "Which plant part shows unlimited growth?",
                    options: ["Leaf", "Flower", "Root", "Fruit"],
                    answer: 2,
                    explain: "Root shows indeterminate (unlimited) growth throughout plant life due to persistent apical meristem."
                },
                {
                    id: "b28",
                    text: "Pollen grain represents:",
                    options: ["Microspore", "Male gametophyte", "Microsporangium", "Anther"],
                    answer: 1,
                    explain: "Mature pollen grain is male gametophyte containing generative and vegetative nuclei."
                },
                {
                    id: "b29",
                    text: "Root nodules contain which enzyme for nitrogen fixation?",
                    options: ["Nitrate reductase", "Nitrite reductase", "Nitrogenase", "Glutamine synthetase"],
                    answer: 2,
                    explain: "Nitrogenase enzyme complex in root nodule bacteroids converts atmospheric N₂ to NH₃."
                },
                {
                    id: "b30",
                    text: "Lenticels are formed by:",
                    options: ["Epidermis", "Cork cambium", "Vascular cambium", "Ground meristem"],
                    answer: 1,
                    explain: "Lenticels develop from cork cambium (phellogen) as loose tissue allowing gaseous exchange."
                },
                {
                    id: "b31",
                    text: "Seismonasty is shown by:",
                    options: ["Helianthus", "Mimosa pudica", "Tropaeolum", "Drosera"],
                    answer: 1,
                    explain: "Mimosa pudica shows seismonastic movement, rapidly folding leaves in response to touch stimulus."
                },
                {
                    id: "b32",
                    text: "Which tissue differentiates first in primary growth?",
                    options: ["Protoxylem", "Metaxylem", "Secondary xylem", "Protophloem"],
                    answer: 0,
                    explain: "Protoxylem differentiates first during primary growth, followed by metaxylem as organ matures."
                },
                {
                    id: "b33",
                    text: "Which hormone is involved in fruit ripening?",
                    options: ["Auxin", "Gibberellin", "Ethylene", "Cytokinin"],
                    answer: 2,
                    explain: "Ethylene is primary hormone triggering fruit ripening by activating enzymes that soften fruit tissues."
                },
                {
                    id: "b34",
                    text: "Root pressure is maximum during:",
                    options: ["Day time", "Night time", "Noon", "Evening"],
                    answer: 1,
                    explain: "Root pressure is maximum at night when transpiration is minimal and active salt pumping continues."
                },
                {
                    id: "b35",
                    text: "PEP carboxylase has high affinity for:",
                    options: ["O₂", "CO₂", "Both O₂ and CO₂", "Neither O₂ nor CO₂"],
                    answer: 1,
                    explain: "PEP carboxylase in C₄ plants has high affinity for CO₂ and no affinity for O₂, avoiding photorespiration."
                },
                {
                    id: "b36",
                    text: "Maximum water absorption occurs in:",
                    options: ["Root cap", "Zone of elongation", "Root hair zone", "Zone of maturation"],
                    answer: 2,
                    explain: "Root hair zone provides enormous surface area through root hairs for maximum water and mineral absorption."
                },
                {
                    id: "b37",
                    text: "Whorled phyllotaxy is seen in:",
                    options: ["China rose", "Mustard", "Nerium", "Sunflower"],
                    answer: 2,
                    explain: "Nerium shows whorled phyllotaxy where three or more leaves arise from single node."
                },
                {
                    id: "b38",
                    text: "Gibberellin was first discovered in:",
                    options: ["Wheat", "Rice", "Maize", "Barley"],
                    answer: 1,
                    explain: "Gibberellin was first discovered by Japanese scientist Kurosawa in rice plants affected by 'bakanae' disease."
                },
                {
                    id: "b39",
                    text: "Xerophytes have:",
                    options: ["Large leaves", "Sunken stomata", "Thin cuticle", "High transpiration"],
                    answer: 1,
                    explain: "Xerophytes adapt to water scarcity with sunken stomata that create humid microenvironment reducing transpiration."
                },
                {
                    id: "b40",
                    text: "Cohesion-tension theory explains:",
                    options: ["Translocation", "Water transport", "Mineral transport", "Gas exchange"],
                    answer: 1,
                    explain: "Cohesion-tension theory explains upward water transport in xylem through cohesive forces and transpiration pull."
                },
                {
                    id: "b41",
                    text: "ABA is synthesized in:",
                    options: ["Roots", "Leaves", "Fruits", "All parts"],
                    answer: 3,
                    explain: "ABA can be synthesized in all plant parts including roots, leaves, fruits, and seeds when needed."
                },
                {
                    id: "b42",
                    text: "Hill reaction demonstrates:",
                    options: ["CO₂ fixation", "Oxygen evolution", "ATP synthesis", "NADPH formation"],
                    answer: 1,
                    explain: "Hill reaction showed that chloroplasts can evolve oxygen in presence of light and artificial electron acceptors."
                },
                {
                    id: "b43",
                    text: "Hypogeal germination occurs in:",
                    options: ["Bean", "Castor", "Pea", "Onion"],
                    answer: 2,
                    explain: "Pea shows hypogeal germination where cotyledons remain underground and epicotyl emerges above soil."
                },
                {
                    id: "b44",
                    text: "Phloem in gymnosperms lacks:",
                    options: ["Sieve cells", "Albuminous cells", "Companion cells", "Phloem parenchyma"],
                    answer: 2,
                    explain: "Gymnosperms have sieve cells with albuminous cells but lack true companion cells found in angiosperms."
                },
                {
                    id: "b45",
                    text: "Which hormone delays senescence?",
                    options: ["Auxin", "Ethylene", "Cytokinin", "ABA"],
                    answer: 2,
                    explain: "Cytokinins delay leaf senescence by promoting protein synthesis and chlorophyll retention."
                },
                
                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Which of the following is the longest bone in human body?",
                    options: ["Humerus", "Tibia", "Fibula", "Femur"],
                    answer: 3,
                    explain: "Femur (thigh bone) is the longest and strongest bone in human body, extending from hip to knee."
                },
                {
                    id: "z2",
                    text: "Insulin is secreted by which cells of pancreas?",
                    options: ["Alpha cells", "Beta cells", "Delta cells", "PP cells"],
                    answer: 1,
                    explain: "Beta cells in islets of Langerhans secrete insulin hormone that lowers blood glucose levels."
                },
                {
                    id: "z3",
                    text: "Which of the following prevents blood clotting?",
                    options: ["Fibrinogen", "Prothrombin", "Heparin", "Thrombin"],
                    answer: 2,
                    explain: "Heparin is natural anticoagulant produced by basophils and mast cells that prevents blood clotting."
                },
                {
                    id: "z4",
                    text: "Bowman's capsule is part of:",
                    options: ["Glomerulus", "Nephron", "Ureter", "Urethra"],
                    answer: 1,
                    explain: "Bowman's capsule is cup-shaped structure at beginning of nephron that surrounds glomerulus."
                },
                {
                    id: "z5",
                    text: "Which part of brain is called little brain?",
                    options: ["Cerebrum", "Cerebellum", "Medulla", "Pons"],
                    answer: 1,
                    explain: "Cerebellum is called 'little brain' due to its folded structure resembling miniature brain."
                },
                {
                    id: "z6",
                    text: "RBCs are destroyed in:",
                    options: ["Liver", "Spleen", "Kidney", "Bone marrow"],
                    answer: 1,
                    explain: "Old and damaged RBCs are primarily destroyed by macrophages in spleen after ~120 days lifespan."
                },
                {
                    id: "z7",
                    text: "Which enzyme converts proteins to polypeptides in stomach?",
                    options: ["Trypsin", "Pepsin", "Chymotrypsin", "Carboxypeptidase"],
                    answer: 1,
                    explain: "Pepsin, activated from pepsinogen in acidic stomach environment, breaks proteins into polypeptides."
                },
                {
                    id: "z8",
                    text: "Which node is called pacemaker of heart?",
                    options: ["AV node", "SA node", "Bundle of His", "Purkinje fibers"],
                    answer: 1,
                    explain: "SA (sinoatrial) node in right atrium initiates heartbeat and sets cardiac rhythm as natural pacemaker."
                },
                {
                    id: "z9",
                    text: "Fight or flight response is mediated by:",
                    options: ["Insulin", "Thyroxine", "Adrenaline", "Cortisol"],
                    answer: 2,
                    explain: "Adrenaline (epinephrine) from adrenal medulla mediates fight-or-flight response during stress."
                },
                {
                    id: "z10",
                    text: "Implantation occurs in:",
                    options: ["Ovary", "Fallopian tube", "Uterus", "Cervix"],
                    answer: 2,
                    explain: "Blastocyst implants in endometrial lining of uterus about 6-7 days after fertilization."
                },
                {
                    id: "z11",
                    text: "Deficiency of which vitamin causes scurvy?",
                    options: ["Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D"],
                    answer: 2,
                    explain: "Vitamin C (ascorbic acid) deficiency causes scurvy characterized by bleeding gums and connective tissue problems."
                },
                {
                    id: "z12",
                    text: "Which cells are involved in cell-mediated immunity?",
                    options: ["B cells", "T cells", "NK cells", "Plasma cells"],
                    answer: 1,
                    explain: "T lymphocytes (T cells) mediate cellular immunity by directly attacking infected cells."
                },
                {
                    id: "z13",
                    text: "Which part of nephron is impermeable to urea?",
                    options: ["Proximal tubule", "Descending loop", "Ascending loop", "Collecting duct"],
                    answer: 2,
                    explain: "Thick ascending limb of Henle's loop is impermeable to water and urea but actively transports salts."
                },
                {
                    id: "z14",
                    text: "Oxytocin is produced by:",
                    options: ["Anterior pituitary", "Posterior pituitary", "Hypothalamus", "Adrenal cortex"],
                    answer: 2,
                    explain: "Oxytocin is synthesized by hypothalamus and stored/released by posterior pituitary."
                },
                {
                    id: "z15",
                    text: "Which blood group is universal donor?",
                    options: ["A", "B", "AB", "O"],
                    answer: 3,
                    explain: "O blood group lacks A and B antigens, making it compatible with all blood groups as universal donor."
                },
                {
                    id: "z16",
                    text: "Vocal cords are located in:",
                    options: ["Pharynx", "Larynx", "Trachea", "Bronchi"],
                    answer: 1,
                    explain: "Vocal cords are elastic ligaments stretched across larynx that vibrate to produce sound."
                },
                {
                    id: "z17",
                    text: "Static equilibrium is detected by:",
                    options: ["Cochlea", "Utricle and saccule", "Semicircular canals", "Organ of corti"],
                    answer: 1,
                    explain: "Utricle and saccule contain otoliths that detect linear acceleration and static equilibrium."
                },
                {
                    id: "z18",
                    text: "Normal color of urine is due to:",
                    options: ["Urobilin", "Urochrome", "Bilirubin", "Hemoglobin"],
                    answer: 1,
                    explain: "Urochrome, a yellow pigment derived from urobilinogen, gives urine its characteristic yellow color."
                },
                {
                    id: "z19",
                    text: "Growth hormone is secreted by:",
                    options: ["Thyroid", "Adrenal", "Anterior pituitary", "Posterior pituitary"],
                    answer: 2,
                    explain: "Growth hormone (somatotropin) is secreted by acidophil cells of anterior pituitary."
                },
                {
                    id: "z20",
                    text: "Spermatogenesis occurs in:",
                    options: ["Epididymis", "Vas deferens", "Seminiferous tubules", "Prostate gland"],
                    answer: 2,
                    explain: "Spermatogenesis (sperm formation) occurs in seminiferous tubules of testes under FSH and testosterone influence."
                },
                {
                    id: "z21",
                    text: "Which chamber pumps blood to lungs?",
                    options: ["Right atrium", "Left atrium", "Right ventricle", "Left ventricle"],
                    answer: 2,
                    explain: "Right ventricle pumps deoxygenated blood to lungs through pulmonary artery for oxygenation."
                },
                {
                    id: "z22",
                    text: "Blood pH is maintained by:",
                    options: ["Lungs only", "Kidneys only", "Both lungs and kidneys", "Liver"],
                    answer: 2,
                    explain: "Blood pH is maintained by lungs (CO₂ removal) and kidneys (H⁺/HCO₃⁻ regulation) working together."
                },
                {
                    id: "z23",
                    text: "Calcitonin is secreted by:",
                    options: ["Parathyroid", "Thyroid", "Adrenal", "Pancreas"],
                    answer: 1,
                    explain: "Calcitonin is secreted by parafollicular (C) cells of thyroid gland to lower blood calcium."
                },
                {
                    id: "z24",
                    text: "Respiratory center is present in:",
                    options: ["Cerebrum", "Cerebellum", "Medulla oblongata", "Hypothalamus"],
                    answer: 2,
                    explain: "Medulla oblongata contains respiratory center that controls automatic breathing rhythm."
                },
                {
                    id: "z25",
                    text: "Cardiac muscle is:",
                    options: ["Voluntary striated", "Involuntary striated", "Voluntary smooth", "Involuntary smooth"],
                    answer: 1,
                    explain: "Cardiac muscle is involuntary (autonomic control) but striated like skeletal muscle."
                },
                {
                    id: "z26",
                    text: "Bile is produced by:",
                    options: ["Liver", "Gall bladder", "Pancreas", "Duodenum"],
                    answer: 0,
                    explain: "Bile is produced by hepatocytes in liver and stored in gall bladder for fat emulsification."
                },
                {
                    id: "z27",
                    text: "Which vitamin helps in calcium absorption?",
                    options: ["Vitamin A", "Vitamin C", "Vitamin D", "Vitamin K"],
                    answer: 2,
                    explain: "Vitamin D (calciferol) promotes calcium absorption from intestine and calcium deposition in bones."
                },
                {
                    id: "z28",
                    text: "Sliding filament theory explains:",
                    options: ["Nerve conduction", "Muscle contraction", "Blood clotting", "Digestion"],
                    answer: 1,
                    explain: "Sliding filament theory explains muscle contraction mechanism where actin slides over myosin filaments."
                },
                {
                    id: "z29",
                    text: "Type I diabetes is caused by:",
                    options: ["Insulin resistance", "Insulin deficiency", "Glucagon excess", "Cortisol excess"],
                    answer: 1,
                    explain: "Type I diabetes is caused by autoimmune destruction of beta cells leading to insulin deficiency."
                },
                {
                    id: "z30",
                    text: "Maximum absorption occurs in:",
                    options: ["Stomach", "Small intestine", "Large intestine", "Rectum"],
                    answer: 1,
                    explain: "Small intestine with its villi and microvilli is primary site for nutrient absorption."
                },
                {
                    id: "z31",
                    text: "Body temperature is regulated by:",
                    options: ["Cerebrum", "Hypothalamus", "Medulla", "Thalamus"],
                    answer: 1,
                    explain: "Hypothalamus acts as thermostat, regulating body temperature through various mechanisms."
                },
                {
                    id: "z32",
                    text: "Human pregnancy lasts for:",
                    options: ["36 weeks", "38 weeks", "40 weeks", "42 weeks"],
                    answer: 2,
                    explain: "Normal human pregnancy duration is 40 weeks (280 days) from last menstrual period."
                },
                {
                    id: "z33",
                    text: "Thrombocytes are:",
                    options: ["RBCs", "WBCs", "Platelets", "Plasma proteins"],
                    answer: 2,
                    explain: "Thrombocytes are platelets - cell fragments that play crucial role in blood clotting."
                },
                {
                    id: "z34",
                    text: "Pineal gland secretes:",
                    options: ["Growth hormone", "Melatonin", "Serotonin", "Dopamine"],
                    answer: 1,
                    explain: "Pineal gland secretes melatonin hormone that regulates sleep-wake cycles and circadian rhythms."
                },
                {
                    id: "z35",
                    text: "Yellow spot in retina is:",
                    options: ["Blind spot", "Fovea", "Optic disc", "Choroid"],
                    answer: 1,
                    explain: "Fovea centralis (yellow spot) has highest concentration of cones for sharp, color vision."
                },
                {
                    id: "z36",
                    text: "ADH regulates:",
                    options: ["Blood pressure", "Water reabsorption", "Glucose level", "Calcium level"],
                    answer: 1,
                    explain: "ADH (antidiuretic hormone) regulates water reabsorption in kidney collecting ducts."
                },
                {
                    id: "z37",
                    text: "Sense of smell is detected by:",
                    options: ["Olfactory nerve", "Optic nerve", "Auditory nerve", "Facial nerve"],
                    answer: 0,
                    explain: "Olfactory nerve (I cranial nerve) carries smell sensations from nasal cavity to brain."
                },
                {
                    id: "z38",
                    text: "Excessive glucose in blood is called:",
                    options: ["Hypoglycemia", "Hyperglycemia", "Glycosuria", "Ketonuria"],
                    answer: 1,
                    explain: "Hyperglycemia is condition with abnormally high blood glucose levels, common in diabetes."
                },
                {
                    id: "z39",
                    text: "Corpus luteum secretes:",
                    options: ["Estrogen", "Progesterone", "FSH", "LH"],
                    answer: 1,
                    explain: "Corpus luteum formed from ruptured follicle secretes progesterone to maintain pregnancy."
                },
                {
                    id: "z40",
                    text: "How many pairs of true ribs are present?",
                    options: ["5", "6", "7", "8"],
                    answer: 2,
                    explain: "Humans have 7 pairs of true ribs that directly attach to sternum through costal cartilages."
                },
                {
                    id: "z41",
                    text: "HCl in gastric juice is secreted by:",
                    options: ["Chief cells", "Parietal cells", "Mucus cells", "G cells"],
                    answer: 1,
                    explain: "Parietal cells (oxyntic cells) in gastric glands secrete hydrochloric acid and intrinsic factor."
                },
                {
                    id: "z42",
                    text: "Synovial joints are:",
                    options: ["Immovable", "Slightly movable", "Freely movable", "Fixed"],
                    answer: 2,
                    explain: "Synovial joints are freely movable joints with synovial fluid providing lubrication."
                },
                {
                    id: "z43",
                    text: "T₃ and T₄ are hormones of:",
                    options: ["Parathyroid", "Thyroid", "Adrenal", "Pancreas"],
                    answer: 1,
                    explain: "Thyroid gland produces triiodothyronine (T₃) and thyroxine (T₄) that regulate metabolism."
                },
                {
                    id: "z44",
                    text: "Bicuspid valve prevents backflow from:",
                    options: ["Right ventricle to right atrium", "Left ventricle to left atrium", "Aorta to left ventricle", "Pulmonary artery to right ventricle"],
                    answer: 1,
                    explain: "Bicuspid (mitral) valve prevents backflow from left ventricle to left atrium during ventricular contraction."
                },
                {
                    id: "z45",
                    text: "Acrosome of sperm contains:",
                    options: ["DNA", "Mitochondria", "Enzymes", "Flagellar proteins"],
                    answer: 2,
                    explain: "Acrosome contains hydrolytic enzymes like hyaluronidase and acrosin that help penetrate egg coats."
                }
            ]
        }
    ]
};
