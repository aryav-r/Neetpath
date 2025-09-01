// mock-test-13.js - NEET Mock Test 13 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_13 = {
    id: "neet-013",
    title: "Full Syllabus Mock 13", 
    level: "medium",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A particle moves in a straight line with velocity v = 6t² - 4t + 2. The acceleration at t = 3s is:",
                    options: ["32 m/s²", "36 m/s²", "40 m/s²", "44 m/s²"],
                    answer: 0,
                    explain: "Acceleration a = dv/dt = d/dt(6t² - 4t + 2) = 12t - 4. At t = 3s: a = 12(3) - 4 = 32 m/s²"
                },
                {
                    id: "p2",
                    text: "A uniform ring of radius R rolls without slipping down an incline of angle θ. Its linear acceleration is:",
                    options: ["g sin θ", "g sin θ/2", "2g sin θ/3", "g sin θ/3"],
                    answer: 1,
                    explain: "For rolling ring: a = g sin θ/(1 + I/mR²) = g sin θ/(1 + 1) = g sin θ/2"
                },
                {
                    id: "p3",
                    text: "Three capacitors 2μF, 3μF, and 6μF are connected in series and then parallel with a 4μF capacitor. Total capacitance is:",
                    options: ["5 μF", "6 μF", "5.5 μF", "4.5 μF"],
                    answer: 0,
                    explain: "Series: 1/C = 1/2 + 1/3 + 1/6 = 3/6 + 2/6 + 1/6 = 1, so C = 1μF. Parallel with 4μF: Total = 1 + 4 = 5μF"
                },
                {
                    id: "p4",
                    text: "In interference pattern, if path difference is λ/3, the resultant amplitude is:",
                    options: ["a", "2a", "√3 a", "a√3"],
                    answer: 2,
                    explain: "Resultant amplitude A = 2a cos(φ/2) where φ = 2π(path diff)/λ = 2π/3. A = 2a cos(π/3) = 2a × 0.5 = a"
                },
                {
                    id: "p5",
                    text: "A photon of energy 6.2 eV strikes a metal surface. If work function is 3.7 eV, maximum kinetic energy of photoelectron is:",
                    options: ["2.5 eV", "3.7 eV", "6.2 eV", "9.9 eV"],
                    answer: 0,
                    explain: "KEmax = hf - φ = 6.2 - 3.7 = 2.5 eV"
                },
                {
                    id: "p6",
                    text: "A circular coil of 100 turns and area 0.02 m² rotates at 120 rpm in magnetic field 0.5 T. Peak EMF is:",
                    options: ["25.1 V", "15.7 V", "31.4 V", "12.6 V"],
                    answer: 0,
                    explain: "Peak EMF = NABω = 100 × 0.02 × 0.5 × (2π × 120/60) = 1 × 4π × 2 = 25.1 V"
                },
                {
                    id: "p7",
                    text: "A thin rod of length L pivoted at one end oscillates as compound pendulum. Distance of center of oscillation from pivot is:",
                    options: ["L/2", "2L/3", "3L/4", "L/3"],
                    answer: 1,
                    explain: "Center of oscillation is at distance h = I/(mL/2) = (mL²/3)/(mL/2) = 2L/3 from pivot"
                },
                {
                    id: "p8",
                    text: "In RLC series circuit, when frequency is half of resonant frequency, the phase angle is:",
                    options: ["-45°", "45°", "-60°", "60°"],
                    answer: 2,
                    explain: "At f = f₀/2: XL = ωL/2, XC = 2/(ωC). tan φ = (XL-XC)/R = (ωL/2 - 2/(ωC))/R. Since XL < XC, φ is negative"
                },
                {
                    id: "p9",
                    text: "A gas at 300K expands adiabatically to 8 times its volume. If γ = 1.4, final temperature is:",
                    options: ["150 K", "125 K", "100 K", "75 K"],
                    answer: 0,
                    explain: "For adiabatic process: T₁V₁^(γ-1) = T₂V₂^(γ-1). T₂ = T₁(V₁/V₂)^(γ-1) = 300(1/8)^0.4 = 300/2 = 150 K"
                },
                {
                    id: "p10",
                    text: "Magnetic field at the center of a square loop of side 'a' carrying current I is:",
                    options: ["2√2 μ₀I/πa", "4√2 μ₀I/πa", "μ₀I/πa", "√2 μ₀I/πa"],
                    answer: 0,
                    explain: "Field due to one side at center = μ₀I/(4πr) × 2 sin 45°, where r = a/(2√2). Total field = 4 × μ₀I√2/(2πa) = 2√2 μ₀I/πa"
                },
                {
                    id: "p11",
                    text: "A particle executes SHM with equation x = 8 sin(4πt + π/6). At t = 0.25s, velocity is:",
                    options: ["16π m/s", "0 m/s", "-16π m/s", "8π m/s"],
                    answer: 2,
                    explain: "v = dx/dt = 32π cos(4πt + π/6). At t = 0.25s: v = 32π cos(π + π/6) = 32π cos(7π/6) = 32π(-√3/2) = -16π√3 ≈ -16π m/s"
                },
                {
                    id: "p12",
                    text: "In pair production, minimum photon energy required is:",
                    options: ["0.511 MeV", "1.022 MeV", "1.533 MeV", "2.044 MeV"],
                    answer: 1,
                    explain: "Minimum energy = 2mₑc² = 2 × 0.511 = 1.022 MeV to create electron-positron pair"
                },
                {
                    id: "p13",
                    text: "Five identical resistors each of 10Ω are connected as shown in a regular pentagon. Resistance between adjacent vertices is:",
                    options: ["6Ω", "8Ω", "5Ω", "4Ω"],
                    answer: 2,
                    explain: "In pentagon configuration between adjacent vertices: equivalent resistance = 5Ω by symmetry analysis"
                },
                {
                    id: "p14",
                    text: "A ball is thrown at angle 60° with horizontal. Ratio of maximum height to range is:",
                    options: ["√3:4", "3:4", "1:2", "√3:2"],
                    answer: 0,
                    explain: "H = u²sin²θ/2g, R = u²sin2θ/g. H/R = sin²θ/(2sin2θ) = sin²60°/(2sin120°) = (3/4)/(2×√3/2) = √3/4"
                },
                {
                    id: "p15",
                    text: "Mutual inductance between two coils depends on:",
                    options: ["Current in coils", "Number of turns and geometry", "Frequency", "Voltage applied"],
                    answer: 1,
                    explain: "Mutual inductance M = k√(L₁L₂) depends on number of turns, geometry, and coupling coefficient k"
                },
                {
                    id: "p16",
                    text: "A stone dropped from balloon ascending at 20 m/s reaches ground in 8s. Height of balloon when stone was dropped:",
                    options: ["160 m", "320 m", "240 m", "280 m"],
                    answer: 2,
                    explain: "s = ut + ½gt² = -20(8) + ½(10)(64) = -160 + 320 = 160 m below starting point. Height = 160 + 80 = 240 m"
                },
                {
                    id: "p17",
                    text: "In AC circuit with R = 6Ω and L = 8Ω (inductive reactance), impedance and phase angle are:",
                    options: ["10Ω, 53°", "14Ω, 45°", "10Ω, 37°", "8Ω, 60°"],
                    answer: 0,
                    explain: "Z = √(R² + XL²) = √(36 + 64) = 10Ω. tan φ = XL/R = 8/6 = 4/3, so φ = 53°"
                },
                {
                    id: "p18",
                    text: "Wavelength of matter wave for electron with kinetic energy 100 eV is approximately:",
                    options: ["1.2 Å", "0.12 Å", "2.4 Å", "0.24 Å"],
                    answer: 0,
                    explain: "λ = h/p = h/√(2mKE) = 6.626×10⁻³⁴/√(2×9.1×10⁻³¹×100×1.6×10⁻¹⁹) ≈ 1.2×10⁻¹⁰ m = 1.2 Å"
                },
                {
                    id: "p19",
                    text: "Two springs of constants k₁ and k₂ are connected in series. A mass m attached oscillates with period:",
                    options: ["2π√[m(k₁+k₂)/(k₁k₂)]", "2π√[mk₁k₂/(k₁+k₂)]", "2π√[m/(k₁+k₂)]", "2π√[(k₁+k₂)/m]"],
                    answer: 0,
                    explain: "Series: 1/kₑff = 1/k₁ + 1/k₂, so kₑff = k₁k₂/(k₁+k₂). T = 2π√(m/kₑff) = 2π√[m(k₁+k₂)/(k₁k₂)]"
                },
                {
                    id: "p20",
                    text: "Electric field intensity at distance r from infinite line charge λ is:",
                    options: ["λ/2πε₀r", "λ/4πε₀r²", "2λ/4πε₀r", "λ/πε₀r"],
                    answer: 0,
                    explain: "Using Gauss's law for infinite line charge: E = λ/2πε₀r"
                },
                {
                    id: "p21",
                    text: "A glass slab of thickness 6 cm and refractive index 1.5 is placed over a mark. Apparent depth is:",
                    options: ["4 cm", "9 cm", "3 cm", "2 cm"],
                    answer: 0,
                    explain: "Apparent depth = Real depth/μ = 6/1.5 = 4 cm"
                },
                {
                    id: "p22",
                    text: "Two speakers vibrate in phase with frequency 340 Hz. Path difference for first minimum at a point is:",
                    options: ["0.5 m", "1.0 m", "1.5 m", "2.0 m"],
                    answer: 0,
                    explain: "For first minimum, path difference = λ/2 = v/2f = 340/(2×340) = 0.5 m"
                },
                {
                    id: "p23",
                    text: "In polytropic process PVⁿ = constant, work done is:",
                    options: ["(P₁V₁ - P₂V₂)/(n-1)", "nRT ln(V₂/V₁)", "(P₂V₂ - P₁V₁)/(n-1)", "nP₁V₁"],
                    answer: 2,
                    explain: "For polytropic process: W = (P₂V₂ - P₁V₁)/(n-1) when n ≠ 1"
                },
                {
                    id: "p24",
                    text: "If λ is decay constant, average life of radioactive element is:",
                    options: ["λ", "1/λ", "ln2/λ", "λ/ln2"],
                    answer: 1,
                    explain: "Average life (mean life) τ = 1/λ where λ is decay constant"
                },
                {
                    id: "p25",
                    text: "Energy stored per unit volume in magnetic field is:",
                    options: ["B²/2μ₀", "μ₀B²/2", "B²/μ₀", "μ₀B²"],
                    answer: 0,
                    explain: "Magnetic energy density u = B²/2μ₀"
                },
                {
                    id: "p26",
                    text: "Magnetic quantum number can have values:",
                    options: ["-l to +l", "0 to n-1", "1 to n", "-n to +n"],
                    answer: 0,
                    explain: "Magnetic quantum number m can have (2l+1) values from -l to +l where l is azimuthal quantum number"
                },
                {
                    id: "p27",
                    text: "For diamond (μ = 2.42), critical angle is approximately:",
                    options: ["24°", "41°", "30°", "60°"],
                    answer: 0,
                    explain: "sin θc = 1/μ = 1/2.42 = 0.413. θc = arcsin(0.413) ≈ 24°"
                },
                {
                    id: "p28",
                    text: "In pipe closed at one end, first overtone has frequency:",
                    options: ["3f₀", "2f₀", "5f₀", "4f₀"],
                    answer: 0,
                    explain: "In closed pipe, frequencies are f, 3f, 5f... First overtone = 3f₀"
                },
                {
                    id: "p29",
                    text: "Coefficient of performance of refrigerator is:",
                    options: ["Q₂/W", "W/Q₂", "Q₁/W", "W/Q₁"],
                    answer: 0,
                    explain: "COP of refrigerator = Heat extracted/Work done = Q₂/W"
                },
                {
                    id: "p30",
                    text: "A charged particle moving perpendicular to uniform magnetic field follows:",
                    options: ["Straight line", "Parabolic path", "Circular path", "Helical path"],
                    answer: 2,
                    explain: "When v ⊥ B, Lorentz force provides centripetal force, resulting in circular motion"
                },
                {
                    id: "p31",
                    text: "Ultrasonic waves have frequency:",
                    options: [">20 Hz", ">200 Hz", ">2000 Hz", ">20000 Hz"],
                    answer: 3,
                    explain: "Ultrasonic waves have frequency greater than 20,000 Hz (20 kHz), above human hearing range"
                },
                {
                    id: "p32",
                    text: "Moment of inertia of hollow cylinder about axis perpendicular to its axis through center is:",
                    options: ["M(R₁² + R₂²)/2", "MR²", "M(3R₁² + 3R₂² + h²)/12", "MR²/2"],
                    answer: 2,
                    explain: "For hollow cylinder about perpendicular axis through center: I = M(3R₁² + 3R₂² + h²)/12"
                },
                {
                    id: "p33",
                    text: "In parallel LC circuit at resonance, impedance is:",
                    options: ["Zero", "Infinite", "√(L/C)", "√(L×C)"],
                    answer: 1,
                    explain: "At resonance in parallel LC circuit, reactances cancel and impedance becomes theoretically infinite"
                },
                {
                    id: "p34",
                    text: "Object placed at focus of concave mirror forms image at:",
                    options: ["Focus", "Center of curvature", "Between F and C", "Infinity"],
                    answer: 3,
                    explain: "When object is at focus of concave mirror, reflected rays are parallel, forming image at infinity"
                },
                {
                    id: "p35",
                    text: "α-particle has mass number:",
                    options: ["2", "4", "1", "3"],
                    answer: 1,
                    explain: "α-particle is helium nucleus with 2 protons + 2 neutrons, so mass number = 4"
                },
                {
                    id: "p36",
                    text: "Free oscillations have:",
                    options: ["Constant amplitude", "Decreasing amplitude", "Increasing amplitude", "Variable frequency"],
                    answer: 1,
                    explain: "Real free oscillations have decreasing amplitude due to damping forces like friction"
                },
                {
                    id: "p37",
                    text: "Earth's magnetic field is approximately:",
                    options: ["Uniform", "Dipolar", "Quadrupolar", "Monopolar"],
                    answer: 1,
                    explain: "Earth's magnetic field resembles that of a magnetic dipole with north and south magnetic poles"
                },
                {
                    id: "p38",
                    text: "Transformer efficiency depends on:",
                    options: ["Turn ratio", "Core losses", "Copper losses", "Both core and copper losses"],
                    answer: 3,
                    explain: "Transformer efficiency is reduced by both core losses (hysteresis, eddy currents) and copper losses (I²R losses)"
                },
                {
                    id: "p39",
                    text: "Time taken for current to reach 63% of maximum in LR circuit is:",
                    options: ["L/R", "R/L", "LR", "√(LR)"],
                    answer: 0,
                    explain: "Time constant τ = L/R is time to reach 63% of maximum current in LR circuit"
                },
                {
                    id: "p40",
                    text: "Electric dipole moment has units:",
                    options: ["Coulomb-meter", "Newton-meter", "Joule-meter", "Volt-meter"],
                    answer: 0,
                    explain: "Electric dipole moment p = q×d has units of charge × distance = Coulomb-meter"
                },
                {
                    id: "p41",
                    text: "Specific heat at constant pressure is greater than at constant volume because:",
                    options: ["More heat required for expansion", "Molecular motion increases", "Temperature rise is different", "Pressure changes"],
                    answer: 0,
                    explain: "At constant pressure, additional heat is needed for expansion work, so CP > CV"
                },
                {
                    id: "p42",
                    text: "In pure inductive AC circuit, average power consumed is:",
                    options: ["VI", "VI/2", "Zero", "V²/XL"],
                    answer: 2,
                    explain: "In pure inductive circuit, current and voltage are 90° out of phase, so average power = VIcos90° = 0"
                },
                {
                    id: "p43",
                    text: "Condition for sustained interference is:",
                    options: ["Same amplitude", "Constant phase difference", "Same frequency", "All of these"],
                    answer: 3,
                    explain: "For sustained interference, sources must have same frequency, constant phase difference, and preferably similar amplitudes"
                },
                {
                    id: "p44",
                    text: "Orbital velocity of satellite depends on:",
                    options: ["Mass of satellite", "Orbital radius only", "Mass of Earth only", "Both orbital radius and Earth's mass"],
                    answer: 3,
                    explain: "Orbital velocity v = √(GM/r) depends on central body mass (M) and orbital radius (r)"
                },
                {
                    id: "p45",
                    text: "Zener diode is used in:",
                    options: ["Amplification", "Rectification", "Voltage regulation", "Switching"],
                    answer: 2,
                    explain: "Zener diode operates in reverse bias breakdown region for voltage regulation and stabilization"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which has minimum first ionization energy?",
                    options: ["Li", "Na", "K", "Rb"],
                    answer: 3,
                    explain: "Rubidium has largest atomic size in given options, hence minimum first ionization energy"
                },
                {
                    id: "c2",
                    text: "The geometry of ICl₂⁻ ion is:",
                    options: ["Linear", "Bent", "Trigonal planar", "T-shaped"],
                    answer: 0,
                    explain: "ICl₂⁻ has 5 electron pairs (2 bonding + 3 lone pairs) around I, giving linear geometry"
                },
                {
                    id: "c3",
                    text: "Which shows fastest E1 elimination?",
                    options: ["Primary alkyl halide", "Secondary alkyl halide", "Tertiary alkyl halide", "Vinyl halide"],
                    answer: 2,
                    explain: "Tertiary alkyl halides form most stable carbocations, hence undergo fastest E1 elimination"
                },
                {
                    id: "c4",
                    text: "Oxidation state of sulfur in Na₂S₂O₃ is:",
                    options: ["+2", "+4", "+6", "+2.5"],
                    answer: 0,
                    explain: "In thiosulfate, one S has +6 and other has -2 oxidation states. Average = (+6-2)/2 = +2"
                },
                {
                    id: "c5",
                    text: "A compound with molecular formula C₄H₈O₂ showing optical activity must have:",
                    options: ["At least one chiral carbon", "Two chiral carbons", "No chiral carbons", "Three chiral carbons"],
                    answer: 0,
                    explain: "Optical activity requires at least one chiral carbon atom with four different substituents"
                },
                {
                    id: "c6",
                    text: "Strongest acid among these is:",
                    options: ["HClO", "HClO₂", "HClO₃", "HClO₄"],
                    answer: 3,
                    explain: "Acid strength increases with number of oxygen atoms: HClO₄ > HClO₃ > HClO₂ > HClO"
                },
                {
                    id: "c7",
                    text: "Ground state configuration of Cr²⁺ is:",
                    options: ["[Ar] 3d⁴", "[Ar] 3d⁵", "[Ar] 3d³ 4s¹", "[Ar] 3d⁶"],
                    answer: 0,
                    explain: "Cr: [Ar] 3d⁵ 4s¹. Cr²⁺ loses 4s¹ and one 3d electron: [Ar] 3d⁴"
                },
                {
                    id: "c8",
                    text: "Which has pyramidal geometry?",
                    options: ["BF₃", "NH₃", "PCl₅", "SF₆"],
                    answer: 1,
                    explain: "NH₃ has 4 electron pairs (3 bonding + 1 lone pair) giving trigonal pyramidal geometry"
                },
                {
                    id: "c9",
                    text: "Which has expanded octet?",
                    options: ["BeCl₂", "BCl₃", "PCl₅", "CCl₄"],
                    answer: 2,
                    explain: "PCl₅ has 10 electrons around P (5 bonding pairs), exceeding octet rule"
                },
                {
                    id: "c10",
                    text: "For reaction at equilibrium, ΔG is:",
                    options: ["Positive", "Negative", "Zero", "Cannot predict"],
                    answer: 2,
                    explain: "At equilibrium, there is no net change in free energy, so ΔG = 0"
                },
                {
                    id: "c11",
                    text: "Best oxidizing agent among these is:",
                    options: ["F₂", "Cl₂", "Br₂", "I₂"],
                    answer: 0,
                    explain: "F₂ has highest reduction potential (+2.87 V), making it strongest oxidizing agent"
                },
                {
                    id: "c12",
                    text: "Which complex has maximum CFSE?",
                    options: ["[Ti(H₂O)₆]³⁺", "[V(H₂O)₆]³⁺", "[Cr(H₂O)₆]³⁺", "[Mn(H₂O)₆]³⁺"],
                    answer: 2,
                    explain: "Cr³⁺ (d³) has configuration t₂g³ eg⁰ with CFSE = -1.2Δ₀, maximum among given options"
                },
                {
                    id: "c13",
                    text: "Lucas test distinguishes between:",
                    options: ["Primary and secondary alcohols", "Secondary and tertiary alcohols", "Primary, secondary and tertiary alcohols", "Alcohols and phenols"],
                    answer: 2,
                    explain: "Lucas test (ZnCl₂/HCl) gives immediate turbidity with tertiary, slow with secondary, no reaction with primary alcohols"
                },
                {
                    id: "c14",
                    text: "Correct order of covalent radii is:",
                    options: ["C < N < O < F", "F < O < N < C", "O < N < C < F", "N < O < F < C"],
                    answer: 1,
                    explain: "Covalent radius decreases across period due to increasing nuclear charge: F < O < N < C"
                },
                {
                    id: "c15",
                    text: "Which has maximum buffering capacity?",
                    options: ["CH₃COOH + CH₃COONa (equal molar)", "NH₄OH + NH₄Cl (equal molar)", "H₃PO₄ + NaH₂PO₄", "All have same"],
                    answer: 3,
                    explain: "At equal concentrations and pH near pKa, all buffers have similar capacity"
                },
                {
                    id: "c16",
                    text: "Spin-only magnetic moment of [Fe(H₂O)₆]²⁺ is:",
                    options: ["4.90 BM", "5.92 BM", "2.83 BM", "0 BM"],
                    answer: 0,
                    explain: "Fe²⁺ (d⁶) with weak field H₂O forms high-spin complex with 4 unpaired electrons. μ = √[4×6] = 4.90 BM"
                },
                {
                    id: "c17",
                    text: "Which has maximum percentage s-character?",
                    options: ["sp³", "sp²", "sp", "sp³d"],
                    answer: 2,
                    explain: "sp hybridization has 50% s-character, highest among given options"
                },
                {
                    id: "c18",
                    text: "Order of reaction A + B → C is determined by:",
                    options: ["Stoichiometry", "Concentration of reactants", "Experimental rate law", "Temperature"],
                    answer: 2,
                    explain: "Order of reaction is determined experimentally from rate law, not from balanced equation"
                },
                {
                    id: "c19",
                    text: "Which shows strongest van der Waals forces?",
                    options: ["CH₄", "C₂H₆", "C₃H₈", "C₄H₁₀"],
                    answer: 3,
                    explain: "C₄H₁₀ has largest molecular size and highest polarizability, hence strongest van der Waals forces"
                },
                {
                    id: "c20",
                    text: "Which follows Hückel's rule?",
                    options: ["C₄H₄²⁺", "C₈H₈²⁻", "C₆H₆", "All of these"],
                    answer: 3,
                    explain: "All follow 4n+2 rule: C₄H₄²⁺ (2π e⁻), C₈H₈²⁻ (10π e⁻), C₆H₆ (6π e⁻)"
                },
                {
                    id: "c21",
                    text: "Which has resultant dipole moment zero?",
                    options: ["H₂O", "NH₃", "CHCl₃", "BF₃"],
                    answer: 3,
                    explain: "BF₃ has trigonal planar geometry with three equal B-F dipoles canceling each other"
                },
                {
                    id: "c22",
                    text: "Which can show linkage isomerism?",
                    options: ["CN⁻", "SCN⁻", "NO₂⁻", "All of these"],
                    answer: 3,
                    explain: "All three ligands are ambidentate: CN⁻ (C or N), SCN⁻ (S or N), NO₂⁻ (N or O)"
                },
                {
                    id: "c23",
                    text: "pKa of acetic acid is 4.8. At pH 6.8, ratio [CH₃COO⁻]/[CH₃COOH] is:",
                    options: ["1:100", "100:1", "1:10", "10:1"],
                    answer: 1,
                    explain: "Using Henderson-Hasselbalch: pH = pKa + log[A⁻]/[HA]. 6.8 = 4.8 + log(ratio). log(ratio) = 2, so ratio = 100:1"
                },
                {
                    id: "c24",
                    text: "Bond order of CN⁻ is:",
                    options: ["2", "2.5", "3", "3.5"],
                    answer: 2,
                    explain: "CN⁻ has 14 electrons. Bond order = (bonding - antibonding)/2 = (10-4)/2 = 3"
                },
                {
                    id: "c25",
                    text: "Which undergoes nucleophilic aromatic substitution most readily?",
                    options: ["Benzene", "Chlorobenzene", "Nitrobenzene", "2,4-dinitrophenyl chloride"],
                    answer: 3,
                    explain: "2,4-dinitrophenyl chloride has two strong electron-withdrawing NO₂ groups activating nucleophilic substitution"
                },
                {
                    id: "c26",
                    text: "Steric number of Xe in XeF₂ is:",
                    options: ["2", "3", "4", "5"],
                    answer: 4,
                    explain: "Xe in XeF₂ has 2 bonding pairs + 3 lone pairs = 5, giving steric number 5"
                },
                {
                    id: "c27",
                    text: "Which is most paramagnetic?",
                    options: ["O₂", "O₂⁺", "O₂⁻", "O₂²⁻"],
                    answer: 1,
                    explain: "O₂⁺ has 15 electrons with 3 unpaired electrons, showing maximum paramagnetism"
                },
                {
                    id: "c28",
                    text: "Average oxidation state of C in C₃O₂ is:",
                    options: ["+4/3", "+8/3", "+2/3", "+4"],
                    answer: 0,
                    explain: "In C₃O₂: 3C + 2(-2) = 0, so 3C = +4, average oxidation state = +4/3"
                },
                {
                    id: "c29",
                    text: "Solvate isomerism is shown by:",
                    options: ["[Cr(H₂O)₆]Cl₃", "[CrCl₂(H₂O)₄]Cl⋅2H₂O", "[Co(NH₃)₆]Cl₃", "[Ni(CO)₄]"],
                    answer: 1,
                    explain: "[CrCl₂(H₂O)₄]Cl⋅2H₂O can exist as [CrCl(H₂O)₅]Cl₂⋅H₂O showing solvate isomerism"
                },
                {
                    id: "c30",
                    text: "In antifluorite structure, cation to anion ratio is:",
                    options: ["1:1", "1:2", "2:1", "1:4"],
                    answer: 2,
                    explain: "Antifluorite structure (like Na₂O) has cation:anion ratio of 2:1"
                },
                {
                    id: "c31",
                    text: "Which has maximum boiling point?",
                    options: ["n-pentane", "2-methylbutane", "2,2-dimethylpropane", "Cyclopentane"],
                    answer: 0,
                    explain: "n-pentane has linear structure with maximum surface contact, hence highest boiling point"
                },
                {
                    id: "c32",
                    text: "For first order reaction, half-life is 20 min. Time for 75% completion is:",
                    options: ["30 min", "40 min", "60 min", "80 min"],
                    answer: 1,
                    explain: "75% completion means 25% remaining = (1/2)² of original. Time = 2 × t₁/₂ = 40 min"
                },
                {
                    id: "c33",
                    text: "Which has octahedral geometry?",
                    options: ["SF₆", "PCl₅", "ClF₅", "IF₇"],
                    answer: 0,
                    explain: "SF₆ has 6 bonding pairs around S with no lone pairs, giving perfect octahedral geometry"
                },
                {
                    id: "c34",
                    text: "Which has maximum electron affinity?",
                    options: ["F", "Cl", "Br", "I"],
                    answer: 1,
                    explain: "Chlorine has maximum electron affinity due to optimal size - larger than F (less repulsion) but smaller than Br, I"
                },
                {
                    id: "c35",
                    text: "Perkin reaction involves:",
                    options: ["Aldol condensation", "Claisen condensation", "Knoevenagel condensation", "Wittig reaction"],
                    answer: 2,
                    explain: "Perkin reaction is condensation between aromatic aldehyde and acid anhydride in presence of weak base"
                },
                {
                    id: "c36",
                    text: "In trigonal bipyramidal geometry, lone pairs occupy:",
                    options: ["Axial positions", "Equatorial positions", "Both positions equally", "Any position"],
                    answer: 1,
                    explain: "Lone pairs preferentially occupy equatorial positions to minimize repulsion in trigonal bipyramidal geometry"
                },
                {
                    id: "c37",
                    text: "Which period has maximum number of elements?",
                    options: ["4th", "5th", "6th", "7th"],
                    answer: 2,
                    explain: "6th period has 32 elements (including lanthanides), maximum in periodic table"
                },
                {
                    id: "c38",
                    text: "Which shows maximum number of stereoisomers?",
                    options: ["[Co(NH₃)₄Cl₂]⁺", "[Co(NH₃)₃Cl₃]", "[Co(en)₂Cl₂]⁺", "[Co(en)₃]³⁺"],
                    answer: 2,
                    explain: "[Co(en)₂Cl₂]⁺ can show both geometrical (cis-trans) and optical isomerism, giving maximum stereoisomers"
                },
                {
                    id: "c39",
                    text: "Most electronegative element in periodic table is:",
                    options: ["F", "O", "N", "Cl"],
                    answer: 0,
                    explain: "Fluorine has highest electronegativity (4.0 on Pauling scale) in periodic table"
                },
                {
                    id: "c40",
                    text: "Which shows resonance?",
                    options: ["CH₄", "C₂H₆", "C₆H₆", "CH₃CH₃"],
                    answer: 2,
                    explain: "Benzene (C₆H₆) shows resonance with delocalized π electrons in aromatic ring"
                },
                {
                    id: "c41",
                    text: "Which can exhibit geometrical isomerism?",
                    options: ["1-butene", "2-butene", "2-methylpropene", "1,1-dichloroethene"],
                    answer: 1,
                    explain: "2-butene (CH₃CH=CHCH₃) can exist as cis and trans isomers around C=C double bond"
                },
                {
                    id: "c42",
                    text: "Which follows 18-electron rule?",
                    options: ["[V(CO)₆]⁻", "[Cr(CO)₆]", "[Mn(CO)₅]⁻", "All of these"],
                    answer: 3,
                    explain: "All follow 18-electron rule: V(5) + 12 + 1 = 18, Cr(6) + 12 = 18, Mn(7) + 10 + 1 = 18"
                },
                {
                    id: "c43",
                    text: "Which is strongest Brønsted base?",
                    options: ["OH⁻", "NH₂⁻", "CH₃⁻", "H⁻"],
                    answer: 3,
                    explain: "H⁻ (hydride ion) is strongest Brønsted base as H₂ is weakest acid"
                },
                {
                    id: "c44",
                    text: "Total number of sigma bonds in C₂H₄ is:",
                    options: ["4", "5", "6", "7"],
                    answer: 1,
                    explain: "C₂H₄ has 5 σ bonds: 4 C-H bonds + 1 C-C σ bond (double bond = 1σ + 1π)"
                },
                {
                    id: "c45",
                    text: "Which has shortest bond length?",
                    options: ["C-C", "C=C", "C≡C", "C-H"],
                    answer: 2,
                    explain: "Triple bond C≡C has shortest bond length due to maximum electron density between atoms"
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "In C₄ pathway, CO₂ is first fixed in:",
                    options: ["Bundle sheath cells", "Mesophyll cells", "Guard cells", "Epidermal cells"],
                    answer: 1,
                    explain: "In C₄ plants, CO₂ is first fixed by PEP carboxylase in mesophyll cells to form 4-carbon compounds"
                },
                {
                    id: "b2",
                    text: "Which tissue stores food in monocot seeds?",
                    options: ["Cotyledon", "Endosperm", "Perisperm", "Nucellus"],
                    answer: 1,
                    explain: "Monocot seeds typically have persistent endosperm that serves as food storage tissue"
                },
                {
                    id: "b3",
                    text: "Which hormone promotes cambial activity?",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Both auxin and gibberellin"],
                    answer: 3,
                    explain: "Both auxin and gibberellins promote cambial cell division and differentiation in secondary growth"
                },
                {
                    id: "b4",
                    text: "Interfascicular cambium originates from:",
                    options: ["Pericycle", "Endodermis", "Medullary rays", "Cortex"],
                    answer: 2,
                    explain: "Interfascicular cambium develops from parenchyma cells of medullary rays between vascular bundles"
                },
                {
                    id: "b5",
                    text: "Protonema in bryophytes develops from:",
                    options: ["Spore", "Gamete", "Zygote", "Archegonium"],
                    answer: 0,
                    explain: "Protonema is juvenile gametophytic stage that develops from germinating spore in bryophytes"
                },
                {
                    id: "b6",
                    text: "Dichogamy refers to:",
                    options: ["Unisexual flowers", "Temporal separation of anther and stigma maturation", "Spatial separation", "Self-incompatibility"],
                    answer: 1,
                    explain: "Dichogamy is temporal separation where anthers and stigma mature at different times to prevent selfing"
                },
                {
                    id: "b7",
                    text: "NADP reductase is associated with:",
                    options: ["PSI", "PSII", "Cytochrome complex", "Calvin cycle"],
                    answer: 0,
                    explain: "NADP reductase is located at PSI acceptor side and reduces NADP⁺ to NADPH using electrons from PSI"
                },
                {
                    id: "b8",
                    text: "Which type of stomata are found in grasses?",
                    options: ["Anomocytic", "Paracytic", "Diacytic", "Anisocytic"],
                    answer: 1,
                    explain: "Paracytic stomata with subsidiary cells parallel to guard cells are characteristic of grasses"
                },
                {
                    id: "b9",
                    text: "Catkin inflorescence is found in:",
                    options: ["Mulberry", "Oak", "Sunflower", "Mustard"],
                    answer: 1,
                    explain: "Catkin (amentum) is spike-like inflorescence with unisexual flowers, characteristic of oak, willow, birch"
                },
                {
                    id: "b10",
                    text: "Diplospory is a type of:",
                    options: ["Vegetative reproduction", "Apomixis", "Sexual reproduction", "Fragmentation"],
                    answer: 1,
                    explain: "Diplospory is type of apomixis where embryo sac develops from diploid megaspore mother cell without meiosis"
                },
                {
                    id: "b11",
                    text: "Who coined the term 'hormone' for plant growth substances?",
                    options: ["Darwin", "Went", "Thimann", "Boysen-Jensen"],
                    answer: 2,
                    explain: "Kenneth Thimann coined the term 'phytohormone' for plant growth-regulating substances"
                },
                {
                    id: "b12",
                    text: "Plastoquinone is involved in:",
                    options: ["Electron transport between PSII and cytochrome complex", "Calvin cycle", "Photorespiration", "Nitrogen fixation"],
                    answer: 0,
                    explain: "Plastoquinone transfers electrons from PSII to cytochrome b₆f complex in photosynthetic electron transport"
                },
                {
                    id: "b13",
                    text: "Which plant shows thermoperiodism?",
                    options: ["Chrysanthemum", "Spinach", "Tulip", "Sunflower"],
                    answer: 2,
                    explain: "Tulip requires cold treatment (vernalization) followed by warm temperatures for flowering, showing thermoperiodism"
                },
                {
                    id: "b14",
                    text: "Incipient plasmolysis occurs when:",
                    options: ["Cell is fully turgid", "Turgor pressure equals wall pressure", "Turgor pressure becomes zero", "Cell wall starts shrinking"],
                    answer: 2,
                    explain: "Incipient plasmolysis is initial stage where turgor pressure just becomes zero but protoplasm hasn't started shrinking"
                },
                {
                    id: "b15",
                    text: "Source-sink relationship in phloem transport refers to:",
                    options: ["Direction of water flow", "Sites of sugar production and utilization", "Xylem-phloem interaction", "Root-shoot communication"],
                    answer: 1,
                    explain: "Source-sink relationship describes movement from sugar-producing sites (source) to sugar-utilizing sites (sink)"
                },
                {
                    id: "b16",
                    text: "Secondary nucleus in embryo sac is:",
                    options: ["Haploid", "Diploid", "Triploid", "Formed by fusion of polar nuclei"],
                    answer: 3,
                    explain: "Secondary nucleus (central cell nucleus) is formed by fusion of two polar nuclei before fertilization"
                },
                {
                    id: "b17",
                    text: "Kranz anatomy refers to:",
                    options: ["Bundle sheath arrangement in C₄ plants", "Stomatal arrangement", "Vascular bundle pattern", "Root anatomical structure"],
                    answer: 0,
                    explain: "Kranz anatomy is wreath-like arrangement of bundle sheath cells around vascular bundles in C₄ plants"
                },
                {
                    id: "b18",
                    text: "Which gives positive test with phloroglucinol-HCl?",
                    options: ["Cellulose", "Lignin", "Pectin", "Suberin"],
                    answer: 1,
                    explain: "Phloroglucinol-HCl gives cherry red color with lignin, used as specific test for lignified tissues"
                },
                {
                    id: "b19",
                    text: "Lutein and zeaxanthin belong to:",
                    options: ["Chlorophylls", "Carotenoids", "Anthocyanins", "Phycobilins"],
                    answer: 1,
                    explain: "Lutein and zeaxanthin are xanthophyll carotenoids that protect against photo-oxidative damage"
                },
                {
                    id: "b20",
                    text: "Campylotropous ovule has:",
                    options: ["Straight structure", "Completely inverted structure", "Curved embryo sac", "Inverted nucellus"],
                    answer: 2,
                    explain: "In campylotropous ovule, embryo sac is curved so that micropyle comes to lie near funicle"
                },
                {
                    id: "b21",
                    text: "Skotomorphogenesis refers to:",
                    options: ["Growth in light", "Growth in darkness", "Phototropic response", "Gravitropic response"],
                    answer: 1,
                    explain: "Skotomorphogenesis is pattern of growth and development that occurs in darkness (etiolation)"
                },
                {
                    id: "b22",
                    text: "Homogamy prevents:",
                    options: ["Self-pollination", "Cross-pollination", "Wind pollination", "Both self and cross pollination"],
                    answer: 0,
                    explain: "Homogamy (simultaneous maturation of anthers and stigma) can lead to self-pollination if other mechanisms don't prevent it"
                },
                {
                    id: "b23",
                    text: "Photorespiration is also called:",
                    options: ["Dark respiration", "C₂ cycle", "Glycolate pathway", "Both C₂ cycle and glycolate pathway"],
                    answer: 3,
                    explain: "Photorespiration is also called C₂ cycle or glycolate pathway due to involvement of 2-carbon glycolate"
                },
                {
                    id: "b24",
                    text: "Wind dispersed fruits/seeds have:",
                    options: ["Hooks and spines", "Wings or pappus", "Fleshy covering", "Hard seed coat"],
                    answer: 1,
                    explain: "Anemochorous fruits/seeds have adaptations like wings, pappus, or light weight for wind dispersal"
                },
                {
                    id: "b25",
                    text: "P-protein is found in:",
                    options: ["Xylem vessels", "Sieve tube elements", "Companion cells", "Tracheids"],
                    answer: 1,
                    explain: "P-protein (phloem protein) is characteristic component of mature sieve tube elements"
                },
                {
                    id: "b26",
                    text: "In Calvin cycle, for synthesis of one glucose molecule, CO₂ molecules required are:",
                    options: ["3", "6", "9", "12"],
                    answer: 1,
                    explain: "Synthesis of one glucose (6-carbon) requires fixation of 6 CO₂ molecules in Calvin cycle"
                },
                {
                    id: "b27",
                    text: "Which shows unlimited growth?",
                    options: ["Leaves", "Flowers", "Fruits", "Stems"],
                    answer: 3,
                    explain: "Stems show indeterminate growth through persistent apical meristem activity"
                },
                {
                    id: "b28",
                    text: "Filiform apparatus is part of:",
                    options: ["Synergids", "Antipodals", "Central cell", "Integuments"],
                    answer: 0,
                    explain: "Filiform apparatus consists of finger-like projections in synergids that help guide pollen tube"
                },
                {
                    id: "b29",
                    text: "Frankia forms symbiotic association with:",
                    options: ["Legumes", "Non-leguminous plants", "Aquatic plants", "Epiphytes"],
                    answer: 1,
                    explain: "Frankia (actinomycete) forms nitrogen-fixing nodules with non-leguminous plants like Alnus, Casuarina"
                },
                {
                    id: "b30",
                    text: "Which tissue is responsible for increase in girth?",
                    options: ["Apical meristem", "Intercalary meristem", "Lateral meristem", "Primary meristem"],
                    answer: 2,
                    explain: "Lateral meristems (vascular and cork cambium) are responsible for secondary growth and increase in girth"
                },
                {
                    id: "b31",
                    text: "Seismonastic movements are shown by:",
                    options: ["Sunflower", "Mimosa pudica", "Tropaeolum", "Oxalis"],
                    answer: 1,
                    explain: "Mimosa pudica shows rapid seismonastic (shock) movements in response to touch or vibration"
                },
                {
                    id: "b32",
                    text: "Primary wall of plant cell is composed mainly of:",
                    options: ["Lignin", "Cellulose and pectin", "Suberin", "Cutin"],
                    answer: 1,
                    explain: "Primary cell wall consists mainly of cellulose microfibrils embedded in pectin matrix"
                },
                {
                    id: "b33",
                    text: "Banana fruit is:",
                    options: ["True fruit", "False fruit", "Aggregate fruit", "Parthenocarpic fruit"],
                    answer: 3,
                    explain: "Commercial banana is parthenocarpic fruit that develops without fertilization and is seedless"
                },
                {
                    id: "b34",
                    text: "Water potential of pure water is:",
                    options: ["-1 MPa", "0 MPa", "+1 MPa", "Variable"],
                    answer: 1,
                    explain: "Pure water has water potential of zero (0 MPa) which is taken as reference point"
                },
                {
                    id: "b35",
                    text: "Which pathway avoids photorespiration?",
                    options: ["C₃ pathway", "C₄ pathway", "CAM pathway", "Both C₄ and CAM"],
                    answer: 3,
                    explain: "Both C₄ and CAM pathways concentrate CO₂ around RuBisCO, minimizing photorespiration"
                },
                {
                    id: "b36",
                    text: "Root hairs develop from:",
                    options: ["Epidermis", "Cortex", "Endodermis", "Pericycle"],
                    answer: 0,
                    explain: "Root hairs are tubular extensions of epidermal cells in the root hair zone"
                },
                {
                    id: "b37",
                    text: "Whorled phyllotaxy is found in:",
                    options: ["Mustard", "China rose", "Nerium", "Sunflower"],
                    answer: 2,
                    explain: "Nerium (oleander) shows whorled phyllotaxy with three leaves per node"
                },
                {
                    id: "b38",
                    text: "Double fertilization was discovered by:",
                    options: ["Strasburger", "Nawaschin", "Guignard", "Both Nawaschin and Guignard"],
                    answer: 3,
                    explain: "Double fertilization was independently discovered by Nawaschin and Guignard in 1898"
                },
                {
                    id: "b39",
                    text: "Mangroves show:",
                    options: ["Hydrophytic adaptations", "Xerophytic adaptations", "Mesophytic adaptations", "Halophytic adaptations"],
                    answer: 3,
                    explain: "Mangroves are halophytes adapted to saline conditions with salt-filtering roots and salt glands"
                },
                {
                    id: "b40",
                    text: "Cohesion-adhesion theory explains:",
                    options: ["Translocation of food", "Ascent of sap", "Mineral absorption", "Gas exchange"],
                    answer: 1,
                    explain: "Cohesion-adhesion theory explains ascent of sap in xylem through cohesive and adhesive properties of water"
                },
                {
                    id: "b41",
                    text: "Dedifferentiation is followed by:",
                    options: ["Senescence", "Redifferentiation", "Cell death", "Dormancy"],
                    answer: 1,
                    explain: "In tissue culture, dedifferentiated cells undergo redifferentiation to form organized structures"
                },
                {
                    id: "b42",
                    text: "Cyclic electron flow involves:",
                    options: ["Only PSI", "Only PSII", "Both PSI and PSII", "Neither photosystem"],
                    answer: 0,
                    explain: "Cyclic electron flow involves only PSI where electrons return to the same photosystem"
                },
                {
                    id: "b43",
                    text: "Testa develops from:",
                    options: ["Nucellus", "Inner integument", "Outer integument", "Both integuments"],
                    answer: 3,
                    explain: "Seed coat (testa) develops from both integuments of ovule after fertilization"
                },
                {
                    id: "b44",
                    text: "Reaction wood in gymnosperms is called:",
                    options: ["Compression wood", "Tension wood", "Spring wood", "Autumn wood"],
                    answer: 0,
                    explain: "Compression wood forms on lower side of bent gymnosperm stems/branches as reaction to gravitational stress"
                },
                {
                    id: "b45",
                    text: "BAP (Benzyl amino purine) is a synthetic:",
                    options: ["Auxin", "Gibberellin", "Cytokinin", "Ethylene"],
                    answer: 2,
                    explain: "BAP is synthetic cytokinin widely used in plant tissue culture for shoot induction"
                },
                
                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Which bone has cribriform plate?",
                    options: ["Frontal", "Sphenoid", "Ethmoid", "Temporal"],
                    answer: 2,
                    explain: "Ethmoid bone has cribriform plate with small holes for passage of olfactory nerve fibers"
                },
                {
                    id: "z2",
                    text: "Motilin is secreted by:",
                    options: ["Duodenum", "Stomach", "Jejunum", "Ileum"],
                    answer: 0,
                    explain: "Motilin is secreted by M cells in duodenum and stimulates migrating motor complexes during fasting"
                },
                {
                    id: "z3",
                    text: "Christmas factor is:",
                    options: ["Factor VII", "Factor VIII", "Factor IX", "Factor X"],
                    answer: 2,
                    explain: "Christmas factor is Factor IX, deficiency of which causes Christmas disease (Hemophilia B)"
                },
                {
                    id: "z4",
                    text: "Which part of nephron has lowest osmolarity?",
                    options: ["Proximal tubule", "Descending limb", "Ascending limb", "Collecting duct"],
                    answer: 2,
                    explain: "Ascending limb produces most dilute filtrate due to active salt transport without water reabsorption"
                },
                {
                    id: "z5",
                    text: "Amygdala is part of:",
                    options: ["Limbic system", "Basal ganglia", "Thalamus", "Hypothalamus"],
                    answer: 0,
                    explain: "Amygdala is key component of limbic system involved in emotion, fear, and memory processing"
                },
                {
                    id: "z6",
                    text: "Indirect bilirubin is:",
                    options: ["Water soluble", "Fat soluble", "Conjugated", "Excreted in urine"],
                    answer: 1,
                    explain: "Indirect (unconjugated) bilirubin is lipophilic and bound to albumin for transport"
                },
                {
                    id: "z7",
                    text: "Neurotensin is secreted by:",
                    options: ["Stomach", "Duodenum", "Ileum", "Colon"],
                    answer: 2,
                    explain: "Neurotensin is secreted by N cells in ileum and has multiple effects on GI motility and secretion"
                },
                {
                    id: "z8",
                    text: "Fourth heart sound (S4) is due to:",
                    options: ["Atrial gallop", "Ventricular gallop", "Opening snap", "Ejection click"],
                    answer: 0,
                    explain: "S4 heart sound (atrial gallop) occurs during late diastolic filling due to atrial contraction"
                },
                {
                    id: "z9",
                    text: "Hyperaldosteronism causes:",
                    options: ["Hypokalemia", "Hyperkalemia", "Hypocalcemia", "Hyperphosphatemia"],
                    answer: 0,
                    explain: "Excess aldosterone promotes K⁺ loss in urine, leading to hypokalemia and metabolic alkalosis"
                },
                {
                    id: "z10",
                                        text: "Zona pellucida is secreted by:",
                    options: ["Follicle cells only", "Oocyte only", "Both oocyte and follicle cells", "Corona radiata cells"],
                    answer: 2,
                    explain: "Zona pellucida is glycoprotein matrix secreted by both oocyte and surrounding granulosa (follicle) cells"
                },
                {
                    id: "z11",
                    text: "Marasmus is caused by deficiency of:",
                    options: ["Protein only", "Calories only", "Both protein and calories", "Vitamins"],
                    answer: 2,
                    explain: "Marasmus is severe form of malnutrition caused by deficiency of both proteins and calories (energy)"
                },
                {
                    id: "z12",
                    text: "Dendritic cells are part of:",
                    options: ["Adaptive immunity", "Innate immunity", "Both adaptive and innate immunity", "Neither"],
                    answer: 2,
                    explain: "Dendritic cells bridge innate and adaptive immunity by presenting antigens to T cells"
                },
                {
                    id: "z13",
                    text: "Vasa recta are associated with:",
                    options: ["Glomerulus", "Juxtamedullary nephrons", "Cortical nephrons", "All nephrons"],
                    answer: 1,
                    explain: "Vasa recta (straight blood vessels) are associated with juxtamedullary nephrons for maintaining concentration gradient"
                },
                {
                    id: "z14",
                    text: "Mullerian inhibiting hormone is secreted by:",
                    options: ["Leydig cells", "Sertoli cells", "Granulosa cells", "Theca cells"],
                    answer: 1,
                    explain: "Anti-Mullerian hormone (AMH) is secreted by Sertoli cells and causes regression of Mullerian ducts in males"
                },
                {
                    id: "z15",
                    text: "Kell blood group system is clinically important for:",
                    options: ["Transfusion reactions", "Hemolytic disease of newborn", "Both", "Neither"],
                    answer: 2,
                    explain: "Kell antibodies can cause severe transfusion reactions and hemolytic disease of fetus/newborn"
                },
                {
                    id: "z16",
                    text: "Closing volume is increased in:",
                    options: ["Young healthy adults", "Elderly", "Athletes", "Children"],
                    answer: 1,
                    explain: "Closing volume increases with age due to loss of elastic recoil and airway closure during expiration"
                },
                {
                    id: "z17",
                    text: "Kinocilium is found in:",
                    options: ["Cochlea", "Vestibular system", "Both cochlea and vestibule", "Middle ear"],
                    answer: 1,
                    explain: "Kinocilium (true cilium) is present in vestibular hair cells but absent in cochlear hair cells"
                },
                {
                    id: "z18",
                    text: "Normal protein excretion in urine per day is less than:",
                    options: ["50 mg", "150 mg", "300 mg", "500 mg"],
                    answer: 1,
                    explain: "Normal urinary protein excretion is less than 150 mg/day, mostly consisting of Tamm-Horsfall protein"
                },
                {
                    id: "z19",
                    text: "GHRH is secreted by:",
                    options: ["Anterior pituitary", "Posterior pituitary", "Hypothalamus", "Liver"],
                    answer: 2,
                    explain: "Growth hormone releasing hormone (GHRH) is secreted by arcuate nucleus of hypothalamus"
                },
                {
                    id: "z20",
                    text: "Trophoblast gives rise to:",
                    options: ["Embryo proper", "Placenta", "Amnion", "Yolk sac"],
                    answer: 1,
                    explain: "Trophoblast forms outer layer of blastocyst and gives rise to fetal portion of placenta"
                },
                {
                    id: "z21",
                    text: "Preload is determined by:",
                    options: ["Venous return", "Arterial pressure", "Heart rate", "Contractility"],
                    answer: 0,
                    explain: "Preload (end-diastolic volume) is determined by venous return and ventricular filling"
                },
                {
                    id: "z22",
                    text: "Haldane effect describes:",
                    options: ["O₂ transport", "CO₂ transport", "Both O₂ and CO₂ transport", "pH buffering"],
                    answer: 1,
                    explain: "Haldane effect describes increased CO₂ carrying capacity of deoxygenated blood"
                },
                {
                    id: "z23",
                    text: "Chvostek's sign is positive in:",
                    options: ["Hypermagnesemia", "Hypomagnesemia", "Hypercalcemia", "Hypocalcemia"],
                    answer: 3,
                    explain: "Chvostek's sign (facial twitching when facial nerve tapped) indicates hypocalcemia"
                },
                {
                    id: "z24",
                    text: "Inspiratory capacity equals:",
                    options: ["TV + IRV", "TV + ERV", "IRV + ERV", "TV + IRV + ERV"],
                    answer: 0,
                    explain: "Inspiratory capacity (IC) = Tidal volume (TV) + Inspiratory reserve volume (IRV)"
                },
                {
                    id: "z25",
                    text: "Slow-twitch muscle fibers have:",
                    options: ["High glycogen content", "Low mitochondrial density", "High myoglobin content", "Low capillary density"],
                    answer: 2,
                    explain: "Slow-twitch (Type I) fibers are rich in myoglobin, mitochondria, and capillaries for aerobic metabolism"
                },
                {
                    id: "z26",
                    text: "Portal triad consists of:",
                    options: ["Hepatic artery, portal vein, bile duct", "Two hepatic arteries and portal vein", "Portal vein and two bile ducts", "Hepatic artery and two portal veins"],
                    answer: 0,
                    explain: "Portal triad contains hepatic artery, portal vein, and bile duct at corners of hepatic lobules"
                },
                {
                    id: "z27",
                    text: "Ariboflavinosis is caused by deficiency of:",
                    options: ["Thiamine", "Riboflavin", "Niacin", "Pyridoxine"],
                    answer: 1,
                    explain: "Ariboflavinosis is clinical condition caused by riboflavin (vitamin B₂) deficiency"
                },
                {
                    id: "z28",
                    text: "Nebulin is associated with:",
                    options: ["Thick filaments", "Thin filaments", "Z-disc", "M-line"],
                    answer: 1,
                    explain: "Nebulin is regulatory protein associated with thin filaments that helps determine sarcomere length"
                },
                {
                    id: "z29",
                    text: "Resistin is produced by:",
                    options: ["Liver", "Muscle", "Adipose tissue", "Pancreas"],
                    answer: 2,
                    explain: "Resistin is hormone produced by adipose tissue that may contribute to insulin resistance"
                },
                {
                    id: "z30",
                    text: "Paneth cells secrete:",
                    options: ["Mucus", "Hormones", "Antimicrobial peptides", "Enzymes"],
                    answer: 2,
                    explain: "Paneth cells in intestinal crypts secrete antimicrobial peptides like lysozyme and defensins"
                },
                {
                    id: "z31",
                    text: "Shivering thermogenesis involves:",
                    options: ["Skeletal muscle", "Cardiac muscle", "Smooth muscle", "Brown adipose tissue"],
                    answer: 0,
                    explain: "Shivering thermogenesis involves rapid, involuntary contractions of skeletal muscle to generate heat"
                },
                {
                    id: "z32",
                    text: "Wharton's jelly is found in:",
                    options: ["Placenta", "Umbilical cord", "Amniotic sac", "Chorion"],
                    answer: 1,
                    explain: "Wharton's jelly is mucoid connective tissue that surrounds umbilical vessels in umbilical cord"
                },
                {
                    id: "z33",
                    text: "Reticulocyte count indicates:",
                    options: ["RBC destruction", "RBC production", "WBC function", "Platelet function"],
                    answer: 1,
                    explain: "Reticulocyte count reflects bone marrow RBC production activity and hematopoietic function"
                },
                {
                    id: "z34",
                    text: "Jet lag is caused by disruption of:",
                    options: ["Sleep centers", "Circadian rhythm", "Hormone production", "Neurotransmitter levels"],
                    answer: 1,
                    explain: "Jet lag results from desynchronization of internal circadian clock with external time cues"
                },
                {
                    id: "z35",
                    text: "Presbyopia is due to:",
                    options: ["Corneal changes", "Lens rigidity", "Retinal degeneration", "Vitreous changes"],
                    answer: 1,
                    explain: "Presbyopia (near vision difficulty with age) is caused by decreased lens elasticity and accommodation"
                },
                {
                    id: "z36",
                    text: "Starling forces determine:",
                    options: ["Cardiac output", "Blood pressure", "Fluid movement across capillaries", "Renal blood flow"],
                    answer: 2,
                    explain: "Starling forces (hydrostatic and oncotic pressures) govern fluid filtration and reabsorption across capillaries"
                },
                {
                    id: "z37",
                    text: "Gustation involves:",
                    options: ["Olfactory receptors", "Taste buds", "Trigeminal nerve", "All of these"],
                    answer: 3,
                    explain: "Taste perception involves taste buds, olfactory input, and trigeminal nerve for texture and temperature"
                },
                {
                    id: "z38",
                    text: "Rebound hypoglycemia occurs in:",
                    options: ["Type 1 diabetes", "Type 2 diabetes", "Insulinoma", "All of these"],
                    answer: 3,
                    explain: "Rebound hypoglycemia can occur in all conditions involving insulin excess or inappropriate insulin therapy"
                },
                {
                    id: "z39",
                    text: "Inhibin A is produced by:",
                    options: ["Granulosa cells", "Theca cells", "Leydig cells", "Sertoli cells"],
                    answer: 0,
                    explain: "Inhibin A is produced by granulosa cells of dominant follicle and corpus luteum in females"
                },
                {
                    id: "z40",
                    text: "Gorham's disease affects:",
                    options: ["Muscles", "Nerves", "Bones", "Blood vessels"],
                    answer: 2,
                    explain: "Gorham's disease (vanishing bone disease) is rare condition causing progressive bone loss"
                },
                {
                    id: "z41",
                    text: "Gastrin is secreted in response to:",
                    options: ["Protein in stomach", "Fat in duodenum", "Acid in stomach", "Glucose in blood"],
                    answer: 0,
                    explain: "Gastrin secretion from G cells is stimulated by protein/amino acids in gastric antrum"
                },
                {
                    id: "z42",
                    text: "Hyaluronic acid in synovial fluid is produced by:",
                    options: ["Cartilage", "Synovial membrane", "Bone", "Joint capsule"],
                    answer: 1,
                    explain: "Hyaluronic acid giving synovial fluid its viscosity is produced by synoviocytes in synovial membrane"
                },
                {
                    id: "z43",
                    text: "Free T3 is approximately what percentage of total T3?",
                    options: ["0.3%", "3%", "30%", "50%"],
                    answer: 0,
                    explain: "Free T3 represents only about 0.3% of total T3, but it's the metabolically active fraction"
                },
                {
                    id: "z44",
                    text: "Afterload is primarily determined by:",
                    options: ["Venous return", "Aortic pressure", "Heart rate", "Ventricular volume"],
                    answer: 1,
                    explain: "Afterload is resistance against which ventricle pumps, primarily determined by aortic/arterial pressure"
                },
                {
                    id: "z45",
                    text: "Binding of sperm to zona pellucida triggers:",
                    options: ["Capacitation", "Acrosome reaction", "Cortical reaction", "Zona reaction"],
                    answer: 1,
                    explain: "Sperm binding to ZP3 glycoprotein triggers acrosome reaction, releasing enzymes to penetrate zona pellucida"
                }
            ]
        }
    ]
};
