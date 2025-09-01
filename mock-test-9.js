// mock-test-9.js - NEET Mock Test 9 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// All NEW questions - no repeats from previous mock tests

window.MOCK_TEST_9 = {
    id: "neet-009",
    title: "Full Syllabus Mock 9", 
    level: "medium",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A particle moves in a circle of radius R with constant angular velocity ω. The magnitude of change in velocity when particle moves through angle π/3 is:",
                    options: ["Rω", "Rω√3", "2Rω sin(π/6)", "Rω√3/2"],
                    answer: 0,
                    explain: "Change in velocity = |v₂ - v₁| = 2v sin(θ/2) = 2Rω sin(π/6) = 2Rω × 1/2 = Rω"
                },
                {
                    id: "p2",
                    text: "Two masses m₁ and m₂ connected by a light string pass over a frictionless pulley. If m₁ > m₂, the tension in string is:",
                    options: ["m₁g", "m₂g", "2m₁m₂g/(m₁+m₂)", "(m₁-m₂)g"],
                    answer: 2,
                    explain: "For Atwood machine: T = 2m₁m₂g/(m₁+m₂). This accounts for both masses and their acceleration."
                },
                {
                    id: "p3",
                    text: "A parallel plate capacitor has plate area A and separation d. When a dielectric slab of thickness t and dielectric constant K is inserted, capacitance becomes:",
                    options: ["ε₀A/[d-t+t/K]", "ε₀A/[d-t(1-1/K)]", "Kε₀A/d", "ε₀A(K+1)/d"],
                    answer: 1,
                    explain: "With partial dielectric insertion: C = ε₀A/[d-t+t/K] = ε₀A/[d-t(1-1/K)]"
                },
                {
                    id: "p4",
                    text: "In Young's double slit experiment, if one slit is covered with transparent sheet of thickness t and refractive index μ, fringe shift is:",
                    options: ["(μ-1)t/λ", "(μ-1)tD/λd", "(μ+1)t/λ", "μt/λ"],
                    answer: 1,
                    explain: "Path difference introduced = (μ-1)t. Number of fringes shifted = (μ-1)t/λ. Actual shift on screen = (μ-1)tD/λd"
                },
                {
                    id: "p5",
                    text: "In hydrogen atom, electron jumps from n = 5 to n = 1. How many spectral lines are possible?",
                    options: ["4", "6", "8", "10"],
                    answer: 3,
                    explain: "Number of spectral lines = n(n-1)/2 where n is number of levels. From n=5 to n=1: transitions possible are 5→4, 5→3, 5→2, 5→1, 4→3, 4→2, 4→1, 3→2, 3→1, 2→1. Total = 10 lines"
                },
                {
                    id: "p6",
                    text: "A wire of resistance 10 Ω is drawn into wire of double length. Its new resistance becomes:",
                    options: ["20 Ω", "40 Ω", "5 Ω", "2.5 Ω"],
                    answer: 1,
                    explain: "When length doubles at constant volume, area becomes A/2. New resistance R' = ρ(2L)/(A/2) = 4ρL/A = 4R = 40 Ω"
                },
                {
                    id: "p7",
                    text: "A solid sphere rolls down an inclined plane of height h. Its velocity at bottom is:",
                    options: ["√(2gh)", "√(10gh/7)", "√(5gh/3)", "√(gh)"],
                    answer: 1,
                    explain: "For rolling sphere: mgh = ½mv² + ½Iω² = ½mv² + ½(2/5)mr²(v/r)² = 7/10 mv². Therefore v = √(10gh/7)"
                },
                {
                    id: "p8",
                    text: "In LCR series circuit, if frequency is doubled, reactance of inductor becomes:",
                    options: ["Half", "Double", "Four times", "One-fourth"],
                    answer: 1,
                    explain: "Inductive reactance XL = ωL = 2πfL. When frequency doubles, XL becomes 2XL (doubles)"
                },
                {
                    id: "p9",
                    text: "A gas expands adiabatically to twice its volume. If initial pressure is P₀, final pressure is:",
                    options: ["P₀/2", "P₀/2^γ", "P₀/4", "P₀√2"],
                    answer: 1,
                    explain: "For adiabatic process: PV^γ = constant. P₁V₁^γ = P₂V₂^γ. P₀ × V₀^γ = P₂(2V₀)^γ. P₂ = P₀/2^γ"
                },
                {
                    id: "p10",
                    text: "Two parallel wires carry currents I₁ and I₂ in same direction separated by distance d. Force per unit length is:",
                    options: ["μ₀I₁I₂/2πd (attractive)", "μ₀I₁I₂/2πd (repulsive)", "μ₀I₁I₂/4πd (attractive)", "μ₀I₁I₂d/2π"],
                    answer: 0,
                    explain: "Force per unit length between parallel current-carrying wires: F/L = μ₀I₁I₂/2πd. Same direction currents attract each other"
                },
                {
                    id: "p11",
                    text: "A particle in SHM has displacement x = 5 sin(2πt + π/4). Its amplitude and initial phase are:",
                    options: ["5, π/4", "5, π/2", "10, π/4", "2.5, π/4"],
                    answer: 0,
                    explain: "Comparing with x = A sin(ωt + φ): Amplitude A = 5, initial phase φ = π/4"
                },
                {
                    id: "p12",
                    text: "In photoelectric effect, if intensity of light is doubled keeping frequency constant, stopping potential:",
                    options: ["Doubles", "Becomes half", "Remains same", "Becomes four times"],
                    answer: 2,
                    explain: "Stopping potential depends only on frequency of light, not intensity. eV₀ = hf - φ. Intensity affects only number of photoelectrons"
                },
                {
                    id: "p13",
                    text: "Three resistors 2Ω, 3Ω and 6Ω are connected in parallel. Their equivalent resistance is:",
                    options: ["1Ω", "2Ω", "11Ω", "0.5Ω"],
                    answer: 0,
                    explain: "1/R = 1/2 + 1/3 + 1/6 = 3/6 + 2/6 + 1/6 = 6/6 = 1. Therefore R = 1Ω"
                },
                {
                    id: "p14",
                    text: "A projectile has range R on horizontal ground. Maximum height reached is:",
                    options: ["R/2", "R/4", "R tan²θ/4", "R/4 tan θ"],
                    answer: 2,
                    explain: "Range R = u²sin2θ/g, Max height H = u²sin²θ/2g. H/R = sin²θ/2sin2θ = sinθ/4cosθ = tanθ/4. But H = R tan²θ/4 for given range"
                },
                {
                    id: "p15",
                    text: "Current in a wire creates magnetic field B at center of circular loop. If wire is bent into square loop of same perimeter, field at center becomes:",
                    options: ["B", "B√2/π", "2√2B/π", "πB/2√2"],
                    answer: 2,
                    explain: "For circular loop: B₁ = μ₀I/2R. For square: B₂ = 4 × μ₀I√2/4πa = μ₀I√2/πa. Since perimeters equal: 2πR = 4a, R = 2a/π. Substituting: B₂/B₁ = 2√2/π"
                },
                {
                    id: "p16",
                    text: "A stone is dropped from tower of height H. In last second of flight, it covers distance H/4. Height of tower is:",
                    options: ["80 m", "125 m", "100 m", "144 m"],
                    answer: 1,
                    explain: "Let total time be t. Distance in last second = H - distance in (t-1) seconds = H - ½g(t-1)² = H/4. Also H = ½gt². Solving: t = 5s, H = 125m"
                },
                {
                    id: "p17",
                    text: "In series RLC circuit, at resonance the phase difference between voltage and current is:",
                    options: ["0°", "90°", "180°", "45°"],
                    answer: 0,
                    explain: "At resonance XL = XC, so net reactance is zero. Circuit behaves as pure resistive, hence phase difference = 0°"
                },
                {
                    id: "p18",
                    text: "Momentum of photon of wavelength λ is:",
                    options: ["hλ", "h/λ", "hc/λ", "λ/h"],
                    answer: 1,
                    explain: "Photon momentum p = E/c = hf/c = h/λ (using c = fλ)"
                },
                {
                    id: "p19",
                    text: "A spring-mass system oscillates with frequency f. If both mass and spring constant are doubled, new frequency is:",
                    options: ["f", "2f", "f/2", "f√2"],
                    answer: 0,
                    explain: "f = (1/2π)√(k/m). When both k and m are doubled: f' = (1/2π)√(2k/2m) = (1/2π)√(k/m) = f"
                },
                {
                    id: "p20",
                    text: "Electric field inside a conducting sphere of radius R carrying charge Q is:",
                    options: ["kQ/R²", "0", "kQ/r²", "kQ/R³"],
                    answer: 1,
                    explain: "Inside a conductor in electrostatic equilibrium, electric field is always zero due to redistribution of charges on surface"
                },
                {
                    id: "p21",
                    text: "A concave mirror of focal length 15 cm forms virtual image of magnification +3. Object distance is:",
                    options: ["5 cm", "10 cm", "12 cm", "20 cm"],
                    answer: 1,
                    explain: "For virtual image: m = +3 = -v/u, so v = -3u. Using mirror equation: 1/f = 1/u + 1/v = 1/u - 1/3u = 2/3u. 1/15 = 2/3u, so u = 10 cm"
                },
                {
                    id: "p22",
                    text: "Two coherent sources with phase difference π/2 interfere. If individual amplitudes are a, resultant amplitude is:",
                    options: ["2a", "a√2", "a", "0"],
                    answer: 1,
                    explain: "Resultant amplitude A = √(a₁² + a₂² + 2a₁a₂cos φ) = √(a² + a² + 2a²cos(π/2)) = √(2a²) = a√2"
                },
                {
                    id: "p23",
                    text: "In isothermal expansion of ideal gas from volume V to 2V, work done is:",
                    options: ["PV", "PV ln 2", "2PV", "PV/2"],
                    answer: 1,
                    explain: "For isothermal process: W = nRT ln(V₂/V₁) = PV ln(2V/V) = PV ln 2"
                },
                {
                    id: "p24",
                    text: "Half-life of radioactive element is 10 days. After 30 days, fraction of sample decayed is:",
                    options: ["1/8", "7/8", "1/4", "3/4"],
                    answer: 1,
                    explain: "After 3 half-lives (30 days): remaining fraction = (1/2)³ = 1/8. Decayed fraction = 1 - 1/8 = 7/8"
                },
                {
                    id: "p25",
                    text: "Energy stored in inductor carrying current I is:",
                    options: ["LI", "½LI²", "LI²", "2LI²"],
                    answer: 1,
                    explain: "Energy stored in inductor: U = ½LI² where L is inductance and I is current"
                },
                {
                    id: "p26",
                    text: "A car travels first half distance at speed v₁ and second half at speed v₂. Average speed is:",
                    options: ["(v₁+v₂)/2", "2v₁v₂/(v₁+v₂)", "√(v₁v₂)", "v₁v₂/(v₁+v₂)"],
                    answer: 1,
                    explain: "Average speed = Total distance/Total time = 2s/(s/v₁ + s/v₂) = 2v₁v₂/(v₁+v₂)"
                },
                {
                    id: "p27",
                    text: "For total internal reflection, light must travel from:",
                    options: ["Rarer to denser medium", "Denser to rarer medium", "Same density media", "Any direction"],
                    answer: 1,
                    explain: "Total internal reflection occurs only when light travels from optically denser to rarer medium at angles greater than critical angle"
                },
                {
                    id: "p28",
                    text: "In stationary wave, distance between consecutive nodes is:",
                    options: ["λ/4", "λ/2", "λ", "2λ"],
                    answer: 1,
                    explain: "In stationary wave, nodes are separated by λ/2 where λ is wavelength of constituent waves"
                },
                {
                    id: "p29",
                    text: "Efficiency of heat engine working between 600K and 300K is:",
                    options: ["25%", "50%", "75%", "100%"],
                    answer: 1,
                    explain: "Carnot efficiency η = 1 - T₂/T₁ = 1 - 300/600 = 1 - 0.5 = 0.5 = 50%"
                },
                {
                    id: "p30",
                    text: "When a charged particle enters magnetic field perpendicularly, its kinetic energy:",
                    options: ["Increases", "Decreases", "Remains constant", "Becomes zero"],
                    answer: 2,
                    explain: "Magnetic force is always perpendicular to velocity, does no work. Hence kinetic energy remains constant"
                },
                {
                    id: "p31",
                    text: "Sound wave travels from air to water. Its frequency:",
                    options: ["Increases", "Decreases", "Remains same", "Becomes zero"],
                    answer: 2,
                    explain: "Frequency is characteristic of source and remains unchanged when wave enters different medium"
                },
                {
                    id: "p32",
                    text: "Moment of inertia of thin ring about diameter is:",
                    options: ["MR²", "MR²/2", "MR²/4", "2MR²"],
                    answer: 1,
                    explain: "For thin ring about diameter (perpendicular axis theorem): I = MR²/2"
                },
                {
                    id: "p33",
                    text: "Power consumed by resistor in AC circuit with RMS current I is:",
                    options: ["I²R", "2I²R", "I²R/2", "I²R/√2"],
                    answer: 0,
                    explain: "Average power in AC circuit: P = I²rmsR where Irms is RMS current"
                },
                {
                    id: "p34",
                    text: "A lens of power -2.5 D has focal length:",
                    options: ["-0.4 m", "+0.4 m", "-2.5 m", "+2.5 m"],
                    answer: 0,
                    explain: "Power P = 1/f. For P = -2.5 D: f = 1/(-2.5) = -0.4 m (diverging lens)"
                },
                {
                    id: "p35",
                    text: "Binding energy per nucleon is minimum for:",
                    options: ["Light nuclei", "Medium nuclei", "Heavy nuclei", "All same"],
                    answer: 0,
                    explain: "Light nuclei have lower binding energy per nucleon due to surface effects and Coulomb repulsion"
                },
                {
                    id: "p36",
                    text: "A particle executes SHM with period 4s. Time to go from mean position to half amplitude is:",
                    options: ["1/3 s", "2/3 s", "1 s", "4/3 s"],
                    answer: 1,
                    explain: "x = A sin(ωt). For x = A/2: sin(ωt) = 1/2, ωt = π/6. t = π/6ω = π/6 × 2π/4 = T/12 = 4/12 = 1/3 s. Wait, let me recalculate: t = (π/6)/(2π/4) = (π/6) × (4/2π) = 4/12 = 1/3 s. Actually for SHM x = A sin(2πt/T), if x = A/2, then sin(2πt/T) = 1/2, so 2πt/T = π/6, t = T/12 = 4/12 = 1/3 s"
                },
                {
                    id: "p37",
                    text: "Gauss's law for magnetism states that:",
                    options: ["∮B⋅dA = 0", "∮B⋅dA = μ₀I", "∮B⋅dl = 0", "∮B⋅dl = μ₀I"],
                    answer: 0,
                    explain: "Gauss's law for magnetism: ∮B⋅dA = 0, indicating no magnetic monopoles exist"
                },
                {
                    id: "p38",
                    text: "Step-up transformer has turn ratio 1:10. If primary current is 5A, secondary current is:",
                    options: ["0.5 A", "5 A", "50 A", "500 A"],
                    answer: 0,
                    explain: "For ideal transformer: I₁/I₂ = N₂/N₁. Given N₁:N₂ = 1:10, so I₂ = I₁ × N₁/N₂ = 5 × 1/10 = 0.5 A"
                },
                {
                    id: "p39",
                    text: "Time constant of LR circuit is 0.1 s. After 0.2 s, current reaches what fraction of maximum?",
                    options: ["(1-e⁻²)", "(1-e⁻¹)", "e⁻²", "e⁻¹"],
                    answer: 0,
                    explain: "Current in LR circuit: I(t) = I₀(1-e⁻ᵗ/τ). At t = 0.2s, τ = 0.1s: I = I₀(1-e⁻²)"
                },
                {
                    id: "p40",
                    text: "Work done in moving charge q from potential V₁ to V₂ is:",
                    options: ["q(V₂-V₁)", "q(V₁-V₂)", "qV₁V₂", "q/(V₂-V₁)"],
                    answer: 0,
                    explain: "Work done = Change in potential energy = q(V₂-V₁)"
                },
                {
                    id: "p41",
                    text: "Gas bubble at bottom of lake has volume V₀. At surface (pressure halves), its volume becomes:",
                    options: ["V₀/2", "V₀", "2V₀", "4V₀"],
                    answer: 2,
                    explain: "At constant temperature (isothermal): P₁V₁ = P₂V₂. If P₂ = P₁/2, then V₂ = 2V₁ = 2V₀"
                },
                {
                    id: "p42",
                    text: "Impedance of pure capacitor in AC circuit is:",
                    options: ["ωC", "1/ωC", "ωL", "R"],
                    answer: 1,
                    explain: "Capacitive reactance XC = 1/ωC acts as impedance for pure capacitor"
                },
                {
                    id: "p43",
                    text: "In interference, two waves of amplitudes 3 and 4 units interfere. Maximum intensity ratio to minimum is:",
                    options: ["7:1", "49:1", "16:9", "25:1"],
                    answer: 1,
                    explain: "Imax = (a₁+a₂)² = (3+4)² = 49, Imin = (a₁-a₂)² = (4-3)² = 1. Ratio = 49:1"
                },
                {
                    id: "p44",
                    text: "Escape velocity from planet of mass M and radius R is:",
                    options: ["√(GM/R)", "√(2GM/R)", "√(GM/2R)", "√(2GM/R²)"],
                    answer: 1,
                    explain: "Escape velocity ve = √(2GM/R) where G is gravitational constant"
                },
                {
                    id: "p45",
                    text: "In common emitter amplifier, if base current changes by 10 μA and collector current by 1 mA, current gain β is:",
                    options: ["10", "50", "100", "1000"],
                    answer: 2,
                    explain: "Current gain β = ΔIC/ΔIB = 1 mA/10 μA = 1000 μA/10 μA = 100"
                }
            ]
        },
        {
            name: "Chemistry",
            questions: [
                {
                    id: "c1",
                    text: "Which has maximum polarizing power?",
                    options: ["Li⁺", "Na⁺", "Mg²⁺", "Al³⁺"],
                    answer: 3,
                    explain: "Polarizing power ∝ charge/size². Al³⁺ has highest charge (+3) and small size, giving maximum polarizing power"
                },
                {
                    id: "c2",
                    text: "The hybridization of I in IF₇ is:",
                    options: ["sp³d²", "sp³d³", "sp³d⁴", "sp³d⁵"],
                    answer: 1,
                    explain: "IF₇ has 7 bonding pairs around I, requiring sp³d³ hybridization with pentagonal bipyramidal geometry"
                },
                {
                    id: "c3",
                    text: "Which carbocation is most stable?",
                    options: ["CH₃⁺", "CH₃CH₂⁺", "(CH₃)₂CH⁺", "C₆H₅CH₂⁺"],
                    answer: 3,
                    explain: "Benzyl carbocation C₆H₅CH₂⁺ is stabilized by resonance with benzene ring, making it most stable"
                },
                {
                    id: "c4",
                    text: "Oxidation state of Mn in K₂MnO₄ is:",
                    options: ["+4", "+5", "+6", "+7"],
                    answer: 2,
                    explain: "In K₂MnO₄: 2(+1) + Mn + 4(-2) = 0, solving gives Mn = +6"
                },
                {
                    id: "c5",
                    text: "Number of chiral carbons in glucose is:",
                    options: ["2", "3", "4", "5"],
                    answer: 2,
                    explain: "Glucose has 4 chiral carbons at C2, C3, C4, and C5 positions in its chain form"
                },
                {
                    id: "c6",
                    text: "Strongest base among the following is:",
                    options: ["F⁻", "Cl⁻", "Br⁻", "I⁻"],
                    answer: 0,
                    explain: "Base strength is inverse of conjugate acid strength. HF is weakest acid among HX, so F⁻ is strongest base"
                },
                {
                    id: "c7",
                    text: "Electronic configuration of Cr³⁺ is:",
                    options: ["[Ar] 3d³", "[Ar] 3d⁴ 4s⁻¹", "[Ar] 3d⁵", "[Ar] 4s¹ 3d²"],
                    answer: 0,
                    explain: "Cr (24): [Ar] 3d⁵ 4s¹. Cr³⁺ loses 4s¹ and two 3d electrons: [Ar] 3d³"
                },
                {
                    id: "c8",
                    text: "Geometry of SF₆ is:",
                    options: ["Trigonal bipyramidal", "Octahedral", "Square pyramidal", "Pentagonal bipyramidal"],
                    answer: 1,
                    explain: "SF₆ has 6 bonding pairs and no lone pairs around S, giving regular octahedral geometry"
                },
                {
                    id: "c9",
                    text: "Which cannot act as Lewis acid?",
                    options: ["BF₃", "AlCl₃", "NH₃", "FeCl₃"],
                    answer: 2,
                    explain: "NH₃ has lone pair and acts as Lewis base (electron donor), not Lewis acid (electron acceptor)"
                },
                {
                    id: "c10",
                    text: "Entropy is maximum for:",
                    options: ["Diamond", "Graphite", "Ice", "Water vapor"],
                    answer: 3,
                    explain: "Water vapor (gas) has maximum disorder and molecular motion, hence maximum entropy"
                },
                {
                    id: "c11",
                    text: "Weakest reducing agent is:",
                    options: ["Li", "Na", "K", "Au"],
                    answer: 3,
                    explain: "Au has most positive reduction potential (+1.50 V), making it weakest reducing agent (strongest oxidizing agent)"
                },
                {
                    id: "c12",
                    text: "Number of isomers of [Co(NH₃)₄Cl₂]⁺ is:",
                    options: ["2", "3", "4", "5"],
                    answer: 0,
                    explain: "Octahedral complex [Co(NH₃)₄Cl₂]⁺ shows only geometrical isomerism: cis and trans forms (2 isomers)"
                },
                {
                    id: "c13",
                    text: "Which undergoes haloform reaction?",
                    options: ["CH₃CHO", "HCHO", "CH₃COCH₃", "C₂H₅CHO"],
                    answer: 2,
                    explain: "Compounds with CH₃CO- group undergo haloform reaction. CH₃COCH₃ (acetone) has this group"
                },
                {
                    id: "c14",
                    text: "Correct order of hydration enthalpy is:",
                    options: ["Li⁺ < Na⁺ < K⁺", "K⁺ < Na⁺ < Li⁺", "Na⁺ < Li⁺ < K⁺", "Li⁺ = Na⁺ = K⁺"],
                    answer: 1,
                    explain: "Hydration enthalpy is inversely proportional to ionic size. Li⁺ is smallest, so order: K⁺ < Na⁺ < Li⁺"
                },
                {
                    id: "c15",
                    text: "pH of 0.1 M CH₃COOH (Ka = 1.8 × 10⁻⁵) is approximately:",
                    options: ["1.0", "2.9", "4.7", "7.0"],
                    answer: 1,
                    explain: "For weak acid: [H⁺] = √(Ka × C) = √(1.8×10⁻⁵ × 0.1) = √(1.8×10⁻⁶) ≈ 1.34×10⁻³. pH = -log(1.34×10⁻³) ≈ 2.9"
                },
                {
                    id: "c16",
                    text: "Magnetic moment of [Co(NH₃)₆]³⁺ is:",
                    options: ["0 BM", "3.87 BM", "4.90 BM", "5.92 BM"],
                    answer: 0,
                    explain: "Co³⁺ (d⁶) with strong field NH₃ ligands forms low-spin complex with all electrons paired, μ = 0 BM"
                },
                {
                    id: "c17",
                    text: "Most ionic compound is:",
                    options: ["LiCl", "NaCl", "KCl", "CsCl"],
                    answer: 3,
                    explain: "Ionic character increases with size difference between cation and anion. CsCl has maximum size difference, hence most ionic"
                },
                {
                    id: "c18",
                    text: "Rate constant doubles when temperature increases from 300K to 310K. Activation energy is:",
                    options: ["53.6 kJ/mol", "26.8 kJ/mol", "107.2 kJ/mol", "13.4 kJ/mol"],
                    answer: 0,
                    explain: "Using Arrhenius equation: ln(k₂/k₁) = Ea/R(1/T₁ - 1/T₂). ln(2) = Ea/8.314(1/300 - 1/310). Solving: Ea ≈ 53.6 kJ/mol"
                },
                {
                    id: "c19",
                    text: "Strongest hydrogen bonding is in:",
                    options: ["HF", "HCl", "HBr", "HI"],
                    answer: 0,
                    explain: "HF forms strongest hydrogen bonds due to highest electronegativity of F and small size allowing close approach"
                },
                {
                    id: "c20",
                    text: "Aromaticity requires:",
                    options: ["4n electrons", "4n+2 electrons", "Even number of electrons", "Any number of electrons"],
                    answer: 1,
                    explain: "Hückel's rule states aromatic compounds must have 4n+2 π electrons where n = 0, 1, 2..."
                },
                {
                    id: "c21",
                    text: "Which has highest dipole moment?",
                    options: ["BeF₂", "BF₃", "NF₃", "CF₄"],
                    answer: 2,
                    explain: "NF₃ has pyramidal geometry with net dipole moment due to lone pair, while others are symmetrical"
                },
                {
                    id: "c22",
                    text: "EDTA is:",
                    options: ["Monodentate", "Bidentate", "Tridentate", "Hexadentate"],
                    answer: 3,
                    explain: "EDTA (ethylenediaminetetraacetic acid) can donate 6 electron pairs, making it hexadentate ligand"
                },
                {
                    id: "c23",
                    text: "Strongest acid is:",
                    options: ["HClO₄", "HClO₃", "HClO₂", "HClO"],
                    answer: 0,
                    explain: "Acid strength in oxoacids increases with number of oxygen atoms. HClO₄ has maximum O atoms, hence strongest"
                },
                {
                    id: "c24",
                    text: "Bond order in N₂⁺ is:",
                    options: ["2", "2.5", "3", "3.5"],
                    answer: 1,
                    explain: "N₂⁺ has 13 electrons. Bond order = (bonding - antibonding)/2 = (9 - 4)/2 = 2.5"
                },
                {
                    id: "c25",
                    text: "Electrophilic aromatic substitution is fastest in:",
                    options: ["Benzene", "Toluene", "Nitrobenzene", "Chlorobenzene"],
                    answer: 1,
                    explain: "Toluene has electron-donating methyl group that activates benzene ring toward electrophilic substitution"
                },
                {
                    id: "c26",
                    text: "Number of lone pairs on central atom in ClF₅ is:",
                    options: ["0", "1", "2", "3"],
                    answer: 1,
                    explain: "Cl in ClF₅ has 7 valence electrons, uses 5 for bonding with F atoms, leaving 2 electrons = 1 lone pair"
                },
                {
                    id: "c27",
                    text: "Which is diamagnetic?",
                    options: ["O₂", "O₂⁺", "O₂⁻", "O₂²⁻"],
                    answer: 3,
                    explain: "O₂²⁻ (peroxide ion) has 18 electrons with all paired in molecular orbitals, making it diamagnetic"
                },
                {
                    id: "c28",
                    text: "Oxidation number of S in Na₂S₄O₆ is:",
                    options: ["+2", "+2.5", "+4", "+6"],
                    answer: 1,
                    explain: "In tetrathionate ion S₄O₆²⁻: 4S + 6(-2) = -2, so 4S = +10, average oxidation state = +2.5"
                },
                {
                    id: "c29",
                    text: "Chelation effect increases:",
                    options: ["Entropy", "Enthalpy", "Both entropy and enthalpy", "Neither"],
                    answer: 0,
                    explain: "Chelation increases entropy due to release of more solvent molecules when multidentate ligand replaces monodentate ligands"
                },
                {
                    id: "c30",
                    text: "Packing efficiency of body-centered cubic is:",
                    options: ["52.4%", "68%", "74%", "90%"],
                    answer: 1,
                    explain: "BCC packing efficiency = (√3/8)π × 100% ≈ 68%"
                },
                {
                    id: "c31",
                    text: "Lowest boiling point is of:",
                    options: ["n-butane", "iso-butane", "sec-butanol", "tert-butanol"],
                    answer: 1,
                    explain: "iso-butane has most branched structure, least surface contact and van der Waals forces, hence lowest boiling point"
                },
                {
                    id: "c32",
                    text: "Unit of rate constant for zero-order reaction is:",
                    options: ["s⁻¹", "mol L⁻¹ s⁻¹", "mol⁻¹ L s⁻¹", "mol L⁻¹"],
                    answer: 1,
                    explain: "Zero order: rate = k[A]⁰ = k. Rate has units mol L⁻¹ s⁻¹, so k has same units"
                },
                {
                    id: "c33",
                    text: "Which molecule is planar?",
                    options: ["NH₃", "PCl₃", "BF₃", "NF₃"],
                    answer: 2,
                    explain: "BF₃ has trigonal planar geometry with all atoms in same plane, while others have pyramidal geometry"
                },
                {
                    id: "c34",
                    text: "Electron gain enthalpy is most negative for:",
                    options: ["F", "Cl", "Br", "I"],
                    answer: 1,
                    explain: "Cl has most negative electron gain enthalpy due to optimal size - not too small like F (electron repulsion) or too large like Br, I"
                },
                {
                    id: "c35",
                    text: "Benzoin condensation involves:",
                    options: ["Two benzaldehyde molecules", "Benzaldehyde and acetone", "Two acetone molecules", "Benzaldehyde and formaldehyde"],
                    answer: 0,
                    explain: "Benzoin condensation is self-condensation of two benzaldehyde molecules in presence of CN⁻ catalyst"
                },
                {
                    id: "c36",
                    text: "Bond angle in NH₃ is less than CH₄ due to:",
                    options: ["Lone pair repulsion", "Electronegativity difference", "Hybridization difference", "Size difference"],
                    answer: 0,
                    explain: "NH₃ has one lone pair that repels bonding pairs more than bonding pairs repel each other, reducing bond angle"
                },
                {
                    id: "c37",
                    text: "Smallest atom in third period is:",
                    options: ["Na", "Mg", "Al", "Cl"],
                    answer: 3,
                    explain: "Atomic size decreases across period due to increasing nuclear charge. Cl is rightmost, hence smallest"
                },
                {
                    id: "c38",
                    text: "mer-[Co(NH₃)₃Cl₃] has how many isomers?",
                    options: ["2", "3", "4", "5"],
                    answer: 0,
                    explain: "mer isomer has three identical ligands in meridional positions. It can show facial (fac) and meridional (mer) isomerism, but mer specifically refers to one form. The complex can show 2 geometrical isomers total"
                },
                {
                    id: "c39",
                    text: "Standard reduction potential is highest for:",
                    options: ["Li⁺/Li", "Na⁺/Na", "K⁺/K", "Au³⁺/Au"],
                    answer: 3,
                    explain: "Au³⁺/Au has highest (most positive) standard reduction potential (+1.50 V), indicating strong tendency to get reduced"
                },
                {
                    id: "c40",
                    text: "Carbon-carbon bond length is shortest in:",
                    options: ["Diamond", "Graphite", "Fullerene", "Graphene"],
                    answer: 1,
                    explain: "Graphite has partial double bond character due to delocalized π electrons, giving shorter C-C bonds than pure single bonds"
                },
                {
                    id: "c41",
                    text: "E-Z nomenclature is used for:",
                    options: ["Optical isomerism", "Geometrical isomerism", "Conformational isomerism", "Structural isomerism"],
                    answer: 1,
                    explain: "E-Z nomenclature (Entgegen-Zusammen) is systematic way to name geometrical isomers based on priority rules"
                },
                {
                    id: "c42",
                    text: "18-electron rule is obeyed by:",
                    options: ["[Ni(CO)₄]", "[Fe(CO)₅]", "[Cr(CO)₆]", "All of these"],
                    answer: 3,
                    explain: "All given carbonyls follow 18-electron rule: Ni(10) + 4×2 = 18, Fe(8) + 5×2 = 18, Cr(6) + 6×2 = 18"
                },
                {
                    id: "c43",
                    text: "Conjugate base of H₂SO₄ is:",
                    options: ["SO₄²⁻", "HSO₄⁻", "H₃SO₄⁺", "H₂SO₃"],
                    answer: 1,
                    explain: "Conjugate base is formed by removing one H⁺: H₂SO₄ → HSO₄⁻ + H⁺"
                },
                {
                    id: "c44",
                    text: "Triple bond has:",
                    options: ["1σ + 2π", "3σ", "2σ + 1π", "3π"],
                    answer: 0,
                    explain: "Triple bond consists of one sigma (σ) bond and two pi (π) bonds"
                },
                {
                    id: "c45",
                    text: "Bond dissociation energy is highest for:",
                    options: ["H-H", "F-F", "Cl-Cl", "N≡N"],
                    answer: 3,
                    explain: "N≡N triple bond has highest bond dissociation energy (~945 kJ/mol) due to strong bonding interaction"
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY (45 questions)
                {
                    id: "b1",
                    text: "Bundle sheath cells in C₄ plants contain:",
                    options: ["PEP carboxylase", "RuBisCO", "Pyruvate phosphate dikinase", "Malic enzyme"],
                    answer: 1,
                    explain: "Bundle sheath cells in C₄ plants contain RuBisCO enzyme for Calvin cycle, while mesophyll cells contain PEP carboxylase"
                },
                {
                    id: "b2",
                    text: "Scutellum in monocot seed represents:",
                    options: ["Cotyledon", "Endosperm", "Plumule", "Radicle"],
                    answer: 0,
                    explain: "Scutellum is modified cotyledon in monocot seeds that helps in absorption of nutrients from endosperm"
                },
                {
                    id: "b3",
                    text: "Which hormone causes leaf abscission?",
                    options: ["Auxin", "Cytokinin", "Ethylene", "Gibberellin"],
                    answer: 2,
                    explain: "Ethylene promotes leaf abscission by weakening abscission zone through enzymatic breakdown of cell wall"
                },
                {
                    id: "b4",
                    text: "Interfascicular cambium develops from:",
                    options: ["Pericycle", "Medullary rays", "Cortex", "Epidermis"],
                    answer: 1,
                    explain: "Interfascicular cambium develops from parenchymatous cells of medullary rays between vascular bundles"
                },
                {
                    id: "b5",
                    text: "Prothallus represents:",
                    options: ["Sporophyte", "Gametophyte", "Sporangium", "Archegonium"],
                    answer: 1,
                    explain: "Prothallus is heart-shaped gametophyte generation in ferns that bears archegonia and antheridia"
                },
                {
                    id: "b6",
                    text: "Pollination by water is called:",
                    options: ["Anemophily", "Hydrophily", "Entomophily", "Ornithophily"],
                    answer: 1,
                    explain: "Hydrophily is pollination by water, seen in aquatic plants like Vallisneria and Zostera"
                },
                {
                    id: "b7",
                    text: "C₄ cycle was discovered by:",
                    options: ["Calvin", "Hatch and Slack", "Blackman", "Arnon"],
                    answer: 1,
                    explain: "C₄ pathway (Hatch-Slack pathway) was discovered by Hatch and Slack in 1966"
                },
                {
                    id: "b8",
                    text: "Casparian strips are made of:",
                    options: ["Lignin", "Suberin", "Cellulose", "Pectin"],
                    answer: 1,
                    explain: "Casparian strips in endodermis are made of suberin, which is impermeable to water and minerals"
                },
                {
                    id: "b9",
                    text: "Verticillaster inflorescence is found in:",
                    options: ["Tulsi", "Mustard", "Sunflower", "Marigold"],
                    answer: 0,
                    explain: "Verticillaster (false whorl) inflorescence is characteristic of Lamiaceae family plants like Tulsi (Ocimum)"
                },
                {
                    id: "b10",
                    text: "Parthenocarpy is induced by:",
                    options: ["Auxin", "Gibberellin", "Both auxin and gibberellin", "Cytokinin"],
                    answer: 2,
                    explain: "Both auxin and gibberellins can induce parthenocarpy (development of fruit without fertilization)"
                },
                {
                    id: "b11",
                    text: "Richmond-Lang effect is related to:",
                    options: ["Auxin", "Cytokinin", "Gibberellin", "Ethylene"],
                    answer: 1,
                    explain: "Richmond-Lang effect demonstrates cytokinin's role in preventing senescence and maintaining chlorophyll"
                },
                {
                    id: "b12",
                    text: "Non-cyclic photophosphorylation produces:",
                    options: ["Only ATP", "Only NADPH", "Both ATP and NADPH", "Neither ATP nor NADPH"],
                    answer: 2,
                    explain: "Non-cyclic photophosphorylation involves both PSI and PSII, producing both ATP and NADPH"
                },
                {
                    id: "b13",
                    text: "Critical photoperiod is:",
                    options: ["Maximum day length", "Minimum day length", "Duration that determines flowering", "Night length"],
                    answer: 2,
                    explain: "Critical photoperiod is specific day/night length that triggers flowering response in photoperiodic plants"
                },
                {
                    id: "b14",
                    text: "DPD (Diffusion Pressure Deficit) is equal to:",
                    options: ["OP - TP", "OP + TP", "TP - OP", "OP × TP"],
                    answer: 0,
                    explain: "DPD = OP - TP where OP is osmotic pressure and TP is turgor pressure. DPD represents water absorbing capacity"
                },
                {
                    id: "b15",
                    text: "Loading of phloem is:",
                    options: ["Passive process", "Active process", "Osmotic process", "Diffusion"],
                    answer: 1,
                    explain: "Phloem loading requires energy (ATP) for active transport of sucrose into sieve tubes against concentration gradient"
                },
                {
                    id: "b16",
                    text: "Aleurone layer in cereal grains is:",
                    options: ["Diploid", "Triploid", "Haploid", "Tetraploid"],
                    answer: 1,
                    explain: "Aleurone layer is part of endosperm, which is triploid (3n) tissue formed by triple fusion"
                },
                {
                    id: "b17",
                    text: "Phyllotaxy 2/5 means:",
                    options: ["2 leaves at 5 nodes", "5 leaves at 2 nodes", "2 complete turns in 5 leaves", "5 complete turns in 2 leaves"],
                    answer: 2,
                    explain: "In 2/5 phyllotaxy, genetic spiral makes 2 complete turns while passing through 5 leaves"
                },
                {
                    id: "b18",
                    text: "Transfusion tissue is found in:",
                    options: ["Angiosperms", "Gymnosperms", "Pteridophytes", "Bryophytes"],
                    answer: 1,
                    explain: "Transfusion tissue is specialized conducting tissue found in gymnosperm leaves for lateral transport"
                },
                {
                    id: "b19",
                    text: "Bacteriochlorophyll differs from chlorophyll in:",
                    options: ["Mg content", "Phytol tail", "Ring structure", "All of these"],
                    answer: 2,
                    explain: "Bacteriochlorophyll has different ring structure with additional double bonds and side chains compared to chlorophyll"
                },
                {
                    id: "b20",
                    text: "Chalazogamy occurs in:",
                    options: ["Normal plants", "Casuarina", "Cucurbita", "Both Casuarina and Cucurbita"],
                    answer: 3,
                    explain: "Chalazogamy (pollen tube entry through chalaza) occurs in plants like Casuarina and some Cucurbitaceae members"
                },
                {
                    id: "b21",
                    text: "Blue light receptors in plants are:",
                    options: ["Phytochromes", "Cryptochromes", "Phototropins", "Both cryptochromes and phototropins"],
                    answer: 3,
                    explain: "Both cryptochromes and phototropins are blue light receptors controlling various photomorphogenic responses"
                },
                {
                    id: "b22",
                    text: "Pollen-pistil interaction involves:",
                    options: ["Recognition", "Rejection", "Both recognition and rejection", "Neither"],
                    answer: 2,
                    explain: "Pollen-pistil interaction involves molecular recognition of compatible pollen and rejection of incompatible pollen"
                },
                {
                    id: "b23",
                    text: "Photorespiration occurs in:",
                    options: ["Peroxisomes only", "Chloroplasts only", "Mitochondria only", "All three organelles"],
                    answer: 3,
                    explain: "Photorespiration involves chloroplasts (oxygenation), peroxisomes (glycolate metabolism), and mitochondria (glycine conversion)"
                },
                {
                    id: "b24",
                    text: "Florigen is:",
                    options: ["Gibberellin", "Cytokinin", "Hypothetical flowering hormone", "Auxin"],
                    answer: 2,
                    explain: "Florigen is hypothetical flowering hormone proposed by Chailakhyan, now identified as FT protein"
                },
                {
                    id: "b25",
                    text: "Bordered pits are characteristic of:",
                    options: ["Parenchyma", "Collenchyma", "Sclerenchyma", "Tracheary elements"],
                    answer: 3,
                    explain: "Bordered pits are specialized pits found in tracheids and vessel elements for controlled water flow"
                },
                {
                    id: "b26",
                    text: "RuBP regeneration requires:",
                    options: ["3 ATP", "6 ATP", "9 ATP", "12 ATP"],
                    answer: 2,
                    explain: "Regeneration of 6 RuBP molecules in Calvin cycle requires 9 ATP molecules"
                },
                {
                    id: "b27",
                    text: "Secondary growth in monocots is:",
                    options: ["Always absent", "Always present", "Present in some", "Present in all"],
                    answer: 2,
                    explain: "Some monocots like Dracaena, Yucca show secondary growth through anomalous cambial activity"
                },
                {
                    id: "b28",
                    text: "Megaspore mother cell is:",
                    options: ["Haploid", "Diploid", "Triploid", "Tetraploid"],
                    answer: 1,
                    explain: "Megaspore mother cell (megasporocyte) is diploid cell that undergoes meiosis to form four haploid megaspores"
                },
                {
                    id: "b29",
                    text: "Leg-haemoglobin in root nodules:",
                    options: ["Transports oxygen", "Scavenges oxygen", "Stores oxygen", "Produces oxygen"],
                    answer: 1,
                    explain: "Leg-haemoglobin maintains low oxygen concentration around nitrogenase by scavenging free oxygen"
                },
                {
                    id: "b30",
                    text: "Periderm includes:",
                    options: ["Cork only", "Cork cambium only", "Cork and cork cambium", "Cork, cork cambium and secondary cortex"],
                    answer: 3,
                    explain: "Periderm consists of cork (phellem), cork cambium (phellogen), and secondary cortex (phelloderm)"
                },
                {
                    id: "b31",
                    text: "Turgor movements are shown by:",
                    options: ["Mimosa", "Drosera", "Both Mimosa and Drosera", "Sunflower"],
                    answer: 2,
                    explain: "Both Mimosa (seismonasty) and Drosera (thigmonasty) show turgor movements in response to stimuli"
                },
                {
                    id: "b32",
                    text: "Spiral thickening is found in:",
                    options: ["Protoxylem", "Metaxylem", "Secondary xylem", "Phloem"],
                    answer: 0,
                    explain: "Protoxylem vessels have spiral (helical) thickening that allows stretching during organ growth"
                },
                {
                    id: "b33",
                    text: "Climacteric fruits are characterized by:",
                    options: ["High respiration during ripening", "Low respiration during ripening", "No ethylene production", "Slow ripening"],
                    answer: 0,
                    explain: "Climacteric fruits show sharp increase in respiration rate and ethylene production during ripening"
                },
                {
                    id: "b34",
                    text: "Xylem loading is:",
                    options: ["Active", "Passive", "Semi-active", "Osmotic"],
                    answer: 1,
                    explain: "Water enters xylem passively through apoplastic and symplastic pathways driven by transpiration pull"
                },
                {
                    id: "b35",
                    text: "Malate valve operates in:",
                    options: ["C₃ plants", "C₄ plants", "CAM plants", "All plants"],
                    answer: 2,
                    explain: "Malate valve operates in CAM plants, allowing CO₂ storage as malate at night and release during day"
                },
                {
                    id: "b36",
                    text: "Root pressure can be demonstrated by:",
                    options: ["Guttation", "Bleeding", "Both guttation and bleeding", "Transpiration"],
                    answer: 2,
                    explain: "Both guttation (water droplets from hydathodes) and bleeding (xylem sap from cut stems) demonstrate root pressure"
                },
                {
                    id: "b37",
                    text: "Distichous phyllotaxy is:",
                    options: ["1/2", "1/3", "2/5", "3/8"],
                    answer: 0,
                    explain: "Distichous phyllotaxy has 1/2 arrangement where leaves alternate in two rows"
                },
                {
                    id: "b38",
                    text: "Bakane disease of rice is caused by:",
                    options: ["Bacteria", "Virus", "Fungus", "Nematode"],
                    answer: 2,
                    explain: "Bakane (foolish seedling) disease is caused by fungus Gibberella fujikuroi, which produces gibberellins"
                },
                {
                    id: "b39",
                    text: "CAM plants keep stomata:",
                    options: ["Open during day", "Open during night", "Always closed", "Always open"],
                    answer: 1,
                    explain: "CAM plants open stomata at night to minimize water loss while fixing CO₂ as organic acids"
                },
                {
                    id: "b40",
                    text: "Pressure flow hypothesis was proposed by:",
                    options: ["Münch", "Dixon", "Curtis", "Crafts"],
                    answer: 0,
                    explain: "Mass flow or pressure flow hypothesis for phloem transport was proposed by Ernst Münch in 1930"
                },
                {
                    id: "b41",
                    text: "Quiescent center is found in:",
                    options: ["Shoot apex", "Root apex", "Cambium", "Cork cambium"],
                    answer: 1,
                    explain: "Quiescent center is region of slowly dividing cells in root apical meristem that acts as reserve"
                },
                {
                    id: "b42",
                    text: "Z-scheme represents:",
                    options: ["Calvin cycle", "Light reaction", "Photorespiration", "Respiration"],
                    answer: 1,
                    explain: "Z-scheme depicts electron transport during light reactions showing energy levels of PSI and PSII"
                },
                {
                    id: "b43",
                    text: "Viviparous germination occurs in:",
                    options: ["Desert plants", "Mangroves", "Alpine plants", "Aquatic plants"],
                    answer: 1,
                    explain: "Viviparous germination (seed germination on parent plant) is adaptation found in mangrove plants"
                },
                {
                    id: "b44",
                    text: "Companion cells are connected to sieve tubes by:",
                    options: ["Plasmodesmata", "Pits", "Perforations", "Channels"],
                    answer: 0,
                    explain: "Companion cells are connected to sieve tube elements through numerous plasmodesmata for metabolic support"
                },
                {
                    id: "b45",
                    text: "Kinetin belongs to:",
                    options: ["Auxins", "Gibberellins", "Cytokinins", "Abscisic acid"],
                    answer: 2,
                    explain: "Kinetin (6-furfurylaminopurine) is first discovered cytokinin, promoting cell division and preventing senescence"
                },
                
                // ZOOLOGY (45 questions)
                {
                    id: "z1",
                    text: "Which bone forms the posterior part of hard palate?",
                    options: ["Maxilla", "Palatine", "Sphenoid", "Vomer"],
                    answer: 1,
                    explain: "Palatine bones form the posterior part of hard palate, while maxilla forms the anterior part"
                },
                {
                    id: "z2",
                    text: "Glucagon is secreted by:",
                    options: ["Alpha cells", "Beta cells", "Delta cells", "PP cells"],
                    answer: 0,
                    explain: "Alpha cells in pancreatic islets secrete glucagon hormone that raises blood glucose levels"
                },
                {
                    id: "z3",
                    text: "Prothrombin is converted to thrombin by:",
                    options: ["Thromboplastin", "Fibrinogen", "Calcium", "Thrombokinase"],
                    answer: 3,
                    explain: "Thrombokinase (thromboplastin + calcium) converts prothrombin to thrombin in blood clotting cascade"
                },
                {
                    id: "z4",
                    text: "Henle's loop is part of:",
                    options: ["Glomerulus", "Bowman's capsule", "Nephron tubule", "Collecting duct"],
                    answer: 2,
                    explain: "Loop of Henle is U-shaped part of nephron tubule consisting of descending and ascending limbs"
                },
                {
                    id: "z5",
                    text: "Arbor vitae is found in:",
                    options: ["Cerebrum", "Cerebellum", "Medulla", "Spinal cord"],
                    answer: 1,
                    explain: "Arbor vitae (tree of life) is white matter in cerebellum showing characteristic branched pattern"
                },
                {
                    id: "z6",
                    text: "Hemoglobin is broken down in:",
                    options: ["Liver", "Spleen", "Kidney", "Bone marrow"],
                    answer: 1,
                    explain: "Old RBCs are hemolyzed in spleen where hemoglobin is broken down into heme and globin"
                },
                {
                    id: "z7",
                    text: "Trypsinogen is activated by:",
                    options: ["Pepsin", "Enterokinase", "HCl", "Chymotrypsin"],
                    answer: 1,
                    explain: "Enterokinase (enteropeptidase) secreted by duodenal mucosa activates trypsinogen to trypsin"
                },
                {
                    id: "z8",
                    text: "Purkinje fibers are located in:",
                    options: ["Atria", "Ventricles", "AV node", "SA node"],
                    answer: 1,
                    explain: "Purkinje fibers form conduction system in ventricular walls for rapid impulse transmission"
                },
                {
                    id: "z9",
                    text: "Addison's disease is caused by hypoactivity of:",
                    options: ["Thyroid", "Parathyroid", "Adrenal cortex", "Adrenal medulla"],
                    answer: 2,
                    explain: "Addison's disease results from insufficient cortisol and aldosterone production by adrenal cortex"
                },
                {
                    id: "z10",
                    text: "Corpus luteum is formed from:",
                    options: ["Primary follicle", "Secondary follicle", "Graafian follicle", "Atretic follicle"],
                    answer: 2,
                    explain: "Corpus luteum develops from ruptured Graafian follicle after ovulation under LH influence"
                },
                {
                    id: "z11",
                    text: "Night blindness is caused by deficiency of:",
                    options: ["Vitamin A", "Vitamin B₁", "Vitamin C", "Vitamin D"],
                    answer: 0,
                    explain: "Vitamin A deficiency affects rhodopsin synthesis in rod cells, causing night blindness (nyctalopia)"
                },
                {
                    id: "z12",
                    text: "Memory cells are derived from:",
                    options: ["T cells", "B cells", "Both T and B cells", "NK cells"],
                    answer: 2,
                    explain: "Both T and B lymphocytes can differentiate into memory cells for long-term immunity"
                },
                {
                    id: "z13",
                    text: "Countercurrent mechanism operates in:",
                    options: ["Bowman's capsule", "Proximal tubule", "Loop of Henle", "Distal tubule"],
                    answer: 2,
                    explain: "Countercurrent multiplier mechanism in loop of Henle concentrates urine by creating osmotic gradient"
                },
                {
                    id: "z14",
                    text: "Milk ejection reflex involves:",
                    options: ["Prolactin", "Oxytocin", "Estrogen", "Progesterone"],
                    answer: 1,
                    explain: "Oxytocin causes contraction of myoepithelial cells around mammary alveoli for milk ejection"
                },
                {
                    id: "z15",
                    text: "Rh factor was discovered in:",
                    options: ["Humans", "Rhesus monkey", "Chimpanzee", "Rats"],
                    answer: 1,
                    explain: "Rh factor was first discovered in red blood cells of Rhesus macaque monkeys, hence the name"
                },
                {
                    id: "z16",
                    text: "Adam's apple is formed by:",
                    options: ["Thyroid cartilage", "Cricoid cartilage", "Arytenoid cartilage", "Epiglottis"],
                    answer: 0,
                    explain: "Adam's apple is prominence of thyroid cartilage, more prominent in males due to testosterone"
                },
                {
                    id: "z17",
                    text: "Organ of Corti is located in:",
                    options: ["External ear", "Middle ear", "Cochlea", "Semicircular canals"],
                    answer: 2,
                    explain: "Organ of Corti is auditory receptor located on basilar membrane inside cochlea"
                },
                {
                    id: "z18",
                    text: "Creatinine clearance test measures:",
                    options: ["Liver function", "Kidney function", "Heart function", "Lung function"],
                    answer: 1,
                    explain: "Creatinine clearance is standard test for measuring glomerular filtration rate and kidney function"
                },
                {
                    id: "z19",
                    text: "Acromegaly is caused by excess of:",
                    options: ["Thyroxine", "Insulin", "Growth hormone", "Cortisol"],
                    answer: 2,
                    explain: "Acromegaly results from excess growth hormone secretion in adults, causing enlarged extremities"
                },
                {
                    id: "z20",
                    text: "Mittelschmerz is associated with:",
                    options: ["Menstruation", "Ovulation", "Fertilization", "Implantation"],
                    answer: 1,
                    explain: "Mittelschmerz is ovulation pain experienced by some women mid-cycle during follicle rupture"
                },
                {
                    id: "z21",
                    text: "Tricuspid valve is located between:",
                    options: ["Right atrium and right ventricle", "Left atrium and left ventricle", "Right ventricle and pulmonary artery", "Left ventricle and aorta"],
                    answer: 0,
                    explain: "Tricuspid valve prevents backflow from right ventricle to right atrium during ventricular systole"
                },
                {
                    id: "z22",
                    text: "Carbonic anhydrase is present in:",
                    options: ["Plasma", "RBCs", "WBCs", "Platelets"],
                    answer: 1,
                    explain: "Carbonic anhydrase in RBCs catalyzes reversible reaction between CO₂ and H₂O for CO₂ transport"
                },
                {
                    id: "z23",
                    text: "Osteoporosis is related to deficiency of:",
                    options: ["Calcitonin", "Parathormone", "Estrogen", "Testosterone"],
                    answer: 2,
                    explain: "Estrogen deficiency (especially post-menopause) leads to increased bone resorption causing osteoporosis"
                },
                {
                    id: "z24",
                    text: "Chemoreceptors for respiration are located in:",
                    options: ["Medulla oblongata", "Carotid body", "Aortic body", "All of these"],
                    answer: 3,
                    explain: "Central chemoreceptors in medulla and peripheral chemoreceptors in carotid and aortic bodies monitor CO₂, O₂, and pH"
                },
                {
                    id: "z25",
                    text: "Intercalated discs are found in:",
                    options: ["Skeletal muscle", "Cardiac muscle", "Smooth muscle", "All muscles"],
                    answer: 1,
                    explain: "Intercalated discs are specialized cell junctions between cardiac muscle cells for synchronized contraction"
                },
                {
                    id: "z26",
                    text: "Kupffer cells are found in:",
                    options: ["Liver", "Spleen", "Kidney", "Lung"],
                    answer: 0,
                    explain: "Kupffer cells are specialized macrophages in liver sinusoids that remove old RBCs and foreign particles"
                },
                {
                    id: "z27",
                    text: "Rickets is caused by deficiency of:",
                    options: ["Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D"],
                    answer: 3,
                    explain: "Vitamin D deficiency causes rickets in children due to impaired calcium absorption and bone mineralization"
                },
                {
                    id: "z28",
                    text: "Troponin is associated with:",
                    options: ["Thick filaments", "Thin filaments", "Z-line", "M-line"],
                    answer: 1,
                    explain: "Troponin complex is located on thin (actin) filaments and regulates muscle contraction with tropomyosin"
                },
                {
                    id: "z29",
                    text: "Ketone bodies are produced in:",
                    options: ["Muscle", "Liver", "Kidney", "Brain"],
                    answer: 1,
                    explain: "Ketone bodies are produced in liver during prolonged fasting or diabetes as alternative fuel source"
                },
                {
                    id: "z30",
                    text: "Peyer's patches are found in:",
                    options: ["Stomach", "Duodenum", "Jejunum", "Ileum"],
                    answer: 3,
                    explain: "Peyer's patches are lymphoid aggregations in ileum that form part of gut-associated lymphoid tissue"
                },
                {
                    id: "z31",
                    text: "Shivering is controlled by:",
                    options: ["Cerebrum", "Cerebellum", "Hypothalamus", "Medulla"],
                    answer: 2,
                    explain: "Hypothalamus detects temperature changes and initiates shivering thermogenesis to generate heat"
                },
                {
                    id: "z32",
                    text: "Implantation occurs on day:",
                    options: ["7th", "14th", "21st", "28th"],
                    answer: 0,
                    explain: "Blastocyst implants in endometrium around 7th day after fertilization in human pregnancy"
                },
                {
                    id: "z33",
                    text: "Thrombopoietin stimulates production of:",
                    options: ["RBCs", "WBCs", "Platelets", "Plasma proteins"],
                    answer: 2,
                    explain: "Thrombopoietin produced by liver and kidneys stimulates megakaryocyte differentiation and platelet production"
                },
                {
                    id: "z34",
                    text: "Circadian rhythm is controlled by:",
                    options: ["Pituitary", "Pineal", "Thyroid", "Adrenal"],
                    answer: 1,
                    explain: "Pineal gland secretes melatonin in response to light-dark cycles, controlling circadian rhythms"
                },
                {
                    id: "z35",
                    text: "Aqueous humor is secreted by:",
                    options: ["Cornea", "Iris", "Ciliary body", "Choroid"],
                    answer: 2,
                    explain: "Ciliary processes of ciliary body secrete aqueous humor that maintains intraocular pressure"
                },
                {
                    id: "z36",
                    text: "Normal urine output per day is:",
                    options: ["500 ml", "1-1.5 L", "2-3 L", "3-4 L"],
                    answer: 1,
                    explain: "Normal daily urine output is 1-1.5 liters, varying with fluid intake and physiological conditions"
                },
                {
                    id: "z37",
                    text: "Taste buds are innervated by:",
                    options: ["Facial nerve only", "Glossopharyngeal nerve only", "Both facial and glossopharyngeal", "Trigeminal nerve"],
                    answer: 2,
                    explain: "Taste buds are innervated by facial nerve (anterior 2/3 tongue) and glossopharyngeal nerve (posterior 1/3)"
                },
                {
                    id: "z38",
                    text: "Insulin resistance is characteristic of:",
                    options: ["Type 1 diabetes", "Type 2 diabetes", "Both types", "Neither type"],
                    answer: 1,
                    explain: "Type 2 diabetes is characterized by insulin resistance where cells don't respond properly to insulin"
                },
                {
                    id: "z39",
                    text: "HCG is secreted by:",
                    options: ["Corpus luteum", "Placenta", "Endometrium", "Ovary"],
                    answer: 1,
                    explain: "Human chorionic gonadotropin (hCG) is secreted by trophoblast cells of developing placenta"
                },
                {
                    id: "z40",
                    text: "Floating ribs are:",
                    options: ["11th and 12th", "9th and 10th", "8th and 9th", "10th and 11th"],
                    answer: 0,
                    explain: "11th and 12th pairs are floating ribs that don't attach to sternum or other ribs anteriorly"
                },
                {
                    id: "z41",
                    text: "Intrinsic factor is secreted by:",
                    options: ["Chief cells", "Parietal cells", "Mucus cells", "G cells"],
                    answer: 1,
                    explain: "Parietal cells secrete intrinsic factor necessary for vitamin B₁₂ absorption in terminal ileum"
                },
                {
                    id: "z42",
                    text: "Synovial membrane lines:",
                    options: ["Cartilage", "Joint capsule", "Bone", "Ligaments"],
                    answer: 1,
                    explain: "Synovial membrane lines inner surface of joint capsule and secretes synovial fluid"
                },
                {
                    id: "z43",
                    text: "Reverse T₃ is:",
                    options: ["Active hormone", "Inactive metabolite", "Storage form", "Precursor"],
                    answer: 1,
                    explain: "Reverse T₃ (rT₃) is inactive metabolite of thyroxine formed during stress or illness"
                },
                {
                    id: "z44",
                    text: "Chordae tendineae prevent:",
                    options: ["Valve prolapse", "Valve stenosis", "Arrhythmia", "Heart block"],
                    answer: 0,
                    explain: "Chordae tendineae (heartstrings) connect AV valve cusps to papillary muscles, preventing valve prolapse"
                },
                {
                    id: "z45",
                    text: "Cortical reaction prevents:",
                    options: ["Fertilization", "Polyspermy", "Implantation", "Cleavage"],
                    answer: 1,
                    explain: "Cortical reaction after sperm entry causes zona pellucida hardening to prevent polyspermy"
                }
            ]
        }
    ]
};
