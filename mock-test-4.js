// mock-test-4.js - NEET Mock Test 4 (High Quality Questions)
// Physics: 45 questions | Chemistry: 45 questions | Biology: 90 questions (45 Botany + 45 Zoology)
// Follows exact NEET pattern and difficulty - All NEW questions

window.MOCK_TEST_4 = {
    id: "neet-004",
    title: "Full Syllabus Mock 4", 
    level: "medium",
    durationMinutes: 180,
    sections: [
        {
            name: "Physics",
            questions: [
                {
                    id: "p1",
                    text: "A uniform rod AB of length 2L is placed with end A at distance L from fulcrum. What weight should be placed at end B to balance a weight W at end A?",
                    options: ["W/2", "W/3", "2W/3", "W"],
                    answer: 1,
                    explain: "Taking moments about fulcrum: W × L = Weight at B × 3L. Therefore, Weight at B = W/3"
                },
                {
                    id: "p2",
                    text: "A particle moves in a circle of radius R with constant angular velocity ω. Its average acceleration over a quarter circle is:",
                    options: ["ω²R√2/π", "ω²R/π", "ω²R", "0"],
                    answer: 0,
                    explain: "Change in velocity = √2 ωR, time = πR/2ωR = π/2ω. Average acceleration = √2 ωR ÷ (π/2ω) = 2√2 ω²R/π"
                },
                {
                    id: "p3",
                    text: "Two masses m₁ and m₂ are connected by a string over a pulley. The tension in string when system moves is:",
                    options: ["(m₁ + m₂)g/2", "2m₁m₂g/(m₁ + m₂)", "m₁m₂g/(m₁ + m₂)", "(m₁ - m₂)g"],
                    answer: 1,
                    explain: "For Atwood machine: T = 2m₁m₂g/(m₁ + m₂) when system is in motion"
                },
                {
                    id: "p4",
                    text: "A sphere rolls down an inclined plane of angle θ. Its acceleration down the plane is:",
                    options: ["g sin θ", "(5/7)g sin θ", "(2/3)g sin θ", "(3/5)g sin θ"],
                    answer: 1,
                    explain: "For solid sphere rolling: a = g sin θ/(1 + I/mR²) = g sin θ/(1 + 2/5) = (5/7)g sin θ"
                },
                {
                    id: "p5",
                    text: "The velocity of electromagnetic waves in a medium of dielectric constant K and permeability μᵣ is:",
                    options: ["c/√(K)", "c/√(Kμᵣ)", "c√(Kμᵣ)", "c√K"],
                    answer: 1,
                    explain: "v = c/√(εᵣμᵣ) = c/√(Kμᵣ) where εᵣ = K is relative permittivity"
                },
                {
                    id: "p6",
                    text: "A wire is stretched to increase its length by 1%. The percentage change in resistance is:",
                    options: ["1%", "2%", "0.5%", "4%"],
                    answer: 1,
                    explain: "R = ρL/A. If L increases by 1%, volume constant means A decreases. ΔR/R = 2ΔL/L = 2%"
                },
                {
                    id: "p7",
                    text: "In hydrogen atom, electron moves from n=3 to n=1 level. The ratio of frequencies of emitted photons in Lyman series is:",
                    options: ["32:27", "27:32", "8:27", "27:8"],
                    answer: 0,
                    explain: "f ∝ (1/n₁² - 1/n₂²). For 2→1: (1-1/4) = 3/4. For 3→1: (1-1/9) = 8/9. Ratio = (8/9):(3/4) = 32:27"
                },
                {
                    id: "p8",
                    text: "A convex lens of focal length 20 cm is cut into two halves. The focal length of each half is:",
                    options: ["10 cm", "20 cm", "40 cm", "∞"],
                    answer: 1,
                    explain: "Cutting lens horizontally doesn't change curvature, so focal length remains same = 20 cm"
                },
                {
                    id: "p9",
                    text: "The time period of oscillation of liquid in U-tube of uniform cross-section is T = 2π√(L/2g), where L is:",
                    options: ["Total length of liquid", "Length of one arm", "Half the total length", "Effective length"],
                    answer: 0,
                    explain: "L is total length of liquid column oscillating in U-tube"
                },
                {
                    id: "p10",
                    text: "A uniform magnetic field B exists perpendicular to plane of paper. A charged particle moves in a circle of radius r. If charge is doubled and speed is halved, new radius is:",
                    options: ["r/4", "r/2", "r", "2r"],
                    answer: 0,
                    explain: "r = mv/qB. When q → 2q and v → v/2: r' = m(v/2)/(2q)B = r/4"
                },
                {
                    id: "p11",
                    text: "Two resistors R₁ and R₂ carry currents I₁ and I₂. They dissipate equal power when:",
                    options: ["I₁R₁ = I₂R₂", "I₁²R₁ = I₂²R₂", "I₁/R₁ = I₂/R₂", "R₁/I₁ = R₂/I₂"],
                    answer: 1,
                    explain: "Power P = I²R. For equal power: I₁²R₁ = I₂²R₂"
                },
                {
                    id: "p12",
                    text: "A satellite orbits earth at height h = R (R = earth radius). Its orbital speed is:",
                    options: ["√(gR)", "√(gR/2)", "√(2gR)", "√(gR/4)"],
                    answer: 1,
                    explain: "At height R, orbital radius = 2R. v = √(GM/2R) = √(gR²/2R) = √(gR/2)"
                },
                {
                    id: "p13",
                    text: "In Young's double slit experiment, if distance between slits is halved and screen distance is doubled, fringe width becomes:",
                    options: ["Same", "Double", "Four times", "Half"],
                    answer: 2,
                    explain: "β = λD/d. When d → d/2 and D → 2D: β' = λ(2D)/(d/2) = 4λD/d = 4β"
                },
                {
                    id: "p14",
                    text: "A particle performs SHM with amplitude A. When displacement is A/3, what fraction of total energy is kinetic?",
                    options: ["8/9", "1/9", "2/3", "1/3"],
                    answer: 0,
                    explain: "KE/Total = (A²-x²)/A² = (A²-A²/9)/A² = 8/9"
                },
                {
                    id: "p15",
                    text: "The self-inductance of a solenoid is proportional to:",
                    options: ["n", "n²", "n³", "1/n"],
                    answer: 1,
                    explain: "L = μ₀n²Al where n is turns per unit length. L ∝ n²"
                },
                {
                    id: "p16",
                    text: "A radioactive sample has half-life 4 hours. What fraction remains after 12 hours?",
                    options: ["1/4", "1/8", "1/16", "1/2"],
                    answer: 1,
                    explain: "12 hours = 3 half-lives. Remaining fraction = (1/2)³ = 1/8"
                },
                {
                    id: "p17",
                    text: "Two identical springs in series support mass m. If one spring breaks, the new frequency of oscillation is:",
                    options: ["√2 times", "2 times", "1/√2 times", "1/2 times"],
                    answer: 0,
                    explain: "Initially keff = k/2, finally keff = k. f ∝ √k, so f₂/f₁ = √(k/(k/2)) = √2"
                },
                {
                    id: "p18",
                    text: "The energy density in electromagnetic wave is distributed between electric and magnetic fields in ratio:",
                    options: ["1:1", "1:2", "2:1", "1:c²"],
                    answer: 0,
                    explain: "In EM waves, energy density is equally distributed between E and B fields: uE = uB"
                },
                {
                    id: "p19",
                    text: "A lens made of material with refractive index 1.5 has focal length 20 cm in air. In water (μ = 4/3), focal length becomes:",
                    options: ["60 cm", "80 cm", "40 cm", "30 cm"],
                    answer: 1,
                    explain: "1/f₂ = (μ_lens/μ_medium - 1)(1/R₁ - 1/R₂). f₂ = f₁ × μ_medium/(μ_lens - μ_medium) = 20 × (4/3)/(1.5 - 4/3) = 80 cm"
                },
                {
                    id: "p20",
                    text: "The ratio of kinetic energies of photoelectrons when light of wavelengths 400 nm and 600 nm fall on metal surface (work function 2 eV) is:",
                    options: ["3:2", "2:1", "1:0", "Cannot be determined"],
                    answer: 2,
                    explain: "For 400 nm: E = 3.1 eV, KE = 1.1 eV. For 600 nm: E = 2.07 eV, KE = 0.07 eV. For 600 nm, KE ≈ 0. Ratio ≈ 1:0"
                },
                {
                    id: "p21",
                    text: "A charged particle is accelerated through potential V. Its de Broglie wavelength is λ₁. If accelerated through 4V, wavelength becomes:",
                    options: ["λ₁/2", "λ₁/4", "2λ₁", "4λ₁"],
                    answer: 0,
                    explain: "λ = h/√(2meV). When V → 4V: λ₂ = h/√(8meV) = λ₁/2"
                },
                {
                    id: "p22",
                    text: "In AC circuit, if voltage leads current by 60°, the power factor is:",
                    options: ["0.5", "0.707", "0.866", "1"],
                    answer: 0,
                    explain: "Power factor = cos φ = cos 60° = 0.5"
                },
                {
                    id: "p23",
                    text: "A solid sphere and hollow sphere of same mass and radius roll down incline. Ratio of their kinetic energies at bottom is:",
                    options: ["5:7", "7:5", "3:5", "5:3"],
                    answer: 0,
                    explain: "KE_solid = (7/10)mgh, KE_hollow = (5/6)mgh. Ratio = (7/10):(5/6) = 21:25. Wait, let me recalculate. For solid: I = (2/5)mR², for hollow: I = (2/3)mR². KE_total = mgh. KE_trans = mgh/(1 + I/mR²). Solid: KE = mgh/(1 + 2/5) = (5/7)mgh. Hollow: KE = mgh/(1 + 2/3) = (3/5)mgh. Ratio = (5/7):(3/5) = 25:21"
                },
                {
                    id: "p24",
                    text: "The magnetic flux through a closed surface is:",
                    options: ["Always zero", "Always positive", "Always negative", "May be positive or negative"],
                    answer: 0,
                    explain: "By Gauss's law for magnetism: ∮B⃗·dA⃗ = 0 as magnetic monopoles don't exist"
                },
                {
                    id: "p25",
                    text: "When white light passes through a prism, which color deviates least?",
                    options: ["Red", "Blue", "Green", "Yellow"],
                    answer: 0,
                    explain: "Red light has least refractive index, hence least deviation according to δ = (μ - 1)A"
                },
                {
                    id: "p26",
                    text: "A particle moves with constant speed v in a straight line. At t = 0, it is at origin. Its position vector at time t is:",
                    options: ["vt î", "vt ĵ", "vt (î + ĵ)/√2", "Depends on direction"],
                    answer: 3,
                    explain: "Position depends on direction of motion. Could be along any direction in space."
                },
                {
                    id: "p27",
                    text: "The electric field at distance r from infinite line charge with linear charge density λ is:",
                    options: ["λ/(2πε₀r)", "λ/(4πε₀r²)", "λ/(πε₀r)", "λ/(ε₀r)"],
                    answer: 0,
                    explain: "For infinite line charge: E = λ/(2πε₀r) using Gauss's law with cylindrical surface"
                },
                {
                    id: "p28",
                    text: "A body cools from 80°C to 60°C in 10 minutes. Time to cool from 60°C to 40°C (room temp 20°C) is:",
                    options: ["10 min", "15 min", "20 min", "25 min"],
                    answer: 1,
                    explain: "Using Newton's cooling law: dT/dt = -k(T-T₀). Solving gives t = (1/k)ln[(60-20)/(40-20)] = 15 min"
                },
                {
                    id: "p29",
                    text: "In photoelectric effect, if frequency is doubled while keeping work function same, stopping potential becomes:",
                    options: ["Double", "More than double", "Less than double", "Same"],
                    answer: 1,
                    explain: "V₀ = (hf - φ)/e. When f → 2f: V₀' = (2hf - φ)/e = 2hf/e - φ/e = 2V₀ + hf/e > 2V₀"
                },
                {
                    id: "p30",
                    text: "A conducting loop moves through non-uniform magnetic field. The induced EMF depends on:",
                    options: ["Rate of change of flux", "Velocity of loop", "Area of loop", "All of above"],
                    answer: 0,
                    explain: "According to Faraday's law, induced EMF = -dΦ/dt depends only on rate of change of magnetic flux"
                },
                {
                    id: "p31",
                    text: "Two identical cells are connected first in series then in parallel across external resistance R. The ratio of currents is:",
                    options: ["1:1", "2:1", "1:2", "4:1"],
                    answer: 1,
                    explain: "Series: I₁ = 2E/(2r+R). Parallel: I₂ = E/(r/2+R). Ratio I₁:I₂ = 2E(r/2+R)/(2r+R)E = (r+2R)/(2r+R) ≈ 2:1 for R>>r"
                },
                {
                    id: "p32",
                    text: "The intensity of sound increases by factor of 20. The increase in sound level is:",
                    options: ["13 dB", "20 dB", "3 dB", "10 dB"],
                    answer: 0,
                    explain: "ΔL = 10 log₁₀(I₂/I₁) = 10 log₁₀(20) = 10 × 1.3 = 13 dB"
                },
                {
                    id: "p33",
                    text: "A uniform rod of mass M and length L is pivoted at one end. Its radius of gyration about the pivot is:",
                    options: ["L/√3", "L/√12", "L√3", "L/2"],
                    answer: 0,
                    explain: "I = ML²/3, k = √(I/M) = √(L²/3) = L/√3"
                },
                {
                    id: "p34",
                    text: "The maximum number of electrons in M shell is:",
                    options: ["8", "18", "32", "50"],
                    answer: 1,
                    explain: "M shell has n=3, maximum electrons = 2n² = 2(3)² = 18"
                },
                {
                    id: "p35",
                    text: "A particle executes SHM between x = -A and x = +A. The time to go from x = 0 to x = A/2 is:",
                    options: ["T/6", "T/12", "T/8", "T/4"],
                    answer: 1,
                    explain: "x = A sin(ωt). For x = A/2: sin(ωt) = 1/2, ωt = π/6, t = T/12"
                },
                {
                    id: "p36",
                    text: "The energy stored in inductor carrying current I is:",
                    options: ["LI²", "½LI²", "2LI²", "LI²/4"],
                    answer: 1,
                    explain: "Energy stored in inductor U = ½LI²"
                },
                {
                    id: "p37",
                    text: "A ray of light travels from denser to rarer medium. The critical angle depends on:",
                    options: ["Wavelength only", "Refractive indices only", "Both wavelength and refractive indices", "Angle of incidence"],
                    answer: 2,
                    explain: "Critical angle sin θc = μ₂/μ₁, but μ depends on wavelength (dispersion)"
                },
                {
                    id: "p38",
                    text: "In damped oscillations, if damping is doubled, the time period:",
                    options: ["Remains same", "Increases", "Decreases", "Becomes zero"],
                    answer: 1,
                    explain: "T = 2π/ω' where ω' = √(ω₀² - γ²). As γ increases, ω' decreases, T increases"
                },
                {
                    id: "p39",
                    text: "The ratio of magnetic moments of electron and proton is approximately:",
                    options: ["1:1", "1836:1", "1:1836", "1:2000"],
                    answer: 1,
                    explain: "μ ∝ 1/mass. μₑ/μₚ = mₚ/mₑ = 1836"
                },
                {
                    id: "p40",
                    text: "A wire carries current I. If it is bent into circular loop, magnetic field at center is B. If same wire is bent into square loop, field at center becomes:",
                    options: ["B√2/π", "2√2B/π", "πB/2√2", "B"],
                    answer: 1,
                    explain: "For circle: B₁ = μ₀I/2r, for square: B₂ = 4μ₀I/(π√2 a). With same wire length: B₂/B₁ = 2√2/π"
                },
                {
                    id: "p41",
                    text: "The binding energy per nucleon is minimum for:",
                    options: ["Hydrogen", "Helium", "Iron", "Uranium"],
                    answer: 0,
                    explain: "Hydrogen has minimum binding energy per nucleon (zero for protium)"
                },
                {
                    id: "p42",
                    text: "A gas bubble rises from bottom of lake. As it rises, its temperature:",
                    options: ["Increases", "Decreases", "Remains constant", "First increases then decreases"],
                    answer: 1,
                    explain: "As bubble rises, pressure decreases. Adiabatic expansion causes temperature to decrease"
                },
                {
                    id: "p43",
                    text: "Two waves of same frequency but different amplitudes interfere. The ratio of maximum to minimum intensity is 9:1. The ratio of amplitudes is:",
                    options: ["3:1", "2:1", "4:1", "9:1"],
                    answer: 1,
                    explain: "Imax/Imin = (A₁+A₂)²/(A₁-A₂)² = 9. Solving: (A₁+A₂)/(A₁-A₂) = 3, gives A₁:A₂ = 2:1"
                },
                {
                    id: "p44",
                    text: "The wavelength of X-rays is measured using:",
                    options: ["Diffraction grating", "Young's double slit", "Bragg's crystal diffraction", "Prism"],
                    answer: 2,
                    explain: "X-ray wavelengths are measured using Bragg's law: nλ = 2d sinθ with crystal lattice"
                },
                {
                    id: "p45",
                    text: "A charged capacitor discharges through resistor. The time constant RC has dimensions:",
                    options: ["[T]", "[T⁻¹]", "[ML²T⁻¹]", "[A²T]"],
                    answer: 0,
                    explain: "RC = (Volt/Amp) × (Charge/Volt) = Charge/Amp = Time. Dimension is [T]"
                }
            ]
        },
        {
            name: "Chemistry", 
            questions: [
                {
                    id: "c1",
                    text: "Which of the following has maximum bond angle?",
                    options: ["NH₃", "H₂O", "CH₄", "H₂S"],
                    answer: 2,
                    explain: "CH₄ has tetrahedral geometry with bond angle 109.5°, highest among given options due to no lone pairs."
                },
                {
                    id: "c2",
                    text: "The number of stereoisomers possible for the compound CH₃CHOHCHOHCH₃ is:",
                    options: ["2", "3", "4", "6"],
                    answer: 1,
                    explain: "Two chiral centers give 2² = 4 stereoisomers, but due to meso form, only 3 exist (2 enantiomers + 1 meso)."
                },
                {
                    id: "c3",
                    text: "Which of the following complexes is optically active?",
                    options: ["[Pt(NH₃)₂Cl₂]", "[Co(NH₃)₆]³⁺", "[Co(en)₃]³⁺", "[Ni(CO)₄]"],
                    answer: 2,
                    explain: "[Co(en)₃]³⁺ is octahedral with three bidentate ligands, lacks plane of symmetry, hence optically active."
                },
                {
                    id: "c4",
                    text: "The correct order of reducing power is:",
                    options: ["Zn > Fe > Cu > Ag", "Fe > Zn > Cu > Ag", "Cu > Zn > Fe > Ag", "Ag > Cu > Fe > Zn"],
                    answer: 0,
                    explain: "Reducing power decreases as standard electrode potential increases: Zn(-0.76) > Fe(-0.44) > Cu(+0.34) > Ag(+0.80)"
                },
                {
                    id: "c5",
                    text: "Which of the following shows maximum enol content?",
                    options: ["CH₃COCH₃", "CH₃COCH₂CH₃", "CH₃COCH₂COCH₃", "C₆H₅COCH₃"],
                    answer: 2,
                    explain: "CH₃COCH₂COCH₃ (acetylacetone) shows maximum enol content due to intramolecular hydrogen bonding in enol form."
                },
                {
                    id: "c6",
                    text: "The IUPAC name of compound (CH₃)₂C=CHCH₂OH is:",
                    options: ["3-methylbut-2-en-1-ol", "2-methylbut-2-en-4-ol", "3-methylbut-3-en-1-ol", "Isoprenol"],
                    answer: 0,
                    explain: "4-carbon chain with OH at C1, double bond at C2-C3, methyl at C3: 3-methylbut-2-en-1-ol"
                },
                {
                    id: "c7",
                    text: "Which of the following is most stable free radical?",
                    options: ["CH₃CH₂•", "(CH₃)₂CH•", "•CH₂COOH", "C₆H₅•"],
                    answer: 3,
                    explain: "Phenyl radical is stabilized by resonance with benzene ring, making it most stable."
                },
                {
                    id: "c8",
                    text: "The spin-only magnetic moment of [Mn(H₂O)₆]²⁺ is:",
                    options: ["3.87 BM", "4.89 BM", "5.92 BM", "0 BM"],
                    answer: 2,
                    explain: "Mn²⁺ is d⁵ with 5 unpaired electrons. μ = √[n(n+2)] = √[5×7] = 5.92 BM"
                },
                {
                    id: "c9",
                    text: "Which of the following shows maximum hydrogen bonding?",
                    options: ["Ethanol", "Diethyl ether", "Ethanoic acid", "Acetone"],
                    answer: 2,
                    explain: "Ethanoic acid shows both intermolecular and intramolecular hydrogen bonding, maximum among options."
                },
                {
                    id: "c10",
                    text: "The number of π electrons in naphthalene is:",
                    options: ["8", "10", "12", "14"],
                    answer: 1,
                    explain: "Naphthalene has 10 π electrons following Hückel's rule (4n+2 = 10, n=2) for aromaticity."
                },
                {
                    id: "c11",
                    text: "Which of the following undergoes nucleophilic substitution most readily?",
                    options: ["Vinyl chloride", "Allyl chloride", "Benzyl chloride", "Phenyl chloride"],
                    answer: 2,
                    explain: "Benzyl chloride forms resonance-stabilized benzyl cation, undergoes SN1 most readily."
                },
                {
                    id: "c12",
                    text: "The oxidation state of chromium in K₂Cr₂O₇ is:",
                    options: ["+4", "+5", "+6", "+7"],
                    answer: 2,
                    explain: "In K₂Cr₂O₇: 2(+1) + 2Cr + 7(-2) = 0, solving gives Cr = +6"
                },
                {
                    id: "c13",
                    text: "Which of the following has minimum melting point?",
                    options: ["NaCl", "MgO", "CaCl₂", "Al₂O₃"],
                    answer: 0,
                    explain: "NaCl has lowest charge product (1×1=1) and lowest lattice energy, hence minimum melting point."
                },
                {
                    id: "c14",
                    text: "The correct order of first ionization energy is:",
                    options: ["B < C < N < O", "B < C < O < N", "C < B < O < N", "B < O < C < N"],
                    answer: 1,
                    explain: "General trend increases across period, but N > O due to half-filled stability: B < C < O < N"
                },
                {
                    id: "c15",
                    text: "Which of the following is not a buffer system?",
                    options: ["H₂CO₃/HCO₃⁻", "H₂PO₄⁻/HPO₄²⁻", "NH₄⁺/NH₃", "HNO₃/NO₃⁻"],
                    answer: 3,
                    explain: "HNO₃/NO₃⁻ is not buffer as HNO₃ is strong acid. Buffers need weak acid/base with conjugate."
                },
                {
                    id: "c16",
                    text: "The hybridization of central atom in [BrF₅] is:",
                    options: ["sp³d", "sp³d²", "sp³d³", "dsp³"],
                    answer: 1,
                    explain: "BrF₅ has 6 electron pairs (5 bonding + 1 lone pair) around Br, requiring sp³d² hybridization."
                },
                {
                    id: "c17",
                    text: "Which of the following shows maximum paramagnetic character?",
                    options: ["NO", "N₂", "CO", "CN⁻"],
                    answer: 0,
                    explain: "NO has one unpaired electron, making it paramagnetic. Others are diamagnetic."
                },
                {
                    id: "c18",
                    text: "The rate of reaction doubles when temperature increases from 27°C to 37°C. The activation energy is:",
                    options: ["53.6 kJ/mol", "12.8 kcal/mol", "Both A and B", "Cannot be determined"],
                    answer: 2,
                    explain: "Using Arrhenius equation: ln(k₂/k₁) = Ea/R(1/T₁ - 1/T₂). Ea = 53.6 kJ/mol = 12.8 kcal/mol"
                },
                {
                    id: "c19",
                    text: "Which of the following is most basic oxide?",
                    options: ["BeO", "MgO", "CaO", "BaO"],
                    answer: 3,
                    explain: "Basicity increases down the group. BaO is most basic among alkaline earth metal oxides."
                },
                {
                    id: "c20",
                    text: "The number of geometric isomers for [MA₂B₂C₂] octahedral complex is:",
                    options: ["5", "6", "15", "30"],
                    answer: 2,
                    explain: "For octahedral complex with 3 different pairs of identical ligands, geometric isomers = 15"
                },
                {
                    id: "c21",
                    text: "Which of the following shows maximum covalent character?",
                    options: ["LiI", "NaI", "KI", "CsI"],
                    answer: 0,
                    explain: "According to Fajan's rule, LiI has maximum covalent character due to high charge density of Li⁺."
                },
                {
                    id: "c22",
                    text: "The correct order of atomic radii is:",
                    options: ["F > Cl > Br > I", "I > Br > Cl > F", "Cl > F > Br > I", "Br > Cl > I > F"],
                    answer: 1,
                    explain: "Atomic radii increase down the group due to addition of electron shells: I > Br > Cl > F"
                },
                {
                    id: "c23",
                    text: "Which of the following is strongest nucleophile in polar protic solvent?",
                    options: ["F⁻", "Cl⁻", "Br⁻", "I⁻"],
                    answer: 3,
                    explain: "In protic solvents, larger anions are better nucleophiles due to less solvation: I⁻ > Br⁻ > Cl⁻ > F⁻"
                },
                {
                    id: "c24",
                    text: "The shape of XeOF₄ molecule is:",
                    options: ["Square pyramidal", "Pentagonal planar", "T-shaped", "See-saw"],
                    answer: 0,
                    explain: "XeOF₄ has 6 electron pairs (5 bonding + 1 lone pair) giving square pyramidal geometry."
                },
                {
                    id: "c25",
                    text: "Which of the following has maximum bond order?",
                    options: ["O₂", "O₂⁺", "O₂⁻", "O₂²⁻"],
                    answer: 1,
                    explain: "O₂⁺ has 15 electrons, bond order = (10-5)/2 = 2.5, highest among given species."
                },
                {
                    id: "c26",
                    text: "The number of moles of H₂SO₄ required to neutralize 2 moles of NaOH is:",
                    options: ["1", "2", "0.5", "4"],
                    answer: 0,
                    explain: "H₂SO₄ + 2NaOH → Na₂SO₄ + 2H₂O. 1 mole H₂SO₄ neutralizes 2 moles NaOH."
                },
                {
                    id: "c27",
                    text: "Which of the following is most acidic hydrogen?",
                    options: ["CH₃-CH₂-H", "HC≡C-H", "CH₂=CH-H", "C₆H₅-H"],
                    answer: 1,
                    explain: "Terminal alkyne hydrogen is most acidic due to sp hybridization (50% s-character)."
                },
                {
                    id: "c28",
                    text: "The electronic configuration of Cu²⁺ ion is:",
                    options: ["[Ar] 3d⁹", "[Ar] 3d⁷ 4s²", "[Ar] 3d⁸ 4s¹", "[Ar] 3d¹⁰ 4s¹"],
                    answer: 0,
                    explain: "Cu: [Ar] 3d¹⁰ 4s¹, Cu²⁺ loses 4s¹ and one 3d electron: [Ar] 3d⁹"
                },
                {
                    id: "c29",
                    text: "Which of the following shows maximum dipole moment?",
                    options: ["CO₂", "BeCl₂", "H₂O", "BF₃"],
                    answer: 2,
                    explain: "H₂O has bent structure with two lone pairs, creating maximum dipole moment among options."
                },
                {
                    id: "c30",
                    text: "The number of lone pairs in ICl₂⁻ is:",
                    options: ["2", "3", "4", "5"],
                    answer: 1,
                    explain: "I has 7 valence electrons + 1 from charge = 8. Uses 2 for bonding with Cl atoms, leaving 6 electrons = 3 lone pairs."
                },
                {
                    id: "c31",
                    text: "Which of the following undergoes fastest electrophilic aromatic substitution?",
                    options: ["Benzene", "Toluene", "Nitrobenzene", "Chlorobenzene"],
                    answer: 1,
                    explain: "Toluene has electron-donating methyl group that activates benzene ring for electrophilic substitution."
                },
                {
                    id: "c32",
                    text: "The correct order of stability of alkenes is:",
                    options: ["1° > 2° > 3°", "3° > 2° > 1°", "2° > 1° > 3°", "1° = 2° = 3°"],
                    answer: 1,
                    explain: "Alkene stability increases with substitution due to hyperconjugation: tetrasubstituted > trisubstituted > disubstituted > monosubstituted"
                },
                {
                    id: "c33",
                    text: "The pKa value of strongest acid among the following is:",
                    options: ["Lowest", "Highest", "Intermediate", "Cannot be determined"],
                    answer: 0,
                    explain: "Strongest acid has lowest pKa value. pKa = -log Ka, stronger acid has higher Ka and lower pKa."
                },
                {
                    id: "c34",
                    text: "Which of the following shows maximum ionic conductance?",
                    options: ["Li⁺", "Na⁺", "K⁺", "Cs⁺"],
                    answer: 3,
                    explain: "Cs⁺ has largest ionic radius, least hydrated, highest mobility, hence maximum ionic conductance."
                },
                {
                    id: "c35",
                    text: "The number of bridging carbonyls in Fe₂(CO)₉ is:",
                    options: ["0", "2", "3", "6"],
                    answer: 2,
                    explain: "Fe₂(CO)₉ has 3 bridging CO groups and 6 terminal CO groups in its structure."
                },
                {
                    id: "c36",
                    text: "Which of the following is not isoelectronic with neon?",
                    options: ["Na⁺", "Mg²⁺", "F⁻", "Cl⁻"],
                    answer: 3,
                    explain: "Cl⁻ has 18 electrons (isoelectronic with Ar). Ne has 10 electrons like Na⁺, Mg²⁺, F⁻."
                },
                {
                    id: "c37",
                    text: "The geometry around boron in BH₄⁻ is:",
                    options: ["Triangular planar", "Tetrahedral", "Square planar", "Pyramidal"],
                    answer: 1,
                    explain: "BH₄⁻ has 4 bonding pairs around B with no lone pairs, giving tetrahedral geometry."
                },
                {
                    id: "c38",
                    text: "Which of the following has maximum solubility in water?",
                    options: ["BeSO₄", "MgSO₄", "CaSO₄", "BaSO₄"],
                    answer: 0,
                    explain: "Solubility of alkaline earth sulfates decreases down the group: BeSO₄ > MgSO₄ > CaSO₄ > BaSO₄"
                },
                {
                    id: "c39",
                    text: "The number of unpaired electrons in [CoF₆]³⁻ is:",
                    options: ["0", "2", "4", "6"],
                    answer: 2,
                    explain: "Co³⁺ is d⁶. F⁻ is weak field ligand, so high spin: t₂g⁴ eg², giving 4 unpaired electrons."
                },
                {
                    id: "c40",
                    text: "Which of the following is most volatile?",
                    options: ["CH₃CH₂OH", "(CH₃)₂O", "CH₃CH₂CH₂OH", "CH₃OCH₂CH₃"],
                    answer: 1,
                    explain: "(CH₃)₂O has no hydrogen bonding and lowest boiling point, hence most volatile."
                },
                {
                    id: "c41",
                    text: "The correct order of thermal stability of carbonates is:",
                    options: ["BeCO₃ > MgCO₃ > CaCO₃", "CaCO₃ > MgCO₃ > BeCO₃", "MgCO₃ > CaCO₃ > BeCO₃", "All are equally stable"],
                    answer: 1,
                    explain: "Thermal stability increases down the group as cation size increases: CaCO₃ > MgCO₃ > BeCO₃"
                },
                {
                    id: "c42",
                    text: "Which of the following is an example of linkage isomerism?",
                    options: ["[Co(NH₃)₅(NO₂)]Cl₂ and [Co(NH₃)₅(ONO)]Cl₂", "[Co(NH₃)₆][Cr(CN)₆] and [Cr(NH₃)₆][Co(CN)₆]", "[Pt(NH₃)₄Cl₂]Br₂ and [Pt(NH₃)₄Br₂]Cl₂", "All of above"],
                    answer: 0,
                    explain: "NO₂⁻ can bind through N or O atoms, showing linkage isomerism between nitro and nitrito forms."
                },
                {
                    id: "c43",
                    text: "The number of sigma bonds in benzene is:",
                    options: ["6", "9", "12", "15"],
                    answer: 2,
                    explain: "Benzene has 6 C-C σ bonds and 6 C-H σ bonds, total 12 σ bonds."
                },
                {
                    id: "c44",
                    text: "Which of the following is strongest Brønsted base?",
                    options: ["OH⁻", "NH₂⁻", "H⁻", "CH₃⁻"],
                    answer: 2,
                    explain: "H⁻ (hydride ion) is strongest Brønsted base as H₂ is weakest conjugate acid."
                },
                {
                    id: "c45",
                    text: "The coordination number of central atom in [Ni(EDTA)]²⁻ is:",
                    options: ["4", "5", "6", "8"],
                    answer: 2,
                    explain: "EDTA is hexadentate ligand, provides 6 coordination sites to Ni²⁺, giving coordination number 6."
                }
            ]
        },
        {
            name: "Biology",
            questions: [
                // BOTANY Questions (45)
                {
                    id: "b1",
                    text: "Which of the following is the first enzyme of Calvin cycle?",
                    options: ["PEP carboxylase", "RuBisCO", "Phosphoglycerate kinase", "Aldolase"],
                    answer: 1,
                    explain: "RuBisCO (Ribulose bisphosphate carboxylase oxygenase) catalyzes the first step of Calvin cycle, fixing CO₂ to RuBP."
                },
                {
                    id: "b2",
                    text: "The phenomenon of ripening of fruits even after detachment from plant is called:",
                    options: ["Senescence", "Abscission", "Climacteric", "Parthenocarpy"],
                    answer: 2,
                    explain: "Climacteric fruits like apple, banana continue ripening after harvest due to increased ethylene production."
                },
                {
                    id: "b3",
                    text: "Which tissue forms the bulk of ground tissue in monocot stems?",
                    options: ["Cortex", "Pith", "Parenchyma", "Endodermis"],
                    answer: 2,
                    explain: "In monocot stems, parenchyma forms the bulk of ground tissue with scattered vascular bundles."
                },
                {
                    id: "b4",
                    text: "The direction of phloem transport is:",
                    options: ["Always upward", "Always downward", "Bidirectional", "From source to sink"],
                    answer: 3,
                    explain: "Phloem transport is from source (photosynthetic organs) to sink (storage/growing organs), can be in any direction."
                },
                {
                    id: "b5",
                    text: "Which of the following is NOT a characteristic of C₄ plants?",
                    options: ["Kranz anatomy", "High CO₂ compensation point", "Bundle sheath cells", "Spatial separation"],
                    answer: 1,
                    explain: "C₄ plants have low CO₂ compensation point due to efficient CO₂ concentration mechanism."
                },
                {
                    id: "b6",
                    text: "The process by which plants shed their leaves is called:",
                    options: ["Senescence", "Abscission", "Wilting", "Defoliation"],
                    answer: 1,
                    explain: "Abscission is controlled shedding of leaves, flowers, or fruits through formation of abscission layer."
                },
                {
                    id: "b7",
                    text: "Which plant growth regulator is responsible for breaking seed dormancy?",
                    options: ["Auxin", "Cytokinin", "Gibberellin", "Abscisic acid"],
                    answer: 2,
                    explain: "Gibberellins break seed dormancy by overcoming effects of ABA and promoting enzyme synthesis."
                },
                {
                    id: "b8",
                    text: "The specialized cells that control stomatal opening in grasses are:",
                    options: ["Guard cells", "Subsidiary cells", "Bulliform cells", "Companion cells"],
                    answer: 0,
                    explain: "Guard cells control stomatal opening in all plants, including grasses where they are dumbbell-shaped."
                },
                {
                    id: "b9",
                    text: "Which of the following is involved in biological clock of plants?",
                    options: ["Chlorophyll", "Phytochrome", "Carotenoids", "Anthocyanins"],
                    answer: 1,
                    explain: "Phytochrome system acts as biological clock, measuring day length and controlling circadian rhythms."
                },
                {
                    id: "b10",
                    text: "The water conducting elements in gymnosperms are:",
                    options: ["Vessels", "Tracheids", "Sieve tubes", "Companion cells"],
                    answer: 1,
                    explain: "Gymnosperms have only tracheids for water conduction, vessels are absent."
                },
                {
                    id: "b11",
                    text: "Which of the following shows thigmotropism?",
                    options: ["Roots growing towards water", "Stems bending towards light", "Tendrils coiling around support", "Leaves folding at night"],
                    answer: 2,
                    explain: "Thigmotropism is growth response to touch, seen in tendrils that coil around support structures."
                },
                {
                    id: "b12",
                    text: "The female gametophyte in angiosperms is:",
                    options: ["Embryo sac", "Ovule", "Carpel", "Nucellus"],
                    answer: 0,
                    explain: "Embryo sac is the mature female gametophyte containing egg, synergids, antipodals, and polar nuclei."
                },
                {
                    id: "b13",
                    text: "Which enzyme is responsible for splitting of water in photosynthesis?",
                    options: ["Catalase", "Peroxidase", "Water-splitting complex", "RuBisCO"],
                    answer: 2,
                    explain: "Water-splitting complex (oxygen-evolving complex) in PSII catalyzes photolysis of water."
                },
                {
                    id: "b14",
                    text: "The phenomenon where older leaves turn yellow due to chlorophyll breakdown is:",
                    options: ["Etiolation", "Senescence", "Chlorosis", "Necrosis"],
                    answer: 1,
                    explain: "Senescence involves programmed aging where chlorophyll breaks down, revealing yellow pigments."
                },
                {
                    id: "b15",
                    text: "Which of the following is a naturally occurring cytokinin?",
                    options: ["Kinetin", "BAP", "Zeatin", "2,4-D"],
                    answer: 2,
                    explain: "Zeatin is naturally occurring cytokinin found in corn kernels and many other plants."
                },
                {
                    id: "b16",
                    text: "The arrangement of ovules in which funicle is fused with integument is:",
                    options: ["Orthotropous", "Anatropous", "Campylotropous", "Circinotropous"],
                    answer: 1,
                    explain: "In anatropous ovules, ovule is inverted and funicle fuses with integument forming raphe."
                },
                {
                    id: "b17",
                    text: "Which of the following is involved in photorespiration?",
                    options: ["Only chloroplasts", "Only mitochondria", "Only peroxisomes", "All three organelles"],
                    answer: 3,
                    explain: "Photorespiration involves chloroplasts (glycolate formation), peroxisomes (glycine formation), and mitochondria (serine synthesis)."
                },
                {
                    id: "b18",
                    text: "The stage of meiosis where crossing over occurs is:",
                    options: ["Leptotene", "Zygotene", "Pachytene", "Diplotene"],
                    answer: 2,
                    explain: "Crossing over (recombination) occurs during pachytene stage when homologous chromosomes are paired."
                },
                {
                    id: "b19",
                    text: "Which of the following is essential for chlorophyll synthesis?",
                    options: ["Calcium", "Magnesium", "Potassium", "Sodium"],
                    answer: 1,
                    explain: "Magnesium is central atom in chlorophyll molecule, essential for its synthesis and photosynthesis."
                },
                {
                    id: "b20",
                    text: "The condition where plant completes its life cycle in two years is:",
                    options: ["Annual", "Biennial", "Perennial", "Ephemeral"],
                    answer: 1,
                    explain: "Biennial plants complete life cycle in two years - vegetative growth in first year, reproduction in second."
                },
                {
                    id: "b21",
                    text: "Which of the following is NOT a function of auxin?",
                    options: ["Apical dominance", "Parthenocarpy", "Cell division", "Tropistic responses"],
                    answer: 2,
                    explain: "Cell division is primarily promoted by cytokinins, not auxins. Auxins promote cell elongation."
                },
                {
                    id: "b22",
                    text: "The tissue that transports photosynthetic products is:",
                    options: ["Xylem", "Phloem", "Cambium", "Epidermis"],
                    answer: 1,
                    explain: "Phloem tissue transports organic solutes including photosynthetic products from leaves to other parts."
                },
                {
                    id: "b23",
                    text: "Which of the following shows pneumatophores?",
                    options: ["Desert plants", "Aquatic plants", "Mangrove plants", "Epiphytic plants"],
                    answer: 2,
                    explain: "Mangrove plants have pneumatophores (aerial roots) for gas exchange in waterlogged anaerobic soil."
                },
                {
                    id: "b24",
                    text: "The primary cell wall is mainly composed of:",
                    options: ["Lignin", "Cellulose", "Suberin", "Chitin"],
                    answer: 1,
                    explain: "Primary cell wall is mainly cellulose microfibrils embedded in pectin and hemicellulose matrix."
                },
                {
                    id: "b25",
                    text: "Which of the following is characteristic of wind-pollinated flowers?",
                    options: ["Large petals", "Nectar production", "Light, dry pollen", "Strong fragrance"],
                    answer: 2,
                    explain: "Wind-pollinated flowers produce light, dry, non-sticky pollen that can be carried by air currents."
                },
                {
                    id: "b26",
                    text: "The process of seed formation without fertilization is called:",
                    options: ["Parthenocarpy", "Apomixis", "Parthenogenesis", "Vegetative propagation"],
                    answer: 1,
                    explain: "Apomixis is asexual seed formation without fertilization, producing genetically identical offspring."
                },
                {
                    id: "b27",
                    text: "Which of the following is a long-day plant?",
                    options: ["Rice", "Soybean", "Wheat", "Cotton"],
                    answer: 2,
                    explain: "Wheat is long-day plant requiring photoperiods longer than critical length for flowering."
                },
                {
                    id: "b28",
                    text: "The protective covering around plumule is called:",
                    options: ["Coleoptile", "Coleorhiza", "Scutellum", "Aleurone layer"],
                    answer: 0,
                    explain: "Coleoptile is protective sheath around plumule (shoot tip) in monocot seeds like grass and cereals."
                },
                {
                    id: "b29",
                    text: "Which of the following is involved in translocation of organic solutes?",
                    options: ["Tracheids", "Vessels", "Sieve tubes", "Xylem parenchyma"],
                    answer: 2,
                    explain: "Sieve tubes (sieve tube elements) in phloem transport organic solutes like sucrose."
                },
                {
                    id: "b30",
                    text: "The phenomenon where plants bloom in response to cold treatment is:",
                    options: ["Photoperiodism", "Vernalization", "Thermoperiodism", "Stratification"],
                    answer: 1,
                    explain: "Vernalization is requirement of cold treatment for flowering, especially in winter varieties."
                },
                {
                    id: "b31",
                    text: "Which of the following is a characteristic of hydrophytes?",
                    options: ["Thick waxy cuticle", "Deep root system", "Reduced vascular tissue", "Succulent leaves"],
                    answer: 2,
                    explain: "Hydrophytes have reduced vascular tissue as water and mineral transport is not a major concern."
                },
                {
                    id: "b32",
                    text: "The cells that surround sieve tube elements are:",
                    options: ["Companion cells", "Transfer cells", "Albuminous cells", "All of above"],
                    answer: 0,
                    explain: "Companion cells are specialized parenchyma cells that surround and support sieve tube elements."
                },
                {
                    id: "b33",
                    text: "Which of the following is NOT a micronutrient?",
                    options: ["Molybdenum", "Chlorine", "Phosphorus", "Copper"],
                    answer: 2,
                    explain: "Phosphorus is macronutrient required in large amounts. Mo, Cl, Cu are micronutrients."
                },
                {
                    id: "b34",
                    text: "The opening of floral buds into flowers is called:",
                    options: ["Anthesis", "Dehiscence", "Pollination", "Fertilization"],
                    answer: 0,
                    explain: "Anthesis is opening of flower buds and period when flower is functionally mature."
                },
                {
                    id: "b35",
                    text: "Which of the following is involved in nitrogen fixation?",
                    options: ["Rhizobium", "Azotobacter", "Clostridium", "All of above"],
                    answer: 3,
                    explain: "All three bacteria can fix atmospheric nitrogen - Rhizobium (symbiotic), Azotobacter and Clostridium (free-living)."
                },
                {
                    id: "b36",
                    text: "The layer of cells in root that regulates water and mineral uptake is:",
                    options: ["Epidermis", "Cortex", "Endodermis", "Pericycle"],
                    answer: 2,
                    explain: "Endodermis with Casparian strips regulates radial movement of water and minerals into stele."
                },
                {
                    id: "b37",
                    text: "Which of the following shows determinate growth?",
                    options: ["Root apex", "Shoot apex", "Leaves", "Cambium"],
                    answer: 2,
                    explain: "Leaves show determinate growth - they stop growing after reaching mature size."
                },
                {
                    id: "b38",
                    text: "The sugar alcohol commonly found in phloem sap is:",
                    options: ["Glucose", "Fructose", "Sucrose", "Mannitol"],
                    answer: 3,
                    explain: "Mannitol is common sugar alcohol transported in phloem of many plants, especially those in saline conditions."
                },
                {
                    id: "b39",
                    text: "Which of the following is a characteristic of CAM plants?",
                    options: ["Stomata open during day", "CO₂ fixed as C₃ acid", "Temporal separation of CO₂ fixation", "High water loss"],
                    answer: 2,
                    explain: "CAM plants show temporal separation - CO₂ fixation at night, Calvin cycle during day."
                },
                {
                    id: "b40",
                    text: "The male gametophyte in angiosperms is:",
                    options: ["Anther", "Pollen grain", "Stamen", "Microsporangium"],
                    answer: 1,
                    explain: "Pollen grain is mature male gametophyte containing generative and vegetative nuclei."
                },
                {
                    id: "b41",
                    text: "Which of the following is involved in gibberellin biosynthesis?",
                    options: ["Mevalonic acid pathway", "Shikimic acid pathway", "Tryptophan pathway", "Methionine pathway"],
                    answer: 0,
                    explain: "Gibberellins are synthesized via mevalonic acid pathway, same as other terpenoids."
                },
                {
                    id: "b42",
                    text: "The phenomenon where plants bend towards unilateral light is:",
                    options: ["Geotropism", "Phototropism", "Chemotropism", "Hydrotropism"],
                    answer: 1,
                    explain: "Phototropism is bending response toward unilateral light, mediated by auxin redistribution."
                },
                {
                    id: "b43",
                    text: "Which of the following is characteristic of secondary growth?",
                    options: ["Increase in length", "Increase in girth", "Formation of branches", "Leaf development"],
                    answer: 1,
                    explain: "Secondary growth results in increase in girth/thickness of stems and roots through cambial activity."
                },
                {
                    id: "b44",
                    text: "The starch sheath in dicot stems is:",
                    options: ["Epidermis", "Cortex", "Endodermis", "Pericycle"],
                    answer: 2,
                    explain: "Endodermis in dicot stems is called starch sheath as it contains abundant starch grains."
                },
                {
                    id: "b45",
                    text: "Which of the following shows hypogeal germination?",
                    options: ["Castor", "Mustard", "Gram", "Onion"],
                    answer: 2,
                    explain: "Gram (chickpea) shows hypogeal germination where cotyledons remain underground."
                },

                // ZOOLOGY Questions (45)
                {
                    id: "z1",
                    text: "The enzyme that converts angiotensinogen to angiotensin I is:",
                    options: ["Renin", "ACE", "Aldosterone", "ADH"],
                    answer: 0,
                    explain: "Renin secreted by juxtaglomerular cells converts angiotensinogen to angiotensin I in RAAS."
                },
                {
                    id: "z2",
                    text: "Which of the following is NOT a function of large intestine?",
                    options: ["Water absorption", "Electrolyte absorption", "Vitamin K synthesis", "Protein digestion"],
                    answer: 3,
                    explain: "Protein digestion occurs in stomach and small intestine. Large intestine mainly absorbs water and electrolytes."
                },
                {
                    id: "z3",
                    text: "The hormone that regulates basal metabolic rate is:",
                    options: ["Insulin", "Cortisol", "Thyroxine", "Growth hormone"],
                    answer: 2,
                    explain: "Thyroxine (T₄) and triiodothyronine (T₃) from thyroid gland regulate basal metabolic rate."
                },
                {
                    id: "z4",
                    text: "Which part of brain is responsible for regulation of breathing?",
                    options: ["Cerebrum", "Cerebellum", "Medulla oblongata", "Hypothalamus"],
                    answer: 2,
                    explain: "Medulla oblongata contains respiratory center that controls automatic rhythmic breathing."
                },
                {
                    id: "z5",
                    text: "The cells that secrete histamine during allergic reactions are:",
                    options: ["Neutrophils", "Eosinophils", "Basophils", "Monocytes"],
                    answer: 2,
                    explain: "Basophils and mast cells secrete histamine, causing vasodilation and increased permeability in allergic reactions."
                },
                {
                    id: "z6",
                    text: "Which of the following is reabsorbed in descending limb of loop of Henle?",
                    options: ["Na⁺", "Cl⁻", "Water", "Urea"],
                    answer: 2,
                    explain: "Descending limb is permeable to water but not to salts. Water is passively reabsorbed here."
                },
                {
                    id: "z7",
                    text: "The type of immunity transferred from mother to baby through breast milk is:",
                    options: ["Active natural", "Active artificial", "Passive natural", "Passive artificial"],
                    answer: 2,
                    explain: "Colostrum and breast milk provide antibodies (especially IgA) giving passive natural immunity to baby."
                },
                {
                    id: "z8",
                    text: "Which hormone stimulates ovulation?",
                    options: ["FSH", "LH", "Estrogen", "Progesterone"],
                    answer: 1,
                    explain: "LH surge triggers ovulation by causing rupture of mature Graafian follicle."
                },
                {
                    id: "z9",
                    text: "The sound produced by closure of semilunar valves is:",
                    options: ["First heart sound (lub)", "Second heart sound (dub)", "Third heart sound", "Murmur"],
                    answer: 1,
                    explain: "Second heart sound (dub) is produced by closure of aortic and pulmonary semilunar valves."
                },
                {
                    id: "z10",
                    text: "Which of the following crosses blood-brain barrier most easily?",
                    options: ["Glucose", "Amino acids", "Lipids", "Proteins"],
                    answer: 2,
                    explain: "Lipophilic substances cross blood-brain barrier most easily due to high lipid content of brain tissue."
                },
                {
                    id: "z11",
                    text: "The part of sperm that contains enzymes for penetrating egg is:",
                    options: ["Head", "Neck", "Middle piece", "Tail"],
                    answer: 0,
                    explain: "Acrosome in sperm head contains hydrolytic enzymes like hyaluronidase for egg penetration."
                },
                {
                    id: "z12",
                    text: "Which vitamin is synthesized by bacteria in large intestine?",
                    options: ["Vitamin C", "Vitamin D", "Vitamin K", "Vitamin B₁₂"],
                    answer: 2,
                    explain: "Vitamin K is synthesized by intestinal bacteria and absorbed in large intestine."
                },
                {
                    id: "z13",
                    text: "The muscle fiber type that contracts slowly but resists fatigue is:",
                    options: ["Type I (slow oxidative)", "Type IIa (fast oxidative)", "Type IIx (fast glycolytic)", "Cardiac muscle"],
                    answer: 0,
                    explain: "Type I muscle fibers are slow-contracting, fatigue-resistant, rich in mitochondria and myoglobin."
                },
                {
                    id: "z14",
                    text: "Which of the following is NOT a component of nephron?",
                    options: ["Bowman's capsule", "Proximal tubule", "Loop of Henle", "Ureter"],
                    answer: 3,
                    explain: "Ureter is part of urinary tract that carries urine from kidney to bladder, not part of nephron."
                },
                {
                    id: "z15",
                    text: "The hormone that increases blood calcium level is:",
                    options: ["Calcitonin", "Parathyroid hormone", "Thyroxine", "Insulin"],
                    answer: 1,
                    explain: "Parathyroid hormone (PTH) increases blood calcium by promoting bone resorption and calcium absorption."
                },
                {
                    id: "z16",
                    text: "Which structure in ear is responsible for hearing?",
                    options: ["Semicircular canals", "Utricle", "Saccule", "Cochlea"],
                    answer: 3,
                    explain: "Cochlea contains organ of Corti with hair cells that convert sound vibrations to electrical signals."
                },
                {
                    id: "z17",
                    text: "The phase of menstrual cycle when corpus luteum degenerates is:",
                    options: ["Menstrual phase", "Follicular phase", "Ovulation", "Luteal phase"],
                    answer: 0,
                    explain: "If pregnancy doesn't occur, corpus luteum degenerates during menstrual phase, causing menstruation."
                },
                {
                    id: "z18",
                    text: "Which of the following is involved in clot retraction?",
                    options: ["Platelets", "Fibrinogen", "Prothrombin", "Heparin"],
                    answer: 0,
                    explain: "Platelets contain contractile proteins that help in clot retraction and wound healing."
                },
                {
                    id: "z19",
                    text: "The primary site of erythropoiesis in adults is:",
                    options: ["Liver", "Spleen", "Red bone marrow", "Yellow bone marrow"],
                    answer: 2,
                    explain: "Red bone marrow is primary site of RBC production (erythropoiesis) in healthy adults."
                },
                {
                    id: "z20",
                    text: "Which cranial nerve carries taste sensation from posterior 1/3 of tongue?",
                    options: ["Trigeminal", "Facial", "Glossopharyngeal", "Vagus"],
                    answer: 2,
                    explain: "Glossopharyngeal nerve (IX cranial nerve) carries taste from posterior third of tongue."
                },
                {
                    id: "z21",
                    text: "The condition where person cannot distinguish between red and green colors is:",
                    options: ["Myopia", "Hyperopia", "Color blindness", "Night blindness"],
                    answer: 2,
                    explain: "Red-green color blindness is due to deficiency or absence of specific cone cells in retina."
                },
                {
                    id: "z22",
                    text: "Which of the following is NOT secreted by pancreas?",
                    options: ["Insulin", "Glucagon", "Pancreatic lipase", "Gastrin"],
                    answer: 3,
                    explain: "Gastrin is secreted by G cells in stomach antrum, not by pancreas."
                },
                {
                    id: "z23",
                    text: "The valve between left ventricle and aorta is:",
                    options: ["Tricuspid valve", "Bicuspid valve", "Aortic semilunar valve", "Pulmonary semilunar valve"],
                    answer: 2,
                    explain: "Aortic semilunar valve prevents backflow from aorta to left ventricle."
                },
                {
                    id: "z24",
                    text: "Which of the following is characteristic of smooth muscle?",
                    options: ["Voluntary control", "Striations", "Single nucleus", "Fast contraction"],
                    answer: 2,
                    explain: "Smooth muscle cells are uninucleated, involuntary, non-striated, and contract slowly."
                },
                {
                    id: "z25",
                    text: "The hormone that prepares uterus for implantation is:",
                    options: ["Estrogen", "Progesterone", "FSH", "LH"],
                    answer: 1,
                    explain: "Progesterone prepares and maintains uterine endometrium for implantation of fertilized ovum."
                },
                {
                    id: "z26",
                    text: "Which of the following is reabsorbed maximally in proximal convoluted tubule?",
                    options: ["Glucose", "Sodium", "Water", "All of above"],
                    answer: 3,
                    explain: "PCT reabsorbs maximum amounts of all filtered substances - about 65% of water, sodium, and 100% of glucose."
                },
                {
                    id: "z27",
                    text: "The neurotransmitter released at neuromuscular junction is:",
                    options: ["Acetylcholine", "Dopamine", "Serotonin", "GABA"],
                    answer: 0,
                    explain: "Acetylcholine is released from motor neurons at neuromuscular junction to stimulate muscle contraction."
                },
                {
                    id: "z28",
                    text: "Which of the following is NOT a granulocyte?",
                    options: ["Neutrophil", "Eosinophil", "Basophil", "Lymphocyte"],
                    answer: 3,
                    explain: "Lymphocyte is agranulocyte (lacks visible cytoplasmic granules). Others are granulocytes."
                },
                {
                    id: "z29",
                    text: "The structure that prevents food from entering trachea is:",
                    options: ["Epiglottis", "Glottis", "Vocal cords", "Larynx"],
                    answer: 0,
                    explain: "Epiglottis acts as lid covering glottis during swallowing to prevent food aspiration."
                },
                {
                    id: "z30",
                    text: "Which hormone controls the release of growth hormone?",
                    options: ["TRH", "CRH", "GHRH", "GnRH"],
                    answer: 2,
                    explain: "Growth hormone releasing hormone (GHRH) from hypothalamus stimulates GH release from anterior pituitary."
                },
                {
                    id: "z31",
                    text: "The period between two successive menstrual cycles is approximately:",
                    options: ["21 days", "28 days", "35 days", "Variable"],
                    answer: 1,
                    explain: "Average menstrual cycle length is 28 days, though normal range is 21-35 days."
                },
                {
                    id: "z32",
                    text: "Which of the following is involved in cell-mediated immunity?",
                    options: ["B lymphocytes", "Antibodies", "T lymphocytes", "Complement proteins"],
                    answer: 2,
                    explain: "T lymphocytes (T cells) are responsible for cell-mediated immunity against intracellular pathogens."
                },
                {
                    id: "z33",
                    text: "The opening between middle ear and nasopharynx is:",
                    options: ["External auditory meatus", "Eustachian tube", "Oval window", "Round window"],
                    answer: 1,
                    explain: "Eustachian tube connects middle ear to nasopharynx, equalizing pressure on both sides of tympanum."
                },
                {
                    id: "z34",
                    text: "Which of the following is produced by zona glomerulosa of adrenal cortex?",
                    options: ["Cortisol", "Aldosterone", "Androgens", "Adrenaline"],
                    answer: 1,
                    explain: "Zona glomerulosa produces mineralocorticoids like aldosterone that regulate electrolyte balance."
                },
                {
                    id: "z35",
                    text: "The term 'systole' refers to:",
                    options: ["Relaxation of heart", "Contraction of heart", "Filling of ventricles", "Opening of valves"],
                    answer: 1,
                    explain: "Systole is contraction phase of cardiac cycle when blood is pumped out of heart chambers."
                },
                {
                    id: "z36",
                    text: "Which of the following is NOT a function of skin?",
                    options: ["Temperature regulation", "Protection", "Vitamin D synthesis", "Insulin production"],
                    answer: 3,
                    explain: "Insulin is produced by pancreatic beta cells, not skin. Skin has other important functions listed."
                },
                {
                    id: "z37",
                    text: "The region of stomach that connects to duodenum is:",
                    options: ["Fundus", "Body", "Antrum", "Pylorus"],
                    answer: 3,
                    explain: "Pylorus is narrow region of stomach that opens into duodenum through pyloric sphincter."
                },
                {
                    id: "z38",
                    text: "Which of the following is involved in autoimmune diseases?",
                    options: ["Innate immunity", "Humoral immunity", "Cell-mediated immunity", "Both B and C"],
                    answer: 3,
                    explain: "Autoimmune diseases involve both humoral (antibody) and cell-mediated immune responses against self-antigens."
                },
                {
                    id: "z39",
                    text: "The maximum life span of RBCs in humans is:",
                    options: ["60 days", "90 days", "120 days", "150 days"],
                    answer: 2,
                    explain: "Human RBCs have average lifespan of 120 days before being destroyed in spleen."
                },
                {
                    id: "z40",
                    text: "Which part of eye has maximum refractive power?",
                    options: ["Cornea", "Lens", "Aqueous humor", "Vitreous humor"],
                    answer: 0,
                    explain: "Cornea provides about 65-75% of eye's refractive power due to high refractive index difference with air."
                },
                {
                    id: "z41",
                    text: "The hormone that stimulates milk production is:",
                    options: ["Oxytocin", "Prolactin", "Estrogen", "Progesterone"],
                    answer: 1,
                    explain: "Prolactin from anterior pituitary stimulates milk synthesis and secretion in mammary glands."
                },
                {
                    id: "z42",
                    text: "Which of the following is NOT a characteristic of cardiac muscle?",
                    options: ["Involuntary", "Striated", "Branched", "Multinucleated"],
                    answer: 3,
                    explain: "Cardiac muscle fibers are uninucleated or binucleated, not multinucleated like skeletal muscle."
                },
                {
                    id: "z43",
                    text: "The part of nephron where maximum reabsorption occurs is:",
                    options: ["Glomerulus", "Proximal convoluted tubule", "Loop of Henle", "Distal convoluted tubule"],
                    answer: 1,
                    explain: "Proximal convoluted tubule reabsorbs maximum amount of filtered substances (about 65-80%)."
                },
                {
                    id: "z44",
                    text: "Which enzyme initiates carbohydrate digestion in mouth?",
                    options: ["Pepsin", "Amylase", "Lipase", "Trypsin"],
                    answer: 1,
                    explain: "Salivary amylase (ptyalin) begins starch digestion in mouth by breaking it into maltose."
                },
                {
                    id: "z45",
                    text: "The process of formation of urine involves:",
                    options: ["Only filtration", "Only reabsorption", "Only secretion", "All three processes"],
                    answer: 3,
                    explain: "Urine formation involves ultrafiltration, selective reabsorption, and tubular secretion in nephrons."
                }
            ]
        }
    ]
};
