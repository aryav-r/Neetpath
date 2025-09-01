// mock-test-5.js – NEET Mock Test 5 (High-Quality, All New Questions)
// Physics: 45 | Chemistry: 45 | Biology (Botany 45 + Zoology 45) = 180 questions total

window.MOCK_TEST_5 = {
  id: "neet-005",
  title: "Full Syllabus Mock 5",
  level: "hard",
  durationMinutes: 180,
  sections: [
    {
      name: "Physics",
      questions: [
        {
          id: "p1",
          text: "A block slides without friction on a wedge of angle θ. If block starts from rest at height H, its speed at bottom is:",
          options: ["√(2gH)", "√(gH)", "√(2gH cosθ)", "√(2gH sinθ)"],
          answer: 0,
          explain: "Energy conservation: potential mgh → kinetic ½mv², so v = √(2gH). Angle irrelevant when no friction."
        },
        {
          id: "p2",
          text: "A mass on a spring oscillates. The maximum speed occurs when displacement is:",
          options: ["0", "A", "A/2", "A/√2"],
          answer: 0,
          explain: "Max speed v_max = ωA occurs at equilibrium (x=0) where all energy is kinetic."
        },
        {
          id: "p3",
          text: "Two capacitors C and 4C are charged to same voltage V. They are disconnected and connected in parallel positive to positive. Final common voltage is:",
          options: ["V", "2V/5", "4V/5", "3V/5"],
          answer: 2,
          explain: "Total initial charge = V(C+4C)=5CV. Equivalent capacitance = 5C. Vf=Qtotal/5C=V"
        },
        {
          id: "p4",
          text: "A thin lens of focal length f produces three times magnified real image. Object distance is:",
          options: ["f/2", "f",  "3f/2", "2f/3"],
          answer: 2,
          explain: "Magnification m = v/u =3 and lens eqn:1/f=1/u+1/v=1/u+1/(3u)=4/(3u)⇒u=4f/3⇒v=3u=4f"
        },
        {
          id: "p5",
          text: "In a cyclotron, the mass of particle is corrected for relativistic effect when v approaches c. This limits maximum energy because:",
          options: ["Period increases", "Period decreases", "Radius decreases", "Magnetic field fails"],
          answer: 0,
          explain: "As mass increases relativistically, cyclotron frequency drops and particle falls out of phase."
        },
        {
          id: "p6",
          text: "A point charge q moves with constant velocity v. The magnetic field at distance r perpendicular to v is:",
          options: ["μ₀qv/(2πr)", "μ₀qv/(4πr²)", "μ₀qv²/(2πr)", "μ₀qv/(2r)"],
          answer: 0,
          explain: "Magnetic field of moving charge B = μ₀qv/(4πr²)×(×2πr) if integrated over circle gives μ₀qv/(2πr)."
        },
        {
          id: "p7",
          text: "A current I flows in a loop and produces magnetic dipole moment μ. The torque in uniform field B is:",
          options: ["μB sinθ", "μB cosθ", "μB tanθ", "μB"],
          answer: 0,
          explain: "Torque τ = μ×B = μB sinθ where θ is angle between μ and B."
        },
        {
          id: "p8",
          text: "An inductor of inductance L stores energy U. Current is:",
          options: ["√(2U/L)", "U/L", "2UL", "√(U/2L)"],
          answer: 0,
          explain: "Energy in inductor U = ½LI² ⇒ I = √(2U/L)."
        },
        {
          id: "p9",
          text: "A DC current produces steady magnetic field. The field lines inside current‐carrying conductor are:",
          options: ["Concentric circles", "Radial lines", "Straight lines", "Helical"],
          answer: 0,
          explain: "Inside DC‐carrying wire, magnetic field lines form concentric circles around axis (Ampère’s law)."
        },
        {
          id: "p10",
          text: "A coil with N turns encloses area A. If flux changes at rate dΦ/dt, induced emf is:",
          options: ["NdΦ/dt", "-NdΦ/dt", "Nd²Φ/dt²", "Φ/(Ndt)"],
          answer: 1,
          explain: "Faraday’s law: emf = -N dΦ/dt (negative sign from Lenz’s law)."
        },
        {
          id: "p11",
          text: "In double slit experiment, slit separation is d. To double fringe width β, screen distance D must be:",
          options: ["D/2", "D", "2D", "4D"],
          answer: 2,
          explain: "Fringe width β = λD/d. To double β, double D."
        },
        {
          id: "p12",
          text: "A photon and electron have equal momentum p. Ratio of their kinetic energies is:",
          options: ["m_ec²/p : p²/(2m_e)", "pc : p²/(2m_e)", "p²/(2m_e c²)", "pc²"],
          answer: 1,
          explain: "Electron KE = p²/(2m), photon E = pc. Ratio photon/electron = pc / [p²/(2m)] = 2m c / p."
        },
        {
          id: "p13",
          text: "A thin rod rotates about center with angular speed ω. Its moment of inertia about axis at one end is:",
          options: ["1/3 mL²", "1/12 mL²", "2/3 mL²", "1/4 mL²"],
          answer: 2,
          explain: "Rod about center I_c = 1/12mL². Use parallel axis: I_end = I_c + m(L/2)² = 1/12mL² + 1/4mL² = 1/3mL²."
        },
        {
          id: "p14",
          text: "In Young’s modulus Y = stress/strain, units are:",
          options: ["N/m²", "Pa", "Both A and B", "N m"],
          answer: 2,
          explain: "Stress/strain units = (N/m²)/(dimensionless) = N/m² = Pascal."
        },
        {
          id: "p15",
          text: "A wire of resistivity ρ is bent into coil of radius R. Its inductance is L. If wire is straightened to same length, inductance:",
          options: ["Increases", "Decreases", "Same", "Zero"],
          answer: 1,
          explain: "Coil has self-inductance due to turns. Straight wire has negligible inductance relative to coil."
        },
        {
          id: "p16",
          text: "The relativistic kinetic energy of particle is KE = (γ-1)mc². At v = 0.6c, γ = 1.25. KE/mc² is:",
          options: ["0.25", "0.5", "1", "1.25"],
          answer: 0,
          explain: "γ - 1 = 0.25, so KE = 0.25 mc²."
        },
        {
          id: "p17",
          text: "Which of the following is NOT a Maxwell equation in vacuum?",
          options: ["∇·E = ρ/ε₀", "∇·B = 0", "∇×E = -∂B/∂t", "∇×B = μ₀ε₀∂E/∂t"],
          answer: 0,
          explain: "In vacuum (no charge) ∇·E = 0, not ρ/ε₀ which applies in presence of charge."
        },
        {
          id: "p18",
          text: "A conductor with free charge distribution reaches electrostatic equilibrium when internal field is:",
          options: ["Maximum", "Minimum", "Zero", "Uniform"],
          answer: 2,
          explain: "Internal electric field in conductor at electrostatic equilibrium is zero."
        },
        {
          id: "p19",
          text: "A particle of charge q and mass m moves in uniform B field with radius r. Cyclotron frequency ω is:",
          options: ["qB/m", "m/qB", "qB²/m", "qB/m²"],
          answer: 0,
          explain: "Cyclotron angular frequency ω = qB/m (independent of speed)."
        },
        {
          id: "p20",
          text: "The torque on magnetic dipole μ in field B is:",
          options: ["μB cosθ", "μB sinθ", "μ×B", "0"],
          answer: 2,
          explain: "Torque vector τ = μ × B (magnitude μB sinθ)."
        },
        {
          id: "p21",
          text: "The electric field inside a hollow conducting shell carrying net charge is:",
          options: ["Zero", "Uniform", "Radial decreasing", "Indeterminate"],
          answer: 0,
          explain: "Inside hollow conductor, field is zero by electrostatic shielding."
        },
        {
          id: "p22",
          text: "A solenoid length l and turns N carries current I. Magnetic field inside is B. If N doubles and I halves, B becomes:",
          options: ["2B", "B", "B/2", "0"],
          answer: 1,
          explain: "B ∝ nI = (N/l)I. Doubling N and halving I leaves nI unchanged."
        },
        {
          id: "p23",
          text: "The speed of sound in air at STP is ≈ 340 m/s. In helium it is:",
          options: ["340 m/s", "786 m/s", "1270 m/s", "450 m/s"],
          answer: 1,
          explain: "v ∝ √(T/M). Helium has lower molecular mass, so speed ≈ 786 m/s at same T."
        },
        {
          id: "p24",
          text: "In a damped oscillator, amplitude decays as e^{-bt/2m}. Damping force is:",
          options: ["-bv", "-bv²", "-bv^0", "bv"],
          answer: 0,
          explain: "Linear damping force F = -bv yields amplitude decay e^{-b t/(2m)}."
        },
        {
          id: "p25",
          text: "Two waves of amplitude A and 2A interfere perfectly in phase. Resultant intensity is:",
          options: ["9A²", "A²", "4A²", "3A²"],
          answer: 0,
          explain: "Resultant amplitude = A + 2A = 3A, intensity ∝ amplitude² = 9A²."
        },
        {
          id: "p26",
          text: "The first quantum number to be filled in multi-electron atom is:",
          options: ["1s", "2s", "2p", "3s"],
          answer: 0,
          explain: "Aufbau principle: fill lowest energy orbital 1s first."
        },
        {
          id: "p27",
          text: "The energy difference between two levels ΔE = hf. If λ=500 nm, frequency f is:",
          options: ["6×10¹⁴ Hz", "3×10¹⁴ Hz", "1.5×10¹⁴ Hz", "2×10¹⁴ Hz"],
          answer: 0,
          explain: "f = c/λ=3×10⁸/5×10⁻⁷=6×10¹⁴ Hz."
        },
        {
          id: "p28",
          text: "Time for half-life of an RC circuit is (ln2)RC. If R doubles and C halves, new half time is:",
          options: ["Same", "Double", "Half", "Quadruple"],
          answer: 0,
          explain: "t_{1/2}=RC ln2; doubling R and halving C leaves RC unchanged."
        },
        {
          id: "p29",
          text: "The gravitational potential energy at distance r from Earth surface is:",
          options: ["-GMm/r", "-GMm/(R+r)", "-mgh", "-GMm/R²"],
          answer: 0,
          explain: "U = -GMm/r for point mass or spherically symmetric Earth."
        },
        {
          id: "p30",
          text: "A particle moves such that velocity v → -v under time reversal. Angular momentum transforms as:",
          options: ["L → L", "L → -L", "L → 0", "Indeterminate"],
          answer: 1,
          explain: "L = r×p; p→-p so L→-L under time reversal."
        },
        {
          id: "p31",
          text: "The principle of superposition in waves requires:",
          options: ["Linear medium", "Nonlinear medium", "Dissipative medium", "Dispersive medium"],
          answer: 0,
          explain: "Superposition holds when wave equation is linear (no interaction between waves)."
        },
        {
          id: "p32",
          text: "In diffraction, resolving power of a circular aperture is:",
          options: ["1.22λ/D", "0.61λ/D", "2λ/D", "λ/2D"],
          answer: 1,
          explain: "Rayleigh criterion: angular resolution = 1.22λ/D in radians for circular aperture. But resolving power ~1/(θ). So θ=1.22λ/D; smallest resolvable ~0.61λ/D for two sources."
        },
        {
          id: "p33",
          text: "The power radiated by accelerating charge is given by Larmor’s formula P ∝:",
          options: ["a", "a²", "v²", "v"],
          answer: 1,
          explain: "P ∝ a² (square of acceleration) in Larmor formula."
        },
        {
          id: "p34",
          text: "A particle is projected horizontally off a cliff. The trajectory is:",
          options: ["Parabola", "Ellipse", "Hyperbola", "Straight line"],
          answer: 0,
          explain: "Projectile motion under uniform gravity yields parabolic trajectory."
        },
        {
          id: "p35",
          text: "In a rotating frame, fictitious force mω²r is called:",
          options: ["Coriolis force", "Centrifugal force", "Euler force", "Gravitational force"],
          answer: 1,
          explain: "mω²r directed outward in rotating frame is centrifugal force."
        },
        {
          id: "p36",
          text: "The energy density of electric field E is:",
          options: ["½ε₀E²", "ε₀E", "ε₀E²", "2ε₀E²"],
          answer: 0,
          explain: "Electric energy density u = ½ε₀E²."
        },
        {
          id: "p37",
          text: "The induced current in a loop to oppose change in flux is per law of:",
          options: ["Faraday", "Lenz", "Ampère", "Ohm"],
          answer: 1,
          explain: "Lenz’s law states induced current opposes change in flux direction."
        },
        {
          id: "p38",
          text: "The moment of inertia of thin spherical shell is:",
          options: ["2/3 mR²", "2/5 mR²", "3/5 mR²", "1/2 mR²"],
          answer: 0,
          explain: "I = 2/3 mR² for thin spherical shell about diameter."
        },
        {
          id: "p39",
          text: "Critical angle for glass (n=1.5) to air is:",
          options: ["42°", "49°", "60°", "30°"],
          answer: 1,
          explain: "sinθ_c = n₂/n₁ = 1/1.5 = 0.6667 ⇒ θ_c = 41.8°. Approx 42°."
        },
        {
          id: "p40",
          text: "The group velocity in dispersive media can exceed c but information velocity:",
          options: ["Also exceeds c", "Remains ≤ c", "Is zero", "Is infinite"],
          answer: 1,
          explain: "Signal speed or front velocity remains at or below c; group velocity can exceed c."
        },
        {
          id: "p41",
          text: "A charge moving with constant velocity does not radiate if motion is:",
          options: ["Uniform", "Accelerated", "Circular", "Oscillatory"],
          answer: 0,
          explain: "Only accelerated charges radiate; uniform motion produces no radiation."
        },
        {
          id: "p42",
          text: "A LC circuit self-resonates at frequency f. If L is halved, new frequency is:",
          options: ["f/√2", "√2 f", "2f", "f/2"],
          answer: 1,
          explain: "f ∝ 1/√L. Halving L gives f' = f √(1/(1/2)) = f√2."
        },
        {
          id: "p43",
          text: "Electric susceptibility χₑ is related to permittivity ε by:",
          options: ["ε = ε₀(1 + χₑ)", "ε = ε₀χₑ", "χₑ = ε/ε₀", "χₑ = ε₀/ε"],
          answer: 0,
          explain: "ε = ε₀(1 + χₑ)."
        },
        {
          id: "p44",
          text: "The diffraction pattern from single slit: central maximum width is:",
          options: ["2λD/a", "λD/a", "λD/2a", "4λD/a"],
          answer: 0,
          explain: "Width of central maximum = 2λD/a where a= slit width, D= screen distance."
        },
        {
          id: "p45",
          text: "A metal cube of side a is charged. Potential at center is kQ/a. If side doubles, potential at center:",
          options: ["kQ/(2a)", "2kQ/a", "kQ/a", "kQ/(4a)"],
          answer: 0,
          explain: "Potential of uniformly charged cube scales inversely with characteristic length. Doubling side halves potential."
        }
      ]
    },
    {
      name: "Chemistry",
      questions: [
        {
          id: "c1",
          text: "Which of the following is not a colligative property?",
          options: ["Vapor pressure lowering", "Boiling point elevation", "Melting point depression", "pH change"],
          answer: 3,
          explain: "pH change is not colligative; colligative properties depend on number of solute particles."
        },
        {
          id: "c2",
          text: "The correct order of acidity in carboxylic acids is:",
          options: ["CH₃COOH > CCl₃COOH > CF₃COOH", "CF₃COOH > CCl₃COOH > CH₃COOH", "CCl₃COOH > CF₃COOH > CH₃COOH", "CH₃COOH = CCl₃COOH = CF₃COOH"],
          answer: 1,
          explain: "Electronegative substituents stabilize conjugate base: CF₃COOH > CCl₃COOH > CH₃COOH"
        },
        {
          id: "c3",
          text: "Which of the following exhibits keto–enol tautomerism most readily?",
          options: ["Acetone", "Acetoacetic ester", "Benzaldehyde", "Ethanal"],
          answer: 1,
          explain: "Acetoacetic ester (β-keto ester) has active methylene and stabilized enol form by H-bond and conjugation."
        },
        {
          id: "c4",
          text: "The coordination number of Ag⁺ in Ag(NH₃)₂⁺ complex is:",
          options: ["2", "4", "6", "8"],
          answer: 0,
          explain: "[Ag(NH₃)₂]⁺ is linear with 2 ammonia ligands, coordination number = 2."
        },
        {
          id: "c5",
          text: "Which of the following has the highest boiling point?",
          options: ["CH₃Cl", "CH₃Br", "CH₃I", "CH₃F"],
          answer: 2,
          explain: "CH₃I has highest molar mass and strongest van der Waals forces, highest boiling point."
        },
        {
          id: "c6",
          text: "The oxidation number of S in Na₂S₂O₃ (thiosulfate) average is:",
          options: ["+5", "+4", "+2", "+0"],
          answer: 1,
          explain: "Thiosulfate has net oxidation state: 2(+1)+2x+3(-2)=0⇒2x=4⇒x=+2 average. But actual are +5 and +3."
        },
        {
          id: "c7",
          text: "Which element shows +2 oxidation state but no +1?",
          options: ["Be", "Mg", "Ca", "Zn"],
          answer: 3,
          explain: "Zn commonly exhibits +2 oxidation state; +1 is not stable due to filled d¹⁰ shell."
        },
        {
          id: "c8",
          text: "The unit cell of NaCl has how many formula units?",
          options: ["2", "4", "6", "8"],
          answer: 1,
          explain: "NaCl (rock salt) has 4 formula units per cubic unit cell."
        },
        {
          id: "c9",
          text: "Which of the following is most basic oxide?",
          options: ["CO₂", "SO₂", "CaO", "SiO₂"],
          answer: 2,
          explain: "CaO is strongly basic oxide among listed options."
        },
        {
          id: "c10",
          text: "In periodic table, the smallest atomic radius occurs at:",
          options: ["Top left", "Top right", "Bottom left", "Bottom right"],
          answer: 1,
          explain: "Atomic radius decreases across period and increases down group, smallest at top right."
        },
        {
          id: "c11",
          text: "Which of the following compounds will not exhibit geometrical isomerism?",
          options: ["CH₃CH=CHCH₃", "CH₂ClCH=CHCl", "CCl₂=CHCl", "CCl₃=CH₂"],
          answer: 3,
          explain: "For cis/trans, both carbons of double bond must have two different substituents—CCl₃=CH₂ fails."
        },
        {
          id: "c12",
          text: "The rate constant of first order reaction has units:",
          options: ["M⁻¹ s⁻¹", "s⁻¹", "M s⁻¹", "M² s⁻¹"],
          answer: 1,
          explain: "First order rate constant k has units time⁻¹ (s⁻¹)."
        },
        {
          id: "c13",
          text: "Which of the following is amphoteric?",
          options: ["ZnO", "Al₂O₃", "MnO₂", "All of above"],
          answer: 3,
          explain: "ZnO, Al₂O₃ and MnO₂ all react both with acids and bases (amphoteric)."
        },
        {
          id: "c14",
          text: "In aqueous solution, the species with highest hydration enthalpy is:",
          options: ["Li⁺", "Na⁺", "K⁺", "Rb⁺"],
          answer: 0,
          explain: "Smallest ion Li⁺ has highest charge density and thus highest hydration enthalpy."
        },
        {
          id: "c15",
          text: "Which of the following is not aromatic?",
          options: ["Cyclobutadiene", "Benzene", "Naphthalene", "Pyridine"],
          answer: 0,
          explain: "Cyclobutadiene has 4 π electrons, fails Hückel’s 4n+2 rule, antiaromatic."
        },
        {
          id: "c16",
          text: "Which of the following compounds shows keto–enol tautomerism only in presence of base?",
          options: ["Acetone", "Acetaldehyde", "Benzaldehyde", "Formaldehyde"],
          answer: 2,
          explain: "Benzaldehyde lacks α hydrogen and requires base‐promoted enolate formation to tautomerize."
        },
        {
          id: "c17",
          text: "The number of sigma bonds in ethylene (C₂H₄) is:",
          options: ["5", "6", "7", "8"],
          answer: 0,
          explain: "Ethylene has 5 σ bonds: 1 C–C and 4 C–H bonds."
        },
        {
          id: "c18",
          text: "Which of the following is strongest oxidizing agent?",
          options: ["KMnO₄", "K₂Cr₂O₇", "H₂O₂", "O₃"],
          answer: 3,
          explain: "Ozone has highest oxidation potential among listed oxidizers."
        },
        {
          id: "c19",
          text: "The correct order of thermal stability of nitrates is:",
          options: ["LiNO₃ > NaNO₃ > KNO₃", "KNO₃ > NaNO₃ > LiNO₃", "NaNO₃ > LiNO₃ > KNO₃", "KNO₃ > LiNO₃ > NaNO₃"],
          answer: 1,
          explain: "Thermal stability increases down group for nitrates of alkali metals."
        },
        {
          id: "c20",
          text: "Which oxide dissolves in water to give acidic solution?",
          options: ["CO₂", "SO₂", "NO₂", "All of above"],
          answer: 3,
          explain: "All listed non-metal oxides form acids (carbonic, sulfurous, nitric)."
        },
        {
          id: "c21",
          text: "The ionization energy of which species is highest?",
          options: ["Na", "Mg", "Al", "Si"],
          answer: 3,
          explain: "Ionization energy generally increases across period; Si highest among given."
        },
        {
          id: "c22",
          text: "Which of the following is a reducing sugar?",
          options: ["Sucrose", "Starch", "Glucose", "Cellulose"],
          answer: 2,
          explain: "Glucose has free aldehyde group that can reduce Cu²⁺ in Benedict’s test."
        },
        {
          id: "c23",
          text: "The pKa of acetic acid is 4.76. The pH of 0.1 M CH₃COOH is approximately:",
          options: ["2.4", "3.4", "4.4", "5.4"],
          answer: 1,
          explain: "pH ≈ ½(pKa - logC) = ½(4.76 + 1) ≈ 2.88. Closest answer 3.4."
        },
        {
          id: "c24",
          text: "Which of the following is linking sugar to protein in glycoproteins?",
          options: ["O-glycosidic bond", "N-glycosidic bond", "Peptide bond", "Phosphodiester bond"],
          answer: 1,
          explain: "N-glycosidic linkages attach sugars to asparagine residues in glycoproteins."
        },
        {
          id: "c25",
          text: "Which of the following coordination compounds exhibits optical isomerism?",
          options: ["[Co(en)₃]³⁺", "[PtCl₂(NH₃)₂]", "[Co(NH₃)₆]³⁺", "[Ni(CO)₄]"],
          answer: 0,
          explain: "Chiral [Co(en)₃]³⁺ has non-superimposable mirror images due to three bidentate ligands."
        },
        {
          id: "c26",
          text: "The rate of reaction doubles for every 10°C rise. Activation energy is approximately:",
          options: ["53 kJ/mol", "100 kJ/mol", "10 kJ/mol", "5 kJ/mol"],
          answer: 0,
          explain: "Rule of thumb: rate doubles per 10°C ⇒ Ea ≈ 53 kJ/mol."
        },
        {
          id: "c27",
          text: "Which of the following is NOT a buffer?",
          options: ["NH₄Cl + NH₄OH", "CH₃COOH + CH₃COONa", "NaH₂PO₄ + Na₂HPO₄", "HCl + KCl"],
          answer: 3,
          explain: "Strong acid HCl with its salt provides no buffering action."
        },
        {
          id: "c28",
          text: "The hybridization in SF₆ is:",
          options: ["sp³", "sp³d", "sp³d²", "sp³d³"],
          answer: 2,
          explain: "SF₆ has six bonding pairs around S, requiring sp³d² hybridization."
        },
        {
          id: "c29",
          text: "Which halogen has highest electron affinity?",
          options: ["F", "Cl", "Br", "I"],
          answer: 1,
          explain: "Chlorine has highest electron affinity due to compact size and electron-electron repulsion in F."
        },
        {
          id: "c30",
          text: "In an octahedral complex, t₂g orbitals are lower in energy because:",
          options: ["Direct overlap with ligands", "No overlap", "Partial overlap", "Symmetric orientation"],
          answer: 1,
          explain: "t₂g orbitals lie between axes so less repulsion, hence lower energy in octahedral field."
        },
        {
          id: "c31",
          text: "Which of the following is strongest oxidizing agent in aqueous solution?",
          options: ["Cl₂", "Br₂", "I₂", "F₂"],
          answer: 3,
          explain: "F₂ has highest standard reduction potential, strongest oxidizer."
        },
        {
          id: "c32",
          text: "Which of the following is paramagnetic?",
          options: ["V²⁺", "V³⁺", "V⁵⁺", "V⁴⁺"],
          answer: 3,
          explain: "V⁴⁺ is d¹, one unpaired electron makes it paramagnetic."
        },
        {
          id: "c33",
          text: "The bond dissociation energy of H–H is 436 kJ/mol. Energy of photon to break one H–H molecule is:",
          options: ["7.25×10⁻¹⁹ J", "4.36×10⁻¹⁹ J", "1.24×10⁻¹⁸ J", "6.37×10⁻¹⁹ J"],
          answer: 0,
          explain: "436 kJ/mol per molecule: (436×10³ J)/(6.02×10²³)≈7.25×10⁻¹⁹ J."
        },
        {
          id: "c34",
          text: "Which of the following alcohols is most easily dehydrated?",
          options: ["1° alcohol", "2° alcohol", "3° alcohol", "All equal"],
          answer: 2,
          explain: "Tertiary alcohols dehydrate most easily due to stable carbocation formation."
        },
        {
          id: "c35",
          text: "The number of chiral carbons in glucose is:",
          options: ["3", "4", "5", "2"],
          answer: 1,
          explain: "Glucose has 4 chiral centers at C2, C3, C4, and C5."
        },
        {
          id: "c36",
          text: "Which of the following is NOT a reducing sugar?",
          options: ["Glucose", "Fructose", "Sucrose", "Galactose"],
          answer: 2,
          explain: "Sucrose has no free hemiacetal group, so not reducing; others are reducing sugars."
        },
        {
          id: "c37",
          text: "Which amino acid is optically inactive due to symmetry?",
          options: ["Alanine", "Glycine", "Proline", "Valine"],
          answer: 1,
          explain: "Glycine has two identical hydrogen substituents, no chiral center, optically inactive."
        },
        {
          id: "c38",
          text: "Which nucleotide pairs via two hydrogen bonds?",
          options: ["A–T", "G–C", "A–U", "C–U"],
          answer: 0,
          explain: "A–T (or A–U in RNA) pairs by two hydrogen bonds. G–C pairs by three."
        },
        {
          id: "c39",
          text: "The enzyme that synthesizes DNA from RNA template in retroviruses is:",
          options: ["DNA polymerase", "RNA polymerase", "Reverse transcriptase", "Ligase"],
          answer: 2,
          explain: "Reverse transcriptase synthesizes complementary DNA (cDNA) using RNA template."
        },
        {
          id: "c40",
          text: "Which amino acid has imidazole side chain?",
          options: ["Lysine", "Histidine", "Arginine", "Proline"],
          answer: 1,
          explain: "Histidine contains imidazole ring in side chain."
        },
        {
          id: "c41",
          text: "Which of the following base pairs is most susceptible to UV damage?",
          options: ["A–T", "G–C", "C–G", "T–A"],
          answer: 3,
          explain: "Thymine dimers form under UV, making T–A pairs most UV-sensitive."
        },
        {
          id: "c42",
          text: "Which vitamin is cofactor for transamination reactions?",
          options: ["Vitamin B₁", "Vitamin B₆", "Vitamin B₁₂", "Vitamin B₂"],
          answer: 1,
          explain: "Pyridoxal phosphate (B₆) acts as coenzyme in transamination."
        },
        {
          id: "c43",
          text: "Which type of enzyme inhibition cannot be overcome by increasing substrate concentration?",
          options: ["Competitive", "Noncompetitive", "Uncompetitive", "Allosteric"],
          answer: 1,
          explain: "Noncompetitive inhibitors bind to allosteric site, reducing Vmax regardless of [S]."
        },
        {
          id: "c44",
          text: "The Michaelis–Menten constant Km equals:",
          options: ["Substrate concentration at Vmax", "Substrate concentration at Vmax/2", "Half of Vmax", "Cannot be determined"],
          answer: 1,
          explain: "Km is substrate concentration at which reaction velocity is half of Vmax."
        },
        {
          id: "c45",
          text: "Which of the following reactions is zero-order?",
          options: ["Enzyme saturation kinetics", "Radioactive decay", "Second order reaction", "Photochemical reaction under intense light"],
          answer: 3,
          explain: "At high light intensities, photochemical rates can become independent of substrate, appearing zero-order."
        }
      ]
    },
    {
      name: "Biology",
      questions: [
        // BOTANY – 45 questions
        {
          id: "b1",
          text: "Which plastid type stores starch in roots?",
          options: ["Chloroplast", "Chromo- plast", "Leucoplast", "Elaioplast"],
          answer: 2,
          explain: "Leucoplasts store starch (amyloplasts) in roots and tubers."
        },
        {
          id: "b2",
          text: "Which enzyme fixes CO₂ in dark reaction of CAM plants at night?",
          options: ["Rubisco", "PEP carboxylase", "Malate dehydrogenase", "NADP reductase"],
          answer: 1,
          explain: "PEP carboxylase fixes CO₂ to PEP forming oxaloacetate at night in CAM plants."
        },
        {
          id: "b3",
          text: "The guard cell walls are thicker on which side?",
          options: ["Inner", "Outer", "Both equal", "One side indeterminate"],
          answer: 0,
          explain: "Guard cells have thicker inner walls facing stomatal pore, causing shape change on turgor."
        },
        {
          id: "b4",
          text: "Which of the following is example of apogamy in plants?",
          options: ["Bryophyllum budding", "Ferns producing sporophyte from gametophyte", "Bulbils in onion", "Runners in strawberry"],
          answer: 1,
          explain: "Apogamy in ferns: sporophyte develops from gametophyte without fertilization."
        },
        {
          id: "b5",
          text: "Which of the following exhibits raphe in seeds?",
          options: ["Legume", "Drupe", "Achene", "Cypsela"],
          answer: 2,
          explain: "Achene (e.g., sunflower seed) develops raphe from funicle fusion to integument."
        },
        {
          id: "b6",
          text: "Which root modification stores water in desert plants?",
          options: ["Pneumatophore", "Storage root", "Tap root", "Aerial root"],
          answer: 1,
          explain: "Storage roots (e.g., sweet potato) store water and nutrients; in desert, roots may be succulent."
        },
        {
          id: "b7",
          text: "Which leaf movement occurs daily in some plants (e.g., leguminous leaves)?",
          options: ["Ephemeral", "Circadian", "Nyctinastic", "Thigmonastic"],
          answer: 2,
          explain: "Nyctinastic movements are sleep movements (night closure, day opening) in response to circadian rhythm."
        },
        {
          id: "b8",
          text: "Which tissue conducts water in gymnosperms?",
          options: ["Vessels", "Tracheids", "Sieve tubes", "Phloem fibers"],
          answer: 1,
          explain: "Gymnosperms have only tracheids for water conduction; vessels are absent."
        },
        {
          id: "b9",
          text: "Which hormone promotes seed dormancy?",
          options: ["Auxin", "Abscisic acid", "Gibberellin", "Cytokinin"],
          answer: 1,
          explain: "Abscisic acid maintains seed dormancy by inhibiting germination."
        },
        {
          id: "b10",
          text: "Vernation refers to arrangement of:",
          options: ["Leaves in bud", "Vascular bundles", "Pollen grains", "Root hairs"],
          answer: 0,
          explain: "Vernation is pattern of unfolding of leaf primordia in bud."
        },
        {
          id: "b11",
          text: "Which of the following is a non‐vascular plant?",
          options: ["Ferns", "Liverworts", "Gymnosperms", "Angiosperms"],
          answer: 1,
          explain: "Liverworts are bryophytes, non‐vascular lower plants."
        },
        {
          id: "b12",
          text: "Which part of flower develops into seed coat?",
          options: ["Ovary wall", "Integument", "Ovule", "Funicle"],
          answer: 1,
          explain: "Integuments harden to form testa (seed coat) around embryo and endosperm."
        },
        {
          id: "b13",
          text: "Which of the following is C₃ plant?",
          options: ["Maize", "Sugarcane", "Wheat", "Sorghum"],
          answer: 2,
          explain: "Wheat is C₃ plant; sugarcane and maize are C₄, sorghum is C₄."
        },
        {
          id: "b14",
          text: "Kranz anatomy refers to:",
          options: ["C₃ leaf structure", "C₄ leaf structure", "CAM leaf structure", "Hydrophyte leaf structure"],
          answer: 1,
          explain: "Kranz anatomy is leaf structure in C₄ plants with bundle sheath cells surrounding veins."
        },
        {
          id: "b15",
          text: "Which of the following cells are living at maturity?",
          options: ["Tracheids", "Vessel elements", "Xylem fibers", "Xylem parenchyma"],
          answer: 3,
          explain: "Xylem parenchyma are only living cells in xylem; others are dead at maturity."
        },
        {
          id: "b16",
          text: "Which structure protects developing ovule from mechanical damage?",
          options: ["Integument", "Funicle", "Nucellus", "Perisperm"],
          answer: 0,
          explain: "Integument develops into protective seed coat around embryo and endosperm."
        },
        {
          id: "b17",
          text: "Which type of placentation is common in Solanaceae (e.g., tomato)?",
          options: ["Marginal", "Axile", "Parietal", "Basal"],
          answer: 1,
          explain: "Tomato has axile placentation: ovules attached to central axis in multilocular ovary."
        },
        {
          id: "b18",
          text: "Which of the following is xerophytic leaf adaptation?",
          options: ["Thick cuticle", "Sunken stomata", "Reduced leaf area", "All of above"],
          answer: 3,
          explain: "Xerophytes show multiple adaptations: thick cuticle, sunken stomata, reduced leaf area."
        },
        {
          id: "b19",
          text: "Which enzyme catalyzes the conversion of nitrate to nitrite in plants?",
          options: ["Nitrate reductase", "Nitrite reductase", "Nitrase", "Nitrogenase"],
          answer: 0,
          explain: "Nitrate reductase converts nitrate (NO₃⁻) to nitrite (NO₂⁻) in nitrogen assimilation."
        },
        {
          id: "b20",
          text: "Which of the following shows mycorrhizal association?",
          options: ["Pea", "Orchid", "Wheat", "Rose"],
          answer: 1,
          explain: "Orchids have specialized mycorrhizal associations for seed germination and nutrient uptake."
        },
        {
          id: "b21",
          text: "Which of the following develops into perisperm in seeds?",
          options: ["Nucellus", "Integument", "Endosperm", "Funicle"],
          answer: 0,
          explain: "Perisperm is nutritive tissue derived from nucellus in some seeds (e.g., coffee)."
        },
        {
          id: "b22",
          text: "Which photoreceptor mediates seed germination in response to red light?",
          options: ["Cryptochrome", "Phototropin", "Phytochrome", "UVR8"],
          answer: 2,
          explain: "Phytochrome senses red/far-red light and triggers germination under red light."
        },
        {
          id: "b23",
          text: "Which cell organelle produces ATP in plant cells?",
          options: ["Chloroplast", "Mitochondrion", "Leucoplast", "Golgi body"],
          answer: 1,
          explain: "Mitochondria perform oxidative phosphorylation to produce ATP in both plant and animal cells."
        },
        {
          id: "b24",
          text: "Which of the following is not part of stoma complex?",
          options: ["Guard cells", "Subsidiary cells", "Hypodermis", "None of above"],
          answer: 2,
          explain: "Hypodermis is a tissue layer under epidermis, not part of stomatal complex."
        },
        {
          id: "b25",
          text: "The movement of water from soil into root xylem is:",
          options: ["Symplastic", "Apoplastic", "Transmembrane", "All of above"],
          answer: 3,
          explain: "Water moves via apoplastic, symplastic, and transmembrane pathways to reach xylem."
        },
        {
          id: "b26",
          text: "Which enzyme splits water into electrons and protons in photosynthesis?",
          options: ["ATP synthase", "Rubisco", "Water‐splitting complex", "NADP reductase"],
          answer: 2,
          explain: "PSII’s oxygen‐evolving complex splits water to replenish electrons to reaction center."
        },
        {
          id: "b27",
          text: "Which of the following is a dioecious plant?",
          options: ["Papaya", "Lemon", "Tomato", "Rice"],
          answer: 0,
          explain: "Papaya has male and female flowers on separate plants (dioecious)."
        },
        {
          id: "b28",
          text: "Which stage of flower development is called anthesis?",
          options: ["Bud formation", "Petal formation", "Flower opening", "Fruit set"],
          answer: 2,
          explain: "Anthesis is the period when flower is fully open and functional for pollination."
        },
        {
          id: "b29",
          text: "Which tissue is responsible for secondary growth in dicot stems?",
          options: ["Apical meristem", "Intercalary meristem", "Vascular cambium", "Cork cambium"],
          answer: 2,
          explain: "Vascular cambium produces secondary xylem and phloem, enabling thickening of stems."
        },
        {
          id: "b30",
          text: "Which of the following is an example of C₃ plant?",
          options: ["Maize", "Sugarcane", "Rice", "Sorghum"],
          answer: 2,
          explain: "Rice fixes CO₂ directly via Calvin cycle (C₃) in mesophyll cells without specialized anatomy."
        },
        {
          id: "b31",
          text: "Which leaf venation pattern is characteristic of dicots?",
          options: ["Parallel", "Reticulate", "Dichotomous", "Pinnate"],
          answer: 1,
          explain: "Dicots have reticulate (net‐like) venation; monocots have parallel venation."
        },
        {
          id: "b32",
          text: "Which of the following structures is NOT diploid in angiosperms?",
          options: ["Zygote", "Embryo sac", "Endosperm", "Embryo"],
          answer: 2,
          explain: "Endosperm is typically triploid following double fertilization."
        },
        {
          id: "b33",
          text: "The phenomenon of reversible leaf folding in response to touch is called:",
          options: ["Thigmotropism", "Thigmonasty", "Phototropism", "Hydrotropism"],
          answer: 1,
          explain: "Thigmonastic movement is rapid reversible response to touch (e.g., Mimosa pudica)."
        },
        {
          id: "b34",
          text: "Which plant hormone is involved in fruit ripening?",
          options: ["Auxin", "Ethylene", "Gibberellin", "Abscisic acid"],
          answer: 1,
          explain: "Ethylene gas triggers ripening in climacteric fruits."
        },
        {
          id: "b35",
          text: "Which type of vascular bundles is found in dicot stem?",
          options: ["Closed", "Open", "Scattered", "Collateral"],
          answer: 1,
          explain: "Dicot stems have open vascular bundles with cambium allowing secondary growth."
        },
        {
          id: "b36",
          text: "Which of the following is involved in gravitropism?",
          options: ["Auxin", "Cytokinin", "Gibberellin", "Ethylene"],
          answer: 0,
          explain: "Auxin redistribution causes differential growth and bending in response to gravity."
        },
        {
          id: "b37",
          text: "Which of the following is a C₄ plant?",
          options: ["Wheat", "Rice", "Maize", "Barley"],
          answer: 2,
          explain: "Maize uses C₄ pathway with Kranz anatomy to fix CO₂."
        },
        {
          id: "b38",
          text: "Which cells store starch in leaves?",
          options: ["Palisade parenchyma", "Spongy parenchyma", "Bundle sheath cells", "Guard cells"],
          answer: 2,
          explain: "Bundle sheath cells in C₄ plants store starch during Calvin cycle."
        },
        {
          id: "b39",
          text: "Which of the following is a connective tissue?",
          options: ["Epidermis", "Phloem", "Xylem", "Parenchyma"],
          answer: 3,
          explain: "Parenchyma is simple plant tissue, not connective. In animals, connective tissue binds other tissues together."
        },
        {
          id: "b40",
          text: "Which of the following produces secondary metabolites as defense?",
          options: ["Alkaloids", "Proteins", "Carbohydrates", "Lipids"],
          answer: 0,
          explain: "Alkaloids are nitrogenous secondary metabolites with defensive roles."
        },
        {
          id: "b41",
          text: "Which of the following exhibits heterospory?",
          options: ["Selaginella", "Pteris", "Cycas", "Equisetum"],
          answer: 0,
          explain: "Selaginella is a heterosporous pteridophyte producing microspores and megaspores."
        },
        {
          id: "b42",
          text: "Which of the following is the water impermeable tissue in roots?",
          options: ["Epidermis", "Exodermis", "Endodermis", "Casparian strip"],
          answer: 3,
          explain: "Casparian strip in endodermis blocks apoplastic flow, forcing symplastic uptake."
        },
        {
          id: "b43",
          text: "Which pigment is primary in C₄ bundle sheath chloroplasts?",
          options: ["Chlorophyll a", "Chlorophyll b", "Xanthophyll", "Anthocyanin"],
          answer: 0,
          explain: "Bundle sheath chloroplasts contain predominantly chlorophyll a for Calvin cycle."
        },
        {
          id: "b44",
          text: "Which part of embryo develops into shoot in dicots?",
          options: ["Radicle", "Hypocotyl", "Epicotyl", "Cotyledon"],
          answer: 2,
          explain: "Epicotyl is region above cotyledons that develops into shoot."
        },
        {
          id: "b45",
          text: "Which plant group has unprotected seeds?",
          options: ["Angiosperms", "Bryophytes", "Gymnosperms", "Pteridophytes"],
          answer: 2,
          explain: "Gymnosperms have naked seeds not enclosed in ovary."
        },

        // ZOOLOGY – 45 questions
        {
          id: "z1",
          text: "Which of the following cells are phagocytic in blood?",
          options: ["Neutrophils", "Lymphocytes", "Erythrocytes", "Platelets"],
          answer: 0,
          explain: "Neutrophils are primary phagocytes in innate immune response."
        },
        {
          id: "z2",
          text: "The hepatic portal vein carries blood from:",
          options: ["Liver to heart", "Gut to liver", "Heart to liver", "Kidney to liver"],
          answer: 1,
          explain: "Hepatic portal vein carries nutrient-rich blood from gastrointestinal tract to liver."
        },
        {
          id: "z3",
          text: "Which vitamin deficiency causes scurvy?",
          options: ["Vit A", "Vit B₁", "Vit C", "Vit D"],
          answer: 2,
          explain: "Vitamin C deficiency impairs collagen synthesis causing scurvy."
        },
        {
          id: "z4",
          text: "The conduction tissue between atria and ventricles is:",
          options: ["SA node", "AV node", "Bundle of His", "Purkinje fibers"],
          answer: 1,
          explain: "AV node delays conduction before passing impulse to ventricles via Bundle of His."
        },
        {
          id: "z5",
          text: "Which hormone increases blood calcium by stimulating bone resorption?",
          options: ["Calcitonin", "Parathyroid hormone", "Insulin", "Glucagon"],
          answer: 1,
          explain: "Parathyroid hormone raises blood Ca²⁺ by promoting osteoclast activity."
        },
        {
          id: "z6",
          text: "The functional unit of kidney is:",
          options: ["Glomerulus", "Nephron", "Bowman’s capsule", "Collecting duct"],
          answer: 1,
          explain: "Nephron is the structural and functional unit of kidney."
        },
        {
          id: "z7",
          text: "Which part of ear detects linear acceleration?",
          options: ["Cochlea", "Utricle", "Semicircular canals", "Eustachian tube"],
          answer: 1,
          explain: "Utricle (in vestibule) detects linear acceleration horizontally."
        },
        {
          id: "z8",
          text: "Which blood cell lacks nucleus?",
          options: ["Neutrophil", "Eosinophil", "Monocyte", "Erythrocyte"],
          answer: 3,
          explain: "Mature mammalian erythrocytes lack nucleus to maximize hemoglobin content."
        },
        {
          id: "z9",
          text: "Which organ produces bile?",
          options: ["Pancreas", "Liver", "Gallbladder", "Stomach"],
          answer: 1,
          explain: "Liver synthesizes bile, stored and concentrated in gallbladder."
        },
        {
          id: "z10",
          text: "Which enzyme breaks down proteins in stomach?",
          options: ["Amylase", "Pepsin", "Lipase", "Nuclease"],
          answer: 1,
          explain: "Pepsin is secreted as pepsinogen and activates in acidic stomach to digest proteins."
        },
        {
          id: "z11",
          text: "The first heart sound (lub) is due to closure of:",
          options: ["AV valves", "Semilunar valves", "Both", "None"],
          answer: 0,
          explain: "Lub is caused by closure of atrioventricular valves (mitral and tricuspid)."
        },
        {
          id: "z12",
          text: "Which hormone regulates circadian rhythm?",
          options: ["Melatonin", "Serotonin", "Dopamine", "Oxytocin"],
          answer: 0,
          explain: "Melatonin from pineal gland regulates sleep–wake cycles."
        },
        {
          id: "z13",
          text: "Which part of neuron conducts impulses away from cell body?",
          options: ["Dendrite", "Axon", "Axon hillock", "Synapse"],
          answer: 1,
          explain: "Axon transmits nerve impulses away from neuron’s cell body."
        },
        {
          id: "z14",
          text: "Which vitamin is essential for clotting?",
          options: ["Vit A", "Vit B₁₂", "Vit C", "Vit K"],
          answer: 3,
          explain: "Vitamin K is required for synthesis of clotting factors II, VII, IX, X."
        },
        {
          id: "z15",
          text: "Which of the following is NOT part of innate immunity?",
          options: ["Skin barrier", "Complement", "Antibodies", "Macrophages"],
          answer: 2,
          explain: "Antibodies are part of adaptive immunity; innate includes barriers, complement, phagocytes."
        },
        {
          id: "z16",
          text: "The largest organ in human body is:",
          options: ["Liver", "Skin", "Lungs", "Brain"],
          answer: 1,
          explain: "Skin is the largest organ by surface area and weight."
        },
        {
          id: "z17",
          text: "Which hormone regulates blood sugar by increasing it?",
          options: ["Insulin", "Glucagon", "Somatostatin", "ADH"],
          answer: 1,
          explain: "Glucagon from α cells of pancreas raises blood glucose by promoting glycogenolysis."
        },
        {
          id: "z18",
          text: "Which valve prevents backflow into left atrium?",
          options: ["Tricuspid", "Mitral", "Aortic", "Pulmonary"],
          answer: 1,
          explain: "Mitral (bicuspid) valve prevents blood returning to left atrium during ventricular systole."
        },
        {
          id: "z19",
          text: "The pacemaker region of heart is:",
          options: ["AV node", "SA node", "Bundle of His", "Purkinje fibers"],
          answer: 1,
          explain: "SA node sets heart rate by initiating electrical impulses."
        },
        {
          id: "z20",
          text: "Which blood group is universal donor?",
          options: ["A", "B", "AB", "O"],
          answer: 3,
          explain: "Type O negative lacks A and B antigens, can donate to all."
        },
        {
          id: "z21",
          text: "Which part of brain controls voluntary movement?",
          options: ["Cerebellum", "Cerebrum", "Medulla", "Pons"],
          answer: 1,
          explain: "Cerebrum (motor cortex) controls voluntary muscular activity."
        },
        {
          id: "z22",
          text: "Which part of nephron is impermeable to water?",
          options: ["Proximal tubule", "Descending loop", "Ascending loop", "Collecting duct"],
          answer: 2,
          explain: "Ascending limb is impermeable to water but actively transports ions."
        },
                        {
                    id: "z23",
                    text: "The condition of high blood sugar is called:",
                    options: ["Hypoglycemia", "Hyperglycemia", "Glycosuria", "Ketonuria"],
                    answer: 1,
                    explain: "Hyperglycemia refers to abnormally high blood glucose levels."
                },
                {
                    id: "z24",
                    text: "Which of the following carries oxygen in blood?",
                    options: ["Platelets", "Plasma", "Hemoglobin", "Leukocytes"],
                    answer: 2,
                    explain: "Hemoglobin in red blood cells binds and transports oxygen."
                },
                {
                    id: "z25",
                    text: "Which part of ear is responsible for dynamic equilibrium?",
                    options: ["Cochlea", "Vestibule", "Semicircular canals", "Eustachian tube"],
                    answer: 2,
                    explain: "Semicircular canals detect rotational movements and maintain dynamic balance."
                },
                {
                    id: "z26",
                    text: "Which gland produces insulin?",
                    options: ["Thyroid", "Adrenal", "Pancreas", "Pituitary"],
                    answer: 2,
                    explain: "Insulin is secreted by beta cells in the islets of Langerhans in the pancreas."
                },
                {
                    id: "z27",
                    text: "Which vitamin deficiency causes beriberi?",
                    options: ["B₁", "B₂", "B₃", "B₆"],
                    answer: 0,
                    explain: "Thiamine (vitamin B₁) deficiency leads to beriberi affecting cardiovascular and nervous systems."
                },
                {
                    id: "z28",
                    text: "The hormone that stimulates cell division and growth is:",
                    options: ["Thyroxine", "Insulin", "Growth hormone", "Glucagon"],
                    answer: 2,
                    explain: "Growth hormone from anterior pituitary stimulates growth and cell reproduction."
                },
                {
                    id: "z29",
                    text: "Which of the following is not a leukocyte?",
                    options: ["Neutrophil", "Platelet", "Lymphocyte", "Monocyte"],
                    answer: 1,
                    explain: "Platelets are cell fragments involved in clotting, not white blood cells."
                },
                {
                    id: "z30",
                    text: "Which organ produces bile acids?",
                    options: ["Pancreas", "Liver", "Gallbladder", "Stomach"],
                    answer: 1,
                    explain: "Bile acids are synthesized in the liver and stored in the gallbladder."
                },
                {
                    id: "z31",
                    text: "Which part of nephron fine-tunes salt balance?",
                    options: ["Proximal tubule", "Loop of Henle", "Distal tubule", "Bowman's capsule"],
                    answer: 2,
                    explain: "Distal tubule adjusts ions under hormonal control for salt balance."
                },
                {
                    id: "z32",
                    text: "Which of the following cells are involved in antibody production?",
                    options: ["T cells", "B cells", "Macrophages", "Neutrophils"],
                    answer: 1,
                    explain: "B lymphocytes differentiate into plasma cells that secrete antibodies."
                },
                {
                    id: "z33",
                    text: "The functional residual capacity is the volume left in lungs after:",
                    options: ["Maximal inhale", "Normal exhale", "Maximal exhale", "Tidal breath"],
                    answer: 1,
                    explain: "Functional residual capacity = residual volume + expiratory reserve after normal exhalation."
                },
                {
                    id: "z34",
                    text: "Which hormone is released in response to low blood calcium?",
                    options: ["Calcitonin", "Parathyroid hormone", "Calcitriol", "Insulin"],
                    answer: 1,
                    explain: "Parathyroid hormone increases blood calcium by mobilizing bone stores."
                },
                {
                    id: "z35",
                    text: "Which blood vessel has valves to prevent backflow?",
                    options: ["Arteries", "Veins", "Capillaries", "Arterioles"],
                    answer: 1,
                    explain: "Veins contain valves to maintain unidirectional blood flow back to the heart."
                },
                {
                    id: "z36",
                    text: "Which cranial nerve mediates vision?",
                    options: ["Optic", "Oculomotor", "Trigeminal", "Facial"],
                    answer: 0,
                    explain: "Optic nerve (II) transmits visual information from retina to brain."
                },
                {
                    id: "z37",
                    text: "Which of the following is not part of the respiratory zone?",
                    options: ["Respiratory bronchioles", "Alveolar ducts", "Bronchioles", "Alveoli"],
                    answer: 2,
                    explain: "Bronchioles are part of conducting zone, not gas exchange region."
                },
                {
                    id: "z38",
                    text: "The process of blood cell formation is called:",
                    options: ["Hematopoiesis", "Erythropoiesis", "Leukopoiesis", "Thrombopoiesis"],
                    answer: 0,
                    explain: "Hematopoiesis is formation of all blood cells in bone marrow."
                },
                {
                    id: "z39",
                    text: "Which of the following is not a component of innate immunity?",
                    options: ["Skin barrier", "Lysozyme", "Complement", "Antibodies"],
                    answer: 3,
                    explain: "Antibodies are adaptive immunity components; others are innate."
                },
                {
                    id: "z40",
                    text: "The most abundant protein in blood plasma is:",
                    options: ["Albumin", "Globulin", "Fibrinogen", "Hemoglobin"],
                    answer: 0,
                    explain: "Albumin maintains oncotic pressure and is most abundant plasma protein."
                },
                {
                    id: "z41",
                    text: "Which part of the conduction system delays impulse briefly before ventricles?",
                    options: ["SA node", "AV node", "Bundle branches", "Purkinje fibers"],
                    answer: 1,
                    explain: "AV node delays conduction allowing ventricles to fill before contraction."
                },
                {
                    id: "z42",
                    text: "Which of the following secretes glucagon?",
                    options: ["Alpha cells", "Beta cells", "Delta cells", "PP cells"],
                    answer: 0,
                    explain: "Alpha cells of pancreatic islets secrete glucagon to raise blood glucose."
                },
                {
                    id: "z43",
                    text: "Which organ stores bile but does not produce it?",
                    options: ["Pancreas", "Liver", "Gallbladder", "Small intestine"],
                    answer: 2,
                    explain: "Gallbladder stores and concentrates bile produced by liver."
                },
                {
                    id: "z44",
                    text: "Which part of brain regulates thirst?",
                    options: ["Hypothalamus", "Cerebrum", "Medulla", "Pons"],
                    answer: 0,
                    explain: "Hypothalamic osmoreceptors sense plasma osmolarity to regulate thirst."
                },
                {
                    id: "z45",
                    text: "The first line of defense in immunity is:",
                    options: ["Phagocytes", "Complement", "Skin and mucous membranes", "Fever"],
                    answer: 2,
                    explain: "Physical barriers like skin and mucosa provide initial protection against pathogens."
                }
      ]
    }
  ]
};
