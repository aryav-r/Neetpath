// neetpath-data.js
// 20 tests, each 180 questions (45 Physics, 45 Chemistry, 90 Biology), with comprehensive explanations
// Exposes window.NEETPATH_TESTS for neetpath.html to consume.

(function(){
    function RNG(seed){ let s = seed>>>0; return function(){ s^=s<<13; s^=s>>>17; s^=s<<5; return (s>>>0)/4294967296 } }
    const pick = (r,a)=>a[Math.floor(r()*a.length)];
    const shuffle=(r,a)=>{const x=a.slice();for(let i=x.length-1;i>0;i--){const j=Math.floor(r()*(i+1));[x[i],x[j]]=[x[j],x[i]]}return x};
    const fmt = n => { const s = Math.abs(n)<1e-3 ? n.toExponential(2) : (Math.round(n*100)/100).toString(); return s.replace(/\.0+$/,'') }

    // PHYSICS QUESTIONS
    const physicsQuestions = {
        easy: [
            {
                id: "p1",
                text: "A car accelerates from rest at 2 m/s² for 5 seconds. What is its final velocity?",
                options: ["8 m/s", "10 m/s", "12 m/s", "15 m/s"],
                answer: 1,
                explain: "Using v = u + at, where u = 0, a = 2 m/s², t = 5s. v = 0 + 2×5 = 10 m/s"
            },
            {
                id: "p2",
                text: "What is the SI unit of electric current?",
                options: ["Volt", "Ampere", "Ohm", "Watt"],
                answer: 1,
                explain: "The SI unit of electric current is Ampere (A), named after André-Marie Ampère."
            },
            {
                id: "p3",
                text: "A force of 10 N acts on a mass of 2 kg. What is the acceleration?",
                options: ["5 m/s²", "8 m/s²", "12 m/s²", "20 m/s²"],
                answer: 0,
                explain: "Using Newton's second law F = ma, a = F/m = 10/2 = 5 m/s²"
            },
            {
                id: "p4",
                text: "The frequency of a wave is 50 Hz. What is its time period?",
                options: ["0.01 s", "0.02 s", "0.05 s", "0.1 s"],
                answer: 1,
                explain: "Time period T = 1/f = 1/50 = 0.02 seconds"
            },
            {
                id: "p5",
                text: "Two resistors 4Ω and 6Ω are connected in series. Total resistance is:",
                options: ["2.4 Ω", "5 Ω", "10 Ω", "24 Ω"],
                answer: 2,
                explain: "In series connection, total resistance = R₁ + R₂ = 4 + 6 = 10 Ω"
            },
            {
                id: "p6",
                text: "The acceleration due to gravity on Earth is approximately:",
                options: ["9.8 m/s²", "8.9 m/s²", "10.8 m/s²", "11.2 m/s²"],
                answer: 0,
                explain: "The standard acceleration due to gravity is g = 9.8 m/s² (or 9.81 m/s²)"
            },
            {
                id: "p7",
                text: "A body travels 100 m in 10 s. Its average speed is:",
                options: ["5 m/s", "10 m/s", "15 m/s", "20 m/s"],
                answer: 1,
                explain: "Average speed = distance/time = 100/10 = 10 m/s"
            },
            {
                id: "p8",
                text: "The unit of electric potential is:",
                options: ["Ampere", "Volt", "Ohm", "Coulomb"],
                answer: 1,
                explain: "Electric potential is measured in Volts (V), named after Alessandro Volta."
            }
        ],
        medium: [
            {
                id: "p9",
                text: "A projectile is launched at 45° with initial velocity 20 m/s. Maximum height reached is:",
                options: ["5.1 m", "10.2 m", "15.3 m", "20.4 m"],
                answer: 1,
                explain: "Maximum height H = (v₀sinθ)²/2g = (20×sin45°)²/(2×9.8) = (20×0.707)²/19.6 ≈ 10.2 m"
            },
            {
                id: "p10",
                text: "Three capacitors 2μF, 4μF and 6μF are connected in parallel. Total capacitance is:",
                options: ["0.92 μF", "12 μF", "2.4 μF", "4.8 μF"],
                answer: 1,
                explain: "In parallel, total capacitance C = C₁ + C₂ + C₃ = 2 + 4 + 6 = 12 μF"
            },
            {
                id: "p11",
                text: "A wire of resistance 10Ω is stretched to double its length. New resistance is:",
                options: ["20 Ω", "30 Ω", "40 Ω", "50 Ω"],
                answer: 2,
                explain: "When length doubles, area becomes half. R = ρL/A, so new R = ρ(2L)/(A/2) = 4ρL/A = 4×10 = 40 Ω"
            },
            {
                id: "p12",
                text: "A ball is dropped from height 20m. Time to reach ground is:",
                options: ["1.8 s", "2.0 s", "2.2 s", "2.4 s"],
                answer: 1,
                explain: "Using s = ut + ½gt², where s = 20m, u = 0, g = 9.8 m/s². 20 = 0 + ½×9.8×t². t = √(40/9.8) ≈ 2.0 s"
            }
        ],
        hard: [
            {
                id: "p13",
                text: "A uniform rod of length L pivots about its end. Time period of small oscillations is:",
                options: ["2π√(L/3g)", "2π√(2L/3g)", "2π√(3L/2g)", "2π√(L/g)"],
                answer: 1,
                explain: "For a rod pivoting about its end, moment of inertia I = ⅓mL² and torque τ = mg(L/2)sinθ ≈ mg(L/2)θ. Using T = 2π√(I/mgd), T = 2π√(2L/3g)"
            },
            {
                id: "p14",
                text: "In LC circuit with L = 2mH and C = 8μF, angular frequency of oscillation is:",
                options: ["2.5×10⁴ rad/s", "5×10⁴ rad/s", "7.5×10⁴ rad/s", "10×10⁴ rad/s"],
                answer: 0,
                explain: "Angular frequency ω = 1/√(LC) = 1/√(2×10⁻³ × 8×10⁻⁶) = 1/√(16×10⁻⁹) = 1/(4×10⁻⁵) = 2.5×10⁴ rad/s"
            },
            {
                id: "p15",
                text: "Photo-electric work function of metal is 2eV. Threshold wavelength is:",
                options: ["620 nm", "680 nm", "720 nm", "780 nm"],
                answer: 0,
                explain: "Work function φ = hc/λ₀. λ₀ = hc/φ = (6.63×10⁻³⁴ × 3×10⁸)/(2 × 1.6×10⁻¹⁹) ≈ 620 nm"
            }
        ]
    };

    // CHEMISTRY QUESTIONS
    const chemistryQuestions = {
        easy: [
            {
                id: "c1",
                text: "The atomic number of carbon is:",
                options: ["4", "6", "8", "12"],
                answer: 1,
                explain: "Carbon has 6 protons in its nucleus, so its atomic number is 6."
            },
            {
                id: "c2",
                text: "Which is the most electronegative element?",
                options: ["Oxygen", "Nitrogen", "Fluorine", "Chlorine"],
                answer: 2,
                explain: "Fluorine has the highest electronegativity (4.0 on Pauling scale) among all elements."
            },
            {
                id: "c3",
                text: "pH of pure water at 25°C is:",
                options: ["6", "7", "8", "9"],
                answer: 1,
                explain: "Pure water is neutral with pH = 7 at 25°C due to equal H⁺ and OH⁻ concentrations."
            },
            {
                id: "c4",
                text: "Number of electrons in Mg²⁺ ion:",
                options: ["10", "12", "14", "24"],
                answer: 0,
                explain: "Mg has 12 electrons. Mg²⁺ has lost 2 electrons, so it has 12 - 2 = 10 electrons."
            },
            {
                id: "c5",
                text: "Which gas is evolved when zinc reacts with HCl?",
                options: ["O₂", "CO₂", "H₂", "Cl₂"],
                answer: 2,
                explain: "Zn + 2HCl → ZnCl₂ + H₂. Hydrogen gas is evolved when active metals react with acids."
            }
        ],
        medium: [
            {
                id: "c6",
                text: "Hybridization of carbon in CH₄ is:",
                options: ["sp", "sp²", "sp³", "sp³d"],
                answer: 2,
                explain: "In CH₄, carbon forms 4 sigma bonds with hydrogen atoms, requiring sp³ hybridization (tetrahedral geometry)."
            },
            {
                id: "c7",
                text: "For reaction N₂ + 3H₂ ⇌ 2NH₃, if Kc = 4, then Kp at 500K is:",
                options: ["4×(RT)²", "4/(RT)²", "4×RT", "4"],
                answer: 1,
                explain: "Δn = 2 - (1 + 3) = -2. Kp = Kc(RT)^Δn = 4×(RT)^(-2) = 4/(RT)²"
            },
            {
                id: "c8",
                text: "Which shows maximum covalent character?",
                options: ["NaCl", "MgO", "AlF₃", "CaF₂"],
                answer: 2,
                explain: "According to Fajan's rules, smaller cation (Al³⁺) and larger anion (F⁻) leads to maximum covalent character."
            }
        ],
        hard: [
            {
                id: "c9",
                text: "Rate constant of first order reaction is 0.693 min⁻¹. Half-life is:",
                options: ["0.5 min", "1 min", "1.5 min", "2 min"],
                answer: 1,
                explain: "For first order reaction, t₁/₂ = 0.693/k = 0.693/0.693 = 1 minute"
            },
            {
                id: "c10",
                text: "Spin-only magnetic moment of Fe³⁺ (Z=26) is:",
                options: ["3.87 BM", "4.90 BM", "5.92 BM", "6.93 BM"],
                answer: 2,
                explain: "Fe³⁺: [Ar] 3d⁵. Number of unpaired electrons n = 5. μ = √[n(n+2)] = √[5×7] = √35 = 5.92 BM"
            }
        ]
    };

    // BIOLOGY QUESTIONS
    const biologyQuestions = {
        easy: [
            {
                id: "b1",
                text: "The powerhouse of the cell is:",
                options: ["Nucleus", "Mitochondria", "Ribosome", "Golgi body"],
                answer: 1,
                explain: "Mitochondria are called powerhouse of cell because they produce ATP through cellular respiration."
            },
            {
                id: "b2",
                text: "DNA is mainly found in:",
                options: ["Cytoplasm", "Nucleus", "Ribosome", "Vacuole"],
                answer: 1,
                explain: "DNA (genetic material) is primarily located in the nucleus of eukaryotic cells."
            },
            {
                id: "b3",
                text: "The functional unit of kidney is:",
                options: ["Neuron", "Nephron", "Alveoli", "Hepatocyte"],
                answer: 1,
                explain: "Nephron is the basic functional unit of kidney, responsible for filtration and urine formation."
            },
            {
                id: "b4",
                text: "Photosynthesis occurs in:",
                options: ["Mitochondria", "Nucleus", "Chloroplast", "Ribosome"],
                answer: 2,
                explain: "Photosynthesis occurs in chloroplasts, which contain chlorophyll pigments."
            },
            {
                id: "b5",
                text: "Human heart has how many chambers?",
                options: ["2", "3", "4", "5"],
                answer: 2,
                explain: "Human heart has 4 chambers: 2 atria (left and right) and 2 ventricles (left and right)."
            }
        ],
        medium: [
            {
                id: "b6",
                text: "Insulin is secreted by:",
                options: ["α-cells", "β-cells", "δ-cells", "Acinar cells"],
                answer: 1,
                explain: "Insulin is secreted by β-cells (beta cells) of islets of Langerhans in pancreas."
            },
            {
                id: "b7",
                text: "First stable product in C₄ plants is:",
                options: ["3-PGA", "Oxaloacetate", "RuBP", "Glucose"],
                answer: 1,
                explain: "In C₄ pathway, CO₂ first combines with PEP to form 4-carbon oxaloacetate as the first stable product."
            },
            {
                id: "b8",
                text: "Double fertilization is characteristic of:",
                options: ["Gymnosperms", "Angiosperms", "Pteridophytes", "Bryophytes"],
                answer: 1,
                explain: "Double fertilization is unique to angiosperms, where one sperm fuses with egg and other with polar nuclei."
            }
        ],
        hard: [
            {
                id: "b9",
                text: "In non-cyclic photophosphorylation, final electron acceptor is:",
                options: ["O₂", "NADP⁺", "ADP", "FAD"],
                answer: 1,
                explain: "In non-cyclic photophosphorylation, NADP⁺ is the final electron acceptor, forming NADPH."
            },
            {
                id: "b10",
                text: "Linkage maps are based on:",
                options: ["Recombination frequency", "Physical distance", "DNA sequence", "Chromosome structure"],
                answer: 0,
                explain: "Linkage maps show genetic distances measured by recombination frequency between genes."
            }
        ]
    };

    function generateQuestions(subject, level, count, seed) {
        const r = RNG(seed);
        const questions = [];
        const pool = subject === 'Physics' ? physicsQuestions[level] :
                     subject === 'Chemistry' ? chemistryQuestions[level] :
                     biologyQuestions[level];
        
        for (let i = 0; i < count; i++) {
            const baseQuestion = pick(r, pool);
            questions.push({
                ...baseQuestion,
                id: `${baseQuestion.id}_${i}_${level}`
            });
        }
        
        return questions;
    }

    function makeSection(subject, testLevel, seed) {
        const r = RNG(seed);
        const questions = [];
        
        // Question distribution based on test level and subject
        let distribution;
        if (subject === 'Physics' || subject === 'Chemistry') {
            // 45 questions each for Physics and Chemistry
            distribution = testLevel === 'easy' ? {easy: 25, medium: 15, hard: 5} :
                         testLevel === 'medium' ? {easy: 15, medium: 20, hard: 10} :
                         {easy: 10, medium: 15, hard: 20};
        } else {
            // 90 questions for Biology
            distribution = testLevel === 'easy' ? {easy: 50, medium: 30, hard: 10} :
                         testLevel === 'medium' ? {easy: 30, medium: 40, hard: 20} :
                         {easy: 20, medium: 30, hard: 40};
        }

        // Generate questions for each difficulty
        Object.entries(distribution).forEach(([level, count]) => {
            const levelQuestions = generateQuestions(subject, level, count, seed + level.length);
            questions.push(...levelQuestions);
        });

        return shuffle(r, questions);
    }

    function makeTest(testNumber) {
        const testLevel = testNumber <= 7 ? 'easy' : testNumber <= 14 ? 'medium' : 'hard';
        const id = `mock-test-${testNumber}`;
        
        const sections = [
            {
                name: "Physics",
                questions: makeSection("Physics", testLevel, (testNumber * 131 + 17 + 1013904223) >>> 0)
            },
            {
                name: "Chemistry", 
                questions: makeSection("Chemistry", testLevel, (testNumber * 131 + 34 + 1013904223) >>> 0)
            },
            {
                name: "Biology",
                questions: makeSection("Biology", testLevel, (testNumber * 131 + 51 + 1013904223) >>> 0)
            }
        ];

        return {
            id,
            title: `Mock Test ${testNumber}`,
            level: testLevel,
            durationMinutes: 180,
            sections
        };
    }

    // Generate exactly 20 tests
    const tests = [];
    for (let i = 1; i <= 20; i++) {
        tests.push(makeTest(i));
    }

    // Export for the main app
    window.NEETPATH_TESTS = tests;
    console.log(`Generated ${tests.length} NEET mock tests with explanations`);

})();
