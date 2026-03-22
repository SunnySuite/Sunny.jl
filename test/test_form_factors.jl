@testitem "Form factors" begin

    # Check parsing of spectroscopic term symbols

    @test Sunny.parse_term_symbol("6S5/2") == (; s=5/2, l=0, j=5/2)
    @test Sunny.parse_term_symbol("10D0") == (; s=9/2, l=2, j=0)
    @test Sunny.parse_term_symbol("2F5/2") == (; s=1/2, l=3, j=5/2)

    # Check display

    @test repr(FormFactor("Fe2")) == "FormFactor(\"Fe2\"; config=\"3d⁶\", j2_weight=0)"
    msg = "Consider a custom j2_weight to account for orbital quenching in 3d transition metals"
    ff = @test_logs (:warn, msg) FormFactor("Fe2"; j2_weight=:free_ion)
    @test repr(ff) == "FormFactor(\"Fe2\"; config=\"3d⁶\", j2_weight=1/3)"
    @test repr(one(FormFactor)) == "one(FormFactor)"
    @test repr(zero(FormFactor)) == "zero(FormFactor)"

    # Compute form factors for all ion types, and on arbitrary arguments

    q2 = 0.57
    function check_ffs(keys, refs, configs=fill(nothing, length(keys)); j2_weight)
        for (key, config, ref) in zip(keys, configs, refs)
            @test ref ≈ Sunny.compute_form_factor(FormFactor(key; config, j2_weight), q2)
        end
    end

    brown_3d_keys = ["Sc0", "Sc1", "Sc2", "Ti0", "Ti1", "Ti2", "Ti3", "V0", "V2", "V3", "V4", "Cr0", "Cr1", "Cr3", "Cr4", "Mn0", "Mn1", "Mn2", "Mn4", "Mn5", "Fe0", "Fe1", "Fe2", "Fe3", "Co0", "Co1", "Co2", "Co3", "Co4", "Ni0", "Ni1", "Ni2", "Ni3", "Ni4", "Cu0", "Cu2", "Cu3", "Cu4", "V1", "Cr2", "Mn3", "Fe4", "Cu1"]
    brown_3d_refs = [0.8697961854, 0.8948823223, 0.9280679196, 0.9284788524, 0.9192421914, 0.9411439174, 0.9528748452, 0.9425756230, 0.9500344548, 0.9595439920, 0.9655051866, 0.9470527864, 0.9442752111, 0.9643926517, 0.9696910698, 0.9575192927, 0.9506780639, 0.9618982113, 0.9725302343, 0.9748460372, 0.9622332230, 0.9564842003, 0.9654220904, 0.9713328208, 0.9655125121, 0.9607398355, 0.9683228814, 0.9735555925, 0.9774696305, 0.9750228207, 0.9647047030, 0.9713015478, 0.9758901268, 0.9788421653, 0.9670133332, 0.9737428228, 0.9774373048, 0.9802693376, 0.9340088191, 0.9566815141, 0.9679234542, 0.9747800566, 0.9685835750]
    check_ffs(brown_3d_keys, brown_3d_refs; j2_weight=0.1)

    brown_keys = ["Y0", "Zr0", "Zr1", "Nb0", "Mo0", "Mo1", "Tc0", "Tc1", "Ru0", "Ru1", "Rh0", "Rh1", "Pd1", "Ce2", "Nd2", "Nd3", "Sm3", "Eu2", "Gd2", "Gd3", "Tb2", "Tb3", "Dy2", "Dy3", "Ho2", "Ho3", "Er2", "Er3", "Tm2", "Tm3", "Yb3", "Pr3", "U3", "U4", "U5", "Np3", "Np4", "Np5", "Np6", "Pu3", "Pu4", "Pu5", "Pu6", "Am2", "Am4", "Am5", "Am6", "Am7", "Ce3"]
    brown_refs = [0.9237751814, 0.9665886930, 1.0693536967, 0.8755127829, 0.9188644657, 0.9017249567, 0.9284121234, 0.9120487339, 0.9469756019, 0.9374084434, 0.9530682750, 0.9514340450, 0.9529680222, 0.9859644599, 0.9965196236, 0.9911020302, 1.0353711653, 0.9716200199, 0.9707082381, 0.9772613581, 0.9796652823, 0.9813662594, 0.9817183183, 0.9836009799, 0.9831343698, 0.9851486780, 0.9846378588, 0.9863509727, 0.9855360814, 0.9871936270, 0.9878390389, 0.9889365244, 0.9796871913, 0.9776921513, 0.9772919614, 0.9936497864, 0.9848881540, 0.9822462046, 0.9802356766, 1.0568699926, 0.9956884083, 0.9873548394, 0.9847976634, 0.9462331063, 1.0577468360, 0.9971133930, 0.9889680030, 0.9858637400, 0.9840897557]
    check_ffs(brown_keys, brown_refs; j2_weight=:free_ion)

    kobayashi_keys = ["Hf2", "Hf3", "Ta2", "Ta3", "Ta4", "W0", "W0", "W1", "W1", "W3", "W4", "W5", "Re0", "Re0", "Re0", "Re1", "Re1", "Re2", "Re4", "Re5", "Re6", "Os0", "Os0", "Os0", "Os1", "Os1", "Os2", "Os3", "Os5", "Os6", "Os7", "Ir0", "Ir0", "Ir0", "Ir1", "Ir1", "Ir2", "Ir3", "Ir4", "Ir6", "Pt1", "Pt2", "Pt3", "Pt4", "Pt5", "Au1", "Au2", "Au3", "Au4", "Au5"]
    kobayashi_configs = [nothing, nothing, nothing, nothing, nothing, "6s⁰5d⁶", "6s¹5d⁵", "6s⁰5d⁵", "6s¹5d⁴", nothing, nothing, nothing, "6s⁰5d⁷", "6s¹5d⁶", "6s²5d⁵", "6s⁰5d⁶", "6s¹5d⁵", nothing, nothing, nothing, nothing, "6s⁰5d⁸", "6s¹5d⁷", "6s²5d⁶", "6s⁰5d⁷", "6s¹5d⁶", nothing, nothing, nothing, nothing, nothing, "6s⁰5d⁹", "6s¹5d⁸", "6s²5d⁷", "6s⁰5d⁸", "6s¹5d⁷", nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing]
    kobayashi_refs = [0.9697898796, 0.9545865128, 1.0603393256, 0.9780727435, 0.9637145024, 0.8700261400, 0.8685191095, 0.8811491124, 0.8740740070, 1.0494421964, 0.9825141085, 0.9698399600, 0.8961489116, 0.8979473918, 0.8961276144, 0.9080587312, 0.9024973689, 0.9085204449, 1.0420620141, 0.9852956945, 0.9739488132, 0.9131917824, 0.9150620109, 0.9190058841, 0.9235770719, 0.9217021429, 0.9275322677, 0.9249061720, 1.0365922415, 0.9874125646, 0.9771334021, 0.9252574144, 0.9263532230, 0.9322720419, 0.9341928215, 0.9331499710, 0.9385768259, 0.9397049516, 0.9360282199, 1.0323550076, 0.9418560047, 0.9461298470, 0.9482647373, 0.9481314065, 0.9442348506, 0.9453306935, 0.9518189600, 0.9542743103, 0.9550841600, 0.9545380045]
    check_ffs(kobayashi_keys, kobayashi_refs, kobayashi_configs; j2_weight=:free_ion)

    brown_singlet_keys = ["Nb1", "Pd0", "Sm2", "Eu3", "Yb2", "Am3"]
    brown_singlet_refs = [0.8886670644, 0.9493697730, 0.9706780518, 0.9769639669, 0.9810535685, 0.9550229409]
    check_ffs(brown_singlet_keys, brown_singlet_refs; j2_weight=0.1)

    kobayashi_singlet_keys = ["W0", "W2", "Re3", "Os4", "Ir5", "Pt6"]
    kobayashi_singlet_configs = ["6s²5d⁴", "6s⁰5d⁴", nothing, nothing, nothing, nothing]
    kobayashi_singlet_refs = [0.8869759320, 0.9027027582, 0.9215934489, 0.9339654429, 0.9427037649, 0.9495372842]
    check_ffs(kobayashi_singlet_keys, kobayashi_singlet_refs, kobayashi_singlet_configs; j2_weight=0.1)

    # Check error messages

    @test_throws "Provide element and ionization state, e.g. \"Fe2\" for Fe²⁺" FormFactor("Fe")
    @test_throws "Select electronic `config` from \"6s⁰5d⁹\" or \"6s¹5d⁸\" or \"6s²5d⁷\"" FormFactor("Ir0")
    @test_logs (:warn, r"Suffix `a` is deprecated and will be ignored") begin
        @test_throws "Select electronic `config` from \"6s⁰5d⁹\" or \"6s¹5d⁸\" or \"6s²5d⁷\"" FormFactor("Ir0a")
    end
    @test_throws "No form factor data for \"Pr5\"; contact us if you can provide it" FormFactor("Pr5")
    @test_throws "Unrecognized magnetic element" FormFactor("H0")
    @test_throws "Free-ion Landé factor is ambiguous when J=0 (config 3d¹⁰)" Sunny.compute_form_factor(FormFactor("Cu1"; j2_weight=:free_ion), q2)

    # Check length units

    units1 = Units(:meV, :angstrom)
    units2 = Units(:meV, :nm)
    x1 = Sunny.compute_form_factor(FormFactor("Sc0"; units1.length), q2 * units1.angstrom^(-2))
    x2 = Sunny.compute_form_factor(FormFactor("Sc0"; units2.length), q2 * units2.angstrom^(-2))
    @test isapprox(x1, x2; atol=1e-14)
end
