@testitem "Form factors" begin

    # Check parsing of spectroscopic term symbols

    @test Sunny.parse_term_symbol("6S5/2") == (; s=5/2, l=0, j=5/2)
    @test Sunny.parse_term_symbol("10D0") == (; s=9/2, l=2, j=0)
    @test Sunny.parse_term_symbol("2F5/2") == (; s=1/2, l=3, j=5/2)

    # Check display

    @test sprint(show, FormFactor("Fe2")) == "FormFactor(\"Fe2\"; config=\"3d⁶\", j2_weight=0)"
    @test sprint(show, FormFactor("Fe2"; j2_weight=:free_ion)) == "FormFactor(\"Fe2\"; config=\"3d⁶\", j2_weight=1/3)"
    @test sprint(show, one(FormFactor)) == "one(FormFactor)"
    @test sprint(show, zero(FormFactor)) == "zero(FormFactor)"

    # Compute form factors for all ion types, and on arbitrary arguments

    q2 = 0.57
    function check_ffs(keys, refs, configs=fill(nothing, length(keys)); j2_weight)
        for (key, config, ref) in zip(keys, configs, refs)
            @test ref ≈ Sunny.compute_form_factor(FormFactor(key; config, j2_weight), q2)
        end
    end

    brown_keys = ["Sc0", "Sc1", "Sc2", "Ti0", "Ti1", "Ti2", "Ti3", "V0", "V2", "V3", "V4", "Cr0", "Cr1", "Cr3", "Cr4", "Mn0", "Mn1", "Mn2", "Mn4", "Mn5", "Fe0", "Fe1", "Fe2", "Fe3", "Co0", "Co1", "Co2", "Co3", "Co4", "Ni0", "Ni1", "Ni2", "Ni3", "Ni4", "Cu0", "Cu2", "Cu3", "Cu4", "Y0", "Zr0", "Zr1", "Nb0", "Mo0", "Mo1", "Tc0", "Tc1", "Ru0", "Ru1", "Rh0", "Rh1", "Pd1", "Ce2", "Nd2", "Nd3", "Sm3", "Eu2", "Gd2", "Gd3", "Tb2", "Tb3", "Dy2", "Dy3", "Ho2", "Ho3", "Er2", "Er3", "Tm2", "Tm3", "Yb3", "Pr3", "U3", "U4", "U5", "Np3", "Np4", "Np5", "Np6", "Pu3", "Pu4", "Pu5", "Pu6", "Am2", "Am4", "Am5", "Am6", "Am7", "Ce3"]
    brown_refs = [0.9404337568, 1.0102252421, 0.9691031535, 0.9845639150, 1.0400364232, 0.9848718135, 0.9793265320, 1.0306567026, 1.0268730999, 0.9909639213, 0.9854311960, 0.9450154960, 0.9421240065, 1.0218504918, 0.9938592553, 0.9558570749, 0.9487615985, 0.9603604067, 1.0172888858, 0.9749140365, 0.9658076229, 0.9596712711, 0.9687388754, 0.9701542597, 0.9708821117, 0.9683743194, 0.9732842981, 0.9760181032, 0.9765403351, 0.9802214709, 0.9724819212, 0.9769926098, 0.9797787893, 0.9808579486, 0.9657441307, 0.9796541349, 0.9819513156, 0.9834544599, 0.9237751814, 0.9665886930, 1.0693536967, 0.8755127829, 0.9188644657, 0.9017249567, 0.9284121234, 0.9120487339, 0.9469756019, 0.9374084434, 0.9530682750, 0.9514340450, 0.9529680222, 0.9859644599, 0.9965196236, 0.9911020302, 1.0353711653, 0.9716200199, 0.9707082381, 0.9772613581, 0.9796652823, 0.9813662594, 0.9817183183, 0.9836009799, 0.9831343698, 0.9851486780, 0.9846378588, 0.9863509727, 0.9855360814, 0.9871936270, 0.9878390389, 0.9889365244, 0.9796871913, 0.9776921513, 0.9772919614, 0.9936497864, 0.9848881540, 0.9822462046, 0.9802356766, 1.0568699926, 0.9956884083, 0.9873548394, 0.9847976634, 0.9462331063, 1.0577468360, 0.9971133930, 0.9889680030, 0.9858637400, 0.9840897557]
    check_ffs(brown_keys, brown_refs; j2_weight=:free_ion)

    kobayashi_keys = ["Hf2", "Hf3", "Ta2", "Ta3", "Ta4", "W0", "W0", "W1", "W1", "W3", "W4", "W5", "Re0", "Re0", "Re0", "Re1", "Re1", "Re2", "Re4", "Re5", "Re6", "Os0", "Os0", "Os0", "Os1", "Os1", "Os2", "Os3", "Os5", "Os6", "Os7", "Ir0", "Ir0", "Ir0", "Ir1", "Ir1", "Ir2", "Ir3", "Ir4", "Ir6", "Pt1", "Pt2", "Pt3", "Pt4", "Pt5", "Au1", "Au2", "Au3", "Au4", "Au5"]
    kobayashi_configs = [nothing, nothing, nothing, nothing, nothing, "6s⁰5d⁶", "6s¹5d⁵", "6s⁰5d⁵", "6s¹5d⁴", nothing, nothing, nothing, "6s⁰5d⁷", "6s¹5d⁶", "6s²5d⁵", "6s⁰5d⁶", "6s¹5d⁵", nothing, nothing, nothing, nothing, "6s⁰5d⁸", "6s¹5d⁷", "6s²5d⁶", "6s⁰5d⁷", "6s¹5d⁶", nothing, nothing, nothing, nothing, nothing, "6s⁰5d⁹", "6s¹5d⁸", "6s²5d⁷", "6s⁰5d⁸", "6s¹5d⁷", nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing]
    kobayashi_refs = [0.9697898796, 0.9545865128, 1.0603393256, 0.9780727435, 0.9637145024, 0.8700261400, 0.8685191095, 0.8811491124, 0.8740740070, 1.0494421964, 0.9825141085, 0.9698399600, 0.8961489116, 0.8979473918, 0.8961276144, 0.9080587312, 0.9024973689, 0.9085204449, 1.0420620141, 0.9852956945, 0.9739488132, 0.9131917824, 0.9150620109, 0.9190058841, 0.9235770719, 0.9217021429, 0.9275322677, 0.9249061720, 1.0365922415, 0.9874125646, 0.9771334021, 0.9252574144, 0.9263532230, 0.9322720419, 0.9341928215, 0.9331499710, 0.9385768259, 0.9397049516, 0.9360282199, 1.0323550076, 0.9418560047, 0.9461298470, 0.9482647373, 0.9481314065, 0.9442348506, 0.9453306935, 0.9518189600, 0.9542743103, 0.9550841600, 0.9545380045]
    check_ffs(kobayashi_keys, kobayashi_refs, kobayashi_configs; j2_weight=:free_ion)

    brown_singlet_keys = ["V1", "Cr2", "Mn3", "Fe4", "Cu1", "Nb1", "Pd0", "Sm2", "Eu3", "Yb2", "Am3"]
    brown_singlet_refs = [0.9296780690, 0.9537320763, 0.9656789199, 0.9730528909, 0.9664889420, 0.8811771880, 0.9458897742, 0.9686572579, 0.9753447860, 0.9797343518, 0.9519159578]
    check_ffs(brown_singlet_keys, brown_singlet_refs; j2_weight=-3/43)

    kobayashi_singlet_keys = ["W0", "W2", "Re3", "Os4", "Ir5", "Pt6"]
    kobayashi_singlet_configs = ["6s²5d⁴", "6s⁰5d⁴", nothing, nothing, nothing, nothing]
    kobayashi_singlet_refs = [0.8792679089, 0.8960938494, 0.9161468433, 0.9293869807, 0.9387216007, 0.9460141761]
    check_ffs(kobayashi_singlet_keys, kobayashi_singlet_refs, kobayashi_singlet_configs; j2_weight=-3/43)

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
