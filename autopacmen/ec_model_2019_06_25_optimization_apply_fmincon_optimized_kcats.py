import ec_model_2019_06_25_data_set_up_model
import cobra

model = ec_model_2019_06_25_data_set_up_model.set_up_ec_model_with_sbml("ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES.xml", .095)

reactions_to_change = [
        # Acetate
        "FLDR2",
        "ACKr_TG_forward",
        "ACtex_TG_forward",
        "POR5_TG_reverse",
        "PTAr_TG_reverse",
        # Glycerol
        "GLYK",
        "TRDR",
        "G3PD2_TG_forward",
        "GLYCtex_TG_forward",
        "GTHOr_TG_forward",
        # Oxoglutarate
        "AKGt2rpp_TG_forward",
        "AKGtex_TG_forward",
        # L-Alanine
        "ASPT",
        "DAAD",
        "PROD2",
        "ALAR_TG_forward",
        "ALATA_L_TG_forward",
        "ALAtex_TG_forward",
        "GLUDy_TG_forward",
        "VALTA_TG_forward",
        # Pyruvate
        "PYRtex_TG_forward",
        # Fructose
        "FRUK",
        "FRUpts2pp",
        "FRUptspp",
        "FRUtex_TG_forward",
        # Guanosine
        "ALLTAMH",
        "ALLTN",
        "GMPR",
        "GSNK",
        "GSNt2pp",
        "GUAD",
        "UGLYCH",
        "XAND",
        "GSNtex_TG_forward",
        "PPM_TG_forward",
        "PRPPS_TG_reverse",
        "PUNP3_TG_forward",
        # N-acetylglucosamine
        "ACGAptspp",
        "AGDC",
        "ACGAtex_TG_forward",
        # Fumarate
        "FUMt2_2pp",
        "FUMt2_3pp",
        "FUMtex_TG_forward",
        # Ribose
        "PGCD",
        "PSP_L",
        "RBK",
        "RIBabcpp",
        "GHMT2r_TG_forward",
        "MTHFD_TG_forward",
        "RIBtex_TG_forward",
        # L-Lactate
        "L_LACt2rpp_TG_forward",
        "L_LACtex_TG_forward",
        # Gluconate
        "GLCNt2rpp_TG_forward",
        "GLCNtex_TG_forward",
        # Sorbitol
        "SBTptspp",
        "SBTPD_TG_forward",
        "SBTtex_TG_forward",
        # Glucosamine
        "GAMptspp",
        "GAMtex_TG_forward",
        # L-Malate
        "MALt2_2pp",
        "MALt2_3pp",
        "MALtex_TG_forward",
        # Succinate
        "SUCCt2_2pp",
        "SUCCt2_3pp",
        "SUCCtex_TG_forward",
        # Maltose
        "AMALT2",
        "AMALT3",
        "AMALT4",
        "MALTabcpp",
        "MALTtexi",
        "MLTP1_TG_forward",
        "MLTP2_TG_forward",
        "MLTP3_TG_forward",
        # Glucose 6-Phosphate
        "G6Pt6_2pp",
        "G6Ptex_TG_forward",
        "PItex_TG_reverse",
        # Mannitol
        "MNLptspp",
        "M1PD_TG_forward",
        "MNLtex_TG_forward",
        # Trehalose
        "TRE6PH",
        "TREHpp",
        "TREtex_TG_forward",
        # Xylose
        "XYLK",
        "XYLt2pp",
        "RPE_TG_reverse",
        "RPI_TG_reverse",
        "XYLI1_TG_forward",
        "XYLtex_TG_forward",
        # Mannose
        "MANptspp",
        "MAN6PI_TG_forward",
        "MANtex_TG_forward",
        # Galactose
        "GALt2pp",
        "GALKr_TG_forward",
        "GALtex_TG_forward",
        "UDPG4E_TG_reverse",
        "UGLT_TG_forward"
]

new_kcats_fmincon = [
    -0.00293949861858511,
    -0.000364498508355306,
    -0.00108736088770387,
    -0.00927694491848188,
    -0.000145966917645298,
    -0.000162931631405698,
    -0.0184389617062742,
    -0.00227906005491414,
    -0.00108084724295836,
    -0.00128343670616397,
    -0.00567514632053575,
    -0.000956917571413286,
    -0.00416478943107872,
    -0.00618529134697917,
    -0.0673956858219002,
    -0.00782786783442272,
    -0.00757635248788353,
    -0.000958094817375905,
    -0.00496102959479366,
    -0.00229272399439915,
    -0.00106848064454288,
    -0.00177175124737780,
    -0.0161736073907213,
    -0.000540425128991031,
    -0.000599190025988858,
    -0.000930867498313126,
    -0.00166118564427948,
    -0.633451862402422,
    -0.00184392221561183,
    -0.00238221847093628,
    -0.00227791977951661,
    -0.000535614206953250,
    -0.00118595030335169,
    -0.000571640484791942,
    -0.00476487084146292,
    -0.00612776660267211,
    -0.00422412201879175,
    -0.000519325049686670,
    -0.000718293421325220,
    -0.000459456094731076,
    -0.00291338905048220,
    -0.00148909066502485,
    -0.000566076567520122,
    -0.000544517462599971,
    -0.00267920300222965,
    -0.00161150967300756,
    -0.00286483670647769,
    -0.00509588723420483,
    -0.00123130856233523,
    -0.000579126582572730,
    -0.00231413578253416,
    -0.000646038585346168,
    -6.81942920753279e-05,
    -0.000181912708955017,
    -0.00311162248178719,
    -0.00161270454268964,
    -0.000567664480012105,
    -0.00727135249708002,
    -0.000722378126689139,
    -0.00289232986547737,
    -0.000382893503423326,
    -0.000679319235452166,
    -0.00132179444920235,
    -0.00278218948357429,
    -0.000570102942871970,
    -0.00420877488273741,
    -0.00420878225364306,
    -0.00420877878346769,
    -0.00617553009249243,
    -4.12580450138011e-05,
    -0.0103360073764172,
    -0.0103360380110049,
    -0.0103360123772782,
    -8.59037846459746e-05,
    -0.000210938598935829,
    -0.000224562919825066,
    -0.000455180928918913,
    -0.000828999913397899,
    -0.000568371588500959,
    -0.0375448658233615,
    -0.00184842308729540,
    -0.000874710033884501,
    -0.000412697200847156,
    -0.00334028747308522,
    -0.000137167092070952,
    -8.04772105655222e-05,
    -0.000840854355253846,
    -0.000622381255125240,
    -0.00996779542515868,
    -2.47580329602882e-05,
    -0.000643518851061915,
    -0.00324712270485216,
    -0.0161922306874758,
    -0.000632509288906700,
    -0.00122234825750356,
    -0.00106923064693653,
]

i = 0
prot_pool_metabolite = model.metabolites.get_by_id("prot_pool")
for reaction_to_change in reactions_to_change[:len(new_kcats_fmincon)]:
    reaction = model.reactions.get_by_id(reaction_to_change)
    old_stoichiometry = reaction.metabolites[prot_pool_metabolite]
    reaction.add_metabolites({
        prot_pool_metabolite: -old_stoichiometry
    })
    reaction.add_metabolites({
        prot_pool_metabolite: new_kcats_fmincon[i]
    })

    i += 1

# cobra.io.write_sbml_model(model, "ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES_FMINCON_CHANGE_FACTOR_50.xml")
cobra.io.write_sbml_model(model, "ec_model_2019_06_25_output_optimization/iJO1366star.xml")
