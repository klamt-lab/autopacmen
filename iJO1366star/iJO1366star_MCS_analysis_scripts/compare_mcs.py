import xlsxwriter


def _get_mcs_sets(mcs_lines):
    mcs_sets = []
    for mcs_line in mcs_lines:
        mcs_line_splitted = mcs_line.split("\t")
        mcs_line_splitted = [x for x in mcs_line_splitted if len(x) > 0]
        mcs_line_tuple = set(mcs_line_splitted)
        mcs_sets.append(mcs_line_tuple)
    return mcs_sets


def _get_subset_mapping(mcs_one, mcs_two):
    subset_mapping = {}
    for reaction_one in mcs_one:
        is_subset = False
        for reaction_two in mcs_two:
            if reaction_one.issubset(reaction_two):
                if reaction_one != reaction_two:
                    is_subset = True
                    break
        subset_mapping[tuple(reaction_one)] = is_subset
    return subset_mapping


def _get_superset_mapping(mcs_one, mcs_two):
    superset_mapping = {}
    for reaction_one in mcs_one:
        is_superset = False
        for reaction_two in mcs_two:
            if reaction_one.issuperset(reaction_two):
                if reaction_one != reaction_two:
                    is_superset = True
                    break
        superset_mapping[tuple(reaction_one)] = is_superset
    return superset_mapping


max_mcs = 6
secretion_products = ["ethanol", "leucine", "succinate", "valine"]

workbook = xlsxwriter.Workbook(
    "./autopacmen/iJOstar_MCS_analysis_scripts/mcs_analysis.xlsx"
)
for secretion_product in secretion_products:
    worksheet = workbook.add_worksheet(secretion_product)

    basepath = f"./autopacmen/iJOstar_MCS_analysis_scripts/results_{secretion_product}/"
    with open(basepath + "original_6.txt", "r") as f:
        original_mcs = f.readlines()[1:]
    original_mcs_lines = [x.replace("R_", "").replace("\n", "") for x in original_mcs]
    with open(basepath + "geckoed_6.txt", "r") as f:
        geckoed_mcs = f.readlines()[1:]
    geckoed_mcs_lines = [
        x.replace("R_", "").replace("_TG_", "TG").replace("\n", "") for x in geckoed_mcs
    ]

    original_mcs = _get_mcs_sets(original_mcs_lines)
    geckoed_mcs = _get_mcs_sets(geckoed_mcs_lines)
    mcs_only_in_original = [x for x in original_mcs if x not in geckoed_mcs]
    mcs_only_in_geckoed = [x for x in geckoed_mcs if x not in original_mcs]

    original_subset_mapping = _get_subset_mapping(original_mcs, geckoed_mcs)
    original_superset_mapping = _get_superset_mapping(original_mcs, geckoed_mcs)
    geckoed_subset_mapping = _get_subset_mapping(geckoed_mcs, original_mcs)
    geckoed_superset_mapping = _get_superset_mapping(geckoed_mcs, original_mcs)

    original_num_identicals = len(
        set([tuple(x) for x in original_mcs]).intersection(
            [tuple(x) for x in geckoed_mcs]
        )
    )
    geckoed_num_identicals = len(
        set([tuple(x) for x in geckoed_mcs]).intersection(
            [tuple(x) for x in original_mcs]
        )
    )

    current_start = 0
    # Original MCS
    worksheet.write(0, current_start, "MCS - Original")
    i = 1
    num_subsets = 0
    num_supersets = 0
    num_identicals = 0
    for mcs in original_mcs:
        j = 0
        for reaction in mcs:
            worksheet.write(i, j, reaction)
            j += 1
        worksheet.write(i, current_start + max_mcs, original_subset_mapping[tuple(mcs)])
        worksheet.write(
            i, current_start + max_mcs + 1, original_superset_mapping[tuple(mcs)]
        )

        if original_subset_mapping[tuple(mcs)]:
            num_subsets += 1
        if original_superset_mapping[tuple(mcs)]:
            num_supersets += 1

        i += 1
    worksheet.write(0, current_start + max_mcs, "Is subset?")
    worksheet.write(0, current_start + max_mcs + 1, "Is superset?")
    worksheet.write(0, current_start + max_mcs + 2, "Length")
    worksheet.write(0, current_start + max_mcs + 3, "#")
    for current_length in range(1, max_mcs + 1):
        worksheet.write(current_length, current_start + max_mcs + 2, current_length)
        worksheet.write(
            current_length,
            current_start + max_mcs + 3,
            len([x for x in original_mcs if len(x) == (current_length)]),
        )
    worksheet.write(max_mcs + 1, current_start + max_mcs + 2, "SUM")
    worksheet.write(max_mcs + 1, current_start + max_mcs + 3, len(original_mcs))
    worksheet.write(max_mcs + 2, current_start + max_mcs + 2, "Num subsets")
    worksheet.write(max_mcs + 2, current_start + max_mcs + 3, num_subsets)
    worksheet.write(max_mcs + 3, current_start + max_mcs + 2, "Num supersets")
    worksheet.write(max_mcs + 3, current_start + max_mcs + 3, num_supersets)
    worksheet.write(max_mcs + 4, current_start + max_mcs + 2, "Num identicals")
    worksheet.write(max_mcs + 4, current_start + max_mcs + 3, original_num_identicals)
    worksheet.write(max_mcs + 5, current_start + max_mcs + 2, "Neither of any")
    worksheet.write(
        max_mcs + 5,
        current_start + max_mcs + 3,
        len(original_mcs) - num_subsets - num_supersets - original_num_identicals,
    )

    current_start = 1 * max_mcs + 5
    # Geckoed MCS
    worksheet.write(0, current_start, "MCS - Geckoed")
    i = 1
    num_subsets = 0
    num_supersets = 0
    num_identicals = 0
    for mcs in geckoed_mcs:
        j = 0
        for reaction in mcs:
            worksheet.write(i, current_start + j, reaction)
            j += 1
        worksheet.write(i, current_start + max_mcs, geckoed_subset_mapping[tuple(mcs)])
        worksheet.write(
            i, current_start + max_mcs + 1, geckoed_superset_mapping[tuple(mcs)]
        )

        if geckoed_subset_mapping[tuple(mcs)]:
            num_subsets += 1
        if geckoed_superset_mapping[tuple(mcs)]:
            num_supersets += 1

        i += 1
    worksheet.write(0, current_start + max_mcs, "Is subset?")
    worksheet.write(0, current_start + max_mcs + 1, "Is superset?")
    worksheet.write(0, current_start + max_mcs + 2, "Length")
    worksheet.write(0, current_start + max_mcs + 3, "#")
    for current_length in range(1, max_mcs + 1):
        worksheet.write(current_length, current_start + max_mcs + 2, current_length)
        worksheet.write(
            current_length,
            current_start + max_mcs + 3,
            len([x for x in geckoed_mcs if len(x) == (current_length)]),
        )
    worksheet.write(max_mcs + 1, current_start + max_mcs + 2, "SUM")
    worksheet.write(max_mcs + 1, current_start + max_mcs + 3, len(geckoed_mcs))

    worksheet.write(max_mcs + 2, current_start + max_mcs + 2, "Num subset")
    worksheet.write(max_mcs + 2, current_start + max_mcs + 3, num_subsets)
    worksheet.write(max_mcs + 3, current_start + max_mcs + 2, "Num superset")
    worksheet.write(max_mcs + 3, current_start + max_mcs + 3, num_supersets)
    worksheet.write(max_mcs + 4, current_start + max_mcs + 2, "Num identicals")
    worksheet.write(max_mcs + 4, current_start + max_mcs + 3, geckoed_num_identicals)
    worksheet.write(max_mcs + 5, current_start + max_mcs + 2, "Neither of any")
    worksheet.write(
        max_mcs + 5,
        current_start + max_mcs + 3,
        len(geckoed_mcs) - num_subsets - num_supersets - geckoed_num_identicals,
    )

    current_start = 2 * max_mcs + 10
    # MCS only in original model
    worksheet.write(0, current_start, "MCS - Original only")
    i = 1
    for mcs in mcs_only_in_original:
        j = 0
        for reaction in mcs:
            worksheet.write(i, current_start + j, reaction)
            j += 1
        worksheet.write(i, current_start + max_mcs, original_subset_mapping[tuple(mcs)])
        worksheet.write(
            i, current_start + max_mcs + 1, original_superset_mapping[tuple(mcs)]
        )

        i += 1
    worksheet.write(0, current_start + max_mcs, "Is subset?")
    worksheet.write(0, current_start + max_mcs + 1, "Is superset?")
    worksheet.write(0, current_start + max_mcs + 2, "Length")
    worksheet.write(0, current_start + max_mcs + 3, "#")
    for current_length in range(1, max_mcs + 1):
        worksheet.write(current_length, current_start + max_mcs + 2, current_length)
        worksheet.write(
            current_length,
            current_start + max_mcs + 3,
            len([x for x in mcs_only_in_original if len(x) == (current_length)]),
        )
    worksheet.write(max_mcs + 1, current_start + max_mcs + 2, "SUM")
    worksheet.write(max_mcs + 1, current_start + max_mcs + 3, len(mcs_only_in_original))

    current_start = 3 * max_mcs + 15
    # MCS only in geckoed model
    worksheet.write(0, current_start, "MCS - Geckoed only")
    i = 1
    for mcs in mcs_only_in_geckoed:
        j = 0
        for reaction in mcs:
            worksheet.write(i, current_start + j, reaction)
            j += 1
        worksheet.write(i, current_start + max_mcs, geckoed_subset_mapping[tuple(mcs)])
        worksheet.write(
            i, current_start + max_mcs + 1, geckoed_superset_mapping[tuple(mcs)]
        )
        i += 1
    worksheet.write(0, current_start + max_mcs, "Is subset?")
    worksheet.write(0, current_start + max_mcs + 1, "Is superset?")
    worksheet.write(0, current_start + max_mcs + 2, "Length")
    worksheet.write(0, current_start + max_mcs + 3, "#")
    for current_length in range(1, max_mcs + 1):
        worksheet.write(current_length, current_start + max_mcs + 2, current_length)
        worksheet.write(
            current_length,
            current_start + max_mcs + 3,
            len([x for x in mcs_only_in_geckoed if len(x) == (current_length)]),
        )
    worksheet.write(max_mcs + 1, current_start + max_mcs + 2, "SUM")
    worksheet.write(max_mcs + 1, current_start + max_mcs + 3, len(mcs_only_in_geckoed))


workbook.close()
