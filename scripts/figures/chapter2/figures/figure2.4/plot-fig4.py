def generate_predictions_and_enzymes():
    """
    Generate the intermediate 'predictions_and_enzymes.tsv' table from three
    supplementary Excel tables (S1, S3, S4) and the bacteria metadata table.

    Uses global variables:
        S1_Table_xlsx           – Path to S1_Table.xlsx  (lytic depolymerases)
        S3_Table_xlsx           – Path to S3_Table.xlsx  (GWAS predictions)
        S4_Table_xlsx           – Path to S4_Table.xlsx  (lysogenic depolymerases)
        bacteria_table          – Path to bacteria.tsv
        predictions_and_enzymes – Path where the output TSV will be written
    """

    import pandas as pd
    import numpy as np

    global S1_Table_xlsx
    global S3_Table_xlsx
    global S4_Table_xlsx
    global bacteria_table
    global predictions_and_enzymes

    # ------------------------------------------------------------------
    # STEP 1 – Process recombinant depolymerases (lysogenic + lytic)
    # ------------------------------------------------------------------
    print('Processing recombinant depolymerases...', end='\t')

    # --- lysogenic (S4) ---
    rename_lysogenic_dict = {
        'protein seq': 'seq',
        'expression level': 'expression',
        'K_locus_specificity': 'specificity'
    }
    rename_lytic_dict = {
        'K_locus_specificity': 'specificity',
        'protein_seq': 'seq'
    }
    genscript_list = [
        '0367_12', '0574_17', '0391_11', '0496_72', '021_14',
        '391_03', '184_04', '145_08', '174_38', '164_08'
    ]
    recombinant_cols_step1 = [
        'proteinID', 'report', 'group', 'expression',
        'specificity', 'K_locus_host', 'seq'
    ]

    lysogenic_df = pd.read_excel(
        S4_Table_xlsx,
        sheet_name='S4A_Produced_depolymerases',
        usecols=[0, 3, 5, 6, 18],
        header=0,
    ).rename(rename_lysogenic_dict, axis=1)

    lytic_df = pd.read_excel(
        S1_Table_xlsx,
        sheet_name='2025-03-18_LITERATURE_SEARCH'
    ).rename(rename_lytic_dict, axis=1)

    # clean lysogenic
    remove_na = ~((lysogenic_df['proteinID'].isna()) |
                   (lysogenic_df['proteinID'] == 'GWAS proteins:'))
    lysogenic_df = lysogenic_df.loc[remove_na].reset_index(drop=True)
    ugly_string = lysogenic_df['specificity'].str.contains('-')
    lysogenic_df.loc[ugly_string, 'specificity'] = 'NO_KL'

    clean_expression_dict = {
        'high': 'PRODUCED', 'medium': 'PRODUCED',
        'low': 'PRODUCED', 'none': 'NOT_PRODUCED'
    }
    lysogenic_df['expression'] = lysogenic_df['expression'].replace(clean_expression_dict)

    # assign groups
    is_produced = (lysogenic_df['expression'] == 'PRODUCED')
    is_active = ~(lysogenic_df['specificity'].str.contains('NO_KL'))
    is_genscript = lysogenic_df['proteinID'].isin(genscript_list)
    is_zdk = ~(lysogenic_df['proteinID'].isin(genscript_list))

    zdk_not_produced = is_zdk & ~is_produced
    zdk_produced_active = is_zdk & is_produced & is_active
    zdk_produced_inactive = is_zdk & is_produced & ~is_active

    genscript_not_produced = is_genscript & ~is_produced
    genscript_produced_active = is_genscript & is_produced & is_active
    genscript_produced_inactive = is_genscript & is_produced & ~is_active

    lysogenic_df.loc[zdk_not_produced, 'group'] = 'lysogenic_zdk_not_produced'
    lysogenic_df.loc[zdk_produced_active, 'group'] = 'lysogenic_zdk_produced_active'
    lysogenic_df.loc[zdk_produced_inactive, 'group'] = 'lysogenic_zdk_produced_inactive'

    lysogenic_df.loc[genscript_not_produced, 'group'] = 'lysogenic_genscript_not_produced'
    lysogenic_df.loc[genscript_produced_active, 'group'] = 'lysogenic_genscript_produced_active'
    lysogenic_df.loc[genscript_produced_inactive, 'group'] = 'lysogenic_genscript_produced_inactive'

    lysogenic_df['report'] = (lysogenic_df['proteinID'] + '_' +
                              lysogenic_df['specificity'] + '_' +
                              lysogenic_df['expression'] + '_HOST_' +
                              lysogenic_df['K_locus_host'])

    # --- lytic (S1) ---
    lytic_df['group'] = 'lytic'
    lytic_df['report'] = lytic_df['proteinID']

    # combine
    recombinant_df = pd.concat([lysogenic_df, lytic_df])
    recombinant_df = recombinant_df[recombinant_cols_step1]
    print('Done!')

    # ------------------------------------------------------------------
    # STEP 2 – Combine recombinant enzymes with GWAS predictions (S3)
    # ------------------------------------------------------------------
    print('Combining predictions and enzymes...', end='\t')

    recombinant_cols_step2 = [
        'proteinID', 'group', 'expression',
        'specificity', 'K_locus_host', 'seq'
    ]
    prediction_cols = ['locus', 'PC', 'prediction_strength', 'seq']

    bacteria_df = pd.read_csv(bacteria_table, sep='\t')
    prediction_df = pd.read_excel(S3_Table_xlsx, sheet_name='GWAS predictions')

    # clean recombinant
    recombinant_df = recombinant_df[recombinant_cols_step2].rename({'group': 'source'}, axis=1)
    recombinant_df['expression'] = recombinant_df['expression'].fillna('PRODUCED')
    recombinant_df.loc[recombinant_df['source'].str.contains('zdk'), 'source'] = 'PROPHAGE_ZDKLAB'
    recombinant_df.loc[recombinant_df['source'].str.contains('genscript'), 'source'] = 'GENSCRIPT'
    recombinant_df.loc[recombinant_df['source'].str.contains('lytic'), 'source'] = 'LITERATURE_SEARCH'

    recombinant_df['activity'] = 'ACTIVE'
    recombinant_df.loc[recombinant_df['specificity'].str.contains('NO_KL'), 'activity'] = 'INACTIVE'
    recombinant_df.loc[recombinant_df['expression'].str.contains('NOT_PRODUCED'), 'activity'] = np.nan
    recombinant_df.loc[recombinant_df['expression'].str.contains('NOT_PRODUCED'), 'specificity'] = np.nan
    recombinant_df.loc[recombinant_df['activity'] == 'INACTIVE', 'specificity'] = np.nan

    # fix literature search specificity naming: K → KL, KLL → KL, KLN → KN
    recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'] = \
        recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'].str.replace('K', 'KL')

    recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'] = \
        recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'].str.replace('KLL', 'KL')

    recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'] = \
        recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'].str.replace('KLN', 'KN')

    # clean literature search proteinID (remove suffix after last underscore)
    recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'proteinID'] = \
        recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'proteinID'].str.rsplit('_', n=1).str[0]

    # clean prediction
    prediction_df = prediction_df[prediction_cols].rename(
        {'locus': 'specificity', 'PC': 'proteinID'}, axis=1)
    prediction_df = prediction_df.dropna()
    prediction_df['source'] = 'PREDICTION'

    # concat
    depos_df = pd.concat([recombinant_df, prediction_df])
    depos_df['assing_K_locus_host_when_specificity_missing'] = \
        depos_df['specificity'].fillna(depos_df['K_locus_host'])

    # strains per K locus
    kaspah_ref_df = (bacteria_df
                     .query('collection == "kaspah_complete"')
                     .groupby("MGG_K_locus").size().reset_index()
                     .rename({'MGG_K_locus': 'specificity', 0: '# KASPAH-REF'}, axis=1))
    ksc_gwas_df = (bacteria_df
                   .groupby("MGG_K_locus").size().reset_index()
                   .rename({'MGG_K_locus': 'specificity', 0: '# GWAS_KSC'}, axis=1))
    strains_availability = kaspah_ref_df.merge(ksc_gwas_df, on='specificity', how='outer')
    strains_availability['assing_K_locus_host_when_specificity_missing'] = strains_availability['specificity']
    strains_availability = strains_availability.drop('specificity', axis=1)

    # merge
    depos_df = depos_df.merge(strains_availability,
                              on='assing_K_locus_host_when_specificity_missing', how='left')

    # save
    final_cols = [
        'proteinID', 'source', 'expression',
        'assing_K_locus_host_when_specificity_missing',
        'specificity', 'K_locus_host', '# KASPAH-REF', '# GWAS_KSC',
        'activity', 'prediction_strength', 'seq'
    ]
    depos_df[final_cols].to_csv(predictions_and_enzymes, sep='\t', index=False)
    print('Done!')


def FIGURE4_PANELA(cell_number_color="#000000"):
    """
    Revised FIGURE4_PANELA visualization that:
      - Uses the provided layout configuration,
      - Plots only the union of K loci coming from:
           (i) unique KASPAH_REF,
           (ii) custom GWAS K loci,
           (iii) unique PROPHAGE active/inactive (after splitting dual specificities),
           (iv) unique lytic phage depolymerases (after splitting dual specificities),
           (v) unique GENSCRIPT active/inactive (after splitting dual specificities).
      - Skips NOT_PRODUCED proteins,
      - Draws split cells only when both active and inactive counts are present,
      - And for "# Kp complex SCs" colors only cells with values ≥ 10 (values less than 10 are white).
      
    New features controlled by parameters:
      1. remove_kn_kloci: Remove any K locus that starts with "KN" from the plot.
      2. remove_inactive: Remove inactive protein counts for PROPHAGE and GENSCRIPT rows.
      
    Row order (merged version with subcells) is:
      Row 0: "KASPAH-REF isolates" 
      Row 1: "Sequence Clusters (SCs)"
      Row 2: "active proteins based on GWAS prediction"  (prediction row)
      Row 3: "active proteins based on manual prediction"  (split: PROPHAGE active/inactive, but inactive removed if flag is True)
      Row 4: "genscript proteins"                        (split: GENSCRIPT active/inactive, but inactive removed if flag is True)
      Row 5: "lytic proteins"
    """

    import os
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import re


    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42



    # ------------------------------------------------------------
    # Feature control flags (set via parameters)
    # ------------------------------------------------------------
    # Set these in your params (default: True)
    remove_kn_kloci = True
    remove_inactive = True

    # ------------------------------------------------------------
    # Correction function: fix naming issues for K locus strings.
    # ------------------------------------------------------------
    def correct_k_locus(k):
        corrections = {'K14': 'KL14', 'K3': 'KL3'}
        return corrections.get(k, k)

    # ------------------------------------------------------------
    # Layout configuration (do not change these values)
    # ------------------------------------------------------------
    layout_config = {
        "num_rows": 6,
        "margin": 0.3,
        "header_space": 0.1,
        # "labels_fontsize": 5.25,
        "cells_fontsize": 4.5,
        "target_fig_width": 7.5,
        "target_fig_height": 1.4,
        "cell_line_width": 0.35,
        "X_labels_fontsize": 6.25,
        "Y_labels_fontsize": 5.25
    }

    # ------------------------------------------------------------
    # Color configuration and row definitions.
    # ------------------------------------------------------------
    color_config = {
        "KASPAH_REF": "#f06c6c",          # remains unchanged
        "GWAS_KSC": "#FFBD4D",             # used if value ≥ 10
        "PROPHAGE_ACTIVE": "#fa61bd",      # dark shade for active (PROPHAGE_ZDKLAB)
        "PROPHAGE_INACTIVE": "#f5bfe0",     # light shade for inactive (PROPHAGE_ZDKLAB)
        "PREDICTION_PERFECT": "#85C187",    # remains the same
        "PREDICTION_GOOD": "#DBEFDC",       # remains the same
        "LITERATURE_SEARCH": "#bc96d9",     # for lytic phage depolymerases
        "GENSCRIPT_ACTIVE": "#6caed7",      # lighter blue for genscript active
        "GENSCRIPT_INACTIVE": "#dfebf8"      # light blue for genscript inactive
    }

    row_config = {
        0: {"data_key": "KASPAH_REF", "label": "KASPAH-REF isolates", "plot_type": "single"},
        1: {"data_key": "GWAS_KSC", "label": "Sequence Clusters (SCs)", "plot_type": "single"},
        2: {"data_key": "PREDICTION", "label": "GWAS prediction", "plot_type": "prediction",
            "categories": ["PREDICTION_PERFECT", "PREDICTION_GOOD"]},
        3: {"data_key": "PROPHAGE", "label": "Active based on manual prediction", "plot_type": "split",
            "categories": ["ACTIVE", "INACTIVE"], "color_keys": ["PROPHAGE_ACTIVE", "PROPHAGE_INACTIVE"]},
        4: {"data_key": "GENSCRIPT", "label": "Active based on GWAS prediction", "plot_type": "split",
            "categories": ["ACTIVE", "INACTIVE"], "color_keys": ["GENSCRIPT_ACTIVE", "GENSCRIPT_INACTIVE"]},
        5: {"data_key": "LITERATURE_SEARCH", "label": "Lytic proteins", "plot_type": "single"}
    }
    
    # ------------------------------------------------------------
    # Data Loading
    # ------------------------------------------------------------

    global custom_gwas_kloci
    global bacteria_table
    global predictions_and_enzymes
    global panelA_pdf_path

    bacteria_df = pd.read_csv(bacteria_table, sep='\t')
    predictions_and_enzymes_df = pd.read_csv(predictions_and_enzymes, sep='\t')
    
    # ------------------------------------------------------------
    # Pre-processing: apply correction and dual-specificity splitting
    # ------------------------------------------------------------
    def process_dual_specificity(df):
        if 'specificity' in df.columns:
            df = df.copy()
            df['assing_specificity'] = df.apply(
                lambda row: [correct_k_locus(s.strip()) for s in row['specificity'].split('/')]
                    if pd.notnull(row['specificity']) and "/" in row['specificity']
                    else [correct_k_locus(row['assing_K_locus_host_when_specificity_missing'])],
                axis=1)
            df = df.explode('assing_specificity')
        else:
            df['assing_specificity'] = df['assing_K_locus_host_when_specificity_missing'].apply(correct_k_locus)
        return df

    # Define predicted_df, literature_df, prophage_df, and genscript_df properly.
    predicted_df = predictions_and_enzymes_df.query('source == "PREDICTION"')
    predicted_df = process_dual_specificity(predicted_df)
    
    literature_df = predictions_and_enzymes_df.query('source == "LITERATURE_SEARCH"')
    literature_df = process_dual_specificity(literature_df)
    
    prophage_df = predictions_and_enzymes_df.query('source == "PROPHAGE_ZDKLAB"')
    prophage_df = process_dual_specificity(prophage_df)
    
    genscript_df = predictions_and_enzymes_df.query('source == "GENSCRIPT"')
    genscript_df = process_dual_specificity(genscript_df)
    
    # correct the K locus names in genscript if needed.
    genscript_df["specificity"] = genscript_df["specificity"].apply(correct_k_locus)
    
    # ------------------------------------------------------------
    # Compute counts for each category.
    # ------------------------------------------------------------
    # KASPAH_REF counts from PROPHAGE and PREDICTION datasets.
    kaspa_ref_series = pd.concat([prophage_df, predicted_df])\
                          .groupby('assing_specificity')["# KASPAH-REF"].first().fillna(0)
    
    # GWAS counts.
    gwas_series = bacteria_df.groupby("MGG_K_locus")["MGG_SC"].nunique().fillna(0)
    
    # PROPHAGE overexpressed proteins: active and inactive.
    prophage_active_series = prophage_df.query('expression == "PRODUCED" and activity == "ACTIVE"')\
                                        .groupby('assing_specificity').size().fillna(0)
    prophage_inactive_series = prophage_df.query('expression == "PRODUCED" and activity == "INACTIVE"')\
                                          .groupby('assing_specificity').size().fillna(0)
    
    # Predicted proteins from predicted_df.
    predicted_group = predicted_df.groupby('assing_specificity')
    predicted_dict = {}
    for kl, group in predicted_group:
        counts = group['prediction_strength'].value_counts().to_dict()
        predicted_dict[kl] = counts
    
    # Lytic phage depolymerases counts.
    literature_series = literature_df.groupby('assing_specificity').size().fillna(0)
    
    # GENSCRIPT proteins: active and inactive.
    genscript_active_series = genscript_df.query('expression == "PRODUCED" and activity == "ACTIVE"')\
                                          .groupby('assing_specificity').size().fillna(0)
    genscript_inactive_series = genscript_df.query('expression == "PRODUCED" and activity == "INACTIVE"')\
                                            .groupby('assing_specificity').size().fillna(0)
    
    # ------------------------------------------------------------
    # Build union of all K loci to plot.
    # ------------------------------------------------------------
    union_kloci = set(kaspa_ref_series.index) | set(custom_gwas_kloci) | \
                  set(prophage_active_series.index) | set(prophage_inactive_series.index) | \
                  set(literature_series.index) | \
                  set(genscript_active_series.index) | set(genscript_inactive_series.index)
    
    # Feature 1: Remove any K locus that starts with "KN" if flag is True.
    if remove_kn_kloci:
        used_kloci = {kl for kl in union_kloci if not kl.startswith("KN")}
    else:
        used_kloci = union_kloci

    # Build input data dictionary for each K locus.
    input_dict = {}
    for kl in used_kloci:
        input_dict[kl] = {
            "KASPAH_REF": int(kaspa_ref_series.get(kl, 0)),
            "GWAS_KSC": int(gwas_series.get(kl, 0)),
            "PROPHAGE": {
                "ACTIVE": int(prophage_active_series.get(kl, 0)),
                "INACTIVE": int(prophage_inactive_series.get(kl, 0))
            },
            "PREDICTION": predicted_dict.get(kl, {}),
            "LITERATURE_SEARCH": int(literature_series.get(kl, 0)),
            "GENSCRIPT": {
                "ACTIVE": int(genscript_active_series.get(kl, 0)),
                "INACTIVE": int(genscript_inactive_series.get(kl, 0))
            }
        }
    
    # Ensure every K locus from custom GWAS is present.
    for kl in custom_gwas_kloci:
        if kl not in input_dict:
            input_dict[kl] = {
                "KASPAH_REF": 0,
                "GWAS_KSC": 0,
                "PROPHAGE": {"ACTIVE": 0, "INACTIVE": 0},
                "PREDICTION": {},
                "LITERATURE_SEARCH": 0,
                "GENSCRIPT": {"ACTIVE": 0, "INACTIVE": 0}
            }
    
    # Feature 2: Remove inactive proteins for PROPHAGE and GENSCRIPT rows if flag is True.
    if remove_inactive:
        for kl in input_dict:
            input_dict[kl]["PROPHAGE"]["INACTIVE"] = 0
            input_dict[kl]["GENSCRIPT"]["INACTIVE"] = 0

    # ------------------------------------------------------------
    # Sort the K loci keys.
    # ------------------------------------------------------------
    def sort_key(key):
        prefix = key[:2]
        if prefix == "KL":
            prefix_order = 0
        elif prefix == "KN":
            prefix_order = 1
        else:
            prefix_order = 2
        match = re.search(r'\d+', key)
        num = int(match.group()) if match else float('inf')
        return (prefix_order, num)
    
    sorted_keys = sorted(input_dict.keys(), key=sort_key)
    num_cols = len(sorted_keys)
    
    # ------------------------------------------------------------
    # Figure layout settings from layout_config.
    # ------------------------------------------------------------
    margin = layout_config["margin"]
    header_space = layout_config["header_space"]
    target_fig_width = layout_config["target_fig_width"]
    target_fig_height = layout_config["target_fig_height"]
    X_labels_fontsize = layout_config["X_labels_fontsize"]
    Y_labels_fontsize = layout_config["Y_labels_fontsize"]
    cells_fontsize = layout_config["cells_fontsize"]
    cell_line_width = layout_config["cell_line_width"]
    num_rows = layout_config["num_rows"]
    
    cell_width = (target_fig_width - 2 * margin) / num_cols
    cell_height = (target_fig_height - 2 * margin) / num_rows
    
    fig, ax = plt.subplots(figsize=(target_fig_width, target_fig_height))
    ax.set_xlim(0, target_fig_width)
    ax.set_ylim(0, target_fig_height)
    ax.axis('off')
    
    # ------------------------------------------------------------
    # Draw column headers (centered above each cell).
    # ------------------------------------------------------------
    header_y = target_fig_height - margin + header_space
    for i, k_locus in enumerate(sorted_keys):
        x = margin + i * cell_width + cell_width / 2
        ax.text(x, header_y, k_locus, ha='center', va='bottom',
                fontsize=X_labels_fontsize, fontweight='bold', rotation=60)
    
    # ------------------------------------------------------------
    # Loop through rows (using row_config order) and columns.
    # ------------------------------------------------------------
    for row_index in range(num_rows):
        config = row_config[row_index]
        data_key = config["data_key"]
        plot_type = config["plot_type"]
        for col, k_locus in enumerate(sorted_keys):
            x = margin + col * cell_width
            y = margin + (num_rows - 1 - row_index) * cell_height
            
            cell_value = input_dict[k_locus].get(data_key, 0)
            
            if plot_type == "single":
                count = cell_value
                # Special handling for "# Kp complex SCs" (GWAS_KSC)
                if data_key == "GWAS_KSC":
                    cell_color = color_config.get(data_key, "white") if count >= 10 else "white"
                else:
                    cell_color = color_config.get(data_key, "white") if count > 0 else "white"
                rect = patches.Rectangle((x, y), cell_width, cell_height,
                                         facecolor=cell_color, edgecolor="black", linewidth=cell_line_width)
                ax.add_patch(rect)

                if count > 0:
                    ax.text(x + cell_width/2, y + cell_height/2, str(count),
                        ha='center', va='center', color=cell_number_color,
                        fontsize=cells_fontsize, fontweight='bold')
            
            elif plot_type == "split":
                cat_counts = cell_value
                count_active = cat_counts.get("ACTIVE", 0)
                count_inactive = cat_counts.get("INACTIVE", 0)
                total = count_active + count_inactive
                
                if total == 0:
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor="white", edgecolor="black", linewidth=cell_line_width)
                    ax.add_patch(rect)
                elif (count_active > 0 and count_inactive == 0) or (count_inactive > 0 and count_active == 0):
                    if count_active > 0:
                        fill_color = color_config[config["color_keys"][0]]
                        count = count_active
                    else:
                        fill_color = color_config[config["color_keys"][1]]
                        count = count_inactive
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor=fill_color, edgecolor="black", linewidth=cell_line_width)
                    ax.add_patch(rect)
                    ax.text(x + cell_width/2, y + cell_height/2, str(count),
                            ha='center', va='center', color=cell_number_color,
                            fontsize=cells_fontsize, fontweight='bold')
                else:
                    subcell_width = cell_width / 2
                    # Left subcell for ACTIVE
                    sub_rect = patches.Rectangle((x, y), subcell_width, cell_height,
                                                 facecolor=color_config[config["color_keys"][0]],
                                                 edgecolor="black", linewidth=cell_line_width)
                    ax.add_patch(sub_rect)
                    ax.text(x + subcell_width/2, y + cell_height/2, str(count_active),
                            ha='center', va='center', color=cell_number_color,
                            fontsize=cells_fontsize, fontweight='bold')
                    # Right subcell for INACTIVE
                    sub_x = x + subcell_width
                    sub_rect = patches.Rectangle((sub_x, y), subcell_width, cell_height,
                                                 facecolor=color_config[config["color_keys"][1]],
                                                 edgecolor="black", linewidth=cell_line_width)
                    ax.add_patch(sub_rect)
                    ax.text(sub_x + subcell_width/2, y + cell_height/2, str(count_inactive),
                            ha='center', va='center', color=cell_number_color,
                            fontsize=cells_fontsize, fontweight='bold')
            
            elif plot_type == "prediction":
                pred_counts = cell_value
                if not pred_counts:
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor="white", edgecolor="black", linewidth=cell_line_width)
                    ax.add_patch(rect)
                elif len(pred_counts) == 1:
                    pred_category, count = list(pred_counts.items())[0]
                    cell_color = color_config["PREDICTION_PERFECT"] if pred_category.lower() == "strong" else color_config["PREDICTION_GOOD"]
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor=cell_color, edgecolor="black", linewidth=cell_line_width)
                    ax.add_patch(rect)
                    ax.text(x + cell_width/2, y + cell_height/2, str(count),
                            ha='center', va='center', color=cell_number_color,
                            fontsize=cells_fontsize, fontweight='bold')
                else:
                    num_subcells = len(pred_counts)
                    subcell_width = cell_width / num_subcells
                    sorted_preds = sorted(pred_counts.items(), key=lambda x: 0 if x[0].lower() == "strong" else 1)
                    for i, (pred_category, count) in enumerate(sorted_preds):
                        sub_x = x + i * subcell_width
                        mapped_category = "PREDICTION_PERFECT" if pred_category.lower() == "strong" else "PREDICTION_GOOD"
                        sub_rect = patches.Rectangle((sub_x, y), subcell_width, cell_height,
                                                     facecolor=color_config[mapped_category],
                                                     edgecolor="black", linewidth=cell_line_width)
                        ax.add_patch(sub_rect)
                        ax.text(sub_x + subcell_width/2, y + cell_height/2, str(count),
                                ha='center', va='center', color=cell_number_color,
                                fontsize=cells_fontsize, fontweight='bold')
    
    # ------------------------------------------------------------
    # Draw row labels (centered in each row).
    # ------------------------------------------------------------
    for row_index in range(num_rows):
        config = row_config[row_index]
        y = margin + (num_rows - 1 - row_index) * cell_height + cell_height/2
        ax.text(margin - 0.075, y, config["label"], ha='right', va='center',
                fontsize=Y_labels_fontsize, fontweight='bold')
    
    plt.savefig(panelA_pdf_path, format='pdf', bbox_inches='tight', pad_inches=0)
    plt.close(fig)


def FIGURE4_PANELB(min_cov=0.5, min_pident=0, max_evalue=1e-3, use_coverage=False):
    import pandas as pd
    import numpy as np
    import subprocess
    from pathlib import Path

    def generate_edges(
        nodes_df: pd.DataFrame,
        tmp_dir: str,
        edges_file: str,
        min_cov: float = 0.5,
        min_pident: float = 20,
        max_evalue: float = 1e-3,
        seqidcol='proteinID',
        use_coverage: bool = False,
        intervals: list = None  # Expected as [bottom, mid, top]
    ) -> pd.DataFrame:
        """
        Generate a DataFrame of edges between proteins based on BLASTp similarity.

        Parameters
        ----------
        nodes_df : pd.DataFrame
            Must include at least columns: ['proteinID', 'seq'].
        tmp_dir : str
            Directory to store intermediate BLAST files.
        edges_file : str
            Path to final edges TSV output file.
        min_cov : float, optional
            Minimum coverage threshold. Default=0.5.
        min_pident : float, optional
            Minimum percent-identity threshold (ignored if use_coverage is True). Default=20.
        max_evalue : float, optional
            Maximum e-value threshold. Default=1e-3.
        seqidcol : str, optional
            Column name for sequence IDs.
        use_coverage : bool, optional
            If True, filtering and interval assignment use coverage instead of percent identity.
        intervals : list, optional
            List of three numbers defining the interval boundaries.
            For pident: [bottom, mid, top] (e.g., [min_pident, 50, 80])
            For coverage: [bottom, mid, top] in percentage (e.g., [min_cov*100, 70, 90])
        
        Returns
        -------
        edges_df : pd.DataFrame
            The filtered edges, ready for downstream usage.
        """
        from pathlib import Path
        import subprocess
        import numpy as np

        # Ensure tmp_dir is a Path object.
        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(exist_ok=True, parents=True)

        # Define working file paths.
        sequences_file = tmp_dir / 'sequences.fasta'
        blastp_raw     = tmp_dir / 'blastp_raw.tsv'
        blast_not_filtered = tmp_dir / 'blast_not_filtered.tsv'
        blast_processed    = tmp_dir / 'blast_processed.tsv'
        
        # Build FASTA from nodes_df.
        fasta_strs = []
        for _, row in nodes_df.iterrows():
            header = f">{row[seqidcol]}\n"
            seq    = f"{row['seq']}\n"
            fasta_strs.append(header + seq)
        with open(sequences_file, 'w') as f:
            f.write(''.join(fasta_strs))

        # Create BLAST database and run BLASTp.
        makeblastdb_cmd = [
            'makeblastdb',
            '-in', str(sequences_file),
            '-dbtype', 'prot'
        ]
        blastp_cmd = [
            'blastp',
            '-query', str(sequences_file),
            '-db', str(sequences_file),
            '-out', str(blastp_raw),
            '-outfmt', '6 qseqid sseqid evalue pident bitscore qlen qstart qend slen sstart send'
        ]
        subprocess.run(makeblastdb_cmd, check=True)
        subprocess.run(blastp_cmd, check=True)

        # Load raw BLAST results.
        edges_df = pd.read_csv(blastp_raw, sep='\t', header=None)
        cols = [
            'query', 'target', 'evalue', 'pident', 'bitscore',
            'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send'
        ]
        edges_df.columns = cols

        # Calculate coverage for query and subject.
        edges_df['qcov'] = np.round((edges_df['qend'] - edges_df['qstart'] + 1) / edges_df['qlen'], 2)
        edges_df['scov'] = np.round((edges_df['send'] - edges_df['sstart'] + 1) / edges_df['slen'], 2)

        # Deduplicate self-hits and symmetrical hits.
        edges_df['query-target'] = edges_df.apply(
            lambda row: '-'.join(sorted([row['query'], row['target']])), axis=1
        )
        edges_df = edges_df.drop_duplicates(subset='query-target')

        # Save unfiltered BLAST results.
        edges_df.to_csv(blast_not_filtered, sep='\t', index=False)

        # Apply filters.
        is_self = edges_df['query'] == edges_df['target']
        filt_cov = (edges_df['scov'] >= min_cov) & (edges_df['qcov'] >= min_cov)
        filt_evalue = edges_df['evalue'] <= max_evalue

        if use_coverage:
            # Use coverage as the metric.
            filt = filt_cov & filt_evalue & (~is_self)
        else:
            # Filter on both percent identity and coverage.
            filt_pident = edges_df['pident'] >= min_pident
            filt = filt_cov & filt_pident & filt_evalue & (~is_self)
        edges_df = edges_df.loc[filt]

        # Assign intervals based on the chosen metric.
        if use_coverage:
            # Calculate combined coverage (minimum of query and subject coverage, in percentage).
            edges_df['coverage'] = edges_df[['qcov', 'scov']].min(axis=1) * 100
            # Set default intervals if not provided.
            if intervals is None:
                intervals = [min_cov * 100, 70, 90]
            bottom, mid, top = intervals[0], intervals[1], intervals[2]
            remove_low_cov = edges_df['coverage'] >= bottom
            filt_low = (edges_df['coverage'] >= bottom) & (edges_df['coverage'] < mid)
            filt_med = (edges_df['coverage'] >= mid) & (edges_df['coverage'] < top)
            filt_high = edges_df['coverage'] >= top

            edges_df = edges_df.loc[remove_low_cov]
            edges_df.loc[filt_low, 'coverage_interval'] = f'low ({bottom}-{mid})'
            edges_df.loc[filt_med, 'coverage_interval'] = f'medium ({mid}-{top})'
            edges_df.loc[filt_high, 'coverage_interval'] = f'high ({top}-100)'
        else:
            # Use percent identity.
            if intervals is None:
                intervals = [min_pident, 50, 80]
            bottom, mid, top = intervals[0], intervals[1], intervals[2]
            remove_low_pident = edges_df['pident'] >= bottom
            filt_low = (edges_df['pident'] >= bottom) & (edges_df['pident'] < mid)
            filt_med = (edges_df['pident'] >= mid) & (edges_df['pident'] < top)
            filt_high = edges_df['pident'] >= top

            edges_df = edges_df.loc[remove_low_pident]
            edges_df.loc[filt_low, 'pident_interval'] = f'low ({bottom}-{mid})'
            edges_df.loc[filt_med, 'pident_interval'] = f'medium ({mid}-{top})'
            edges_df.loc[filt_high, 'pident_interval'] = f'high ({top}-100)'

        # Save processed BLAST results.
        edges_df.to_csv(blast_processed, sep='\t', index=False)

        # Add singletons (nodes with no edges).
        all_nodes = list(nodes_df[seqidcol].unique())
        blast_nodes = set(edges_df['query'].tolist() + edges_df['target'].tolist())
        missing_nodes = [node for node in all_nodes if node not in blast_nodes]
        if missing_nodes:
            singletons_df = pd.DataFrame({
                'query': missing_nodes,
                'target': missing_nodes,
                'pident_interval': ['high ({0}-100)'.format(top)] * len(missing_nodes) if not use_coverage 
                                    else ['high ({}-100)'.format(top)] * len(missing_nodes)
            })
            edges_df = pd.concat([edges_df, singletons_df], ignore_index=True)

        edges_df['interaction'] = 'A'
        edges_df.to_csv(edges_file, sep='\t', index=False)
        return edges_df

    

    def assign_categories_and_labels(df):
        """
        Assigns a category and final label to each row in the DataFrame based on:
          - source (e.g., PROPHAGE_ZDKLAB, LITERATURE_SEARCH, GENSCRIPT, PREDICTION)
          - expression and activity (for PROPHAGE_ZDKLAB and GENSCRIPT)
          - prediction_strength (for PREDICTION)
        
        Final label construction:
          - For categories that map to PROTEINID_K_LOCUS_HOST:
                label = proteinID + "_" + K_locus_host
          - For categories that map to PROTEINID_SPECIFICITY:
                label = proteinID + "_" + specificity
          - For LITERATURE_SEARCH, label = proteinID directly
        """
        
        # Function to determine category for each row
        def get_category(row):
            source = row.get('source', '')
            if source == "PROPHAGE_ZDKLAB":
                if row.get('expression', '') == "NOT_PRODUCED":
                    return "ZDK_NOT_PRODUCED"
                elif row.get('expression', '') == "PRODUCED":
                    if row.get('activity', '') == "ACTIVE":
                        return "ZDK_ACTIVE"
                    elif row.get('activity', '') == "INACTIVE":
                        return "ZDK_INACTIVE"
            elif source == "LITERATURE_SEARCH":
                return "LITERATURE_SEARCH"
            elif source == "GENSCRIPT":
                if row.get('activity', '') == "ACTIVE":
                    return "GENSCRIPT_ACTIVE"
                elif row.get('activity', '') == "INACTIVE":
                    return "GENSCRIPT_INACTIVE"
            elif source == "PREDICTION":
                ps = str(row.get('prediction_strength', ''))
                if ps == "strong":
                    return "PREDICTION_PERFECT"
                elif ps == "likely":
                    return "PREDICTION_GOOD"
            return None

        # Apply category assignment
        df['category'] = df.apply(get_category, axis=1)

        # Function to build the final label by concatenating columns.
        def get_final_label(row):
            cat = row.get('category', None)
            protein_id = str(row.get('proteinID', ''))
            if cat == "LITERATURE_SEARCH":
                # This category already has specificity in the proteinID.
                return f"{protein_id}_{row.get('specificity', '')}"
            elif cat in ["ZDK_NOT_PRODUCED", "ZDK_INACTIVE", "GENSCRIPT_INACTIVE"]:
                # Use K_locus_host column
                return f"{protein_id}_HOST_{row.get('K_locus_host', '')}"
            elif cat in ["ZDK_ACTIVE", "GENSCRIPT_ACTIVE", "PREDICTION_PERFECT", "PREDICTION_GOOD"]:
                # Use specificity column
                return f"{protein_id}_{row.get('specificity', '')}"
            return None

        # Apply final label construction
        df['label1'] = df.apply(get_final_label, axis=1)
        return df

    def get_label2(label):
        return label.split('_')[-1]
        # if 'HOST' in label:
        #     return '_'.join(label.split('_')[-2:])
        # else: 
        #     return label.split('_')[-1]

    # paths
    global predictions_and_enzymes
    global nodes_outfile
    global edges_outfile
    global tmp_blast_dir

    # read
    predictions_and_enzymes_df = pd.read_csv(predictions_and_enzymes, sep='\t')

    # nodes: generate and save
    cols = ['query', 'proteinID', 'source', 'expression', 'category', 'label1', 'label2', 'seq']
    predictions_and_enzymes_df = assign_categories_and_labels(predictions_and_enzymes_df)

    predictions_and_enzymes_df.index = predictions_and_enzymes_df.index + 1
    predictions_and_enzymes_df['query'] = predictions_and_enzymes_df.index.astype(str) + '_' + predictions_and_enzymes_df['label1']
    predictions_and_enzymes_df['label2'] = predictions_and_enzymes_df['label1'].apply(get_label2)

    predictions_and_enzymes_df[cols].to_csv(nodes_outfile, sep='\t', index=False)

    # edges: generate and save
    edges_df = generate_edges(predictions_and_enzymes_df, tmp_blast_dir, edges_outfile,
                            min_cov=min_cov, min_pident=min_pident, max_evalue=max_evalue,
                            seqidcol='query', use_coverage=True,
                            intervals=[00, 50, 80])


def FIGURE4_PANELC():
    """
    Create a visualization with one protein per row and K loci as columns.
    
    Steps:
      1. Extract all unique K loci from both 'K_locus_host' and 'specificity' columns.
      2. Strip the 'KL' prefix and sort them numerically in ascending order.
      3. Create a matrix with one row per protein (labeled on the left) and each sorted K locus as columns.
         For each protein:
           - If the protein is active against the same K locus as its host, color that cell green.
           - If the protein is active against a different K locus than its host, mark two cells:
             the active one in green and the host cell in red.
      4. Sort the protein rows by the lowest active K locus number.
      5. Move (append) all columns that are “host‐only” (i.e. no protein is active on that column) to the right.
      6. If swap_axes is True, swap the X and Y axes (proteins become columns and K loci become rows).
      
    Colors used:
       Green: #0E470E
       Red:   #f06c6c
    """
    
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    # Use TrueType fonts in PDF for easier editing in Illustrator.
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42

    # ---- Feature Switch: swap_axes ----
    # If set True (in params), proteins become columns and K loci become rows.
    swap_axes = False

    # --- Data Loading ---
    global predictions_and_enzymes
    global panelC_pdf_path
    
    df = pd.read_csv(predictions_and_enzymes, sep='\t')
    # Use only proteins that are active
    df = df.query('source == "PROPHAGE_ZDKLAB" and activity == "ACTIVE"')
    df = df[['proteinID', 'K_locus_host', 'specificity']]
    
    # --- Step 1: Extract Unique K Loci ---
    unique_kloci = set()
    for host in df['K_locus_host']:
        unique_kloci.add(host.strip())
    for spec in df['specificity']:
        for k in spec.split('/'):
            unique_kloci.add(k.strip())
    
    # --- Step 2: Strip 'KL' and sort ascending ---
    def get_num(k):
        k = k.strip()
        if k.startswith("KL"):
            try:
                return int(k[2:])
            except:
                return float('inf')
        else:
            try:
                return int(k)
            except:
                return float('inf')
    
    sorted_kloci = sorted(unique_kloci, key=get_num)
    # Use the numeric part as string for internal computations.
    sorted_kloci_str = [str(get_num(k)) for k in sorted_kloci]
    
    # --- Step 3: Build Data for Each Protein ---
    protein_data = []
    for idx, row in df.iterrows():
        protein = row['proteinID']
        host = row['K_locus_host'].strip()
        host_num = str(get_num(host))
        active_list = [k.strip() for k in row['specificity'].split('/')]
        active_nums = [str(get_num(k)) for k in active_list]
        protein_data.append((protein, host_num, active_nums))
    
    # --- Step 3a/3b: Determine Cell Colors per Protein Row ---
    GREEN = "#0E470E"
    RED = "#f06c6c"
    cell_colors = {}
    for protein, host_num, active_nums in protein_data:
        if not active_nums:
            continue
        if host_num in active_nums:
            cell_colors[(protein, host_num)] = GREEN
        else:
            for act in active_nums:
                cell_colors[(protein, act)] = GREEN
            cell_colors[(protein, host_num)] = RED

    # --- Step 4: Sort Protein Rows ---
    def protein_sort_key(item):
        if item[2]:
            return min(int(x) for x in item[2])
        return float('inf')
    sorted_protein_data = sorted(protein_data, key=protein_sort_key)
    sorted_proteins = [protein for (protein, _, _) in sorted_protein_data]
    
    # --- Plot Configuration Dictionary ---
    # Controls overall figure size, cell dimensions, and label fonts/offsets.
    plot_config = {
        "figsize": (len(sorted_kloci_str) * 0.5 + 3, len(sorted_proteins) * 0.5 + 3),
        "cell_width": 1,
        "cell_height": 1,
        "protein_label_fontsize": 12,
        "protein_label_offset": -0.2,  # for protein labels in non-swap mode (left side)
        "k_label_fontsize": 12,
        "k_label_rotation": 60,
        "k_label_offset": 0.1         # vertical offset for K loci headers in non-swap mode (top)
    }
    
    # --- New Feature: Sort Columns ---
    # Separate sorted K loci columns into those that receive an active (GREEN) and those that are host-only.
    active_set = set()
    for (protein, col), color in cell_colors.items():
        if color == GREEN:
            active_set.add(col)
    active_columns = [k for k in sorted_kloci_str if k in active_set]
    host_only_columns = [k for k in sorted_kloci_str if k not in active_set]
    # Final column order: active columns first, then host-only columns.
    final_columns = active_columns + host_only_columns
    
    # --- Establish Grid Dimensions & Label Order Based on swap_axes ---
    if not swap_axes:
        # Normal mode: rows = proteins, columns = K loci.
        row_labels = sorted_proteins
        col_labels = final_columns
        # For headers, protein labels remain as is; K loci get "KL" prefix.
        row_prefix = ""          # protein labels: no prefix.
        col_prefix = "KL"        # add "KL" to column labels.
        # Offsets and font sizes come from plot_config as defined.
    else:
        # Swapped mode: rows = final_columns, columns = proteins.
        row_labels = final_columns
        col_labels = sorted_proteins
        # Now, the row labels (left side) are the K loci; add "KL" prefix there.
        row_prefix = "KL"
        col_prefix = ""          # protein labels as column headers without prefix.
        # Optionally swap the label sizes and offsets; here I simply swap the values.
        # For example, using protein_label_fontsize for column headers and k_label_fontsize for row headers:
        temp_fontsize = plot_config["protein_label_fontsize"]
        plot_config["protein_label_fontsize"] = plot_config["k_label_fontsize"]
        plot_config["k_label_fontsize"] = temp_fontsize
        # Similarly swap the offsets if desired:
        temp_offset = plot_config["protein_label_offset"]
        plot_config["protein_label_offset"] = plot_config["k_label_offset"]
        plot_config["k_label_offset"] = temp_offset

    num_rows_grid = len(row_labels)
    num_cols_grid = len(col_labels)
    
    # --- Create Visualization ---
    fig, ax = plt.subplots(figsize=plot_config["figsize"])
    
    # Draw grid cells.
    # Loop over rows and columns based on the current mode.
    for i, row_lab in enumerate(row_labels):
        # Compute y position: top row is index 0 -> y = num_rows_grid - 1 - i.
        y = num_rows_grid - 1 - i
        # Draw row label on left.
        ax.text(plot_config["protein_label_offset"], y + plot_config["cell_height"]/2,
                row_prefix + row_lab,
                ha='right', va='center', fontsize=plot_config["protein_label_fontsize"],
                fontweight='bold')
        for j, col_lab in enumerate(col_labels):
            x = j  # x position for column j.
            # In cell, if we did not swap: cell corresponds to (protein, K locus).
            # If swapped: cell corresponds to (protein, K locus) with protein = col_lab and K locus = row_lab.
            if not swap_axes:
                cell_color = cell_colors.get((row_lab, col_lab), "white")
            else:
                cell_color = cell_colors.get((col_lab, row_lab), "white")
            rect = plt.Rectangle((x, y), plot_config["cell_width"], plot_config["cell_height"],
                                 facecolor=cell_color, edgecolor="black")
            ax.add_patch(rect)
    
    # Draw column headers.
    for j, col_lab in enumerate(col_labels):
        x = j + plot_config["cell_width"]/2
        # Place header above the grid.
        ax.text(x, num_rows_grid + plot_config["k_label_offset"], 
                col_prefix + col_lab,
                ha='center', va='bottom', fontsize=plot_config["k_label_fontsize"],
                rotation=plot_config["k_label_rotation"], fontweight='bold')
    
    # Set axis limits.
    ax.set_xlim(plot_config["protein_label_offset"] - 0.5, num_cols_grid)
    ax.set_ylim(0, num_rows_grid)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.savefig(panelC_pdf_path, format='pdf', bbox_inches='tight')
    plt.close(fig)

def remove_file(path):
    from pathlib import Path
    Path(path).unlink(missing_ok=True)



if __name__ == "__main__":

    from pathlib import Path
    import yaml


    # input
    config_yaml = '../../config/config.yaml'

    # Load the entire YAML document into a Python object
    with open(config_yaml, 'r') as f:
        config_dict = yaml.safe_load(f)

    # params
    custom_gwas_kloci = config_dict['params']['custom_kl_list'].split(' ')
    custom_gwas_kloci.append('KL146')

    # I/O paths
    
    # input directories
    user_path = config_dict['paths']['janusz']['main']
    figshare_dir = Path(user_path, config_dict['paths']['janusz']['figshare_dir'])
    supplement_dir = Path(user_path, config_dict['paths']['janusz']['supplement_dir'])

    # input paths
    bacteria_table = Path(figshare_dir, 'GWAS/1_INTERMEDIATE/1_PROCESSED_INPUT/PROTMINLEN300_SPECIES-KPN-KVV-KQQ-KQS_MINNCONTIGS500/bacteria.tsv')
    S1_Table_xlsx = Path(supplement_dir, 'S1_Table.xlsx')
    S3_Table_xlsx = Path(supplement_dir, 'S3_Table.xlsx')
    S4_Table_xlsx = Path(supplement_dir, 'S4_Table.xlsx')


    # output directories
    output_dir = Path().cwd()
    suppl_figs_dir = Path(Path().cwd(), 'suppl-figs')
    tmp_blast_dir = Path(user_path, 'tmp_blast_dir')
    
    # output paths
    panelA_pdf_path = Path(output_dir, 'Figure4A.pdf')
    nodes_outfile = Path(output_dir, 'Figure4B_NODES.tsv')
    edges_outfile = Path(output_dir, 'Figure4B_EDGES.tsv')
    panelC_pdf_path = Path(suppl_figs_dir, 'recombinant-host.pdf')
    predictions_and_enzymes = Path(output_dir, 'predictions_and_enzymes.tsv')

    # create dirs
    suppl_figs_dir.mkdir(exist_ok=True, parents=True)

    # run
    generate_predictions_and_enzymes()
    FIGURE4_PANELA()
    FIGURE4_PANELB()
    FIGURE4_PANELC()

    # clean
    remove_file(predictions_and_enzymes)