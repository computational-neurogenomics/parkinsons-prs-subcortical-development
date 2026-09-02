# ==============================================================================
# OLS REGRESSION, SENSITIVITY ANALYSES, AND VISUALISATION
# FOR ABCD COHORT 
# ==============================================================================

# 1. LOAD LIBRARIES
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

# ==============================================================================
# 2. USER CONFIGURATION: FILE PATHS & VARIABLES
# Note for users: Ensure your working directory is set to the root of the repository.
# ==============================================================================

# Input Data
input_file = "./data/output/abcd_input_for_ols.csv"

# Output Directories
output_dir = "./data/output/"
figures_dir = "./data/output/figures/"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

# Define Regions of Interest and Predictors
RoiCodes = ['Accumbens', 'Amygdala', 'Brainstem', 'Caudate', 'Hippocampus', 
            'Pallidum', 'Putamen', 'Thalamus', 'ventralDC', 'ICV']

# Principal Components to include (1 to 10)
PCcovarStr = '+'.join([f'PC{i}' for i in range(1, 11)])
thresholds = [f'S{i}' for i in range(1, 9)]

# Dictionary to format region names for the final plot
rename_dict = {
    'ventralDC': 'Ventral Diencephalon', 
    'ICV': 'Intracranial Volume', 
    'Accumbens': 'Nucleus accumbens', 
    'Caudate': 'Caudate Nucleus', 
    'Pallidum': 'Globus Pallidus'
}

# Threshold mapping for final supplementary tables
threshold_mapping = {
    'S1': 'P < 5E-08', 'S2': 'P < 1E-05', 'S3': 'P < 0.001', 'S4': 'P < 0.01',
    'S5': 'P < 0.05', 'S6': 'P < 0.1', 'S7': 'P < 0.5', 'S8': 'P < 1.0'
}

# ==============================================================================
# 3. HELPER FUNCTIONS FOR LINEAR MODELS
# ==============================================================================

def build_formula(roi, prs_col=None, adjust_icv=False, pc_covars="", scanner_covars=""):
    """ Safely builds the Patsy formula for statsmodels. """
    covars = f"Sex + Age + Sex:Age + I(Age**2) + Sex:I(Age**2) + {pc_covars}"
    if adjust_icv:
        covars = f"ICV + {covars}"
    if scanner_covars:
        covars = f"{covars} + {scanner_covars}"
    
    if prs_col:
        return f"{roi} ~ {prs_col} + {covars}"
    return f"{roi} ~ {covars}"

def run_sbayes(data, roi_codes, pc_covars, prs_col, scanner_covars="", adjust_icv=False):
    """ Runs OLS regressions for SBayesRC scores. """
    results = []
    for roi in roi_codes:
        red_fe = build_formula(roi, None, adjust_icv, pc_covars, scanner_covars)
        full_fe = build_formula(roi, prs_col, adjust_icv, pc_covars, scanner_covars)

        m_red = smf.ols(red_fe, data=data).fit()
        m_full = smf.ols(full_fe, data=data).fit()

        r2_marg = m_full.rsquared - m_red.rsquared

        results.append({
            'ROI': roi,
            'Beta': m_full.params[prs_col],
            'SE': m_full.bse[prs_col],
            'T-value': m_full.tvalues[prs_col],
            'P': m_full.pvalues[prs_col],
            'Rsq': r2_marg
        })
        
    df = pd.DataFrame(results).set_index('ROI')
    df['Beta (SE) [t]'] = df.apply(lambda x: f'{x["Beta"]:.2f} ({x["SE"]:.2f}) [{x["T-value"]:.2f}]', axis=1)
    df['Variance Explained (%)'] = df['Rsq'] * 100
    
    return df

def run_ct(data, roi_codes, thresholds, pc_covars, scanner_covars="", adjust_icv=True):
    """ Runs OLS regressions for C+T Polygenic Risk Scores. """
    results = []
    for thr in thresholds:
        for roi in roi_codes:
            red_fe = build_formula(roi, None, adjust_icv, pc_covars, scanner_covars)
            full_fe = build_formula(roi, thr, adjust_icv, pc_covars, scanner_covars)

            m_red = smf.ols(red_fe, data=data).fit()
            m_full = smf.ols(full_fe, data=data).fit()

            r2_marg = m_full.rsquared - m_red.rsquared

            results.append({
                'Threshold': thr,
                'ROI': roi,
                'Beta': m_full.params[thr],
                'SE': m_full.bse[thr],
                'T-value': m_full.tvalues[thr],
                'P-value': m_full.pvalues[thr],
                'Variance Explained (%)': r2_marg * 100
            })
            
    df = pd.DataFrame(results)
    df['Beta (SE) [t]'] = df.apply(lambda x: f'{x["Beta"]:.2f} ({x["SE"]:.2f}) [{x["T-value"]:.2f}]', axis=1)
    return df

# ==============================================================================
# 4. SENSITIVITY ANALYSES FUNCTIONS
# ==============================================================================

def compute_icv_improvement(data, roi_codes, pc_covars, scanner_covars):
    """ Computes the R-squared improvement when adding ICV to the base demographic model. """
    results = []
    for roi in roi_codes:
        f_wo_icv = build_formula(roi, None, False, pc_covars, scanner_covars)
        f_w_icv = build_formula(roi, None, True, pc_covars, scanner_covars)

        m_wo_icv = smf.ols(f_wo_icv, data=data).fit()
        m_w_icv = smf.ols(f_w_icv, data=data).fit()

        delta_r2 = m_w_icv.rsquared - m_wo_icv.rsquared

        results.append({
            'ROI': roi,
            'R2_without_ICV': m_wo_icv.rsquared,
            'R2_with_ICV': m_w_icv.rsquared,
            'Delta_R2': delta_r2,
            'Delta_R2 (%)': delta_r2 * 100
        })
    return pd.DataFrame(results)

def compute_scanner_sensitivity(data, roi_codes, pc_covars, prs_col, safe_scanner_cols):
    """ Evaluates the R-squared improvement and F-statistic when adding scanner IDs. """
    results = []
    scanner_term = " + ".join(safe_scanner_cols) if safe_scanner_cols else ""

    for roi in roi_codes:
        f_base = build_formula(roi, prs_col, True, pc_covars, "")
        f_scan = build_formula(roi, prs_col, True, pc_covars, scanner_term)

        m1 = smf.ols(f_base, data=data).fit()
        m2 = smf.ols(f_scan, data=data).fit()

        delta_r2 = m2.rsquared - m1.rsquared
        f_stat, p_val = np.nan, np.nan

        try:
            anova_res = anova_lm(m1, m2)
            row_idx = anova_res.index[-1]
            if 'F' in anova_res.columns: f_stat = anova_res.loc[row_idx, 'F']
            if 'Pr(>F)' in anova_res.columns: p_val = anova_res.loc[row_idx, 'Pr(>F)']
        except Exception:
            pass

        results.append({
            'ROI': roi,
            'R2_without_scanner': m1.rsquared,
            'R2_with_scanner': m2.rsquared,
            'Delta_R2': delta_r2,
            'Delta_R2 (%)': delta_r2 * 100.0,
            'F_stat': f_stat,
            'P_value': p_val
        })
        
    df = pd.DataFrame(results)
    df['P_value'] = df['P_value'].apply(lambda x: "{:.2e}".format(x) if pd.notna(x) else x)
    return df

# ==============================================================================
# 5. HEATMAP VISUALISATION FUNCTIONS
# ==============================================================================

def create_heatmap_data(df_sbayes, df_ct, roi_codes, n_corrections):
    """ Formats statistical output for heatmap rendering aligning by ROI index. """
    R2_data = pd.DataFrame(columns=roi_codes)
    P_data = pd.DataFrame(columns=roi_codes)

    R2_data.loc['SbayesRC'] = df_sbayes['Variance Explained (%)']
    P_data.loc['SbayesRC'] = df_sbayes['P']

    for thr in df_ct['Threshold'].unique():
        df_thr = df_ct[df_ct['Threshold'] == thr]
        R2_data.loc[thr] = df_thr.set_index('ROI')['Variance Explained (%)'].astype(float)
        P_data.loc[thr] = df_thr.set_index('ROI')['P-value'].astype(float)

    annotations = P_data.copy().astype(str)
    annotations[:] = ''
    
    for col in P_data.columns:
        for idx in P_data.index:
            p_val = P_data.loc[idx, col]
            if pd.notna(p_val):
                if p_val <= 0.05 / n_corrections:
                    annotations.loc[idx, col] = '**'
                elif p_val <= 0.05:
                    annotations.loc[idx, col] = '*'

    return R2_data, annotations

def plot_heatmap(data, annotations, filename_base, vmin=0, vmax=2):
    """ Renders and saves heatmap visuals. """
    sns.set_theme(font_scale=1.4, style='whitegrid')
    fig, ax = plt.subplots(figsize=(10, 10))
    cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)

    sns.heatmap(
        data, cmap=cmap, ax=ax, vmin=vmin, vmax=vmax, annot=annotations, fmt="", 
        cbar_kws={'label': 'Variance explained (%)', 'orientation': 'horizontal', 'location': 'top'}
    )

    y_labels = ['SbayesRC'] + [f'$p \\leq {thr}$' for thr in ['5e-08', '1e-05', '0.001', '0.01', '0.05', '0.1', '0.5', '1.0']]
    ax.set_yticklabels(y_labels, rotation=0)
    ax.set_xlabel('')
    ax.set_ylabel('')
    
    plt.tight_layout()
    for ext in ['svg', 'pdf', 'png']:
        fig.savefig(f"{figures_dir}{filename_base}.{ext}", transparent=True, bbox_inches='tight')
    plt.close()

# ==============================================================================
# 6. EXECUTION PIPELINE
# ==============================================================================

if __name__ == "__main__":
    
    print("Loading ABCD dataset...")
    df_abcd = pd.read_csv(input_file)

    # Extract scanner covariates dynamically and drop the first one to avoid Singular Matrix Warning
    scanner_cols = [col for col in df_abcd.columns if col.startswith('ScannerID_')]
    safe_scanner_cols = scanner_cols[1:] if len(scanner_cols) > 1 else scanner_cols
    safe_scanner_term = '+'.join(safe_scanner_cols) if safe_scanner_cols else ""

    # ---------------------------------------------------------
    # RUN MAIN REGRESSIONS (Incorporating Scanner Covariates)
    # ---------------------------------------------------------
    print("Running Model 1 (without ICV adjustment, with scanner covariates)...")
    df_sbayes1 = run_sbayes(df_abcd, RoiCodes, PCcovarStr, 'PDPRS', safe_scanner_term, adjust_icv=False)
    df_ct1     = run_ct(df_abcd, RoiCodes, thresholds, PCcovarStr, safe_scanner_term, adjust_icv=False)

    print("Running Model 2 (with ICV adjustment, with scanner covariates)...")
    df_sbayes2 = run_sbayes(df_abcd, RoiCodes, PCcovarStr, 'PDPRS', safe_scanner_term, adjust_icv=True)
    df_ct2     = run_ct(df_abcd, RoiCodes, thresholds, PCcovarStr, safe_scanner_term, adjust_icv=True)

    # ---------------------------------------------------------
    # SENSITIVITY ANALYSES
    # ---------------------------------------------------------
    print("Running ICV Sensitivity Analysis...")
    df_icv_effect = compute_icv_improvement(df_abcd, RoiCodes, PCcovarStr, safe_scanner_term)

    print("Running Scanner Sensitivity ANOVA...")
    df_scanner_effects = compute_scanner_sensitivity(df_abcd, RoiCodes, PCcovarStr, 'PDPRS', safe_scanner_cols)

    # ---------------------------------------------------------
    # GENERATE HEATMAPS
    # ---------------------------------------------------------
    print("Generating Heatmap visualisations...")
    R2_data1, DFstar1 = create_heatmap_data(df_sbayes1, df_ct1, RoiCodes, n_corrections=90)
    plot_heatmap(R2_data1.rename(columns=rename_dict), DFstar1.rename(columns=rename_dict), 'Model1_ABCD_PRS_Heatmap')

    RoiCodes_noICV = [roi for roi in RoiCodes if roi != 'ICV']
    R2_data2, DFstar2 = create_heatmap_data(df_sbayes2, df_ct2, RoiCodes_noICV, n_corrections=81)
    plot_heatmap(R2_data2.rename(columns=rename_dict), DFstar2.rename(columns=rename_dict), 'Model2_ABCD_PRS_Heatmap_noICV')

    # ---------------------------------------------------------
    # EXPORT MASTER RESULTS (EXCEL & CSV)
    # ---------------------------------------------------------
    print("Exporting supplementary tables...")
    
    # Supplementary Table 1 (SBayesRC) - Separated columns for direct validation
    supplementary_table_1_export = pd.DataFrame({
        'Subcortical Structure': df_sbayes2.index,
        'Beta (SE) [t]': [f'{df_sbayes2["Beta"].iloc[i]:.2f} ({df_sbayes2["SE"].iloc[i]:.2f}) [{df_sbayes2["T-value"].iloc[i]:.2f}]' for i in range(len(df_sbayes2))],
        'P-value': [f'{df_sbayes2["P"].iloc[i]:.2e}' for i in range(len(df_sbayes2))],
        'Variance Explained (%)': [f'{df_sbayes2["Variance Explained (%)"].iloc[i]:.2f}' for i in range(len(df_sbayes2))]
    })
    supplementary_table_1_export.to_csv(f'{output_dir}Supplementary_Table1_ABCD_SBayesRC_T-value_added.csv', index=False)

    # Supplementary Table 2 (Clumping & Thresholding) - Separated columns for direct validation
    supplementary_table_2_export = pd.DataFrame({
        'Threshold': df_ct2['Threshold'].map(threshold_mapping).fillna(df_ct2['Threshold']),
        'ROI': df_ct2['ROI'],
        'Beta (SE) [t]': [f'{df_ct2["Beta"].iloc[i]:.2f} ({df_ct2["SE"].iloc[i]:.2f}) [{df_ct2["T-value"].iloc[i]:.2f}]' for i in range(len(df_ct2))],
        'P-value': [f'{df_ct2["P-value"].iloc[i]:.2e}' for i in range(len(df_ct2))],
        'Variance Explained (%)': [f'{df_ct2["Variance Explained (%)"].iloc[i]:.2f}' for i in range(len(df_ct2))]
    })
    supplementary_table_2_export.to_csv(f'{output_dir}Supplementary_Table2_ABCD_CT_T-value_added.csv', index=False)

    # Master Excel Workbook
    with pd.ExcelWriter(f"{output_dir}ABCD_PRS_results_master.xlsx", engine='openpyxl') as writer:
        df_sbayes1.to_excel(writer, sheet_name="Sbayes_Model1")
        df_ct1.to_excel(writer, sheet_name="CT_Model1", index=False)
        df_sbayes2.to_excel(writer, sheet_name="Sbayes_Model2")
        df_ct2.to_excel(writer, sheet_name="CT_Model2", index=False)
        df_icv_effect.to_excel(writer, sheet_name="ICV_Sensitivity", index=False)
        df_scanner_effects.to_excel(writer, sheet_name="Scanner_Sensitivity", index=False)

    print("\nABCD OLS regression pipeline and visualisations completed successfully.")