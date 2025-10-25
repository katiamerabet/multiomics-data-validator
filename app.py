import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from io import StringIO, BytesIO
import base64
from typing import Optional, Dict, Any

st.set_page_config(layout="wide", page_title="Multi-Omics Data Validator")

# -----------------------
# Utilities & Loaders
# -----------------------
@st.cache_data
def load_data(file) -> Optional[pd.DataFrame]:
    try:
        name = getattr(file, "name", "")
        if name.endswith('.csv'):
            return pd.read_csv(file)
        elif name.endswith('.tsv'):
            return pd.read_csv(file, sep='\t')
        elif name.endswith('.xlsx') or name.endswith('.xls'):
            return pd.read_excel(file)
        else:
            # last attempt: try csv read
            try:
                file.seek(0)
            except Exception:
                pass
            raise ValueError("Unsupported file format. Please upload CSV, TSV, or Excel.")
    except Exception as e:
        st.error(f"Error reading file: {e}")
        return None

def df_to_bytes_csv(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode('utf-8')

def dict_to_bytes_json(d: Dict[str, Any]) -> bytes:
    return json.dumps(d, indent=2, default=str).encode('utf-8')

# -----------------------
# Detection & Validation
# -----------------------
def check_completeness(df: pd.DataFrame) -> pd.Series:
    missing_data = df.isnull().sum()
    return missing_data[missing_data > 0].sort_values(ascending=False)

def find_duplicates(df: pd.DataFrame, subset: Optional[list] = None, keep: str = 'first') -> pd.DataFrame:
    """
    If subset provided, detect duplicates with respect to subset columns.
    """
    if subset:
        valid_subset = [c for c in subset if c in df.columns]
        return df[df.duplicated(subset=valid_subset, keep=keep)]
    return df[df.duplicated(keep=keep)]

def detect_outliers(df: pd.DataFrame, method: str = 'zscore', z_thresh: float = 3.0, iqr_multiplier: float = 1.5):
    numeric = df.select_dtypes(include=[np.number]).copy()
    if numeric.empty:
        return pd.DataFrame(columns=df.columns)  # no numeric -> no outliers

    outlier_mask = pd.DataFrame(False, index=df.index, columns=numeric.columns)

    if method == 'zscore':
        # robustly compute z-scores; handle zero std
        for col in numeric.columns:
            col_data = numeric[col].dropna()
            mean = col_data.mean()
            std = col_data.std()
            if pd.isna(std) or std == 0:
                outlier_mask[col] = False
                continue
            z = (numeric[col] - mean) / std
            outlier_mask[col] = z.abs() > z_thresh
    elif method == 'iqr':
        for col in numeric.columns:
            col_data = numeric[col].dropna()
            q1 = col_data.quantile(0.25)
            q3 = col_data.quantile(0.75)
            iqr = q3 - q1
            if pd.isna(iqr) or iqr == 0:
                outlier_mask[col] = False
                continue
            lower = q1 - iqr_multiplier * iqr
            upper = q3 + iqr_multiplier * iqr
            outlier_mask[col] = (numeric[col] < lower) | (numeric[col] > upper)
    else:
        raise ValueError("Unsupported outlier method")

    flagged = df[outlier_mask.any(axis=1)].copy()
    # annotate which columns flagged
    if not flagged.empty:
        flagged['_outlier_columns'] = outlier_mask.loc[flagged.index].apply(lambda r: list(r[r].index), axis=1)
    return flagged

def check_data_consistency(df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    """
    Heuristics:
      - If column name suggests non-negative (count, abundance, level, expression, cell, reads), flag negative values.
      - If column name suggests binary/mutation, expect only 0/1 (or 0/1/NaN)
      - If column name suggests 'log' or 'fold' allow negatives
    Returns a dict of issues per column with counts and examples.
    """
    issues = {}
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        name = col.lower()
        col_series = df[col]
        col_issues = {}

        # binary / mutation checks
        if any(k in name for k in ['mutation', 'mut', 'binary', 'is_', 'has_']) or col_series.dropna().isin([0,1]).all():
            # expect 0/1 mostly
            unique_vals = sorted(col_series.dropna().unique())
            non_binary = [v for v in unique_vals if v not in (0,1)]
            if non_binary:
                col_issues['unexpected_values'] = non_binary[:10]

        # non-negative expectation heuristics
        non_negative_tokens = ['count', 'abundance', 'level', 'expression', 'cell', 'reads', 'amount']
        if any(t in name for t in non_negative_tokens) and not any(k in name for k in ['log', 'fold', 'delta', 'diff']):
            negs = col_series[col_series < 0]
            if not negs.empty:
                col_issues['negative_values_count'] = int(negs.count())
                col_issues['negative_examples'] = negs.head(3).tolist()

        if col_issues:
            issues[col] = col_issues
    return issues

def validate_categorical_data(df: pd.DataFrame, top_n: int = 10) -> Dict[str, Any]:
    # include object dtype, category dtype, and numeric columns with few unique values
    cat_cols = []
    for col in df.columns:
        if pd.api.types.is_categorical_dtype(df[col]) or df[col].dtype == object:
            cat_cols.append(col)
        else:
            # numeric but small number of unique values -> treat as categorical
            if pd.api.types.is_integer_dtype(df[col]) or pd.api.types.is_float_dtype(df[col]):
                nunique = df[col].nunique(dropna=True)
                if nunique > 0 and nunique <= min(50, max(5, int(0.05 * len(df)))):
                    cat_cols.append(col)
    issues = {}
    for col in cat_cols:
        vals = df[col].dropna().astype(str)
        freqs = vals.value_counts().head(top_n)
        issues[col] = {
            "n_unique": int(df[col].nunique(dropna=True)),
            "top_values": freqs.to_dict(),
            "rare_levels": df[col].value_counts()[df[col].value_counts() < max(2, int(0.01 * len(df)))].index.tolist()[:10]
        }
    return issues

# -----------------------
# Report
# -----------------------
def generate_report(df: pd.DataFrame,
                    missing_data: pd.Series,
                    duplicates: pd.DataFrame,
                    outliers: pd.DataFrame,
                    consistency_issues: dict,
                    categorical_issues: dict) -> str:
    report_lines = []
    report_lines.append("## Data Validation Report\n")
    report_lines.append(f"- **Number of Rows**: {df.shape[0]}")
    report_lines.append(f"- **Number of Columns**: {df.shape[1]}\n")

    report_lines.append("### 1) Missing Data")
    if not missing_data.empty:
        for col, cnt in missing_data.items():
            report_lines.append(f"- **{col}**: {int(cnt)} missing")
    else:
        report_lines.append("- No missing data detected.")

    report_lines.append("\n### 2) Duplicates")
    if not duplicates.empty:
        report_lines.append(f"- Duplicate rows found: {len(duplicates)}")
    else:
        report_lines.append("- No duplicate rows found.")

    report_lines.append("\n### 3) Outliers")
    if not outliers.empty:
        report_lines.append(f"- Outlier rows detected: {len(outliers)}")
        # show sample
        sample_out = outliers.head(5).to_dict(orient='records')
        report_lines.append(f"- Example outlier rows (first 5):\n```\n{json.dumps(sample_out, indent=2, default=str)}\n```")
    else:
        report_lines.append("- No outliers detected.")

    report_lines.append("\n### 4) Consistency Issues")
    if consistency_issues:
        for col, iss in consistency_issues.items():
            report_lines.append(f"- **{col}**: {json.dumps(iss, default=str)}")
    else:
        report_lines.append("- No consistency issues detected.")

    report_lines.append("\n### 5) Categorical Data Summary")
    if categorical_issues:
        for col, info in categorical_issues.items():
            report_lines.append(f"- **{col}**: {info['n_unique']} unique values. Top: {list(info['top_values'].keys())[:5]}")
    else:
        report_lines.append("- No categorical columns detected or no issues.")

    return "\n".join(report_lines)

# -----------------------
# App UI
# -----------------------
def main():
    st.title("🧬 Multi-Omics Data Validator")

    st.sidebar.header("Upload Dataset")
    st.sidebar.markdown("This app helps you validate the completeness, consistency, and structure of multi-omics datasets.")
    uploaded_file = st.sidebar.file_uploader("Upload your dataset (CSV, TSV, or Excel)", type=["csv", "tsv", "xlsx", "xls"])

    # Sidebar controls for thresholds and behavior
    st.sidebar.header("Validation Settings")
    outlier_method = st.sidebar.selectbox("Outlier detection method", options=['zscore', 'iqr'], index=0)
    z_thresh = st.sidebar.number_input("Z-score threshold", value=3.0, step=0.5)
    iqr_multiplier = st.sidebar.number_input("IQR multiplier", value=1.5, step=0.1, format="%.2f")
    plot_limit = st.sidebar.number_input("Max numeric columns to plot", min_value=1, max_value=50, value=10)
    duplicate_subset = st.sidebar.multiselect("Columns to consider for duplicate check (optional)", options=[], default=[])
    keep_dup_choice = st.sidebar.selectbox("Keep duplicate behavior", options=['first', 'last', False], index=0)
    if uploaded_file:
        # update duplicate_subset options once dataframe loaded (we need to be careful, so load first)
        df = load_data(uploaded_file)
        if df is None:
            st.stop()

        # Refresh duplicate options in sidebar (to show actual columns)
        # Note: this doesn't re-render the multiselect choices automatically in the same run,
        # but we will allow the user to re-run by changing a control. For usability, show columns.
        st.sidebar.write("Columns detected:")
        st.sidebar.write(list(df.columns[:50]))

        st.write("### Dataset Preview")
        st.dataframe(df.head(100))

        # Main validation flow
        has_issues = False

        with st.spinner("Running completeness check..."):
            missing_data = check_completeness(df)
            with st.expander("📋 Data Completeness Check", expanded=True):
                if not missing_data.empty:
                    has_issues = True
                    st.error("❌ Missing Data Detected!")
                    st.write(missing_data)
                    # show missingness per row summary
                    missing_per_row = df.isnull().sum(axis=1)
                    st.write("Missing values per row (top 10 rows with most missing):")
                    st.write(missing_per_row.sort_values(ascending=False).head(10))
                else:
                    st.success("✅ No missing data detected!")

        with st.spinner("Checking duplicates..."):
            # attempt to use chosen subset if it exists and columns are present
            subset_cols = [c for c in duplicate_subset if c in df.columns] if duplicate_subset else None
            duplicates = find_duplicates(df, subset=subset_cols, keep=keep_dup_choice if keep_dup_choice != False else False)
            with st.expander("🔁 Duplicate Data Identification", expanded=False):
                if not duplicates.empty:
                    has_issues = True
                    st.error("❌ Duplicate Data Found!")
                    st.write(duplicates.head(20))
                    st.write(f"Total duplicates: {len(duplicates)}")
                else:
                    st.success("✅ No duplicates detected!")

        with st.spinner("Detecting outliers..."):
            outliers = detect_outliers(df, method=outlier_method, z_thresh=z_thresh, iqr_multiplier=iqr_multiplier)
            with st.expander("📉 Outlier Detection", expanded=False):
                if not outliers.empty:
                    has_issues = True
                    st.error("❌ Outliers Detected!")
                    st.write(outliers.head(10))
                else:
                    st.success("✅ No outliers detected!")

        with st.spinner("Creating visualizations..."):
            with st.expander("📊 Data Visualizations", expanded=False):
                numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
                if numeric_cols:
                    # limit by variance ranking
                    variances = df[numeric_cols].var().sort_values(ascending=False)
                    cols_to_plot = variances.head(plot_limit).index.tolist()
                    st.write(f"Plotting top {len(cols_to_plot)} numeric columns by variance.")
                    for col in cols_to_plot:
                        fig, ax = plt.subplots(figsize=(8, 3.5))
                        sns.histplot(df[col].dropna(), kde=True, ax=ax, bins=30)
                        ax.set_title(f"Distribution of {col}", fontsize=12)
                        fig.tight_layout()
                        st.pyplot(fig)
                else:
                    st.warning("No numeric columns available for distribution plots.")

        with st.spinner("Checking consistency..."):
            consistency_issues = check_data_consistency(df)
            with st.expander("🔍 Data Consistency Check", expanded=False):
                if consistency_issues:
                    has_issues = True
                    st.error("❌ Data Consistency Issues Found!")
                    st.write(consistency_issues)
                else:
                    st.success("✅ No consistency issues detected!")

        with st.spinner("Validating categorical data..."):
            categorical_issues = validate_categorical_data(df)
            with st.expander("📊 Categorical Data Validation", expanded=False):
                if categorical_issues:
                    has_issues = True
                    st.error("❌ Issues with Categorical Data Detected!")
                    # show a compact summary table
                    for col, info in categorical_issues.items():
                        st.write(f"**{col}** — unique: {info['n_unique']}")
                        st.write(pd.DataFrame(list(info['top_values'].items()), columns=["value", "count"]).head(10))
                else:
                    st.success("✅ Categorical data is valid or none detected!")

        # Generate and show report + downloads
        with st.expander("📋 Report", expanded=True):
            report_md = generate_report(df, missing_data, duplicates, outliers, consistency_issues, categorical_issues)
            st.markdown(report_md)

            # JSON summary for programmatic use
            summary = {
                "n_rows": int(df.shape[0]),
                "n_cols": int(df.shape[1]),
                "missing": missing_data.to_dict(),
                "n_duplicates": int(len(duplicates)),
                "n_outliers": int(len(outliers)),
                "consistency_issues": consistency_issues,
                "categorical_summary": categorical_issues
            }

            # downloadable files
            st.download_button("Download report (Markdown)", data=report_md, file_name="data_validation_report.md", mime="text/markdown")
            st.download_button("Download summary (JSON)", data=dict_to_bytes_json(summary), file_name="data_validation_summary.json", mime="application/json")

            # problematic rows CSV: union of missing-rows, duplicates, outliers
            problematic_idx = set()
            problematic_idx.update(df[df.isnull().any(axis=1)].index.tolist())
            problematic_idx.update(duplicates.index.tolist())
            problematic_idx.update(outliers.index.tolist())
            if problematic_idx:
                problem_df = df.loc[sorted(problematic_idx)]
                st.write("Sample of problematic rows (first 20):")
                st.write(problem_df.head(20))
                st.download_button("Download problematic rows (CSV)", data=df_to_bytes_csv(problem_df), file_name="problematic_rows.csv", mime="text/csv")
            else:
                st.info("No problematic rows to download.")

        # Final validation result
        st.write("### Final Validation Result")
        if has_issues:
            st.error("❌ Your data is NOT valid based on the validation checks!")
        else:
            st.success("✅ Your data is VALID based on the validation checks!")

    else:
        st.info("Upload a CSV / TSV / or Excel file to begin validation.")

if __name__ == "__main__":
    main()
