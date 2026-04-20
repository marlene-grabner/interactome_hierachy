import pandas as pd


def summarize_level_enrichment(
    cluster_enrichment_dict: dict, level_name: str, output_folder: str = None
):
    """
    Takes a dictionary where keys are cluster IDs and values are Enrichr DataFrames.
    Returns a single consolidated DataFrame.
    """
    summary_dfs = []

    for cluster_id, df in cluster_enrichment_dict.items():
        # 1. Sort by significance
        df_sorted = df.sort_values("Adjusted P-value")

        # 2. Take the top 5 rows
        top_5 = df_sorted.head(5).copy()

        # 3. Keep only essential columns
        cols_to_keep = ["Term", "Adjusted P-value", "Combined Score"]
        top_5 = top_5[cols_to_keep]

        # 4. Add the cluster ID
        top_5.insert(0, "Cluster_ID", cluster_id)

        summary_dfs.append(top_5)

    # Combine all clusters for this level into one DataFrame
    level_summary = pd.concat(summary_dfs, ignore_index=True)

    # Save this single file if specified
    if output_folder is not None:
        level_summary.to_csv(
            f"{output_folder}/{level_name}_enrichment_summary.csv", index=False
        )

    return level_summary
