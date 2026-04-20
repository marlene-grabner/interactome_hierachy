import mygene


def dict_ncbi_to_hgnc(ncbi_names: list):
    """
    This function takes a list of NCBI gene IDs and returns a dictionary
    mapping those IDs to HGNC symbols using the MyGeneInfo API.
    """

    # One query to MyGene
    mg = mygene.MyGeneInfo()
    query_results = mg.querymany(
        ncbi_names,
        scopes="entrezgene",
        fields="symbol",
        species="human",
        verbose=True,
    )

    # Constrcut lookup dictionary
    ncbi_to_hgnc = {}
    for match in query_results:
        if "query" in match and "symbol" in match:
            ncbi_to_hgnc[match["query"]] = match["symbol"]

    # Inform user of how successful the mapping was
    mapped_count = len(ncbi_to_hgnc)
    total_count = len(ncbi_names)
    print(
        f"Successfully mapped {mapped_count} out of {total_count} NCBI IDs to HGNC symbols."
    )

    return ncbi_to_hgnc
