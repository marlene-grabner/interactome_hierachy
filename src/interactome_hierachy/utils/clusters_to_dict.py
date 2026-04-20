def clusters_to_dict(filepath):
    cluster_dict = {}
    current_cluster = None

    with open(filepath, "r") as file:
        for line in file:
            # Remove any leading/trailing whitespace and newline characters
            line = line.strip()

            # Skip empty lines if there are any
            if not line:
                continue

            # Check if it's a header line
            if line.startswith(">"):
                # Split the line at " | " and take the first part
                # Then slice off the ">" character
                current_cluster = line.split(" | ")[0][1:]

                # Initialize a new empty list for this cluster in the dictionary
                cluster_dict[current_cluster] = []

            # If it's not a header, and we have established a current cluster, it's an ID
            elif current_cluster is not None:
                cluster_dict[current_cluster].append(line)

    return cluster_dict
