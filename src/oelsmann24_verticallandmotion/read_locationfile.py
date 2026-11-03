import numpy as np

""" read_locationfile.py

Reads in the location file in order to get site names, site ids, lats, and lons

Parameters:
location_file = Background Rate file

"""


def ReadLocationFile(location_file):
    # Initialize variables to hold data and site information
    names = []
    ids = []
    lats = []
    lons = []

    # Open the rate file
    with open(location_file, "r") as f:
        # Loop over the lines of the file
        for number, line in enumerate(f, 1):
            # Skip commented lines
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            # Split the line into components
            try:
                (this_name, this_id, this_lat, this_lon) = line.strip().split("\t")
                this_id, this_lat, this_lon = (
                    int(this_id),
                    float(this_lat),
                    float(this_lon),
                )
                if not this_name or this_id in ids:
                    raise ValueError("Site names must be nonempty and IDs unique")
                if (
                    not np.isfinite(this_lat)
                    or abs(this_lat) > 90
                    or not np.isfinite(this_lon)
                ):
                    raise ValueError(
                        "Coordinates must be finite and latitude in -90..90"
                    )
                if not -180 <= this_lon <= 360:
                    raise ValueError("Longitude must use -180..180 or 0..360 degrees")
            except ValueError as exc:
                raise ValueError(f"{location_file}: line {number}: {exc}") from exc

            # Store the information
            names.append(this_name)
            ids.append(int(this_id))
            lats.append(float(this_lat))
            lons.append(float(this_lon))

    if not ids:
        raise ValueError(f"{location_file}: no locations found")
    # Cast everything as numpy arrays
    names = np.array(names)
    ids = np.array(ids)
    lats = np.array(lats)
    lons = np.array(lons)

    # Return variables
    return (names, ids, lats, lons)


if __name__ == "__main__":
    import argparse

    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(
        description="Reads in a location file used for localization in FACTS",
        epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)",
    )

    # Define the command line arguments to be expected
    parser.add_argument("infile", help="File from which to read input data")

    # Parse the arguments
    args = parser.parse_args()

    # Load the locations from the file
    (names, ids, lats, lons) = ReadLocationFile(args.infile)

    # Print the info to the screen
    print("--- infile: {} ---".format(args.infile))
    for i in np.arange(len(ids)):
        print(
            "Name: {}\tID: {}\tLat: {}\tLon: {}".format(
                names[i], ids[i], lats[i], lons[i]
            )
        )

    # Done
    exit()
