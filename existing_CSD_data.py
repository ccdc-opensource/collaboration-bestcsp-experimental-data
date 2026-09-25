import ccdc
import sys
import os

entries = ccdc.io.EntryReader('CSD')
MAX_REFCODE_SUFFIX = 100


def retrieve_crystal_data(family):
    entries = ccdc.io.EntryReader('CSD')
    volumes_by_polymorph = {}
    temperature_by_polymorph = {}
    dois_by_polymorph = {}
    r_factor_by_polymorph = {}
    for i in range(MAX_REFCODE_SUFFIX):
        if i == 0:
            num = ""
        else:
            num = "%02i" % i
        full_name = f"{family}{num}"
#        print(full_name)
        try:
            entry = entries.entry(full_name)
        except Exception as e:
#            print(f"Could not retrieve crystal for {full_name}: {e}")
            continue
        polymorph = entry.polymorph if entry.polymorph else "Unknown"
        polymorph = polymorph.strip()
        polymorph = polymorph.replace("polymorph","")
        polymorph = polymorph.replace("/","")
        volume = entry.crystal.cell_volume
        if volume is not None:
            if polymorph not in volumes_by_polymorph:
                volumes_by_polymorph[polymorph] = []
            volumes_by_polymorph[polymorph].append(volume)
        temperature = entry.temperature if entry.temperature else "Unknown"
        if temperature is not None:
            if polymorph not in temperature_by_polymorph:
                temperature_by_polymorph[polymorph] = []
            temperature_by_polymorph[polymorph].append(temperature)
        doi = entry.publication.doi if entry.publication and entry.publication.doi else full_name
        if polymorph not in dois_by_polymorph:
            dois_by_polymorph[polymorph] = []
        dois_by_polymorph[polymorph].append(doi)
        r_factor = entry.r_factor if entry.r_factor else "Unknown"
        if polymorph not in r_factor_by_polymorph:
            r_factor_by_polymorph[polymorph] = []
        r_factor_by_polymorph[polymorph].append([full_name,r_factor])
    return volumes_by_polymorph, temperature_by_polymorph, dois_by_polymorph, r_factor_by_polymorph


def write_volume_csv(common_name, volumes_by_polymorph, temperature_by_polymorph, dois_by_polymorph, r_factor_by_polymorph):
    for polymorph, volumes in volumes_by_polymorph.items():
        output_file = open(os.path.join(common_name, f"{polymorph.strip(' ')}_volume.csv"), 'w')
        lowest_r_factor_rep = r_factor_by_polymorph[polymorph][0][0]
        min_r_factor = 100.0
        for full_name, r in r_factor_by_polymorph[polymorph]:
            if r != 'Unknown':
                if float(r) < min_r_factor:
                    min_r_factor = float(r)
                    lowest_r_factor_rep = full_name
        print(f'{polymorph} {lowest_r_factor_rep} {min_r_factor}')
        output_file.write("%s\n" % (lowest_r_factor_rep))  # should be best R factor refcode rather than first really
        output_file.write("Identifier, Property, Value (Angstrom^3), Std, N, Name, Reference, Comment\n")
        for i, volume in enumerate(volumes):
            output_file.write(
                f"{i+1}, Unit Cell Volume @  {temperature_by_polymorph[polymorph][i]}, {volume:.2f}, , 1, ,{dois_by_polymorph[polymorph][i]}, generated automatically; check temperature etc.\n")


def main():
    family = sys.argv[1]
    common_name = sys.argv[2]
    volumes_by_polymorph, temperature_by_polymorph, dois_by_polymorph, r_factor_by_polymorph = retrieve_crystal_data(family)
    if os.path.exists(common_name):
        print(f"folder {common_name} exists.")
    else:
        print(f"making folder {common_name}")
        os.mkdir(common_name)
    write_volume_csv(common_name, volumes_by_polymorph, temperature_by_polymorph, dois_by_polymorph, r_factor_by_polymorph)
#    for polymorph, r_factors in r_factor_by_polymorph.items():
#        for full_name, r in r_factors:
#            print(f"{polymorph}: {full_name} - {r}")

if __name__ == "__main__":
    main()
