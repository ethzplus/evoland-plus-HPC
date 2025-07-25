"""Script to show existing NCP results by scenario

Run script from conda 'ncps' environment:
(ncps) $ python show_files.py

1. Read numbers from files to give back histogram of counts.

File format:
```
1000/: 216
1001/: 243
100/: 219
1002/: 216
1003/: 216
1004/: 216
1005/: 216
```

Extraxt the numbers and give back histogram of counts.

2. Check if each expected file is present in each scenario.

Files to expect:
```
{simID}/ (216 files)
├── run_ncp_params_{year}.yml (9 files)
├── CAR (9 files)
│   ├── tot_c_cur_2020.tif
│   ├── ...
│   └── tot_c_cur_2060.tif
├── FF (9 files)
│   ├── FF_S_CH_2020.tif
│   ├── ...
│   └── FF_S_CH_2060.tif
├── HAB (9 folders, 18 files)
│   ├── 2020
│   │   ├── deg_sum_c.tif
│   │   └── quality_c.tif
│   ├── ...
│   └── 2060
├── NDR (9 folders, 9 files)
│   ├── 2020
│   │   └── watershed_results_ndr.gpkg
│   ├── ...
│   └── 2060
├── POL (9 files)
│   ├── POL_S_CH_2020.tif
│   ├── ...
│   └── POL_S_CH_2060.tif
├── SED (9 folders, 36 files)
│   ├── 2020
│   │   ├── watershed_results_sdr.shp
│   │   ├── watershed_results_sdr.shx
│   │   ├── watershed_results_sdr.dbf
│   │   └── watershed_results_sdr.prj
│   ├── ...
│   └── 2060
└── WY (9 folders, 13*9=177 files)
    ├── 2020
    │   ├── watershed_results_wyield.shp, .shx, .dbf, .prj, .csv
    │   ├── subwatershed_results_wyield.shp, .shx, .dbf, .prj, .csv
    │   └── per_pixel
    │       ├── aet.tif
    │       ├── fractp.tif
    │       └── wyield.tif
    ├── ...
    └── 2060
```

"""

import sys
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

output_files = [
    "run_ncp_params_{year}.yml",
    "CAR/tot_c_cur_{year}.tif",
    "FF/FF_S_CH_{year}.tif",
    "HAB/{year}/deg_sum_c.tif",
    "HAB/{year}/quality_c.tif",
    "NDR/{year}/n_subsurface_export.tif",
    "NDR/{year}/n_surface_export.tif",
    "NDR/{year}/n_total_export.tif",
    "NDR/{year}/p_surface_export.tif",
    "POL/POL_S_CH_{year}.tif",
    "REC/REC_S_CH_{year}.tif",
    "SDR/{year}/avoided_erosion.tif",
    "SDR/{year}/avoided_export.tif",
    "SDR/{year}/sed_deposition.tif",
    "SDR/{year}/sed_export.tif",
    "SDR/{year}/stream.tif",
    "SDR/{year}/rkls.tif",
    "SDR/{year}/usle.tif",
    "WY/{year}/watershed_results_wyield.shp",
    "WY/{year}/watershed_results_wyield.shx",
    "WY/{year}/watershed_results_wyield.dbf",
    "WY/{year}/watershed_results_wyield.prj",
    "WY/{year}/watershed_results_wyield.csv",
    "WY/{year}/subwatershed_results_wyield.shp",
    "WY/{year}/subwatershed_results_wyield.shx",
    "WY/{year}/subwatershed_results_wyield.dbf",
    "WY/{year}/subwatershed_results_wyield.prj",
    "WY/{year}/subwatershed_results_wyield.csv",
    "WY/{year}/per_pixel/aet.tif",
    "WY/{year}/per_pixel/fractp.tif",
    "WY/{year}/per_pixel/wyield.tif",
]

output_files_complete = {
    f: [f.format(year=year) for year in range(2020, 2061, 5)]
    for f in output_files
}


def print_loading_bar(progress, total, bar_length=40):
    percent = float(progress) / total
    arrow = '-' * int(round(percent * bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f"\r[{arrow}{spaces}] {percent * 100:.2f}%")
    sys.stdout.flush()


def process_scenario(scenario_path):
    found_files = set()
    with os.scandir(scenario_path) as it:
        for entry in it:
            if entry.is_file():
                found_files.add(os.path.relpath(entry.path, scenario_path))
            elif entry.is_dir():
                for root, _, files in os.walk(entry.path):
                    for file in files:
                        found_files.add(
                            os.path.relpath(os.path.join(root, file),
                                            scenario_path))
    return found_files


def get_expected_files():
    expected_files = set()
    for paths in output_files_complete.values():
        expected_files.update(paths)
    return expected_files


def collect_unexpected_files(scenario_path, expected_files):
    unexpected_files = []
    with os.scandir(scenario_path) as it:
        for entry in it:
            if entry.is_file():
                rel_path = os.path.relpath(entry.path, scenario_path)
                if rel_path not in expected_files:
                    unexpected_files.append(entry.path)
            elif entry.is_dir():
                for root, _, files in os.walk(entry.path):
                    for file in files:
                        rel_path = os.path.relpath(os.path.join(root, file),
                                                   scenario_path)
                        if rel_path not in expected_files:
                            unexpected_files.append(os.path.join(root, file))
    return unexpected_files


def check_files(scenario_dir):
    coverage = {file: 0 for file in output_files}
    missing_sim_ids = {file: [] for file in output_files}
    unexpected_files = []
    expected_files = get_expected_files()
    scenarios = list(os.scandir(scenario_dir))
    scenarios = [scenario for scenario in scenarios if scenario.is_dir()]
    total_scenarios = len(scenarios)

    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_scenario, os.path.join(scenario_dir,
                                                                  scenario.name)): scenario
                   for scenario in scenarios}
        for i, future in enumerate(as_completed(futures)):
            scenario = futures[future]
            found_files = future.result()
            unexpected_files.extend(collect_unexpected_files(
                os.path.join(scenario_dir, scenario.name), expected_files))
            for file, paths in output_files_complete.items():
                for path in paths:
                    if path in found_files:
                        coverage[file] += 1
                    else:
                        missing_sim_ids[file].append(scenario.name)
            print_loading_bar(i + 1, total_scenarios)

    # Convert to a fraction
    for file in coverage:
        coverage[file] /= total_scenarios * len(range(2020, 2061, 5))

    # Print results in a table
    print("\n\nCoverage:")
    print("File".ljust(41), "Coverage".ljust(10), "Missing in")
    print("-" * 60)
    for file, cov in coverage.items():
        # show maximal 10 missing scenarios, indicate more
        print(
            file.ljust(41),
            f"{cov * 100:.2f}%".ljust(10),
            ", ".join(list(set(missing_sim_ids[file]))[:10]) + (
                ", ..." if len(set(missing_sim_ids[file])) > 10 else "")
        )
    print("-" * 60)
    print("Total".ljust(41),
          f"{sum(coverage.values()) / len(coverage) * 100:.2f}%".ljust(10))
    print("-" * 60)

    # write paths to unexpected files to file
    if unexpected_files:
        with open("unexpected_files.txt", "w") as f:
            for path in unexpected_files:
                f.write(path + "\n")
        print(f"Paths to unexpected files written to 'unexpected_files.txt'.")
    else:
        print("No unexpected files found.")
    # write scenarios with missing files to file
    if any(missing_sim_ids.values()):
        with open("scenarios_with_missing_files.txt", "w") as f:
            for file, sim_ids in missing_sim_ids.items():
                f.write(f"{file}:\n")
                for sim_id in set(sim_ids):
                    f.write(f"  {sim_id}\n")
        print("Scenarios with missing files written to "
              "'scenarios_with_missing_files.txt'.")
    else:
        print("No scenarios with missing files found.")


def main():
    if len(sys.argv) > 2:
        print(f"Usage: {sys.argv[0]} <ncp_output_dir>")
        sys.exit(1)
    if len(sys.argv) == 2:
        check_files(sys.argv[1])
    else:
        print(f"Checking files in current directory {os.getcwd()}, "
              "should be .../ncp_output")
        check_files(os.getcwd())


if __name__ == '__main__':
    main()

# unexpected_files.txt can be used to delete unexpected files
# the straight forward method would be (sensitive to spaces in file names):
#   cat unexpected_files.txt | xargs -d '\n' rm
# This can take long, as each file is deleted sequentially.
# A faster way is to use parallel:
#   cat unexpected_files.txt | parallel rm
# This will delete multiple files in parallel, which can be faster.
# Additionally showing a loading bar can be done with:
#   cat unexpected_files.txt | parallel --bar rm
# Including -d '\n' in the parallel command is not necessary, as it is the default.
# One can specify the number of jobs to run in parallel with -j <number>. In a
# slurm job with 8 cpu cores, one could use -j 8 to use all cores.
# Choosing a higher number of jobs can be faster, but can also lead to a
# significant load on the system, which can slow down other processes.

# After deleting, there are empty directories left.
# These can be removed with:
#   find . -type d -empty -delete
# To see that the directories are deleted, one can add -print to the command:
#   find . -type d -empty -delete -print
