"""Script to check whether the focal window output files are complete.

Run script with python:
$ python src/steps/20_FocalLULC/show_files.py path/to/focal_output base_name

The output folder should have the following structure,
first it goes down the base_name, but each _ is a new folder:

```
path/to/focal_output/base/name/ (9 year folders, 2020, 2025, ..., 2060)
├── {year}/ (1080 sim folders)
│   ├── {simID}/ (11 reg folders)
│   │   ├── reg1/ (5 files each)
│   │   │   ├── base_name_{year}_{simID}_reg1_100.rds
│   │   │   ├── base_name_{year}_{simID}_reg1_200.rds
│   │   │   ├── base_name_{year}_{simID}_reg1_500.rds
│   │   │   ├── base_name_{year}_{simID}_reg1_1500.rds
│   │   │   └── base_name_{year}_{simID}_reg1_3000.rds
│   │   ├── reg10/ (5 files each)
│   │   ├── ...
│   │   └── reg19/
│   ├── ...
│   └── {simID}/
├── ...
└── {year}/
```

In total, there should be {n_years} * {n_sims} * {n_regs} * {n_files} files.
9 * 1080 * 11 * 5 = 534600 files (59400 for each year).

The script will return the percentage of finished files for each year
in a table.
Additionally, it will return the filenames of the missing files and
save them to a text file called missing_files.txt.
"""

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import itertools


class ProgressCounter:
    def __init__(self, total):
        self.counter = itertools.count()
        self.total = total

    def increment(self):
        return next(self.counter)


def print_loading_bar(progress, total, bar_length=40):
    percent = float(progress) / total
    arrow = '-' * int(round(percent * bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f"\r[{arrow}{spaces}] {percent * 100:.2f}%")
    sys.stdout.flush()


def check_year(output_dir, base_name, year, scenario_id_list, reg_list,
               window_size_list, progress_counter):
    """Check the output files for a single year."""
    num_files = 0
    missing_files = []
    for sim_id in scenario_id_list:
        for reg in reg_list:
            for window_size in window_size_list:
                file_name = f"{base_name}_{year}_{sim_id}_{reg}_{window_size}.rds"
                file_path = os.path.join(output_dir,
                                         base_name.replace("_", "/"), str(year),
                                         str(sim_id), reg, file_name)
                if os.path.isfile(file_path):
                    num_files += 1
                else:
                    missing_files.append(file_path)
        progress = progress_counter.increment()
        print_loading_bar(progress, progress_counter.total)

    return {
        "year": year,
        "covered_files": num_files,
        "covered_percent": num_files / (
                    len(scenario_id_list) * len(reg_list) * len(
                window_size_list)) * 100,
        "missing_files": missing_files,
    }


def check_files(output_dir, base_name, year_list, scenario_id_list, reg_list,
                window_size_list):
    """Check the output files for the focal LULC model.

    Show table of coverage/completion for each year.
    Save missing files to text files.
    """
    results_by_year = {}
    total_tasks = len(year_list) * len(scenario_id_list)
    progress_counter = ProgressCounter(total_tasks)

    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(check_year, output_dir, base_name, year,
                                   scenario_id_list, reg_list, window_size_list,
                                   progress_counter): year for year in
                   year_list}
        for future in as_completed(futures):
            year = futures[future]
            try:
                results = future.result()
                results_by_year[year] = results
            except Exception as exc:
                print(f"Year {year} generated an exception: {exc}")

    # Print table of coverage for each year
    print("\nYearly coverage:")
    print("Year".ljust(6), "Files".ljust(10), "Percent")
    print("-" * 25)
    for year, results in results_by_year.items():
        print(str(year).ljust(6), str(results["covered_files"]).ljust(10),
              f"{results['covered_percent']:.2f}%")
    print("-" * 25)
    total_files = len(year_list) * len(scenario_id_list) * len(reg_list) * len(
        window_size_list)
    total_covered = sum(
        results["covered_files"] for results in results_by_year.values())
    print("Total".ljust(6), str(total_covered).ljust(10),
          f"{total_covered / total_files * 100:.2f}%")
    print("-" * 25)

    # Save missing files to text file
    missing_files = []
    for year, results in results_by_year.items():
        missing_files.extend(results["missing_files"])

    if missing_files:
        with open("missing_files.txt", "w") as f:
            f.write("\n".join(missing_files))
        print(f"Missing files written to missing_files.txt.")
    else:
        print("No missing files found.")


def main():
    if len(sys.argv) != 3:
        print("Usage: python show_files.py path/to/focal_output base_name")
        sys.exit(1)
    output_dir = sys.argv[1]
    if not os.path.isdir(output_dir):
        print(f"Error: {output_dir} is not a directory.")
        sys.exit(1)
    base_name = sys.argv[2]
    year_list = list(range(2020, 2061, 5))
    scenario_id_list = list(range(1, 1081))
    reg_list = [f"reg{i}" for i in [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]
    window_size_list = [100, 200, 500, 1500, 3000]

    check_files(output_dir, base_name, year_list, scenario_id_list, reg_list,
                window_size_list)


if __name__ == "__main__":
    main()
