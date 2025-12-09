# PHY 105N Engineering Physics II Laboratory
## Fall 2025

This repository contains lab materials, analysis scripts, and figures for the introductory physics laboratory course **PHY 105N**. The code is primarily MATLAB and is used to analyze data collected during weekly labs and for the final project.

## Repository Structure

- `Lab 0/` – Introductory lab on measurement and uncertainty (report in `.pages` format).
- `Lab 1/` – Focus on uncertainty and basic data analysis (report in `.pages` format).
- `Lab 2/` – Lab 2 materials (folder placeholder; contents may be external or on Canvas).
- `Lab 3/` – Lab 3 materials (folder placeholder; contents may be external or on Canvas).
- `Lab 4/` – Diffraction and single‑slit analysis scripts and figures (e.g. `lab4PartI.m`, `lab4PartII.m`, and saved `.fig` plots).
- `Lab 5/` – Linear and quadratic model fits for lightbulb and rheostat experiments (`linearModel_Lightbulb.m`, `quadModelA_Lightbulb.m`, `newLinModel_Rheostat.m`, etc.).
- `Lab 6/` – Modeling scripts for log–log and linearized relationships (`PartI_model.m`, `PartII_model.m`) and associated figures.
- `Lab 7/` – Statistical testing scripts (e.g. `testingStats.m`).
- `Lab 8/` – Lab 8 materials (folder placeholder; contents may be external or on Canvas).
- `Final Project/` – Final project analysis and figures for the capacitor experiment (`graphsCapacitance.m`, `capVsPlatArea.fig`, `capVsPlatSep.fig`).
- `Prelab 1/` – Pre‑lab materials and figures (`prelab1.fig`).
- Top‑level `.pages` files – Written lab reports and project documentation in Apple Pages format.

## Requirements

- **MATLAB** (R2020b or later is recommended, but most scripts should work in older versions that support `readtable`, plotting, and basic statistics functions).
- Ability to open **MATLAB `.fig` files** (MATLAB desktop recommended).
- For the final project script `graphsCapacitance.m`, the Excel file `CapacitorFinalProjectDataSet.xlsx` must be available in the same directory as the script.

## How to Use the MATLAB Scripts

1. Open MATLAB and set the current folder to the appropriate lab directory (e.g. `Lab 4/`, `Lab 5/`, or `Final Project/`).
2. Make sure any required data files (e.g. Excel spreadsheets or CSV files referenced in the scripts) are present in that folder.
3. Run the desired script from the MATLAB Command Window, for example:
   ```matlab
   lab4PartI
   linearModel_Lightbulb
   graphsCapacitance
   ```
4. Review the generated figures and command‑window output (e.g. best‑fit parameters, chi‑square statistics, and uncertainty estimates) as part of your lab analysis.

## Notes

- Scripts are intended as analysis tools for specific PHY 105N lab activities; they are not general‑purpose libraries.
- Many scripts assume data values are hard‑coded inside the file (as given in the lab handouts). To reuse them with new data, edit the vectors at the top of each script accordingly.
- This repository does not include grading information or instructor solutions; it only contains the student work and analysis code used in the course.
