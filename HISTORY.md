Changelog
=========
0.1.5 (2024-05-08)
- Add "monviso analysis" tool
- Move standard monviso worklof to "monviso isoform"
- Update Docs

0.1.4 (2024-04-10)
- Fix: Correct residue number at chain breaks (for more cases)

0.1.3 (2024-04-08)
- Fix: Removing sorting of models in Modeller script
  makes the software compatible with the Modeller's
  version in conda

0.1.2 (2024-04-08)
-----------------
- Added Containerfile for Docker
- Fix: Removal of maximum 10 models reported


0.1.1 (2024-04-02)
------------------
- Fix: The model now is assigned to chain A
- Fix: Correct residue number at chain breaks
- Fix: Alignment, scoring, and filtering is repeated as long as
       the list of templates keeps changing
- Fix: If mutation can be mapped on sequence, but it is not
       covered by templates, it gets reported in the report.md file


