# Development Notes

This file logs internal design decisions, implementation rationale, deprecated functionality, and architectural changes.

It complements the `NEWS.md` file, which is limited to user-facing changes. The entries here may reference implementation-level changes, decisions taken during refactors, or notes for future development.

---

## [2025-05-06] [Replaced Function to Compute Gradient Vector for Interactions]

### Summary
Replaced the former C++ function responsible for updating threshold parameters in the Gibbs sampler with a new, optimized version.

### Rationale
The old implementation was functionally correct but could be optimized. The new function improves numerical efficiency, as it computes the pseudolikelihoods only once.

### Action Taken
- Integrated the new gradient function in active source files.
- Archived the old function in: src/archived code/deprecated_gibbs_functions.cpp

---

## [2025-05-03] Replaced Threshold Update Function in Gibbs Sampler

### Summary
Replaced the former C++ function responsible for updating threshold parameters in the Gibbs sampler with a new, optimized version.

### Rationale
The old implementation was functionally correct but lacked performance and modularity needed for recent algorithmic improvements. The new function improves numerical stability and better aligns with the sampler structure.

### Action Taken
- Integrated the new threshold update function in active source files.
- Archived the old function in: src/archived code/deprecated_gibbs_functions.cpp




### Notes
- Make sure any new function names or interfaces are documented in header files if reused elsewhere.
- When future updates are made to this logic, review the archived version to avoid reintroducing old assumptions.

---

## [YYYY-MM-DD] [Title of Change]
<!-- TEMPLATE FOR FUTURE ENTRIES -->

### Summary
<!-- One-sentence summary of the internal change -->

### Rationale
<!-- Why was the change made? Performance? Correctness? Maintainability? -->

### Action Taken
<!-- Describe the actual change: renamed files, swapped algorithms, archived code, etc. -->

### Notes
<!-- Optional: mention implications, caveats, things to revisit, or links to discussions/issues -->
<!-- You can also list TODOs, known limitations, or future improvements here -->

---

## Known Technical Debt / Refactor Targets
<!-- Optional section you can update over time -->

- [ ] 
- [ ] 
- [ ] 

