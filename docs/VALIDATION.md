# Validation

On 23 November 2026, every one of the release's 50 workflow stages completed successfully in a validation checkout initialized without generated data, numerical results, figures, or tables. Previously downloaded public inputs were reused after checksum verification. The clinical API request and response were retained for that run. Affected stages were rerun after final packaging and figure corrections.

The output-coverage check passed for seven main figures, eight supplementary figures, and Tables S1–S3. The three workbooks contain 72 labelled scientific sections, with numerical source paths and SHA-256 checksums. Coverage means the declared reconstructed outputs were produced; it does not assert equality with an earlier manuscript version.

Thirteen implementation tests passed, covering:

- independent Breslow Cox estimates and standard errors;
- event-code validation and BH multiplicity families;
- count-library/annotation identity;
- TRACERx patient-score scaling and independent HC3 variance calculation;
- explicit joint-model Rest3 units;
- matched patient masks for proliferation models;
- checksum verification and explicit primary-diagnosis selection.

The reconstructed ClearCode34 classifier agreed with 376 of 380 published benchmark calls (98.9%). Its inputs and benchmark labels are downloaded from public sources, not retained model outputs. EPIC convergence diagnostics are exported, and only converged samples enter the adjusted models.

Validation used Linux, Python 3.12, R 4.4, and the package versions in `requirements.txt`. Figures were visually checked for labels, legends, and layout. No test hard-codes a desired manuscript association. The reconstructed spatial and phosphosite definitions and their interpretation limits are documented in [Methods](METHODS.md).

To repeat the implementation and output checks after computation:

```bash
python -m unittest discover -s tests -v
python scripts/reproduce.py --check-coverage
```
