# Changelog

## [v0.2.0] - Initial beta release

This release marks the first public changelog entry for DePSI and captures the current state of the Python implementation of the DePSI workflow. It follows the early v0.1.x prototype phase and introduces a more structured, maintainable package for Persistent Scatterer Interferometric SAR (PS-InSAR) processing.

### Notes
- DePSI is still a beta project and is under active development.
- The codebase is intended to be tested and documented further as development continues.
- This release focuses on establishing the Python implementation as a usable, extensible foundation rather than a complete production-grade feature freeze.

### Added
- Initial Python package release of DePSI for PS-InSAR processing.
- End-to-end workflow support for:
  - PS candidate selection
  - network generation and deformation estimation
  - atmospheric phase screen (APS) estimation
  - densification of PS points
  - geocoding of deformation estimates
- Modernized Python-based implementation built around the scientific Python ecosystem.
- PyPI packaging and installation workflow for end users.
- Documentation and usage guides covering installation, workflow steps, and API reference.
- Developer-oriented project structure for future maintenance and contribution.

### Changed
- Ported the original MATLAB-based DePSI workflow into a Python software package to improve maintainability and extensibility.
- Reorganized the project around a reusable package structure and documentation-first workflow.
- Improved accessibility for a broader InSAR and research community through modern tooling and packaging.


### Previous versions
- v0.1.x: early prototype / pre-release work; not yet a stable package baseline.