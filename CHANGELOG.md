## v0.3.0 - 30 Sep 2026
- Extensive remodeling of the class layout
    - **This is a breaking change**: old workflows in v0.2 basis will not work without modifications. The backward compatibility will be enhanced in the future, as the pickle serialization is removed.
- Full refactorization of the file I/O
    - Remove 'dill' dependency in favor of minimal ASCII serialization
    - Implement I/O caching features for all central classes, like ParameterSet, ParameterHessian, etc
    - Consistent, simplified syntax to treat I/O using 'path' keyword
    - Use pathlib.Path for all file paths
    - Streamline, tweak code and examples
- Full refactorization of the PesFunction
    - Use FunctionCaller class for treating all callable operations with arguments
    - Derived classes (e.g. FilesFunction, NexusPes) only override simple hooks
    - Add syntactic sugar by making PesFunction callable
    - Streamline class layout
    - Updated treatment of warnings, exceptions
- Comphrehensive renaming and reorganization of package files
- Revise std I/O and add logger features
- Add TransitionPathWay features, including example
- Update Nexus integration to support new packaging
- Update Nexus features, e.g.
    - Bundling of jobs
    - Dependent jobs
    - Robustness updates to specific loaders
- Enhance Parameter-derived classes, printout
    - Add BondAngle, PhaseAngle
    - Update examples in this regard
- Add EffectiveVarianceMap and analysis of apparent error to characterize 'black-box' statistical properties
- Add new fitting classes: MorseFit, SplineFit
- Add qiskit VQE examples
- Revised documentation and added first tutorials + examples
- Revise usability functions, like plotting, printing
- Addition of type hints, numerous fixes to typing checks
- Numerous minor improvements, tweaks and bug fixes

## v0.2.1 - 28 May 2025
- Error surface optimization: Performance update and minor bug fixes
- Plotting and printing updates for line-searches
- Refactor examples and add README for documentation
    - Add Nexus examples: carbon dimer, carbon diamond
    - Add PySCF examples: benzene, H2O
- Add Nexus job bundle feature and improved prompt
- Minor tweaks and fixes

## v0.2.0 - 15 Apr 2025

- Comprehensive code refactorization, as overviewed in the following:
- PyPI support
- Isolation and streamlining of Nexus-related functionalities
    - Nexus support is strictly in stalk.nexus module and only available with Nexus
    - The code also works without nexus unless Nexus features are explicitly requested
- Added wrapper classes for core functionalities, to enable checks and enforce API:
    - e.g. PesFunction, NexusGenerator, FittinfFunction, GeometryLoader etc
- Functionality updates in line-search and optimizer, including:
    - Automatic ordering of grid, and omission of grid points
    - Target bias bracketing for improved optimizer performance
- Streamlining of API and script usage, e.g.
    - Object pickling/unpickling
    - Nexus job control; no more need to specify mode of operation
    - Surrogate optimization
- Increased use of @property for better property control in classes and to avoid redundancy
of numerical properties
- More diverse class inheritance chains and type hinting
- Increased and more straightforward unit test coverage
- Revised printing and plotting features and their inheritance
- Additional examples and PySCF support
- Documentation updates

## Rebranding – 19 Dec 2024

- The repository was renamed 'surrogate_hessian_relax'->'stalk' on 19 Dec 2024, and
the code usage has changed substantially upon python packaging. To complete projects in the
old code base, keep using [v0.1](https://github.com/QMCPACK/stalk/releases/tag/v0.1) or
reach out for help in migration.

## v0.1 - 11 Nov 2024 and earlier

- A development version used in original or modified variants in selected research projects.