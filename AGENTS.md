This repo is an implementation of the PISM (https://github.com/pism/pism/) model in .esm files.

PISM is made up of process components that are composed together into the large-scale model. Implementation will take place in three stages.

In the first stage, the original fortran codes will be instrumented to output inputs and outputs from each model component during simulations. Outputs can either be instantaneous derivatives (preferred) or integrated trajectories over a time period. The code will also be instrumented to output inputs and outputs of subassemblies of components as appropriate. A select subset of these inputs and outputs will be used as tests in .esm stub files for each component or subassembly which will be added as pull requests to https://github.com/EarthSciML/EarthSciModels.

In stage 2, The stub .esm files will be filled in with the actual physics of each process so that the tests pass, then each PR will be reviewed by a human and merged. Then, the indvidual .esm components will be imported by reference into the subassembly .esm files and coupled together so that the subassembly .esm tests pass.

In stage 3, the components that have been merged into EarthSciModels will be imported by reference and coupled into top-level eqice.esm model and tested against simulations of the full fortran model.

All model logic should be contained in .esm files. Tests and examples should be in the "tests" and "analysis" sections of .esm files, respectively. All local tests should be conducted using the EarthSciAST rust CLI binary, a copy of which should be kept untracked in the root directory of this repo.

The .esm files should be authored compositionally, with separate components, subcomponents, and expression templates used liberally and imported by reference into other components to keep the .esm code succint, simple, and human-interpretable. Avoid repeating logic or calculations, instead factor reused pieces into their own files/components/expression templates and import them by reference. Avoid the use of scripts to mechanically generate .esm files, as they often end up producing expressios that are not properly factored or succinct. (Scripts MAY be used to extract numeric test tuples from Fortran dumps into an existing hand-authored `tests` block, and to translate tabular data such as KPP `.eqn` reaction lists or lookup tables into `reaction_systems`, `function_tables`, or `data_sources` entries; the equations, structure, and factoring of a component are always authored by hand.)

Always, the equations in the .esm files should like like math rather than einsum operations whenever possible. For example, if a component includes a aggregate operation that is a numerical representation of a mathematical integral, the .esm equivalent should be written as an integral instead of a aggregate operation. Likewise, PDEs that are discretized in the fortran could should be written as PDEs, with discretization templates either used from or submitted as pull requests to https://github.com/EarthSciML/EarthSciDiscretizations.

Some components may already have implementations in EarthSciModels. If so, add additional tests from the instrumented FORTRAN to ensure the existing .esm implementations match their fortran counterparts.

Keep a centralized list of any bugs or errors that are found in the fortran code.

Do not store anything large in /tmp, as it is backed by RAM rather than hard drive and can cause an OOM.

Do not put "Co-Authored-By: Claude ..." in any commit messages. Claude is an LLM and LLMs cannot take responsibility for outputs, therefore they cannot be authors.

